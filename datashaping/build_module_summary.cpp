#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

// ============================================================
// build_module_summary.cpp
//
// 目的:
//  - NaIパルス(4ch) + PSパルス(2ch)を統合し、runごとの最終出力を作る
//
// 出力1 (txt):
//  run event module integral t_fall_start_ns particle nai_ch ps_ch
//   - integral は互換のため列名を維持するが、実体は NaI 線形変換後の
//     エネルギー相当量（MeVスケール）である
//
// 出力2 (csv):
//  run,event,module,ch,time_ns,nearest_dt_same_ch_ns
//  ※ NaIに対応しなかったPSのみ
//  ※ nearest_dt_same_ch は同一run/event/ch内で最短の時間差（なければ-1）
// ============================================================

// 許容時間差 [sample]（実数サンプル差）
static const int kNaiPairWindowSamples = 5;
static const double kPsCoincidenceWindowSamples = 2.5; // 2.5 sample
static const double kSampleToTimeNs = 4.0;

// CHごとの一括時刻オフセット [sample]
static const double kOffset_NaI_A1 = 0.0;
static const double kOffset_NaI_A2 = 0.0;
static const double kOffset_NaI_B1 = 0.0;
static const double kOffset_NaI_B2 = 0.5;   // +2 ns
static const double kOffset_PS_A = -1.25;   // -5 ns
static const double kOffset_PS_B = 0.0;

// NaI積分値 -> エネルギー線形変換の既定値
// 旧 step2 の初期値を採用
static const double kNaI_A1_Coeff = 0.352608;
static const double kNaI_A2_Coeff = 0.344478;
static const double kNaI_B1_Coeff = 0.548069;
static const double kNaI_B2_Coeff = 0.342241;
static const double kNaI_A1_Offset = -558.571;
static const double kNaI_A2_Offset = -674.076;
static const double kNaI_B1_Offset = -580.215;
static const double kNaI_B2_Offset = -649.299;

struct NaiPulse {
    int run;
    long long event;
    std::string ch;   // NaI_A1, NaI_A2, ...
    std::string module; // A or B
    double integral;
    double time;
};

struct PsPulse {
    int run;
    long long event;
    std::string ch;      // PS_A or PS_B
    std::string module;  // A or B
    double time;
    bool matched_to_nai;
};

struct FinalRec {
    int run;
    long long event;
    std::string module;
    double energy_like; // 互換列名 integral に書くが実体はエネルギー相当量
    double time;        // 補正後時刻 [sample]
    std::string particle;   // positron or gamma
    std::string nai_ch;
    std::string ps_ch;      // "-" if gamma
};

static double NaiEnergyFromIntegral(const std::string& ch, double integral);

struct TimeIdx {
    double t;
    int idx;
};

struct TimeIdxLess {
    bool operator()(const TimeIdx& l, const TimeIdx& r) const {
        if (l.t < r.t) return true;
        if (l.t > r.t) return false;
        return l.idx < r.idx;
    }
};

static bool FindNearestAvailable(const std::multiset<TimeIdx, TimeIdxLess>& avail,
                                 double tref,
                                 int& out_idx,
                                 double& out_dt) {
    out_idx = -1;
    out_dt = 1.0e30;
    if (avail.empty()) return false;

    const TimeIdx key{tref, -1};
    auto it = avail.lower_bound(key);
    if (it != avail.end()) {
        out_idx = it->idx;
        out_dt = std::fabs(it->t - tref);
    }
    if (it != avail.begin()) {
        auto it_prev = it;
        --it_prev;
        const double dt = std::fabs(it_prev->t - tref);
        if (dt < out_dt ||
            (std::fabs(dt - out_dt) <= 1.0e-12 &&
             (it_prev->t < (it != avail.end() ? it->t : 1.0e300) ||
              (it != avail.end() && it_prev->t == it->t && it_prev->idx < out_idx)))) {
            out_idx = it_prev->idx;
            out_dt = dt;
        }
    }
    return out_idx >= 0;
}

static int FindNearestByTimeIndex(const std::vector<size_t>& sorted_idx,
                                  const std::vector<double>& sorted_time,
                                  double tref,
                                  double& out_dt) {
    out_dt = 1.0e30;
    if (sorted_idx.empty()) return -1;
    auto it = std::lower_bound(sorted_time.begin(), sorted_time.end(), tref);
    int best_pos = -1;
    if (it != sorted_time.end()) {
        best_pos = static_cast<int>(it - sorted_time.begin());
        out_dt = std::fabs(sorted_time[static_cast<size_t>(best_pos)] - tref);
    }
    if (it != sorted_time.begin()) {
        const int pos_prev = static_cast<int>((it - sorted_time.begin()) - 1);
        const double dt_prev = std::fabs(sorted_time[static_cast<size_t>(pos_prev)] - tref);
        if (dt_prev < out_dt ||
            (std::fabs(dt_prev - out_dt) <= 1.0e-12 &&
             sorted_time[static_cast<size_t>(pos_prev)] < sorted_time[static_cast<size_t>(best_pos)])) {
            best_pos = pos_prev;
            out_dt = dt_prev;
        }
    }
    if (best_pos < 0) return -1;
    return static_cast<int>(sorted_idx[static_cast<size_t>(best_pos)]);
}

struct PsMatchScore {
    int n_pair;
    double sum_dt;
};

static bool BetterPsMatchScore(const PsMatchScore& a, const PsMatchScore& b) {
    if (a.n_pair != b.n_pair) return a.n_pair > b.n_pair;
    if (std::fabs(a.sum_dt - b.sum_dt) > 1.0e-12) return a.sum_dt < b.sum_dt;
    return false;
}

static std::vector<std::pair<size_t, size_t>> FindOptimalPsMatches(
        const std::vector<size_t>& fi,
        const std::vector<size_t>& pi,
        const std::vector<FinalRec>& finals,
        const std::vector<PsPulse>& ps_all,
        double window_sample) {
    const int nf = static_cast<int>(fi.size());
    const int np = static_cast<int>(pi.size());
    std::vector<std::vector<PsMatchScore>> dp(static_cast<size_t>(nf + 1),
                                              std::vector<PsMatchScore>(static_cast<size_t>(np + 1),
                                                                        PsMatchScore{0, 0.0}));
    std::vector<std::vector<char>> act(static_cast<size_t>(nf + 1),
                                       std::vector<char>(static_cast<size_t>(np + 1), 0));

    for (int i = 1; i <= nf; ++i) act[static_cast<size_t>(i)][0] = 1; // skip final
    for (int j = 1; j <= np; ++j) act[0][static_cast<size_t>(j)] = 2; // skip ps

    for (int i = 1; i <= nf; ++i) {
        for (int j = 1; j <= np; ++j) {
            PsMatchScore best = dp[static_cast<size_t>(i - 1)][static_cast<size_t>(j)];
            char best_act = 1;

            const PsMatchScore skip_ps = dp[static_cast<size_t>(i)][static_cast<size_t>(j - 1)];
            if (BetterPsMatchScore(skip_ps, best)) {
                best = skip_ps;
                best_act = 2;
            }

            const double dt = std::fabs(finals[fi[static_cast<size_t>(i - 1)]].time -
                                        ps_all[pi[static_cast<size_t>(j - 1)]].time);
            if (dt <= window_sample) {
                PsMatchScore m = dp[static_cast<size_t>(i - 1)][static_cast<size_t>(j - 1)];
                m.n_pair += 1;
                m.sum_dt += dt;
                if (BetterPsMatchScore(m, best)) {
                    best = m;
                    best_act = 3; // match
                }
            }

            dp[static_cast<size_t>(i)][static_cast<size_t>(j)] = best;
            act[static_cast<size_t>(i)][static_cast<size_t>(j)] = best_act;
        }
    }

    std::vector<std::pair<size_t, size_t>> rev;
    int i = nf;
    int j = np;
    while (i > 0 && j > 0) {
        const char a = act[static_cast<size_t>(i)][static_cast<size_t>(j)];
        if (a == 3) {
            rev.push_back({fi[static_cast<size_t>(i - 1)], pi[static_cast<size_t>(j - 1)]});
            --i;
            --j;
        } else if (a == 1) {
            --i;
        } else if (a == 2) {
            --j;
        } else {
            break;
        }
    }
    std::reverse(rev.begin(), rev.end());
    return rev;
}

struct NaiMatchScore {
    int n_pair;
    double sum_dt;
};

static bool BetterNaiMatchScore(const NaiMatchScore& a, const NaiMatchScore& b) {
    if (a.n_pair != b.n_pair) return a.n_pair > b.n_pair;
    if (std::fabs(a.sum_dt - b.sum_dt) > 1.0e-12) return a.sum_dt < b.sum_dt;
    return false;
}

// 同一event/module内の NaI 2ch を 1対1で全体最適化して対応付ける。
// 目的関数は「ペア数最大化」→「Σ|dt|最小化」。
// pair_window_sample を超える組は採用しない。
static std::vector<std::pair<size_t, size_t>> FindOptimalNaiMatches(
        const std::vector<NaiPulse>& a,
        const std::vector<NaiPulse>& b,
        double pair_window_sample) {
    const int na = static_cast<int>(a.size());
    const int nb = static_cast<int>(b.size());
    std::vector<std::vector<NaiMatchScore>> dp(static_cast<size_t>(na + 1),
                                               std::vector<NaiMatchScore>(static_cast<size_t>(nb + 1),
                                                                          NaiMatchScore{0, 0.0}));
    std::vector<std::vector<char>> act(static_cast<size_t>(na + 1),
                                       std::vector<char>(static_cast<size_t>(nb + 1), 0));

    for (int i = 1; i <= na; ++i) act[static_cast<size_t>(i)][0] = 1; // skip a
    for (int j = 1; j <= nb; ++j) act[0][static_cast<size_t>(j)] = 2; // skip b

    for (int i = 1; i <= na; ++i) {
        for (int j = 1; j <= nb; ++j) {
            NaiMatchScore best = dp[static_cast<size_t>(i - 1)][static_cast<size_t>(j)];
            char best_act = 1;

            const NaiMatchScore skip_b = dp[static_cast<size_t>(i)][static_cast<size_t>(j - 1)];
            if (BetterNaiMatchScore(skip_b, best)) {
                best = skip_b;
                best_act = 2;
            }

            const double dt = std::fabs(a[static_cast<size_t>(i - 1)].time -
                                        b[static_cast<size_t>(j - 1)].time);
            if (dt <= pair_window_sample) {
                NaiMatchScore m = dp[static_cast<size_t>(i - 1)][static_cast<size_t>(j - 1)];
                m.n_pair += 1;
                m.sum_dt += dt;
                if (BetterNaiMatchScore(m, best)) {
                    best = m;
                    best_act = 3; // match
                }
            }

            dp[static_cast<size_t>(i)][static_cast<size_t>(j)] = best;
            act[static_cast<size_t>(i)][static_cast<size_t>(j)] = best_act;
        }
    }

    std::vector<std::pair<size_t, size_t>> rev;
    int i = na;
    int j = nb;
    while (i > 0 && j > 0) {
        const char a_act = act[static_cast<size_t>(i)][static_cast<size_t>(j)];
        if (a_act == 3) {
            rev.push_back({static_cast<size_t>(i - 1), static_cast<size_t>(j - 1)});
            --i;
            --j;
        } else if (a_act == 1) {
            --i;
        } else if (a_act == 2) {
            --j;
        } else {
            break;
        }
    }
    std::reverse(rev.begin(), rev.end());
    return rev;
}

static std::string Trim(const std::string& s) {
    size_t a = 0;
    while (a < s.size() && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    size_t b = s.size();
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}

static std::vector<std::string> SplitCsv(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, ',')) out.push_back(Trim(tok));
    return out;
}

static bool ParseLongLong(const std::string& s, long long& out) {
    char* endp = nullptr;
    errno = 0;
    const long long v = std::strtoll(s.c_str(), &endp, 10);
    if (errno != 0 || endp == s.c_str() || *endp != '\0') return false;
    out = v;
    return true;
}

static bool ParseDouble(const std::string& s, double& out) {
    char* endp = nullptr;
    errno = 0;
    const double v = std::strtod(s.c_str(), &endp);
    if (errno != 0 || endp == s.c_str() || *endp != '\0') return false;
    out = v;
    return true;
}

static std::string ModuleFromCh(const std::string& ch) {
    if (ch.find("_A") != std::string::npos) return "A";
    if (ch.find("_B") != std::string::npos) return "B";
    return "?";
}

static double OffsetFromCh(const std::string& ch) {
    if (ch == "NaI_A1") return kOffset_NaI_A1;
    if (ch == "NaI_A2") return kOffset_NaI_A2;
    if (ch == "NaI_B1") return kOffset_NaI_B1;
    if (ch == "NaI_B2") return kOffset_NaI_B2;
    if (ch == "PS_A") return kOffset_PS_A;
    if (ch == "PS_B") return kOffset_PS_B;
    return 0.0;
}

static double NaiEnergyFromIntegral(const std::string& ch, double integral) {
    if (ch == "NaI_A1") return kNaI_A1_Coeff * integral + kNaI_A1_Offset;
    if (ch == "NaI_A2") return kNaI_A2_Coeff * integral + kNaI_A2_Offset;
    if (ch == "NaI_B1") return kNaI_B1_Coeff * integral + kNaI_B1_Offset;
    if (ch == "NaI_B2") return kNaI_B2_Coeff * integral + kNaI_B2_Offset;
    return integral;
}

static void LoadNaiPulses(const std::string& path,
                          const std::string& ch,
                          std::vector<NaiPulse>& out) {
    std::ifstream ifs(path.c_str());
    if (!ifs) throw std::runtime_error("Failed to open NaI pulse file: " + path);

    std::string line;
    while (std::getline(ifs, line)) {
        line = Trim(line);
        if (line.empty()) continue;
        std::vector<std::string> cols = SplitCsv(line);
        if (cols.size() < 5) continue;

        long long run_ll = 0;
        long long event_id = 0;
        double area = 0.0;

        if (!ParseLongLong(cols[0], run_ll)) continue;
        if (!ParseLongLong(cols[1], event_id)) continue;
        if (!ParseDouble(cols[3], area)) continue;
        double t_raw = 0.0;
        if (!ParseDouble(cols[4], t_raw)) continue;
        const double t = t_raw + OffsetFromCh(ch);

        NaiPulse p;
        p.run = static_cast<int>(run_ll);
        p.event = event_id;
        p.ch = ch;
        p.module = ModuleFromCh(ch);
        p.integral = area;
        p.time = t;
        out.push_back(p);
    }
}

static void LoadPsPulses(const std::string& path,
                         const std::string& ch,
                         std::vector<PsPulse>& out) {
    std::ifstream ifs(path.c_str());
    if (!ifs) throw std::runtime_error("Failed to open PS signal file: " + path);

    std::string line;
    while (std::getline(ifs, line)) {
        line = Trim(line);
        if (line.empty()) continue;
        std::vector<std::string> cols = SplitCsv(line);
        if (cols.size() < 4) continue;

        long long run_ll = 0;
        long long event_id = 0;
        if (!ParseLongLong(cols[0], run_ll)) continue;
        if (!ParseLongLong(cols[1], event_id)) continue;

        double t_raw = 0.0;
        if (!ParseDouble(cols[3], t_raw)) continue;
        const double t = t_raw + OffsetFromCh(ch);

        PsPulse p;
        p.run = static_cast<int>(run_ll);
        p.event = event_id;
        p.ch = ch;
        p.module = ModuleFromCh(ch);
        p.time = t;
        p.matched_to_nai = false;
        out.push_back(p);
    }
}

static void MergeNaiModule(const std::vector<NaiPulse>& x,
                           const std::vector<NaiPulse>& y,
                           std::vector<FinalRec>& out) {
    std::map<std::pair<int, long long>, std::vector<NaiPulse>> groups_a;
    std::map<std::pair<int, long long>, std::vector<NaiPulse>> groups_b;
    for (const auto& p : x) groups_a[{p.run, p.event}].push_back(p);
    for (const auto& p : y) groups_b[{p.run, p.event}].push_back(p);

    std::set<std::pair<int, long long>> keys;
    for (const auto& kv : groups_a) keys.insert(kv.first);
    for (const auto& kv : groups_b) keys.insert(kv.first);

    for (const auto& key : keys) {
        std::vector<NaiPulse> a = groups_a[key];
        std::vector<NaiPulse> b = groups_b[key];
        std::sort(a.begin(), a.end(), [](const NaiPulse& l, const NaiPulse& r) { return l.time < r.time; });
        std::sort(b.begin(), b.end(), [](const NaiPulse& l, const NaiPulse& r) { return l.time < r.time; });
        const std::vector<std::pair<size_t, size_t>> pairs =
            FindOptimalNaiMatches(a, b, static_cast<double>(kNaiPairWindowSamples));
        std::vector<int> matched_b_for_a(a.size(), -1);
        std::vector<bool> matched_b(b.size(), false);
        for (const auto& p : pairs) {
            matched_b_for_a[p.first] = static_cast<int>(p.second);
            matched_b[p.second] = true;
        }

        for (size_t i = 0; i < a.size(); ++i) {
            FinalRec r;
            r.run = a[i].run;
            r.event = a[i].event;
            r.module = a[i].module;
            r.particle = "gamma";
            r.ps_ch = "-";

            const int best_j = matched_b_for_a[i];
            if (best_j >= 0) {
                const double ea = NaiEnergyFromIntegral(a[i].ch, a[i].integral);
                const double eb = NaiEnergyFromIntegral(b[static_cast<size_t>(best_j)].ch,
                                                        b[static_cast<size_t>(best_j)].integral);
                r.energy_like = ea + eb;
                r.time = std::min(a[i].time, b[static_cast<size_t>(best_j)].time);
                r.nai_ch = a[i].ch + "+" + b[static_cast<size_t>(best_j)].ch;
            } else {
                r.energy_like = NaiEnergyFromIntegral(a[i].ch, a[i].integral);
                r.time = a[i].time;
                r.nai_ch = a[i].ch;
            }
            out.push_back(r);
        }

        for (size_t j = 0; j < b.size(); ++j) {
            if (matched_b[j]) continue;
            const NaiPulse& bj = b[j];
            FinalRec r;
            r.run = bj.run;
            r.event = bj.event;
            r.module = bj.module;
            r.energy_like = NaiEnergyFromIntegral(bj.ch, bj.integral);
            r.time = bj.time;
            r.particle = "gamma";
            r.nai_ch = bj.ch;
            r.ps_ch = "-";
            out.push_back(r);
        }
    }
}

static void MatchPsAndClassify(std::vector<FinalRec>& finals,
                               std::vector<PsPulse>& ps_all) {
    std::map<std::tuple<int, long long, std::string>, std::vector<size_t>> finals_by_key;
    std::map<std::tuple<int, long long, std::string>, std::vector<size_t>> ps_by_key;

    for (size_t i = 0; i < finals.size(); ++i) {
        finals_by_key[{finals[i].run, finals[i].event, finals[i].module}].push_back(i);
    }
    for (size_t i = 0; i < ps_all.size(); ++i) {
        ps_by_key[{ps_all[i].run, ps_all[i].event, ps_all[i].module}].push_back(i);
    }

    for (auto& kv : finals_by_key) {
        const auto& key = kv.first;
        std::vector<size_t>& fi = kv.second;
        std::vector<size_t> pi = ps_by_key[key];

        std::sort(fi.begin(), fi.end(), [&](size_t l, size_t r) {
            return finals[l].time < finals[r].time;
        });
        std::sort(pi.begin(), pi.end(), [&](size_t l, size_t r) {
            return ps_all[l].time < ps_all[r].time;
        });

        // 初期化
        for (size_t fidx : fi) {
            FinalRec& r = finals[fidx];
            r.particle = "gamma";
            r.ps_ch = "-";
        }
        for (size_t pidx : pi) ps_all[pidx].matched_to_nai = false;

        // 1対1で、ペア数最大化 + Σ|dt|最小化
        const std::vector<std::pair<size_t, size_t>> pairs =
            FindOptimalPsMatches(fi, pi, finals, ps_all, static_cast<double>(kPsCoincidenceWindowSamples));

        for (const auto& m : pairs) {
            FinalRec& r = finals[m.first];
            PsPulse& p = ps_all[m.second];
            r.particle = "positron";
            r.ps_ch = p.ch;
            p.matched_to_nai = true;
        }
    }
}

static void WriteFinalTxt(const std::string& path, const std::vector<FinalRec>& finals) {
    std::filesystem::path p(path);
    std::error_code ec;
    if (!p.parent_path().empty()) {
        std::filesystem::create_directories(p.parent_path(), ec);
        if (ec) throw std::runtime_error("Failed to create output dir: " + ec.message());
    }

    std::ofstream ofs(path.c_str());
    if (!ofs) throw std::runtime_error("Failed to open output: " + path);

    ofs << "run event module integral t_fall_start_ns particle nai_ch ps_ch\n";
    for (size_t i = 0; i < finals.size(); ++i) {
        const FinalRec& r = finals[i];
        ofs << r.run << " " << r.event << " " << r.module << " "
            << r.energy_like << " " << (r.time * kSampleToTimeNs) << " "
            << r.particle << " " << r.nai_ch << " " << r.ps_ch << "\n";
    }
}

static void WriteUnmatchedPsCsv(const std::string& path, const std::vector<PsPulse>& ps_all) {
    std::filesystem::path p(path);
    std::error_code ec;
    if (!p.parent_path().empty()) {
        std::filesystem::create_directories(p.parent_path(), ec);
        if (ec) throw std::runtime_error("Failed to create output dir: " + ec.message());
    }

    std::ofstream ofs(path.c_str());
    if (!ofs) throw std::runtime_error("Failed to open output: " + path);

    ofs << "run,event,module,ch,time_ns,nearest_dt_same_ch_ns\n";

    for (size_t i = 0; i < ps_all.size(); ++i) {
        const PsPulse& p = ps_all[i];
        if (p.matched_to_nai) continue;

        double nearest = -1.0;
        for (size_t j = 0; j < ps_all.size(); ++j) {
            if (i == j) continue;
            const PsPulse& q = ps_all[j];
            if (q.run != p.run || q.event != p.event || q.ch != p.ch) continue;
            const double dt = std::fabs(q.time - p.time);
            if (nearest < 0 || dt < nearest) nearest = dt;
        }

        ofs << p.run << "," << p.event << "," << p.module << "," << p.ch
            << "," << (p.time * kSampleToTimeNs) << ","
            << ((nearest >= 0.0) ? (nearest * kSampleToTimeNs) : -1.0) << "\n";
    }
}

static void Usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " --nai-a1 FILE --nai-a2 FILE --nai-b1 FILE --nai-b2 FILE \\\n"
        << "            --ps-a FILE --ps-b FILE --out-main FILE --out-unmatched FILE [--debug]\n";
}

int main(int argc, char** argv) {
    try {
        std::string nai_a1, nai_a2, nai_b1, nai_b2;
        std::string ps_a, ps_b;
        std::string out_main, out_unmatched;
        bool debug = false;

        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            auto need = [&](const std::string& opt) -> std::string {
                if (i + 1 >= argc) throw std::runtime_error("Missing value for " + opt);
                return std::string(argv[++i]);
            };

            if (a == "--nai-a1") nai_a1 = need(a);
            else if (a == "--nai-a2") nai_a2 = need(a);
            else if (a == "--nai-b1") nai_b1 = need(a);
            else if (a == "--nai-b2") nai_b2 = need(a);
            else if (a == "--ps-a") ps_a = need(a);
            else if (a == "--ps-b") ps_b = need(a);
            else if (a == "--out-main") out_main = need(a);
            else if (a == "--out-unmatched") out_unmatched = need(a);
            else if (a == "--debug") debug = true;
            else if (a == "-h" || a == "--help") {
                Usage(argv[0]);
                return 0;
            } else {
                throw std::runtime_error("Unknown option: " + a);
            }
        }

        if (nai_a1.empty() || nai_a2.empty() || nai_b1.empty() || nai_b2.empty() ||
            ps_a.empty() || ps_b.empty() || out_main.empty() || out_unmatched.empty()) {
            Usage(argv[0]);
            return 1;
        }

        std::vector<NaiPulse> v_a1, v_a2, v_b1, v_b2;
        std::vector<PsPulse> v_ps;

        LoadNaiPulses(nai_a1, "NaI_A1", v_a1);
        LoadNaiPulses(nai_a2, "NaI_A2", v_a2);
        LoadNaiPulses(nai_b1, "NaI_B1", v_b1);
        LoadNaiPulses(nai_b2, "NaI_B2", v_b2);
        LoadPsPulses(ps_a, "PS_A", v_ps);
        LoadPsPulses(ps_b, "PS_B", v_ps);

        std::vector<FinalRec> finals;
        finals.reserve(v_a1.size() + v_a2.size() + v_b1.size() + v_b2.size());
        MergeNaiModule(v_a1, v_a2, finals);
        MergeNaiModule(v_b1, v_b2, finals);

        MatchPsAndClassify(finals, v_ps);

        std::sort(finals.begin(), finals.end(), [](const FinalRec& l, const FinalRec& r) {
            return std::tie(l.run, l.event, l.module, l.time) <
                   std::tie(r.run, r.event, r.module, r.time);
        });

        WriteFinalTxt(out_main, finals);
        WriteUnmatchedPsCsv(out_unmatched, v_ps);

        if (debug) {
            size_t n_pos = 0;
            size_t n_unmatched_ps = 0;
            for (size_t i = 0; i < finals.size(); ++i) {
                if (finals[i].particle == "positron") ++n_pos;
            }
            for (size_t i = 0; i < v_ps.size(); ++i) {
                if (!v_ps[i].matched_to_nai) ++n_unmatched_ps;
            }
            std::cout << "final_records=" << finals.size()
                      << " positron=" << n_pos
                      << " gamma=" << (finals.size() - n_pos)
                      << " unmatched_ps=" << n_unmatched_ps << "\n";
        }

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
