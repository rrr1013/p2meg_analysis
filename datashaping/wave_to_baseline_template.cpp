#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// ============================================================
// wave_to_baseline_template.cpp
//
// 目的:
//  - wave_*.txt から event_baseline_noise.csv を作る
//  - NaI の場合は template.json も作る（PS は baseline のみ）
//
// 注意:
//  - 既定パラメータは fit_singlepulse_waveforms.cc と同一
//  - 図の作成は行わない
// ============================================================

// ------------------------------
// 既定パラメータ（fit_singlepulse_waveforms.cc と同一）
// ------------------------------
static const int kBaselineWinSamples = 100;
static const double kSnrThreshold = 8.0;
static const int kAlignWinPreSamples = 15;
static const int kAlignWinPostSamples = 15;
static const int kMaxShift = 15;
static const double kDerivThreshold = 2.5;
static const double kRecoveryDip = 2.5;
static const double kSaturationMargin = 2.0;
static const bool kUseSaturationCheck = false;
static const double kFitFrac = 0.20;
static const double kSnrBottomFrac = 0.15;
static const int kMinPrePulseSeparation = 200;
static const double kPrePulseThrSigma = 6.0;
static const int kBaselineRecoverLen = 15;
static const double kBaselineRecoverSigma = 2.5;
static const int kPreIgnoreBeforeMin = 50;
static const int kMinPreSamples = 50;
static const int kMinPostSamples = 200;
static const double kCoarseStartSigma = 2.5;
static const int kCoarseLookback = 120;
static const int kCoarseTargetIndex = 120;
static const int kCoarseStartConsec = 3;
static const double kLeadReboundSigma = 1.0;
static const double kLeadDipSigma = 1.0;
static const double kSampleDt = 1.0;
static const char* kNormDescription = "min=-1 (event-wise, then averaged, final renorm)";

struct SinglePulseDecisionRow {
    size_t event_id = 0;
    double baseline = 0.0;
    double noise_rms = 0.0;
    double noise_mad = 0.0;
    double snr = 0.0;
    int fail_code = 0;
    std::string reason;
};

static double Median(std::vector<double> v) {
    std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
    return v[v.size() / 2];
}

static double ComputeRMS(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double mean = 0.0;
    for (double x : v) mean += x;
    mean /= static_cast<double>(v.size());

    double s2 = 0.0;
    for (double x : v) {
        const double d = x - mean;
        s2 += d * d;
    }
    return std::sqrt(s2 / static_cast<double>(v.size()));
}

static double ComputeMadSigma(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    std::vector<double> tmp = v;
    const double med = Median(tmp);
    std::vector<double> dev;
    dev.reserve(v.size());
    for (double x : v) dev.push_back(std::abs(x - med));
    const double mad = Median(dev);
    return 1.4826 * mad;
}

static const char* FailCodeToReason(int fail_code) {
    switch (fail_code) {
        case 0: return "accepted";
        case 1: return "snr_rms_below_threshold";
        case 2: return "saturation_low";
        case 3: return "recovery_has_negative_derivative";
        case 4: return "recovery_has_secondary_dip";
        case 5: return "saturation_high";
        case 6: return "pre_pulse_too_close";
        case 7: return "baseline_not_recovered_after_pre_pulse";
        case 8: return "insufficient_post_samples";
        case 9: return "lead_rebound";
        case 10: return "insufficient_pre_samples";
        case 12: return "lead_rebound_then_redip";
        case 101: return "snr_mad_below_threshold";
        case 102: return "snr_bottom_fraction_rejected";
        default: return "unknown_rejection";
    }
}

static std::string MakeReasonFileTag(const std::string& reason) {
    std::string out;
    out.reserve(reason.size());
    for (char c : reason) {
        if (std::isalnum(static_cast<unsigned char>(c))) out.push_back(c);
        else out.push_back('_');
    }
    return out;
}

static void WriteSinglePulseDecisionCsv(const std::string& path,
                                        const std::vector<SinglePulseDecisionRow>& rows,
                                        int run_no,
                                        const std::string& ch_name) {
    std::filesystem::path p(path);
    std::error_code ec;
    if (!p.parent_path().empty()) {
        std::filesystem::create_directories(p.parent_path(), ec);
        if (ec) throw std::runtime_error("Failed to create output dir: " + ec.message());
    }

    std::ofstream ofs(path.c_str());
    if (!ofs) throw std::runtime_error("Failed to open output file: " + path);

    ofs << std::setprecision(12);
    ofs << "run,event_id,CH,baseline,noise_rms,noise_mad,snr,fail_code,reason\n";
    for (const auto& row : rows) {
        ofs << run_no << "," << row.event_id << "," << ch_name << ","
            << row.baseline << "," << row.noise_rms << "," << row.noise_mad << ","
            << row.snr << "," << row.fail_code << "," << row.reason << "\n";
    }
}

static void WriteSinglePulseDecisionReports(const std::filesystem::path& out_dir,
                                            const std::vector<SinglePulseDecisionRow>& accepted,
                                            const std::vector<SinglePulseDecisionRow>& rejected,
                                            int run_no,
                                            const std::string& ch_name) {
    const std::string accepted_path = (out_dir / "singlepulse_selected_events.csv").string();
    const std::string rejected_all_path = (out_dir / "singlepulse_rejected_all.csv").string();
    WriteSinglePulseDecisionCsv(accepted_path, accepted, run_no, ch_name);
    WriteSinglePulseDecisionCsv(rejected_all_path, rejected, run_no, ch_name);

    std::vector<std::pair<int, std::string>> reasons;
    for (const auto& row : rejected) {
        bool exists = false;
        for (const auto& kv : reasons) {
            if (kv.first == row.fail_code) {
                exists = true;
                break;
            }
        }
        if (!exists) reasons.emplace_back(row.fail_code, row.reason);
    }
    std::sort(reasons.begin(), reasons.end(),
              [](const std::pair<int, std::string>& a, const std::pair<int, std::string>& b) {
                  return a.first < b.first;
              });

    for (const auto& kv : reasons) {
        std::vector<SinglePulseDecisionRow> rows;
        for (const auto& row : rejected) {
            if (row.fail_code == kv.first) rows.push_back(row);
        }
        const std::string name =
            "singlepulse_rejected_reason_" + std::to_string(kv.first) + "_" +
            MakeReasonFileTag(kv.second) + ".csv";
        WriteSinglePulseDecisionCsv((out_dir / name).string(), rows, run_no, ch_name);
    }
}

static bool ParseInputMetaFromPath(const std::string& path, int& run_no, std::string& ch_name) {
    run_no = -1;
    ch_name = "unknown";
    const std::string base = std::filesystem::path(path).filename().string();

    size_t prun = base.rfind("run");
    if (prun != std::string::npos) {
        size_t p = prun + 3;
        size_t q = p;
        while (q < base.size() && std::isdigit(static_cast<unsigned char>(base[q]))) ++q;
        if (q > p) run_no = std::atoi(base.substr(p, q - p).c_str());
    }

    size_t pmod = std::string::npos;
    if (base.rfind("wave_", 0) == 0) pmod = 5;
    else if (base.rfind("Michel_", 0) == 0) pmod = 7;

    if (pmod != std::string::npos) {
        const size_t end = base.find("_run", pmod);
        if (end != std::string::npos && end > pmod) ch_name = base.substr(pmod, end - pmod);
    }

    return (run_no >= 0);
}

static std::vector<double> LoadBins(const std::string& path) {
    std::ifstream ifs(path.c_str());
    if (!ifs) throw std::runtime_error("Failed to open input file: " + path);

    std::vector<double> bins;
    bins.reserve(1 << 20);

    std::string line;
    long long line_no = 0;
    while (std::getline(ifs, line)) {
        ++line_no;
        auto is_space = [](unsigned char c) { return std::isspace(c); };
        line.erase(line.begin(), std::find_if(line.begin(), line.end(),
                    [&](unsigned char c) { return !is_space(c); }));
        line.erase(std::find_if(line.rbegin(), line.rend(),
                    [&](unsigned char c) { return !is_space(c); }).base(), line.end());
        if (line.empty()) continue;

        char* endptr = nullptr;
        errno = 0;
        const double v = std::strtod(line.c_str(), &endptr);
        if (errno != 0 || endptr == line.c_str() || *endptr != '\0') {
            throw std::runtime_error("Non-numeric line in " + path + " at line " +
                                     std::to_string(line_no) + ": " + line);
        }
        bins.push_back(v);
    }
    return bins;
}

static std::vector<std::vector<double>> SplitEvents(const std::vector<double>& all_bins,
                                                     int bins_per_event) {
    if (bins_per_event <= 0) throw std::runtime_error("bins_per_event must be > 0");
    const size_t n = all_bins.size();
    const size_t N = static_cast<size_t>(bins_per_event);
    if (n % N != 0) {
        throw std::runtime_error("Total bins (" + std::to_string(n) +
                                 ") is not divisible by bins_per_event (" +
                                 std::to_string(bins_per_event) + ").");
    }

    const size_t nevt = n / N;
    std::vector<std::vector<double>> events;
    events.reserve(nevt);
    for (size_t e = 0; e < nevt; ++e) {
        auto start = all_bins.begin() + static_cast<long long>(e * N);
        auto end = start + static_cast<long long>(N);
        events.emplace_back(start, end);
    }
    return events;
}

static double EstimateBaseline(const std::vector<double>& p, double* noise_rms, double* noise_mad) {
    const int N = static_cast<int>(p.size());
    const int nwin = std::min(N, std::max(1, kBaselineWinSamples));

    int best_start = 0;
    double best_mad = 1e99;
    for (int start = 0; start + nwin <= N; ++start) {
        std::vector<double> w(p.begin() + start, p.begin() + start + nwin);
        const double mad = ComputeMadSigma(w);
        if (mad < best_mad) {
            best_mad = mad;
            best_start = start;
        }
    }

    std::vector<double> head(p.begin() + best_start, p.begin() + best_start + nwin);
    const double rms = ComputeRMS(head);
    const double mad = ComputeMadSigma(head);
    if (noise_rms) *noise_rms = rms;
    if (noise_mad) *noise_mad = mad;
    return Median(head);
}

static int FindShiftByCorrelation(const std::vector<double>& y,
                                  const std::vector<double>& ref,
                                  int win_start,
                                  int win_end,
                                  int max_shift) {
    double best = -1e99;
    int best_shift = 0;
    for (int s = -max_shift; s <= max_shift; ++s) {
        double c = 0.0;
        for (int i = win_start; i < win_end; ++i) {
            const int j = i + s;
            if (j < 0 || j >= static_cast<int>(ref.size())) continue;
            c += y[i] * ref[j];
        }
        if (c > best) {
            best = c;
            best_shift = s;
        }
    }
    return best_shift;
}

static std::vector<double> ShiftWave(const std::vector<double>& y, int shift) {
    std::vector<double> out(y.size(), 0.0);
    for (size_t i = 0; i < y.size(); ++i) {
        const int j = static_cast<int>(i) + shift;
        if (j >= 0 && j < static_cast<int>(y.size())) out[static_cast<size_t>(j)] = y[i];
    }
    return out;
}

static void BuildAlignWindowFromT0(int t0_index, int n_samples, int& win_start, int& win_end) {
    win_start = t0_index - kAlignWinPreSamples;
    win_end = t0_index + kAlignWinPostSamples;
    if (win_start < 0) win_start = 0;
    if (win_end > n_samples) win_end = n_samples;
    if (win_end <= win_start + 1) {
        win_start = 0;
        win_end = std::min(n_samples, 2);
    }
}

static int FindTemplateT0ByMaxNegSlope(const std::vector<double>& templ) {
    if (templ.size() < 2) return 0;
    double min_dy = 1e99;
    int idx = 0;
    for (int i = 0; i + 1 < static_cast<int>(templ.size()); ++i) {
        const double dy = templ[static_cast<size_t>(i + 1)] - templ[static_cast<size_t>(i)];
        if (dy < min_dy) {
            min_dy = dy;
            idx = i + 1;
        }
    }
    return idx;
}

static void WriteTemplateJson(const std::string& path,
                              const std::vector<double>& templ,
                              double dt,
                              const char* norm_desc,
                              int t0_index,
                              int run_no,
                              const std::string& ch_name) {
    std::filesystem::path p(path);
    std::error_code ec;
    if (!p.parent_path().empty()) {
        std::filesystem::create_directories(p.parent_path(), ec);
        if (ec) throw std::runtime_error("Failed to create output dir: " + ec.message());
    }

    std::ofstream ofs(path.c_str());
    if (!ofs) throw std::runtime_error("Failed to open output file: " + path);

    double sum = 0.0;
    double sum2 = 0.0;
    for (double v : templ) {
        sum += v;
        sum2 += v * v;
    }
    const double areaT = sum * dt;
    const double sumT2 = sum2;

    const double tmin = *std::min_element(templ.begin(), templ.end());
    const char* polarity = (tmin < 0.0) ? "negative" : "positive";

    ofs << std::setprecision(12);
    ofs << "{\n";
    ofs << "  \"dt\": " << dt << ",\n";
    ofs << "  \"run\": " << run_no << ",\n";
    ofs << "  \"CH\": \"" << ch_name << "\",\n";
    ofs << "  \"polarity\": \"" << polarity << "\",\n";
    ofs << "  \"norm\": \"" << norm_desc << "\",\n";
    ofs << "  \"t0_index\": " << t0_index << ",\n";
    ofs << "  \"template\": [";
    for (size_t i = 0; i < templ.size(); ++i) {
        if (i != 0) ofs << ", ";
        ofs << templ[i];
    }
    ofs << "],\n";
    ofs << "  \"areaT\": " << areaT << ",\n";
    ofs << "  \"sumT2\": " << sumT2 << "\n";
    ofs << "}\n";
}

static void WriteBaselineCsv(const std::string& path,
                             const std::vector<double>& baseline,
                             const std::vector<double>& noise_mad,
                             int run_no,
                             const std::string& ch_name) {
    std::filesystem::path p(path);
    std::error_code ec;
    if (!p.parent_path().empty()) {
        std::filesystem::create_directories(p.parent_path(), ec);
        if (ec) throw std::runtime_error("Failed to create output dir: " + ec.message());
    }

    std::ofstream ofs(path.c_str());
    if (!ofs) throw std::runtime_error("Failed to open output file: " + path);

    ofs << std::setprecision(12);
    ofs << "run,event_id,CH,baseline,noise\n";
    const size_t n = std::min(baseline.size(), noise_mad.size());
    for (size_t i = 0; i < n; ++i) {
        ofs << run_no << "," << i << "," << ch_name << ","
            << baseline[i] << "," << noise_mad[i] << "\n";
    }
}

static int FindFallingStartIndex(const std::vector<double>& y, double noise_mad, int idx_min) {
    if (noise_mad <= 0) return -1;
    if (idx_min <= 0) return -1;

    const double thr = -kCoarseStartSigma * noise_mad;
    const int N = static_cast<int>(y.size());
    const int end = std::min(N - 1, idx_min);
    const int start = std::max(0, end - kCoarseLookback);

    for (int i = end - 1; i >= start; --i) {
        bool ok = true;
        for (int j = 1; j <= kCoarseStartConsec; ++j) {
            if (i + j > end) {
                ok = false;
                break;
            }
            if (!(y[static_cast<size_t>(i + j)] < thr)) {
                ok = false;
                break;
            }
        }
        if (ok && y[static_cast<size_t>(i)] >= thr) return i + 1;
    }
    return -1;
}

static bool IsCorrelationWindowValid(int shift, int n_samples, int t0_index) {
    int win_start = 0;
    int win_end = 0;
    BuildAlignWindowFromT0(t0_index, n_samples, win_start, win_end);
    const int win_last = win_end - 1;
    const int src_start = win_start - shift;
    const int src_end = win_last - shift;
    return (src_start >= 0 && src_end < n_samples);
}

static bool IsSinglePulse(const std::vector<double>& y,
                          double noise_rms,
                          double noise_mad,
                          double adc_min,
                          double adc_max,
                          int* fail_code = nullptr) {
    if (y.empty()) return false;

    auto it_min = std::min_element(y.begin(), y.end());
    const double minv = *it_min;
    const int idx_min = static_cast<int>(it_min - y.begin());

    if (idx_min < kMinPreSamples) {
        if (fail_code) *fail_code = 10;
        return false;
    }
    if (static_cast<int>(y.size()) - 1 - idx_min < kMinPostSamples) {
        if (fail_code) *fail_code = 8;
        return false;
    }
    if (std::abs(minv) <= kSnrThreshold * noise_rms) {
        if (fail_code) *fail_code = 1;
        return false;
    }

    const double thr = kFitFrac * minv;
    int fit_start = 0;
    int fit_end = static_cast<int>(y.size()) - 1;
    for (int i = 0; i < static_cast<int>(y.size()); ++i) {
        if (y[static_cast<size_t>(i)] < thr) {
            fit_start = i;
            break;
        }
    }
    for (int i = static_cast<int>(y.size()) - 1; i >= 0; --i) {
        if (y[static_cast<size_t>(i)] < thr) {
            fit_end = i;
            break;
        }
    }

    if (noise_mad > 0.0) {
        int last_pre_idx = -1;
        const double pre_thr = -kPrePulseThrSigma * noise_mad;
        const int pre_end = std::max(0, idx_min - kPreIgnoreBeforeMin);
        for (int i = 0; i < pre_end; ++i) {
            if (y[static_cast<size_t>(i)] < pre_thr) last_pre_idx = i;
        }

        if (last_pre_idx >= 0) {
            const int dt = idx_min - last_pre_idx;
            if (dt < kMinPrePulseSeparation) {
                if (fail_code) *fail_code = 6;
                return false;
            }

            int consec = 0;
            const double rec_thr = kBaselineRecoverSigma * noise_mad;
            for (int i = last_pre_idx + 1; i < pre_end; ++i) {
                if (std::abs(y[static_cast<size_t>(i)]) < rec_thr) {
                    ++consec;
                    if (consec >= kBaselineRecoverLen) break;
                } else {
                    consec = 0;
                }
            }
            if (consec < kBaselineRecoverLen) {
                if (fail_code) *fail_code = 7;
                return false;
            }
        }
    }

    if (noise_mad > 0.0) {
        const int t_le = FindFallingStartIndex(y, noise_mad, idx_min);
        if (t_le >= 0 && idx_min > t_le + 2) {
            double running_min = y[static_cast<size_t>(t_le)];
            double max_rebound = 0.0;
            bool rebound_seen = false;
            double rebound_peak = y[static_cast<size_t>(t_le)];
            for (int i = t_le + 1; i <= idx_min; ++i) {
                const double yi = y[static_cast<size_t>(i)];
                if (yi < running_min) running_min = yi;
                const double rebound = yi - running_min;
                if (rebound > max_rebound) max_rebound = rebound;
                if (rebound > kLeadReboundSigma * noise_mad) {
                    rebound_seen = true;
                    if (yi > rebound_peak) rebound_peak = yi;
                }
                if (rebound_seen && (rebound_peak - yi) > kLeadDipSigma * noise_mad) {
                    if (fail_code) *fail_code = 12;
                    return false;
                }
            }
            if (max_rebound > kLeadReboundSigma * noise_mad) {
                if (fail_code) *fail_code = 9;
                return false;
            }
        }
    }

    if (kUseSaturationCheck && minv <= adc_min + kSaturationMargin) {
        if (fail_code) *fail_code = 2;
        return false;
    }

    for (int i = idx_min + 1; i + 1 <= fit_end; ++i) {
        const double dy = y[static_cast<size_t>(i + 1)] - y[static_cast<size_t>(i)];
        if (dy < -kDerivThreshold * noise_mad) {
            if (fail_code) *fail_code = 3;
            return false;
        }
    }

    double max_after = y[static_cast<size_t>(idx_min)];
    for (int i = idx_min + 1; i <= fit_end; ++i) {
        const double yi = y[static_cast<size_t>(i)];
        if (yi > max_after) max_after = yi;
        if (yi < max_after - kRecoveryDip * noise_mad) {
            if (fail_code) *fail_code = 4;
            return false;
        }
    }

    const double maxv = *std::max_element(y.begin(), y.end());
    if (kUseSaturationCheck && maxv >= adc_max - kSaturationMargin) {
        if (fail_code) *fail_code = 5;
        return false;
    }

    return true;
}

static void BuildTemplateNoPlot(const std::vector<std::vector<double>>& pulses,
                                int run_no,
                                const std::string& ch_name,
                                const std::string& out_baseline,
                                const std::string& out_template,
                                bool baseline_only,
                                bool debug) {
    if (pulses.empty()) throw std::runtime_error("No events in input waveform.");
    const int N = static_cast<int>(pulses[0].size());

    double adc_min = 1e99;
    double adc_max = -1e99;
    for (const auto& p : pulses) {
        const auto mm = std::minmax_element(p.begin(), p.end());
        if (*mm.first < adc_min) adc_min = *mm.first;
        if (*mm.second > adc_max) adc_max = *mm.second;
    }

    std::vector<std::vector<double>> y_all;
    std::vector<double> noise_rms_all;
    std::vector<double> noise_mad_all;
    std::vector<double> snr_all;
    std::vector<double> baseline_all;
    y_all.reserve(pulses.size());
    noise_rms_all.reserve(pulses.size());
    noise_mad_all.reserve(pulses.size());
    snr_all.reserve(pulses.size());
    baseline_all.reserve(pulses.size());

    for (const auto& p : pulses) {
        double noise_rms = 0.0;
        double noise_mad = 0.0;
        const double base = EstimateBaseline(p, &noise_rms, &noise_mad);

        std::vector<double> y = p;
        for (double& v : y) v -= base;

        const double minv = *std::min_element(y.begin(), y.end());
        const double snr = (noise_mad > 0.0) ? std::abs(minv) / noise_mad : 0.0;

        y_all.push_back(std::move(y));
        noise_rms_all.push_back(noise_rms);
        noise_mad_all.push_back(noise_mad);
        snr_all.push_back(snr);
        baseline_all.push_back(base);
    }

    WriteBaselineCsv(out_baseline, baseline_all, noise_mad_all, run_no, ch_name);

    if (baseline_only) return;

    std::vector<double> snr_pass;
    snr_pass.reserve(snr_all.size());
    for (double s : snr_all) {
        if (s > kSnrThreshold) snr_pass.push_back(s);
    }

    double snr_bottom_cut = -1.0;
    if (!snr_pass.empty()) {
        std::sort(snr_pass.begin(), snr_pass.end());
        const size_t idx = static_cast<size_t>(
            std::floor(kSnrBottomFrac * static_cast<double>(snr_pass.size() - 1)));
        snr_bottom_cut = snr_pass[idx];
    }

    std::vector<std::vector<double>> cleaned;
    std::vector<double> cleaned_mad;
    std::vector<SinglePulseDecisionRow> accepted_rows;
    std::vector<SinglePulseDecisionRow> rejected_rows;
    cleaned.reserve(y_all.size());
    cleaned_mad.reserve(y_all.size());
    accepted_rows.reserve(y_all.size());
    rejected_rows.reserve(y_all.size());

    for (size_t i = 0; i < y_all.size(); ++i) {
        const std::vector<double>& y = y_all[i];
        const double noise_rms = noise_rms_all[i];
        const double noise_mad = noise_mad_all[i];
        const double snr = snr_all[i];

        SinglePulseDecisionRow row;
        row.event_id = i;
        row.baseline = baseline_all[i];
        row.noise_rms = noise_rms;
        row.noise_mad = noise_mad;
        row.snr = snr;

        if (snr <= kSnrThreshold) {
            row.fail_code = 101;
            row.reason = FailCodeToReason(row.fail_code);
            rejected_rows.push_back(row);
            continue;
        }
        if (snr_bottom_cut >= 0 && snr < snr_bottom_cut) {
            row.fail_code = 102;
            row.reason = FailCodeToReason(row.fail_code);
            rejected_rows.push_back(row);
            continue;
        }

        int fail_code = 0;
        if (IsSinglePulse(y, noise_rms, noise_mad, adc_min, adc_max, &fail_code)) {
            cleaned.push_back(y);
            cleaned_mad.push_back(noise_mad);
            row.fail_code = 0;
            row.reason = FailCodeToReason(0);
            accepted_rows.push_back(row);
        } else {
            row.fail_code = fail_code;
            row.reason = FailCodeToReason(fail_code);
            rejected_rows.push_back(row);
        }
    }

    WriteSinglePulseDecisionReports(std::filesystem::path(out_baseline).parent_path(),
                                    accepted_rows, rejected_rows, run_no, ch_name);

    if (cleaned.empty()) {
        throw std::runtime_error("No single-pulse events passed selection for template.");
    }

    std::vector<std::vector<double>> coarse_aligned;
    coarse_aligned.reserve(cleaned.size());
    for (size_t ie = 0; ie < cleaned.size(); ++ie) {
        const std::vector<double>& y = cleaned[ie];
        const double noise_mad = cleaned_mad[ie];
        std::vector<double> ycoarse = y;

        const int idx_min = static_cast<int>(std::min_element(y.begin(), y.end()) - y.begin());
        const int fall_idx = FindFallingStartIndex(y, noise_mad, idx_min);
        if (fall_idx >= 0) {
            const int s_coarse = kCoarseTargetIndex - fall_idx;
            if (IsCorrelationWindowValid(s_coarse, N, kCoarseTargetIndex)) {
                ycoarse = ShiftWave(y, s_coarse);
            }
        }
        coarse_aligned.push_back(ycoarse);
    }

    std::vector<double> ref(static_cast<size_t>(N), 0.0);
    for (const auto& y : coarse_aligned) {
        for (int i = 0; i < N; ++i) ref[static_cast<size_t>(i)] += y[static_cast<size_t>(i)];
    }
    for (double& v : ref) v /= static_cast<double>(coarse_aligned.size());

    int align_win_start = 0;
    int align_win_end = 0;
    BuildAlignWindowFromT0(kCoarseTargetIndex, N, align_win_start, align_win_end);

    std::vector<std::vector<double>> aligned_norm;
    aligned_norm.reserve(coarse_aligned.size());
    for (const auto& ycoarse : coarse_aligned) {
        const int s = FindShiftByCorrelation(ycoarse, ref, align_win_start, align_win_end, kMaxShift);
        std::vector<double> yshift = ShiftWave(ycoarse, -s);

        double ymin = *std::min_element(yshift.begin(), yshift.end());
        if (ymin < 0.0) {
            for (double& v : yshift) v /= -ymin;
        }
        aligned_norm.push_back(std::move(yshift));
    }

    std::vector<double> templ(static_cast<size_t>(N), 0.0);
    for (const auto& y : aligned_norm) {
        for (int i = 0; i < N; ++i) templ[static_cast<size_t>(i)] += y[static_cast<size_t>(i)];
    }
    for (double& v : templ) v /= static_cast<double>(aligned_norm.size());

    const double tmin = *std::min_element(templ.begin(), templ.end());
    if (tmin < 0.0) {
        for (double& v : templ) v /= -tmin;
    }

    const int t0_slope = FindTemplateT0ByMaxNegSlope(templ);
    WriteTemplateJson(out_template, templ, kSampleDt, kNormDescription, t0_slope, run_no, ch_name);

    if (debug) {
        std::cout << "template: selected=" << cleaned.size()
                  << " total=" << pulses.size()
                  << " t0_index=" << t0_slope
                  << " out=" << out_template << "\n";
    }
}

static void Usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " <input_wave.txt> --out-baseline <event_baseline_noise.csv> [options]\n"
        << "Options:\n"
        << "  --out-template <template.json>   NaI用テンプレート出力先\n"
        << "  --bins-per-event N               1イベントサンプル数 (default: 1000)\n"
        << "  --baseline-only                  baseline/noise のみ出力\n"
        << "  --debug                          進捗表示\n";
}

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            Usage(argv[0]);
            return 1;
        }

        std::string input_path;
        std::string out_baseline;
        std::string out_template;
        int bins_per_event = 1000;
        bool baseline_only = false;
        bool debug = false;

        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            if (a == "--out-baseline" && i + 1 < argc) {
                out_baseline = argv[++i];
            } else if (a == "--out-template" && i + 1 < argc) {
                out_template = argv[++i];
            } else if (a == "--bins-per-event" && i + 1 < argc) {
                bins_per_event = std::atoi(argv[++i]);
            } else if (a == "--baseline-only") {
                baseline_only = true;
            } else if (a == "--debug") {
                debug = true;
            } else if (a == "-h" || a == "--help") {
                Usage(argv[0]);
                return 0;
            } else if (!a.empty() && a[0] == '-') {
                std::cerr << "Unknown option: " << a << "\n";
                Usage(argv[0]);
                return 1;
            } else if (input_path.empty()) {
                input_path = a;
            } else {
                std::cerr << "Too many positional arguments.\n";
                Usage(argv[0]);
                return 1;
            }
        }

        if (input_path.empty() || out_baseline.empty()) {
            Usage(argv[0]);
            return 1;
        }

        if (!baseline_only && out_template.empty()) {
            std::cerr << "error: --out-template is required unless --baseline-only is set.\n";
            return 1;
        }

        int run_no = -1;
        std::string ch_name = "unknown";
        ParseInputMetaFromPath(input_path, run_no, ch_name);

        const std::vector<double> all_bins = LoadBins(input_path);
        const std::vector<std::vector<double>> events = SplitEvents(all_bins, bins_per_event);

        BuildTemplateNoPlot(events, run_no, ch_name, out_baseline, out_template, baseline_only, debug);

        std::cout << "Wrote baseline: " << out_baseline << "\n";
        if (!baseline_only) std::cout << "Wrote template: " << out_template << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
