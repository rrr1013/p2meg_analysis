// ============================================================
// signal_mock.dat と rmd_mock.dat を行単位でミックスして新しい模擬データを作る
//
// - 1行=1イベントとして扱う（列フォーマットに依存しない）
// - 先頭が '#' の行や空行は「ヘッダ扱い」で入力から読み飛ばし、出力には生成情報を付与
// - N_sig, N_rmd を直接指定
// - 出力ファイル名に Ns/Nr/total が反映される
// - 保存先は data/mockdata/（デフォルト）
//
// 例:
//   mkdir -p build
//   g++ -O2 -std=c++17 -Wall -Wextra -pedantic -o build/mix_mockdata scripts/mix_mockdata.cc
//   ./build/mix_mockdata --sig data/signal_mock.dat --rmd data/rmd_mock.dat --nsig 3000 --nrmd 7000
// ============================================================

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

static bool IsBlankLine(const std::string& s) {
    for (unsigned char c : s) {
        if (!std::isspace(c)) return false;
    }
    return true;
}

static bool IsHeaderLine(const std::string& s) {
    if (IsBlankLine(s)) return true;
    std::size_t i = 0;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    return (i < s.size() && s[i] == '#');
}

static bool ReadDataLines(const std::string& path,
                          std::vector<std::string>& header_out,
                          std::vector<std::string>& data_out) {
    std::ifstream fin(path);
    if (!fin) {
        std::cerr << "ERROR: failed to open: " << path << "\n";
        return false;
    }
    std::string line;
    while (std::getline(fin, line)) {
        if (IsHeaderLine(line)) header_out.push_back(line);
        else data_out.push_back(line);
    }
    return true;
}

static bool ParseInt64(const std::string& s, long long& out) {
    errno = 0;
    char* endp = nullptr;
    long long v = std::strtoll(s.c_str(), &endp, 10);
    if (errno != 0) return false;
    if (endp == s.c_str() || *endp != '\0') return false;
    out = v;
    return true;
}

static void PrintUsage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " --sig <signal_mock.dat> --rmd <rmd_mock.dat> --nsig <Nsig> --nrmd <Nrmd> [options]\n"
        << "\n"
        << "Required:\n"
        << "  --sig   PATH   input signal mock data\n"
        << "  --rmd   PATH   input rmd mock data\n"
        << "  --nsig  INT    number of signal events in output\n"
        << "  --nrmd  INT    number of rmd events in output\n"
        << "\n"
        << "Options:\n"
        << "  --seed  INT    RNG seed (default: 12345)\n"
        << "  --outdir PATH  output directory (default: data/mockdata)\n"
        << "  --prefix STR   output filename prefix (default: mixed)\n"
        << "  --strict       if input has fewer than requested, exit with error (no replacement)\n"
        << "\n";
}

int main(int argc, char** argv) {
    std::string sig_path;
    std::string rmd_path;
    std::string outdir = "data/mockdata";
    std::string prefix = "mixed";
    long long nsig_ll = -1;
    long long nrmd_ll = -1;
    long long seed_ll = 12345;
    bool strict = false;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto need_value = [&](const char* opt) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "ERROR: missing value for " << opt << "\n";
                std::exit(2);
            }
            return std::string(argv[++i]);
        };

        if (a == "--sig") {
            sig_path = need_value("--sig");
        } else if (a == "--rmd") {
            rmd_path = need_value("--rmd");
        } else if (a == "--nsig") {
            std::string v = need_value("--nsig");
            if (!ParseInt64(v, nsig_ll) || nsig_ll < 0) {
                std::cerr << "ERROR: invalid --nsig: " << v << "\n";
                return 2;
            }
        } else if (a == "--nrmd") {
            std::string v = need_value("--nrmd");
            if (!ParseInt64(v, nrmd_ll) || nrmd_ll < 0) {
                std::cerr << "ERROR: invalid --nrmd: " << v << "\n";
                return 2;
            }
        } else if (a == "--seed") {
            std::string v = need_value("--seed");
            if (!ParseInt64(v, seed_ll)) {
                std::cerr << "ERROR: invalid --seed: " << v << "\n";
                return 2;
            }
        } else if (a == "--outdir") {
            outdir = need_value("--outdir");
        } else if (a == "--prefix") {
            prefix = need_value("--prefix");
        } else if (a == "--strict") {
            strict = true;
        } else if (a == "-h" || a == "--help") {
            PrintUsage(argv[0]);
            return 0;
        } else {
            std::cerr << "ERROR: unknown option: " << a << "\n";
            PrintUsage(argv[0]);
            return 2;
        }
    }

    if (sig_path.empty() || rmd_path.empty() || nsig_ll < 0 || nrmd_ll < 0) {
        PrintUsage(argv[0]);
        return 2;
    }

    const std::size_t n_sig = static_cast<std::size_t>(nsig_ll);
    const std::size_t n_rmd = static_cast<std::size_t>(nrmd_ll);
    const std::size_t n_total = n_sig + n_rmd;

    // 入力読み込み
    std::vector<std::string> sig_header, sig_data;
    std::vector<std::string> rmd_header, rmd_data;

    if (!ReadDataLines(sig_path, sig_header, sig_data)) return 1;
    if (!ReadDataLines(rmd_path, rmd_header, rmd_data)) return 1;

    if (sig_data.empty() && n_sig > 0) {
        std::cerr << "ERROR: no data lines in signal file: " << sig_path << "\n";
        return 1;
    }
    if (rmd_data.empty() && n_rmd > 0) {
        std::cerr << "ERROR: no data lines in rmd file: " << rmd_path << "\n";
        return 1;
    }

    std::mt19937_64 rng(static_cast<std::uint64_t>(seed_ll));

    auto sample_lines = [&](const std::vector<std::string>& src, std::size_t k, const char* label) {
        std::vector<std::string> out;
        out.reserve(k);
        if (k == 0) return out;

        if (k <= src.size()) {
            // without replacement
            std::vector<std::size_t> idx(src.size());
            for (std::size_t i = 0; i < src.size(); ++i) idx[i] = i;
            std::shuffle(idx.begin(), idx.end(), rng);
            for (std::size_t i = 0; i < k; ++i) out.push_back(src[idx[i]]);
        } else {
            if (strict) {
                std::cerr << "ERROR: not enough events in " << label << " (requested " << k
                          << ", available " << src.size() << ") with --strict\n";
                std::exit(1);
            }
            // with replacement
            std::uniform_int_distribution<std::size_t> dist(0, src.size() - 1);
            std::cerr << "WARNING: requested " << k << " events but only " << src.size()
                      << " available in " << label << "; sampling WITH replacement.\n";
            for (std::size_t i = 0; i < k; ++i) out.push_back(src[dist(rng)]);
        }
        return out;
    };

    std::vector<std::string> mixed;
    mixed.reserve(n_total);

    auto sig_part = sample_lines(sig_data, n_sig, "signal");
    auto rmd_part = sample_lines(rmd_data, n_rmd, "rmd");

    mixed.insert(mixed.end(), sig_part.begin(), sig_part.end());
    mixed.insert(mixed.end(), rmd_part.begin(), rmd_part.end());

    // 全体シャッフル
    std::shuffle(mixed.begin(), mixed.end(), rng);

    // 出力ディレクトリ作成
    std::filesystem::create_directories(outdir);

    std::ostringstream outname;
    outname << prefix
            << "_Ns" << n_sig
            << "_Nr" << n_rmd
            << "_N" << n_total
            << ".dat";

    const std::filesystem::path outpath = std::filesystem::path(outdir) / outname.str();

    std::ofstream fout(outpath);
    if (!fout) {
        std::cerr << "ERROR: failed to open output: " << outpath.string() << "\n";
        return 1;
    }

    // 生成情報ヘッダ
    fout << "# mixed mockdata generated by mix_mockdata\n";
    fout << "# sig_path: " << sig_path << "\n";
    fout << "# rmd_path: " << rmd_path << "\n";
    fout << "# N_sig: " << n_sig << "\n";
    fout << "# N_rmd: " << n_rmd << "\n";
    fout << "# N_total: " << n_total << "\n";
    fout << "# seed: " << seed_ll << "\n";
    fout << "# strict: " << (strict ? 1 : 0) << "\n";

    for (const auto& line : mixed) {
        fout << line << "\n";
    }

    std::cout << outpath.string() << "\n";
    return 0;
}
