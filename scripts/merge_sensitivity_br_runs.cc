// scripts/merge_sensitivity_br_runs.cc
//
// 複数の run_sensitivity_br_toymc 出力を結合し、
// toy ごとの BR90 リストをまとめた summary を作る。

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct RunMeta {
    std::string path;
    int n_calib_toys;
    int local_exact_inner_toys;
    int br_scan_points;
    double br_min;
    double br_max;
    double br_step;
    double observed_local_exact_br90;
    std::vector<double> toy_br90;
    bool ok;
};

static bool StartsWith(const std::string& s, const char* prefix)
{
    const std::string p(prefix);
    return s.rfind(p, 0) == 0;
}

static double ParseValueAfterEquals(const std::string& line, bool& ok)
{
    ok = false;
    const std::size_t pos = line.find('=');
    if (pos == std::string::npos) return 0.0;
    std::istringstream iss(line.substr(pos + 1U));
    double v = 0.0;
    if (!(iss >> v)) return 0.0;
    ok = std::isfinite(v);
    return v;
}

static int ParseIntAfterEquals(const std::string& line, bool& ok)
{
    ok = false;
    const std::size_t pos = line.find('=');
    if (pos == std::string::npos) return 0;
    std::istringstream iss(line.substr(pos + 1U));
    int v = 0;
    if (!(iss >> v)) return 0;
    ok = true;
    return v;
}

static bool LoadRun(const char* path, RunMeta& out)
{
    out = RunMeta{};
    out.path = path ? path : "";
    out.n_calib_toys = -1;
    out.local_exact_inner_toys = -1;
    out.br_scan_points = -1;
    out.br_min = std::numeric_limits<double>::quiet_NaN();
    out.br_max = std::numeric_limits<double>::quiet_NaN();
    out.br_step = std::numeric_limits<double>::quiet_NaN();
    out.observed_local_exact_br90 = std::numeric_limits<double>::quiet_NaN();
    out.ok = false;

    std::ifstream fin(path);
    if (!fin) return false;

    std::string line;
    bool in_obs_block = false;
    bool in_toy_block = false;
    double prev_obs_br = std::numeric_limits<double>::quiet_NaN();
    while (std::getline(fin, line)) {
        bool ok = false;
        if (StartsWith(line, "n_calib_toys")) {
            out.n_calib_toys = ParseIntAfterEquals(line, ok);
            continue;
        }
        if (StartsWith(line, "local_exact_inner_toys")) {
            out.local_exact_inner_toys = ParseIntAfterEquals(line, ok);
            continue;
        }
        if (StartsWith(line, "br_scan_points")) {
            out.br_scan_points = ParseIntAfterEquals(line, ok);
            continue;
        }
        if (StartsWith(line, "observed_local_exact_BR90")) {
            out.observed_local_exact_br90 = ParseValueAfterEquals(line, ok);
            continue;
        }

        if (line.find("=================== Observed") != std::string::npos) {
            in_obs_block = true;
            continue;
        }
        if (line.find("=============== Sensitivity Toys") != std::string::npos) {
            in_obs_block = false;
            in_toy_block = true;
            continue;
        }
        if (in_toy_block && line.find("===============================================") != std::string::npos) {
            in_toy_block = false;
            continue;
        }

        if (in_obs_block) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);
            double br = 0.0;
            double q_obs = 0.0;
            double p_value = 0.0;
            int accepted = 0;
            if (!(iss >> br >> q_obs >> p_value >> accepted)) continue;
            if (!std::isfinite(br)) continue;
            if (!std::isfinite(prev_obs_br)) {
                out.br_min = br;
            }
            if (std::isfinite(prev_obs_br) && !std::isfinite(out.br_step)) {
                const double step = br - prev_obs_br;
                if (step > 0.0 && std::isfinite(step)) out.br_step = step;
            }
            prev_obs_br = br;
            out.br_max = br;
            continue;
        }

        if (in_toy_block) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);
            int idx = -1;
            double br90 = 0.0;
            if (!(iss >> idx >> br90)) continue;
            if (!std::isfinite(br90)) continue;
            out.toy_br90.push_back(br90);
        }
    }

    out.ok = (out.n_calib_toys > 0 &&
              out.local_exact_inner_toys > 0 &&
              out.br_scan_points > 0 &&
              !out.toy_br90.empty());
    return out.ok;
}

static double QuantileFromSorted(const std::vector<double>& sorted, double prob)
{
    if (sorted.empty()) return 0.0;
    if (prob <= 0.0) return sorted.front();
    if (prob >= 1.0) return sorted.back();
    const double x = prob * static_cast<double>(sorted.size() - 1U);
    const std::size_t i0 = static_cast<std::size_t>(std::floor(x));
    const std::size_t i1 = static_cast<std::size_t>(std::ceil(x));
    const double t = x - static_cast<double>(i0);
    return (1.0 - t) * sorted[i0] + t * sorted[i1];
}

int main(int argc, char** argv)
{
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <output.txt> <input1.txt> <input2.txt> [input3.txt ...]\n";
        return 1;
    }

    const char* output = argv[1];
    std::vector<RunMeta> runs;
    runs.reserve(static_cast<std::size_t>(argc - 2));
    for (int i = 2; i < argc; ++i) {
        RunMeta run;
        if (!LoadRun(argv[i], run)) {
            std::cerr << "[merge_sensitivity_br_runs] failed to load " << argv[i] << "\n";
            return 2;
        }
        runs.push_back(run);
    }

    const RunMeta& ref = runs.front();
    for (std::size_t i = 1; i < runs.size(); ++i) {
        const RunMeta& run = runs[i];
        if (run.n_calib_toys != ref.n_calib_toys ||
            run.local_exact_inner_toys != ref.local_exact_inner_toys ||
            run.br_scan_points != ref.br_scan_points) {
            std::cerr << "[merge_sensitivity_br_runs] incompatible config in "
                      << run.path << "\n";
            return 3;
        }
        if (std::fabs(run.br_min - ref.br_min) > 1e-12 ||
            std::fabs(run.br_max - ref.br_max) > 1e-12 ||
            std::fabs(run.br_step - ref.br_step) > 1e-12) {
            std::cerr << "[merge_sensitivity_br_runs] incompatible BR grid in "
                      << run.path << "\n";
            return 4;
        }
    }

    std::vector<double> merged;
    for (const auto& run : runs) {
        merged.insert(merged.end(), run.toy_br90.begin(), run.toy_br90.end());
    }
    std::sort(merged.begin(), merged.end());

    std::ofstream fout(output);
    if (!fout) {
        std::cerr << "[merge_sensitivity_br_runs] cannot open output " << output << "\n";
        return 5;
    }

    fout << std::setprecision(10);
    fout << "==================== Merged Sensitivity ====================\n";
    fout << "merged_input_files          = " << (argc - 2) << "\n";
    fout << "n_calib_toys                = " << ref.n_calib_toys << "\n";
    fout << "local_exact_inner_toys      = " << ref.local_exact_inner_toys << "\n";
    fout << "br_scan_points              = " << ref.br_scan_points << "\n";
    fout << "br_min                      = " << ref.br_min << "\n";
    fout << "br_max                      = " << ref.br_max << "\n";
    fout << "br_step                     = " << ref.br_step << "\n";
    fout << "observed_local_exact_BR90   = " << ref.observed_local_exact_br90 << "\n";
    fout << "sensitivity_valid_toys      = " << merged.size() << "\n";
    fout << "sensitivity_BR50            = " << QuantileFromSorted(merged, 0.50) << "\n";
    fout << "sensitivity_BR16            = " << QuantileFromSorted(merged, 0.16) << "\n";
    fout << "sensitivity_BR84            = " << QuantileFromSorted(merged, 0.84) << "\n";
    fout << "sensitivity_BR05            = " << QuantileFromSorted(merged, 0.05) << "\n";
    fout << "sensitivity_BR95            = " << QuantileFromSorted(merged, 0.95) << "\n";
    fout << "===============================================\n";
    fout << "=============== Sensitivity Toys ==============\n";
    fout << "# toy_index  BR90\n";
    for (std::size_t i = 0; i < merged.size(); ++i) {
        fout << i << " " << merged[i] << "\n";
    }
    fout << "===============================================\n";
    fout << "# source files\n";
    for (const auto& run : runs) {
        fout << "# " << run.path << "\n";
    }

    return 0;
}
