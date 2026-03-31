#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// ============================================================
// wave_to_pileup_csv.cpp
//
// 目的:
//   1サンプル/行の波形テキストを、1行1イベントのCSVに変換する。
//
// 入力:
//   wave_*.txt (1サンプル/行)
//
// 出力:
//   pileup_waveforms.csv
//   形式: event_id,y0,y1,...,y(N-1)
// ============================================================

static std::string Trim(const std::string& s) {
    size_t a = 0;
    while (a < s.size() && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    size_t b = s.size();
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
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
        line = Trim(line);
        if (line.empty()) continue;

        char* endptr = nullptr;
        errno = 0;
        double v = std::strtod(line.c_str(), &endptr);
        if (errno != 0 || endptr == line.c_str() || *endptr != '\0') {
            throw std::runtime_error(
                "Non-numeric line in " + path + " at line " + std::to_string(line_no) + ": " + line
            );
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
        throw std::runtime_error(
            "Total bins (" + std::to_string(n) + ") is not divisible by bins_per_event (" +
            std::to_string(bins_per_event) + ")."
        );
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

static void WriteWaveformsCsv(const std::string& path,
                              const std::vector<std::vector<double>>& events) {
    std::filesystem::path p(path);
    std::error_code ec;
    if (!p.parent_path().empty()) {
        std::filesystem::create_directories(p.parent_path(), ec);
        if (ec) {
            throw std::runtime_error("Failed to create output dir: " + p.parent_path().string() +
                                     " (" + ec.message() + ")");
        }
    }

    std::ofstream ofs(path.c_str());
    if (!ofs) throw std::runtime_error("Failed to open output file: " + path);

    for (size_t i = 0; i < events.size(); ++i) {
        ofs << i;
        const std::vector<double>& y = events[i];
        for (size_t j = 0; j < y.size(); ++j) {
            ofs << "," << y[j];
        }
        ofs << "\n";
    }
}

static void Usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " <input_wave.txt> <output_pileup.csv> [--bins-per-event N]\n";
}

int main(int argc, char** argv) {
    try {
        if (argc < 3) {
            Usage(argv[0]);
            return 1;
        }

        std::string input_path;
        std::string output_path;
        int bins_per_event = 1000;

        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            if (a == "--bins-per-event" && i + 1 < argc) {
                bins_per_event = std::atoi(argv[++i]);
            } else if (a == "-h" || a == "--help") {
                Usage(argv[0]);
                return 0;
            } else if (!a.empty() && a[0] == '-') {
                std::cerr << "Unknown option: " << a << "\n";
                Usage(argv[0]);
                return 1;
            } else if (input_path.empty()) {
                input_path = a;
            } else if (output_path.empty()) {
                output_path = a;
            } else {
                std::cerr << "Too many positional arguments.\n";
                Usage(argv[0]);
                return 1;
            }
        }

        if (input_path.empty() || output_path.empty()) {
            Usage(argv[0]);
            return 1;
        }

        const std::vector<double> all_bins = LoadBins(input_path);
        const std::vector<std::vector<double>> events = SplitEvents(all_bins, bins_per_event);
        WriteWaveformsCsv(output_path, events);

        std::cout << "Wrote " << events.size() << " events to " << output_path << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
