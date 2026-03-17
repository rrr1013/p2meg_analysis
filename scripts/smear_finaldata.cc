#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "TRandom3.h"

#include "p2meg/DetectorResolution.h"

// 5列テキスト:
//   Ee Eg t phi_detector_e phi_detector_g
// のうち、Ee と Eg のみに DetectorResolution.h の応答を乱数でかけて
// 新しい 5列テキストとして保存する。
// 先頭のヘッダ行や '#' コメント行はそのままコピーする。

static bool ParseEventLine(const std::string& line,
                           double& Ee,
                           double& Eg,
                           double& t,
                           double& phi_e,
                           double& phi_g) {
    if (line.empty()) return false;

    std::istringstream iss(line);
    if (!(iss >> Ee >> Eg >> t >> phi_e >> phi_g)) return false;
    if (!std::isfinite(Ee) || !std::isfinite(Eg) || !std::isfinite(t) ||
        !std::isfinite(phi_e) || !std::isfinite(phi_g)) {
        return false;
    }
    return true;
}

static void PrintUsage(const char* argv0) {
    std::cerr
        << "usage: " << argv0 << " INPUT OUTPUT [SEED]\n"
        << "  INPUT/OUTPUT: 5-column text data\n"
        << "  SEED        : TRandom3 seed (default: 0)\n";
}

int main(int argc, char** argv) {
    if (argc < 3 || argc > 4) {
        PrintUsage(argv[0]);
        return 1;
    }

    const std::string infile = argv[1];
    const std::string outfile = argv[2];
    const unsigned int seed = (argc >= 4)
        ? static_cast<unsigned int>(std::stoul(argv[3]))
        : 0u;

    std::ifstream fin(infile);
    if (!fin) {
        std::cerr << "[smear_finaldata] failed to open input: " << infile << "\n";
        return 1;
    }

    std::ofstream fout(outfile);
    if (!fout) {
        std::cerr << "[smear_finaldata] failed to open output: " << outfile << "\n";
        return 1;
    }

    TRandom3 rng(seed);

    long long n_lines = 0;
    long long n_events = 0;
    long long n_passthrough = 0;

    fout << std::setprecision(10);

    std::string line;
    while (std::getline(fin, line)) {
        ++n_lines;

        double Ee = 0.0;
        double Eg = 0.0;
        double t = 0.0;
        double phi_e = 0.0;
        double phi_g = 0.0;
        if (!ParseEventLine(line, Ee, Eg, t, phi_e, phi_g)) {
            fout << line << "\n";
            ++n_passthrough;
            continue;
        }

        const double Ee_smeared = smear_energy_trandom3_e(rng, Ee);
        const double Eg_smeared = smear_energy_trandom3_g(rng, Eg);

        fout << Ee_smeared << "\t"
             << Eg_smeared << "\t"
             << t << "\t"
             << phi_e << "\t"
             << phi_g << "\n";
        ++n_events;
    }

    std::cout << "[smear_finaldata] input       : " << infile << "\n";
    std::cout << "[smear_finaldata] output      : " << outfile << "\n";
    std::cout << "[smear_finaldata] seed        : " << seed << "\n";
    std::cout << "[smear_finaldata] lines       : " << n_lines << "\n";
    std::cout << "[smear_finaldata] smeared evs : " << n_events << "\n";
    std::cout << "[smear_finaldata] passthrough : " << n_passthrough << "\n";
    return 0;
}
