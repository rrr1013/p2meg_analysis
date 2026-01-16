#include <algorithm>
#include <chrono>
#include <cctype>
#include <fstream>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <vector>

// Compile: g++ -std=c++17 -O2 scripts/merge_mockdata.cc -o build/merge_mockdata
// Run example: build/merge_mockdata data/mockdata/pdfmix_10_300000.dat data/mockdata/acc_700000.dat -o data/mockdata/testdata1.dat --seed 42

namespace {

std::string Trim(const std::string& text) {
    size_t begin = 0;
    while (begin < text.size() && std::isspace(static_cast<unsigned char>(text[begin]))) {
        ++begin;
    }
    size_t end = text.size();
    while (end > begin && std::isspace(static_cast<unsigned char>(text[end - 1]))) {
        --end;
    }
    return text.substr(begin, end - begin);
}

std::vector<std::string> SplitTokens(const std::string& text) {
    std::istringstream iss(text);
    std::vector<std::string> tokens;
    std::string token;
    while (iss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

bool IsNumericToken(const std::string& token) {
    char* end_ptr = nullptr;
    const char* cstr = token.c_str();
    std::strtod(cstr, &end_ptr);
    return end_ptr != cstr && *end_ptr == '\0';
}

std::vector<std::string> ReadFiveColumnData(const std::string& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("failed to open: " + path);
    }
    std::vector<std::string> cleaned;
    std::string raw_line;
    while (std::getline(input, raw_line)) {
        const std::string line = Trim(raw_line);
        if (line.empty()) {
            continue;
        }
        auto tokens = SplitTokens(line);
        if (tokens.size() != 5) {
            continue;
        }
        bool all_numeric = true;
        for (const auto& t : tokens) {
            if (!IsNumericToken(t)) {
                all_numeric = false;
                break;
            }
        }
        if (!all_numeric) {
            continue;
        }
        cleaned.push_back(tokens[0] + "\t" + tokens[1] + "\t" + tokens[2] + "\t" + tokens[3] + "\t" + tokens[4]);
    }
    return cleaned;
}

struct Options {
    std::vector<std::string> inputs;
    std::string output;
    bool has_output = false;
    std::optional<unsigned long long> seed;
};

Options ParseArguments(int argc, char* argv[]) {
    if (argc < 3) {
        throw std::runtime_error("usage: merge_mockdata <input1> <input2> [-o output] [--seed N]");
    }
    Options options;
    int idx = 1;
    while (idx < argc) {
        std::string arg = argv[idx];
        if (arg == "-o" || arg == "--output") {
            if (idx + 1 >= argc) {
                throw std::runtime_error(arg + " requires an argument");
            }
            options.output = argv[++idx];
            options.has_output = true;
        } else if (arg == "--seed") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--seed requires a numeric argument");
            }
            options.seed = std::stoull(argv[++idx]);
        } else {
            options.inputs.push_back(arg);
        }
        ++idx;
    }
    if (options.inputs.size() != 2) {
        throw std::runtime_error("exactly two input files must be provided");
    }
    return options;
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        const Options options = ParseArguments(argc, argv);

        std::vector<std::string> merged;
        for (const auto& input : options.inputs) {
            const auto lines = ReadFiveColumnData(input);
            merged.insert(merged.end(), lines.begin(), lines.end());
        }

        unsigned long long seed_value = options.seed.has_value()
                                            ? options.seed.value()
                                            : static_cast<unsigned long long>(
                                                  std::chrono::high_resolution_clock::now()
                                                      .time_since_epoch()
                                                      .count());
        std::mt19937_64 rng(seed_value);
        std::shuffle(merged.begin(), merged.end(), rng);

        if (options.has_output) {
            std::ofstream out(options.output);
            if (!out) {
                throw std::runtime_error("cannot open output: " + options.output);
            }
            for (const auto& line : merged) {
                out << line << '\n';
            }
        } else {
            for (const auto& line : merged) {
                std::cout << line << '\n';
            }
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << '\n';
        return 1;
    }
}
