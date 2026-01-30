// src/MichelEData.cc
#include "p2meg/MichelEData.h"

#include <fstream>
#include <sstream>
#include <string>
#include <cctype>

static inline bool IsBlankLine(const std::string& s) {
    for (char c : s) {
        if (!std::isspace(static_cast<unsigned char>(c))) return false;
    }
    return true;
}

static inline bool IsCommentLine(const std::string& s) {
    for (char c : s) {
        if (std::isspace(static_cast<unsigned char>(c))) continue;
        return (c == '#');
    }
    return false;
}

std::vector<MichelEEvent> ReadMichelEData(const std::string& path,
                                         long long* n_read,
                                         long long* n_skipped) {
    if (n_read) *n_read = 0;
    if (n_skipped) *n_skipped = 0;

    std::vector<MichelEEvent> out;

    std::ifstream ifs(path);
    if (!ifs) {
        // 致命: ファイルが開けない → 空で返す（例外は投げない）
        return out;
    }

    std::string line;
    while (std::getline(ifs, line)) {
        if (IsBlankLine(line) || IsCommentLine(line)) {
            if (n_skipped) (*n_skipped)++;
            continue;
        }

        std::istringstream iss(line);
        double Ee = 0.0;
        if (!(iss >> Ee)) {
            if (n_skipped) (*n_skipped)++;
            continue;
        }

        // 余分なトークンがあっても無視（Ee-only として最初の数値だけ採用）
        MichelEEvent ev;
        ev.Ee = Ee;
        out.push_back(ev);
        if (n_read) (*n_read)++;
    }

    return out;
}
