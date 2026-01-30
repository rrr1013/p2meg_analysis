// include/p2meg/MichelEData.h
#ifndef P2MEG_MICHEL_E_DATA_H
#define P2MEG_MICHEL_E_DATA_H

#include <string>
#include <vector>

// ============================================================
// Michel 偏極測定（別実験）用データ形式: Ee-only
//
// ファイル形式:
//  - 1行に1事象
//  - 1列: Ee [MeV]
//  - 空行は無視
//  - 先頭が '#' の行はコメントとして無視
//
// エラーハンドリング方針:
//  - 読めない行・数値化できない行はスキップ
//  - ファイルが開けない等の致命的エラー時は空 vector を返す
// ============================================================

struct MichelEEvent {
    double Ee; // [MeV]
};

// Ee-only データを読み込む。
// n_read / n_skipped は nullptr なら更新しない。
std::vector<MichelEEvent> ReadMichelEData(const std::string& path,
                                         long long* n_read = nullptr,
                                         long long* n_skipped = nullptr);

#endif // P2MEG_MICHEL_E_DATA_H
