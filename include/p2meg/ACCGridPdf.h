#ifndef P2MEG_ACC_GRID_PDF_H
#define P2MEG_ACC_GRID_PDF_H

// ============================================================
// ACCGridPdf
//
// 目的
//  - オフラインで作成した ACC 4D 格子 PDF を読み込み、
//    任意の観測値 (Ee, Eg, t, phi_detector_e, phi_detector_g) に対して
//    ACC の PDF 値を返す。
//  - 格子は 4D（Ee, Eg, phi_detector_e, phi_detector_g）で保存しておき、
//    評価時に theta_eg = |phi_e - phi_g| を作って解析窓カットを行う。
//  - 時間因子は解析窓内一様として解析的に掛ける。
// ============================================================

// 初期化：root ファイルから格子PDF（4D）を読み込む
//  - filepath : 例 "data/pdf_cache/acc_grid.root"
//  - key      : 例 "acc_grid"
// 成功したら true
bool ACCGridPdf_Load(const char* filepath, const char* key);

// 現在ロード済みかどうか（デバッグ・安全用）
bool ACCGridPdf_IsLoaded();

// PDF 評価（ロード済みが前提）
// 解析窓外は 0 を返す
// 時間因子は一様 (1/(t_max - t_min)) を掛ける
// phi は DetectorResolution の範囲にクリップし離散化して評価する
double ACCGridPdf(double Ee, double Eg, double t,
                  double phi_detector_e, double phi_detector_g);

#endif // P2MEG_ACC_GRID_PDF_H
