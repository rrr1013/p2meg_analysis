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
//  - 時間因子は root 内の Eg 条件付き key_tshape_egbin* テンプレートを優先し、
//    無い場合は単一 key_tshape、さらに無ければ解析窓内一様へフォールバックする。
// ============================================================

// 初期化：root ファイルから格子PDF（4D）を読み込む
//  - filepath : 例 "data/pdf_cache/acc_grid.root"
//  - key      : 例 "acc_grid"
// 成功したら true
bool ACCGridPdf_Load(const char* filepath, const char* key);

// 現在ロード済みかどうか（デバッグ・安全用）
bool ACCGridPdf_IsLoaded();

// 現在ロード中の root ファイルと key を返す
// 未ロード時は nullptr を返す
const char* ACCGridPdf_LoadedFilepath();
const char* ACCGridPdf_LoadedKey();

// PDF 評価（ロード済みが前提）
// 解析窓外は 0 を返す
// 時間因子は Eg 条件付き key_tshape_egbin* から作る p_t(t|Eg) を掛ける
// （無ければ単一 key_tshape、さらに無ければ一様）
// phi は DetectorResolution の範囲にクリップし離散化して評価する
double ACCGridPdf(double Ee, double Eg, double t,
                  double phi_detector_e, double phi_detector_g);

#endif // P2MEG_ACC_GRID_PDF_H
