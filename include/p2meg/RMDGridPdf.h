// include/p2meg/RMDGridPdf.h
#ifndef P2MEG_RMD_GRID_PDF_H
#define P2MEG_RMD_GRID_PDF_H

// ============================================================
// RMDGridPdf
//
// 目的
//  - オフラインで作成した「RMD 畳み込み済み 格子PDF」を読み込み、
//    任意の観測値に対して PDF 値を返す。
//
// 入力（Eval）
//  - Ee            : 陽電子エネルギー [MeV]
//  - Eg            : ガンマ線エネルギー [MeV]
//  - t             : 到達時間差 Δt [ns]
//  - theta         : （互換用。内部では使用しない想定）
//  - cos_detector_e: cosθ_e = P̂ · d̂_e（検出器代表方向）
//  - cos_detector_g: cosθ_γ = P̂ · d̂_γ（検出器代表方向）
//
// 角度の扱い（重要）
//  - cos_detector_e/g は「角度の離散化」に対応して最近傍の格子点に丸める。
//    DetectorResolutionConst::N_theta により 0..pi を N_theta 分割し、
//      theta_i = i*pi/N_theta (i=0..N_theta)
//    として、ie, ig を決める（cos は acos→丸め→cos に対応）。
//  - 装置幾何の仮定：偏極軸・e検出器方向・γ検出器方向が同一平面上。
//    さらに cosΔφ=+1（Δφ=0）で固定。
//    よって相対角は
//      cosθ_eγ = cosθ_e cosθ_γ + sinθ_e sinθ_γ = cos(theta_e - theta_γ)
//    とし、解析窓の theta カットはこの θ_eγ で判定する。
//
// 実装方針（重要）
//  - root から読み込む格子は 4D（Ee, Eg, ie, ig）の THnD。
//    ここで ie, ig は離散インデックス（0..N_theta）であり、補間しない。
//  - Ee,Eg 方向のみ 2D 双一次補間（4点）で p2(Ee,Eg|ie,ig) を求める。
//  - 時間因子は解析的に掛ける（窓内正規化ガウシアン）：
//      p(t) = N(t_mean, sigma_t) / ∫_{tmin}^{tmax} N dt
// ============================================================

// 初期化：root ファイルから格子PDF（4D）を読み込む
//  - filepath : 例 "data/pdf_cache/rmd_grid.root"
//  - key      : 例 "rmd_grid"
// 成功したら true
bool RMDGridPdf_Load(const char* filepath, const char* key);

// 現在ロード済みかどうか（デバッグ・安全用）
bool RMDGridPdf_IsLoaded();

// PDF 評価（ロード済みが前提）
// 解析窓外は 0 を返す
double RMDGridPdf(double Ee, double Eg, double t, double theta,
                  double cos_detector_e, double cos_detector_g);

#endif // P2MEG_RMD_GRID_PDF_H
