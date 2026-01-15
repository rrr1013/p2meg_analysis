// include/p2meg/PdfWrappers.h
#ifndef P2MEG_PDF_WRAPPERS_H
#define P2MEG_PDF_WRAPPERS_H

#include "p2meg/Event.h"
#include "p2meg/Likelihood.h" // PdfEval, PdfComponent

#include "p2meg/AnalysisWindow.h"
#include "p2meg/DetectorResolution.h"
#include "p2meg/Constants.h"

// SignalPdf / RMDGridPdf の本体ヘッダ
#include "p2meg/SignalPdf.h"
#include "p2meg/RMDGridPdf.h"

// ============================================================
// p2MEG PDF ラッパ
//
// 目的:
//  - 既存の PDF 関数（SignalPdf, RMDGridPdf）を、尤度側の共通I/F
//    PdfComponent (PdfEval + ctx) に接続する。
//  - 尤度側は「成分の配列」を足し算するだけにする。
//
// 注意:
//  - ここで返す p(x) は「解析窓内で正規化された密度」を想定。
//    SignalPdf / RMDGridPdf はその前提を満たしている。
// ============================================================

// ---- Signal 用 ctx ----
struct SignalPdfContext {
    AnalysisWindow4D win;          // 解析窓
    DetectorResolutionConst res;   // 分解能
    ParticleMasses ms;             // 質量（Ee0=Eg0=m_mu/2 の計算に使用）
};

// Signal 成分評価（PdfEval 互換）
double SignalPdfEval(const Event& ev, const void* ctx);

// Signal 成分の PdfComponent を作る（ctx は呼び出し側が生存管理）
PdfComponent MakeSignalComponent(const SignalPdfContext* ctx);

// ---- RMD 用 ----
// RMDGridPdf は内部でロード済み格子を参照するため、通常 ctx は不要。
// 角度情報は Event の (phi_detector_e, phi_detector_g) をそのまま渡す。
double RMDGridPdfEval(const Event& ev, const void* ctx);

// RMD 成分の PdfComponent を作る（ctx は NULL でよい）
PdfComponent MakeRMDComponent();

#endif // P2MEG_PDF_WRAPPERS_H
