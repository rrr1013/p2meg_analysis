// src/PdfWrappers.cc
#include "p2meg/PdfWrappers.h"

#include <cmath>

// ============================================================
// Signal ラッパ
// ============================================================

double SignalPdfEval(const Event& ev, const void* ctx)
{
    if (!ctx) return 0.0;
    const auto* c = static_cast<const SignalPdfContext*>(ctx);

    // SignalPdf は窓外や不正入力で 0 を返す仕様
    const double p = SignalPdf(
        ev.Ee, ev.Eg, ev.t,
        ev.phi_detector_e, ev.phi_detector_g,
        c->win, c->res, c->ms
    );
    if (!std::isfinite(p) || p < 0.0) return 0.0;
    return p;
}

PdfComponent MakeSignalComponent(const SignalPdfContext* ctx)
{
    PdfComponent comp;
    comp.name = "sig";
    comp.eval = &SignalPdfEval;
    comp.ctx  = ctx;
    return comp;
}

// ============================================================
// RMD ラッパ
// ============================================================

double RMDGridPdfEval(const Event& ev, const void* /*ctx*/)
{
    // RMDGridPdf は未ロード・窓外・不正入力で 0 を返す仕様
    const double p = RMDGridPdf(
        ev.Ee, ev.Eg, ev.t,
        ev.phi_detector_e, ev.phi_detector_g
    );
    if (!std::isfinite(p) || p < 0.0) return 0.0;
    return p;
}

PdfComponent MakeRMDComponent()
{
    PdfComponent comp;
    comp.name = "rmd";
    comp.eval = &RMDGridPdfEval;
    comp.ctx  = nullptr;
    return comp;
}

// ============================================================
// ACC ラッパ
// ============================================================

double ACCGridPdfEval(const Event& ev, const void* /*ctx*/)
{
    // ACCGridPdf は未ロード・窓外・不正入力で 0 を返す仕様
    const double p = ACCGridPdf(
        ev.Ee, ev.Eg, ev.t,
        ev.phi_detector_e, ev.phi_detector_g
    );
    if (!std::isfinite(p) || p < 0.0) return 0.0;
    return p;
}

PdfComponent MakeACCComponent()
{
    PdfComponent comp;
    comp.name = "acc";
    comp.eval = &ACCGridPdfEval;
    comp.ctx  = nullptr;
    return comp;
}
