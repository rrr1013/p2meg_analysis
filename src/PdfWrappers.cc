// src/PdfWrappers.cc
#include "p2meg/PdfWrappers.h"

#include <cmath>
#include <limits>

// ============================================================
// 角度ヘルパ
// ============================================================

static double ThetaFromPhi(double phi_e, double phi_g)
{
    // phi は有限値を想定（不正入力は NaN を返す）
    if (!std::isfinite(phi_e) || !std::isfinite(phi_g)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double dphi_raw = std::fabs(phi_e - phi_g);
    const double two_pi = 2.0 * pi;
    double dphi = std::fmod(dphi_raw, two_pi);
    if (dphi > pi) dphi = two_pi - dphi;
    return dphi;
}

// ============================================================
// Signal ラッパ
// ============================================================

double SignalPdfEval(const Event& ev, const void* ctx)
{
    if (!ctx) return 0.0;
    const auto* c = static_cast<const SignalPdfContext*>(ctx);

    // SignalPdf は窓外や不正入力で 0 を返す仕様
    const double theta_eg = ThetaFromPhi(ev.phi_detector_e, ev.phi_detector_g);
    const double p = SignalPdf(ev.Ee, ev.Eg, ev.t, theta_eg, c->win, c->res, c->ms);
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
    const double theta_eg = ThetaFromPhi(ev.phi_detector_e, ev.phi_detector_g);
    const double cos_detector_e = std::cos(ev.phi_detector_e);
    const double cos_detector_g = std::cos(ev.phi_detector_g);
    const double p = RMDGridPdf(
        ev.Ee, ev.Eg, ev.t, theta_eg,
        cos_detector_e, cos_detector_g
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
