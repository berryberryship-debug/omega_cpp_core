
from pathlib import Path

content = r'''#include "integrated_network.hpp"

namespace omega::v1 {

IntegratedDerivative calculate_integrated_derivative(
    const IntegratedState& s,
    const IntegratedParameters& p
) {
    const auto d_camkii =
        calculate_camkii_derivative(s.camkii, p.camkii);

    const DephosphorylationState dephos_state{
        s.PP1,
        s.camkii.KpC,
        s.camkii.KC,
        s.PP1KpC
    };

    const auto d_dephos =
        calculate_dephosphorylation_derivative(
            dephos_state,
            p.dephosphorylation
        );

    const PP1I1State pp1_i1_state{
        s.PP1,
        s.I1,
        s.I1p,
        s.PP1I1p
    };

    const auto d_pp1_i1 =
        calculate_pp1_i1_derivative(
            pp1_i1_state,
            p.pp1_i1,
            p.pka_activity
        );

    IntegratedDerivative out{
        d_camkii,
        0.0,
        d_pp1_i1.I1,
        d_pp1_i1.I1p,
        d_pp1_i1.PP1I1p,
        d_dephos.PP1KpC
    };

    out.camkii.KC += d_dephos.KC;
    out.camkii.KpC += d_dephos.KpC;

    // Un seul pool libre PP1, partagé entre les deux réactions.
    out.PP1 = d_dephos.PP1 + d_pp1_i1.PP1;

    return out;
}

double total_integrated_camkii(const IntegratedState& s) {
    return s.camkii.K0
         + s.camkii.KC
         + 2.0 * s.camkii.D
         + s.camkii.KpC
         + s.camkii.Kp
         + s.PP1KpC;
}

double total_integrated_pp1(const IntegratedState& s) {
    return s.PP1 + s.PP1KpC + s.PP1I1p;
}

double total_integrated_i1(const IntegratedState& s) {
    return s.I1 + s.I1p + s.PP1I1p;
}

} // namespace omega::v1
'''

Path("integrated_network.cpp").write_text(content, encoding="utf-8")
