
from pathlib import Path

content = r'''#include "omega_full.hpp"

namespace omega::v1 {

OmegaFullDerivative calculate_omega_full_derivative(
    const OmegaFullState& s,
    const OmegaFullParameters& p
) {
    const auto d_cam = calculate_derivative(
        s.cam, p.calcium, p.calmodulin
    );

    // M4 libre provient exclusivement du pool M1.
    const CaMKIIState kinase_state{
        s.camkii.K0,
        s.camkii.KC,
        s.camkii.D,
        s.camkii.KpC,
        s.camkii.Kp,
        s.cam.M4
    };

    const auto d_kinase =
        calculate_camkii_derivative(kinase_state, p.camkii);

    const DephosphorylationState dephos_state{
        s.PP1,
        s.camkii.KpC,
        s.camkii.KC,
        s.PP1KpC
    };

    const auto d_dephos =
        calculate_dephosphorylation_derivative(
            dephos_state, p.dephosphorylation
        );

    const PP1I1State pp1_i1_state{
        s.PP1, s.I1, s.I1p, s.PP1I1p
    };

    const auto d_pp1_i1 =
        calculate_pp1_i1_derivative(
            pp1_i1_state, p.pp1_i1, p.pka_activity
        );

    OmegaFullDerivative out{
        d_cam,
        {
            d_kinase.K0,
            d_kinase.KC + d_dephos.KC,
            d_kinase.D,
            d_kinase.KpC + d_dephos.KpC,
            d_kinase.Kp
        },
        d_dephos.PP1 + d_pp1_i1.PP1,
        d_pp1_i1.I1,
        d_pp1_i1.I1p,
        d_pp1_i1.PP1I1p,
        d_dephos.PP1KpC
    };

    // Échange de CaM libre avec CaMKII :
    // liaison K0 + M4 -> KC et dissociation KpC -> Kp + M4.
    out.cam.M4 += d_kinase.M4;

    return out;
}

double total_omega_cam(const OmegaFullState& s) {
    return s.cam.M0 + s.cam.M2N + s.cam.M2C + s.cam.M4
         + s.camkii.KC + 2.0 * s.camkii.D
         + s.camkii.KpC + s.PP1KpC;
}

double total_omega_camkii(const OmegaFullState& s) {
    return s.camkii.K0 + s.camkii.KC + 2.0 * s.camkii.D
         + s.camkii.KpC + s.camkii.Kp + s.PP1KpC;
}

double total_omega_pp1(const OmegaFullState& s) {
    return s.PP1 + s.PP1I1p + s.PP1KpC;
}

double total_omega_i1(const OmegaFullState& s) {
    return s.I1 + s.I1p + s.PP1I1p;
}

} // namespace omega::v1
'''

Path("omega_full.cpp").write_text(content, encoding="utf-8")
