
from pathlib import Path

content = r'''#include "coupled_dephosphorylation.hpp"

namespace omega::v1 {

M2M3Derivative calculate_m2_m3_derivative(
    const M2M3State& s,
    const M2M3Parameters& p
) {
    const auto d2 =
        calculate_camkii_derivative(s.camkii, p.camkii);

    const DephosphorylationState dstate{
        s.PP1,
        s.camkii.KpC,
        s.camkii.KC,
        s.PP1KpC
    };

    const auto d3 =
        calculate_dephosphorylation_derivative(
            dstate,
            p.dephosphorylation
        );

    M2M3Derivative out{d2, 0.0, 0.0};

    out.camkii.KC += d3.KC;
    out.camkii.KpC += d3.KpC;

    out.PP1 = d3.PP1;
    out.PP1KpC = d3.PP1KpC;

    return out;
}

double total_m2_m3_camkii(const M2M3State& s) {
    return s.camkii.K0
         + s.camkii.KC
         + 2.0 * s.camkii.D
         + s.camkii.KpC
         + s.camkii.Kp
         + s.PP1KpC;
}

double total_m2_m3_pp1(const M2M3State& s) {
    return s.PP1 + s.PP1KpC;
}

} // namespace omega::v1
'''

Path("coupled_dephosphorylation.cpp").write_text(content, encoding="utf-8")
