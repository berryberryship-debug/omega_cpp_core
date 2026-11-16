#ifndef OMEGA_V1_COUPLED_DEPHOSPHORYLATION_HPP
#define OMEGA_V1_COUPLED_DEPHOSPHORYLATION_HPP

#include "../m2_camkii/camkii.hpp"
#include "../m3_camkii_dephosphorylation/camkii_dephosphorylation.hpp"

namespace omega::v1 {

struct M2M3State {
    CaMKIIState camkii;
    double PP1;
    double PP1KpC;
};

struct M2M3Derivative {
    CaMKIIDerivative camkii;
    double PP1;
    double PP1KpC;
};

struct M2M3Parameters {
    CaMKIIParameters camkii;
    DephosphorylationParameters dephosphorylation;
};

M2M3Derivative calculate_m2_m3_derivative(
    const M2M3State& state,
    const M2M3Parameters& parameters
);

double total_m2_m3_camkii(
    const M2M3State& state
);

double total_m2_m3_pp1(
    const M2M3State& state
);

} // namespace omega::v1

#endif // OMEGA_V1_COUPLED_DEPHOSPHORYLATION_HPP
