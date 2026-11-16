#ifndef OMEGA_V1_INTEGRATED_NETWORK_HPP
#define OMEGA_V1_INTEGRATED_NETWORK_HPP

#include "../m2_camkii/camkii.hpp"
#include "../m3_pp1_i1/pp1_i1.hpp"
#include "../m3_camkii_dephosphorylation/camkii_dephosphorylation.hpp"

namespace omega::v1 {

struct IntegratedState {
    CaMKIIState camkii;
    double PP1;
    double I1;
    double I1p;
    double PP1I1p;
    double PP1KpC;
};

struct IntegratedParameters {
    CaMKIIParameters camkii;
    PP1I1Parameters pp1_i1;
    DephosphorylationParameters dephosphorylation;
    double pka_activity;
};

struct IntegratedDerivative {
    CaMKIIDerivative camkii;
    double PP1;
    double I1;
    double I1p;
    double PP1I1p;
    double PP1KpC;
};

IntegratedDerivative calculate_integrated_derivative(
    const IntegratedState& state,
    const IntegratedParameters& parameters
);

double total_integrated_camkii(
    const IntegratedState& state
);

double total_integrated_pp1(
    const IntegratedState& state
);

double total_integrated_i1(
    const IntegratedState& state
);

} // namespace omega::v1

#endif // OMEGA_V1_INTEGRATED_NETWORK_HPP
