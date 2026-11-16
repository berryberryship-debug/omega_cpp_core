#ifndef OMEGA_V1_FULL_HPP
#define OMEGA_V1_FULL_HPP

#include "../m1_calmodulin/calmodulin.hpp"
#include "../m2_camkii/camkii.hpp"
#include "../m3_pp1_i1/pp1_i1.hpp"
#include "../m3_camkii_dephosphorylation/camkii_dephosphorylation.hpp"

namespace omega::v1 {

struct CaMKIIComplexState {
    double K0;
    double KC;
    double D;
    double KpC;
    double Kp;
};

struct CaMKIIComplexDerivative {
    double K0;
    double KC;
    double D;
    double KpC;
    double Kp;
};

struct OmegaFullState {
    CaMState cam;
    CaMKIIComplexState camkii;
    double PP1;
    double I1;
    double I1p;
    double PP1I1p;
    double PP1KpC;
};

struct OmegaFullParameters {
    CaMParameters calmodulin;
    double calcium;
    CaMKIIParameters camkii;
    PP1I1Parameters pp1_i1;
    DephosphorylationParameters dephosphorylation;
    double pka_activity;
};

struct OmegaFullDerivative {
    CaMDerivative cam;
    CaMKIIComplexDerivative camkii;
    double PP1;
    double I1;
    double I1p;
    double PP1I1p;
    double PP1KpC;
};

OmegaFullDerivative calculate_omega_full_derivative(
    const OmegaFullState& state,
    const OmegaFullParameters& parameters
);

double total_omega_cam(
    const OmegaFullState& state
);

double total_omega_camkii(
    const OmegaFullState& state
);

double total_omega_pp1(
    const OmegaFullState& state
);

double total_omega_i1(
    const OmegaFullState& state
);

} // namespace omega::v1

#endif // OMEGA_V1_FULL_HPP
