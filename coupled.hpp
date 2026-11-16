
from pathlib import Path

content = r'''#ifndef OMEGA_V1_COUPLED_HPP
#define OMEGA_V1_COUPLED_HPP

#include "../m1_calmodulin/calmodulin.hpp"
#include "../m2_camkii/camkii.hpp"

namespace omega::v1 {

struct CoupledState {
    CaMState cam;

    // M4 n'est pas dupliqué ici.
    double K0;
    double KC;
    double D;
    double KpC;
    double Kp;
};

struct CoupledParameters {
    CaMParameters cam;
    CaMKIIParameters camkii;
};

struct CoupledDerivative {
    CaMDerivative cam;

    double K0;
    double KC;
    double D;
    double KpC;
    double Kp;

    // Dérivée totale de M4 : M1 + contribution M2.
    double M4_total;
};

CoupledDerivative calculate_coupled_derivative(
    const CoupledState& state,
    double calcium,
    const CoupledParameters& parameters
);

double total_coupled_camkii(const CoupledState& state);
double total_coupled_calmodulin(const CoupledState& state);

} // namespace omega::v1

#endif
'''

Path("coupled.hpp").write_text(content, encoding="utf-8")
