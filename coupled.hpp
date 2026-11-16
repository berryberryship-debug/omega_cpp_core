#ifndef OMEGA_V1_COUPLED_HPP
#define OMEGA_V1_COUPLED_HPP

#include "../m1_calmodulin/calmodulin.hpp"
#include "../m2_camkii/camkii.hpp"

namespace omega::v1 {

struct CoupledState {
    CaMState cam;

    // La calmoduline libre M4 est stockée uniquement dans cam.M4.
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
    // Contribution du module M1 aux états de calmoduline libre.
    CaMDerivative cam;

    double K0;
    double KC;
    double D;
    double KpC;
    double Kp;

    // Dérivée totale de M4 : contribution M1 + contribution M2.
    // Pour mettre à jour state.cam.M4, utiliser M4_total,
    // et non cam.M4 ; ne pas additionner les deux.
    double M4_total;
};

CoupledDerivative calculate_coupled_derivative(
    const CoupledState& state,
    double calcium,
    const CoupledParameters& parameters
);

double total_coupled_camkii(
    const CoupledState& state
);

double total_coupled_calmodulin(
    const CoupledState& state
);

} // namespace omega::v1

#endif // OMEGA_V1_COUPLED_HPP
