#ifndef OMEGA_V1_CALMODULIN_HPP
#define OMEGA_V1_CALMODULIN_HPP

namespace omega::v1 {

struct CaMState {
    double M0;    // CaM sans calcium lié
    double M2N;   // Deux Ca2+ liés au domaine N
    double M2C;   // Deux Ca2+ liés au domaine C
    double M4;    // Quatre Ca2+ liés
};

struct CaMParameters {
    // Association : concentration^-2 * temps^-1
    // Dissociation : temps^-1
    double kon_0N;
    double koff_0N;

    double kon_0C;
    double koff_0C;

    double kon_N4;
    double koff_N4;

    double kon_C4;
    double koff_C4;
};

struct CaMFluxes {
    double v_0N_forward;
    double v_0N_reverse;

    double v_0C_forward;
    double v_0C_reverse;

    double v_N4_forward;
    double v_N4_reverse;

    double v_C4_forward;
    double v_C4_reverse;
};

struct CaMDerivative {
    double M0;
    double M2N;
    double M2C;
    double M4;
};

void validate_calmodulin(
    const CaMState& state,
    double calcium,
    const CaMParameters& parameters
);

CaMFluxes calculate_fluxes(
    const CaMState& state,
    double calcium,
    const CaMParameters& parameters
);

CaMDerivative calculate_derivative(
    const CaMState& state,
    double calcium,
    const CaMParameters& parameters
);

double total_calmodulin(const CaMState& state);

void advance_calmodulin(
    CaMState& state,
    double calcium,
    const CaMParameters& parameters,
    double duration
);

} // namespace omega::v1

#endif
