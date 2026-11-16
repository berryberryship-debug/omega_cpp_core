#include "integrated_network.hpp"

namespace omega::v1 {

IntegratedDerivative calculate_integrated_derivative(
    const IntegratedState& state,
    const IntegratedParameters& parameters
) {
    const CaMKIIDerivative camkii_rate =
        calculate_camkii_derivative(
            state.camkii,
            parameters.camkii
        );

    const DephosphorylationState dephosphorylation_state{
        state.PP1,
        state.camkii.KpC,
        state.camkii.KC,
        state.PP1KpC
    };

    const DephosphorylationDerivative dephosphorylation_rate =
        calculate_dephosphorylation_derivative(
            dephosphorylation_state,
            parameters.dephosphorylation
        );

    const PP1I1State pp1_i1_state{
        state.PP1,
        state.I1,
        state.I1p,
        state.PP1I1p
    };

    const PP1I1Derivative pp1_i1_rate =
        calculate_pp1_i1_derivative(
            pp1_i1_state,
            parameters.pp1_i1,
            parameters.pka_activity
        );

    IntegratedDerivative rate{
        camkii_rate,
        dephosphorylation_rate.PP1 + pp1_i1_rate.PP1,
        pp1_i1_rate.I1,
        pp1_i1_rate.I1p,
        pp1_i1_rate.PP1I1p,
        dephosphorylation_rate.PP1KpC
    };

    // La déphosphorylation s’ajoute à la dynamique CaMKII pour les
    // sous-unités KpC et KC.
    rate.camkii.KpC += dephosphorylation_rate.KpC;
    rate.camkii.KC += dephosphorylation_rate.KC;

    return rate;
}

double total_integrated_camkii(const IntegratedState& state) {
    const CaMKIIState& camkii = state.camkii;

    // K0, KC, KpC et Kp comptent pour une sous-unité chacun.
    // D représente un dimère et compte donc pour deux sous-unités.
    // PP1KpC séquestre une sous-unité KpC.
    return camkii.K0
         + camkii.KC
         + 2.0 * camkii.D
         + camkii.KpC
         + camkii.Kp
         + state.PP1KpC;
}

double total_integrated_pp1(const IntegratedState& state) {
    // PP1 est libre, lié à KpC dans PP1KpC, ou lié à I1p dans PP1I1p.
    return state.PP1
         + state.PP1KpC
         + state.PP1I1p;
}

double total_integrated_i1(const IntegratedState& state) {
    // I1 est libre, phosphorylé sous I1p, ou lié à PP1 dans PP1I1p.
    return state.I1
         + state.I1p
         + state.PP1I1p;
}

} // namespace omega::v1
