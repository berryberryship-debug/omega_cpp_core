#include "coupled.hpp"

namespace omega::v1 {

CoupledDerivative calculate_coupled_derivative(
    const CoupledState& state,
    double calcium,
    const CoupledParameters& parameters
) {
    // M1 : évolution de la calmoduline libre.
    const CaMDerivative calmodulin_rate =
        calculate_derivative(
            state.cam,
            calcium,
            parameters.cam
        );

    // M2 : M4 est lu depuis M1, sans seconde copie persistante.
    const CaMKIIState kinase{
        state.K0,
        state.KC,
        state.D,
        state.KpC,
        state.Kp,
        state.cam.M4
    };

    const CaMKIIDerivative kinase_rate =
        calculate_camkii_derivative(
            kinase,
            parameters.camkii
        );

    return {
        calmodulin_rate,
        kinase_rate.K0,
        kinase_rate.KC,
        kinase_rate.D,
        kinase_rate.KpC,
        kinase_rate.Kp,
        calmodulin_rate.M4 + kinase_rate.M4
    };
}

double total_coupled_camkii(const CoupledState& state) {
    const CaMKIIState kinase{
        state.K0,
        state.KC,
        state.D,
        state.KpC,
        state.Kp,
        state.cam.M4
    };

    return total_camkii(kinase);
}

double total_coupled_calmodulin(const CoupledState& state) {
    // M0, M2N, M2C et M4 sont des formes libres de calmoduline.
    // KC et KpC contiennent chacun une calmoduline liée.
    // D est un dimère et contient donc deux calmodulines liées.
    return state.cam.M0
         + state.cam.M2N
         + state.cam.M2C
         + state.cam.M4
         + state.KC
         + 2.0 * state.D
         + state.KpC;
}

} // namespace omega::v1
