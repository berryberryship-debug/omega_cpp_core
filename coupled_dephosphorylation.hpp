#ifndef OMEGA_V1_COUPLED_DEPHOSPHORYLATION_HPP
#define OMEGA_V1_COUPLED_DEPHOSPHORYLATION_HPP

#include "../m2_camkii/camkii.hpp"
#include "../m3_camkii_dephosphorylation/camkii_dephosphorylation.hpp"

namespace omega::v1 {

/// État du modèle couplé CaMKII–PP1.
struct M2M3State {
    CaMKIIState camkii;
    double PP1;
    double PP1KpC;
};

/// Dérivée temporelle de l’état du modèle couplé CaMKII–PP1.
struct M2M3Derivative {
    CaMKIIDerivative camkii;
    double PP1;
    double PP1KpC;
};

/// Paramètres des sous-modèles CaMKII et de déphosphorylation.
struct M2M3Parameters {
    CaMKIIParameters camkii;
    DephosphorylationParameters dephosphorylation;
};

/// Calcule la dérivée temporelle du système couplé CaMKII–PP1.
[[nodiscard]] M2M3Derivative calculate_m2_m3_derivative(
    const M2M3State& state,
    const M2M3Parameters& parameters
);

/// Calcule la quantité totale des espèces associées à CaMKII.
[[nodiscard]] double total_m2_m3_camkii(
    const M2M3State& state
);

/// Calcule la quantité totale de PP1.
[[nodiscard]] double total_m2_m3_pp1(
    const M2M3State& state
);

} // namespace omega::v1

#endif // OMEGA_V1_COUPLED_DEPHOSPHORYLATION_HPP
