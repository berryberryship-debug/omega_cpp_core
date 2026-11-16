#ifndef OMEGA_V1_COUPLED_DEPHOSPHORYLATION_HPP
#define OMEGA_V1_COUPLED_DEPHOSPHORYLATION_HPP

#include "../m2_camkii/camkii.hpp"
#include "../m3_camkii_dephosphorylation/camkii_dephosphorylation.hpp"

namespace omega::v1 {

/// État du modèle couplé CaMKII–PP1.
///
/// Les concentrations sont supposées exprimées dans la même unité que
/// les sous-modèles CaMKII et déphosphorylation, généralement en µM.
struct M2M3State {
    /// État CaMKII : formes inactives, actives, phosphorylées et complexes.
    CaMKIIState camkii;

    /// PP1 libre.
    double PP1;

    /// Complexe PP1–KpC.
    double PP1KpC;
};

/// Dérivée temporelle de l’état du modèle couplé.
struct M2M3Derivative {
    /// Dérivée des variables CaMKII.
    CaMKIIDerivative camkii;

    /// Dérivée de la concentration de PP1 libre.
    double PP1;

    /// Dérivée de la concentration du complexe PP1–KpC.
    double PP1KpC;
};

/// Paramètres des deux sous-modèles couplés.
struct M2M3Parameters {
    /// Paramètres de la dynamique CaMKII.
    CaMKIIParameters camkii;

    /// Paramètres de la dynamique de déphosphorylation.
    DephosphorylationParameters dephosphorylation;
};

/// Calcule la dérivée temporelle du système couplé CaMKII–PP1.
///
/// La dynamique CaMKII fournit les variations de KpC et KC. La
/// déphosphorylation modifie ensuite KpC, KC, PP1 et PP1KpC. Les
/// contributions des deux sous-modèles sont donc additionnées pour KpC
/// et KC.
[[nodiscard]] M2M3Derivative calculate_m2_m3_derivative(
    const M2M3State& state,
    const M2M3Parameters& parameters
);

/// Calcule la quantité totale conservée des espèces associées à CaMKII.
///
/// Le complexe PP1KpC est inclus car il séquestre une sous-unité KpC.
[[nodiscard]] double total_m2_m3_camkii(const M2M3State& state);

/// Calcule la quantité totale conservée de PP1.
///
/// PP1 existe sous forme libre ou liée à KpC dans le complexe PP1KpC.
[[nodiscard]] double total_m2_m3_pp1(const M2M3State& state);

} // namespace omega::v1

#endif // OMEGA_V1_COUPLED_DEPHOSPHORYLATION_HPP
