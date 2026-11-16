#ifndef OMEGA_V1_INTEGRATED_NETWORK_HPP
#define OMEGA_V1_INTEGRATED_NETWORK_HPP

#include "../m2_camkii/camkii.hpp"
#include "../m3_pp1_i1/pp1_i1.hpp"
#include "../m3_camkii_dephosphorylation/camkii_dephosphorylation.hpp"

namespace omega::v1 {

/// État du réseau intégré CaMKII–PP1–I1.
///
/// Le pool de PP1 libre est partagé entre la déphosphorylation de CaMKII
/// et l’interaction PP1–I1. Les concentrations sont exprimées dans la
/// même unité que les sous-modèles.
struct IntegratedState {
    /// État CaMKII : sous-unités inactives, actives, phosphorylées
    /// et complexes.
    CaMKIIState camkii;

    /// Concentration de PP1 libre, partagée entre les réactions couplées.
    double PP1;

    /// Concentration de I1 non phosphorylé.
    double I1;

    /// Concentration de I1 phosphorylé.
    double I1p;

    /// Concentration du complexe PP1–I1p.
    double PP1I1p;

    /// Concentration du complexe PP1–KpC.
    double PP1KpC;
};

/// Paramètres du réseau intégré CaMKII–PP1–I1.
struct IntegratedParameters {
    /// Paramètres de la dynamique CaMKII.
    CaMKIIParameters camkii;

    /// Paramètres de la dynamique PP1–I1.
    PP1I1Parameters pp1_i1;

    /// Paramètres de la dynamique de déphosphorylation de CaMKII.
    DephosphorylationParameters dephosphorylation;

    /// Activité de la PKA, utilisée par le sous-modèle PP1–I1.
    double pka_activity;
};

/// Dérivée temporelle de l’état du réseau intégré.
struct IntegratedDerivative {
    /// Dérivée des variables d’état CaMKII.
    CaMKIIDerivative camkii;

    /// Dérivée de la concentration de PP1 libre.
    double PP1;

    /// Dérivée de la concentration de I1 non phosphorylé.
    double I1;

    /// Dérivée de la concentration de I1 phosphorylé.
    double I1p;

    /// Dérivée de la concentration du complexe PP1–I1p.
    double PP1I1p;

    /// Dérivée de la concentration du complexe PP1–KpC.
    double PP1KpC;
};

/// Calcule la dérivée temporelle du réseau intégré CaMKII–PP1–I1.
///
/// La dynamique CaMKII produit les variations de KpC et KC. La
/// déphosphorylation ajoute ses contributions à KpC et KC, tandis que
/// le sous-modèle PP1–I1 détermine les variations de I1, I1p et
/// PP1I1p. Les deux réactions utilisent le même pool de PP1 libre :
/// leurs contributions sont donc additionnées pour PP1.
[[nodiscard]] IntegratedDerivative calculate_integrated_derivative(
    const IntegratedState& state,
    const IntegratedParameters& parameters
);

/// Calcule la quantité totale des espèces associées à CaMKII.
///
/// Le dimère D compte pour deux sous-unités. Le complexe PP1KpC
/// séquestre une sous-unité KpC et est donc inclus dans le total.
[[nodiscard]] double total_integrated_camkii(
    const IntegratedState& state
);

/// Calcule la quantité totale de PP1.
///
/// PP1 existe sous forme libre, lié à KpC dans PP1KpC, ou lié à I1p
/// dans PP1I1p.
[[nodiscard]] double total_integrated_pp1(
    const IntegratedState& state
);

/// Calcule la quantité totale de I1.
///
/// I1 existe sous forme libre, phosphorylée sous I1p, ou liée à PP1
/// dans le complexe PP1I1p.
[[nodiscard]] double total_integrated_i1(
    const IntegratedState& state
);

} // namespace omega::v1

#endif // OMEGA_V1_INTEGRATED_NETWORK_HPP
