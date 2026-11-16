#ifndef OMEGA_V1_COUPLED_HPP
#define OMEGA_V1_COUPLED_HPP

#include "../m1_calmodulin/calmodulin.hpp"
#include "../m2_camkii/camkii.hpp"

namespace omega::v1 {

/// État du système couplé calmoduline–CaMKII.
///
/// La calmoduline libre M4 est stockée uniquement dans cam.M4. Elle est
/// lue par le sous-modèle CaMKII, mais n’est pas dupliquée dans l’état
/// couplé.
struct CoupledState {
    /// État de la calmoduline, incluant la forme libre M4.
    CaMState cam;

    /// Concentration de CaMKII inactive.
    double K0;

    /// Concentration de CaMKII liée à la calmoduline.
    double KC;

    /// Concentration du dimère de CaMKII.
    double D;

    /// Concentration de CaMKII phosphorylée liée à la calmoduline.
    double KpC;

    /// Concentration de CaMKII phosphorylée.
    double Kp;
};

/// Paramètres des sous-modèles calmoduline et CaMKII.
struct CoupledParameters {
    /// Paramètres de la dynamique de la calmoduline.
    CaMParameters cam;

    /// Paramètres de la dynamique de CaMKII.
    CaMKIIParameters camkii;
};

/// Dérivée temporelle de l’état du système couplé.
struct CoupledDerivative {
    /// Dérivée des états de calmoduline produits par le module M1.
    CaMDerivative cam;

    /// Dérivée de la concentration de CaMKII inactive.
    double K0;

    /// Dérivée de la concentration de CaMKII liée à la calmoduline.
    double KC;

    /// Dérivée de la concentration du dimère de CaMKII.
    double D;

    /// Dérivée de la concentration de CaMKII phosphorylée liée à la
    /// calmoduline.
    double KpC;

    /// Dérivée de la concentration de CaMKII phosphorylée.
    double Kp;

    /// Dérivée totale de la calmoduline libre M4.
    ///
    /// Elle combine la contribution du module M1 et celle du module M2.
    /// Pour mettre à jour state.cam.M4, utiliser M4_total ; ne pas
    /// additionner cam.M4 et M4_total.
    double M4_total;
};

/// Calcule la dérivée temporelle du système couplé calmoduline–CaMKII.
///
/// Le module M1 détermine l’évolution de la calmoduline, dont la forme
/// libre M4 est utilisée par le module M2. La dérivée totale de M4 est
/// la somme des contributions de M1 et M2.
[[nodiscard]] CoupledDerivative calculate_coupled_derivative(
    const CoupledState& state,
    double calcium,
    const CoupledParameters& parameters
);

/// Calcule la quantité totale des espèces associées à CaMKII.
///
/// Le dimère D compte pour deux sous-unités. Les complexes KC et KpC
/// séquestrent chacun une sous-unité et sont donc inclus dans le total.
[[nodiscard]] double total_coupled_camkii(
    const CoupledState& state
);

/// Calcule la quantité totale de calmoduline.
///
/// La calmoduline libre M4, les formes M0, M2N et M2C, ainsi que la
/// calmoduline liée aux espèces CaMKII KC, D et KpC sont incluses.
[[nodiscard]] double total_coupled_calmodulin(
    const CoupledState& state
);

} // namespace omega::v1

#endif // OMEGA_V1_COUPLED_HPP
