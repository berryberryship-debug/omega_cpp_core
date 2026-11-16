#ifndef OMEGA_V1_CAMKII_HPP
#define OMEGA_V1_CAMKII_HPP

namespace omega::v1 {

/// État de CaMKII couplé à la calmoduline libre M4.
struct CaMKIIState {
    /// CaMKII inactive.
    double K0;

    /// CaMKII liée à la calmoduline.
    double KC;

    /// Dimère de CaMKII.
    double D;

    /// CaMKII phosphorylée liée à la calmoduline.
    double KpC;

    /// CaMKII phosphorylée.
    double Kp;

    /// Calmoduline libre M4 disponible pour la liaison.
    double M4;
};

/// Paramètres cinétiques du sous-modèle CaMKII.
struct CaMKIIParameters {
    /// Constante de liaison CaMKII–calmoduline.
    double kon_bind;

    /// Constante de dissociation du complexe CaMKII–calmoduline.
    double koff_bind;

    /// Constante de formation du dimère.
    double kon_dimer;

    /// Constante de dissociation du dimère.
    double koff_dimer;

    /// Constante de phosphorylation.
    double k_phosphorylation;

    /// Constante de dissociation de la calmoduline liée à CaMKII phosphorylée.
    double koff_phospho_cam;

    /// Constante de liaison de la calmoduline à CaMKII phosphorylée.
    double kon_phospho_cam;
};

/// Dérivée temporelle de l’état CaMKII.
struct CaMKIIDerivative {
    double K0;
    double KC;
    double D;
    double KpC;
    double Kp;
    double M4;
};

/// Flux élémentaires du sous-modèle CaMKII.
struct CaMKIIRates {
    /// Flux de liaison CaMKII–calmoduline.
    double bind_forward;

    /// Flux de dissociation du complexe CaMKII–calmoduline.
    double bind_reverse;

    /// Flux de formation du dimère.
    double dimer_forward;

    /// Flux de dissociation du dimère.
    double dimer_reverse;

    /// Flux de phosphorylation.
    double phosphorylation;

    /// Flux de dissociation de la calmoduline liée à CaMKII phosphorylée.
    double phospho_cam_forward;

    /// Flux de liaison de la calmoduline à CaMKII phosphorylée.
    double phospho_cam_reverse;
};

/// Calcule les flux élémentaires de CaMKII.
[[nodiscard]] CaMKIIRates calculate_camkii_rates(
    const CaMKIIState& state,
    const CaMKIIParameters& parameters
);

/// Calcule la dérivée temporelle de l’état CaMKII.
[[nodiscard]] CaMKIIDerivative calculate_camkii_derivative(
    const CaMKIIState& state,
    const CaMKIIParameters& parameters
);

/// Calcule la quantité totale de sous-unités CaMKII.
///
/// Le dimère D compte pour deux sous-unités.
[[nodiscard]] double total_camkii(const CaMKIIState& state);

/// Valide l’état et les paramètres de CaMKII.
void validate_camkii(
    const CaMKIIState& state,
    const CaMKIIParameters& parameters
);

} // namespace omega::v1

#endif // OMEGA_V1_CAMKII_HPP
