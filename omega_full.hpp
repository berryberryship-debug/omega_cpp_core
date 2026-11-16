#include "omega_full.hpp"
#include <algorithm> // Pour std::max

namespace omega::v1 {

OmegaFullDerivative calculate_omega_full_derivative(
    const OmegaFullState& s,
    const OmegaFullParameters& p
) {
    const CaMDerivative d_cam = calculate_derivative(s.cam, p.calcium, p.calmodulin);

    // M4 libre provient exclusivement du pool M1.
    const CaMKIIState kinase_state{
        s.camkii.K0,
        s.camkii.KC,
        s.camkii.D,
        s.camkii.KpC,
        s.camkii.Kp,
        s.cam.M4
    };

    const CaMKIIDerivative d_kinase = calculate_camkii_derivative(kinase_state, p.camkii);

    const DephosphorylationState dephos_state{
        s.PP1,
        s.camkii.KpC,
        s.camkii.KC,
        s.PP1KpC
    };

    const DephosphorylationDerivative d_dephos = calculate_dephosphorylation_derivative(dephos_state, p.dephosphorylation);

    const PP1I1State pp1_i1_state{
        s.PP1,
        s.I1,
        s.I1p,
        s.PP1I1p
    };

    const PP1I1Derivative d_pp1_i1 = calculate_pp1_i1_derivative(pp1_i1_state, p.pp1_i1, p.pka_activity);

    // --- INITIALISATION EXPLICITE (Plus robuste) ---
    OmegaFullDerivative out{};

    out.cam = d_cam;

    out.camkii.K0  = d_kinase.K0;
    out.camkii.KC  = d_kinase.KC + d_dephos.KC;
    out.camkii.D   = d_kinase.D;
    out.camkii.KpC = d_kinase.KpC + d_dephos.KpC;
    out.camkii.Kp  = d_kinase.Kp;

    out.PP1     = d_dephos.PP1 + d_pp1_i1.PP1;
    out.I1      = d_pp1_i1.I1;
    out.I1p     = d_pp1_i1.I1p;
    out.PP1I1p  = d_pp1_i1.PP1I1p;
    out.PP1KpC  = d_dephos.PP1KpC;

    // Contribution nette de M2 au pool de M4 libre.
    out.cam.M4 += d_kinase.M4;

    return out;
}

/**
 * @brief Calcule la concentration totale de Calmoduline dans le système.
 * 
 * Invariant de conservation : Cette valeur DOIT rester constante 
 * au cours du temps pour valider la stabilité de l'intégration numérique.
 * 
 * Note sur la stœchiométrie : Le dimère de CaMKII (D) lie deux 
 * sous-unités de calmoduline, d'où le facteur 2.0.
 */
double total_omega_cam(const OmegaFullState& s) {
    return s.cam.M0
         + s.cam.M2N
         + s.cam.M2C
         + s.cam.M4
         + s.camkii.KC
         + 2.0 * s.camkii.D
         + s.camkii.KpC; // On doit aussi inclure la CaM liée à la forme phosphorylée !
         // ATTENTION : s.ca a été retiré. Le calcium libre n'est pas de la calmoduline.
}
