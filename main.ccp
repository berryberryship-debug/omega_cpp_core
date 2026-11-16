#include <iostream>
#include <fstream>
#include <iomanip>
#include "omega_full.hpp"

using namespace omega::v1;

// --- FONCTION UTILITAIRE POUR RK4 ---
// Ajoute un multiple d'une dérivée à un état.
// Permet d'éviter de répéter 14 fois les additions pour chaque sous-étape.
OmegaFullState add_scaled(const OmegaFullState& s, const OmegaFullDerivative& d, double factor) {
    OmegaFullState out = s;
    out.cam.M0 += d.cam.M0 * factor;
    out.cam.M2N += d.cam.M2N * factor;
    out.cam.M2C += d.cam.M2C * factor;
    out.cam.M4 += d.cam.M4 * factor;
    
    out.camkii.K0 += d.camkii.K0 * factor;
    out.camkii.KC += d.camkii.KC * factor;
    out.camkii.D += d.camkii.D * factor;
    out.camkii.KpC += d.camkii.KpC * factor;
    out.camkii.Kp += d.camkii.Kp * factor;
    
    out.PP1 += d.PP1 * factor;
    out.I1 += d.I1 * factor;
    out.I1p += d.I1p * factor;
    out.PP1I1p += d.PP1I1p * factor;
    out.PP1KpC += d.PP1KpC * factor;
    return out;
}

// --- VRAI INTÉGRATEUR RUNGE-KUTTA 4 (RK4) ---
OmegaFullState rk4_step(const OmegaFullState& s, const OmegaFullParameters& p, double dt) {
    // Étape 1 : k1
    OmegaFullDerivative k1 = calculate_omega_full_derivative(s, p);
    
    // Étape 2 : k2 (évaluation à mi-pas avec k1)
    OmegaFullState s2 = add_scaled(s, k1, dt * 0.5);
    OmegaFullDerivative k2 = calculate_omega_full_derivative(s2, p);
    
    // Étape 3 : k3 (évaluation à mi-pas avec k2)
    OmegaFullState s3 = add_scaled(s, k2, dt * 0.5);
    OmegaFullDerivative k3 = calculate_omega_full_derivative(s3, p);
    
    // Étape 4 : k4 (évaluation au pas complet avec k3)
    OmegaFullState s4 = add_scaled(s, k3, dt);
    OmegaFullDerivative k4 = calculate_omega_full_derivative(s4, p);
    
    // Combinaison finale : S_next = S + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
    OmegaFullState s_next = s;
    double dt6 = dt / 6.0;
    
    s_next.cam.M0 += dt6 * (k1.cam.M0 + 2.0 * k2.cam.M0 + 2.0 * k3.cam.M0 + k4.cam.M0);
    s_next.cam.M2N += dt6 * (k1.cam.M2N + 2.0 * k2.cam.M2N + 2.0 * k3.cam.M2N + k4.cam.M2N);
    s_next.cam.M2C += dt6 * (k1.cam.M2C + 2.0 * k2.cam.M2C + 2.0 * k3.cam.M2C + k4.cam.M2C);
    s_next.cam.M4 += dt6 * (k1.cam.M4 + 2.0 * k2.cam.M4 + 2.0 * k3.cam.M4 + k4.cam.M4);

    s_next.camkii.K0 += dt6 * (k1.camkii.K0 + 2.0 * k2.camkii.K0 + 2.0 * k3.camkii.K0 + k4.camkii.K0);
    s_next.camkii.KC += dt6 * (k1.camkii.KC + 2.0 * k2.camkii.KC + 2.0 * k3.camkii.KC + k4.camkii.KC);
    s_next.camkii.D += dt6 * (k1.camkii.D + 2.0 * k2.camkii.D + 2.0 * k3.camkii.D + k4.camkii.D);
    s_next.camkii.KpC += dt6 * (k1.camkii.KpC + 2.0 * k2.camkii.KpC + 2.0 * k3.camkii.KpC + k4.camkii.KpC);
    s_next.camkii.Kp += dt6 * (k1.camkii.Kp + 2.0 * k2.camkii.Kp + 2.0 * k3.camkii.Kp + k4.camkii.Kp);

    s_next.PP1 += dt6 * (k1.PP1 + 2.0 * k2.PP1 + 2.0 * k3.PP1 + k4.PP1);
    s_next.I1 += dt6 * (k1.I1 + 2.0 * k2.I1 + 2.0 * k3.I1 + k4.I1);
    s_next.I1p += dt6 * (k1.I1p + 2.0 * k2.I1p + 2.0 * k3.I1p + k4.I1p);
    s_next.PP1I1p += dt6 * (k1.PP1I1p + 2.0 * k2.PP1I1p + 2.0 * k3.PP1I1p + k4.PP1I1p);
    s_next.PP1KpC += dt6 * (k1.PP1KpC + 2.0 * k2.PP1KpC + 2.0 * k3.PP1KpC + k4.PP1KpC);

    return s_next;
}

int main() {
    std::cout << "=== BANC D'ESSAI OMEGA V1 : BALAYAGE CALCIUM (RK4) ===\n";
    std::cout << "Génération du fichier 'courbe_calcium.csv'...\n";

    std::ofstream csv("courbe_calcium.csv");
    csv << "Calcium_uM,CaMKII_Active_Kp,PP1_Libre,Invariant_CaM\n";

    // Balayage du calcium de 0.01 à 10.0 µM
    for (double ca = 0.01; ca <= 10.0; ca *= 1.5) {
        
        OmegaFullParameters p{};
        // Paramètres Bhalla/Hayer
        p.calmodulin.kon1 = 5.0; p.calmodulin.koff1 = 50.0;
        p.calmodulin.kon2 = 100.0; p.calmodulin.koff2 = 2000.0;
        p.camkii.kon_bind = 1.27e-3; p.camkii.koff_bind = 0.07; p.camkii.k_phosphorylation = 0.5;
        p.dephosphorylation.kon = 1.0; p.dephosphorylation.koff = 0.001; p.dephosphorylation.kcat = 0.41;
        p.pp1_i1.kon = 1.0; p.pp1_i1.koff = 0.001;

        p.calcium = ca;
        p.pka_activity = 0.05;

        // État initial
        OmegaFullState s{};
        s.cam.M0 = 10.0; s.camkii.K0 = 50.0; s.PP1 = 2.0; s.I1 = 1.0;

        // Simulation pour atteindre l'équilibre
        double dt = 0.01; 
        int steps = 2000;

        for (int i = 0; i < steps; ++i) {
            s = rk4_step(s, p, dt); // VRAI RK4
        }

        // Écriture des résultats dans le CSV
        double inv = total_omega_cam(s);
        csv << ca << "," << s.camkii.Kp << "," << s.PP1 << "," << inv << "\n";
        
        std::cout << "Calcium: " << std::setw(6) << ca 
                  << " uM | CaMKII: " << std::setw(8) << s.camkii.Kp 
                  << " | PP1: " << std::setw(8) << s.PP1 
                  << " | Invariant: " << inv << "\n";
    }

    csv.close();
    std::cout << "\nTermine. Ouvre 'courbe_calcium.csv' pour voir la courbe.\n";
    return 0;
}
