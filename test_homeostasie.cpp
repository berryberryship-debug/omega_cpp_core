#include <iostream>
#include <iomanip>
#include <cmath>
#include "biophysical_fann_adapter.hpp"

int main() {
    BiophysicalFannAdapter syn_A(1.0);
    BiophysicalFannAdapter syn_B(1.0);

    const double dt = 0.001;
    const int steps = 30000; // 30 secondes

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Test de l'Oubli Catastrophique (avec Homeostasie)\n";
    std::cout << "t (s) | Phase      | Stim_A | Stim_B | Post | W_A   | W_B\n";
    std::cout << "---------------------------------------------------------------\n";

    for (int i = 0; i < steps; ++i) {
        double t = i * dt;
        
        double stim_A = 0.0;
        double stim_B = 0.0;
        std::string phase = "Repos";

        if (t >= 2.0 && t < 2.2) {
            stim_A = 1.5; phase = "Train A";
        } else if (t >= 8.0 && t < 8.2) {
            stim_B = 1.5; phase = "Train B";
        } else if (t >= 15.0 && t < 18.0) {
            stim_B = 1.5; phase = "Interference"; // B force le neurone à décharger
        } else if (t >= 22.0 && t < 22.1) {
            stim_A = 0.5; phase = "Test A";
        }

        double total_input = (stim_A * syn_A.get_weight()) + (stim_B * syn_B.get_weight());
        double post_activity = (total_input > 0.8) ? 1.0 : 0.0; 
        
        syn_A.update_synapse(stim_A, post_activity, dt);
        syn_B.update_synapse(stim_B, post_activity, dt);

        if ((t >= 1.9 && t <= 2.3 && i % 50 == 0) || 
            (t >= 7.9 && t <= 8.3 && i % 50 == 0) || 
            (t >= 14.9 && t <= 18.1 && i % 500 == 0) ||
            (t >= 21.9 && t <= 22.2 && i % 50 == 0) ||
            (i % 5000 == 0)) {
            std::cout << std::setw(5) << t << " | "
                      << std::setw(10) << phase << " | "
                      << std::setw(6) << stim_A << " | "
                      << std::setw(6) << stim_B << " | "
                      << std::setw(4) << post_activity << " | "
                      << std::setw(5) << syn_A.get_weight() << " | "
                      << std::setw(5) << syn_B.get_weight() << "\n";
        }
    }

    return 0;
}
