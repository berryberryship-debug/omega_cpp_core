#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "biophysical_fann_adapter.hpp"

bool run_trial(double c_homeo, int trial_id) {
    // LA CORRECTION : Une graine unique pour chaque essai
    unsigned int seed = 1000 + trial_id; 
    
    BiophysicalFannAdapter syn_A(1.0, c_homeo, seed);
    BiophysicalFannAdapter syn_B(1.0, c_homeo, seed);
    
    const double dt = 0.001;
    const int steps = 25000; 
    
    for (int i = 0; i < steps; ++i) {
        double t = i * dt;
        double stim_A = 0.0, stim_B = 0.0;
        
        if (t >= 2.0 && t < 2.2) stim_A = 1.5; // Train A
        else if (t >= 8.0 && t < 8.2) stim_B = 1.5; // Train B
        else if (t >= 15.0 && t < 20.0) stim_B = 1.5; // INTERFERENCE
        
        double total_input = (stim_A * syn_A.get_weight()) + (stim_B * syn_B.get_weight());
        double post = (total_input > 0.8) ? 1.0 : 0.0;
        
        syn_A.update_synapse(stim_A, post, dt);
        syn_B.update_synapse(stim_B, post, dt);
    }
    
    // Test final de survie
    double stim_A = 0.5;
    double total_input = (stim_A * syn_A.get_weight());
    double post = (total_input > 0.8) ? 1.0 : 0.0;
    
    return (post > 0.5);
}

int main() {
    std::cout << "Test de Rupture (avec graines independantes)\n";
    std::cout << "Force  | Taux de survie (sur 20 essais)\n";
    std::cout << "------------------------------------------------------\n";
    
    std::vector<double> c_values = {0.1, 0.5, 1.0, 2.0, 4.0, 8.0};
    
    for (double c : c_values) {
        int survie = 0;
        for (int trial = 0; trial < 20; ++trial) {
            if (run_trial(c, trial)) survie++; // On passe le numéro de l'essai
        }
        std::cout << "  " << std::setw(3) << c << "  | " 
                  << std::setw(2) << survie << "/20 (" 
                  << (survie * 5) << "%)\n";
    }
    
    return 0;
}
