#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "biophysical_fann_adapter.hpp"

bool run_trial(double c_homeo) {
    BiophysicalFannAdapter syn_A(1.0, c_homeo);
    BiophysicalFannAdapter syn_B(1.0, c_homeo);
    
    const double dt = 0.001;
    const int steps = 25000; 
    
    for (int i = 0; i < steps; ++i) {
        double t = i * dt;
        double stim_A = 0.0, stim_B = 0.0;
        
        if (t >= 2.0 && t < 2.2) stim_A = 1.5; // Train A
        else if (t >= 8.0 && t < 8.2) stim_B = 1.5; // Train B
        else if (t >= 15.0 && t < 20.0) stim_B = 1.5; // INTERFERENCE PROLONGEE (5 secondes)
        
        double total_input = (stim_A * syn_A.get_weight()) + (stim_B * syn_B.get_weight());
        double post = (total_input > 0.8) ? 1.0 : 0.0;
        
        syn_A.update_synapse(stim_A, post, dt);
        syn_B.update_synapse(stim_B, post, dt);
    }
    
    // Test de survie à t=22s
    double stim_A = 0.5;
    double total_input = (stim_A * syn_A.get_weight());
    double post = (total_input > 0.8) ? 1.0 : 0.0;
    
    return (post > 0.5);
}

int main() {
    std::cout << "Test de Rupture du Verrou (Balayage Homéostatique)\n";
    std::cout << "Force  | Taux de survie de la mémoire (sur 20 essais)\n";
    std::cout << "------------------------------------------------------\n";
    
    std::vector<double> c_values = {0.1, 0.5, 1.0, 2.0, 4.0, 8.0};
    
    for (double c : c_values) {
        int survie = 0;
        for (int trial = 0; trial < 20; ++trial) {
            if (run_trial(c)) survie++;
        }
        std::cout << "  " << std::setw(3) << c << "  | " 
                  << std::setw(2) << survie << "/20 (" 
                  << (survie * 5) << "%)\n";
    }
    
    return 0;
}
