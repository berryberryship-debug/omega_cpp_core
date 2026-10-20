#pragma once
#include <random>

class BiophysicalSynapse {
private:
    double w_min = 0.05;
    double w_max = 2.50;
    double k_d = 1.5;   
    double k_p = 2.0;   
    double eta_k = 0.08;  
    double eta_p = 0.08;  
    double tau_k = 10.0;
    double tau_p = 10.0;    
    double sigma_noise = 0.05; 

    double k_auto = 2.0; 
    double k_half = 0.25; 
    double n_hill = 4.0;

    double k_I = 2.0;    
    double K_I_half = 0.4; 
    double n_I = 4.0;

    double homeo_force; // <--- NOUVEAU

    double W, K, P, I;
    std::mt19937 gen;
    std::normal_distribution<double> dist;

public:
    BiophysicalSynapse(double init_W, double init_K, double init_P, double homeo_f, unsigned int seed_val = 12345);
    void update(double Ca, double post_activity, double dt);
    double get_weight() const;
    double get_K() const;
    double get_P() const;
    double get_I() const;
};
