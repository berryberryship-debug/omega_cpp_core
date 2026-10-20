#include "biophysical_synapse.hpp"
#include <cmath>
#include <algorithm>

BiophysicalSynapse::BiophysicalSynapse(double init_W, double init_K, double init_P, double homeo_f, unsigned int seed_val)
    : W(init_W), K(init_K), P(init_P), homeo_force(homeo_f), gen(seed_val), dist(0.0, 1.0) {}

void BiophysicalSynapse::update(double Ca, double post_activity, double dt) {
    double Ca_norm_K = std::pow(Ca, 3) / (std::pow(2.0, 3) + std::pow(Ca, 3)); 
    double act_K = Ca_norm_K; 
    double act_P = 1.5 / (1.5 + Ca); 

    double auto_K = k_auto * (std::pow(K, n_hill) / (std::pow(k_half, n_hill) + std::pow(K, n_hill)));
    double auto_P = k_auto * (std::pow(P, n_hill) / (std::pow(k_half, n_hill) + std::pow(P, n_hill)));

    double act_I_Ca = k_I * (std::pow(Ca, 2) / (std::pow(1.0, 2) + std::pow(Ca, 2)));
    double act_I_K  = k_I * (std::pow(K, 2) / (std::pow(0.3, 2) + std::pow(K, 2)));
    double target_I = std::max(act_I_Ca, act_I_K);

    if (target_I > I) {
        I += (target_I - I) * 10.0 * dt; 
    } else {
        I -= I * (1.0 / 5.0) * dt; 
    }
    I = std::clamp(I, 0.0, 1.0);

    double inh_P_on_K = 2.0 * P;
    double inh_I_on_P = 5.0 * I; 

    double dK_dt = (act_K + auto_K - inh_P_on_K - K) * tau_k;
    double dP_dt = (act_P + auto_P - inh_I_on_P - P) * tau_p;
    
    // L'homéostasie utilise maintenant homeo_force
    if (post_activity > 0.5 && Ca < 0.2) { 
        K -= homeo_force * dt; 
        K = std::clamp(K, 0.0, 1.0);
    }

    double noise_K = sigma_noise * std::sqrt(K + 0.1) * dist(gen) * std::sqrt(dt);
    double noise_P = sigma_noise * std::sqrt(P + 0.1) * dist(gen) * std::sqrt(dt);

    K = std::clamp(K + dK_dt * dt + noise_K, 0.0, 1.0);
    P = std::clamp(P + dP_dt * dt + noise_P, 0.0, 1.0);

    double drive_K = eta_k * K * (w_max - W);
    double drive_P = eta_p * P * (W - w_min);
    double dW_dt = drive_K - drive_P;
    double noise_W = sigma_noise * std::sqrt(K + P + 0.1) * dist(gen) * std::sqrt(dt);

    W = std::clamp(W + dW_dt * dt + noise_W, w_min, w_max);
}

double BiophysicalSynapse::get_weight() const { return W; }
double BiophysicalSynapse::get_K() const { return K; }
double BiophysicalSynapse::get_P() const { return P; }
double BiophysicalSynapse::get_I() const { return I; }
