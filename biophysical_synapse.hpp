#ifndef BIOPHYSICAL_SYNAPSE_HPP
#define BIOPHYSICAL_SYNAPSE_HPP

#include <random>
#include <algorithm>
#include <cmath>

class BiophysicalSynapse {
private:
    // =========================
    // États dynamiques
    // =========================
    double W;              // poids synaptique
    double K;              // activité CaMKII
    double P;              // activité PP1
    double I = 0.0;        // état de I-1

    // =========================
    // Paramètres homéostatiques
    // =========================
    double homeo_force;

    // =========================
    // Générateur stochastique
    // =========================
    std::mt19937 gen;

    // Bruit centré sur zéro.
    // Important : un bruit uniforme [0,1] introduirait
    // artificiellement une dérive positive.
    std::normal_distribution<double> dist;

    // =========================
    // Cinétique CaMKII / PP1
    // =========================
    static constexpr double k_auto = 0.30;
    static constexpr double n_hill = 2.0;
    static constexpr double k_half = 0.50;

    static constexpr double k_I = 1.00;

    // Coefficients cinétiques utilisés actuellement
    // multiplicativement dans biophysical_synapse.cpp.
    static constexpr double tau_k = 1.00;
    static constexpr double tau_p = 1.00;

    // =========================
    // Bruit
    // =========================
    static constexpr double sigma_noise = 0.01;

    // =========================
    // Évolution du poids
    // =========================
    static constexpr double eta_k = 0.10;
    static constexpr double eta_p = 0.10;

    static constexpr double w_min = 0.0;
    static constexpr double w_max = 1.0;

public:

    BiophysicalSynapse(
        double init_W,
        double init_K,
        double init_P,
        double homeo_f,
        unsigned int seed_val
    );

    void update(
        double Ca,
        double post_activity,
        double dt
    );

    double get_weight() const;
    double get_K() const;
    double get_P() const;
    double get_I() const;
};

#endif // BIOPHYSICAL_SYNAPSE_HPP
