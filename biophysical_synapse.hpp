cat > biophysical_synapse.hpp <<'EOF'
#ifndef BIOPHYSICAL_SYNAPSE_HPP
#define BIOPHYSICAL_SYNAPSE_HPP

#include <random>

class BiophysicalSynapse {
private:
    double W;
    double K;
    double P;
    double I;
    double homeo_force;
    bool use_noise;

    std::mt19937 gen;
    std::normal_distribution<double> dist;

    static constexpr double k_auto = 0.30;
    static constexpr double n_hill = 2.0;
    static constexpr double k_half = 0.50;
    static constexpr double k_I = 1.00;

    static constexpr double tau_k = 1.00;
    static constexpr double tau_p = 1.00;

    static constexpr double sigma_noise = 0.01;
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
        unsigned int seed_val,
        bool noise_enabled = true
    );

    void update(double Ca, double post_activity, double dt);

    double get_weight() const;
    double get_K() const;
    double get_P() const;
    double get_I() const;
};

#endif
EOF
