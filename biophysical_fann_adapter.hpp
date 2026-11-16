#ifndef BIOPHYSICAL_FANN_ADAPTER_HPP
#define BIOPHYSICAL_FANN_ADAPTER_HPP

#include "biophysical_synapse.hpp"

class BiophysicalFannAdapter {
private:
    BiophysicalSynapse synapse;

    double current_Ca;

    // Dynamique calcique de l'interface FANN → synapse
    static constexpr double Ca_basal = 0.10;
    static constexpr double alpha_Ca = 1.00;
    static constexpr double tau_Ca = 0.050;

public:
    BiophysicalFannAdapter(
        double initial_weight,
        double homeo_force,
        unsigned int seed
    );

    void update_synapse(
        double pre_activity,
        double post_activity,
        double dt
    );

    double get_weight() const;
    double get_Ca() const;
    double get_K() const;
    double get_P() const;
};

#endif // BIOPHYSICAL_FANN_ADAPTER_HPP
