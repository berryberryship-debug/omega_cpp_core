#pragma once
#include "biophysical_synapse.hpp"

class BiophysicalFannAdapter {
private:
    BiophysicalSynapse synapse;
    double Ca_basal = 0.1;
    double alpha_Ca = 3.0;
    double tau_Ca = 0.050;   
    double current_Ca;

public:
    explicit BiophysicalFannAdapter(double initial_weight = 1.0, double homeo_force = 0.3, unsigned int seed = 12345);
    void update_synapse(double pre_activity, double post_activity, double dt);
    double get_weight() const;
    double get_Ca() const;
    double get_K() const;
    double get_P() const;
};
