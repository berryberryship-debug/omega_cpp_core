#include "biophysical_fann_adapter.hpp"
#include <algorithm>

BiophysicalFannAdapter::BiophysicalFannAdapter(double initial_weight, double homeo_force, unsigned int seed)
    : synapse(initial_weight, 0.3, 0.3, homeo_force, seed), current_Ca(Ca_basal) {}

void BiophysicalFannAdapter::update_synapse(double pre_activity, double post_activity, double dt) {
    pre_activity = std::clamp(pre_activity, 0.0, 1.0);
    post_activity = std::clamp(post_activity, 0.0, 1.0);

    double coincidence = pre_activity * post_activity;
    double target_Ca = Ca_basal + alpha_Ca * coincidence;

    double calcium_rate = std::clamp(dt / tau_Ca, 0.0, 1.0);
    current_Ca += (target_Ca - current_Ca) * calcium_rate;
    current_Ca = std::max(current_Ca, Ca_basal);

    synapse.update(current_Ca, post_activity, dt);
}

double BiophysicalFannAdapter::get_weight() const { return synapse.get_weight(); }
double BiophysicalFannAdapter::get_Ca() const { return current_Ca; }
double BiophysicalFannAdapter::get_K() const { return synapse.get_K(); }
double BiophysicalFannAdapter::get_P() const { return synapse.get_P(); }
