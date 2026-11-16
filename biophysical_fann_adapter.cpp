cat > biophysical_fann_adapter.cpp <<'EOF'
#include "biophysical_fann_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

BiophysicalFannAdapter::BiophysicalFannAdapter(
    double initial_weight,
    double homeo_force,
    unsigned int seed,
    bool use_noise
)
    : synapse(
        initial_weight,
        0.3,
        0.3,
        homeo_force,
        seed,
        use_noise
      ),
      current_Ca(Ca_basal)
{}

void BiophysicalFannAdapter::update_synapse(
    double pre_activity,
    double post_activity,
    double dt
) {
    if (!std::isfinite(pre_activity) ||
        !std::isfinite(post_activity) ||
        !std::isfinite(dt) ||
        dt <= 0.0) {
        throw std::invalid_argument("Invalid adapter inputs");
    }

    pre_activity = std::clamp(pre_activity, 0.0, 1.0);
    post_activity = std::clamp(post_activity, 0.0, 1.0);

    const double coincidence = pre_activity * post_activity;
    const double target_Ca = Ca_basal + alpha_Ca * coincidence;
    const double rate = std::min(1.0, dt / tau_Ca);

    current_Ca += (target_Ca - current_Ca) * rate;
    current_Ca = std::max(current_Ca, Ca_basal);

    synapse.update(current_Ca, post_activity, dt);
}

double BiophysicalFannAdapter::get_weight() const {
    return synapse.get_weight();
}

double BiophysicalFannAdapter::get_Ca() const {
    return current_Ca;
}

double BiophysicalFannAdapter::get_K() const {
    return synapse.get_K();
}

double BiophysicalFannAdapter::get_P() const {
    return synapse.get_P();
}

double BiophysicalFannAdapter::get_I() const {
    return synapse.get_I();
}
EOF
