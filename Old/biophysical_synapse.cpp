cat > biophysical_synapse.cpp <<'EOF'
#include "biophysical_synapse.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

BiophysicalSynapse::BiophysicalSynapse(
    double init_W,
    double init_K,
    double init_P,
    double homeo_f,
    unsigned int seed_val,
    bool noise_enabled
)
    : W(init_W),
      K(init_K),
      P(init_P),
      I(0.0),
      homeo_force(homeo_f),
      use_noise(noise_enabled),
      gen(seed_val),
      dist(0.0, 1.0)
{
    if (!std::isfinite(init_W) ||
        !std::isfinite(init_K) ||
        !std::isfinite(init_P) ||
        !std::isfinite(homeo_f) ||
        init_W < 0.0 || init_W > 1.0 ||
        init_K < 0.0 || init_K > 1.0 ||
        init_P < 0.0 || init_P > 1.0 ||
        homeo_f < 0.0) {
        throw std::invalid_argument("Invalid synapse parameters");
    }
}

void BiophysicalSynapse::update(
    double Ca,
    double post_activity,
    double dt
) {
    if (!std::isfinite(Ca) ||
        !std::isfinite(post_activity) ||
        !std::isfinite(dt) ||
        Ca < 0.0 || dt <= 0.0) {
        throw std::invalid_argument("Invalid update inputs");
    }

    post_activity = std::clamp(post_activity, 0.0, 1.0);

    const double ca3 = std::pow(Ca, 3);
    const double act_K = ca3 / (8.0 + ca3);

    const double act_P = 1.5 / (1.5 + Ca);

    const double k2 = K * K;
    const double p2 = P * P;

    const double auto_K =
        k_auto * k2 / (k_half * k_half + k2);

    const double auto_P =
        k_auto * p2 / (k_half * k_half + p2);

    const double ca2 = Ca * Ca;
    const double act_I_Ca = ca2 / (1.0 + ca2);

    const double act_I_K =
        k2 / (0.09 + k2);

    const double target_I = std::max(act_I_Ca, act_I_K);

    if (target_I > I) {
        I += (target_I - I) * std::min(1.0, 10.0 * dt);
    } else {
        I -= I * std::min(1.0, 0.2 * dt);
    }

    I = std::clamp(I, 0.0, 1.0);

    const double dK_dt =
        act_K + auto_K - 2.0 * P - K;

    const double dP_dt =
        act_P + auto_P - 5.0 * I - P;

    double noise_K = 0.0;
    double noise_P = 0.0;

    if (use_noise) {
        noise_K =
            sigma_noise * std::sqrt(K + 0.1) *
            dist(gen) * std::sqrt(dt);

        noise_P =
            sigma_noise * std::sqrt(P + 0.1) *
            dist(gen) * std::sqrt(dt);
    }

    double next_K = K + dK_dt * dt + noise_K;
    double next_P = P + dP_dt * dt + noise_P;

    if (post_activity > 0.5 && Ca < 0.2) {
        next_K -= homeo_force * dt;
    }

    K = std::clamp(next_K, 0.0, 1.0);
    P = std::clamp(next_P, 0.0, 1.0);

    const double drive_K =
        eta_k * K * (w_max - W);

    const double drive_P =
        eta_p * P * (W - w_min);

    const double dW_dt = drive_K - drive_P;

    double noise_W = 0.0;

    if (use_noise) {
        noise_W =
            sigma_noise * std::sqrt(K + P + 0.1) *
            dist(gen) * std::sqrt(dt);
    }

    W = std::clamp(
        W + dW_dt * dt + noise_W,
        w_min,
        w_max
    );
}

double BiophysicalSynapse::get_weight() const { return W; }
double BiophysicalSynapse::get_K() const { return K; }
double BiophysicalSynapse::get_P() const { return P; }
double BiophysicalSynapse::get_I() const { return I; }
EOF
