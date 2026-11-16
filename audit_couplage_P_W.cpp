cat > audit_couplage_P_W.cpp <<'EOF'
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include "biophysical_fann_adapter.hpp"

struct Stats {
    double mean = 0.0;
    double sd = 0.0;
};

Stats calculate_stats(const std::vector<double>& values) {
    Stats s;

    if (values.empty()) return s;

    for (double x : values)
        s.mean += x;

    s.mean /= static_cast<double>(values.size());

    double variance = 0.0;

    for (double x : values)
        variance += (x - s.mean) * (x - s.mean);

    variance /= static_cast<double>(values.size());
    s.sd = std::sqrt(variance);

    return s;
}

int main() {
    const double dt = 0.001;
    const int steps = 1000;
    const int seeds = 100;
    const int interval = 100;

    std::cout << std::fixed << std::setprecision(8);

    std::cout << "t\tW_det_moy\tW_stoch_moy"
              << "\tW_det_sd\tW_stoch_sd"
              << "\tDelta_moy\n";

    for (int step = 0; step <= steps; step += interval) {
        std::vector<double> W_det;
        std::vector<double> W_stoch;

        for (int seed = 1; seed <= seeds; ++seed) {
            BiophysicalFannAdapter det(0.6, 0.0, seed, false);
            BiophysicalFannAdapter stoch(0.6, 0.0, seed, true);

            for (int i = 0; i < step; ++i) {
                det.update_synapse(0.0, 0.0, dt);
                stoch.update_synapse(0.0, 0.0, dt);
            }

            W_det.push_back(det.get_weight());
            W_stoch.push_back(stoch.get_weight());
        }

        const Stats d = calculate_stats(W_det);
        const Stats s = calculate_stats(W_stoch);

        std::cout
            << step * dt << '\t'
            << d.mean << '\t'
            << s.mean << '\t'
            << d.sd << '\t'
            << s.sd << '\t'
            << s.mean - d.mean << '\n';
    }

    return 0;
}
EOF
