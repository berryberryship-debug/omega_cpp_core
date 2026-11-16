cat > audit_extinction_P.cpp <<'EOF'
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <limits>

#include "biophysical_fann_adapter.hpp"

struct Stats {
    double mean = 0.0;
    double sd = 0.0;
    double min = 0.0;
    double max = 0.0;
};

Stats calculate_stats(const std::vector<double>& values) {
    Stats s;

    if (values.empty()) return s;

    s.min = values[0];
    s.max = values[0];

    for (double x : values) {
        s.mean += x;
        if (x < s.min) s.min = x;
        if (x > s.max) s.max = x;
    }

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
    const double threshold = 1e-12;

    std::vector<double> times_det;
    std::vector<double> times_stoch;

    int not_extinct_det = 0;
    int not_extinct_stoch = 0;

    for (int seed = 1; seed <= seeds; ++seed) {
        BiophysicalFannAdapter det(0.6, 0.0, seed, false);
        BiophysicalFannAdapter stoch(0.6, 0.0, seed, true);

        double t_det = std::numeric_limits<double>::quiet_NaN();
        double t_stoch = std::numeric_limits<double>::quiet_NaN();

        for (int i = 1; i <= steps; ++i) {
            det.update_synapse(0.0, 0.0, dt);
            stoch.update_synapse(0.0, 0.0, dt);

            const double t = i * dt;

            if (std::isnan(t_det) && det.get_P() <= threshold)
                t_det = t;

            if (std::isnan(t_stoch) && stoch.get_P() <= threshold)
                t_stoch = t;
        }

        if (std::isnan(t_det))
            ++not_extinct_det;
        else
            times_det.push_back(t_det);

        if (std::isnan(t_stoch))
            ++not_extinct_stoch;
        else
            times_stoch.push_back(t_stoch);
    }

    const Stats det = calculate_stats(times_det);
    const Stats stoch = calculate_stats(times_stoch);

    std::cout << std::fixed << std::setprecision(6);

    std::cout << "========== TEMPS D'EXTINCTION DE P ==========\n";

    std::cout << "Deterministe :\n";
    std::cout << "Extinctions : " << times_det.size() << "/100\n";
    std::cout << "Non eteints : " << not_extinct_det << "/100\n";

    if (!times_det.empty()) {
        std::cout << "Moyenne : " << det.mean << " s\n";
        std::cout << "Ecart-type : " << det.sd << " s\n";
        std::cout << "Minimum : " << det.min << " s\n";
        std::cout << "Maximum : " << det.max << " s\n";
    }

    std::cout << "\nStochastique :\n";
    std::cout << "Extinctions : " << times_stoch.size() << "/100\n";
    std::cout << "Non eteints : " << not_extinct_stoch << "/100\n";

    if (!times_stoch.empty()) {
        std::cout << "Moyenne : " << stoch.mean << " s\n";
        std::cout << "Ecart-type : " << stoch.sd << " s\n";
        std::cout << "Minimum : " << stoch.min << " s\n";
        std::cout << "Maximum : " << stoch.max << " s\n";
    }

    return 0;
}
EOF
