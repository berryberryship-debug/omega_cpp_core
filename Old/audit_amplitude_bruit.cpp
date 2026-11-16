cat > audit_amplitude_bruit.cpp <<'EOF'
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "biophysical_fann_adapter.hpp"

struct Statistics {
    double mean;
    double sd;
};

Statistics calculate(const std::vector<double>& values) {
    double sum = 0.0;

    for (double x : values) {
        sum += x;
    }

    const double mean = sum / values.size();

    double variance = 0.0;

    for (double x : values) {
        const double d = x - mean;
        variance += d * d;
    }

    variance /= values.size();

    return {mean, std::sqrt(variance)};
}

int main() {
    constexpr int N = 100;
    constexpr double dt = 0.001;
    constexpr int steps = 1000;

    const std::vector<double> amplitudes = {
        0.00, 0.01, 0.02, 0.05
    };

    std::cout << std::fixed << std::setprecision(8);

    std::cout
        << "Sigma   Variable   Moyenne       Ecart-type\n"
        << "--------------------------------------------\n";

    for (double sigma : amplitudes) {
        std::vector<double> W_values;
        std::vector<double> K_values;
        std::vector<double> P_values;
        std::vector<double> I_values;

        for (int seed = 1; seed <= N; ++seed) {
            BiophysicalFannAdapter model(
                0.6,
                0.0,
                static_cast<unsigned int>(seed),
                true,
                sigma
            );

            for (int step = 0; step < steps; ++step) {
                model.update_synapse(0.0, 0.0, dt);
            }

            W_values.push_back(model.get_weight());
            K_values.push_back(model.get_K());
            P_values.push_back(model.get_P());
            I_values.push_back(model.get_I());
        }

        const Statistics W = calculate(W_values);
        const Statistics K = calculate(K_values);
        const Statistics P = calculate(P_values);
        const Statistics I = calculate(I_values);

        auto print = [sigma](const char* name,
                             const Statistics& s) {
            std::cout
                << std::setw(5) << sigma << "   "
                << std::setw(8) << name << "   "
                << std::setw(12) << s.mean << "   "
                << std::setw(12) << s.sd << '\n';
        };

        print("W", W);
        print("K", K);
        print("P", P);
        print("I", I);

        std::cout << "--------------------------------------------\n";
    }

    std::cout << "Graines par amplitude : " << N << '\n';
    std::cout << "Duree : " << dt * steps << " s\n";

    return 0;
}
EOF
