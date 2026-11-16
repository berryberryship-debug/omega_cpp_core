cat > audit_stationnarite.cpp <<'EOF'
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

    for (double x : values)
        sum += x;

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
    constexpr int steps_per_second = 1000;
    constexpr int duration = 10;

    const std::vector<double> amplitudes = {
        0.00, 0.01, 0.02, 0.05
    };

    std::cout << std::fixed << std::setprecision(8);

    std::cout
        << "Sigma  Temps  Variable  Moyenne       Ecart-type\n"
        << "---------------------------------------------------\n";

    for (double sigma : amplitudes) {
        std::vector<double> W, K, P, I;

        W.reserve(N);
        K.reserve(N);
        P.reserve(N);
        I.reserve(N);

        std::vector<BiophysicalFannAdapter> models;
        models.reserve(N);

        for (int seed = 1; seed <= N; ++seed) {
            models.emplace_back(
                0.6,
                0.0,
                static_cast<unsigned int>(seed),
                true,
                sigma
            );
        }

        for (int second = 1; second <= duration; ++second) {
            W.clear();
            K.clear();
            P.clear();
            I.clear();

            for (auto& model : models) {
                for (int step = 0; step < steps_per_second; ++step) {
                    model.update_synapse(0.0, 0.0, dt);
                }

                W.push_back(model.get_weight());
                K.push_back(model.get_K());
                P.push_back(model.get_P());
                I.push_back(model.get_I());
            }

            const Statistics sw = calculate(W);
            const Statistics sk = calculate(K);
            const Statistics sp = calculate(P);
            const Statistics si = calculate(I);

            auto print = [sigma, second](
                const char* name,
                const Statistics& s
            ) {
                std::cout
                    << std::setw(5) << sigma << "  "
                    << std::setw(5) << second << "  "
                    << std::setw(8) << name << "  "
                    << std::setw(12) << s.mean << "  "
                    << std::setw(12) << s.sd << '\n';
            };

            print("W", sw);
            print("K", sk);
            print("P", sp);
            print("I", si);

            std::cout << "---------------------------------------------------\n";
        }
    }

    std::cout << "Graines : " << N << '\n';
    std::cout << "Duree : " << duration << " s\n";
    std::cout << "Pas : " << dt << " s\n";

    return 0;
}
EOF
