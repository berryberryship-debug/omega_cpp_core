cat > audit_convergence_dt.cpp <<'EOF'
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

struct Result {
    double W, K, P, I;
};

Result simulate(double dt, int steps, unsigned int seed, bool noise) {
    BiophysicalFannAdapter model(0.6, 0.0, seed, noise);

    for (int i = 0; i < steps; ++i)
        model.update_synapse(0.0, 0.0, dt);

    return {
        model.get_weight(),
        model.get_K(),
        model.get_P(),
        model.get_I()
    };
}

int main() {
    const double duration = 1.0;
    const int seeds = 100;

    const std::vector<double> time_steps = {
        0.001,
        0.0005,
        0.0001
    };

    std::vector<Result> deterministic_results;

    std::cout << std::fixed << std::setprecision(8);

    std::cout << "========== CONVERGENCE DETERMINISTE ==========\n";
    std::cout << "dt\t\tW\t\tK\t\tP\t\tI\n";

    for (double dt : time_steps) {
        const int steps = static_cast<int>(std::round(duration / dt));

        Result r = simulate(dt, steps, 42, false);
        deterministic_results.push_back(r);

        std::cout << dt << '\t'
                  << r.W << '\t'
                  << r.K << '\t'
                  << r.P << '\t'
                  << r.I << '\n';
    }

    const Result reference = deterministic_results.back();

    std::cout << "\n========== ECART AU PLUS PETIT dt ==========\n";
    std::cout << "dt\t\tDelta W\t\tDelta K\t\tDelta P\t\tDelta I\n";

    for (std::size_t i = 0; i < time_steps.size(); ++i) {
        const Result& r = deterministic_results[i];

        std::cout << time_steps[i] << '\t'
                  << r.W - reference.W << '\t'
                  << r.K - reference.K << '\t'
                  << r.P - reference.P << '\t'
                  << r.I - reference.I << '\n';
    }

    std::cout << "\n========== ENSEMBLES STOCHASTIQUES ==========\n";
    std::cout << "dt\t\tW_moy\t\tW_sd\t\tK_moy\t\tK_sd"
              << "\tP_moy\t\tP_sd\t\tI_moy\t\tI_sd\n";

    for (double dt : time_steps) {
        const int steps = static_cast<int>(std::round(duration / dt));

        std::vector<double> W, K, P, I;

        for (int seed = 1; seed <= seeds; ++seed) {
            Result r = simulate(dt, steps, seed, true);

            W.push_back(r.W);
            K.push_back(r.K);
            P.push_back(r.P);
            I.push_back(r.I);
        }

        Stats sw = calculate_stats(W);
        Stats sk = calculate_stats(K);
        Stats sp = calculate_stats(P);
        Stats si = calculate_stats(I);

        std::cout << dt << '\t'
                  << sw.mean << '\t' << sw.sd << '\t'
                  << sk.mean << '\t' << sk.sd << '\t'
                  << sp.mean << '\t' << sp.sd << '\t'
                  << si.mean << '\t' << si.sd << '\n';
    }

    return 0;
}
EOF
