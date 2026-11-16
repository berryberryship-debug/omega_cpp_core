cat > audit_multi_graines.cpp <<'EOF'
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

bool at_boundary(double x) {
    return x <= 1e-12 || x >= 1.0 - 1e-12;
}

int main() {
    const double dt = 0.001;
    const int steps = 1000;
    const int seeds = 100;

    std::vector<double> dW, dK, dP, dI;

    int hits_det_W = 0, hits_stoch_W = 0;
    int hits_det_K = 0, hits_stoch_K = 0;
    int hits_det_P = 0, hits_stoch_P = 0;

    for (int seed = 1; seed <= seeds; ++seed) {
        BiophysicalFannAdapter det(0.6, 0.0, seed, false);
        BiophysicalFannAdapter stoch(0.6, 0.0, seed, true);

        bool det_W = false, stoch_W = false;
        bool det_K = false, stoch_K = false;
        bool det_P = false, stoch_P = false;

        for (int i = 0; i < steps; ++i) {
            det.update_synapse(0.0, 0.0, dt);
            stoch.update_synapse(0.0, 0.0, dt);

            det_W |= at_boundary(det.get_weight());
            stoch_W |= at_boundary(stoch.get_weight());

            det_K |= at_boundary(det.get_K());
            stoch_K |= at_boundary(stoch.get_K());

            det_P |= at_boundary(det.get_P());
            stoch_P |= at_boundary(stoch.get_P());
        }

        dW.push_back(stoch.get_weight() - det.get_weight());
        dK.push_back(stoch.get_K() - det.get_K());
        dP.push_back(stoch.get_P() - det.get_P());
        dI.push_back(stoch.get_I() - det.get_I());

        hits_det_W += det_W;
        hits_stoch_W += stoch_W;

        hits_det_K += det_K;
        hits_stoch_K += stoch_K;

        hits_det_P += det_P;
        hits_stoch_P += stoch_P;
    }

    const Stats W = calculate_stats(dW);
    const Stats K = calculate_stats(dK);
    const Stats P = calculate_stats(dP);
    const Stats I = calculate_stats(dI);

    std::cout << std::fixed << std::setprecision(8);

    std::cout << "========== DIFFERENCES FINALES ==========\n";
    std::cout << "Variable\tMoyenne\t\tEcart-type\n";

    std::cout << "Delta W\t\t" << W.mean << "\t" << W.sd << '\n';
    std::cout << "Delta K\t\t" << K.mean << "\t" << K.sd << '\n';
    std::cout << "Delta P\t\t" << P.mean << "\t" << P.sd << '\n';
    std::cout << "Delta I\t\t" << I.mean << "\t" << I.sd << '\n';

    std::cout << "\n========== ATTEINTE DES BORNES ==========\n";
    std::cout << "Variable\tDeterministe\tStochastique\n";

    std::cout << "W\t\t" << hits_det_W << "/100\t\t"
              << hits_stoch_W << "/100\n";

    std::cout << "K\t\t" << hits_det_K << "/100\t\t"
              << hits_stoch_K << "/100\n";

    std::cout << "P\t\t" << hits_det_P << "/100\t\t"
              << hits_stoch_P << "/100\n";

    return 0;
}
EOF
