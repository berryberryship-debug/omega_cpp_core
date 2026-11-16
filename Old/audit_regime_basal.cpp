
#include <iostream>
#include <iomanip>
#include <cmath>

#include "biophysical_fann_adapter.hpp"

int main() {
    const double dt = 0.001;
    const double duration = 1.0;
    const int steps = static_cast<int>(duration / dt);

    BiophysicalFannAdapter synapse(0.6, 0.0, 42);

    std::cout << std::fixed << std::setprecision(8);

    std::cout
        << "t\tCa\tK\tP\tI"
        << "\tactK\tautoK\t-perteK"
        << "\tactP\tautoP\t-inhibP\t-perteP\n";

    for (int i = 0; i <= steps; ++i) {
        const double t = i * dt;

        const double Ca = synapse.get_Ca();
        const double K = synapse.get_K();
        const double P = synapse.get_P();
        const double I = synapse.get_I();

        const double ca3 = std::pow(Ca, 3);
        const double actK = ca3 / (8.0 + ca3);

        const double k2 = K * K;
        const double autoK =
            0.30 * k2 / (0.25 + k2);

        const double actP = 1.5 / (1.5 + Ca);

        const double p2 = P * P;
        const double autoP =
            0.30 * p2 / (0.25 + p2);

        const double lossK = -K;
        const double lossP = -P;
        const double inhibP = -5.0 * I;

        if (i % 100 == 0) {
            std::cout
                << t << '\t'
                << Ca << '\t'
                << K << '\t'
                << P << '\t'
                << I << '\t'
                << actK << '\t'
                << autoK << '\t'
                << lossK << '\t'
                << actP << '\t'
                << autoP << '\t'
                << inhibP << '\t'
                << lossP << '\n';
        }

        if (i < steps) {
            synapse.update_synapse(0.0, 0.0, dt);
        }
    }

    std::cout << "\n========== ETAT FINAL ==========\n";
    std::cout << "W  : " << synapse.get_weight() << '\n';
    std::cout << "Ca : " << synapse.get_Ca() << '\n';
    std::cout << "K  : " << synapse.get_K() << '\n';
    std::cout << "P  : " << synapse.get_P() << '\n';
    std::cout << "I  : " << synapse.get_I() << '\n';

    return 0;
}
