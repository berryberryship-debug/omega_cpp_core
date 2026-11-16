
#include <iostream>
#include <iomanip>
#include "biophysical_fann_adapter.hpp"

int main() {
    BiophysicalFannAdapter synapse(1.0, 0.5, 42);

    const double dt = 0.001;
    const int steps = 20000;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "t (s) | Ca    | K     | P     | I     | W\n";
    std::cout << "----------------------------------------------------\n";

    for (int i = 0; i < steps; ++i) {
        double t = i * dt;
        double pre = 0.0;
        double post = 0.0;

        if (t >= 2.0 && t < 2.1) {
            pre = 1.0;
            post = 1.0;
        }

        synapse.update_synapse(pre, post, dt);

        if ((t >= 1.9 && t <= 2.5 && i % 50 == 0)
            || (i % 2000 == 0)) {

            std::cout << std::setw(5) << t << " | "
                      << std::setw(5) << synapse.get_Ca() << " | "
                      << std::setw(5) << synapse.get_K() << " | "
                      << std::setw(5) << synapse.get_P() << " | "
                      << std::setw(5) << synapse.get_I() << " | "
                      << std::setw(5) << synapse.get_weight()
                      << "\n";
        }
    }

    return 0;
}
