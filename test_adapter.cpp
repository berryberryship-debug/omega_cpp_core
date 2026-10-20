#include <iostream>
#include <iomanip>
#include "biophysical_fann_adapter.hpp"

int main() {
    BiophysicalFannAdapter synapse(1.0);

    const double dt = 0.001;       // 1 ms
    const int steps = 20000;       // 20 secondes

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "t (s) | Ca    | K     | P     | I     | W\n";
    std::cout << "----------------------------------------------------\n";

    for (int i = 0; i < steps; ++i) {
        double t = i * dt;
        double pre = 0.0;
        double post = 0.0;

        // Impulsion de 100 ms à t=2.0s (Haute Fréquence)
        if (t >= 2.0 && t < 2.1) {
            pre = 1.0;
            post = 1.0;
        }

        synapse.update_synapse(pre, post, dt);

        // Affichage autour du pic et ensuite toutes les 2 secondes
        if ((t >= 1.9 && t <= 2.5 && i % 50 == 0) || (i % 2000 == 0)) {
            std::cout << std::setw(5) << t << " | " 
                      << std::setw(5) << synapse.get_Ca() << " | "
                      << std::setw(5) << synapse.get_K() << " | "
                      << std::setw(5) << synapse.get_P() << " | "
                      // On affiche pas I pour l'instant car pas d'accesseur, 
                      // mais on peut le déduire de P
                      << std::setw(5) << synapse.get_weight() << "\n";
        }
    }

    return 0;
}
