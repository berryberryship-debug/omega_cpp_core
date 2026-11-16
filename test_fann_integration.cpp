
#include <iostream>
#include <iomanip>
#include <fann.h>
#include "biophysical_fann_adapter.hpp"

int main() {
    struct fann *ann = fann_create_standard(2, 1, 1);

    if (ann == nullptr) {
        std::cerr << "Erreur : création du réseau FANN impossible.\n";
        return 1;
    }

    // Réseau : entrée 0, biais 1, sortie 2
    fann_set_weight(ann, 0, 2, 1.0);
    fann_set_weight(ann, 1, 2, 0.0);

    BiophysicalFannAdapter synapse(1.0, 0.5, 42);

    const double dt = 0.001;
    const int steps = 20000;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Simulation FANN + BiophysicalSynapse\n";
    std::cout << "t (s) | Stimulus | Post | Ca    | W (FANN)\n";
    std::cout << "--------------------------------------------------\n";

    for (int i = 0; i < steps; ++i) {
        double t = i * dt;

        fann_type input[1] = {0.0};

        if (t >= 2.0 && t < 2.1) {
            input[0] = 1.0;
        }

        fann_type *output = fann_run(ann, input);

        double pre_activity = input[0];
        double post_activity =
            (output[0] > 0.5) ? 1.0 : 0.0;

        synapse.update_synapse(pre_activity, post_activity, dt);

        // Réinjection du poids biophysique dans FANN
        fann_set_weight(ann, 0, 2, synapse.get_weight());

        if ((t >= 1.9 && t <= 2.5 && i % 50 == 0)
            || (i % 2000 == 0)) {

            std::cout << std::setw(5) << t << " | "
                      << std::setw(8) << pre_activity << " | "
                      << std::setw(4) << post_activity << " | "
                      << std::setw(5) << synapse.get_Ca() << " | "
                      << std::setw(5) << synapse.get_weight()
                      << "\n";
        }
    }

    fann_destroy(ann);
    return 0;
}
