#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fann.h>
#include "biophysical_fann_adapter.hpp"

int main() {
    struct fann *ann = fann_create_standard(2, 1, 1);
    fann_set_weight(ann, 1, 2, 1.0); // Poids initial de 1.0

    BiophysicalFannAdapter synapse(1.0);

    const double dt = 0.001; // 1 ms
    const int steps = 20000; // 20 secondes

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Simulation FANN + BiophysicalSynapse\n";
    std::cout << "t (s) | Stimulus | Post | Ca    | W (FANN)\n";
    std::cout << "--------------------------------------------------\n";

    for (int i = 0; i < steps; ++i) {
        double t = i * dt;
        
        fann_type input[1] = {0.0};
        if (t >= 2.0 && t < 2.1) { // Rafale de 100 ms
            input[0] = 1.0;
        }

        fann_type *output = fann_run(ann, input);
        
        double pre_activity = input[0];
        double post_activity = output[0];

        synapse.update_synapse(pre_activity, post_activity, dt);
        fann_set_weight(ann, 1, 2, synapse.get_weight());

        if ((t >= 1.9 && t <= 2.5 && i % 50 == 0) || (i % 2000 == 0)) {
            std::cout << std::setw(5) << t << " | "
                      << std::setw(8) << pre_activity << " | "
                      << std::setw(4) << post_activity << " | "
                      << std::setw(5) << synapse.get_Ca() << " | "
                      << std::setw(5) << synapse.get_weight() << "\n";
        }
    }

    fann_destroy(ann);
    return 0;
}
