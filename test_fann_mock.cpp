
#include <iostream>
#include <iomanip>
#include "biophysical_fann_adapter.hpp"

// Simulation minimaliste de FANN
struct MockFann {
    double weight;
};

MockFann* fann_create_standard(int num_layers, int input, int output) {
    return new MockFann{1.0};
}

double* fann_run(MockFann* ann, double* input) {
    static double output[1];
    output[0] = input[0] * ann->weight;
    return output;
}

void fann_set_weight(MockFann* ann, int from, int to, double w) {
    ann->weight = w;
}

void fann_destroy(MockFann* ann) {
    delete ann;
}

int main() {
    MockFann* ann = fann_create_standard(2, 1, 1);
    BiophysicalFannAdapter synapse(1.0, 0.5, 42);

    const double dt = 0.001;
    const int steps = 20000;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Simulation Mock-FANN + BiophysicalSynapse\n";
    std::cout << "t (s) | Stim | Post | Ca    | W (FANN)\n";
    std::cout << "--------------------------------------------------\n";

    for (int i = 0; i < steps; ++i) {
        double t = i * dt;

        double input[1] = {0.0};

        if (t >= 2.0 && t < 2.1) {
            input[0] = 1.0;
        }

        double* output = fann_run(ann, input);

        double pre_activity = input[0];
        double post_activity =
            (output[0] > 0.5) ? 1.0 : 0.0;

        synapse.update_synapse(pre_activity, post_activity, dt);

        fann_set_weight(ann, 1, 2, synapse.get_weight());

        if ((t >= 1.9 && t <= 2.5 && i % 50 == 0)
            || (i % 2000 == 0)) {

            std::cout << std::setw(5) << t << " | "
                      << std::setw(4) << pre_activity << " | "
                      << std::setw(4) << post_activity << " | "
                      << std::setw(5) << synapse.get_Ca() << " | "
                      << std::setw(5) << ann->weight << "\n";
        }
    }

    fann_destroy(ann);
    return 0;
}
