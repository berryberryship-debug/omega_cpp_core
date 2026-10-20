#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "biophysical_fann_adapter.hpp"

// Simulation minimaliste de FANN pour le test
struct MockFann {
    double weight;
};

MockFann* fann_create_standard(int num_layers, int input, int output) {
    return new MockFann{1.0}; // Poids initial de 1.0
}

double* fann_run(MockFann* ann, double* input) {
    static double output[1];
    // L'activité de sortie est simplement l'entrée multipliée par le poids
    output[0] = input[0] * ann->weight;
    // FANN utilise généralement une fonction sigmoïde, mais pour ce test, 
    // une simple multiplication suffit à démontrer l'adaptation du poids.
    return output;
}

void fann_set_weight(MockFann* ann, int from, int to, double w) {
    ann->weight = w;
}

void fann_destroy(MockFann* ann) {
    delete ann;
}

int main() {
    MockFann *ann = fann_create_standard(2, 1, 1);
    BiophysicalFannAdapter synapse(1.0);

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
            input[0] = 1.0; // Rafale de 100 ms
        }

        double *output = fann_run(ann, input);
        
        double pre_activity = input[0];
        double post_activity = output[0];

        // Mise à jour biophysique
        synapse.update_synapse(pre_activity, post_activity, dt);
        
        // On réinjecte le poids dans le réseau
        fann_set_weight(ann, 1, 2, synapse.get_weight());

        if ((t >= 1.9 && t <= 2.5 && i % 50 == 0) || (i % 2000 == 0)) {
            std::cout << std::setw(5) << t << " | "
                      << std::setw(4) << pre_activity << " | "
                      << std::setw(4) << post_activity << " | "
                      << std::setw(5) << synapse.get_Ca() << " | "
                      << std::setw(5) << synapse.get_weight() << "\n";
        }
    }

    fann_destroy(ann);
    return 0;
}
