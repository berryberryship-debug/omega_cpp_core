
#include <iostream>
#include <iomanip>
#include <string>
#include "biophysical_fann_adapter.hpp"

struct Result {
    double W_A_after_training;
    double W_A_before_B;
    double W_A_after_interference;
    double W_A_final;
    double W_B_final;
    bool recall_A;
};

Result run_experiment(
    bool train_B,
    double homeo_A,
    unsigned int seed_A,
    unsigned int seed_B
) {
    // Même état initial dans les quatre groupes.
    // Homéostasie de B désactivée pour isoler celle de A.
    BiophysicalFannAdapter syn_A(0.6, homeo_A, seed_A);
    BiophysicalFannAdapter syn_B(0.6, 0.0, seed_B);

    const double dt = 0.001;
    const int steps = 20000;

    Result result{0.0, 0.0, 0.0, 0.0, 0.0, false};

    for (int i = 0; i < steps; ++i) {
        double t = i * dt;

        double stim_A = 0.0;
        double stim_B = 0.0;

        // Acquisition de A
        if (t >= 2.0 && t < 2.2) {
            stim_A = 1.5;
        }

        // Acquisition de B et interférence
        if (train_B) {
            if (t >= 8.0 && t < 8.2) {
                stim_B = 1.5;
            }

            if (t >= 15.0 && t < 18.0) {
                stim_B = 1.5;
            }
        }

        // Activité postsynaptique commune
        double total_input =
            stim_A * syn_A.get_weight()
            + stim_B * syn_B.get_weight();

        double post_activity =
            (total_input > 0.8) ? 1.0 : 0.0;

        // Mesures aux frontières des phases
        if (i == 2200) {
            result.W_A_after_training = syn_A.get_weight();
        }

        if (i == 8000) {
            result.W_A_before_B = syn_A.get_weight();
        }

        if (i == 18000) {
            result.W_A_after_interference = syn_A.get_weight();
        }

        // Mise à jour des deux synapses
        syn_A.update_synapse(stim_A, post_activity, dt);
        syn_B.update_synapse(stim_B, post_activity, dt);
    }

    result.W_A_final = syn_A.get_weight();
    result.W_B_final = syn_B.get_weight();

    // Rappel virtuel, sans modifier les synapses.
    const double recall_stimulus = 1.5;

    result.recall_A =
        (recall_stimulus * result.W_A_final > 0.8);

    return result;
}

void print_result(
    const std::string& group,
    const Result& r
) {
    std::cout
        << std::setw(5) << group << " | "
        << std::setw(9) << r.W_A_after_training << " | "
        << std::setw(9) << r.W_A_before_B << " | "
        << std::setw(10) << r.W_A_after_interference << " | "
        << std::setw(8) << r.W_A_final << " | "
        << std::setw(8) << r.W_B_final << " | "
        << (r.recall_A ? "OUI" : "NON")
        << "\n";
}

int main() {
    const unsigned int seed_A = 1001;
    const unsigned int seed_B = 9001;

    std::cout << std::fixed << std::setprecision(4);

    std::cout << "\n";
    std::cout << "============================================\n";
    std::cout << "  EXPERIMENT : OUBLI CATASTROPHIQUE\n";
    std::cout << "============================================\n";

    std::cout
        << "Groupe | A apres  | Avant B  | Apres B   | "
        << "A final | B final | Rappel\n";

    std::cout
        << "------------------------------------------------"
        << "-------------------\n";

    // G1 : pas d'apprentissage B, homéostasie A active
    Result G1 = run_experiment(
        false, 0.5, seed_A, seed_B
    );

    // G2 : apprentissage B, homéostasie A active
    Result G2 = run_experiment(
        true, 0.5, seed_A, seed_B
    );

    // G3 : apprentissage B, homéostasie A désactivée
    Result G3 = run_experiment(
        true, 0.0, seed_A, seed_B
    );

    // G4 : pas d'apprentissage B, homéostasie A désactivée
    Result G4 = run_experiment(
        false, 0.0, seed_A, seed_B
    );

    print_result("G1", G1);
    print_result("G2", G2);
    print_result("G3", G3);
    print_result("G4", G4);

    std::cout << "\n";
    std::cout << "Seuil de rappel : W_A > "
              << (0.8 / 1.5) << "\n";

    std::cout << "Fin de l'experience.\n";

    return 0;
}
