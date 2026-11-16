#pragma once

namespace omega::v1 {

// Ordre des membres important : omega_full.cpp initialise par accolades
// dans l'ordre PP1, I1, I1p, PP1I1p.
struct PP1I1State {
    double PP1 = 0.0;
    double I1 = 0.0;
    double I1p = 0.0;
    double PP1I1p = 0.0;
};

struct PP1I1Parameters {
    double k_phosphorylation = 0.0;
    double kon_inhibition = 0.0;
    double koff_inhibition = 0.0;
};

struct PP1I1Rates {
    double phosphorylation = 0.0;
    double inhibition_forward = 0.0;
    double inhibition_reverse = 0.0;
};

struct PP1I1Derivative {
    double PP1 = 0.0;
    double I1 = 0.0;
    double I1p = 0.0;
    double PP1I1p = 0.0;
};

void validate_pp1_i1(
    const PP1I1State& state,
    const PP1I1Parameters& parameters,
    double pka_activity
);

PP1I1Rates calculate_pp1_i1_rates(
    const PP1I1State& state,
    const PP1I1Parameters& parameters,
    double pka_activity
);

PP1I1Derivative calculate_pp1_i1_derivative(
    const PP1I1State& state,
    const PP1I1Parameters& parameters,
    double pka_activity
);

double total_pp1(const PP1I1State& state);
double total_i1(const PP1I1State& state);

} // namespace omega::v1
