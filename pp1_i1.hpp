#ifndef OMEGA_V1_PP1_I1_HPP
#define OMEGA_V1_PP1_I1_HPP

namespace omega::v1 {

struct PP1I1State {
    double PP1;
    double I1;
    double I1p;
    double PP1I1p;
};

struct PP1I1Parameters {
    double k_phosphorylation;
    double kon_inhibition;
    double koff_inhibition;
};

struct PP1I1Rates {
    double phosphorylation;
    double inhibition_forward;
    double inhibition_reverse;
};

struct PP1I1Derivative {
    double PP1;
    double I1;
    double I1p;
    double PP1I1p;
};

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

void validate_pp1_i1(
    const PP1I1State& state,
    const PP1I1Parameters& parameters,
    double pka_activity
);

} // namespace omega::v1

#endif
