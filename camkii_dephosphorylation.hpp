
from pathlib import Path

content = r'''#ifndef OMEGA_V1_CAMKII_DEPHOSPHORYLATION_HPP
#define OMEGA_V1_CAMKII_DEPHOSPHORYLATION_HPP

namespace omega::v1 {

struct DephosphorylationState {
    double PP1;
    double KpC;
    double KC;
    double PP1KpC;
};

struct DephosphorylationParameters {
    double kon;
    double koff;
    double kcat;
};

struct DephosphorylationRates {
    double association;
    double dissociation;
    double catalysis;
};

struct DephosphorylationDerivative {
    double PP1;
    double KpC;
    double KC;
    double PP1KpC;
};

DephosphorylationRates calculate_dephosphorylation_rates(
    const DephosphorylationState& state,
    const DephosphorylationParameters& parameters
);

DephosphorylationDerivative calculate_dephosphorylation_derivative(
    const DephosphorylationState& state,
    const DephosphorylationParameters& parameters
);

double total_dephosphorylation_pp1(
    const DephosphorylationState& state
);

double total_dephosphorylation_camkii(
    const DephosphorylationState& state
);

void validate_dephosphorylation(
    const DephosphorylationState& state,
    const DephosphorylationParameters& parameters
);

} // namespace omega::v1

#endif
'''

Path("camkii_dephosphorylation.hpp").write_text(content, encoding="utf-8")
