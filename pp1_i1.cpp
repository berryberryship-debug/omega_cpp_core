#include "pp1_i1.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace omega::v1 {
namespace {

void check_nonnegative(
    double value,
    const std::string& name
) {
    if (!std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument(
            "Invalid " + name + ": " + std::to_string(value)
        );
    }
}

void check_finite_rate(
    double value,
    const std::string& name
) {
    if (!std::isfinite(value) || value < 0.0) {
        throw std::overflow_error(
            "Invalid calculated rate " + name + ": "
            + std::to_string(value)
        );
    }
}

} // namespace

void validate_pp1_i1(
    const PP1I1State& state,
    const PP1I1Parameters& parameters,
    double pka_activity
) {
    check_nonnegative(state.PP1, "PP1");
    check_nonnegative(state.I1, "I1");
    check_nonnegative(state.I1p, "I1p");
    check_nonnegative(state.PP1I1p, "PP1I1p");

    check_nonnegative(
        parameters.k_phosphorylation,
        "k_phosphorylation"
    );

    check_nonnegative(
        parameters.kon_inhibition,
        "kon_inhibition"
    );

    check_nonnegative(
        parameters.koff_inhibition,
        "koff_inhibition"
    );

    check_nonnegative(
        parameters.k_dephosphorylation,
        "k_dephosphorylation"
    );

    check_nonnegative(
        pka_activity,
        "PKA activity"
    );
}

PP1I1Rates calculate_pp1_i1_rates(
    const PP1I1State& state,
    const PP1I1Parameters& parameters,
    double pka_activity
) {
    validate_pp1_i1(
        state,
        parameters,
        pka_activity
    );

    const double phosphorylation =
        parameters.k_phosphorylation
        * pka_activity
        * state.I1;

    const double dephosphorylation =
        parameters.k_dephosphorylation
        * state.I1p;

    const double inhibition_forward =
        parameters.kon_inhibition
        * state.PP1
        * state.I1p;

    const double inhibition_reverse =
        parameters.koff_inhibition
        * state.PP1I1p;

    check_finite_rate(
        phosphorylation,
        "phosphorylation"
    );

    check_finite_rate(
        dephosphorylation,
        "dephosphorylation"
    );

    check_finite_rate(
        inhibition_forward,
        "inhibition_forward"
    );

    check_finite_rate(
        inhibition_reverse,
        "inhibition_reverse"
    );

    return {
        phosphorylation,
        dephosphorylation,
        inhibition_forward,
        inhibition_reverse
    };
}

PP1I1Derivative calculate_pp1_i1_derivative(
    const PP1I1State& state,
    const PP1I1Parameters& parameters,
    double pka_activity
) {
    const PP1I1Rates rates =
        calculate_pp1_i1_rates(
            state,
            parameters,
            pka_activity
        );

    const double Jphosphorylation =
        rates.phosphorylation;

    const double Jdephosphorylation =
        rates.dephosphorylation;

    const double Jnet_inhibition =
        rates.inhibition_forward
        - rates.inhibition_reverse;

    PP1I1Derivative derivative;

    derivative.PP1 =
        -Jnet_inhibition;

    derivative.I1 =
        -Jphosphorylation
        + Jdephosphorylation;

    derivative.I1p =
        Jphosphorylation
        - Jdephosphorylation
        - Jnet_inhibition;

    derivative.PP1I1p =
        Jnet_inhibition;

    return derivative;
}

double total_pp1(
    const PP1I1State& state
) {
    check_nonnegative(state.PP1, "PP1");
    check_nonnegative(state.PP1I1p, "PP1I1p");

    return state.PP1 + state.PP1I1p;
}

double total_i1(
    const PP1I1State& state
) {
    check_nonnegative(state.I1, "I1");
    check_nonnegative(state.I1p, "I1p");
    check_nonnegative(state.PP1I1p, "PP1I1p");

    return state.I1
        + state.I1p
        + state.PP1I1p;
}

} // namespace omega::v1
