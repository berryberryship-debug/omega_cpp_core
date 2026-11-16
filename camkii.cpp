#include "calmodulin.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace omega::v1 {

namespace {

void require_finite_nonnegative(
    double value,
    const char* name
) {
    if (!std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument(name);
    }
}

constexpr double integration_safety_factor = 0.90;
constexpr double maximum_substeps = 10'000'000.0;

} // namespace

void validate_calmodulin(
    const CaMState& state,
    double calcium,
    const CaMParameters& parameters
) {
    require_finite_nonnegative(state.M0, "Invalid M0");
    require_finite_nonnegative(state.M2N, "Invalid M2N");
    require_finite_nonnegative(state.M2C, "Invalid M2C");
    require_finite_nonnegative(state.M4, "Invalid M4");

    require_finite_nonnegative(calcium, "Invalid calcium");

    require_finite_nonnegative(parameters.kon_0N, "Invalid kon_0N");
    require_finite_nonnegative(parameters.koff_0N, "Invalid koff_0N");
    require_finite_nonnegative(parameters.kon_0C, "Invalid kon_0C");
    require_finite_nonnegative(parameters.koff_0C, "Invalid koff_0C");
    require_finite_nonnegative(parameters.kon_N4, "Invalid kon_N4");
    require_finite_nonnegative(parameters.koff_N4, "Invalid koff_N4");
    require_finite_nonnegative(parameters.kon_C4, "Invalid kon_C4");
    require_finite_nonnegative(parameters.koff_C4, "Invalid koff_C4");
}

CaMFluxes calculate_fluxes(
    const CaMState& state,
    double calcium,
    const CaMParameters& parameters
) {
    validate_calmodulin(state, calcium, parameters);

    const double calcium_squared = calcium * calcium;

    return {
        parameters.kon_0N * calcium_squared * state.M0,
        parameters.koff_0N * state.M2N,

        parameters.kon_0C * calcium_squared * state.M0,
        parameters.koff_0C * state.M2C,

        parameters.kon_N4 * calcium_squared * state.M2N,
        parameters.koff_N4 * state.M4,

        parameters.kon_C4 * calcium_squared * state.M2C,
        parameters.koff_C4 * state.M4
    };
}

CaMDerivative calculate_derivative(
    const CaMState& state,
    double calcium,
    const CaMParameters& parameters
) {
    const CaMFluxes fluxes =
        calculate_fluxes(state, calcium, parameters);

    return {
        -fluxes.v_0N_forward + fluxes.v_0N_reverse
        -fluxes.v_0C_forward + fluxes.v_0C_reverse,

        fluxes.v_0N_forward - fluxes.v_0N_reverse
        -fluxes.v_N4_forward + fluxes.v_N4_reverse,

        fluxes.v_0C_forward - fluxes.v_0C_reverse
        -fluxes.v_C4_forward + fluxes.v_C4_reverse,

        fluxes.v_N4_forward - fluxes.v_N4_reverse
        +fluxes.v_C4_forward - fluxes.v_C4_reverse
    };
}

double total_calmodulin(const CaMState& state) {
    require_finite_nonnegative(state.M0, "Invalid M0");
    require_finite_nonnegative(state.M2N, "Invalid M2N");
    require_finite_nonnegative(state.M2C, "Invalid M2C");
    require_finite_nonnegative(state.M4, "Invalid M4");

    return state.M0
         + state.M2N
         + state.M2C
         + state.M4;
}

void advance_calmodulin(
    CaMState& state,
    double calcium,
    const CaMParameters& parameters,
    double duration
) {
    validate_calmodulin(state, calcium, parameters);

    if (!std::isfinite(duration) || duration <= 0.0) {
        throw std::invalid_argument("Invalid integration duration");
    }

    const double calcium_squared = calcium * calcium;

    const double exit_rate_M0 =
        (parameters.kon_0N + parameters.kon_0C) * calcium_squared;

    const double exit_rate_M2N =
        parameters.koff_0N
        + parameters.kon_N4 * calcium_squared;

    const double exit_rate_M2C =
        parameters.koff_0C
        + parameters.kon_C4 * calcium_squared;

    const double exit_rate_M4 =
        parameters.koff_N4 + parameters.koff_C4;

    const double maximum_exit_rate = std::max(
        std::max(exit_rate_M0, exit_rate_M2N),
        std::max(exit_rate_M2C, exit_rate_M4)
    );

    const double required_substeps =
        std::ceil(
            duration * maximum_exit_rate / integration_safety_factor
        );

    if (!std::isfinite(required_substeps) ||
        required_substeps > maximum_substeps) {
        throw std::runtime_error("Too many integration substeps");
    }

    const std::size_t substep_count =
        static_cast<std::size_t>(std::max(1.0, required_substeps));

    const double substep_duration =
        duration / static_cast<double>(substep_count);

    CaMState next_state = state;

    for (std::size_t substep = 0;
         substep < substep_count;
         ++substep) {
        const CaMDerivative derivative =
            calculate_derivative(
                next_state,
                calcium,
                parameters
            );

        const CaMState candidate_state{
            next_state.M0 + substep_duration * derivative.M0,
            next_state.M2N + substep_duration * derivative.M2N,
            next_state.M2C + substep_duration * derivative.M2C,
            next_state.M4 + substep_duration * derivative.M4
        };

        const bool is_valid =
            std::isfinite(candidate_state.M0) &&
            std::isfinite(candidate_state.M2N) &&
            std::isfinite(candidate_state.M2C) &&
            std::isfinite(candidate_state.M4) &&
            candidate_state.M0 >= 0.0 &&
            candidate_state.M2N >= 0.0 &&
            candidate_state.M2C >= 0.0 &&
            candidate_state.M4 >= 0.0;

        if (!is_valid) {
            throw std::runtime_error(
                "Positivity failure during integration"
            );
        }

        next_state = candidate_state;
    }

    state = next_state;
}

} // namespace omega::v1
