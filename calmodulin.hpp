
from pathlib import Path

content = r'''#include "calmodulin.hpp"
#include <algorithm>

#include <cmath>
#include <stdexcept>

namespace omega::v1 {

namespace {

void require_finite_nonnegative(double value, const char* name) {
    if (!std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument(name);
    }
}

}

void validate_calmodulin(
    const CaMState& s,
    double calcium,
    const CaMParameters& p
) {
    require_finite_nonnegative(s.M0, "Invalid M0");
    require_finite_nonnegative(s.M2N, "Invalid M2N");
    require_finite_nonnegative(s.M2C, "Invalid M2C");
    require_finite_nonnegative(s.M4, "Invalid M4");

    require_finite_nonnegative(calcium, "Invalid calcium");

    require_finite_nonnegative(p.kon_0N, "Invalid kon_0N");
    require_finite_nonnegative(p.koff_0N, "Invalid koff_0N");
    require_finite_nonnegative(p.kon_0C, "Invalid kon_0C");
    require_finite_nonnegative(p.koff_0C, "Invalid koff_0C");
    require_finite_nonnegative(p.kon_N4, "Invalid kon_N4");
    require_finite_nonnegative(p.koff_N4, "Invalid koff_N4");
    require_finite_nonnegative(p.kon_C4, "Invalid kon_C4");
    require_finite_nonnegative(p.koff_C4, "Invalid koff_C4");
}

CaMFluxes calculate_fluxes(
    const CaMState& s,
    double calcium,
    const CaMParameters& p
) {
    validate_calmodulin(s, calcium, p);

    const double calcium_squared = calcium * calcium;

    return {
        p.kon_0N * calcium_squared * s.M0,
        p.koff_0N * s.M2N,

        p.kon_0C * calcium_squared * s.M0,
        p.koff_0C * s.M2C,

        p.kon_N4 * calcium_squared * s.M2N,
        p.koff_N4 * s.M4,

        p.kon_C4 * calcium_squared * s.M2C,
        p.koff_C4 * s.M4
    };
}

CaMDerivative calculate_derivative(
    const CaMState& s,
    double calcium,
    const CaMParameters& p
) {
    const CaMFluxes f = calculate_fluxes(s, calcium, p);

    return {
        -f.v_0N_forward + f.v_0N_reverse
        -f.v_0C_forward + f.v_0C_reverse,

        f.v_0N_forward - f.v_0N_reverse
        -f.v_N4_forward + f.v_N4_reverse,

        f.v_0C_forward - f.v_0C_reverse
        -f.v_C4_forward + f.v_C4_reverse,

        f.v_N4_forward - f.v_N4_reverse
        +f.v_C4_forward - f.v_C4_reverse
    };
}

double total_calmodulin(const CaMState& s) {
    require_finite_nonnegative(s.M0, "Invalid M0");
    require_finite_nonnegative(s.M2N, "Invalid M2N");
    require_finite_nonnegative(s.M2C, "Invalid M2C");
    require_finite_nonnegative(s.M4, "Invalid M4");

    return s.M0 + s.M2N + s.M2C + s.M4;
}

} // namespace omega::v1

namespace omega::v1 {

void advance_calmodulin(
    CaMState& state,
    double calcium,
    const CaMParameters& p,
    double duration
) {
    validate_calmodulin(state, calcium, p);

    if (!std::isfinite(duration) || duration <= 0.0) {
        throw std::invalid_argument("Invalid integration duration");
    }

    const double c2 = calcium * calcium;

    const double exit_M0 = (p.kon_0N + p.kon_0C) * c2;
    const double exit_M2N = p.koff_0N + p.kon_N4 * c2;
    const double exit_M2C = p.koff_0C + p.kon_C4 * c2;
    const double exit_M4 = p.koff_N4 + p.koff_C4;

    const double max_exit = std::max(
        std::max(exit_M0, exit_M2N),
        std::max(exit_M2C, exit_M4)
    );

    constexpr double safety = 0.90;
    constexpr double max_substeps = 10000000.0;

    const double required =
        std::ceil(duration * max_exit / safety);

    if (!std::isfinite(required) || required > max_substeps) {
        throw std::runtime_error("Too many integration substeps");
    }

    const std::size_t steps =
        static_cast<std::size_t>(std::max(1.0, required));

    const double h = duration / static_cast<double>(steps);

    CaMState next = state;

    for (std::size_t i = 0; i < steps; ++i) {
        const CaMDerivative d =
            calculate_derivative(next, calcium, p);

        CaMState candidate{
            next.M0   + h * d.M0,
            next.M2N  + h * d.M2N,
            next.M2C  + h * d.M2C,
            next.M4   + h * d.M4
        };

        if (!std::isfinite(candidate.M0) ||
            !std::isfinite(candidate.M2N) ||
            !std::isfinite(candidate.M2C) ||
            !std::isfinite(candidate.M4) ||
            candidate.M0 < 0.0 ||
            candidate.M2N < 0.0 ||
            candidate.M2C < 0.0 ||
            candidate.M4 < 0.0) {
            throw std::runtime_error(
                "Positivity failure during integration"
            );
        }

        next = candidate;
    }

    state = next;
}

} // namespace omega::v1
'''

Path("calmodulin.cpp").write_text(content, encoding="utf-8")
