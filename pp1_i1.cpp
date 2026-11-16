
from pathlib import Path

content = r'''#include "pp1_i1.hpp"

#include <cmath>
#include <stdexcept>

namespace omega::v1 {
namespace {

void check_nonnegative(double x, const char* name) {
    if (!std::isfinite(x) || x < 0.0)
        throw std::invalid_argument(name);
}

}

void validate_pp1_i1(
    const PP1I1State& s,
    const PP1I1Parameters& p,
    double pka
) {
    check_nonnegative(s.PP1, "Invalid PP1");
    check_nonnegative(s.I1, "Invalid I1");
    check_nonnegative(s.I1p, "Invalid I1p");
    check_nonnegative(s.PP1I1p, "Invalid PP1I1p");

    check_nonnegative(p.k_phosphorylation, "Invalid k_phosphorylation");
    check_nonnegative(p.kon_inhibition, "Invalid kon_inhibition");
    check_nonnegative(p.koff_inhibition, "Invalid koff_inhibition");
    check_nonnegative(pka, "Invalid PKA activity");
}

PP1I1Rates calculate_pp1_i1_rates(
    const PP1I1State& s,
    const PP1I1Parameters& p,
    double pka
) {
    validate_pp1_i1(s, p, pka);

    return {
        p.k_phosphorylation * pka * s.I1,
        p.kon_inhibition * s.PP1 * s.I1p,
        p.koff_inhibition * s.PP1I1p
    };
}

PP1I1Derivative calculate_pp1_i1_derivative(
    const PP1I1State& s,
    const PP1I1Parameters& p,
    double pka
) {
    const auto r = calculate_pp1_i1_rates(s, p, pka);

    const double Jphos = r.phosphorylation;
    const double Jinhib =
        r.inhibition_forward - r.inhibition_reverse;

    return {
        -Jinhib,
        -Jphos,
        Jphos - Jinhib,
        Jinhib
    };
}

double total_pp1(const PP1I1State& s) {
    check_nonnegative(s.PP1, "Invalid PP1");
    check_nonnegative(s.PP1I1p, "Invalid PP1I1p");
    return s.PP1 + s.PP1I1p;
}

double total_i1(const PP1I1State& s) {
    check_nonnegative(s.I1, "Invalid I1");
    check_nonnegative(s.I1p, "Invalid I1p");
    check_nonnegative(s.PP1I1p, "Invalid PP1I1p");
    return s.I1 + s.I1p + s.PP1I1p;
}

} // namespace omega::v1
'''

Path("pp1_i1.cpp").write_text(content, encoding="utf-8")
