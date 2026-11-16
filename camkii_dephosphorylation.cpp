
from pathlib import Path

content = r'''#include "camkii_dephosphorylation.hpp"

#include <cmath>
#include <stdexcept>

namespace omega::v1 {
namespace {

void check_nonnegative(double x, const char* name) {
    if (!std::isfinite(x) || x < 0.0)
        throw std::invalid_argument(name);
}

}

void validate_dephosphorylation(
    const DephosphorylationState& s,
    const DephosphorylationParameters& p
) {
    check_nonnegative(s.PP1, "Invalid PP1");
    check_nonnegative(s.KpC, "Invalid KpC");
    check_nonnegative(s.KC, "Invalid KC");
    check_nonnegative(s.PP1KpC, "Invalid PP1KpC");

    check_nonnegative(p.kon, "Invalid kon");
    check_nonnegative(p.koff, "Invalid koff");
    check_nonnegative(p.kcat, "Invalid kcat");
}

DephosphorylationRates calculate_dephosphorylation_rates(
    const DephosphorylationState& s,
    const DephosphorylationParameters& p
) {
    validate_dephosphorylation(s, p);

    return {
        p.kon * s.PP1 * s.KpC,
        p.koff * s.PP1KpC,
        p.kcat * s.PP1KpC
    };
}

DephosphorylationDerivative calculate_dephosphorylation_derivative(
    const DephosphorylationState& s,
    const DephosphorylationParameters& p
) {
    const auto r = calculate_dephosphorylation_rates(s, p);

    const double Jassociation = r.association;
    const double Jdissociation = r.dissociation;
    const double Jcatalysis = r.catalysis;

    return {
        -Jassociation + Jdissociation + Jcatalysis,
        -Jassociation + Jdissociation,
        Jcatalysis,
        Jassociation - Jdissociation - Jcatalysis
    };
}

double total_dephosphorylation_pp1(
    const DephosphorylationState& s
) {
    check_nonnegative(s.PP1, "Invalid PP1");
    check_nonnegative(s.PP1KpC, "Invalid PP1KpC");
    return s.PP1 + s.PP1KpC;
}

double total_dephosphorylation_camkii(
    const DephosphorylationState& s
) {
    check_nonnegative(s.KpC, "Invalid KpC");
    check_nonnegative(s.KC, "Invalid KC");
    check_nonnegative(s.PP1KpC, "Invalid PP1KpC");
    return s.KpC + s.KC + s.PP1KpC;
}

} // namespace omega::v1
'''

Path("camkii_dephosphorylation.cpp").write_text(content, encoding="utf-8")
