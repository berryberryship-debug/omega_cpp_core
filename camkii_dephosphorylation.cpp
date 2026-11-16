#include "camkii_dephosphorylation.hpp"

#include <cmath>
#include <stdexcept>

namespace omega::v1 {
namespace {

// Fonction utilitaire pour la validation initiale uniquement
void check_nonnegative(double value, const char* name) {
    if (!std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument(name);
    }
}

} // namespace

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
    // ATTENTION : Pas de validate_dephosphorylation() ici !
    // L'intégrateur appelle cette fonction des milliers de fois. 
    // Cela doit rester une opération mathématique pure et rapide.

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
    const DephosphorylationRates r =
        calculate_dephosphorylation_rates(s, p);

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
    return s.PP1 + s.PP1KpC;
}

double total_dephosphorylation_camkii(
    const DephosphorylationState& s
) {
    return s.KpC + s.KC + s.PP1KpC;
}

} // namespace omega::v1
