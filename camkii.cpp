cat > camkii.cpp <<'EOF'
#include "camkii.hpp"

#include <cmath>
#include <stdexcept>

namespace omega::v1 {

namespace {

void check_nonnegative(double value, const char* name) {
    if (!std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument(name);
    }
}

} // namespace

void validate_camkii(
    const CaMKIIState& s,
    const CaMKIIParameters& p
) {
    check_nonnegative(s.K0, "Invalid K0");
    check_nonnegative(s.KC, "Invalid KC");
    check_nonnegative(s.D, "Invalid D");
    check_nonnegative(s.KpC, "Invalid KpC");
    check_nonnegative(s.Kp, "Invalid Kp");
    check_nonnegative(s.M4, "Invalid M4");

    check_nonnegative(p.kon_bind, "Invalid kon_bind");
    check_nonnegative(p.koff_bind, "Invalid koff_bind");
    check_nonnegative(p.kon_dimer, "Invalid kon_dimer");
    check_nonnegative(p.koff_dimer, "Invalid koff_dimer");
    check_nonnegative(
        p.k_phosphorylation,
        "Invalid k_phosphorylation"
    );
    check_nonnegative(
        p.koff_phospho_cam,
        "Invalid koff_phospho_cam"
    );
    check_nonnegative(
        p.kon_phospho_cam,
        "Invalid kon_phospho_cam"
    );
}

CaMKIIRates calculate_camkii_rates(
    const CaMKIIState& s,
    const CaMKIIParameters& p
) {
    validate_camkii(s, p);

    return {
        p.kon_bind * s.K0 * s.M4,
        p.koff_bind * s.KC,

        p.kon_dimer * s.KC * s.KC,
        p.koff_dimer * s.D,

        p.k_phosphorylation * s.D,

        p.koff_phospho_cam * s.KpC,
        p.kon_phospho_cam * s.Kp * s.M4
    };
}

CaMKIIDerivative calculate_camkii_derivative(
    const CaMKIIState& s,
    const CaMKIIParameters& p
) {
    const CaMKIIRates r = calculate_camkii_rates(s, p);

    const double Jbind =
        r.bind_forward - r.bind_reverse;

    const double Jdimer =
        r.dimer_forward - r.dimer_reverse;

    const double Jphospho = r.phosphorylation;

    const double Jphospho_cam =
        r.phospho_cam_forward - r.phospho_cam_reverse;

    return {
        -Jbind,
        Jbind - 2.0 * Jdimer + Jphospho,
        Jdimer - Jphospho,
        Jphospho - Jphospho_cam,
        Jphospho_cam,
        -Jbind + Jphospho_cam
    };
}

double total_camkii(const CaMKIIState& s) {
    check_nonnegative(s.K0, "Invalid K0");
    check_nonnegative(s.KC, "Invalid KC");
    check_nonnegative(s.D, "Invalid D");
    check_nonnegative(s.KpC, "Invalid KpC");
    check_nonnegative(s.Kp, "Invalid Kp");

    return s.K0 + s.KC + 2.0 * s.D + s.KpC + s.Kp;
}

} // namespace omega::v1
EOF
