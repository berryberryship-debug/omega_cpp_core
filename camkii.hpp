cd ~/axiomes_base_labs_v8/omega_v1/m2_camkii || exit 1

cat > camkii.hpp <<'EOF'
#ifndef OMEGA_V1_CAMKII_HPP
#define OMEGA_V1_CAMKII_HPP

namespace omega::v1 {

struct CaMKIIState {
    double K0;
    double KC;
    double D;
    double KpC;
    double Kp;
    double M4;
};

struct CaMKIIParameters {
    double kon_bind;
    double koff_bind;

    double kon_dimer;
    double koff_dimer;

    double k_phosphorylation;

    double koff_phospho_cam;
    double kon_phospho_cam;
};

struct CaMKIIDerivative {
    double K0;
    double KC;
    double D;
    double KpC;
    double Kp;
    double M4;
};

struct CaMKIIRates {
    double bind_forward;
    double bind_reverse;

    double dimer_forward;
    double dimer_reverse;

    double phosphorylation;

    double phospho_cam_forward;
    double phospho_cam_reverse;
};

CaMKIIRates calculate_camkii_rates(
    const CaMKIIState& state,
    const CaMKIIParameters& parameters
);

CaMKIIDerivative calculate_camkii_derivative(
    const CaMKIIState& state,
    const CaMKIIParameters& parameters
);

double total_camkii(const CaMKIIState& state);

void validate_camkii(
    const CaMKIIState& state,
    const CaMKIIParameters& parameters
);

} // namespace omega::v1

#endif
EOF
