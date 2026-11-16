#include "coupled.hpp"

namespace omega::v1 {

CoupledDerivative calculate_coupled_derivative(
    const CoupledState& s,
    double calcium,
    const CoupledParameters& p
) {
    // M1 : évolution de la calmoduline libre.
    const CaMDerivative dm1 =
        calculate_derivative(s.cam, calcium, p.cam);

    // M2 : M4 est lu depuis M1, sans seconde copie persistante.
    const CaMKIIState kinase{
        s.K0,
        s.KC,
        s.D,
        s.KpC,
        s.Kp,
        s.cam.M4
    };

    const CaMKIIDerivative dm2 =
        calculate_camkii_derivative(kinase, p.camkii);

    return {
        dm1,
        dm2.K0,
        dm2.KC,
        dm2.D,
        dm2.KpC,
        dm2.Kp,
        dm1.M4 + dm2.M4
    };
}

double total_coupled_camkii(const CoupledState& s) {
    const CaMKIIState kinase{
        s.K0,
        s.KC,
        s.D,
        s.KpC,
        s.Kp,
        s.cam.M4
    };

    return total_camkii(kinase);
}

double total_coupled_calmodulin(const CoupledState& s) {
    return s.cam.M0
         + s.cam.M2N
         + s.cam.M2C
         + s.cam.M4
         + s.KC
         + 2.0 * s.D
         + s.KpC;
}

} // namespace omega::v1
