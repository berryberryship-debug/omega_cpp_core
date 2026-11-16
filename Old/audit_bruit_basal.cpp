cat > audit_bruit_basal.cpp <<'EOF'
#include <iostream>
#include <iomanip>

#include "biophysical_fann_adapter.hpp"

int main() {
    const double dt = 0.001;
    const double duration = 1.0;
    const int steps = static_cast<int>(duration / dt);

    BiophysicalFannAdapter deterministic(
        0.6, 0.0, 42, false
    );

    BiophysicalFannAdapter stochastic(
        0.6, 0.0, 42, true
    );

    std::cout << std::fixed << std::setprecision(8);

    std::cout
        << "t"
        << "\tW_det\tW_stoch"
        << "\tK_det\tK_stoch"
        << "\tP_det\tP_stoch"
        << "\tI_det\tI_stoch\n";

    for (int i = 0; i <= steps; ++i) {
        const double t = i * dt;

        if (i % 100 == 0) {
            std::cout
                << t << '\t'
                << deterministic.get_weight() << '\t'
                << stochastic.get_weight() << '\t'
                << deterministic.get_K() << '\t'
                << stochastic.get_K() << '\t'
                << deterministic.get_P() << '\t'
                << stochastic.get_P() << '\t'
                << deterministic.get_I() << '\t'
                << stochastic.get_I()
                << '\n';
        }

        if (i < steps) {
            deterministic.update_synapse(0.0, 0.0, dt);
            stochastic.update_synapse(0.0, 0.0, dt);
        }
    }

    std::cout << "\n========== DIFFERENCES FINALES ==========\n";

    std::cout
        << "Delta W (stoch - det) : "
        << stochastic.get_weight()
           - deterministic.get_weight() << '\n';

    std::cout
        << "Delta K (stoch - det) : "
        << stochastic.get_K()
           - deterministic.get_K() << '\n';

    std::cout
        << "Delta P (stoch - det) : "
        << stochastic.get_P()
           - deterministic.get_P() << '\n';

    std::cout
        << "Delta I (stoch - det) : "
        << stochastic.get_I()
           - deterministic.get_I() << '\n';

    return 0;
}
EOF
