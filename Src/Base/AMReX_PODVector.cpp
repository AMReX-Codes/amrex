#include <AMReX_PODVector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

namespace amrex::VectorGrowthStrategy
{
    Real growth_factor = 1.5_rt;

    // clamp user input to reasonable values
    constexpr Real min_factor = 1.001_rt;
    constexpr Real max_factor = 4._rt;

    namespace detail
    {
        void ValidateUserInput() {
            if (growth_factor < min_factor) {
                if (Verbose()) {
                    amrex::Print() << "Warning: user-provided vector growth factor is too small."
                                   << " Clamping to " << min_factor << ". \n";
                }
                growth_factor = min_factor;
            }

            if (growth_factor > max_factor) {
                if (Verbose()) {
                    amrex::Print() << "Warning: user-provided vector growth factor is too large."
                                   << " Clamping to " << max_factor << ". \n";
                }
                growth_factor = max_factor;
            }
        }
    }

    void Initialize () {
        ParmParse pp("amrex");
        pp.queryAdd("vector_growth_factor", growth_factor);

        detail::ValidateUserInput();
    }

    void SetGrowthFactor (Real a_factor) {
        growth_factor = a_factor;
        detail::ValidateUserInput();
    }
}
