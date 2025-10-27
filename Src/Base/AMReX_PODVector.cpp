#include <AMReX_PODVector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <cmath>

namespace amrex::VectorGrowthStrategy
{
    Real growth_factor = 1.5_rt;
    amrex::GrowthStrategy default_strategy = amrex::GrowthStrategy::Exact;
    amrex::GrowthStrategy override_strategy = amrex::GrowthStrategy::Default;

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
        pp.query_enum_case_insensitive("vector_growth_strategy", default_strategy);
        pp.query_enum_case_insensitive("vector_growth_strategy_override", override_strategy);

        detail::ValidateUserInput();
    }

    void SetGrowthFactor (Real a_factor) {
        growth_factor = a_factor;
        detail::ValidateUserInput();
    }

    std::size_t GetNewCapacity ([[maybe_unused]] std::size_t old_size, std::size_t old_capacity,
                                std::size_t new_size, amrex::GrowthStrategy strategy,
                                std::size_t sizeof_T) {

        std::size_t new_capacity = 0;

        if (override_strategy != amrex::GrowthStrategy::Default) {
            return GetNewCapacity(old_size, old_capacity, new_size, override_strategy, sizeof_T);
        }

        switch (strategy) {
            case GrowthStrategy::Default:
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(default_strategy != GrowthStrategy::Default,
                    "Must set default growth strategy to something other than Default");
                return GetNewCapacity(old_size, old_capacity, new_size, default_strategy, sizeof_T);
            case GrowthStrategy::Exact:
                if (new_size > old_capacity) {
                    new_capacity = new_size;
                } else {
                    new_capacity = old_capacity;
                }
                break;
            case GrowthStrategy::Geometric:
                if (new_size > old_capacity) {
                    const std::size_t min_capacity = std::max(64/sizeof_T, new_size);
                    Real const gf = GetGrowthFactor();
                    if (amrex::almostEqual(gf, Real(1.5))) {
                        new_capacity = std::max((old_capacity*3+1)/2, min_capacity);
                    } else {
                        new_capacity = std::max(std::size_t(gf*Real(old_capacity+1)), min_capacity);
                    }
                } else {
                    new_capacity = old_capacity;
                }
                break;
            case GrowthStrategy::Poisson:
                if (new_size > old_capacity) {
                    new_capacity = new_size + static_cast<std::size_t>(3 * std::sqrt(new_size));
                } else {
                    new_capacity = old_capacity;
                }
                break;
        }

        AMREX_ALWAYS_ASSERT(new_capacity >= new_size);
        return new_capacity;
    }
}
