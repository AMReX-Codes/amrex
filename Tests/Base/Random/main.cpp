#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>

#include <cstdint>
#include <limits>
#include <random>

using namespace amrex;

namespace {

// A URBG pinned to its maximum output: the input that drives
// std::uniform_real_distribution hardest towards its excluded upper bound.
// LWG 2524 notes that std::generate_canonical may return exactly 1.0, because
// the quotient it forms can round up. Every standard library AMReX currently
// supports clamps instead (MSVC since it implemented P0952R2 in VS 2022
// 17.12), so on those this returns the largest value below one rather than one
// itself -- but nothing in the standard AMReX targets requires that, which is
// why Random() does not depend on it.
struct MaxEngine
{
    using result_type = std::uint32_t;
    static constexpr result_type min () { return 0; }
    static constexpr result_type max () { return std::numeric_limits<std::uint32_t>::max(); }
    result_type operator() () { return max(); }
};

// The interval conversions that amrex::Random and amrex::RandomPositive apply
// to a raw draw must hold the documented intervals for every input, including
// the endpoints. These checks are deterministic: the endpoints are far too
// rare to reach by sampling, so they are fed in directly.
void test_endpoint_guard ()
{
    // Whatever the standard library hands back for the hardest input, the
    // guard must leave it inside [0,1).
    MaxEngine engine;
    std::uniform_real_distribution<Real> distribution(0.0, 1.0);
    Real const raw = distribution(engine);
    Real const guarded = random_util::clamp_below_one(raw);
    AMREX_ALWAYS_ASSERT(guarded < Real(1));
    AMREX_ALWAYS_ASSERT(guarded >= Real(0));

    // The clamp must be a no-op on every value that is already in range.
    AMREX_ALWAYS_ASSERT(random_util::clamp_below_one(Real(0)) == Real(0));
    AMREX_ALWAYS_ASSERT(random_util::clamp_below_one(Real(0.5)) == Real(0.5));
    Real const below_one = Real(1) - Real(0.5)*std::numeric_limits<Real>::epsilon();
    AMREX_ALWAYS_ASSERT(random_util::clamp_below_one(below_one) == below_one);

    // The clamped upper bound is the largest Real below one, so nothing is
    // needlessly discarded.
    AMREX_ALWAYS_ASSERT(random_util::clamp_below_one(Real(1)) == below_one);

    // The interval conversions relocate exactly one endpoint and leave every
    // other value untouched, so they cannot lose resolution.
    AMREX_ALWAYS_ASSERT(random_util::one_to_zero(Real(1)) == Real(0));   // (0,1] -> [0,1)
    AMREX_ALWAYS_ASSERT(random_util::one_to_zero(below_one) == below_one);
    AMREX_ALWAYS_ASSERT(random_util::one_to_zero(Real(0.5)) == Real(0.5));

    AMREX_ALWAYS_ASSERT(random_util::zero_to_one(Real(0)) == Real(1));   // [0,1) -> (0,1]
    AMREX_ALWAYS_ASSERT(random_util::zero_to_one(below_one) == below_one);
    AMREX_ALWAYS_ASSERT(random_util::zero_to_one(Real(0.5)) == Real(0.5));

    // zero_to_one needs no clamp of its own: even a draw of exactly 1.0, which
    // would violate the generator's own [0,1) contract, is inside (0,1].
    AMREX_ALWAYS_ASSERT(random_util::zero_to_one(Real(1)) > Real(0));
    AMREX_ALWAYS_ASSERT(random_util::zero_to_one(Real(1)) <= Real(1));
    AMREX_ALWAYS_ASSERT(random_util::zero_to_one(raw) > Real(0));
}

// The intervals themselves, over the host generators.
void test_host_intervals ()
{
    for (int i = 0; i < 1000000; ++i) {
        Real const r = amrex::Random();
        AMREX_ALWAYS_ASSERT(r >= Real(0) && r < Real(1));

        Real const p = amrex::RandomPositive();
        AMREX_ALWAYS_ASSERT(p > Real(0) && p <= Real(1));
    }
}

// The same intervals through the RandomEngine overloads, which on a GPU build
// runs the device code paths. Note this part is statistical, not
// deterministic: the endpoints are rare enough that it would take far more
// draws than this to hit one, so it guards against a wholesale breakage of a
// path rather than against a regression at the endpoint itself -- that is what
// test_endpoint_guard above is for. RandomGamma is included because it
// consumes RandomPositive internally and would return a non-finite value if a
// zero ever reached its std::log.
void test_engine_intervals ()
{
    Long const N = 2000000;

    Gpu::DeviceScalar<Long> d_bad_random(0), d_bad_positive(0), d_bad_gamma(0);
    Long* p_bad_random   = d_bad_random.dataPtr();
    Long* p_bad_positive = d_bad_positive.dataPtr();
    Long* p_bad_gamma    = d_bad_gamma.dataPtr();

    amrex::ParallelForRNG(N,
    [=] AMREX_GPU_DEVICE (Long, RandomEngine const& engine) noexcept
    {
        Real const r = amrex::Random(engine);
        if (!(r >= Real(0) && r < Real(1))) {
            Gpu::Atomic::AddNoRet(p_bad_random, Long(1));
        }

        Real const p = amrex::RandomPositive(engine);
        if (!(p > Real(0) && p <= Real(1))) {
            Gpu::Atomic::AddNoRet(p_bad_positive, Long(1));
        }

        Real const g1 = amrex::RandomGamma(Real(2.5), Real(1.0), engine);
        Real const g2 = amrex::RandomGamma(Real(0.4), Real(1.0), engine); // alpha < 1 branch
        if (!(g1 > Real(0)) || !(g2 > Real(0))) {
            Gpu::Atomic::AddNoRet(p_bad_gamma, Long(1));
        }
    });
    Gpu::streamSynchronize();

    AMREX_ALWAYS_ASSERT(d_bad_random.dataValue() == 0);
    AMREX_ALWAYS_ASSERT(d_bad_positive.dataValue() == 0);
    AMREX_ALWAYS_ASSERT(d_bad_gamma.dataValue() == 0);
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    test_endpoint_guard();
    test_host_intervals();
    test_engine_intervals();
    amrex::Print() << "SUCCESS\n";
    amrex::Finalize();
}
