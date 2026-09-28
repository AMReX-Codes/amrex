#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#if defined(AMREX_USE_EB) && (AMREX_SPACEDIM > 1)
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EBFabFactory.H>
#define NAN_TEST_EB 1
#endif
#if __has_include(<AMReX_AlgVector.H>)
#include <AMReX_AlgVector.H>
#define NAN_TEST_ALGVECTOR 1
#endif

#include <bit>
#include <cstdint>
#include <cstring>
#include <limits>
#include <type_traits>

//
// Tests amrex::isnan, isinf and isfinite, which must work under fast math
// (the Small ctest set runs in the fast-math CI jobs), and that norminf
// reports a NaN as infinity.
//
// Under fast math, clang treats a NaN or inf passed to or returned from a
// function as undefined, and drops it if it can prove the value is one.  So
// the special values are passed around as bit patterns that depend on a
// run-time value, and turned into floating-point numbers in memory only.
//

using namespace amrex;

namespace {

int nfailures = 0;

void check (bool ok, char const* what)
{
    if (!ok) {
        ++nfailures;
        amrex::Print() << "FAILED: " << what << "\n";
    }
}

template <typename T>
struct FPBits
{
    using U = std::conditional_t<sizeof(T) == 8, std::uint64_t, std::uint32_t>;
    static constexpr int nman = std::numeric_limits<T>::digits - 1;
    static constexpr U sign = U(1) << (sizeof(T)*8-1);
    static constexpr U expo = ((U(1) << (sizeof(T)*8-1-nman)) - 1) << nman;
    static constexpr U quiet = U(1) << (nman-1);
};

template <typename T>
AMREX_FORCE_INLINE bool is_pos_inf (T x)
{
    return amrex::isinf(x) && x > T(0);
}

// bit 0: isnan, bit 1: isinf, bit 2: isfinite
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int classify (T const* p) noexcept
{
    return int(amrex::isnan(*p)) | (int(amrex::isinf(*p)) << 1)
        | (int(amrex::isfinite(*p)) << 2);
}

template <typename T>
void test_classify (int zero)
{
    using B = FPBits<T>;
    using U = typename B::U;
    constexpr U sign = B::sign;
    constexpr U expo = B::expo;
    U const z = U(zero);

    constexpr int n = 11;
    GpuArray<U,n> const bits{(expo | B::quiet) + z, // qNaN
                             (sign | expo | B::quiet) + z, // -qNaN
                             (expo | U(1)) + z, // NaN, payload 1
                             expo + z, // +inf
                             (sign | expo) + z, // -inf
                             z, // 0
                             sign + z, // -0
                             U(1) + z, // smallest denormal
                             (expo - U(1)) + z, // max
                             (sign | (expo - U(1))) + z, // lowest
                             std::bit_cast<U>(T(1.5)) + z};
    GpuArray<int,n> const expected{1, 1, 1, 2, 2, 4, 4, 4, 4, 4, 4};

    Gpu::HostVector<T> h_x(n);
    std::memcpy(h_x.data(), bits.data(), sizeof(T)*n);
    for (int i = 0; i < n; ++i) {
        T const* p = h_x.data() + i;
        check(classify(p) == expected[i], "host classification");
        check(Math::isnan(*p) == amrex::isnan(*p) &&
              Gpu::isinf(*p) == amrex::isinf(*p) &&
              Math::isfinite(*p) == amrex::isfinite(*p), "aliases");
    }

    Gpu::DeviceVector<T> d_x(n);
    Gpu::copy(Gpu::hostToDevice, h_x.begin(), h_x.end(), d_x.begin());
    Gpu::DeviceVector<int> d_result(n);
    auto const* px = d_x.data();
    auto* pr = d_result.data();
    amrex::ParallelFor(n, [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        pr[i] = classify(px+i);
    });
    Gpu::HostVector<int> h_result(n);
    Gpu::copy(Gpu::deviceToHost, d_result.begin(), d_result.end(), h_result.begin());
    for (int i = 0; i < n; ++i) {
        check(h_result[i] == expected[i], "device classification");
    }
}

using RBits = FPBits<Real>;

// mf(iv,comp) = the number with the given bits
void set_cell (MultiFab& mf, IntVect const& iv, int comp, RBits::U bits)
{
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        if (mfi.fabbox().contains(iv)) {
            auto const& a = mf.array(mfi);
            amrex::ParallelFor(Box(iv,iv), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Gpu::memcpy(a.ptr(i,j,k,comp), &bits, sizeof(Real));
            });
        }
    }
    Gpu::streamSynchronize();
}

void set_mask (iMultiFab& mask, IntVect const& iv, int v)
{
    for (MFIter mfi(mask); mfi.isValid(); ++mfi) {
        if (mfi.fabbox().contains(iv)) {
            auto const& a = mask.array(mfi);
            amrex::ParallelFor(Box(iv,iv), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                a(i,j,k) = v;
            });
        }
    }
    Gpu::streamSynchronize();
}

void test_norminf (int zero)
{
    BoxArray ba(Box(IntVect(0), IntVect(AMREX_D_DECL(15,7,7))));
    ba.maxSize(8); // two boxes
    DistributionMapping dm(ba);
    MultiFab mf(ba, dm, 2, 1);
    mf.setVal(Real(1.5));

    using U = RBits::U;
    U const z = U(zero);
    U const nan = (RBits::expo | RBits::quiet) + z;
    U const minus_inf = (RBits::sign | RBits::expo) + z;
    U const one_half = std::bit_cast<U>(Real(1.5));
    Real const vmax = Real(3.5);
    IntVect const lo(0);
    IntVect const cell(AMREX_D_DECL(12,3,3));
    IntVect const ghost(AMREX_D_DECL(-1,0,0));

    set_cell(mf, lo, 0, std::bit_cast<U>(-vmax));
    check(mf.norminf(0, 2, IntVect(0)) == vmax, "finite norminf");

    set_cell(mf, cell, 1, nan);
    check(is_pos_inf(mf.norminf(0, 2, IntVect(0))), "NaN in valid cell");
    check(is_pos_inf(mf.norminf(0, 2, IntVect(1))), "NaN in valid cell, nghost=1");
    check(mf.norminf(0, 1, IntVect(0)) == vmax, "NaN in another component");
    check(is_pos_inf(mf.norm0(1)), "norm0");
    check(is_pos_inf(amrex::norminf(mf, 0, 2, IntVect(0))), "free norminf");
    check(mf.contains_nan(0, 2, 0), "contains_nan");

    iMultiFab mask(ba, dm, 1, 0);
    mask.setVal(1);
    check(is_pos_inf(mf.norminf(mask, 0, 2, IntVect(0))), "mask norminf");
    set_mask(mask, cell, 0);
    check(mf.norminf(mask, 0, 2, IntVect(0)) == vmax, "NaN masked out");

    set_cell(mf, cell, 1, one_half);
    set_cell(mf, ghost, 1, nan);
    check(mf.norminf(0, 2, IntVect(0)) == vmax, "NaN in ghost cell, nghost=0");
    check(is_pos_inf(mf.norminf(0, 2, IntVect(1))), "NaN in ghost cell, nghost=1");
    check(!mf.contains_nan(0, 2, 0), "contains_nan, nghost=0");

    set_cell(mf, ghost, 1, one_half);
    set_cell(mf, lo, 0, nan);
    check(is_pos_inf(mf.norminf(0, 2, IntVect(0))), "NaN in first cell");

    set_cell(mf, lo, 0, minus_inf);
    check(is_pos_inf(mf.norminf(0, 2, IntVect(0))), "-inf");
}

#ifdef NAN_TEST_EB
void test_norminf_eb (int zero)
{
    Box const domain(IntVect(0), IntVect(15));
    Geometry geom(domain, RealBox(AMREX_D_DECL(0.,0.,0.), AMREX_D_DECL(1.,1.,1.)),
                  CoordSys::cartesian, {AMREX_D_DECL(0,0,0)});
    EB2::PlaneIF plane({AMREX_D_DECL(Real(0.4),Real(0.),Real(0.))},
                       {AMREX_D_DECL(Real(1.),Real(0.),Real(0.))});
    EB2::Build(EB2::makeShop(plane), geom, 0, 0);
    BoxArray ba(domain);
    ba.maxSize(8);
    DistributionMapping dm(ba);
    auto factory = makeEBFabFactory(geom, ba, dm, {1,1,1}, EBSupport::basic);
    MultiFab mf(ba, dm, 1, 0, MFInfo(), *factory);
    mf.setVal(Real(1.5));
    auto const& flags = factory->getMultiEBCellFlagFab();

    // Fill covered or cut cells with NaN
    auto fill_nan = [&] (bool covered) {
        RBits::U const nan = (RBits::expo | RBits::quiet) + RBits::U(zero);
        auto const& ma = mf.arrays();
        auto const& fma = flags.const_arrays();
        amrex::ParallelFor(mf, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
        {
            auto const flag = fma[b](i,j,k);
            if (covered ? flag.isCovered() : flag.isSingleValued()) {
                Gpu::memcpy(ma[b].ptr(i,j,k), &nan, sizeof(Real));
            }
        });
        Gpu::streamSynchronize();
    };

    fill_nan(true);
    check(mf.norminf(0, 1, IntVect(0), false, true) == Real(1.5), "EB: NaN in covered cells ignored");
    check(is_pos_inf(mf.norminf(0, 1, IntVect(0), false, false)), "EB: NaN in covered cells");
    fill_nan(false);
    check(is_pos_inf(mf.norminf(0, 1, IntVect(0), false, true)), "EB: NaN in cut cells");
}
#endif

#ifdef NAN_TEST_ALGVECTOR
void test_norminf_algvector (int zero)
{
    AlgVector<Real> v(Long(100));
    v.setVal(Real(-2.5));
    check(v.norminf() == Real(2.5), "AlgVector: finite norminf");
    RBits::U const nan = (RBits::expo | RBits::quiet) + RBits::U(zero);
    if (v.numLocalRows() > 0) {
        Real* p = v.data();
        amrex::single_task([=] AMREX_GPU_DEVICE () noexcept
        {
            Gpu::memcpy(p, &nan, sizeof(Real));
        });
        Gpu::streamSynchronize();
    }
    check(is_pos_inf(v.norminf()), "AlgVector: NaN");
}
#endif

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        int zero = 0;
        ParmParse pp;
        pp.query("zero", zero);

        test_classify<float>(zero);
        test_classify<double>(zero);
        test_norminf(zero);
#ifdef NAN_TEST_EB
        test_norminf_eb(zero);
#endif
#ifdef NAN_TEST_ALGVECTOR
        test_norminf_algvector(zero);
#endif

        if (nfailures > 0) {
            amrex::Abort("NaN test failed");
        }
        amrex::Print() << "NaN test passed\n";
    }
    amrex::Finalize();
}
