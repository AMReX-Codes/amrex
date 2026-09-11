//
// Unit test for the extdir/hoextrap-aware EB slope kernels
//     amrex_calc_slopes_extdir_eb        and
//     amrex_calc_slopes_extdir_eb_grown
//
// The data is an exactly linear profile q = sum_d coef[d] * x_d, where x_d is
// measured in units of dx and the cell centroids coincide with the cell centers
// (ccent = 0).  In the direction "dir" the domain has an ext_dir boundary on both
// sides, so the value in the first ghost cell lives *on* the domain face rather
// than at the cell center; cells further outside the domain are filled with
// garbage to make sure they never enter the stencil.
//
// Every kernel must return the exact slopes coef[] for such a profile, including
// at the cells one in from the boundary (domlo+1 / domhi-1), which the 4th-order
// one-sided stencil and the grown (n == 2) least squares stencil both reach.
//
#include <AMReX.H>
#include <AMReX_BaseFab.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EB_Slopes_K.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>

#include <cmath>

using namespace amrex;

namespace {

constexpr int ncell = 16;

// Coefficients of the exactly linear test profile
constexpr Real coef_x = Real(1.5);
constexpr Real coef_y = Real(-0.75);
constexpr Real coef_z = Real(0.25);

// A value that must never make it into any stencil
constexpr Real garbage = Real(12345.);

#ifdef BL_USE_FLOAT
constexpr Real tol = Real(1.e-4);
#else
constexpr Real tol = Real(1.e-11);
#endif

int
test_dir (int dir)
{
    const Box domain(IntVect(0), IntVect(AMREX_D_DECL(ncell-1,ncell-1,ncell-1)));
    const Box gbx = amrex::grow(domain, 3);

    const IntVect dlo = domain.smallEnd();
    const IntVect dhi = domain.bigEnd();

    const int domlo = dlo[dir];
    const int domhi = dhi[dir];

    FArrayBox state_fab(gbx, 1, The_Arena());
    FArrayBox ccent_fab(gbx, AMREX_SPACEDIM, The_Arena());
    FArrayBox vfrac_reg_fab(gbx, 1, The_Arena());
    FArrayBox vfrac_cut_fab(gbx, 1, The_Arena());
    BaseFab<EBCellFlag> flag_fab(gbx, 1, The_Arena());

    // Face centroids are all at the center of the face
    AMREX_D_TERM(FArrayBox fcx_fab(amrex::surroundingNodes(gbx,0), AMREX_SPACEDIM-1, The_Arena());,
                 FArrayBox fcy_fab(amrex::surroundingNodes(gbx,1), AMREX_SPACEDIM-1, The_Arena());,
                 FArrayBox fcz_fab(amrex::surroundingNodes(gbx,2), AMREX_SPACEDIM-1, The_Arena()););

    ccent_fab.setVal<RunOn::Device>(Real(0.));
    AMREX_D_TERM(fcx_fab.setVal<RunOn::Device>(Real(0.));,
                 fcy_fab.setVal<RunOn::Device>(Real(0.));,
                 fcz_fab.setVal<RunOn::Device>(Real(0.)););

    // vfrac == 1 everywhere selects the "regular slope" path inside the kernels;
    // vfrac < 1 everywhere forces the EB least squares fit to be used instead.
    vfrac_reg_fab.setVal<RunOn::Device>(Real(1.));
    vfrac_cut_fab.setVal<RunOn::Device>(Real(0.5));

    auto const& state = state_fab.array();
    auto const& flag  = flag_fab.array();

    ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        EBCellFlag f;
        f.setRegular();
        f.setConnected();
        flag(i,j,k) = f;

        int idx[3] = {0,0,0};
        AMREX_D_TERM(idx[0] = i;, idx[1] = j;, idx[2] = k;);

        if (idx[dir] < domlo-1 || idx[dir] > domhi+1) {
            state(i,j,k,0) = garbage;
        } else {
            Real x[3] = {Real(idx[0]), Real(idx[1]), Real(idx[2])};
            // The first ghost cell outside an ext_dir boundary holds the face value
            if (idx[dir] == domlo-1) { x[dir] = Real(domlo) - Real(0.5); }
            if (idx[dir] == domhi+1) { x[dir] = Real(domhi) + Real(0.5); }
            state(i,j,k,0) = coef_x*x[0] + coef_y*x[1] + coef_z*x[2];
        }
    });

    // The cells we evaluate the slopes at: both boundary cells, both cells one in
    // from the boundary, and an interior cell as a control.
    const Vector<int> h_pos = {domlo, domlo+1, ncell/2, domhi-1, domhi};
    const int npos = static_cast<int>(h_pos.size());

    Gpu::DeviceVector<int> d_pos(npos);
    Gpu::copyAsync(Gpu::hostToDevice, h_pos.begin(), h_pos.end(), d_pos.begin());

    // Two kernels tested at each position: 0 = 3^DIM stencil, 1 = grown stencil
    const int ntest = 2*npos;
    Gpu::DeviceVector<Real> d_slopes(std::size_t(ntest)*AMREX_SPACEDIM);

    int const* pos_p = d_pos.data();
    Real* slopes_p = d_slopes.data();

    auto const& s   = state_fab.const_array();
    auto const& cc  = ccent_fab.const_array();
    auto const& vr  = vfrac_reg_fab.const_array();
    auto const& vc  = vfrac_cut_fab.const_array();
    auto const& fl  = flag_fab.const_array();
    AMREX_D_TERM(auto const& fx = fcx_fab.const_array();,
                 auto const& fy = fcy_fab.const_array();,
                 auto const& fz = fcz_fab.const_array(););

    Gpu::streamSynchronize();

    ParallelFor(ntest, [=] AMREX_GPU_DEVICE (int t) noexcept
    {
        const int ip   = t / 2;
        const int kern = t % 2;

        IntVect iv(AMREX_D_DECL(ncell/2,ncell/2,ncell/2));
        iv[dir] = pos_p[ip];

        AMREX_D_TERM(const int i = iv[0];,
                     const int j = iv[1];,
                     const int k = iv[2];);
#if (AMREX_SPACEDIM == 2)
        const int k = 0;
#endif

        bool edlo[3] = {false,false,false};
        bool edhi[3] = {false,false,false};
        edlo[dir] = true;
        edhi[dir] = true;

        // Grow the least squares stencil in the direction being tested
        int nn[3] = {1,1,1};
        nn[dir] = 2;

        GpuArray<Real,AMREX_SPACEDIM> sl;
        if (kern == 0) {
            sl = amrex_calc_slopes_extdir_eb(i,j,k,0,s,cc,vr,
                                             AMREX_D_DECL(fx,fy,fz),fl,
                                             AMREX_D_DECL(edlo[0],edlo[1],edlo[2]),
                                             AMREX_D_DECL(edhi[0],edhi[1],edhi[2]),
                                             AMREX_D_DECL(dlo[0],dlo[1],dlo[2]),
                                             AMREX_D_DECL(dhi[0],dhi[1],dhi[2]),
                                             4);
        } else {
            sl = amrex_calc_slopes_extdir_eb_grown(i,j,k,0,
                                                   AMREX_D_DECL(nn[0],nn[1],nn[2]),
                                                   s,cc,vc,
                                                   AMREX_D_DECL(fx,fy,fz),fl,
                                                   AMREX_D_DECL(edlo[0],edlo[1],edlo[2]),
                                                   AMREX_D_DECL(edhi[0],edhi[1],edhi[2]),
                                                   AMREX_D_DECL(dlo[0],dlo[1],dlo[2]),
                                                   AMREX_D_DECL(dhi[0],dhi[1],dhi[2]),
                                                   4);
        }

        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            slopes_p[t*AMREX_SPACEDIM+d] = sl[d];
        }
    });

    Vector<Real> h_slopes(std::size_t(ntest)*AMREX_SPACEDIM);
    Gpu::copyAsync(Gpu::deviceToHost, d_slopes.begin(), d_slopes.end(), h_slopes.begin());
    Gpu::streamSynchronize();

    const Real expected[3] = {coef_x, coef_y, coef_z};

    int nfail = 0;
    for (int t = 0; t < ntest; ++t) {
        const int ip   = t / 2;
        const int kern = t % 2;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const Real got = h_slopes[std::size_t(t)*AMREX_SPACEDIM+d];
            if (std::abs(got-expected[d]) > tol*(Real(1.)+std::abs(expected[d]))) {
                ++nfail;
                amrex::Print() << "FAIL: dir = " << dir
                               << ", index = " << h_pos[ip]
                               << ", kernel = "
                               << (kern == 0 ? "amrex_calc_slopes_extdir_eb"
                                             : "amrex_calc_slopes_extdir_eb_grown")
                               << ", slope[" << d << "] = " << got
                               << ", expected " << expected[d] << '\n';
            }
        }
    }

    return nfail;
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        int nfail = 0;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            nfail += test_dir(dir);
        }
        if (nfail == 0) {
            amrex::Print() << "EB extdir slopes: PASS\n";
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nfail == 0, "EB extdir slopes test failed");
    }
    amrex::Finalize();
}
