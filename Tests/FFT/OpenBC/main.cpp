#include <AMReX_FFT_Poisson.H> // Put this at the top for testing

#include <AMReX.H>
#include <AMReX_FFT_OpenBCSolver.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Reduce.H>

using namespace amrex;

namespace {

#if (AMREX_SPACEDIM == 3)

void fill_rhs (MultiFab& rho, Geometry const& geom, IndexType ixtype)
{
    auto const& dx = geom.CellSizeArray();
    auto const& problo = geom.ProbLoArray();
    auto const& rhoma = rho.arrays();

    constexpr int nsub = 4;
    Real dxsub = dx[0]/nsub;
    Real dysub = dx[1]/nsub;
    Real dzsub = dx[2]/nsub;

    ParallelFor(rho, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        Real x = (Real(i)+0.5_rt/nsub)*dx[0] + problo[0];
        Real y = (Real(j)+0.5_rt/nsub)*dx[1] + problo[1];
        Real z = (Real(k)+0.5_rt/nsub)*dx[2] + problo[2];
        if (ixtype.nodeCentered()) {
            x -= 0.5_rt*dx[0];
            y -= 0.5_rt*dx[1];
            z -= 0.5_rt*dx[2];
        }
        int n = 0;
        for (int isub = 0; isub < nsub; ++isub) {
        for (int jsub = 0; jsub < nsub; ++jsub) {
        for (int ksub = 0; ksub < nsub; ++ksub) {
            auto xs = x + Real(isub)*dxsub;
            auto ys = y + Real(jsub)*dysub;
            auto zs = z + Real(ksub)*dzsub;
            if ((xs*xs+ys*ys+zs*zs) < 0.25_rt) { ++n; }
        }}}
        rhoma[b](i,j,k) = Real(n) / Real(nsub*nsub*nsub);
    });
}

#endif

// OpenBCSolver computes phi(i) = sum_j G(|i-j|) rho(j), so a direct sum with a
// simple G is an exact reference that does not depend on the FFT machinery.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real test_greens_function (int i, int j, int k)
{
    return Real(1) / (Real(1) + Real(i*i+j*j+k*k));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real test_rhs (int i, int j, int k)
{
    auto h = Real((i*7919 + j*104729 + k*15485863) % 1013);
    return h / Real(1013) - Real(0.5);
}

// Solve with OpenBCSolver and compare against the direct convolution. This is
// the only in-tree coverage of OpenBCSolver in 2D, and of a domain whose
// smallEnd is not zero.
void test_convolution (Box const& domain, int max_grid_size)
{
    amrex::Print() << "\nTesting OpenBCSolver on " << domain
                   << " against a direct convolution\n";

    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    DistributionMapping dm(ba);

    MultiFab rho(ba,dm,1,0), phi(ba,dm,1,0);
    for (MFIter mfi(rho); mfi.isValid(); ++mfi) {
        auto const& a = rho.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            a(i,j,k) = test_rhs(i,j,k);
        });
    }

    auto const lo = amrex::lbound(domain);
    FFT::OpenBCSolver<Real> solver(domain);
    // The Green's function is called with absolute indices.
    solver.setGreensFunction([=] AMREX_GPU_DEVICE (int i, int j, int k) -> Real
    {
        return test_greens_function(i-lo.x, j-lo.y, k-lo.z);
    });
    solver.solve(phi, rho);

    // The reference is analytic, so each rank can check its own cells in place.
    auto const hi = amrex::ubound(domain);
    ReduceOps<ReduceOpMax, ReduceOpMax, ReduceOpSum> reduce_op;
    ReduceData<Real, Real, Long> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        auto const& a = phi.const_array(mfi);
        reduce_op.eval(mfi.validbox(), reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real exact = 0;
            for (int kk = lo.z; kk <= hi.z; ++kk) {
            for (int jj = lo.y; jj <= hi.y; ++jj) {
            for (int ii = lo.x; ii <= hi.x; ++ii) {
                exact += test_greens_function(amrex::Math::abs(i-ii),
                                              amrex::Math::abs(j-jj),
                                              amrex::Math::abs(k-kk)) * test_rhs(ii,jj,kk);
            }}}
            // Math::max would keep the other operand if a(i,j,k) were NaN, so
            // count the non-finite values separately.
            return {amrex::Math::abs(a(i,j,k)-exact), amrex::Math::abs(exact),
                    Long(amrex::isnan(a(i,j,k)) || amrex::isinf(a(i,j,k)))};
        });
    }

    auto hv = reduce_data.value(reduce_op);
    auto errmax = amrex::get<0>(hv);
    auto refmax = amrex::get<1>(hv);
    auto nbad   = amrex::get<2>(hv);
    ParallelDescriptor::ReduceRealMax(errmax);
    ParallelDescriptor::ReduceRealMax(refmax);
    ParallelDescriptor::ReduceLongSum(nbad);

    auto const error = errmax / refmax;
    amrex::Print() << "  relative error " << error
                   << ", non-finite values " << nbad << "\n";
    AMREX_ALWAYS_ASSERT(nbad == 0);
#ifdef AMREX_USE_FLOAT
    constexpr Real eps = 1.e-4;
#else
    constexpr Real eps = 1.e-12;
#endif
    AMREX_ALWAYS_ASSERT(error < eps);
}


#if (AMREX_SPACEDIM == 3)

// In twod_mode each z plane is an independent 2D convolution, so the Green's
// function must not depend on k. A domain whose z range does not start at zero
// is the interesting case here.
void test_twod_mode (Box const& domain, int max_grid_size)
{
    amrex::Print() << "\nTesting OpenBCSolver twod_mode on " << domain
                   << " against a direct convolution\n";

    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    DistributionMapping dm(ba);

    MultiFab rho(ba,dm,1,0), phi(ba,dm,1,0);
    for (MFIter mfi(rho); mfi.isValid(); ++mfi) {
        auto const& a = rho.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            a(i,j,k) = test_rhs(i,j,k);
        });
    }

    auto const lo = amrex::lbound(domain);
    FFT::Info info{};
    info.setTwoDMode(true);
    FFT::OpenBCSolver<Real> solver(domain, info);
    solver.setGreensFunction([=] AMREX_GPU_DEVICE (int i, int j, int) -> Real
    {
        return test_greens_function(i-lo.x, j-lo.y, 0);
    });
    solver.solve(phi, rho);

    auto const hi = amrex::ubound(domain);
    ReduceOps<ReduceOpMax, ReduceOpMax, ReduceOpSum> reduce_op;
    ReduceData<Real, Real, Long> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        auto const& a = phi.const_array(mfi);
        reduce_op.eval(mfi.validbox(), reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real exact = 0;
            for (int jj = lo.y; jj <= hi.y; ++jj) {
            for (int ii = lo.x; ii <= hi.x; ++ii) {
                exact += test_greens_function(amrex::Math::abs(i-ii),
                                              amrex::Math::abs(j-jj), 0) * test_rhs(ii,jj,k);
            }}
            return {amrex::Math::abs(a(i,j,k)-exact), amrex::Math::abs(exact),
                    Long(amrex::isnan(a(i,j,k)) || amrex::isinf(a(i,j,k)))};
        });
    }

    auto hv = reduce_data.value(reduce_op);
    auto errmax = amrex::get<0>(hv);
    auto refmax = amrex::get<1>(hv);
    auto nbad   = amrex::get<2>(hv);
    ParallelDescriptor::ReduceRealMax(errmax);
    ParallelDescriptor::ReduceRealMax(refmax);
    ParallelDescriptor::ReduceLongSum(nbad);

    auto const error = errmax / refmax;
    amrex::Print() << "  relative error " << error
                   << ", non-finite values " << nbad << "\n";
    AMREX_ALWAYS_ASSERT(nbad == 0);
#ifdef AMREX_USE_FLOAT
    constexpr Real eps = 1.e-4;
#else
    constexpr Real eps = 1.e-12;
#endif
    AMREX_ALWAYS_ASSERT(error < eps);
}

#endif

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        BL_PROFILE("main");

#if (AMREX_SPACEDIM == 3)
        int n_cell_x = 128;
        int n_cell_y = 128;
        int n_cell_z = 128;

        int max_grid_size_x = 32;
        int max_grid_size_y = 32;
        int max_grid_size_z = 32;

        {
            ParmParse pp;
            pp.query("n_cell_x", n_cell_x);
            pp.query("n_cell_y", n_cell_y);
            pp.query("n_cell_z", n_cell_z);
            pp.query("max_grid_size_x", max_grid_size_x);
            pp.query("max_grid_size_y", max_grid_size_y);
            pp.query("max_grid_size_z", max_grid_size_z);
        }

        Box domain(IntVect(0), IntVect(n_cell_x-1,n_cell_y-1,n_cell_z-1));
        BoxArray ba(domain);
        ba.maxSize(IntVect(max_grid_size_x, max_grid_size_y, max_grid_size_z));
        DistributionMapping dm(ba);

        Geometry geom(domain, RealBox(-1._rt, -1._rt, -1._rt, 1._rt, 1._rt, 1._rt),
                      CoordSys::cartesian, {AMREX_D_DECL(0,0,0)});

        auto const& dx = geom.CellSizeArray();

        std::array<IndexType,2> ixtypes{IndexType::TheCellType(),
                                        IndexType::TheNodeType()};
        for (auto const ixtype : ixtypes)
        {
            amrex::Print() << "\nTesting " << ixtype << "\n";

            BoxArray const& iba = amrex::convert(ba, ixtype);
            int ng = ixtype.cellCentered() ? 1 : 0;
            MultiFab rho(iba,dm,1,0);
            MultiFab phi(iba,dm,1,ng);
            phi.setVal(std::numeric_limits<Real>::max());

            fill_rhs(rho, geom, ixtype);

            FFT::PoissonOpenBC solver(geom, ixtype, IntVect(ng));
            solver.solve(phi, rho);

            Real mass = rho.sum_unique(0) * dx[0]*dx[1]*dx[2];
            Real offset = ixtype.cellCentered() ? 0.5_rt : 0.0_rt;
            auto x0 = -1._rt + offset*dx[0];
            auto y0 = -1._rt + offset*dx[1];
            auto z0 = -1._rt + offset*dx[2];
            auto r0 = std::sqrt(x0*x0+y0*y0+z0*z0); // radius of the corner cell
            auto expected = -mass/(4._rt*Math::pi<Real>()*r0);
            amrex::Print() << "  Expected phi: " << expected << "\n";

            int iextra = ixtype.cellCentered() ? 1 : 0;

            for (int k = 0; k < 2; ++k) {
            for (int j = 0; j < 2; ++j) {
            for (int i = 0; i < 2; ++i) {
                int ii = (i == 0) ? 0 : n_cell_x-iextra;
                int jj = (j == 0) ? 0 : n_cell_y-iextra;
                int kk = (k == 0) ? 0 : n_cell_z-iextra;
                IntVect corner(ii,jj,kk);
                auto v = amrex::get_cell_data(phi, corner);
                if (!v.empty()) {
                    amrex::AllPrint() << "  phi at " << corner << " is " << v[0] << "\n";
                    auto error = std::abs(expected-v[0])/std::max(std::abs(expected),std::abs(v[0]));
                    amrex::AllPrint() << "  error " << error << "\n";
#ifdef AMREX_USE_FLOAT
                    constexpr Real eps = Real(1.e-5);
#else
                    constexpr Real eps = 1.e-6;
#endif
                    AMREX_ALWAYS_ASSERT(error < eps);
                }
            }}}
        }

        {
            amrex::Print() << "\nTesting OpenBC padding against unpadded solve\n";

            AMREX_ALWAYS_ASSERT(FFT::Info{}.openbc_padding_nfactors ==
                                FFT::FastNumPrimeFactors());
            AMREX_ALWAYS_ASSERT(FFT::nextFastLen(13) ==
                                FFT::nextFastLen(13, FFT::FastNumPrimeFactors()));

            Box domain2(IntVect(0), IntVect(64,66,68));
            BoxArray ba2(domain2);
            ba2.maxSize(IntVect(32, 32, 32));
            DistributionMapping dm2(ba2);

            Geometry geom2(domain2,
                           RealBox(-1._rt, -1._rt, -1._rt, 1._rt, 1._rt, 1._rt),
                           CoordSys::cartesian, {AMREX_D_DECL(0,0,0)});

            MultiFab rho2(ba2, dm2, 1, 0);
            MultiFab phi_padded(ba2, dm2, 1, 0);
            MultiFab phi_unpadded(ba2, dm2, 1, 0);
            fill_rhs(rho2, geom2, IndexType::TheCellType());

            FFT::PoissonOpenBC padded_solver(geom2);
            FFT::Info unpadded_info;
            unpadded_info.setOpenBCPadding(false);
            FFT::PoissonOpenBC unpadded_solver(geom2, IndexType::TheCellType(),
                                               IntVect(0), unpadded_info);

            IntVect expected_padded_length = domain2.length();
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                expected_padded_length[idim] = FFT::nextFastLen(expected_padded_length[idim]);
            }
            AMREX_ALWAYS_ASSERT(padded_solver.PaddedLength() == expected_padded_length);
            AMREX_ALWAYS_ASSERT(unpadded_solver.PaddedLength() == domain2.length());

            padded_solver.solve(phi_padded, rho2);
            unpadded_solver.solve(phi_unpadded, rho2);

            MultiFab diff(ba2, dm2, 1, 0);
            MultiFab::Copy(diff, phi_padded, 0, 0, 1, 0);
            MultiFab::Subtract(diff, phi_unpadded, 0, 0, 1, 0);

            Real const refnorm = phi_unpadded.norm0(0);
            Real const error = diff.norm0(0) / refnorm;
            amrex::Print() << "  relative padded/unpadded error " << error << "\n";
#ifdef AMREX_USE_FLOAT
            constexpr Real eps = Real(1.e-5);
#else
            constexpr Real eps = 1.e-13;
#endif
            AMREX_ALWAYS_ASSERT(error < eps);
        }
#endif

        {
            int n_cell = 32;
            int max_grid_size = 16;
            ParmParse pp;
            pp.query("conv_n_cell", n_cell);
            pp.query("conv_max_grid_size", max_grid_size);

            test_convolution(Box(IntVect(0), IntVect(n_cell-1)), max_grid_size);
            // A domain that does not start at zero.
            test_convolution(Box(IntVect(-3), IntVect(n_cell-4)), max_grid_size);

#if (AMREX_SPACEDIM == 3)
            test_twod_mode(Box(IntVect(0), IntVect(n_cell-1)), max_grid_size);
            // A domain whose z range is entirely negative.
            test_twod_mode(Box(IntVect(0,0,-n_cell), IntVect(n_cell-1,n_cell-1,-1)),
                           max_grid_size);
#endif
        }
    }
    amrex::Finalize();
}
