#include <AMReX_FFT_Poisson.H> // Put this at the top for testing

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

namespace {

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
        Real x = (i+0.5_rt/nsub)*dx[0] + problo[0];
        Real y = (j+0.5_rt/nsub)*dx[1] + problo[1];
        Real z = (k+0.5_rt/nsub)*dx[2] + problo[2];
        if (ixtype.nodeCentered()) {
            x -= 0.5_rt*dx[0];
            y -= 0.5_rt*dx[1];
            z -= 0.5_rt*dx[2];
        }
        int n = 0;
        for (int isub = 0; isub < nsub; ++isub) {
        for (int jsub = 0; jsub < nsub; ++jsub) {
        for (int ksub = 0; ksub < nsub; ++ksub) {
            auto xs = x + isub*dxsub;
            auto ys = y + jsub*dysub;
            auto zs = z + ksub*dzsub;
            if ((xs*xs+ys*ys+zs*zs) < 0.25_rt) { ++n; }
        }}}
        rhoma[b](i,j,k) = Real(n) / Real(nsub*nsub*nsub);
    });
}

}

int main (int argc, char* argv[])
{
    static_assert(AMREX_SPACEDIM == 3);

    amrex::Initialize(argc, argv);
    {
        BL_PROFILE("main");

        int n_cell_x = 128;
        int n_cell_y = 128;
        int n_cell_z = 128;

        int max_grid_size_x = 32;
        int max_grid_size_y = 32;
        int max_grid_size_z = 32;

        ParmParse pp;
        pp.query("n_cell_x", n_cell_x);
        pp.query("n_cell_y", n_cell_y);
        pp.query("n_cell_z", n_cell_z);
        pp.query("max_grid_size_x", max_grid_size_x);
        pp.query("max_grid_size_y", max_grid_size_y);
        pp.query("max_grid_size_z", max_grid_size_z);

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
                    constexpr Real eps = 1.e-5;
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
            constexpr Real eps = 1.e-5;
#else
            constexpr Real eps = 1.e-13;
#endif
            AMREX_ALWAYS_ASSERT(error < eps);
        }


        // Vico-Greengard-Ferrando Green's function, in Fourier space, for a
        // Gaussian source with a closed-form potential: (a) a nearly cubic
        // domain with a resolved source, where the spectral accuracy shows,
        // and (b) an elongated domain with a shifted origin, where the
        // per-direction oversampling rule matters and the fixed 2*N grid
        // of the literature aliases.
        // The grid-comparison tolerances come from the periodic images of
        // the kernel's ringing tail beyond L, which enter through the source
        // at the domain boundary: ~3e-9 of its peak in (a), ~4e-6 in (b).
        struct VGCase {
            Box domain; RealBox rb; Real sigma; Real tol_exact; Real tol_grid; bool check_fixed_2N;
        };
        VGCase const vg_cases[] = {
            {Box(IntVect(0), IntVect(63,47,79)),
             RealBox(-1._rt, -0.75_rt, -1.25_rt, 1._rt, 0.75_rt, 1.25_rt),
#ifdef AMREX_USE_FLOAT
             0.12_rt, 1.e-4_rt, 1.e-5_rt, false},
#else
             0.12_rt, 1.e-9_rt, 1.e-12_rt, false},
#endif
            {Box(IntVect(5,-3,7), IntVect(68,12,22)),
             RealBox(-1._rt, -0.25_rt, -0.25_rt, 1._rt, 0.25_rt, 0.25_rt),
#ifdef AMREX_USE_FLOAT
             0.05_rt, 1.e-4_rt, 1.e-5_rt, true}
#else
             // source not fully resolved: truncation error ~ 3e-7
             0.05_rt, 1.e-5_rt, 1.e-8_rt, true}
#endif
        };
        for (auto const& vgc : vg_cases) {
        for (auto const ixtype : ixtypes)
        {
            amrex::Print() << "\nTesting Vico-Greengard-Ferrando Green's function " << ixtype
                           << " on " << vgc.domain << "\n";

            // Boxes chosen so that the slabs of the solver's doubled domain
            // and the pencils of the DCT differ from the user's layout.
            Box domain3 = vgc.domain;
            BoxArray ba3(domain3);
            ba3.maxSize(IntVect(32, 24, 40));
            DistributionMapping dm3(ba3);
            domain3.convert(ixtype);
            ba3.convert(ixtype);
            Real const offset = ixtype.cellCentered() ? 0.5_rt : 0.0_rt;

            Geometry geom3(amrex::enclosedCells(domain3), vgc.rb,
                           CoordSys::cartesian, {AMREX_D_DECL(0,0,0)});
            auto const& dx3 = geom3.CellSizeArray();
            auto const& problo3 = geom3.ProbLoArray();
            auto const& dlo = domain3.smallEnd();

            Real const sigma = vgc.sigma;
            MultiFab rho3(ba3, dm3, 1, 0);
            MultiFab phi3(ba3, dm3, 1, 0);
            MultiFab phi_exact(ba3, dm3, 1, 0);
            {
                auto const& rhoma = rho3.arrays();
                auto const& phima = phi_exact.arrays();
                Real const q = std::pow(2._rt*Math::pi<Real>(), 1.5_rt) * sigma*sigma*sigma;
                ParallelFor(rho3, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
                {
                    Real x = (i-dlo[0]+offset)*dx3[0] + problo3[0];
                    Real y = (j-dlo[1]+offset)*dx3[1] + problo3[1];
                    Real z = (k-dlo[2]+offset)*dx3[2] + problo3[2];
                    Real r = std::sqrt(x*x+y*y+z*z);
                    rhoma[b](i,j,k) = std::exp(-0.5_rt*r*r/(sigma*sigma));
                    // Solution of laplace(phi) = rho
                    if (r > 0._rt) {
                        phima[b](i,j,k) = -q/(4._rt*Math::pi<Real>()*r)
                            * std::erf(r/(std::sqrt(2._rt)*sigma));
                    } else {
                        phima[b](i,j,k) = -q/(4._rt*Math::pi<Real>())
                            * std::sqrt(2._rt/Math::pi<Real>())/sigma;
                    }
                });
            }
            Real const refnorm = phi_exact.norm0(0);

            // setGreensFunction(MF) must reproduce the functor version. The
            // octant MultiFab is cell-centered and larger than the octant on
            // purpose: any index type and any covering layout must be accepted.
            {
                FFT::OpenBCSolver<Real> mf_solver(domain3);
                Box const octant(domain3.smallEnd(),
                                 domain3.smallEnd()+mf_solver.PaddedLength()-IntVect(1));
                BoxArray oba(Box(octant.smallEnd()-IntVect(2), octant.bigEnd()+IntVect(31,39,23)));
                oba.maxSize(IntVect(40, 64, 48));
                DistributionMapping odm(oba);
                MultiFab goct(oba, odm, 1, 0);
                auto const gfac = 1._rt/std::sqrt(12._rt);
                auto const gpfac = -0.125_rt * (dx3[0]*dx3[1]*dx3[2]) / (4._rt*Math::pi<Real>());
                auto const& gma = goct.arrays();
                ParallelFor(goct, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
                {
                    auto x = (Real(i-dlo[0]) - gfac) * dx3[0];
                    auto y = (Real(j-dlo[1]) - gfac) * dx3[1];
                    auto z = (Real(k-dlo[2]) - gfac) * dx3[2];
                    Real r = 0;
                    for (int gx = 0; gx < 2; ++gx) {
                    for (int gy = 0; gy < 2; ++gy) {
                    for (int gz = 0; gz < 2; ++gz) {
                        auto xg = x + 2*gx*gfac*dx3[0];
                        auto yg = y + 2*gy*gfac*dx3[1];
                        auto zg = z + 2*gz*gfac*dx3[2];
                        r += 1._rt/std::sqrt(xg*xg+yg*yg+zg*zg);
                    }}}
                    gma[b](i,j,k) = gpfac * r;
                });
                mf_solver.setGreensFunction(goct);
                MultiFab phi_mf(ba3, dm3, 1, 0);
                mf_solver.solve(phi_mf, rho3);

                FFT::PoissonOpenBC ref_solver(geom3, ixtype, IntVect(0));
                MultiFab phi_ref(ba3, dm3, 1, 0);
                ref_solver.solve(phi_ref, rho3);
                MultiFab::Subtract(phi_ref, phi_exact, 0, 0, 1, 0);
                amrex::Print() << "  PoissonOpenBC relative error "
                               << phi_ref.norm0(0) / refnorm << "\n";
                MultiFab::Add(phi_ref, phi_exact, 0, 0, 1, 0);

                MultiFab::Subtract(phi_mf, phi_ref, 0, 0, 1, 0);
                Real const diff = phi_mf.norm0(0) / phi_ref.norm0(0);
                amrex::Print() << "  setGreensFunction(MF) vs functor: " << diff << "\n";
#ifdef AMREX_USE_FLOAT
                AMREX_ALWAYS_ASSERT(diff < 1.e-5);
#else
                AMREX_ALWAYS_ASSERT(diff < 1.e-13);
#endif
            }

            // Truncation radius: larger than the largest distance between
            // two points of the domain.
            Real const L = 1.1_rt * std::sqrt(AMREX_D_TERM(  std::pow(domain3.length(0)*dx3[0],2),
                                                           + std::pow(domain3.length(1)*dx3[1],2),
                                                           + std::pow(domain3.length(2)*dx3[2],2)));
            // Fourier transform of -1/(4 pi r) truncated at radius L:
            // -2 sin^2(L|k|/2)/|k|^2, with limit -L^2/2 at k = 0, sampled at
            // k = pi*n/(K*dx) for a half period K. One solver is reused
            // across several K, as a caller with a changing mesh would.
            // In case (b) the FFTs are restricted to one process, so that
            // under MPI the other processes have no local data.
            FFT::Info vg_info{};
            if (vgc.check_fixed_2N) { vg_info.setNumProcs(1); }
            FFT::OpenBCSolver<Real> solver(domain3, vg_info);
            auto solve_vg = [&] (IntVect const& Kh, MultiFab& phi_out)
            {
                GpuArray<Real,3> dk{Math::pi<Real>()/(Kh[0]*dx3[0]),
                                    Math::pi<Real>()/(Kh[1]*dx3[1]),
                                    Math::pi<Real>()/(Kh[2]*dx3[2])};
                solver.setGreensFunctionFromSpectral(
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) -> Real
                {
                    Real const kx = i*dk[0];
                    Real const ky = j*dk[1];
                    Real const kz = k*dk[2];
                    Real const k2 = kx*kx + ky*ky + kz*kz;
                    if (k2 == 0._rt) { return -0.5_rt*L*L; }
                    Real const s = std::sin(0.5_rt*L*std::sqrt(k2));
                    return -2._rt*s*s/k2;
                }, Kh);
                solver.solve(phi_out, rho3);
            };

            // Half period of the oversampled grid, per direction: the
            // periodic images of the truncated kernel must not reach the
            // lags used, which are up to length-1 cells, i.e.
            // 2K > length-1 + L/dx, plus two cells of margin for the ringing
            // of the band-limited kernel beyond L, and K+1 must cover the
            // octant the solver reads.
            IntVect K;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                int const kmin = static_cast<int>(std::floor(
                    (Real(domain3.length(idim)-1) + L/dx3[idim]) / 2._rt)) + 2;
                K[idim] = FFT::nextFastLen(std::max(kmin, solver.PaddedLength()[idim]-1));
            }
            amrex::Print() << "  truncation radius " << L
                           << ", DCT-I grid " << K+IntVect(1) << "\n";

            solve_vg(K, phi3);
            MultiFab::Subtract(phi3, phi_exact, 0, 0, 1, 0);
            Real const vg_error = phi3.norm0(0) / refnorm;
            amrex::Print() << "  Vico-Greengard relative error " << vg_error << "\n";
            AMREX_ALWAYS_ASSERT(vg_error < vgc.tol_exact);

            // Setting the same Green's function again must give the same
            // result bit for bit.
            {
                MultiFab phi_again(ba3, dm3, 1, 0);
                solve_vg(K, phi_again);
                MultiFab::Subtract(phi_again, phi_exact, 0, 0, 1, 0);
                MultiFab::Subtract(phi_again, phi3, 0, 0, 1, 0);
                AMREX_ALWAYS_ASSERT(phi_again.norm0(0) == 0._rt);
            }

            // A larger grid must give the same result to round-off: it
            // differs only in how far the periodic images of the truncated
            // kernel are kept from the lags used.
            {
                MultiFab phi_big(ba3, dm3, 1, 0);
                solve_vg(K+IntVect(16), phi_big);
                MultiFab::Subtract(phi_big, phi_exact, 0, 0, 1, 0);
                MultiFab::Subtract(phi_big, phi3, 0, 0, 1, 0);
                Real const diff = phi_big.norm0(0) / refnorm;
                amrex::Print() << "  relative difference to the K+16 grid " << diff << "\n";
                AMREX_ALWAYS_ASSERT(diff < vgc.tol_grid);
            }

            // The fixed grid of the literature, 4N cells i.e. K = 2N, is
            // insufficient in the short directions of an elongated domain
            // and aliases.
            if (vgc.check_fixed_2N) {
                MultiFab phi_2N(ba3, dm3, 1, 0);
                solve_vg(2*domain3.length(), phi_2N);
                MultiFab::Subtract(phi_2N, phi_exact, 0, 0, 1, 0);
                MultiFab::Subtract(phi_2N, phi3, 0, 0, 1, 0);
                Real const diff = phi_2N.norm0(0) / refnorm;
                amrex::Print() << "  relative difference to the fixed 2*N grid " << diff
                               << " (expected to be large)\n";
                AMREX_ALWAYS_ASSERT(diff > 1.e-2);
            }
        }}
    }
    amrex::Finalize();
}
