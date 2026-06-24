#include <AMReX_FFT_Stokes.H>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        BL_PROFILE("main");

        AMREX_D_TERM(int n_cell_x = 128;,
                     int n_cell_y = 32;,
                     int n_cell_z = 64);

        AMREX_D_TERM(int max_grid_size_x = 16;,
                     int max_grid_size_y = 16;,
                     int max_grid_size_z = 16);

        AMREX_D_TERM(Real prob_lo_x = 0.;,
                     Real prob_lo_y = 0.;,
                     Real prob_lo_z = 0.);
        AMREX_D_TERM(Real prob_hi_x = 1.;,
                     Real prob_hi_y = 1.;,
                     Real prob_hi_z = 1.);

        {
            ParmParse pp;
            AMREX_D_TERM(pp.query("n_cell_x", n_cell_x);,
                         pp.query("n_cell_y", n_cell_y);,
                         pp.query("n_cell_z", n_cell_z));
            AMREX_D_TERM(pp.query("max_grid_size_x", max_grid_size_x);,
                         pp.query("max_grid_size_y", max_grid_size_y);,
                         pp.query("max_grid_size_z", max_grid_size_z));
        }

        auto bcs = std::pair<FFT::Boundary,FFT::Boundary>{FFT::Boundary::periodic,
                                                          FFT::Boundary::periodic};

        Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM>
                fft_bc{AMREX_D_DECL(bcs,bcs,bcs)};

        Box domain(IntVect(0),IntVect(AMREX_D_DECL(n_cell_x-1,n_cell_y-1,n_cell_z-1)));
        BoxArray ba(domain);
        ba.maxSize(IntVect(AMREX_D_DECL(max_grid_size_x,
                                        max_grid_size_y,
                                        max_grid_size_z)));
        DistributionMapping dm(ba);

        Geometry geom;
        {
            geom.define(domain,
                        RealBox(AMREX_D_DECL(prob_lo_x,prob_lo_y,prob_lo_z),
                                AMREX_D_DECL(prob_hi_x,prob_hi_y,prob_hi_z)),
                        CoordSys::cartesian, {AMREX_D_DECL(1,1,1)});
        }
        auto const& dx = geom.CellSizeArray();

        MultiFab pres(ba,dm,1,1);
        auto const& p = pres.arrays();
        ParallelFor(pres, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            AMREX_D_TERM(Real x = (i+0.5_rt) * dx[0] - 0.5_rt;,
                         Real y = (j+0.5_rt) * dx[1] - 0.5_rt;,
                         Real z = (k+0.5_rt) * dx[2] - 0.5_rt);
            p[b](i,j,k) = std::cos(Real(2.)*Math::pi<Real>()*x)
                         + std::cos(Real(2.)*Math::pi<Real>()*y);
#if (BL_SPACEDIM == 3)
            p[b](i,j,k) += std::sin(Real(2.)*Math::pi<Real>()*z);
#endif
        });

        pres.FillBoundary(geom.periodicity());

        MultiFab pres_return(ba,dm,1,0);

        BoxArray fbau = amrex::convert(ba, IntVect::TheDimensionVector(0));
        BoxArray fbav = amrex::convert(ba, IntVect::TheDimensionVector(1));

        MultiFab u_mf(fbau, dm, 1, 1);
        MultiFab u_mf2(fbau, dm, 1, 1);
        MultiFab rhsx(fbau,dm,1,0);
        MultiFab v_mf(fbav, dm, 1, 1);
        MultiFab v_mf2(fbav, dm, 1, 1);
        MultiFab rhsy(fbav,dm,1,0);

#if (BL_SPACEDIM == 3)

        BoxArray fbaw = amrex::convert(ba, IntVect::TheDimensionVector(2));

        MultiFab w_mf(fbaw, dm, 1, 1);
        MultiFab w_mf2(fbaw, dm, 1, 1);
        MultiFab rhsz(fbaw,dm,1,0);

#endif
        AMREX_ALWAYS_ASSERT(pres.boxArray().size() == u_mf.boxArray().size());
        AMREX_ALWAYS_ASSERT(pres.boxArray().size() == v_mf.boxArray().size());
#if (BL_SPACEDIM == 3)
        AMREX_ALWAYS_ASSERT(pres.boxArray().size() == w_mf.boxArray().size());
#endif
        auto const& face_u = u_mf.arrays();
        ParallelFor(u_mf, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            AMREX_D_TERM(int ii = (i == n_cell_x) ? 0 : i;,
                         int jj = (j == n_cell_y) ? 0 : j;,
                         int kk = (k == n_cell_z) ? 0 : k);
            AMREX_D_TERM(Real x = Real(ii) * dx[0] - 0.5_rt;,
                         Real y = Real(jj) * dx[1] - 0.5_rt;,
                         Real z = Real(kk) * dx[2] - 0.5_rt);
            face_u[b](i,j,k) = std::exp(-10._rt*
                        (AMREX_D_TERM(x*x*1.05_rt, + y*y*0.90_rt, + z*z)));
#if (BL_SPACEDIM == 3)
            z = Real(kk) * dx[2];
            face_u[b](i,j,k) = std::cos(Real(2.)*Math::pi<Real>()*z);
#else
            y = Real(jj) * dx[1];
            face_u[b](i,j,k) = std::cos(Real(2.)*Math::pi<Real>()*y);
#endif
        });

        u_mf.FillBoundary(geom.periodicity());

        auto const& face_v = v_mf.arrays();
        ParallelFor(v_mf, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            AMREX_D_TERM(int ii = (i == n_cell_x) ? 0 : i;,
                         int jj = (j == n_cell_y) ? 0 : j;,
                         int kk = (k == n_cell_z) ? 0 : k);
            AMREX_D_TERM(Real x = Real(ii) * dx[0] - 0.5_rt;,
                         Real y = Real(jj) * dx[1] - 0.5_rt;,
                         Real z = Real(kk) * dx[2] - 0.5_rt);
            face_v[b](i,j,k) = std::exp(-10._rt*
                        (AMREX_D_TERM(x*x*1.05_rt, + y*y*0.90_rt, + z*z)));
            face_v[b](i,j,k) = 0.;
            face_v[b](i,j,k) = std::cos(Real(2.)*Math::pi<Real>()*x);
        });

        v_mf.FillBoundary(geom.periodicity());

#if (BL_SPACEDIM == 3)
        auto const& face_w = w_mf.arrays();
        ParallelFor(w_mf, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            AMREX_D_TERM(int ii = (i == n_cell_x) ? 0 : i;,
                         int jj = (j == n_cell_y) ? 0 : j;,
                         int kk = (k == n_cell_z) ? 0 : k);
            AMREX_D_TERM(Real x = Real(ii) * dx[0] - 0.5_rt;,
                         Real y = Real(jj) * dx[1] - 0.5_rt;,
                         Real z = Real(kk) * dx[2] - 0.5_rt);
            face_w[b](i,j,k) = std::exp(-10._rt*
                        (AMREX_D_TERM(x*x*1.05_rt, + y*y*0.90_rt, + z*z)));
            face_w[b](i,j,k) = 0.;
            face_w[b](i,j,k) = std::cos(Real(2.)*Math::pi<Real>()*x) +  std::sin(2+Real(2.)*Math::pi<Real>()*y)  ;
        });

        w_mf.FillBoundary(geom.periodicity());
#endif

        auto alpha = Real(.01);
        auto eta = Real(.01);

        {
            auto const& pres_arr = pres.const_arrays();
            auto const& u = u_mf.const_arrays();
            auto const& rx = rhsx.arrays();

            ParallelFor(u_mf, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept {
                Real lap = (u[b](i-1,j,k)-Real(2.)*u[b](i,j,k)+u[b](i+1,j,k))/(dx[0]*dx[0])
                         + (u[b](i,j-1,k)-Real(2.)*u[b](i,j,k)+u[b](i,j+1,k))/(dx[1]*dx[1]);
#if (BL_SPACEDIM == 3)
                lap += (u[b](i,j,k-1)-Real(2.)*u[b](i,j,k)+u[b](i,j,k+1))/(dx[2]*dx[2]);
#endif
                rx[b](i,j,k) = alpha*u[b](i,j,k) - eta*lap
                             + (pres_arr[b](i,j,k)-pres_arr[b](i-1,j,k))/dx[0];
            });
        }

        {
            auto const& pres_arr = pres.const_arrays();
            auto const& v = v_mf.const_arrays();
            auto const& ry = rhsy.arrays();

            ParallelFor(v_mf, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept {
                Real lap = (v[b](i-1,j,k)-Real(2.)*v[b](i,j,k)+v[b](i+1,j,k))/(dx[0]*dx[0])
                         + (v[b](i,j-1,k)-Real(2.)*v[b](i,j,k)+v[b](i,j+1,k))/(dx[1]*dx[1]);
#if (BL_SPACEDIM == 3)
                lap += (v[b](i,j,k-1)-Real(2.)*v[b](i,j,k)+v[b](i,j,k+1))/(dx[2]*dx[2]);
#endif
                ry[b](i,j,k) = alpha*v[b](i,j,k) - eta*lap
                             + (pres_arr[b](i,j,k)-pres_arr[b](i,j-1,k))/dx[1];
            });
        }
#if (BL_SPACEDIM == 3)

        {
            auto const& pres_arr = pres.const_arrays();
            auto const& w = w_mf.const_arrays();
            auto const& rz = rhsz.arrays();

            ParallelFor(w_mf, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept {
                Real lap = (w[b](i-1,j,k)-Real(2.)*w[b](i,j,k)+w[b](i+1,j,k))/(dx[0]*dx[0])
                         + (w[b](i,j-1,k)-Real(2.)*w[b](i,j,k)+w[b](i,j+1,k))/(dx[1]*dx[1])
                         + (w[b](i,j,k-1)-Real(2.)*w[b](i,j,k)+w[b](i,j,k+1))/(dx[2]*dx[2]);
                rz[b](i,j,k) = alpha*w[b](i,j,k) - eta*lap
                            + (pres_arr[b](i,j,k)-pres_arr[b](i,j,k-1))/dx[2];
            });
        }
#endif
        {
            FFT::Stokes stokes(geom, fft_bc);
            stokes.solve(AMREX_D_DECL(u_mf2, v_mf2, w_mf2), pres_return,
                         AMREX_D_DECL(rhsx, rhsy, rhsz), alpha, eta);

            MultiFab u_err(u_mf.boxArray(), dm, 1, 0);
            MultiFab v_err(v_mf.boxArray(), dm, 1, 0);
            u_err.ParallelCopy(u_mf);
            v_err.ParallelCopy(v_mf);
            MultiFab::Subtract(u_err, u_mf2, 0, 0, 1, 0);
            MultiFab::Subtract(v_err, v_mf2, 0, 0, 1, 0);

            Real u_error = u_err.norminf();
            Real v_error = v_err.norminf();
            amrex::Print() << "  U error expected to be close to zero: " << u_error << "\n";
            amrex::Print() << "  V error expected to be close to zero: " << v_error << "\n";
#ifdef AMREX_USE_FLOAT
            auto eps = 1.e-3f;
#else
            auto eps = 1.e-12;
#endif
            AMREX_ALWAYS_ASSERT(u_error < eps);
            AMREX_ALWAYS_ASSERT(v_error < eps);

#if (BL_SPACEDIM == 3)
            MultiFab w_err(w_mf.boxArray(), dm, 1, 0);
            w_err.ParallelCopy(w_mf);
            MultiFab::Subtract(w_err, w_mf2, 0, 0, 1, 0);
            Real w_error = w_err.norminf();
            amrex::Print() << "  W error expected to be close to zero: " << w_error << "\n";
            AMREX_ALWAYS_ASSERT(w_error < eps);
#endif

            MultiFab p_err(ba, dm, 1, 0);
            p_err.ParallelCopy(pres);
            MultiFab::Subtract(p_err, pres_return, 0, 0, 1, 0);
            Real mean_ref = pres.sum(0, false) / Real(geom.Domain().d_numPts());
            Real mean_sol = pres_return.sum(0, false) / Real(geom.Domain().d_numPts());
            p_err.plus(mean_sol - mean_ref, 0, 1, 0);
            Real p_error = p_err.norminf();
            amrex::Print() << "  Pressure error (mean-adjusted) expected to be close to zero: "
                           << p_error << "\n";
            AMREX_ALWAYS_ASSERT(p_error < eps);
        }
    }
    amrex::Finalize();
}
