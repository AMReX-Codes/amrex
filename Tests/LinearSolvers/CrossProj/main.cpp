#include <AMReX_CrossProj.H>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        BL_PROFILE("main");

        int n_cell_x = 64;
        int n_cell_y = 64;
        int max_grid_size = 16;
        int max_grid_size_2 = 32; // second layout for the layout check
        Real prob_hi_y = 1.0;     // square cells by default
        int verbose = 0;
        {
            ParmParse pp;
            pp.query("n_cell_x", n_cell_x);
            pp.query("n_cell_y", n_cell_y);
            pp.query("max_grid_size", max_grid_size);
            pp.query("max_grid_size_2", max_grid_size_2);
            pp.query("prob_hi_y", prob_hi_y);
            pp.query("verbose", verbose);
        }

        Box domain(IntVect(0), IntVect(n_cell_x-1, n_cell_y-1));
        Geometry geom(domain, RealBox({0.,0.}, {1.,prob_hi_y}), CoordSys::cartesian, {1,1});
        auto const& dx = geom.CellSizeArray();
        GpuArray<Real,AMREX_SPACEDIM> const plen{geom.ProbLength(0), geom.ProbLength(1)};

        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);

        MultiFab vel(ba, dm, AMREX_SPACEDIM, 1);
        auto const& u = vel.arrays();
        ParallelFor(vel, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            constexpr Real tpi = Real(2.)*Math::pi<Real>();
            Real x = (Real(i)+0.5_rt) * dx[0] / plen[0] - 0.5_rt;
            Real y = (Real(j)+0.5_rt) * dx[1] / plen[1] - 0.5_rt;
            // Plus a checkerboard, which the cross stencil barely sees.
            Real const chk = ((i+j) % 2 == 0) ? 0.1_rt : -0.1_rt;
            u[b](i,j,k,0) = std::cos(tpi*x) + std::cos(tpi*y) + chk;
            u[b](i,j,k,1) = std::cos(tpi*x+0.2_rt) + std::cos(tpi*y+0.7_rt);
        });
        vel.FillBoundary(geom.periodicity());

        // Node-centered cross-stencil divergence at (i+1/2,j+1/2)
        auto nodediv = [=] AMREX_GPU_DEVICE (Array4<Real const> const& v,
                                             int i, int j, int k) -> Real
        {
            return (v(i+1,j+1,k,0) + v(i+1,j,k,0) - v(i,j+1,k,0) - v(i,j,k,0))
                * (0.5_rt/dx[0])
                +  (v(i+1,j+1,k,1) + v(i,j+1,k,1) - v(i+1,j,k,1) - v(i,j,k,1))
                * (0.5_rt/dx[1]);
        };

        GMG::CrossProj crossproj(geom);
        crossproj.setVerbose(verbose);

        MultiFab vel_return(ba, dm, AMREX_SPACEDIM, 1);
        crossproj.solve(vel_return, vel); // also fills ghost cells
        amrex::Print() << "  iterations " << crossproj.getNumIters() << "\n";

        MultiFab div_orig(ba, dm, 1, 0);
        MultiFab div_err(ba, dm, 1, 0);
        {
            auto const& ud = vel_return.const_arrays();
            auto const& uo = vel.const_arrays();
            auto const& divinit = div_orig.arrays();
            auto const& div = div_err.arrays();
            ParallelFor(div_err, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                divinit[b](i,j,k) = nodediv(uo[b], i, j, k);
                div[b](i,j,k) = nodediv(ud[b], i, j, k);
            });
        }
        Real const div_initial = div_orig.norminf();
        Real const div_error = div_err.norminf();
        amrex::Print() << "  original divergence " << div_initial << "\n"
                       << "  divergence expected to be close to zero: " << div_error << "\n";

        // A second projection should not change the result.
        MultiFab vel_twice(ba, dm, AMREX_SPACEDIM, 0);
        crossproj.solve(vel_twice, vel_return);
        MultiFab::Subtract(vel_twice, vel_return, 0, 0, AMREX_SPACEDIM, 0);
        Real const idem_error = vel_twice.norminf(0, AMREX_SPACEDIM, IntVect(0));
        Real const vel_max = vel_return.norminf(0, AMREX_SPACEDIM, IntVect(0));
        amrex::Print() << "  idempotency error: " << idem_error << "\n";

        // The result does not depend on the box layout.
        BoxArray ba2(domain);
        ba2.maxSize(max_grid_size_2);
        DistributionMapping dm2(ba2);
        MultiFab vel2(ba2, dm2, AMREX_SPACEDIM, 0);
        vel2.ParallelCopy(vel, 0, 0, AMREX_SPACEDIM);
        GMG::CrossProj crossproj2(geom);
        crossproj2.solve(vel2, vel2);
        MultiFab vel_return2(ba, dm, AMREX_SPACEDIM, 0);
        vel_return2.ParallelCopy(vel2, 0, 0, AMREX_SPACEDIM);
        MultiFab::Subtract(vel_return2, vel_return, 0, 0, AMREX_SPACEDIM, 0);
        Real const layout_error = vel_return2.norminf(0, AMREX_SPACEDIM, IntVect(0));
        amrex::Print() << "  layout difference: " << layout_error
                       << " (iterations " << crossproj2.getNumIters() << ")\n";

        Real const tol = std::is_same_v<Real,float> ? 1.e-3_rt : 1.e-8_rt;
        AMREX_ALWAYS_ASSERT(div_error <= tol*div_initial);
        AMREX_ALWAYS_ASSERT(idem_error <= tol*vel_max);
        AMREX_ALWAYS_ASSERT(layout_error <= tol*vel_max);
    }
    amrex::Finalize();
}
