#include <cuda_device_runtime_api.h>

#include <iostream>
#include <AMReX.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_BaseFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

using namespace amrex;

int main (int argc, char* argv[])
{
    std::cout << "**********************************\n";

    amrex::Initialize(argc, argv);
    {
        int n_cell, max_grid_size, nsteps, plot_int;
        Vector<int> is_periodic(AMREX_SPACEDIM,1);  // periodic in all direction by default

        {
            ParmParse pp;

            pp.get("n_cell",n_cell);
            pp.get("max_grid_size",max_grid_size);
            pp.queryarr("is_periodic", is_periodic);

            plot_int = -1;
            pp.query("plot_int",plot_int);

            nsteps = 10;
            pp.query("nsteps",nsteps);
        }

        // make BoxArray, Geometry & DistributionMapping
        BoxArray ba;
        Geometry geom;
        {
            IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
            IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
            Box domain(dom_lo, dom_hi);
            ba.define(domain);
            ba.maxSize(max_grid_size);

            RealBox real_box({AMREX_D_DECL(-1.0,-1.0,-1.0)},
                             {AMREX_D_DECL( 1.0, 1.0, 1.0)});

            geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
        }

        int Nghost = 1;
        int Ncomp  = 1;
        DistributionMapping dm(ba);

        amrex::Print() << std::endl << std::endl;

        {
            BL_PROFILE("GPU Tests");
            Gpu::LaunchSafeGuard lsg(true);
            amrex::Print() << "GPU Tests" << std::endl;
            MultiFab x(ba, dm, Ncomp, Nghost);
            MultiFab y(ba, dm, Ncomp, Nghost);

            x.setVal(1701.4);
            y.setVal(74.656);

            x.setVal(2.0);
            Real selfdot_begin = amrex::second();
            double selfdot = MultiFab::Dot(x, 0, Ncomp, 0);
            Real selfdot_end = amrex::second();

            Real sum_begin = amrex::second();
            double sum = x.sum(); 
            Real sum_end = amrex::second();

            Real dot_begin = amrex::second();
            double dot = MultiFab::Dot(x, 0, y, 0, Ncomp, Nghost); 
            Real dot_end = amrex::second();

            Real norm0_begin = amrex::second();
            double norm0 = x.norm0(); 
            Real norm0_end = amrex::second();

            Real norm1_begin = amrex::second();
            double norm1 = x.norm1(); 
            Real norm1_end = amrex::second();

            Real norm2_begin = amrex::second();
            double norm2 = x.norm2(); 
            Real norm2_end = amrex::second();

            amrex::Print() << "Self Dot: " << selfdot << ". " << selfdot_end - selfdot_begin << " secs." << std::endl;
            amrex::Print() << "Sum: " << sum << ". " << sum_end - sum_begin << " secs." << std::endl;
            amrex::Print() << "Dot: " << dot << ". " << dot_end - dot_begin << " secs." << std::endl;
            amrex::Print() << "Norm0: " << norm0 << ". " << norm0_end - norm0_begin << " secs." << std::endl;
            amrex::Print() << "Norm1: " << norm1 << ". " << norm1_end - norm1_begin << " secs." << std::endl;
            amrex::Print() << "Norm2: " << norm2 << ". " << norm2_end - norm2_begin << " secs." << std::endl;
        }

        amrex::Print() << std::endl << std::endl;

        {
            BL_PROFILE("GPU Tests B");
            Gpu::LaunchSafeGuard lsg(true);
            amrex::Print() << "GPU Tests 2" << std::endl;
            MultiFab x(ba, dm, Ncomp, Nghost);
            MultiFab y(ba, dm, Ncomp, Nghost);

            x.setVal(1701.4);
            y.setVal(74.656);

            Real sum_begin = amrex::second();
            double sum = x.sum(); 
            Real sum_end = amrex::second();

            Real dot_begin = amrex::second();
            double dot = MultiFab::Dot(x, 0, y, 0, Ncomp, Nghost); 
            Real dot_end = amrex::second();

            Real norm0_begin = amrex::second();
            double norm0 = x.norm0(); 
            Real norm0_end = amrex::second();

            Real norm1_begin = amrex::second();
            double norm1 = x.norm1(); 
            Real norm1_end = amrex::second();

            Real norm2_begin = amrex::second();
            double norm2 = x.norm2(); 
            Real norm2_end = amrex::second();

            amrex::Print() << "Sum: " << sum << ". " << sum_end - sum_begin << " secs." << std::endl;
            amrex::Print() << "Dot: " << dot << ". " << dot_end - dot_begin << " secs." << std::endl;
            amrex::Print() << "Norm0: " << norm0 << ". " << norm0_end - norm0_begin << " secs." << std::endl;
            amrex::Print() << "Norm1: " << norm1 << ". " << norm1_end - norm1_begin << " secs." << std::endl;
            amrex::Print() << "Norm2: " << norm2 << ". " << norm2_end - norm2_begin << " secs." << std::endl;
        }

        amrex::Print() << std::endl << std::endl;

        {
            BL_PROFILE("CPU Tests");
            Gpu::LaunchSafeGuard lsg(false);
            amrex::Print() << "CPU Tests" << std::endl;
            MultiFab x(ba, dm, Ncomp, Nghost);
            MultiFab y(ba, dm, Ncomp, Nghost);

            x.setVal(1701.4);
            y.setVal(74.656);

            Real sum_begin = amrex::second();
            double sum = x.sum(); 
            Real sum_end = amrex::second();

            Real dot_begin = amrex::second();
            double dot = MultiFab::Dot(x, 0, y, 0, Ncomp, Nghost); 
            Real dot_end = amrex::second();

            Real norm0_begin = amrex::second();
            double norm0 = x.norm0(); 
            Real norm0_end = amrex::second();

            Real norm1_begin = amrex::second();
            double norm1 = x.norm1(); 
            Real norm1_end = amrex::second();

            Real norm2_begin = amrex::second();
            double norm2 = x.norm2(); 
            Real norm2_end = amrex::second();

            amrex::Print() << "Sum: " << sum << ". " << sum_end - sum_begin << " secs." << std::endl;
            amrex::Print() << "Dot: " << dot << ". " << dot_end - dot_begin << " secs." << std::endl;
            amrex::Print() << "Norm0: " << norm0 << ". " << norm0_end - norm0_begin << " secs." << std::endl;
            amrex::Print() << "Norm1: " << norm1 << ". " << norm1_end - norm1_begin << " secs." << std::endl;
            amrex::Print() << "Norm2: " << norm2 << ". " << norm2_end - norm2_begin << " secs." << std::endl;
        }


    }
    amrex::Finalize();
}
