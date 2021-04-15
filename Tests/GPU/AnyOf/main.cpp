#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>

#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;
void main_main ();

// ================================================

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    main_main();

    amrex::Finalize();
}

void main_main ()
{
    BL_PROFILE("main");

    int ncell, ncomp, max_grid_size, nitem;
    int nghost = 0;
    {
        ParmParse pp;
        pp.get("ncell", ncell);
        pp.get("ncomp", ncomp);
        pp.get("nitem", nitem);
        pp.get("max_grid_size", max_grid_size);
        pp.query("nghost", nghost);
    }

    MultiFab mf;
    Gpu::DeviceVector<Real> vec(nitem, 7.0);

// ***************************************************************
    // Build the Multifabs and Geometries.
    {
        Box domain(IntVect{0}, IntVect{ncell-1});
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);

        mf.define(ba, dm, ncomp, nghost);
        mf.setVal(3.0);
    }

    {
        BL_PROFILE("Vector AnyOf");

        Gpu::Device::synchronize();

        bool anyof_M = Reduce::AnyOf(nitem, vec.dataPtr(),
                        [=] AMREX_GPU_DEVICE (Real item) noexcept -> int
                        {
                            return ( item < 2.0 ) ;
                        });

        bool anyof_N = Reduce::AnyOf(nitem, vec.dataPtr(),
                        [=] AMREX_GPU_DEVICE (Real item) noexcept -> int
                        {
                            return ( item > 2.0 ) ;
                        });
        amrex::Print() << "Vector, (0, 1): "
                       << anyof_M << ", " << anyof_N << std::endl;
    }

    // Redo, confirming works for a single item.
    vec[int(nitem-1)] = 1.0;

    {
        BL_PROFILE("Vector AnyOf - Just 1");

        bool anyof_M = Reduce::AnyOf(nitem, vec.dataPtr(),
                        [=] AMREX_GPU_DEVICE (Real item) noexcept -> int
                        {
                            return ( item < 2.0 ) ;
                        });

        bool anyof_N = Reduce::AnyOf(nitem, vec.dataPtr(),
                        [=] AMREX_GPU_DEVICE (Real item) noexcept -> int
                        {
                            return ( item > 2.0 ) ;
                        });

        amrex::Print() << "Vector, both true: "
                       << anyof_M << ", " << anyof_N << std::endl;
    }

    amrex::Print() << "\nTesting Box Version.\n One box = 0,1,1,1, all others 0,1,0,1\n";
    {
        BL_PROFILE("Box AnyOf");

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            Array4<Real> arr = mf.array(mfi);

            bool anyof_A = Reduce::AnyOf(bx,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                           {
                               return ( int(arr(i,j,k) < 2.0) ) ;
                           });

            bool anyof_B = Reduce::AnyOf(bx,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                           {
                               return ( int(arr(i,j,k) > 2.0) ) ;
                           });

            // Redo, confirming works for a single value.
            if (bx.contains(IntVect{ncell/3,ncell/2,ncell-1}))
            {
                arr(ncell/3,ncell/2,ncell-1) = 1.0;
            }

            bool anyof_C = Reduce::AnyOf(bx,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                           {
                               return ( int(arr(i,j,k) < 2.0) ) ;
                           });

            bool anyof_D = Reduce::AnyOf(bx,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                           {
                               return ( int(arr(i,j,k) > 2.0) ) ;
                           });

            amrex::Print() << "Box #" << mfi.LocalIndex() << " = "
                           << anyof_A << ", " << anyof_B << ", "
                           << anyof_C << ", " << anyof_D << std::endl;
        }
    }

}
