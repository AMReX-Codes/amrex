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

Real MFdiff(const MultiFab& lhs, const MultiFab& rhs,
            int strt_comp, int num_comp, int nghost, const std::string name = "")
{
    MultiFab temp(rhs.boxArray(), rhs.DistributionMap(), rhs.nComp(), nghost);
    temp.ParallelCopy(lhs);
    temp.minus(rhs, strt_comp, num_comp, nghost);

    if (name != "")
        { amrex::VisMF::Write(temp, std::string("pltfiles/" + name)); }

    Real max_diff = 0;
    for (int i=0; i<num_comp; ++i)
    {
        Real max_i = std::abs(temp.max(i));
        max_diff = (max_diff > max_i) ? max_diff : max_i;
    }

    return max_diff; 
}

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

    int ncell, ncomp, max_grid_size;
    int nghost = 0;
    {
        ParmParse pp;
        pp.get("ncell", ncell);
        pp.get("ncomp", ncomp);
        pp.get("max_grid_size", max_grid_size);
        pp.query("nghost", nghost);
    }

    MultiFab mf_src, mf_dst;
    MultiFab mf_srcnew, mf_dstnew;

// ***************************************************************
    // Build the Multifabs and Geometries.
    {
        Box domain(IntVect{0}, IntVect{ncell-1});
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);

        DistributionMapping dm_src(ba);

        Vector<int> dst_map = dm_src.ProcessorMap();
        for (int& b : dst_map)
        {
           if(b != ParallelDescriptor::NProcs()-1) 
               { b++; } 
           else 
               { b=0; }
        }

        DistributionMapping dm_dst(dst_map);

        mf_src.define(ba, dm_src, ncomp, nghost);
        mf_src.setVal(3.0);

        mf_srcnew.define(ba, dm_src, ncomp, nghost);
        mf_srcnew.setVal(3.0);

        mf_dst.define(ba, dm_dst, ncomp, nghost);
        mf_dst.setVal(5.0);

        mf_dstnew.define(ba, dm_dst, ncomp, nghost);
        mf_dstnew.setVal(5.0);
/*
        amrex::Print() << " Copying from:\n " << dm_src << std::endl
                       << " to:\n "           << dm_dst << std::endl;
*/
    }

    {
        BL_PROFILE("**** ParallelCopy - 1");
        mf_dst.ParallelCopy(mf_src);
    }

    {
        BL_PROFILE("---- ParallelCopy_nowait -- 1");
        mf_dst.ParallelCopy_nowait(mf_src);
        mf_dst.ParallelCopy_finish();
    }


    amrex::Print() << "Error in old PC: " 
                   << MFdiff(mf_src, mf_dst, 0, ncomp, nghost) << std::endl;

    {
        BL_PROFILE("---- ParallelCopy_nowait -- 2");
        mf_dstnew.ParallelCopy_nowait(mf_srcnew);
        mf_dstnew.ParallelCopy_finish();
    }

    {
        BL_PROFILE("**** ParallelCopy - 2");
        mf_dstnew.ParallelCopy(mf_src);
    }


    {
        BL_PROFILE("---- ParallelCopy_nowait -- 3");
        mf_dst.ParallelCopy_nowait(mf_srcnew);
        mf_dst.ParallelCopy_finish();
    }

    {
        BL_PROFILE("**** ParallelCopy - 3");
        mf_src.ParallelCopy(mf_dstnew);
    }



    amrex::Print() << "Error in new PC: " 
                   << MFdiff(mf_srcnew, mf_dstnew, 0, ncomp, nghost) << std::endl;

}
