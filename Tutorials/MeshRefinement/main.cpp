#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>

#include <MyAmr.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    // Add some ParmParse parameters.
    // This could also be done by putting them in inputs file
    {
        ParmParse ppg("geometry");
        ppg.add("coord_sys", 0);
        ppg.addarr("prob_lo", std::vector<Real>{0.0, 0.0, 0.0});
        ppg.addarr("prob_hi", std::vector<Real>{1.0, 1.0, 1.0});
        // By default, geometry.is_periodic is false

        ParmParse ppa("amr");
        ppa.add("max_level", 2);
        ppa.addarr("n_cell", std::vector<int>{64,64,64});
        // amr.max_grid_size = 32
        // amr.blocking_factor = 8
    }

    {
        MyAmr amr;

        Array<BoxArray> grids(amr.maxLevel()+1);
        Array<DistributionMapping> dmap(amr.maxLevel()+1);

        grids[0] = amr.MakeBaseGrids();
        dmap[0] = DistributionMapping{grids[0]};

        // must call these for MakeNewGrids to work
        amr.SetFinestLevel(0);
        amr.SetBoxArray(0, grids[0]);
        amr.SetDistributionMap(0, dmap[0]);

        Real time = 0.0;

        for (int lbase = 0; lbase < amr.maxLevel(); ++lbase)
        {
            int new_finest;

            // Add (at most) one level at a time.
	    amr.MakeNewGrids(lbase, time, new_finest, grids);

            if (new_finest <= amr.finestLevel()) break;

            amr.SetFinestLevel(new_finest);
            amr.SetBoxArray(new_finest, grids[new_finest]);
            amr.SetDistributionMap(new_finest, DistributionMapping{grids[new_finest]});
        }
        
        MultiFab mf(amr.boxArray(0), amr.DistributionMap(0), 1, 0);
        mf.setVal(0.0);

        IntVect ref_ratio = IntVect::TheUnitVector();
        for (int lev = 1; lev <= amr.finestLevel(); ++lev)
        {
            MultiFab fmf(amr.boxArray(lev), amr.DistributionMap(lev), 1, 0);
            fmf.setVal(static_cast<Real>(lev));
            ref_ratio *= amr.refRatio(lev-1);
            amrex::average_down(fmf, mf, 0, 1, ref_ratio);
        }

        VisMF::Write(mf, "mf");
    }

    amrex::Finalize();
}

