#include <iostream>

#include <AMReX.H>
#include <AMReX_RealBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_BLFort.H>

#include "MyAmr.H"

extern "C" void f(BL_FORT_FAB_ARG_3D(fab));

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
        RealBox prob_domain(AMREX_D_DECL(0.,0.,0.), AMREX_D_DECL(1.,1.,1.));
        int max_level = 2;
        Vector<int> n_cell{AMREX_D_DECL(64,64,64)};
        MyAmr amr(&prob_domain, max_level, n_cell);

        amr.MakeNewGrids();

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

        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            FArrayBox& fab = mf[mfi];
            f(BL_TO_FORTRAN_3D(fab));
        }
    }

    amrex::Finalize();
}

