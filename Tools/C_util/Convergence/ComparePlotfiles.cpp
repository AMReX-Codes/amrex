
#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
using std::ios;

#include <unistd.h>

#include <WritePlotFile.H>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>
#include <AMReX_AVGDOWN_F.H>

using namespace amrex;

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    // plotfile names for the coarse, fine, and subtracted output
    std::string iFile1, iFile2, difFile;

    // read in parameters from inputs file
    ParmParse pp;

    // coarse MultiFab
    pp.query("infile1", iFile1);
    if (iFile1.empty())
        amrex::Abort("You must specify `infile1'");

    // fine MultiFab (might have same resolution as coarse)
    pp.query("reffile", iFile2);
    if (iFile2.empty())
        amrex::Abort("You must specify `reffile'");

    // subtracted output (optional)
    pp.query("diffile", difFile);
        
    // single-level for now
    // AMR comes later, where we iterate over each level in isolation
    iFile1 += "/Level_0/Cell";
    iFile2 += "/Level_0/Cell";

    // storage for the input coarse and fine MultiFabs
    MultiFab mf_c, mf_f;
    
    // read in plotfiles, 'coarse' and 'fine' to MultiFabs
    // note: fine could be the same resolution as coarse
    VisMF::Read(mf_c, iFile1);
    VisMF::Read(mf_f, iFile2);

    // check number of components
    Print() << "ncomp_c = " << mf_c.nComp() << std::endl;
    Print() << "ncomp_f = " << mf_f.nComp() << std::endl;
    if (mf_c.nComp() != mf_f.nComp()) {
        Abort("plotfiles do not have the same number of variables");
    }

    // check nodality
    IntVect c_nodality = mf_c.ixType().toIntVect();
    IntVect f_nodality = mf_f.ixType().toIntVect();
    Print() << "c_nodality " << c_nodality << std::endl;
    Print() << "f_nodality " << f_nodality << std::endl;
    if (c_nodality != f_nodality) {
        Abort("plotfiles do not have the same nodality");
    }

    // get coarse and fine boxArrays
    BoxArray ba_c = mf_c.boxArray();
    BoxArray ba_f = mf_f.boxArray();

    Print() << "ba_c " << ba_c << std::endl;

    // minimalBox() computes a single box to enclose all the boxes
    // enclosedCells() converts it to a cell-centered Box
    Box bx_c = ba_c.minimalBox().enclosedCells();
    Box bx_f = ba_f.minimalBox().enclosedCells();

    // assume ref_ratio is the same in each direction
    int ref_ratio = bx_f.length(0)/bx_c.length(0);

    // check to make sure refinement ratio is an integer
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bx_f.length(i)%bx_c.length(i) != 0) {
            Abort("not an integer refinement ratio");
        }
    }

    // check to make sure refinement ration is the same in each direction
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if ( bx_f.length(i)/bx_c.length(i) != ref_ratio ) {
            Abort("ref_ratio not the same in each direction");
        }
    }

    // make a new BoxArray that is a refined version of the coarse BoxArray with the same
    // problem domain as the fine BoxArray
    BoxArray ba_c2 = ba_c;
    ba_c2.refine(ref_ratio);

    // grab the distribtion map from the coarse MultiFab
    DistributionMapping dm = mf_c.DistributionMap();

    // create a fine MultiFab with same distribution mapping as coarse MultiFab
    MultiFab mf_f2(ba_c2,dm,mf_f.nComp(),0);

    // copy fine data into new fine MultiFab
    mf_f2.ParallelCopy(mf_f,0,0,mf_f.nComp(),0,0);

    // storage for averaged-down fine MultiFab
    MultiFab mf_c2(ba_c,dm,mf_c.nComp(),0);

    // now we average down mf_f2 into mf_c2

    int how_many_nodal = 0;
    for (int i=0; i<AMREX_SPACEDIM; ++i ) {
        if (c_nodality[i] == 1) {
            ++how_many_nodal;
        }
    }

    if (how_many_nodal == 0) {
        // cell-centered
        
    } else if (how_many_nodal == 1) {
        // face-centered

    } else if (how_many_nodal == 2) {
        // nodal (2D)
        // edge  (3D)

    } else if (how_many_nodal == 3) {
        // nodal (3D)
        
    }
    
    

    
        // coarsenable
        
    // figure out refinement ratio, make sure it is 1, or an even number,
    // and is the same in all directions

    // get the coarse boxarray

    // refine this boxarray to the resolution of 'fine'

    // create a new multifab, 'fine2' with the same nodality and distribution mapping of 'coarse',
    // but with fine resolution
    
    // parallel copy 'fine' into 'fine2'

    // create a new multifab 'coarse2' (same nodality and distribution mapping of 'coarse')

    // loop over components
    
      // average 'fine2' down to 'coarse2' (check nodality to call the correct averaging routine)

      // subtract 'coarse' and 'coarse2'
    
      // apply the OverlapMask and compute the norm (L0, L1, or L2)

    // write out results

    
    
    
}
