
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

    ParmParse pp;

    std::string iFile1, iFile2, difFile;

    pp.query("infile1", iFile1);
    if (iFile1.empty())
        amrex::Abort("You must specify `infile1'");

    pp.query("reffile", iFile2);
    if (iFile2.empty())
        amrex::Abort("You must specify `reffile'");

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
