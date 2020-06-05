
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

    // read in only level 0 for now
    iFile1 += "/Level_0/Cell";
    iFile2 += "/Level_0/Cell";
    
    MultiFab mf_c, mf_f;
    VisMF::Read(mf_c, iFile1);
    VisMF::Read(mf_f, iFile2);

    Print() << "HACK " << mf_c.nComp();

        
    // single-level for now
    // AMR comes later, where we iterate over each level in isolation

    // read in plotfiles, 'coarse' and 'fine' to MultiFabs
    // note: fine could be the same resolution as coarse

    // get the nodality of the data

    // check to make sure they have the same number of variables and variable names
      // const Vector<std::string>& varNames () const noexcept { return m_impl->varNames(); }
    
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
