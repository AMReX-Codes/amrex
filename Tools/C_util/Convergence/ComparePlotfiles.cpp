
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

    int ncomp = mf_c.nComp();

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
    Print() << "ba_f " << ba_f << std::endl;

    // minimalBox() computes a single box to enclose all the boxes
    // enclosedCells() converts it to a cell-centered Box
    Box bx_c = ba_c.minimalBox().enclosedCells();
    Box bx_f = ba_f.minimalBox().enclosedCells();

    // number of cells in the coarse domain
    long numPts = bx_c.numPts();

    Print() << "numPts in coarse domain = " << numPts << std::endl;
    
    // assume ref_ratio is the same in each direction
    int rr = bx_f.length(0)/bx_c.length(0);

    Print() << "ref_ratio = " << rr << std::endl;
    
    // check to make sure refinement ratio is an integer
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bx_f.length(i)%bx_c.length(i) != 0) {
            Abort("not an integer refinement ratio");
        }
    }

    // check to make sure refinement ratio is the same in each direction
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if ( bx_f.length(i)/bx_c.length(i) != rr ) {
            Abort("ref_ratio not the same in each direction");
        }
    }

    // make a new BoxArray that is a refined version of the coarse BoxArray with the same
    // problem domain as the fine BoxArray
    BoxArray ba_c2 = ba_c;
    ba_c2.refine(rr);

    // grab the distribtion map from the coarse MultiFab
    DistributionMapping dm = mf_c.DistributionMap();

    Print() << "ba_c2 " << ba_c2 << std::endl;
    
    // create a fine MultiFab with same distribution mapping as coarse MultiFab
    MultiFab mf_f2(ba_c2,dm,ncomp,0);

    // copy fine data into new fine MultiFab
    mf_f2.ParallelCopy(mf_f,0,0,ncomp,0,0);

    // storage for averaged-down fine MultiFab
    MultiFab mf_c2(ba_c,dm,ncomp,0);

    // now we average down mf_f2 into mf_c2

    int how_many_nodal = 0;
    for (int i=0; i<AMREX_SPACEDIM; ++i ) {
        if (c_nodality[i] == 1) {
            ++how_many_nodal;            
        }
    }

    if (how_many_nodal == 0) {
        // cell-centered
        // average down ref_ratio^dim fine cells into coarse
        int npts = pow(rr,AMREX_SPACEDIM);
        
        for ( MFIter mfi(mf_c,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        
            const Box& bx = mfi.tilebox();

            const Array4<Real const>& fine   = mf_f2.array(mfi);
            const Array4<Real      >& coarse = mf_c2.array(mfi);

            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                coarse(i,j,k,n) = 0.;

#if (AMREX_SPACEDIM==3)
                for (int kk=0; kk<rr; ++kk) {
#else
                    int kk=0;    
#endif
                    for (int jj=0; jj<rr; ++jj) {
                        for (int ii=0; ii<rr; ++ii) {
                            coarse(i,j,k,n) += fine(rr*i+ii,rr*j+jj,rr*k+kk,n);
                        }
                    }
#if (AMREX_SPACEDIM==3)
                }
#endif
                coarse(i,j,k,n) /= npts;            
            });
            
        } // end MFIter
        
    } else if (how_many_nodal == 1) {
        // face-centered
        // average down ref_ratio^{dim-1} fine faces into coarse
        
    } else if (how_many_nodal == 2 && AMREX_SPACEDIM == 3) {
        // edge (3D)
        // average down ref_ratio^{dim-2} fine edges into coarse

    } else {
        // nodal
        // direct inject fine nodes into coarse nodes
        
    }

    // subtract coarse from coarsened fine
    MultiFab::Subtract(mf_c2,mf_c,0,0,ncomp,0);
    
    // force periodicity so faces/edges/nodes get weighted accordingly for L1 and L2 norms
    Abort("FIXME GENERALIZE");
    IntVect iv(AMREX_D_DECL(16,16,16));
    Periodicity period(iv);

    // compute norms of mf_c2
    for (int i=0; i<ncomp; ++i) {
        Real norm = mf_c2.norm0(i);
        Print() << "L0 comp " << i << " " << norm << std::endl;
    }

    for (int i=0; i<ncomp; ++i) {
        Real norm = mf_c2.norm1(i,period);
        Print() << "L1 comp " << i << " " << norm/numPts << std::endl;
    }

    for (int i=0; i<ncomp; ++i) {
        Real norm = mf_c2.norm2(i,period);
        Print() << "L2 comp " << i << " " << norm/sqrt(numPts) << std::endl;
    }
    
    



    

    // loop over components
    
      // average 'fine2' down to 'coarse2' (check nodality to call the correct averaging routine)

      // subtract 'coarse' and 'coarse2'
    
      // apply the OverlapMask and compute the norm (L0, L1, or L2)

    // write out results

    
    
    
}
