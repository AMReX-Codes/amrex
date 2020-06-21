
#include <fstream>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

static
void
PrintUsage (const char* progName)
{
    Print() << std::endl
            << "This utility performs a diff operation between two"           << std::endl
            << "plotfiles that have the same geometrical domain and nodality" << std::endl
            << "(supports all nodality types; cell, face, edge, node)"        << std::endl
            << "but possibly a factor of refinement between the cells,"       << std::endl
            << "and outputs the L0, L1, and L2 norms"                         << std::endl
            << "L1 = sum(|diff_ijk|)/npts_coarsedomain"                       << std::endl
            << "L2 = sqrt[sum(diff_ijk^2)]/sqrt(npts_coarsedomain)"           << std::endl
            <<  "(only single-level supported)"                               << std::endl << std::endl;
    
    Print() << "Usage:" << '\n';
    Print() << progName << '\n';
    Print() << "    infile1 = inputFileName1" << '\n';
    Print() << "    reffile = refinedPlotFile" << '\n';
    Print() << "    diffile = differenceFileName" << '\n';
    Print() << "              (If not specified no file is written)" << '\n' << '\n';
    
    Print() << "You can either point to the plotfile base directory itself, e.g."      << std::endl
            << "  infile=plt00000"                                                     << std::endl
            << "Or the raw data itself, e.g."                                          << std::endl
            << "  infile=plt00000/Level_0/Cell"                                        << std::endl
            << "the latter is useful for some applications that dump out raw"          << std::endl
            << "nodal data within a plotfile directory."                               << std::endl
            << "The program will first try appending 'Level_0/Cell'"                   << std::endl
            << "onto the specified filenames."                                         << std::endl
            << "If that _H file doesn't exist, it tries using the full specified name" << std::endl << std::endl;
        
    exit(1);
}

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    if (argc == 1) {
        PrintUsage(argv[0]);
    }

    // plotfile names for the coarse, fine, and subtracted output
    std::string iFile1, iFile2, difFile="";

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

    // check to see whether the user pointed to the plotfile base directory
    // or the data itself
    if (amrex::FileExists(iFile1+"/Level_0/Cell_H")) {
       iFile1 += "/Level_0/Cell";
    }
    if (amrex::FileExists(iFile2+"/Level_0/Cell_H")) {
       iFile2 += "/Level_0/Cell";
    }

    // storage for the input coarse and fine MultiFabs
    MultiFab mf_c, mf_f;
    
    // read in plotfiles, 'coarse' and 'fine' to MultiFabs
    // note: fine could be the same resolution as coarse
    VisMF::Read(mf_c, iFile1);
    VisMF::Read(mf_f, iFile2);

    // check number of components
    if (mf_c.nComp() != mf_f.nComp()) {
        Abort("plotfiles do not have the same number of variables");
    }
    int ncomp = mf_c.nComp();
    Print() << "ncomp = " << ncomp << std::endl;

    // check nodality
    IntVect c_nodality = mf_c.ixType().toIntVect();
    IntVect f_nodality = mf_f.ixType().toIntVect();
    if (c_nodality != f_nodality) {
        Abort("plotfiles do not have the same nodality");
    }
    Print() << "nodality " << c_nodality << std::endl;

    // get coarse and fine boxArrays
    BoxArray ba_c = mf_c.boxArray();
    BoxArray ba_f = mf_f.boxArray();

    // minimalBox() computes a single box to enclose all the boxes
    // enclosedCells() converts it to a cell-centered Box
    Box bx_c = ba_c.minimalBox().enclosedCells();
    Box bx_f = ba_f.minimalBox().enclosedCells();

    // number of cells in the coarse domain
    Print() << "npts in coarse domain = " << bx_c.numPts() << std::endl;
    Print() << "npts in fine   domain = " << bx_f.numPts() << std::endl;
    long npts_coarsedomain = bx_c.numPts();
    
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

    int npts_avg = pow(rr,AMREX_SPACEDIM-how_many_nodal);

    int rr_i = (c_nodality[0] == 0) ? rr : 1;
    int rr_j = (c_nodality[1] == 0) ? rr : 1;
#if (AMREX_SPACEDIM == 3)
    int rr_k = (c_nodality[2] == 0) ? rr : 1;
#else
    int rr_k = 0;
#endif
        
    for ( MFIter mfi(mf_c,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        
        const Box& bx = mfi.tilebox();

        const Array4<Real const>& fine   = mf_f2.array(mfi);
        const Array4<Real      >& coarse = mf_c2.array(mfi);

        amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            coarse(i,j,k,n) = 0.;

#if (AMREX_SPACEDIM==3)
            for (int kk=0; kk<rr_k; ++kk) {
#else
                int kk=0;    
#endif
                for (int jj=0; jj<rr_j; ++jj) {
                    for (int ii=0; ii<rr_i; ++ii) {
                        coarse(i,j,k,n) += fine(rr*i+ii,rr*j+jj,rr*k+kk,n);
                    }
                }
#if (AMREX_SPACEDIM==3)
            }
#endif
            coarse(i,j,k,n) /= npts_avg;
        });
            
    } // end MFIter        

    // subtract coarse from coarsened fine
    MultiFab::Subtract(mf_c2,mf_c,0,0,ncomp,0);
    
    // force periodicity so faces/edges/nodes get weighted accordingly for L1 and L2 norms
    IntVect iv(AMREX_D_DECL(bx_c.length(0),
                            bx_c.length(1),
                            bx_c.length(2)));
    Periodicity period(iv);

    // compute norms of mf_c2
    for (int i=0; i<ncomp; ++i) {
        Real norm0 = mf_c2.norm0(i);
        Real norm1 = mf_c2.norm1(i,period);
        Real norm2 = mf_c2.norm2(i,period);
        Print() << "(comp,L0,L1,L2) " << i << " "
                << norm0 << " "
                << norm1/npts_coarsedomain << " "
                << norm2/sqrt(npts_coarsedomain) << " " << std::endl;
    }

    // write out the subtracted plotfile if diffile was specified at the command line
    if (difFile != "") {

        // define the problem domain as (0,1) for now
        RealBox real_box({AMREX_D_DECL(0.,0.,0.)},
                         {AMREX_D_DECL(1.,1.,1.)});

        Vector<int> is_periodic(AMREX_SPACEDIM,1);

        // build a geometry object so we can use WriteSingleLevelPlotfile
        Geometry geom(bx_c,&real_box,CoordSys::cartesian,is_periodic.data());

        // give generic variables names for now
        Vector<std::string> varNames(ncomp);
        for (int i=0; i<ncomp; ++i) {
            varNames[i] = std::to_string(i);
        }

        WriteSingleLevelPlotfile(difFile,mf_c2,varNames,geom,0.,0);
    }
    
}
