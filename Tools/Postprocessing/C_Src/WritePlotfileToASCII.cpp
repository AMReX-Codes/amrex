
#include <fstream>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

using namespace std;
using namespace amrex;

static
void
PrintUsage (const char* progName)
{
    Print() << std::endl
            << "This utility reads in a single-level plotfile, copies it to " << std::endl
            << "a MultiFab with a single box, then writes out all the data in" << std::endl
            << "'i j k comp <value>' format to the screen, for each component one by one." << std::endl
            << "The optional flag [comp_in_line = 1] will plot data " << std::endl
            << "with the format (i,j,k) followed by all components in a row.  " << std::endl
            << "The user can modify this cpp file to write out on certain components," << std::endl
            << "coordinates, row/column formatting, etc." << std::endl << std::endl;

    Print() << "Usage:" << '\n';
    Print() << progName << " infile=inputFileName [comp_in_line = 1]" << '\n' << '\n';

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
    std::string iFile;

    // read in parameters from inputs file
    ParmParse pp;

    // coarse MultiFab
    pp.query("infile", iFile);
    if (iFile.empty())
        amrex::Abort("You must specify `infile'");

    int comp_in_line = 0;
    pp.query("comp_in_line", comp_in_line);

    // single-level for now
    // AMR comes later, where we iterate over each level in isolation

    // for the Header
    std::string iFile2 = iFile;
    iFile2 += "/Header";

    // open header
    ifstream x;
    x.open(iFile2.c_str(), ios::in);

    // read in first line of header
    string str;
    x >> str;

    // read in number of components from header
    int ncomp;
    x >> ncomp;

    // read in variable names from header
    for (int n=0; n<ncomp; ++n) {
        x >> str;
    }

    // read in dimensionality from header
    int dim;
    x >> dim;

    if (dim != AMREX_SPACEDIM) {
        Print() << "\nError: you are using a " << AMREX_SPACEDIM << "D build to open a "
                << dim << "D plotfile\n\n";
        Abort();
    }

    // now read in the plotfile data
    // check to see whether the user pointed to the plotfile base directory
    // or the data itself
    if (amrex::FileExists(iFile+"/Level_0/Cell_H")) {
       iFile += "/Level_0/Cell";
    }
    if (amrex::FileExists(iFile+"/Level_00/Cell_H")) {
       iFile += "/Level_00/Cell";
    }

    // storage for the input coarse and fine MultiFabs
    MultiFab mf;

    // read in plotfiles, 'coarse' and 'fine' to MultiFabs
    // note: fine could be the same resolution as coarse
    VisMF::Read(mf, iFile);

    ncomp = mf.nComp();
    Print() << "ncomp = " << ncomp << std::endl;

    // check nodality
    IntVect c_nodality = mf.ixType().toIntVect();
    Print() << "nodality " << c_nodality << std::endl;

    // get boxArray
    BoxArray ba = mf.boxArray();

    // minimalBox() computes a single box to enclose all the boxes
    // enclosedCells() converts it to a cell-centered Box
    Box bx_onegrid = ba.minimalBox().enclosedCells();

    // number of cells in the coarse domain
    Print() << "npts in coarse domain = " << bx_onegrid.numPts() << std::endl;
    long npts_coarsedomain = bx_onegrid.numPts();

    // BoxArray, DistributionMapping, and MultiFab with one grid
    BoxArray ba_onegrid(bx_onegrid);
    DistributionMapping dmap_onegrid(ba_onegrid);
    MultiFab mf_onegrid(ba_onegrid,dmap_onegrid,ncomp,0);

    // copy data into MultiFab with one grid
    mf_onegrid.ParallelCopy(mf,0,0,ncomp,0,0);

    for ( MFIter mfi(mf_onegrid,false); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.validbox();
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const Array4<Real>& mfdata = mf_onegrid.array(mfi);

        if (comp_in_line == 1){
          std::cout << mf_onegrid[mfi];
        }else{
          for (auto n=0; n<ncomp; ++n) {
            for (auto k = lo.z; k <= hi.z; ++k) {
              for (auto j = lo.y; j <= hi.y; ++j) {
                for (auto i = lo.x; i <= hi.x; ++i) {
                  std::cout << i << " " << j << " " << k << " " << n << " " << mfdata(i,j,k,n) << "\n";
                }
              }
            }
          }
        }
    } // end MFIter

}
