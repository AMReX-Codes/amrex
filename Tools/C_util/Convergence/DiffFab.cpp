
#include <new>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string>
using std::ios;

#include <unistd.h>

#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#ifdef AMREX_DEBUG
#include <TV_TempWrite.H>
#endif

#include <AMReX_AVGDOWN_F.H>

#define GARBAGE 666.e+40
using namespace amrex;
static
void
PrintUsage (const char* progName)
{
    std::cout << '\n';
    std::cout << "Usage:" << '\n';
    std::cout << progName << '\n';
    std::cout << "    infile  = inputFileName" << '\n';
    std::cout << "    exact   = exactFileName" << '\n';
    std::cout << "    outfile = outputFileName" << '\n';
    std::cout << "   [-help]" << '\n';
    std::cout << '\n';
    exit(1);
}

IntVect
getRefRatio(const Box& crse,
	    const Box& fine)
{
    // Compute refinement ratio between crse and fine boxes, return invalid
    // IntVect if there is none suitable
    IntVect ref_ratio;
    for (int i=0; i<BL_SPACEDIM; ++i)
	ref_ratio[i] = fine.length(i) / crse.length(i);

    // Check results
    Box test1 = ::Box(fine).coarsen(ref_ratio);
    Box test2 = ::Box(test1).refine(ref_ratio);
    if (test1 != crse  ||  test2 != fine)
	ref_ratio = IntVect();
    return ref_ratio;
}

int
main (int   argc,
      char* argv[])
{
    if (argc == 1)
        PrintUsage(argv[0]);

    amrex::Initialize(argc,argv);

    ParmParse pp;

    FArrayBox::setFormat(FABio::FAB_IEEE_32);
    //
    // Scan the arguments.
    //
    std::string iFileDir, iFile, eFile, oFile, oFileDir;


    pp.query("infile", iFile);
    if (iFile.empty())
        amrex::Abort("You must specify `infile'");

    pp.query("exact", eFile);
    if (eFile.empty())
        amrex::Abort("You must specify `exact' file");

    pp.query("outfile", oFile);
    bool do_out = !oFile.empty();

    std::ifstream is1(iFile.c_str(),ios::in);
    std::ifstream is2(eFile.c_str(),ios::in);
    std::ofstream os;
    if (do_out) {
      os.open(oFile.c_str(),ios::out);
    }

    FArrayBox dataI, dataE;
    dataI.readFrom(is1);
    dataE.readFrom(is2);

    bool has_nan_I = dataE.contains_nan();
    bool has_nan_E = dataI.contains_nan();
    bool fabs_agree = false;

    if (!has_nan_I && !has_nan_E) {
      BL_ASSERT(dataI.nComp() == dataE.nComp());
	
      //
      // Compute the error
      //
      int nComp = dataI.nComp();
      const Box& domainI = dataI.box();
      const Box& domainE = dataE.box();
      IntVect refine_ratio = getRefRatio(domainI, domainE);

      if (refine_ratio == IntVect())
        amrex::Error("Cannot find refinement ratio from data to exact");
    
      FArrayBox error(domainI,nComp);
      error.setVal(GARBAGE);
 
      FArrayBox exactAvg(domainI,nComp);
      
      FORT_CV_AVGDOWN(exactAvg.dataPtr(),
                      ARLIM(exactAvg.loVect()), 
                      ARLIM(exactAvg.hiVect()), &nComp,
                      dataE.dataPtr(),
                      ARLIM(dataE.loVect()), ARLIM(dataE.hiVect()),
                      domainI.loVect(), domainI.hiVect(),
                      refine_ratio.getVect());

      error.copy(exactAvg);
      error.minus(dataI);
      if (do_out) {
        error.writeOn(os);
      }
      Box bx = error.box();
      int points = bx.numPts();

      Real norm0 = error.norm(0);
      Real norm1 = error.norm(1)/points;
      Real norm2 = error.norm(2)/sqrt(points);
      std::cout << "L0 NORM                           = " << norm0 << std::endl;
      std::cout << "L1 NORM (normalized by 1/N)       = " << norm1 << std::endl;
      std::cout << "L2 NORM (normalized by 1/sqrt(N)) = " << norm2 << std::endl;

      fabs_agree = (norm0 <= 1.e-20) && (norm1 <= 1.e-20) && (norm2 <= 1.e-20);
    }

    if (!has_nan_I && !has_nan_E && fabs_agree) {
      std::cout << "FABS AGREE" << std::endl;
    }

    amrex::Finalize();
}
