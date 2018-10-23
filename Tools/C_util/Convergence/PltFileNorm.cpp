
#include <iostream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <unistd.h>

#include <ComputeAmrDataNorms.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>
#include <AMReX_Utility.H>

#ifdef AMREX_DEBUG
#include <TV_TempWrite.H>
#endif

using std::ios;
using namespace amrex;

static
void
PrintUsage (const char* progName)
{
    std::cout << '\n';
    std::cout << "This routine reads a pltfile and calculates the Linfty," << std::endl
         << "L1 and L2 norms ov every component.                    " << std::endl
         << std::endl;
    std::cout << "Usage:" << '\n';
    std::cout << progName << '\n';
    std::cout << "    infile=inputFileName" << '\n';
    std::cout << "   [outfile=outputFileName]" << '\n';
    std::cout << "   [-help]" << '\n';
    std::cout << "   [-verbose]" << '\n';
    std::cout << '\n';
    std::cout << " Note: outfile required if verbose used" << '\n';
    exit(1);
}

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    if (argc == 1)
        PrintUsage(argv[0]);

    ParmParse pp;

    if (pp.contains("help"))
        PrintUsage(argv[0]);

    FArrayBox::setFormat(FABio::FAB_IEEE_32);
    //
    // Scan the arguments.
    //
    std::string iFile;
    int integrate = 1;
    pp.query("integrate", integrate);
    bool verbose = false;
    if (pp.contains("verbose"))
    {
        verbose = true;
        AmrData::SetVerbose(true);
    }
    pp.query("infile", iFile);
    if (iFile.empty())
        amrex::Abort("You must specify `infile'");

    Vector<Real> norm0, norm1, norm2;

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    
    DataServices dataServices(iFile, fileType);

    AmrData& amrData = dataServices.AmrDataRef();

    
    
    if (integrate)
    {
      ComputeAmrDataInt(amrData, norm1,  verbose);

      // Write norms to screen
      if (ParallelDescriptor::IOProcessor())
      {
	const Vector<std::string>& names = amrData.PlotVarNames();
	int maxl = 0;
	for (int i=0; i<names.size(); ++i)
	    maxl = std::max(maxl,int(names[i].size()));

        std::string maxl_str = amrex::Concatenate("", maxl, 1);

	std::string formatStr =
	    std::string("\t%") + maxl_str + std::string("s |  %10e   \n");
	std::string sformatStr =
	    std::string("\t%") + maxl_str + std::string("s |  %10s   \n");
	
	std::cout << '\n' << "Norms for pltfile = " << iFile << ": " << '\n' << '\n';
	printf(sformatStr.c_str(),"Derived", "Integral");
	std::cout << '\t'
	     << "--------------+------------------------------------------" << '\n';
	
	for (int i=0; i<names.size(); ++i)
	{
	    printf(formatStr.c_str(),names[i].c_str(),norm1[i]);
	}
	std::cout << '\n';
	
      }
    }
    else
    {  
      ComputeAmrDataNorms(amrData, norm0, norm1, norm2, verbose);

      // Write norms to screen
      if (ParallelDescriptor::IOProcessor())
      {
	const Vector<std::string>& names = amrData.PlotVarNames();
	int maxl = 0;
	for (int i=0; i<names.size(); ++i)
	    maxl = std::max(maxl,int(names[i].size()));

        std::string maxl_str = amrex::Concatenate("", maxl, 1);

	std::string formatStr =
	    std::string("\t%") + maxl_str + std::string("s |  %10e   %10e   %10e\n");
	std::string sformatStr =
	    std::string("\t%") + maxl_str + std::string("s |  %10s   %10s   %10s\n");
	
	std::cout << '\n' << "Norms for pltfile = " << iFile << ": " << '\n' << '\n';
	printf(sformatStr.c_str(),"Derived","L-inf","L1","L2");
	std::cout << '\t'
	     << "--------------+------------------------------------------" << '\n';
	
	for (int i=0; i<names.size(); ++i)
	{
	    printf(formatStr.c_str(),names[i].c_str(),norm0[i],norm1[i],norm2[i]);
	}
	std::cout << '\n';
	
      }
    }
    amrex::Finalize();
}


