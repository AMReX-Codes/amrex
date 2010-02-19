
#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
using std::ios;

#include <unistd.h>

#include "ComputeAmrDataNorms.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"

#ifndef NDEBUG
#include "TV_TempWrite.H"
#endif

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
    BoxLib::Initialize(argc,argv);

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

    bool verbose = false;
    if (pp.contains("verbose"))
    {
        verbose = true;
        AmrData::SetVerbose(true);
    }
    pp.query("infile", iFile);
    if (iFile.empty())
        BoxLib::Abort("You must specify `infile'");

    Array<Real> norm0, norm1, norm2;

    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);
    
    DataServices dataServices(iFile, fileType);

    AmrData& amrData = dataServices.AmrDataRef();

    ComputeAmrDataNorms(amrData, norm0, norm1, norm2, verbose);

    // Write norms to screen
    if (ParallelDescriptor::IOProcessor())
    {
	const Array<std::string>& names = amrData.PlotVarNames();
	int maxl = 0;
	for (int i=0; i<names.size(); ++i)
	    maxl = std::max(maxl,int(names[i].size()));
	char sbuf[128];
	sprintf(sbuf,"%d",maxl);
	std::string formatStr =
	    std::string("\t%") + sbuf + std::string("s |  %10e   %10e   %10e\n");
	std::string sformatStr =
	    std::string("\t%") + sbuf + std::string("s |  %10s   %10s   %10s\n");
	
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
    
    BoxLib::Finalize();
}


