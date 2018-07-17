
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
#include <AMReX_AmrData.H>

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
    std::cout << "Usage:" << '\n';
    std::cout << progName << '\n';
    std::cout << "    errorC=CrseErrorFileName" << '\n';
    std::cout << "    errorF=FineErrorFileName" << '\n';
    std::cout << "   [-help]" << '\n';
    std::cout << "   [-verbose]" << '\n';
    std::cout << '\n';
    exit(1);
}

bool
amrDatasHaveSameDerives(const AmrData& amrd1,
                        const AmrData& amrd2);

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
    std::string cFile, fFile;

    bool verbose = false;
    if (pp.contains("verbose"))
    {
        verbose = true;
        AmrData::SetVerbose(true);
    }
    pp.query("errorC", cFile);
    if (cFile.empty())
        amrex::Abort("You must specify `errorC'");

    pp.query("errorF", fFile);
    if (fFile.empty())
        amrex::Abort("You must specify `errorF'");

    Vector<Real> norm0c, norm1c, norm2c;
    Vector<Real> norm0f, norm1f, norm2f;

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    
    DataServices dataServicesC(cFile, fileType);
    DataServices dataServicesF(fFile, fileType);

    AmrData& amrDataC = dataServicesC.AmrDataRef();
    AmrData& amrDataF = dataServicesF.AmrDataRef();
    BL_ASSERT(amrDatasHaveSameDerives(amrDataC,amrDataF));

    ComputeAmrDataNorms(amrDataC, norm0c, norm1c, norm2c, verbose);
    ComputeAmrDataNorms(amrDataF, norm0f, norm1f, norm2f, verbose);

    // Write norms to screen
    if (ParallelDescriptor::IOProcessor())
    {
	const Vector<std::string>& names = amrDataC.PlotVarNames();
	int maxl = 0;
	for (int i=0; i<names.size(); ++i)
	    maxl = std::max(maxl,int(names[i].size()));

        std::string maxl_str = amrex::Concatenate("", maxl, 1);

	std::string formatStr =
	    std::string("\t%") + maxl_str + std::string("s |  %10e   %10e   %10e\n");
	std::string sformatStr =
	    std::string("\t%") + maxl_str + std::string("s |  %10s   %10s   %10s\n");
	
	std::cout << '\n' << "Rates for pltfiles = "
	     << cFile << ", "
	     << fFile << ": " << '\n' << '\n';
	printf(sformatStr.c_str(),"Derived","rate_L-inf","rate_L1","rate_L2");
	std::cout << '\t'
	     << "--------------+------------------------------------------" << '\n';

	Real log2 = log(2.0);
	for (int i=0; i<names.size(); ++i)
	{
	    Real rate0, rate1, rate2;
	    rate0 = (norm0f[i]==0 ? 0.0 : log(norm0c[i]/norm0f[i])/log2);
	    rate1 = (norm1f[i]==0 ? 0.0 : log(norm1c[i]/norm1f[i])/log2);
	    rate2 = (norm2f[i]==0 ? 0.0 : log(norm2c[i]/norm2f[i])/log2);
	    printf(formatStr.c_str(),names[i].c_str(),rate0,rate1,rate2);
	}
	std::cout << '\n';
	
    }
    
    amrex::Finalize();
}


bool
amrDatasHaveSameDerives(const AmrData& amrd1,
                        const AmrData& amrd2)
{
    const Vector<std::string>& derives1 = amrd1.PlotVarNames();
    const Vector<std::string>& derives2 = amrd2.PlotVarNames();
    int length = derives1.size();
    if (length != derives2.size())
        return false;
    for (int i=0; i<length; ++i)
        if (derives1[i] != derives2[i])
            return false;
    return true;
}
    
