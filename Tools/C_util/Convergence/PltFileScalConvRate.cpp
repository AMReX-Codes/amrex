
#ifdef BL_USE_NEW_HFILES
#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
using std::ios;
using std::set_new_handler;
#else
#include <new.h>
#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#endif

#include <unistd.h>

#include "ComputeAmrDataNorms.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"

#ifndef NDEBUG
#include "TV_TempWrite.H"
#endif

//
// This MUST be defined if don't have pubsetbuf() in I/O Streams Library.
//
#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif

static
void
PrintUsage (const char* progName)
{
    cout << '\n';
    cout << "Usage:" << '\n';
    cout << progName << '\n';
    cout << "    errorC=CrseErrorFileName" << '\n';
    cout << "    errorF=FineErrorFileName" << '\n';
    cout << "   [-help]" << '\n';
    cout << "   [-verbose]" << '\n';
    cout << '\n';
    exit(1);
}

bool
amrDatasHaveSameDerives(const AmrData& amrd1,
                        const AmrData& amrd2);

int
main (int   argc,
      char* argv[])
{
    if (argc == 1)
        PrintUsage(argv[0]);
    //
    // Make sure to catch new failures.
    //
    set_new_handler(Utility::OutOfMemory);

    ParallelDescriptor::StartParallel(&argc, &argv);

    ParmParse pp(argc-1,argv+1);

    if (pp.contains("help"))
        PrintUsage(argv[0]);

    FArrayBox::setFormat(FABio::FAB_IEEE_32);
    //
    // Scan the arguments.
    //
    aString cFile, fFile;

    bool verbose = false;
    if (pp.contains("verbose"))
    {
        verbose = true;
        AmrData::SetVerbose(true);
    }
    pp.query("errorC", cFile);
    if (cFile.isNull())
        BoxLib::Abort("You must specify `errorC'");

    pp.query("errorF", fFile);
    if (fFile.isNull())
        BoxLib::Abort("You must specify `errorF'");

    Array<Real> norm0c, norm1c, norm2c;
    Array<Real> norm0f, norm1f, norm2f;

    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);
    
    DataServices dataServicesC(cFile, fileType);
    DataServices dataServicesF(fFile, fileType);

    if (!dataServicesC.AmrDataOk() || !dataServicesF.AmrDataOk())
        //
        // This calls ParallelDescriptor::EndParallel() and exit()
        //
        DataServices::Dispatch(DataServices::ExitRequest, NULL);

    AmrData& amrDataC = dataServicesC.AmrDataRef();
    AmrData& amrDataF = dataServicesF.AmrDataRef();
    BL_ASSERT(amrDatasHaveSameDerives(amrDataC,amrDataF));

    ComputeAmrDataNorms(amrDataC, norm0c, norm1c, norm2c, verbose);
    ComputeAmrDataNorms(amrDataF, norm0f, norm1f, norm2f, verbose);

    // Write norms to screen
    if (ParallelDescriptor::IOProcessor())
    {
	const Array<aString>& names = amrDataC.PlotVarNames();
	int maxl = 0;
	for (int i=0; i<names.length(); ++i)
	    maxl = Max(maxl,names[i].length());
	char sbuf[128];
	sprintf(sbuf,"%d",maxl);
	aString formatStr =
	    aString("\t%") + sbuf + aString("s |  %10e   %10e   %10e\n");
	aString sformatStr =
	    aString("\t%") + sbuf + aString("s |  %10s   %10s   %10s\n");
	
	cout << '\n' << "Rates for pltfiles = "
	     << cFile << ", "
	     << fFile << ": " << '\n' << '\n';
	printf(sformatStr.c_str(),"Derived","rate_L-inf","rate_L1","rate_L2");
	cout << '\t'
	     << "--------------+------------------------------------------" << '\n';

	Real log2 = log(2.0);
	for (int i=0; i<names.length(); ++i)
	{
	    Real rate0, rate1, rate2;
	    rate0 = (norm0f[i]==0 ? 0.0 : log(norm0c[i]/norm0f[i])/log2);
	    rate1 = (norm1f[i]==0 ? 0.0 : log(norm1c[i]/norm1f[i])/log2);
	    rate2 = (norm2f[i]==0 ? 0.0 : log(norm2c[i]/norm2f[i])/log2);
	    printf(formatStr.c_str(),names[i].c_str(),rate0,rate1,rate2);
	}
	cout << '\n';
	
    }
    
    //
    // This calls ParallelDescriptor::EndParallel() and exit()
    //
    DataServices::Dispatch(DataServices::ExitRequest, NULL);
}


bool
amrDatasHaveSameDerives(const AmrData& amrd1,
                        const AmrData& amrd2)
{
    const Array<aString>& derives1 = amrd1.PlotVarNames();
    const Array<aString>& derives2 = amrd2.PlotVarNames();
    int length = derives1.length();
    if (length != derives2.length())
        return false;
    for (int i=0; i<length; ++i)
        if (derives1[i] != derives2[i])
            return false;
    return true;
}
    
