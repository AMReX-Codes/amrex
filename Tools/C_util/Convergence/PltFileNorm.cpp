//BL_COPYRIGHT_NOTICE

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
    cout << "    infile=inputFileName" << '\n';
    cout << "   [outfile=outputFileName]" << '\n';
    cout << "   [-help]" << '\n';
    cout << "   [-verbose]" << '\n';
    cout << '\n';
    cout << " Note: outfile required if verbose used" << '\n';
    exit(1);
}

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
    aString iFile;

    bool verbose = false;
    if (pp.contains("verbose"))
    {
        verbose = true;
        AmrData::SetVerbose(true);
    }
    pp.query("infile", iFile);
    if (iFile.isNull())
        BoxLib::Abort("You must specify `infile'");

    Array<Real> norm0, norm1, norm2;

    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);
    
    DataServices dataServices(iFile, fileType);

    if (!dataServices.AmrDataOk())
        //
        // This calls ParallelDescriptor::EndParallel() and exit()
        //
        DataServices::Dispatch(DataServices::ExitRequest, NULL);

    AmrData& amrData = dataServices.AmrDataRef();

    ComputeAmrDataNorms(amrData, norm0, norm1, norm2, verbose);

    // Write norms to screen
    if (ParallelDescriptor::IOProcessor())
    {
	const Array<aString>& names = amrData.PlotVarNames();
	int maxl = 0;
	for (int i=0; i<names.length(); ++i)
	    maxl = Max(maxl,names[i].length());
	char sbuf[128];
	sprintf(sbuf,"%d",maxl);
	aString formatStr =
	    aString("\t%") + sbuf + aString("s |  %10e   %10e   %10e\n");
	aString sformatStr =
	    aString("\t%") + sbuf + aString("s |  %10s   %10s   %10s\n");
	
	cout << '\n' << "Norms for pltfile = " << iFile << ": " << '\n' << '\n';
	printf(sformatStr.c_str(),"Derived","L-inf","L1","L2");
	cout << '\t'
	     << "--------------+------------------------------------------" << '\n';
	
	for (int i=0; i<names.length(); ++i)
	{
	    printf(formatStr.c_str(),names[i].c_str(),norm0[i],norm1[i],norm2[i]);
	}
	cout << '\n';
	
    }
    
    //
    // This calls ParallelDescriptor::EndParallel() and exit()
    //
    DataServices::Dispatch(DataServices::ExitRequest, NULL);
}


