
//
// $Id: DiffFab.cpp,v 1.5 2000-10-02 21:08:53 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
using std::ios;
using std::set_new_handler;
#else
#include <new.h>
#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#endif

#include <unistd.h>

#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "Utility.H"
#include "VisMF.H"

#ifndef NDEBUG
#include "TV_TempWrite.H"
#endif

#include "AVGDOWN_F.H"

//
// This MUST be defined if don't have pubsetbuf() in I/O Streams Library.
//
#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif

#define GARBAGE 666.e+40

static
void
PrintUsage (const char* progName)
{
    cout << '\n';
    cout << "Usage:" << '\n';
    cout << progName << '\n';
    cout << "    infile  = inputFileName" << '\n';
    cout << "    exact   = exactFileName" << '\n';
    cout << "    outfile = outputFileName" << '\n';
    cout << "   [-help]" << '\n';
    cout << '\n';
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
    //
    // Make sure to catch new failures.
    //
    set_new_handler(Utility::OutOfMemory);

    ParmParse pp(argc-1,argv+1);

    FArrayBox::setFormat(FABio::FAB_IEEE_32);
    //
    // Scan the arguments.
    //
    aString iFileDir, iFile, eFile, oFile, oFileDir;


    pp.query("infile", iFile);
    if (iFile.isNull())
        BoxLib::Abort("You must specify `infile'");

    pp.query("exact", eFile);
    if (eFile.isNull())
        BoxLib::Abort("You must specify `exact' file");

    pp.query("outfile", oFile);
    if (oFile.isNull())
        BoxLib::Abort("You must specify `outfile'");

    ifstream is1(iFile.c_str(),ios::in);
    ifstream is2(eFile.c_str(),ios::in);
    ofstream os (oFile.c_str(),ios::out);

    FARRAYBOX dataI, dataE;
    dataI.readFrom(is1);
    dataE.readFrom(is2);

    BL_ASSERT(dataI.nComp() == dataE.nComp());
	
    //
    // Compute the error
    //
    int nComp = dataI.nComp();
    const BOX& domainI = dataI.box();
    const BOX& domainE = dataE.box();
    IntVect refine_ratio = getRefRatio(domainI, domainE);

    if (refine_ratio == IntVect())
      BoxLib::Error("Cannot find refinement ratio from data to exact");
    
    FARRAYBOX error(domainI,nComp);
    error.setVal(GARBAGE);
 
    FARRAYBOX exactAvg(domainI,nComp);
      
    FORT_CV_AVGDOWN(exactAvg.dataPtr(),
		    ARLIM(exactAvg.loVect()), 
		    ARLIM(exactAvg.hiVect()), &nComp,
		    dataE.dataPtr(),
		    ARLIM(dataE.loVect()), ARLIM(dataE.hiVect()),
		    domainI.loVect(), domainI.hiVect(),
		    refine_ratio.getVect());

    error.copy(exactAvg);
    error.minus(dataI);

    error.writeOn(os);
}
