
//
// $Id: Add.cpp,v 1.3 2000-10-02 21:08:52 lijewski Exp $
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
    cout << "     factor = factor" << '\n';
    cout << "    outfile = outputFileName" << '\n';
    cout << "   [-help]" << '\n';
    cout << '\n';
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

    ParmParse pp(argc-1,argv+1);

    FArrayBox::setFormat(FABio::FAB_IEEE_32);
    //
    // Scan the arguments.
    //
    aString iFileDir, iFile, eFile, oFile, oFileDir;
    REAL factor;

    pp.query("infile", iFile);
    if (iFile.isNull())
        BoxLib::Abort("You must specify `infile'");

    pp.query("factor", factor);

    pp.query("outfile", oFile);
    if (oFile.isNull())
        BoxLib::Abort("You must specify `outfile'");

    ifstream is(iFile.c_str(),ios::in);
    ofstream os(oFile.c_str(),ios::out);

    FARRAYBOX dataI, dataE;
    dataI.readFrom(is);
  
    dataI.plus(factor);   

    dataI.writeOn(os);
}
