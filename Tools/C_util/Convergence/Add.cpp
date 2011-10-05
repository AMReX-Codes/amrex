
#include <new>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string>
using std::ios;

#include <unistd.h>

#include <REAL.H>
#include <Box.H>
#include <FArrayBox.H>
#include <ParmParse.H>
#include <Utility.H>
#include <VisMF.H>

#ifndef NDEBUG
#include <TV_TempWrite.H>
#endif

#include <AVGDOWN_F.H>

#define GARBAGE 666.e+40

static
void
PrintUsage (const char* progName)
{
    std::cout << '\n';
    std::cout << "Usage:" << '\n';
    std::cout << progName << '\n';
    std::cout << "    infile  = inputFileName" << '\n';
    std::cout << "     factor = factor" << '\n';
    std::cout << "    outfile = outputFileName" << '\n';
    std::cout << "   [-help]" << '\n';
    std::cout << '\n';
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

    FArrayBox::setFormat(FABio::FAB_IEEE_32);
    //
    // Scan the arguments.
    //
    std::string iFileDir, iFile, eFile, oFile, oFileDir;
    Real factor;

    pp.query("infile", iFile);
    if (iFile.empty())
        BoxLib::Abort("You must specify `infile'");

    pp.query("factor", factor);

    pp.query("outfile", oFile);
    if (oFile.empty())
        BoxLib::Abort("You must specify `outfile'");

    std::ifstream is(iFile.c_str(),ios::in);
    std::ofstream os(oFile.c_str(),ios::out);

    FArrayBox dataI, dataE;
    dataI.readFrom(is);
  
    dataI.plus(factor);   

    dataI.writeOn(os);

    BoxLib::Finalize();
}
