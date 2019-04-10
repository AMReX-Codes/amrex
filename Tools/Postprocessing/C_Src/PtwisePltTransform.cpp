/*
  A very simple example of reading a plotfile and calling a function to perform a pointwise transformation
  based on a set of components that are specified by name on the command line.  The transformation is done
  in the accompanying fortran routine.  No grow cells are used, so the transformation cannot involve a 
  stencil operation.

  The output is a new plotfile with a single component set to the output of the transform routine.  This 
  new plotfile has metadata (number of levels, boxarray, grid spacing, etc) that is identical to the original
  plotfile.
 */
#include <string>
#include <iostream>

#include "AMReX_ParmParse.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>
#include <WritePlotFile.H>

#include <AMReX_BLFort.H>

extern "C" {
  void transform (const int* lo, const int* hi,
                  const amrex_real* sIn, const int* sInlo, const int* sInhi, const int* ncIn,
                  amrex_real* sOut, const int* sOutlo, const int* sOuthi, const int* ncOut);
}

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile=<plotfilename> varNames=v1 v2 ... \n";
  exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  vector<std::string> tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);
  {
    if (argc < 2)
      print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc,argv);

    std::string infile; pp.get("infile",infile);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int nv = pp.countval("varNames");
    Vector<std::string> varNames(nv); pp.getarr("varNames",varNames,0,nv);

    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    int nCompIn = varNames.size();
    Vector<int> destFillComps(nCompIn);
    for (int i=0; i<nCompIn; ++i) {
      destFillComps[i] = i;
    }


    for (int i=0; i<nCompIn; ++i) {
      int ivar = -1;
      for (int j=0; j<plotVarNames.size(); ++j) {
        if (plotVarNames[j] == varNames[i]) {ivar = j;}
      }
      if (ParallelDescriptor::IOProcessor() && ivar<0) {
        Abort("Cannot find variable="+varNames[i]+" in pltfile");
      }
    }

    const int nCompOut = 1;
    const int nGrow = 0;
    const int nLev = amrData.FinestLevel() + 1;

    Vector<MultiFab*> stateOut(nLev);
    for (int lev=0; lev<nLev; ++lev) {
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dmap(ba);
      MultiFab stateIn(ba,dmap,nCompIn,nGrow);
      stateOut[lev] = new MultiFab(ba,dmap,nCompOut,0);

      // Load input data from pltfile 
      amrData.FillVar(stateIn,lev,varNames,destFillComps);

      // Compute transformation
      for (MFIter mfi(stateIn); mfi.isValid(); ++mfi) {
        const FArrayBox& sIn = stateIn[mfi];
        FArrayBox& sOut = (*stateOut[lev])[mfi];
        const Box& box = mfi.validbox();

        transform(BL_TO_FORTRAN_BOX(box),
                  BL_TO_FORTRAN_ANYD(sIn),&nCompIn,
                  BL_TO_FORTRAN_ANYD(sOut),&nCompOut);
      
      }
    }

    // Write result to new plotfile in local folder
    std::string outfile=getFileRoot(infile) + "_tr";
    Vector<std::string> outNames;
    outNames.push_back("transform");
    WritePlotFile(stateOut,amrData,outfile,false,outNames);
  }
  amrex::Finalize();
  return 0;
}

