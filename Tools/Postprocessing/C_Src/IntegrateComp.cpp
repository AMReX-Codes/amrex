/*
  A very simple example of reading a plotfile and calling a function to an analysis.  Here, we want
  to do a volume integral of a specified component.

  The output is an array of numbers. */
#include <string>
#include <iostream>

#include "AMReX_ParmParse.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile=<plotfilename> varName=v1 v2 ... \n";
  exit(1);
}

int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);

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

  Vector<int> pComp(nCompIn,-1);
  for (int i=0; i<nCompIn; ++i) {
    for (int j=0; j<plotVarNames.size() && pComp[i]<0; ++j) {
      if (plotVarNames[j] == varNames[i]) {pComp[i] = j;}
    }
    if (ParallelDescriptor::IOProcessor() && pComp[i]<0) {
      Abort("Cannot find variable="+varNames[i]+" in pltfile");
    }
  }

  const int nGrow = 0;
  const int finestLevel = amrData.FinestLevel();
  const int nLev = finestLevel + 1;

  Vector<Real> integrals(nCompIn,0);
  for (int lev=0; lev<nLev; ++lev) {
    const BoxArray ba = amrData.boxArray(lev);

    const BoxArray baf = lev < finestLevel
                               ? BoxArray(amrData.boxArray(lev+1)).coarsen(amrData.RefRatio()[lev])
                               : BoxArray();
    Real vol=1.0;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
      vol *= amrData.ProbSize()[d] / amrData.ProbDomain()[lev].length(d);
    }

    for (int n=0; n<nCompIn; ++n) {
      const MultiFab &pfData = amrData.GetGrids(lev,pComp[n]);
      const DistributionMapping dmap(pfData.DistributionMap());
      MultiFab mf(ba,dmap,1,nGrow);

      MultiFab::Copy(mf,pfData,0,0,1,nGrow);

      Real sm = 0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:sm)
#endif
      for (MFIter mfi(mf,true); mfi.isValid(); ++mfi) {
        FArrayBox& myFab = mf[mfi];
        const Box& box = mfi.tilebox();
        if (lev < finestLevel) {
          std::vector< std::pair<int,Box> > isects = baf.intersections(box);                
          for (int ii = 0; ii < isects.size(); ii++) {
            myFab.setVal(0,isects[ii].second,0,1);
          }
        }
        sm += myFab.sum(box,0,1);
      }

      integrals[n] += sm * vol;

    }
  }

  ParallelDescriptor::ReduceRealSum(&(integrals[0]),nCompIn);

  for (int n=0; n<nCompIn; ++n) {
    Print() << integrals[n] << " ";
  }
  Print() << std::endl;

  amrex::Finalize();
  return 0;
}

