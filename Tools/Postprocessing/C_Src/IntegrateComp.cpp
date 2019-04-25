/*
  A very simple example of reading a plotfile and doing a simple analysis.  Here, we want
  to do a volume integral of a component specified by name.

  The twist here is to demonstrate what this might look like if the amr data were coming down a 
  pipe (such as SENSEI, e.g.).  So, we read the plotfile as usual, but mine it for the required
  data, then put that data into a couple of structs that then are queried from that point forward.

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

struct AMReXMeshHierarchy
{
/*
  A minimal class to describe the AMR hierarchy for analysis routines
 */
public:
  AMReXMeshHierarchy() {}
  void define(const AmrData& ad) {
    finestLevel = ad.FinestLevel();
    int nlevs = finestLevel + 1;
    ba.resize(nlevs);
    probSize = ad.ProbSize();
    probDomain = ad.ProbDomain();
    refRatio = ad.RefRatio();
    for (int lev=0; lev<nlevs; ++lev) {
      ba[lev] = &ad.boxArray(lev);
    }
  }
  int FinestLevel() const {return finestLevel;}
  const BoxArray &boxArray(int level) const {return *ba[level];}
  const Vector<int> &RefRatio() const {return refRatio;}
  const Vector<Real> &ProbSize() const {return probSize;}
  const Vector<Box> &ProbDomain() const {return probDomain;}

protected:
  int finestLevel;
  std::vector<const BoxArray*> ba;
  Vector<int> refRatio;
  Vector<Real> probSize;
  Vector<Box> probDomain;
};

struct AMReXDataHierarchy
{
/*
  Data on a AMReXMeshHierarchy, currently pointing to MultiFabs of 
  named variables managed by an AmrData object.
*/
public:
  AMReXDataHierarchy(AmrData& ad, const Vector<std::string>& varNames) {
    mesh.define(ad);
    const Vector<std::string>& plotVarNames = ad.PlotVarNames();
    int nComp = varNames.size();
    int nlevs = mesh.FinestLevel() + 1;
    for (int i=0; i<nComp; ++i) {
      int idx = -1;
      for (int j=0; j<plotVarNames.size() && idx<0; ++j) {
        if (plotVarNames[j] == varNames[i]) {idx = j;}
      }
      if (ParallelDescriptor::IOProcessor() && idx<0) {
        Abort("Cannot find variable="+varNames[i]+" in pltfile");
      }
      std::vector<MultiFab*> mfs(nlevs);
      for (int lev=0; lev<nlevs; ++lev) {
        mfs[lev] = &ad.GetGrids(lev,idx); // Note: This lazily triggers a MultiFab read in the AmrData
      }
      varMap[varNames[i]] = std::make_pair(0,mfs);
    }
  }

  MultiFab &GetGrids(int level, const std::string& name) {
    if (varMap.find(name) == varMap.end()) {
      Abort("Unknown component requested");
    }
    return *(varMap[name].second[level]);
  }

  const AMReXMeshHierarchy& Mesh() const {return mesh;}

protected:
  AMReXMeshHierarchy mesh;
  std::map<std::string,std::pair<int,std::vector<MultiFab*>>> varMap;
};

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

    // Create the AmrData object from a pltfile on disk
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

    // Make a data struct for just the variables needed
    AMReXDataHierarchy data(amrData,varNames);
    const AMReXMeshHierarchy& mesh = data.Mesh();



    // Compute the volume integrals
    const int nGrow = 0;
    const int finestLevel = mesh.FinestLevel();
    const int nLev = finestLevel + 1;
    const int nComp = varNames.size();

    Vector<Real> integrals(nComp,0); // Results, initialized to zero
    for (int lev=0; lev<nLev; ++lev) {
      const BoxArray ba = mesh.boxArray(lev);

      // Make boxes that are projection of finer ones (if exist)
      const BoxArray baf = lev < finestLevel
                                 ? BoxArray(mesh.boxArray(lev+1)).coarsen(mesh.RefRatio()[lev])
                                 : BoxArray();

      // Compute volume of a cell at this level
      Real vol=1.0;
      for (int d=0; d<AMREX_SPACEDIM; ++d) {
        vol *= mesh.ProbSize()[d] / mesh.ProbDomain()[lev].length(d);
      }

      // For each component listed...
      for (int n=0; n<nComp; ++n) {

        // Make copy of original data because we will modify here
        const MultiFab &pfData = data.GetGrids(lev,varNames[n]);
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

          // Zero out covered cells
          if (lev < finestLevel) {
            std::vector< std::pair<int,Box> > isects = baf.intersections(box);                
            for (int ii = 0; ii < isects.size(); ii++) {
              myFab.setVal(0,isects[ii].second,0,1);
            }
          }

          // Get sum over tile, including zeroed covered cells
          sm += myFab.sum(box,0,1);
        }

        // Accumulate to this ranks total
        integrals[n] += sm * vol;

      }
    }

    // Accumulate over ranks
    ParallelDescriptor::ReduceRealSum(&(integrals[0]),nComp);

    // Integrals to stdout
    for (int n=0; n<nComp; ++n) {
      Print() << integrals[n] << " ";
    }
    Print() << std::endl;
  }
  amrex::Finalize();
  return 0;
}

