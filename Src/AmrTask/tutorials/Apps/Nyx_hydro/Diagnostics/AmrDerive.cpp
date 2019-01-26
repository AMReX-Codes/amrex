//BL_COPYRIGHT_NOTICE

#include <new>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <list>
#include <map>
#include <algorithm>

using std::cout;
using std::cerr;
using std::endl;

#ifndef WIN32
#include <unistd.h>
#endif

#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_DataServices.H"
#include "AMReX_Utility.H"
#include "AMReX_FArrayBox.H"

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
    cout << '\n';
    cout << "Usage:" << '\n';
    cout << "    infile  = plotfilename" << '\n';
    cout << "   [-verbose]" << '\n';
    cout << '\n';
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

    bool verbose = false;
    if (pp.contains("verbose"))
    {
        verbose = true;
        AmrData::SetVerbose(true);
    }

    std::string infile; 
    pp.get("infile",infile);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "Reading " << infile << "...";

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if( ! dataServices.AmrDataOk())  
        DataServices::Dispatch(DataServices::ExitRequest, NULL);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "Done reading plot file" << std::endl;

    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel); finestLevel=std::min(finestLevel,amrData.FinestLevel());
    int Nlev = finestLevel + 1;

    if (ParallelDescriptor::IOProcessor())
        std::cout << "... finest level " << finestLevel << std::endl;

    Vector<int> comps;
    if (int nc = pp.countval("comps"))
    {
        comps.resize(nc);
        pp.getarr("comps",comps,0,nc);
    }
    else
    {
        int sComp = 0;
        pp.query("sComp",sComp);
        int nComp = 1;
        pp.query("nComp",nComp);
        BL_ASSERT(sComp+nComp <= amrData.NComp());
        comps.resize(nComp);
        std::cout << "NCOMP NOW " << nComp << std::endl;
        for (int i=0; i<nComp; ++i)
            comps[i] = sComp + i;
    }

    int nComp = comps.size();
    const Vector<string>& plotVarNames=amrData.PlotVarNames();
    Vector<string> inVarNames(nComp);
    Vector<int> destFillComps(nComp);
    for (int i=0; i<nComp; ++i)
    {
        inVarNames[i] = plotVarNames[comps[i]];
        std::cout << "plotVarNames " << plotVarNames[comps[i]] << std::endl;;
        destFillComps[i] = i;
    }

    const int nGrow = 0;

    if (ParallelDescriptor::IOProcessor() && verbose>0)
    {
        cerr << "Will read the following states: ";
        for (int i=0; i<nComp; ++i)
            cerr << " " << amrData.StateNumber(inVarNames[i]) << " (" << inVarNames[i] << ")" ;
        cerr << '\n';
    }

    const Box& probDomain = amrData.ProbDomain()[finestLevel];
    int dir=BL_SPACEDIM-1; pp.query("dir",dir);
    const IntVect lo=probDomain.smallEnd();
    IntVect hi=lo; hi[dir] = probDomain.bigEnd()[dir];
    const Box resBox(lo,hi);
    FArrayBox resFab(resBox,nComp); resFab.setVal(0.);

    int accumFac = 1;
    for (int lev=finestLevel; lev>=0; --lev)
    {
        const BoxArray& ba = amrData.boxArray(lev);
	const DistributionMapping& dm = amrData.DistributionMap(lev);
        MultiFab mf(ba,dm,nComp,nGrow);
        if (ParallelDescriptor::IOProcessor() && verbose>0)
            cerr << "...filling data at level " << lev << endl;
        amrData.FillVar(mf,lev,inVarNames,destFillComps);

        // Zero covered regions
        if (lev < finestLevel)
        {
            const BoxArray baf = BoxArray(amrData.boxArray(lev+1)).coarsen(amrData.RefRatio()[lev]);

	    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
		FArrayBox& myFab = mf[mfi];
                std::vector< std::pair<int,Box> > isects = baf.intersections(ba[mfi.index()]);
                
		for (int ii = 0; ii < isects.size(); ii++)
		    myFab.setVal(0,isects[ii].second,0,nComp);
            }
        }

        // Scale by "volume" (plane, 1 cell thick) of cells, normalized to volume of finest cells
        mf.mult(accumFac*accumFac,0,nComp);

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const FArrayBox& fab = mf[mfi];
            const Box& box = mfi.validbox();
            IntVect ivlo = box.smallEnd();
            IntVect ivhi = box.bigEnd(); ivhi[dir] = ivlo[dir];

            for (int plane=box.smallEnd(dir); plane<=box.bigEnd(dir); ++plane)
            {
                ivlo[dir] = plane;
                ivhi[dir] = plane;

                Box subbox(ivlo,ivhi);
                Vector<Real> thisSum(nComp);
                for (int n=0; n<nComp; ++n)
                    thisSum[n] = fab.sum(subbox,n,1);

                // Now, increment each sum affected by this coarse plane
                for (int r=0; r<accumFac; ++r)
                {
                    const int finePlane = plane*accumFac + r;
                    IntVect ivF=resBox.smallEnd(); ivF[dir] = finePlane;
                    for (int n=0; n<nComp; ++n)
                        resFab(ivF,n) += thisSum[n];
                }
            }
        }

        if (lev>0)
            accumFac *= amrData.RefRatio()[lev-1];
    }

    // Accumulate sums from all processors to IOProc
    ParallelDescriptor::ReduceRealSum(resFab.dataPtr(),nComp*resBox.numPts(),ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
    {
        // Scale result by total "volume" (plane, 1 cell thick) of slab at each plane
        Real volume=1;
        for (int d=0; d<BL_SPACEDIM; ++d)
            if (d!=dir)
                volume*=probDomain.length(d);

        resFab.mult(1/volume);
        std::cout << "RESFAB " << resFab << std::endl;
    }

    
    amrex::Finalize();
    return 0;
}

