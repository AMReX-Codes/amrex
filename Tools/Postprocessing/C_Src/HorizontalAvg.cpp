//
// Do planar averages of a list of variables, optionally weighted with one of the variables in the plotfile.
//

#include <fstream>

#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include "AMReX_ParmParse.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>
#include "AMReX_Utility.H"
#include <AMReX_VisMF.H>

using namespace amrex;

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        ParmParse pp;

        FArrayBox::setFormat(FABio::FAB_IEEE_32);
    
        int verbose; pp.query("verbose",verbose);
        if (verbose > 2)
            AmrData::SetVerbose(true);
    
        std::string infile; pp.get("infile",infile);
        std::string outdir = infile + std::string("_stats"); pp.query("outdir",outdir);
        std::string outfile("horizAvg.dat"); pp.query("outfile",outfile);

        DataServices::SetBatchMode();
        Amrvis::FileType fileType(Amrvis::NEWPLT);
        DataServices dataServices(infile, fileType);
        if (!dataServices.AmrDataOk())
            DataServices::Dispatch(DataServices::ExitRequest, NULL);    
        AmrData& amrData = dataServices.AmrDataRef();

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
            for (int i=0; i<nComp; ++i)
                comps[i] = sComp + i;
        }

        int weightComp = -1; pp.query("weightComp",weightComp);
        int localWtComp = -1;
        if (weightComp>=0) {
            for (int i=0; i<comps.size(); ++i) {
                if (comps[i]==weightComp) localWtComp = i;
            }
            if (localWtComp<0) {
                localWtComp = comps.size();
                comps.resize(localWtComp+1);
                comps[localWtComp] = weightComp;
            }
        }
    
        int nComp = comps.size();
        Vector<string> varNames(nComp);
        for (int i=0; i<nComp; ++i)
            varNames[i] = amrData.PlotVarNames()[comps[i]];

        Vector<int> destFillComps(nComp);
        for (int i=0; i<nComp; ++i)
            destFillComps[i] = i;

        int finestLevel = 0; pp.query("finestLevel",finestLevel);
        int max_grid_size = 16; pp.query("max_grid_size",max_grid_size);

        Box domain = amrData.ProbDomain()[finestLevel];

        Vector<Real> dx(BL_SPACEDIM);
        Vector<Real> plo = amrData.ProbLo();
        for (int j=0; j<BL_SPACEDIM; ++j) {
            dx[j] = amrData.ProbSize()[j]/domain.length(j);
        }

        BoxArray ba = BoxArray(domain).maxSize(max_grid_size);
        DistributionMapping dmap(ba);
        MultiFab mf(ba,dmap,nComp,0);
        amrData.FillVar(mf,finestLevel,varNames,destFillComps);

        if (ParallelDescriptor::IOProcessor() && verbose)
            std::cerr << "Data has been read" << std::endl;
    
        Box probDomain = amrData.ProbDomain()[finestLevel];
        int dir = BL_SPACEDIM-1;
        int ksize(probDomain.length(dir));
        Vector<Real> havg(ksize*nComp,0.);
        Vector<Real> wavg(ksize,0.);
        Vector<Real> area(ksize,0.);
    
        FArrayBox mask;
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.validbox();
            FArrayBox& fab = mf[mfi];
            mask.resize(box,1); mask.setVal(1);
        
            if (weightComp>=0) {
                for (int n=0; n<nComp; ++n) {
                    if (n!=localWtComp) {
                        fab.mult(fab,localWtComp,n,1);
                    }
                }
            }

            for (int n=0; n<nComp; ++n) {
                fab.mult(mask,0,n,1);
            }

            const IntVect& se = box.smallEnd();
            const IntVect& be = box.bigEnd();
            for (int k=se[dir]; k<=be[dir]; ++k)
            {
                Box sub(box); sub.setSmall(dir,k); sub.setBig(dir,k);  // slab box on this fab
                area[k] += mask.sum(sub,0,1);
                if (weightComp>=0) {
                    wavg[k] += fab.sum(sub,localWtComp,1);
                }
                for (int n=0; n<nComp; ++n) {
                    havg[k*nComp+n] += fab.sum(sub,n,1);
                }
            }        
        }

        if (ParallelDescriptor::IOProcessor() && verbose>1)
            std::cerr << "Finished grinding data" << '\n';

        ParallelDescriptor::ReduceRealSum(havg.dataPtr(), havg.size(), ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::ReduceRealSum(area.dataPtr(), area.size(), ParallelDescriptor::IOProcessorNumber());
        if (weightComp>=0) {
            ParallelDescriptor::ReduceRealSum(wavg.dataPtr(), wavg.size(), ParallelDescriptor::IOProcessorNumber());
        }

        if (ParallelDescriptor::IOProcessor())
        {
            if (weightComp>=0) {
                for (int k=0; k<ksize; ++k) {
                    wavg[k] /= area[k];
                    for (int n=0; n<nComp; ++n) {
                        if (n!=localWtComp) {
                            havg[k*nComp+n] /= wavg[k];
                        }
                    }
                }
            }
        
            if (!amrex::UtilCreateDirectory(outdir, 0755))
                amrex::CreateDirectoryFailed(outdir);

            std::string outfileFULL = outdir + '/' + outfile;
            if (verbose) {
                std::cerr << "Writing averages to " << outfileFULL << '\n';
            }

            std::ofstream os;
            os.open(outfileFULL.c_str());

            bool writeTecHdr = true; pp.query("writeTecHdr",writeTecHdr);
            if (writeTecHdr) {
                os << "VARIABLES = \"x\" ";
                for (int i=0; i<nComp; ++i) {
                    os << "\"" << amrData.PlotVarNames()[comps[i]] << "\" ";
                }
                os << "\nZONE F=POINT I=" << ksize << '\n';
            }

            for (int k=0; k<ksize; ++k)
            {
                os << (k+0.5)*dx[dir] << " ";
                for (int n=0; n<nComp; ++n)
                    os << havg[k*nComp+n]/area[k] << " ";
                os << '\n';
            }
            os.close();
        }
    }
    amrex::Finalize();
    return 0;
}
