#include "AMReX_DataServices.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;
using namespace std;

int main ( int argc, char* argv[] )
{
    amrex::Initialize (argc, argv, false);
    
    {
        if (argc < 3)
        {
            amrex::Print() << "Usage: " << argv[0] << " in_file out_file \n";
            amrex::Abort();
        }

        const string in_file = argv[1];
        const string out_file = argv[2];

        const int factor = 2;
        const int lev = 0;
        
        DataServices::SetBatchMode();

        Amrvis::FileType fileType(Amrvis::NEWPLT);
        DataServices dataServices(in_file, fileType);        
        
        if ( !dataServices.AmrDataOk() )
            DataServices::Dispatch(DataServices::ExitRequest, NULL);
        
        AmrData& data = dataServices.AmrDataRef();
        
        int finestLevel = data.FinestLevel();
        AMREX_ALWAYS_ASSERT(finestLevel == 0);
        
        const int ng = data.NGrow();
        AMREX_ALWAYS_ASSERT(ng == 0);
        
        const int ncomp = data.NComp();
        Vector<int> destFillComps;
        for (int i = 0; i < ncomp; ++i)
            destFillComps.push_back(i);
        
        const Vector<Box>& bb = data.ProbDomain();
        const Vector<Real>& plo = data.ProbLo();
        const Vector<Real>& phi = data.ProbHi();
        
        Vector<Real> extracted_phi;
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            Real psize = phi[i] - plo[i];
            Real extracted_psize = psize / factor;
            extracted_phi.push_back(plo[i] + extracted_psize);
        }
        
        const Box extracted_bb(bb[0].smallEnd(), bb[0].smallEnd() + bb[0].size()/factor - 1);
        const RealBox extracted_rb(plo.data(), extracted_phi.data());
        Geometry extracted_geom(extracted_bb, &extracted_rb);
        
        const BoxArray& ba = data.boxArray(lev);
        const IntVect box_size = ba[0].size();
        const DistributionMapping& dm = data.DistributionMap(lev);
        
        BoxArray extracted_ba(extracted_bb);
        extracted_ba.maxSize(box_size);
        DistributionMapping extracted_dm(extracted_ba);
        
        MultiFab extracted_mf(extracted_ba, extracted_dm, ncomp, ng);
        MultiFab mf(ba, dm, ncomp, ng);

        const Vector<string>& varNames = data.PlotVarNames();
        for (int nv = 0; nv < ncomp; ++nv)
        {
            data.FillVar (mf, lev, varNames, destFillComps);
        }
        
        extracted_mf.copy(mf, 0, 0, ncomp);
        
        WriteSingleLevelPlotfile (out_file, extracted_mf, varNames, extracted_geom, 0.0, 0);        
    }
    
    amrex::Finalize ();
}
