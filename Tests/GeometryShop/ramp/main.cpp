
#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include <AMReX_EBLevelGrid.H>
#include <AMReX_GeometryShop.H>
#include <AMReX_SphereIF.H>
#include <AMReX_PlaneIF.H>
#include <AMReX_AllRegularService.H>
#include <AMReX_FlatPlateGeom.H>
#include <AMReX_EBISLayout.H>
#include <AMReX_EBGraph.H>
#include <AMReX_EBDebugOut.H>
#include <AMReX_EBCellFAB.H>
#include <AMReX_EBCellFactory.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_UnionIF.H>
#include <AMReX_TransformIF.H>
#include <AMReX_ComplementIF.H>
#include <AMReX_IntersectionIF.H>
#include <AMReX_LatheIF.H>


using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    int max_grid_size;
    vector<int> n_cell;

    // initialize eb
    {
        ParmParse pp;
        ParmParse ppa("amr");
        ParmParse ppg("geometry");

        int max_level;
        ppa.get("max_level", max_level);

        ppa.getarr("n_cell", n_cell, 0, SpaceDim);

        ppa.get("max_grid_size", max_grid_size);

        vector<int> ref_ratio;
        ppa.getarr("ref_ratio", ref_ratio, 0, max_level);
      
        IntVect ivlo = IntVect::Zero;
        IntVect ivhi;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          ivhi[idir] = n_cell[idir]-1;
        }
        Box finest_domain(ivlo, ivhi);
        // < maxlev because there is one less refinement than number of levels
        for(int ilev = 0; ilev < max_level; ilev++)
        {
          finest_domain.refine(ref_ratio[ilev]);
        }
      
        vector<Real> prob_lo, prob_hi;
        ppg.getarr("prob_lo", prob_lo, 0, SpaceDim);
        ppg.getarr("prob_hi", prob_hi, 0, SpaceDim);
        Real fine_dx = (prob_hi[0]-prob_lo[0])/finest_domain.size()[0];
        RealVect origin = RealVect::Zero;
  

        amrex::Print() << "defining EBIS with..." << endl;
        amrex::Print() << "finest domain = " << finest_domain << endl;
        amrex::Print() << "finest dx     = " << fine_dx       << endl;

        
        std::string geom_type;
        pp.get("geom_type", geom_type);

        if (geom_type != "ramp") {
            amrex::Abort("must be ramp for this test");
        } else {
          std::unique_ptr<BaseIF> impfunc;
          
          amrex::Print() << "ramp geometry\n";
          int upDir;
          int indepVar;
          Real startPt;
          Real slope;
          pp.get("up_dir",upDir);
          pp.get("indep_var",indepVar);
          pp.get("start_pt", startPt);
          pp.get("ramp_slope", slope);
          
          RealVect normal = RealVect::Zero;
          normal[upDir] = 1.0;
          normal[indepVar] = -slope;
          
          RealVect point = RealVect::Zero;
          point[upDir] = -slope*startPt;
          
          bool normalInside = true;
          
          impfunc.reset(static_cast<BaseIF*>(new PlaneIF(normal,point,normalInside)));

          bool eb_verbosity = false;
          pp.query("eb_verbosity", eb_verbosity);
          GeometryShop gshop(*impfunc, eb_verbosity);
          AMReX_EBIS::instance()->define(finest_domain, origin, fine_dx, gshop, max_grid_size, max_level);
        }

        Box geom_domain{Box{IntVect{0,0,0}, IntVect{n_cell[0]-1,n_cell[1]-1,n_cell[2]-1}}};
        BoxArray ba{geom_domain};
        ba.maxSize(max_grid_size);
        DistributionMapping dm{ba};
        const int ng = 0;

        const EBLevelGrid levelgrid(ba,dm,geom_domain,ng);
        const EBISLayout& ebisl = levelgrid.getEBISL();

        //        IntVect debugcell(945,137,7);
        IntVect debugcell(994,213,7);

        for (MFIter mfi(ba,dm); mfi.isValid(); ++mfi)
        {
          const EBISBox& ebisbox = ebisl[mfi];
          const Box& bx = ba[mfi];
            
          if (bx.contains(debugcell))
          {
            int ireg = 0;
            int icov = 0;
            if(ebisbox.isRegular(debugcell))
            {
              ireg = 1;
            }
            if(ebisbox.isCovered(debugcell))
            {
              icov = 1;
            }
            amrex::AllPrint() << "Box " << bx << " on Proc. " << ParallelDescriptor::MyProc()
                              << " contains Cell " << debugcell << "\n"
                              << "  Is regular? " <<  ireg << "\n"
                              << "  Is covered? " <<  icov << "\n";
          }
        }
    }

    amrex::Finalize();
}
