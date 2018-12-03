#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

int main (int argc, char* argv[])
{
  amrex::Initialize(argc, argv);

  {
    int n_cell = 128;
    int max_grid_size = 32;

    // read parameters
    {
      ParmParse pp;
      pp.query("n_cell", n_cell);
      pp.query("max_grid_size", max_grid_size);
    }

    Geometry geom;
    {
      RealBox rb({-5.2,-5.2,-1.0}, {5.2,5.2,2.0}); // physical domain
      Array<int,AMREX_SPACEDIM> is_periodic{false, false, false};
      Geometry::Setup(&rb, 0, is_periodic.data());
      Box domain(IntVect(0), IntVect(n_cell-1));
      geom.define(domain);            
    }

    {
      std::vector<amrex::RealVect> splpts;
      amrex::RealVect p;
      p = amrex::RealVect(D_DECL(36.193*0.1, 7.8583*0.1, 0.0));
      splpts.push_back(p);
      p = amrex::RealVect(D_DECL(35.924*0.1, 7.7881*0.1, 0.0));
      splpts.push_back(p);
      p = amrex::RealVect(D_DECL(35.713*0.1, 7.5773*0.1, 0.0));
      splpts.push_back(p);
      p = amrex::RealVect(D_DECL(35.643*0.1, 7.3083*0.1, 0.0));
      splpts.push_back(p);
      p = amrex::RealVect(D_DECL(35.3*0.1, 7.0281*0.1, 0.0));
      splpts.push_back(p);
      p = amrex::RealVect(D_DECL(35.421*0.1, 6.241*0.1, 0.0));
      splpts.push_back(p);
      p = amrex::RealVect(D_DECL(34.82*0.1, 5.686*0.1, 0.0));
      splpts.push_back(p);
      p = amrex::RealVect(D_DECL(30.539*0.1, 3.5043*0.1, 0.0));
      splpts.push_back(p);
      p = amrex::RealVect(D_DECL(29.677*0.1, 2.6577*0.1, 0.0));
      splpts.push_back(p);
      p = amrex::RealVect(D_DECL(29.457*0.1, 1.47*0.1, 0.0));
      splpts.push_back(p);
      // p = amrex::RealVect(D_DECL(29.38*0.1, -1.1038*0.1, 0.0));
      // splpts.push_back(p);
      // p = amrex::RealVect(D_DECL(29.3*0.1, -2.7262*0.1, 0.0));
      // splpts.push_back(p);
      // p = amrex::RealVect(D_DECL(29.273*0.1, -4.3428*0.1, 0.0));
      // splpts.push_back(p);
      p = amrex::RealVect(D_DECL(28.364*0.1, -5.7632*0.1, 0.0));
      splpts.push_back(p);
      p = amrex::RealVect(D_DECL(27.151*0.1, -6.8407*0.1, 0.0));
      splpts.push_back(p);
      p = amrex::RealVect(D_DECL(25.694*0.1, -7.5555*0.1, 0.0));
      splpts.push_back(p);
      p = amrex::RealVect(D_DECL(24.035*0.1, -7.8586*0.1, 0.0));
      splpts.push_back(p);
      p = amrex::RealVect(D_DECL(22.358*0.1, -7.6902*0.1, 0.0));
      splpts.push_back(p);
      EB2::SplineIF Piston;
      Piston.addSplineElement(splpts);
      std::vector<amrex::RealVect> lnpts;

      p = amrex::RealVect(D_DECL(22.358*0.1, -7.6902*0.1, 0.0));
      lnpts.push_back(p);
      p = amrex::RealVect(D_DECL(1.9934*0.1, 3.464*0.1, 0.0));
      lnpts.push_back(p);
      p = amrex::RealVect(D_DECL(0.0, 3.464*0.1, 0.0));
      lnpts.push_back(p);
      Piston.addLineElement(lnpts);
      lnpts.clear();

      p = amrex::RealVect(D_DECL(49.0*0.1, 7.8583*0.1,  0.0));
      lnpts.push_back(p);
      p = amrex::RealVect(D_DECL(36.193*0.1, 7.8583*0.1, 0.0));
      lnpts.push_back(p);
      Piston.addLineElement(lnpts);
      lnpts.clear();

      EB2::CylinderIF cylinder(4.80, 7.0, 2, {0.0, 0.0, -1.0}, false);

      auto revolvePiston  = EB2::lathe(Piston);
      auto PistonComplement = EB2::makeComplement(revolvePiston);
      auto PistonCylinder = EB2::makeIntersection(PistonComplement, cylinder);
      auto gshop = EB2::makeShop(PistonCylinder);
      EB2::Build(gshop, geom, 0, 0);
    }

    MultiFab mf;
    {
      BoxArray ba(geom.Domain());
      ba.maxSize(max_grid_size);
      DistributionMapping dm{ba};

      std::unique_ptr<EBFArrayBoxFactory> factory
          = amrex::makeEBFabFactory(geom, ba, dm, {2,2,2}, EBSupport::full);

      mf.define(ba, dm, 1, 0, MFInfo(), *factory);
      mf.setVal(1.0);
    }

    EB_WriteSingleLevelPlotfile("plt", mf, {"rho"}, geom, 0.0, 0);
  }

  amrex::Finalize();
}
