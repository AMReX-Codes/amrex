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
        int which_geom = 0;

        // read parameters
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("which_geom", which_geom);
        }

        Geometry geom;
        {
            RealBox rb({-1.3,-1.3,-0.475}, {1.3,1.3,0.5}); // physical domain
            // RealBox rb({-5.2,-5.2,-1.9}, {5.2,5.2,2.0}); // physical domain
            Array<int,AMREX_SPACEDIM> is_periodic{false, false, false};
            Geometry::Setup(&rb, 0, is_periodic.data());
            Box domain(IntVect(0), IntVect(n_cell-1));
            geom.define(domain);            
        }

        if (which_geom == 0) {
            EB2::SphereIF sphere(0.5, {0.0,0.0,0.0}, false);
            EB2::BoxIF cube({-0.4,-0.4,-0.4}, {0.4,0.4,0.4}, false);
            auto cubesphere = EB2::makeIntersection(sphere, cube);
            
            EB2::CylinderIF cylinder_x(0.25, 0, {0.0,0.0,0.0}, false);
            EB2::CylinderIF cylinder_y(0.25, 1, {0.0,0.0,0.0}, false);
            EB2::CylinderIF cylinder_z(0.25, 2, {0.0,0.0,0.0}, false);
            auto three_cylinders = EB2::makeUnion(cylinder_x, cylinder_y, cylinder_z);

            auto csg = EB2::makeDifference(cubesphere, three_cylinders);

            auto gshop = EB2::makeShop(csg);
            EB2::Build(gshop, geom, 0, 0);
        } else if (which_geom == 1) {
          // This geometry uses splines to put an interesting 'piston bowl' shape at the
          // bottom of the cylinder
          Real scaleFact;
          scaleFact = 0.25;

          std::vector<amrex::RealVect> splpts;
          amrex::RealVect p;
          p = amrex::RealVect(D_DECL(36.193*0.1*scaleFact, 7.8583*0.1*scaleFact, 0.0));
          splpts.push_back(p);
          p = amrex::RealVect(D_DECL(35.924*0.1*scaleFact, 7.7881*0.1*scaleFact, 0.0));
          splpts.push_back(p);
          p = amrex::RealVect(D_DECL(35.713*0.1*scaleFact, 7.5773*0.1*scaleFact, 0.0));
          splpts.push_back(p);
          p = amrex::RealVect(D_DECL(35.643*0.1*scaleFact, 7.3083*0.1*scaleFact, 0.0));
          splpts.push_back(p);
          p = amrex::RealVect(D_DECL(35.3*0.1*scaleFact, 7.0281*0.1*scaleFact, 0.0));
          splpts.push_back(p);
          p = amrex::RealVect(D_DECL(35.421*0.1*scaleFact, 6.241*0.1*scaleFact, 0.0));
          splpts.push_back(p);
          p = amrex::RealVect(D_DECL(34.82*0.1*scaleFact, 5.686*0.1*scaleFact, 0.0));
          splpts.push_back(p);
          p = amrex::RealVect(D_DECL(30.539*0.1*scaleFact, 3.5043*0.1*scaleFact, 0.0));
          splpts.push_back(p);
          p = amrex::RealVect(D_DECL(29.677*0.1*scaleFact, 2.6577*0.1*scaleFact, 0.0));
          splpts.push_back(p);
          p = amrex::RealVect(D_DECL(29.457*0.1*scaleFact, 1.47*0.1*scaleFact, 0.0));
          splpts.push_back(p);
          p = amrex::RealVect(D_DECL(28.364*0.1*scaleFact, -5.7632*0.1*scaleFact, 0.0));
          splpts.push_back(p);
          p = amrex::RealVect(D_DECL(27.151*0.1*scaleFact, -6.8407*0.1*scaleFact, 0.0));
          splpts.push_back(p);
          p = amrex::RealVect(D_DECL(25.694*0.1*scaleFact, -7.5555*0.1*scaleFact, 0.0));
          splpts.push_back(p);
          p = amrex::RealVect(D_DECL(24.035*0.1*scaleFact, -7.8586*0.1*scaleFact, 0.0));
          splpts.push_back(p);
          p = amrex::RealVect(D_DECL(22.358*0.1*scaleFact, -7.6902*0.1*scaleFact, 0.0));
          splpts.push_back(p);
          EB2::SplineIF Piston;
          Piston.addSplineElement(splpts);
          
          std::vector<amrex::RealVect> lnpts;

          p = amrex::RealVect(D_DECL(22.358*0.1*scaleFact, -7.6902*0.1*scaleFact, 0.0));
          lnpts.push_back(p);
          p = amrex::RealVect(D_DECL(1.9934*0.1*scaleFact, 3.464*0.1*scaleFact, 0.0));
          lnpts.push_back(p);
          p = amrex::RealVect(D_DECL(0.0, 3.464*0.1*scaleFact, 0.0));
          lnpts.push_back(p);
          Piston.addLineElement(lnpts);
          lnpts.clear();

          p = amrex::RealVect(D_DECL(49.0*0.1*scaleFact, 7.8583*0.1*scaleFact,  0.0));
          lnpts.push_back(p);
          p = amrex::RealVect(D_DECL(36.193*0.1*scaleFact, 7.8583*0.1*scaleFact, 0.0));
          lnpts.push_back(p);
          Piston.addLineElement(lnpts);
          lnpts.clear();

          EB2::CylinderIF cylinder(4.80*scaleFact, 7.0*scaleFact, 2, {0.0, 0.0, -1.0*scaleFact}, true);

          auto revolvePiston  = EB2::lathe(Piston);
          auto PistonCylinder = EB2::makeUnion(revolvePiston, cylinder);
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
