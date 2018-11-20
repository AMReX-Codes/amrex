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
            RealBox rb({-0.5,-0.5,-0.5}, {0.5,0.5,0.5}); // physical domain
            Array<int,AMREX_SPACEDIM> is_periodic{false, false, false};
            Geometry::Setup(&rb, 0, is_periodic.data());
            Box domain(IntVect(0), IntVect(n_cell-1));
            geom.define(domain);            
        }

        {
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
