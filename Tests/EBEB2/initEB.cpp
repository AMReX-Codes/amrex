
#include <AMReX_GeometryShop.H>
#include <AMReX_SphereIF.H>
#include <AMReX_EBISLayout.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_EBTower.H>

#include "MyTest.H"
#include <AMReX_ParmParse.H>
#include <AMReX_VisMF.H>

#include <cmath>
#include <algorithm>

using namespace amrex;

void
MyTest::initializeEB ()
{
    BL_PROFILE("initializeEB");

    ParmParse pp("eb2");
    std::string geom_type;
    pp.get("geom_type", geom_type);

    if (geom_type == "sphere")
    {
        std::vector<Real> vc;
        pp.getarr("sphere_center", vc);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vc.size() >= AMREX_SPACEDIM,
                                         "eb2.sphere_center doesn't have enough items");
        
        Real radius;
        pp.get("sphere_radius", radius);

        bool has_fluid_inside;
        pp.get("sphere_has_fluid_inside", has_fluid_inside);

        SphereIF sif(radius, {AMREX_D_DECL(vc[0],vc[1],vc[2])}, has_fluid_inside);
        GeometryShop gshop(sif, false);
        const Real* dx = geom.CellSize();
        AMReX_EBIS::instance()->define(geom.Domain(), RealVect::Zero, dx[0], gshop,
                                       max_grid_size, max_coarsening_level);
    }
    else
    {
        amrex::Abort("geom_type "+geom_type+ " not supported in Mytest::initializeEB");
    }

    EBTower::Build();
    old_factory.reset(new EBFArrayBoxFactory(geom, grids, dmap, {1, 1, 1}, EBSupport::full));
}
