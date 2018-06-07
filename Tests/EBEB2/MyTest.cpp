
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB2_IF_Lathe.H>
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Sphere.H>
#include <AMReX_EB2_IF_Plane.H>

#include <AMReX_EBIndexSpace.H>
#include <AMReX_EBTower.H>

#include "MyTest.H"
#include <AMReX_ParmParse.H>
#include <AMReX_VisMF.H>

#include <cmath>
#include <algorithm>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();

    initGrids();

//    initializeEB();

    initializeEB2();

    initData();
}

MyTest::~MyTest ()
{
    EB2::Finalize();
    AMReX_EBIS::reset();
    EBTower::Destroy();
}

void
MyTest::readParameters ()
{
    ParmParse pp;
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);
    pp.query("max_coarsening_level", max_coarsening_level);
    max_coarsening_level = std::max(max_coarsening_level, 0);
}

void
MyTest::initGrids ()
{
    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
    std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
    Geometry::Setup(&rb, 0, is_periodic.data());
    Box domain0(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});

    geom.define(domain0);

    grids.define(domain0);
    grids.maxSize(max_grid_size);

    dmap.define(grids);
}

void
MyTest::initData ()
{
}

void
MyTest::test ()
{
}

void
MyTest::initializeEB2 ()
{
    BL_PROFILE("initializeEB2");

    ParmParse pp("eb2");
    std::string geom_type;
    pp.get("geom_type", geom_type);
    
    EB2::Info info;
    info.setMaxCoarseningLevel(max_coarsening_level)
        .setMaxGridSize(max_grid_size);

    if (geom_type == "combustor")
    {
        EB2::PlaneIF farwall({AMREX_D_DECL(0.45,0.,0.)},
                             {AMREX_D_DECL(1.  ,0.,0.)});
        auto ramp = EB2::makeIntersection(EB2::PlaneIF({AMREX_D_DECL(0.25, 0.75, 0.)},
                                                       {AMREX_D_DECL(0.  , -1. , 0.)}),
                                          EB2::PlaneIF({AMREX_D_DECL(0.25, 0.75, 0.)},
                                                       {AMREX_D_DECL(.69 , -.19, 0.)}),
                                          EB2::PlaneIF({AMREX_D_DECL(0.06, 0.  , 0.)},
                                                       {AMREX_D_DECL(1.  , 0.  , 0.)}));
        EB2::BoxIF pipe({AMREX_D_DECL(0.06, -1.0, -100.0)},
                        {AMREX_D_DECL(0.08, 0.5, 100.0)}, false);
        EB2::BoxIF flat_corner({AMREX_D_DECL(0.05999, -1.0, -100.0)},
                               {AMREX_D_DECL(1.0, 0.26, 100.0)}, false);
        auto polys = EB2::makeUnion(farwall, ramp, pipe, flat_corner);
        auto pr = EB2::translate(EB2::lathe(polys), {AMREX_D_DECL(0.5,0.5,0.)});

        auto gshop = EB2::makeShop(pr);
        EB2::Initialize(gshop, geom, info);
    }
    else
    {
        EB2::Initialize(geom, info);
    }

    MultiFab vfrc(grids, dmap, 1, 1);
    const EB2::Level& eb2_level = EB2::getLevel(geom);
    eb2_level.fillVolFrac(vfrc, geom);

    VisMF::Write(vfrc, "vfrc-new");
}

