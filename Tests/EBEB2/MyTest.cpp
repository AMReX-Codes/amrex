
#include <AMReX_GeometryShop.H>
#include <AMReX_SphereIF.H>
#include <AMReX_EBISLayout.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_EBTower.H>

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

#include "MyTest.H"
#include <AMReX_ParmParse.H>
#include <AMReX_VisMF.H>

#include <cmath>

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
    info.setMaxCoarseningLevel(0)
        .setMaxGridSize(max_grid_size);

    if (geom_type == "combustor")
    {
        EB2::PlaneIF farwall({0.45,0.,0.}, {1.,0.,0.});
        auto ramp = EB2::makeIntersection(EB2::PlaneIF({0.25, 0.75, 0.}, {0. , -1. , 0.}),
                                          EB2::PlaneIF({0.25, 0.75, 0.}, {.69, -.19, 0.}),
                                          EB2::PlaneIF({0.06, 0.  , 0.}, {1. ,  0. , 0.}));
        EB2::BoxIF pipe({0.06, -1.0, -100.0}, {0.08, 0.5, 100.0}, false);
        EB2::BoxIF flat_corner({0.05999, -1.0, -100.0}, {1.0, 0.26, 100.0}, false);
        auto polys = EB2::makeUnion(farwall, ramp, pipe, flat_corner);
        auto pr = EB2::translate(EB2::lathe(polys), {0.5,0.5,0.});

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
                                       max_grid_size, 0);
    }
    else
    {
        return;
        amrex::Abort("geom_type "+geom_type+ " not supported in Mytest::initializeEB");
    }


    EBTower::Build();
    EBFArrayBoxFactory factory(geom, grids, dmap, {1, 1, 1}, EBSupport::full);

    const MultiFab& vfrc = factory.getVolFrac();
    VisMF::Write(vfrc, "vfrc-old");
}
