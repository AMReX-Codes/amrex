
#include <AMReX_GeometryShop.H>
#include <AMReX_SphereIF.H>
#include <AMReX_PlaneIF.H>
#include <AMReX_UnionIF.H>
#include <AMReX_TransformIF.H>
#include <AMReX_ComplementIF.H>
#include <AMReX_IntersectionIF.H>
#include <AMReX_LatheIF.H>
#include <AMReX_PolynomialIF.H>
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
    else if (geom_type == "box")
    {
        std::vector<Real> vlo;
        pp.getarr("box_lo", vlo);
        RealArray lo{AMREX_D_DECL(vlo[0],vlo[1],vlo[2])};

        std::vector<Real> vhi;
        pp.getarr("box_hi", vhi);
        RealArray hi{AMREX_D_DECL(vhi[0],vhi[1],vhi[2])};

        bool has_fluid_inside;
        pp.get("box_has_fluid_inside", has_fluid_inside);
        
        Vector<BaseIF*> planes;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
//            planes.push_back(new PlaneIF(BASISREALV(idim));
            
        }

    }
    else if (geom_type == "plane")
    {
        std::vector<Real> vpoint;
        pp.getarr("plane_point", vpoint);

        std::vector<Real> vnormal;
        pp.getarr("plane_normal", vnormal);

        PlaneIF pif({AMREX_D_DECL(vnormal[0],vnormal[1],vnormal[2])},
                    {AMREX_D_DECL(vpoint[0],vpoint[1],vpoint[2])}, false);
        GeometryShop gshop(pif, false);
        const Real* dx = geom.CellSize();
        AMReX_EBIS::instance()->define(geom.Domain(), RealVect::Zero, dx[0], gshop,
                                       max_grid_size, max_coarsening_level);
    }
    else
    {
        amrex::Abort("geom_type "+geom_type+ " not supported in Mytest::initializeEB");
    }

    EBTower::Build();
    old_factory.emplace_back(new EBFArrayBoxFactory(geom, grids, dmap, {ng_eb, ng_eb, ng_eb}, EBSupport::full));

    for (int ilev = 1; ilev <= max_coarsening_level; ++ilev)
    {
        int ratio = std::pow(2,ilev);
        Geometry cgeom(amrex::coarsen(geom.Domain(),ratio));
        old_factory.emplace_back(new EBFArrayBoxFactory(cgeom, amrex::coarsen(grids,ratio),
                                                        dmap, {ng_eb, ng_eb, ng_eb}, EBSupport::full));
    }
}
