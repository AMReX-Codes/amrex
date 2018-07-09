
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Scale.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB2_IF_Rotation.H>
#include <AMReX_EB2_IF_Lathe.H>
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Ellipsoid.H>
#include <AMReX_EB2_IF_Sphere.H>
#include <AMReX_EB2_IF_Plane.H>

#include <AMReX_EBIndexSpace.H>
#include <AMReX_EBTower.H>

#include "MyTest.H"
#include <AMReX_ParmParse.H>
#include <AMReX_VisMF.H>
#include <AMReX_MultiFabUtil.H>

#include <cmath>
#include <algorithm>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();

    initGrids();

    if (test_old_eb) {
        initializeEB();
    }

    initializeEB2();
}

MyTest::~MyTest ()
{
    old_factory.clear();
    new_factory.clear();
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
    pp.query("test_old_eb", test_old_eb);
    pp.query("test_cellflag", test_cellflag);
    pp.query("ng_eb", ng_eb);
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
MyTest::test ()
{
    for (int ilev = 0, nlevs = new_factory.size(); ilev < nlevs; ++ilev)
    {
        amrex::Print() << "\n Testing level " << ilev << "\n";

        const Box& lev_domain_m1 = amrex::grow(new_factory[ilev]->getDomain(), -1);

        const FabArray<EBCellFlagFab>& cellflag_new = new_factory[ilev]->getMultiEBCellFlagFab();

        const MultiFab& vfrc_new = new_factory[ilev]->getVolFrac();
        VisMF::Write(vfrc_new, "new-vfrc-lev"+std::to_string(ilev));
        
        const MultiCutFab& cent_new = new_factory[ilev]->getCentroid();
        VisMF::Write(cent_new.ToMultiFab(0.,0.), "new-cent-lev"+std::to_string(ilev));
        
        const MultiCutFab& bcent_new = new_factory[ilev]->getBndryCent();
        VisMF::Write(bcent_new.ToMultiFab(-1.,-1.), "new-bcent-lev"+std::to_string(ilev));

        const auto& areafrac_new = new_factory[ilev]->getAreaFrac();
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            VisMF::Write(areafrac_new[idim]->ToMultiFab(1.,0.),
                         "new-area"+std::to_string(idim)+"-lev"+std::to_string(ilev));
        }
        
        const auto& facecent_new = new_factory[ilev]->getFaceCent();
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            VisMF::Write(facecent_new[idim]->ToMultiFab(1.,0.),
                         "new-fcent"+std::to_string(idim)+"-lev"+std::to_string(ilev));
        }
        
        if (test_old_eb)
        {
            const FabArray<EBCellFlagFab>& cellflag_old = old_factory[ilev]->getMultiEBCellFlagFab();

            if (test_cellflag)
            {
                for (MFIter mfi(cellflag_new); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.validbox() & lev_domain_m1;
                    const auto& new_fab = cellflag_new[mfi];
                    const auto& old_fab = cellflag_old[mfi];

                    for (BoxIterator bi(bx); bi.ok(); ++bi) {
                        const IntVect& iv = bi();
                        if (new_fab(iv) != old_fab(iv)) {
                            amrex::AllPrint() << "cellflag diff " << iv << ": " << old_fab(iv)
                                              << ", " << new_fab(iv) << "\n";
                        }
                    }
                }
            }
            
            const MultiFab& vfrc_old = old_factory[ilev]->getVolFrac();
            VisMF::Write(vfrc_old, "old-vfrc-lev"+std::to_string(ilev));
            
            const MultiCutFab& cent_old = old_factory[ilev]->getCentroid();
            VisMF::Write(cent_old.ToMultiFab(0.,0.),  "old-cent-lev"+std::to_string(ilev));
            
            const MultiCutFab& bcent_old = old_factory[ilev]->getBndryCent();
            VisMF::Write(bcent_old.ToMultiFab(-1.,-1.),  "old-bcent-lev"+std::to_string(ilev));

            const auto& areafrac_old = old_factory[ilev]->getAreaFrac();
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                VisMF::Write(areafrac_old[idim]->ToMultiFab(1.,0.),
                             "old-area"+std::to_string(idim)+"-lev"+std::to_string(ilev));
            }
            
            const auto& facecent_old = old_factory[ilev]->getFaceCent();
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                VisMF::Write(facecent_old[idim]->ToMultiFab(1.,0.),
                             "old-fcent"+std::to_string(idim)+"-lev"+std::to_string(ilev));
            }
        }
    }
}

void
MyTest::initializeEB2 ()
{
    BL_PROFILE("initializeEB2");

    ParmParse pp("eb2");
    std::string geom_type;
    pp.get("geom_type", geom_type);
    
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
        EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level);
    }
    else if(geom_type == "rotate")
    {
	EB2::PlaneIF myplane({AMREX_D_DECL(0.5,0.5,0.5)},{AMREX_D_DECL(1., 1.,0.)}); 
	auto planerotz = EB2::rotate(myplane, atan(1.)*.25, 2); //rotate plane by pi/16 holding z axis
	auto gshop = EB2::makeShop(planerotz); 
	EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level); 
    } 
    else
    {
        EB2::Build(geom, max_coarsening_level, max_coarsening_level);
    }

    const EB2::IndexSpace& index_space = EB2::IndexSpace::top();
    const EB2::Level& eb_level = index_space.getLevel(geom);

    new_factory.emplace_back(
        new EBFArrayBoxFactory(eb_level, geom, grids, dmap, {ng_eb, ng_eb, ng_eb}, EBSupport::full));

    for (int ilev = 1; ilev <= max_coarsening_level; ++ilev)
    {
        int ratio = std::pow(2,ilev);
        Geometry cgeom(amrex::coarsen(geom.Domain(),ratio));
        const EB2::Level& crse_eb_level = index_space.getLevel(cgeom);
        new_factory.emplace_back(
            new EBFArrayBoxFactory(crse_eb_level, cgeom, amrex::coarsen(grids,ratio),
                                   dmap, {ng_eb, ng_eb, ng_eb}, EBSupport::full));
    }
}

