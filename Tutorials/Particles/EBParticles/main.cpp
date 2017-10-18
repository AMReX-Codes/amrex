#include <iostream>

#include <AMReX.H>
#include <AMReX_GeometryShop.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_SphereIF.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBTower.H>
#include <AMReX_EBFArrayBox.H>

#include "EBParticleContainer.H"

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    
    ParmParse pp;
    
    int size, max_step, max_grid_size;
    bool write_particles, do_nl;
    Real dx, dt;

    pp.get("size", size);
    pp.get("max_step", max_step);
    pp.get("max_grid_size", max_grid_size);
    pp.get("write_particles", write_particles);
    pp.get("dt", dt);
    pp.get("do_nl", do_nl);

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, size);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(size - 1, size - 1, size - 1));
    const Box domain(domain_lo, domain_hi);
    dx = 1.0;
    
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) 
        is_per[i] = 0; 
    Geometry geom(domain, &real_box, CoordSys::cartesian, is_per);

    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    
    DistributionMapping dmap(ba);

    // set up our sphere
    Real radius = 0.25*size;
    RealVect center = 0.5*size*RealVect::Unit;
    bool insideRegular = true;
    SphereIF sphere(radius, center, insideRegular);
    GeometryShop workshop(sphere);
    EBIndexSpace* ebis = AMReX_EBIS::instance();
    ebis->define(domain, RealVect::Zero, dx, workshop);

    // set up ebfactory
    int m_eb_basic_grow_cells = 5;
    int m_eb_volume_grow_cells = 4;
    int m_eb_full_grow_cells = 2;
    EBSupport m_eb_support_level = EBSupport::full;

    EBTower::Build();

    EBFArrayBoxFactory ebfactory(geom, ba, dmap,
                                 {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                                 m_eb_support_level);

    std::array<const MultiCutFab*, AMREX_SPACEDIM> areafrac;
    const MultiFab* volfrac;
    const MultiCutFab* bndrycent;

    areafrac  =  ebfactory.getAreaFrac();
    volfrac   = &ebfactory.getVolFrac();
    bndrycent = &ebfactory.getBndryCent();

    MultiFab dummy(ba, dmap, 1, 0, MFInfo(), ebfactory);

    int num_neighbor_cells = 1;
    EBParticleContainer myPC(geom, dmap, ba);

    myPC.InitParticles(radius, center);

    const int lev = 0;

    for (int i = 0; i < max_step; i++) {
        if (write_particles) myPC.writeParticles(i);
        
        myPC.moveParticles(dt);

        myPC.Redistribute();

        myPC.bounceWalls(dummy, volfrac, bndrycent, areafrac);

        myPC.Redistribute();

    }

    if (write_particles) myPC.writeParticles(max_step);

    EBTower::Destroy();
    
    amrex::Finalize();
}
