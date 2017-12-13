#include <iostream>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_RealBox.H>
#include <AMReX_Particles.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_LoadBalanceKD.H>

using namespace amrex;

typedef ParticleContainer<0> MyParticleContainer;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    int num_cells, max_grid_size, num_procs;
   
    ParmParse pp;    
    pp.get("num_cells", num_cells);
    pp.get("num_procs", num_procs);
    pp.get("max_grid_size", max_grid_size);

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(num_cells - 1, num_cells - 1, num_cells - 1));
    const Box domain(domain_lo, domain_hi);
    
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) 
        is_per[i] = 0; 
    Geometry geom(domain, &real_box, CoordSys::cartesian, is_per);
    
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    DistributionMapping dmap(ba);

    MyParticleContainer myPC(geom, dmap, ba);    
    MyParticleContainer::ParticleInitData pdata = {};
    myPC.InitFromBinaryFile("binary_particle_file.dat", 0);

    BoxArray new_ba;
    Vector<Real> costs;
    loadBalanceKD::balance<MyParticleContainer>(myPC, new_ba, num_procs, 0.0, costs);

    std::cout << new_ba << std::endl;    
    for (int i = 0; i < new_ba.size(); ++i) {
        std::cout << costs[i] << std::endl;
    }

    Vector<int> new_pmap;
    for (int i = 0; i < new_ba.size(); ++i) {
        new_pmap.push_back(0);
    }

    DistributionMapping new_dm(new_pmap);

    myPC.SetParticleBoxArray(0, new_ba);
    myPC.SetParticleDistributionMap(0, new_dm);
    
    myPC.Redistribute();
    MultiFab new_local_cost;
    MultiFab new_global_cost;    
    loadBalanceKD::computeCost<MyParticleContainer>(myPC, new_local_cost, new_global_cost, domain, 0.0);

    WriteSingleLevelPlotfile("plt00000", new_local_cost, {"cost"},
                             geom, 0.0, 0);
    myPC.Checkpoint("plt00000", "particle0", true);
    
    amrex::Finalize();
}
