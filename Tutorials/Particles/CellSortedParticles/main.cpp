#include <iostream>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

#include "CellSortedPC.H"

using namespace amrex;

struct TestParams {
    IntVect ncell;      // num cells in domain
    IntVect nppc;       // number of particles per cell in each dim
    int max_grid_size;
    int nsteps;
};

void test_cell_sorted(const TestParams& parms)
{

    BL_PROFILE("test_em_pic");
    BL_PROFILE_VAR_NS("evolve_time", blp_evolve);
    BL_PROFILE_VAR_NS("init_time", blp_init);

    BL_PROFILE_VAR_START(blp_init);

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }
    
    IntVect domain_lo(AMREX_D_DECL(0, 0, 0)); 
    IntVect domain_hi(AMREX_D_DECL(parms.ncell[0]-1,
                                   parms.ncell[1]-1,
                                   parms.ncell[2]-1)); 
    const Box domain(domain_lo, domain_hi);
    
    int coord = 0;
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) 
        is_per[i] = 1; 
    Geometry geom(domain, &real_box, coord, is_per);
    
    BoxArray ba(domain);
    ba.maxSize(parms.max_grid_size);
    DistributionMapping dm(ba);
    
    amrex::Print() << "Initializing particles... ";

    CellSortedParticleContainer particles(geom, dm, ba);
    particles.InitParticles(parms.nppc);
    
    amrex::Print() << "Done. " << std::endl;

    BL_PROFILE_VAR_STOP(blp_init);
    
    amrex::Print() << "Starting main loop... " << std::endl;

    BL_PROFILE_VAR_START(blp_evolve);

    for (int step = 0; step < parms.nsteps; ++step)
    {
        amrex::Print() << "    Time step: " <<  step << std::endl;
        amrex::Print() << " Number of particles is " << particles.TotalNumberOfParticles() << std::endl;
        particles.MoveParticles();
        particles.Redistribute();
        particles.ReBin();
	amrex::Print() << " Number of particles in cell vectors is " << particles.SumCellVectors() << std::endl;
    }

    amrex::Print() << "Done. " << std::endl;

    BL_PROFILE_VAR_STOP(blp_evolve);
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    
    amrex::InitRandom(451);

    ParmParse pp;

    TestParams parms;

    pp.get("ncell", parms.ncell);
    pp.get("nppc",  parms.nppc);
    pp.get("max_grid_size", parms.max_grid_size);
    pp.get("nsteps", parms.nsteps);

    test_cell_sorted(parms);

    amrex::Finalize();
}
