#include <iostream>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

#include "TestParticleContainer.H"

#include "test_F.H"

using namespace amrex;

struct TestParams
{
    IntVect ncell;      // num cells in domain
    IntVect nppc;       // number of particles per cell in each dim
    int max_grid_size;
    int nsteps;
};

Real compute_dt(const Geometry& geom)
{
    const static Real cfl = 1.0;
    const Real* dx = geom.CellSize();
    const Real dt  = cfl * 1./( std::sqrt(D_TERM(  1./(dx[0]*dx[0]),
                                                 + 1./(dx[1]*dx[1]),
                                                 + 1./(dx[2]*dx[2]))) * 3.0e8 );
    return dt;
}

void run_test(const TestParams& parms)
{
    BL_PROFILE("run_test");
    BL_PROFILE_VAR_NS("evolve_time", blp_evolve);
    BL_PROFILE_VAR_NS("init_time", blp_init);

    BL_PROFILE_VAR_START(blp_init);

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
        real_box.setLo(n, -20e-6);
        real_box.setHi(n,  20e-6);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(parms.ncell[0]-1,parms.ncell[1]-1,parms.ncell[2]-1));
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

    TestParticleContainer mypc(geom, dm, ba);
    RealBox bounds = RealBox(AMREX_D_DECL(-20e-6, -20e-6, -20e-6),
                             AMREX_D_DECL( 20e-6,  20e-6,  20e-6));    
    mypc.InitParticles(parms.nppc, 0.01, 10.0, 1e12, bounds);

    amrex::Print() << "Done. " << std::endl;

    amrex::Print() << "Starting main PIC loop... " << std::endl;

    int nsteps = parms.nsteps;
    const Real dt = compute_dt(geom);

    BL_PROFILE_VAR_STOP(blp_init);

    BL_PROFILE_VAR_START(blp_evolve);

    Real time = 0.0;
    for (int step = 0; step < nsteps; ++step)
    {
        amrex::Print() << "    Time step: " <<  step << std::endl;

        mypc.MoveParticles();
        mypc.Redistribute(0, 0, 0, 1);
        AMREX_ALWAYS_ASSERT(mypc.OK());

        time += dt;
    }

    amrex::Print() << "Done. " << std::endl;
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

    run_test(parms);

    amrex::Finalize();
}
