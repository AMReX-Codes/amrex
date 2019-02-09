
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

#include "MDParticleContainer.H"

using namespace amrex;

struct TestParams
{
    int size;
    int max_grid_size;
    int nsteps;
    bool print_neighbor_list;
};

void main_main();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    main_main();
    amrex::Finalize();
}

void get_test_params(TestParams& params)
{
    ParmParse pp;
    pp.get("size", params.size);
    pp.get("max_grid_size", params.max_grid_size);
    pp.get("nsteps", params.nsteps);
    pp.get("print_neighbor_list", params.print_neighbor_list);
}

void main_main ()
{
    TestParams params;
    get_test_params(params);

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, params.size);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(params.size-1,params.size-1,params.size-1));
    const Box domain(domain_lo, domain_hi);

    int coord = 0;
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per[i] = 1;
    Geometry geom(domain, &real_box, coord, is_per);
    
    BoxArray ba(domain);
    ba.maxSize(params.max_grid_size);
    DistributionMapping dm(ba);

    MDParticleContainer pc(geom, dm, ba);

    IntVect nppc = IntVect(AMREX_D_DECL(1, 1, 1));
    constexpr Real dt = 0.1;

    pc.InitParticles(nppc, 1.0, 0.0);

    for (int step = 0; step < params.nsteps; ++step) {

        pc.BuildNeighborList();

        if (params.print_neighbor_list) pc.printNeighborList();

        pc.computeForces();
        
        pc.moveParticles(dt);
    }
}
