
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

#include "CheckPair.H"

#include "MDParticleContainer.H"

using namespace amrex;

struct TestParams
{
    IntVect size;
    int max_grid_size;
    int nsteps;
    bool print_neighbor_list;
    bool write_particles;
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
    pp.get("write_particles", params.write_particles);
}

void main_main ()
{
    TestParams params;
    get_test_params(params);

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, params.size[n]);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(params.size[0]-1,params.size[1]-1,params.size[2]-1));
    const Box domain(domain_lo, domain_hi);

    int coord = 0;
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per[i] = 0;
    Geometry geom(domain, &real_box, coord, is_per);
    
    BoxArray ba(domain);
    ba.maxSize(params.max_grid_size);
    DistributionMapping dm(ba);

    const int ncells = 1;
    MDParticleContainer pc(geom, dm, ba, ncells);

    IntVect nppc = IntVect(AMREX_D_DECL(1, 1, 1));
    constexpr Real dt = 0.0005;

    pc.InitParticles(nppc, 1.0, 0.0);

    for (int step = 0; step < params.nsteps; ++step) {

        amrex::Print() << "Taking step " << step << "\n";

        pc.fillNeighbors();

        pc.buildNeighborList(CheckPair());

        if (params.print_neighbor_list) pc.printNeighborList();

        pc.computeForces();

        pc.moveParticles(dt);

        pc.RedistributeLocal();
    }

    if (params.write_particles) pc.writeParticles(params.nsteps);
}
