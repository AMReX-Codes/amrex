
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

#include "MDParticleContainer.H"

using namespace amrex;

struct TestParams
{
    IntVect ncell;      // num cells in domain
    IntVect nppc;       // number of particles per cell in each dim
    int max_grid_size;
    int nsteps;
    bool write_plot;
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
    pp.get("ncell", params.ncell);
    pp.get("nppc",  params.nppc);
    pp.get("max_grid_size", params.max_grid_size);
    pp.get("nsteps", params.nsteps);
    pp.get("write_plot", params.write_plot);
}

void main_main ()
{
    TestParams params;
    get_test_params(params);

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(params.ncell[0]-1,params.ncell[1]-1,params.ncell[2]-1));
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

    pc.InitParticles(params.nppc, 0.01, 10.0);

    pc.BuildNeighborList();
}
