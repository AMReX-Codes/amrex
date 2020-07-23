#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

struct TestParams {
  int nx;
  int ny;
  int nz;
  int max_grid_size;
  int nppc;
  int nlevs;
  bool verbose;
};

void test_ghosts_and_virtuals (TestParams& parms)
{
    int nlevs = parms.nlevs;

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    RealBox fine_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
       fine_box.setLo(n,0.25);
       fine_box.setHi(n,0.75);
    }

    IntVect domain_lo(D_DECL(0 , 0, 0));
    IntVect domain_hi(D_DECL(parms.nx - 1, parms.ny - 1, parms.nz-1));
    const Box domain(domain_lo, domain_hi);

    // Define the refinement ratio
    Vector<int> rr(nlevs-1);
    for (int lev = 1; lev < nlevs; lev++)
        rr[lev-1] = 2;

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per[i] = 1;

    // This defines a Geometry object which is useful for writing the plotfiles
    Vector<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < nlevs; lev++) {
	geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
			 &real_box, CoordSys::cartesian, is_per);
    }

    Vector<BoxArray> ba(nlevs);
    ba[0].define(domain);

    if (nlevs > 1) {
        int n_fine = parms.nx*rr[0];
        IntVect refined_lo(D_DECL(n_fine/4,n_fine/4,n_fine/4));
        IntVect refined_hi(D_DECL(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1));

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        ba[1].define(refined_patch);
    }

    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(parms.max_grid_size);
    }

    Vector<DistributionMapping> dmap(nlevs);
    for (int lev = 0; lev < nlevs; lev++) {
        dmap[lev] = DistributionMapping{ba[lev]};
    }

    typedef AmrParticleContainer<0, 0, 1> MyParticleContainer;
    MyParticleContainer myPC(geom, dmap, ba, rr);
    myPC.SetVerbose(false);

    int num_particles = parms.nppc * AMREX_D_TERM(parms.nx, * parms.ny, * parms.nz);
    bool serialize = true;
    int iseed = 451;
    Real mass = 10.0;
    MyParticleContainer::ParticleInitData pdata = {mass};

    myPC.InitRandom(num_particles, iseed, pdata, serialize);

    {
        const int src_lev = 1;
        const int dst_lev = 0;
        MyParticleContainer virtPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType virts;
        myPC.CreateVirtualParticles(1, virts);
        virtPC.AddParticlesAtLevel(virts, 0);
    }

    {
        const int ngrow = 1;
        const int src_lev = 0;
        const int dst_lev = 1;
        MyParticleContainer ghostPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType ghosts;
        myPC.CreateGhostParticles(src_lev, ngrow, ghosts);
        ghostPC.AddParticlesAtLevel(ghosts, dst_lev, ngrow);
    }
}

int main(int argc, char* argv[])
{
  amrex::Initialize(argc,argv);

  ParmParse pp;

  TestParams parms;

  pp.get("nx", parms.nx);
  pp.get("ny", parms.ny);
  pp.get("nz", parms.nz);
  pp.get("max_grid_size", parms.max_grid_size);
  pp.get("nlevs", parms.nlevs);
  pp.get("nppc", parms.nppc);
  if (parms.nppc < 1 && ParallelDescriptor::IOProcessor())
    amrex::Abort("Must specify at least one particle per cell");

  parms.verbose = false;
  pp.query("verbose", parms.verbose);

  if (parms.verbose && ParallelDescriptor::IOProcessor()) {
    std::cout << std::endl;
    std::cout << "Number of particles per cell : ";
    std::cout << parms.nppc  << std::endl;
    std::cout << "Size of domain               : ";
    std::cout << "Num levels: ";
    std::cout << parms.nlevs << std::endl;
    std::cout << parms.nx << " " << parms.ny << " " << parms.nz << std::endl;
  }

  test_ghosts_and_virtuals(parms);

  amrex::Finalize();
}
