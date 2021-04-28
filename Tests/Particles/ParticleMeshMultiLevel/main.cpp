#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"
#include "AMReX_PlotFileUtil.H"
#include <AMReX_AmrParticles.H>

using namespace amrex;

struct TestParams {
    int nx;
    int ny;
    int nz;
    int nlevs;
    int max_grid_size;
    int nppc;
    bool verbose;
};

void testParticleMesh (TestParams& parms)
{
    Vector<IntVect> rr(parms.nlevs-1);
    for (int lev = 1; lev < parms.nlevs; lev++)
        rr[lev-1] = IntVect(AMREX_D_DECL(2,2,2));

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(parms.nx - 1, parms.ny - 1, parms.nz-1));
    const Box base_domain(domain_lo, domain_hi);

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per[i] = 1;

    Vector<Geometry> geom(parms.nlevs);
    geom[0].define(base_domain, &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < parms.nlevs; lev++) {
        geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                         &real_box, CoordSys::cartesian, is_per);
    }

    Vector<BoxArray> ba(parms.nlevs);
    Vector<DistributionMapping> dm(parms.nlevs);

    Box domain = base_domain;
    IntVect size = IntVect(AMREX_D_DECL(parms.nx, parms.ny, parms.nz));
    for (int lev = 0; lev < parms.nlevs; ++lev)
    {
        ba[lev].define(domain);
        ba[lev].maxSize(parms.max_grid_size);
        dm[lev].define(ba[lev]);
        domain.grow(-size/4);   // fine level cover the middle of the coarse domain
        domain.refine(2);
    }

    Vector<MultiFab> density(parms.nlevs);
    for (int lev = 0; lev < parms.nlevs; lev++) {
        density[lev].define(ba[lev], dm[lev], 1, 1);
        density[lev].setVal(0.0);
    }

    typedef ParticleContainer<1> MyParticleContainer;
    MyParticleContainer myPC(geom, dm, ba, rr);
    myPC.SetVerbose(false);

    int num_particles = parms.nppc * parms.nx * parms.ny * parms.nz;
    amrex::Print() << "Total number of particles    : " << num_particles << '\n' << '\n';

    bool serialize = true;
    int iseed = 451;
    Real mass = 10.0;

    MyParticleContainer::ParticleInitData pdata = {{mass}, {}, {}, {}};
    myPC.InitRandom(num_particles, iseed, pdata, serialize);

    amrex::ParticleToMesh(myPC, GetVecOfPtrs(density), 0, parms.nlevs-1,
        [=] AMREX_GPU_DEVICE (const MyParticleContainer::ParticleType& p,
                              amrex::Array4<amrex::Real> const& rho,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi) noexcept
        {
          amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
          amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
          amrex::Real lz = (p.pos(2) - plo[2]) * dxi[2] + 0.5;

          int i = static_cast<int>(amrex::Math::floor(lx));
          int j = static_cast<int>(amrex::Math::floor(ly));
          int k = static_cast<int>(amrex::Math::floor(lz));

          amrex::Real xint = lx - i;
          amrex::Real yint = ly - j;
          amrex::Real zint = lz - k;

          amrex::Real sx[] = {1.-xint, xint};
          amrex::Real sy[] = {1.-yint, yint};
          amrex::Real sz[] = {1.-zint, zint};

          for (int kk = 0; kk <= 1; ++kk) {
              for (int jj = 0; jj <= 1; ++jj) {
                  for (int ii = 0; ii <= 1; ++ii) {
                      amrex::Gpu::Atomic::AddNoRet(&rho(i+ii-1, j+jj-1, k+kk-1, 0),
                                                   sx[ii]*sy[jj]*sz[kk]*p.rdata(0));
                  }
              }
          }
      });

    Vector<std::string> varnames;
    varnames.push_back("density");

    Vector<std::string> particle_varnames;
    particle_varnames.push_back("mass");

    Vector<int> level_steps;
    level_steps.push_back(0);
    level_steps.push_back(0);

    int output_levs = parms.nlevs;

    Vector<const MultiFab*> outputMF(output_levs);
    Vector<IntVect> outputRR(output_levs);
    for (int lev = 0; lev < output_levs; ++lev) {
        outputMF[lev] = &density[lev];
        outputRR[lev] = IntVect(AMREX_D_DECL(2, 2, 2));
    }

    WriteMultiLevelPlotfile("plt00000", output_levs, outputMF,
                            varnames, geom, 0.0, level_steps, outputRR);
    myPC.Checkpoint("plt00000", "particle0", true, particle_varnames);
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
  pp.get("nppc", parms.nppc);
  pp.get("nlevs", parms.nlevs);
  if (parms.nppc < 1 && ParallelDescriptor::IOProcessor())
    amrex::Abort("Must specify at least one particle per cell");

  parms.verbose = false;
  pp.query("verbose", parms.verbose);

  if (parms.verbose && ParallelDescriptor::IOProcessor()) {
    std::cout << std::endl;
    std::cout << "Number of particles per cell : ";
    std::cout << parms.nppc  << std::endl;
    std::cout << "Size of domain               : ";
    std::cout << parms.nx << " " << parms.ny << " " << parms.nz << std::endl;
  }

  testParticleMesh(parms);

  amrex::Finalize();
}
