#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"
#include "AMReX_PlotFileUtil.H"
#include <AMReX_AmrParticles.H>

#include "mypc.H"
#include "trilinear_deposition_K.H"

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
    int is_per[] = {AMREX_D_DECL(1,1,1)};

    Vector<Geometry> geom(parms.nlevs);
    geom[0].define(base_domain, &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < parms.nlevs; lev++) {
        geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                         &real_box, CoordSys::cartesian, is_per);
    }

    Vector<BoxArray> ba(parms.nlevs);
    Vector<DistributionMapping> dm(parms.nlevs);

    Box domain = base_domain;
    IntVect size(AMREX_D_DECL(parms.nx, parms.ny, parms.nz));
    for (int lev = 0; lev < parms.nlevs; ++lev)
    {
        ba[lev].define(domain);
        ba[lev].maxSize(parms.max_grid_size);
        dm[lev].define(ba[lev]);
        domain.grow(-size/4);   // fine level cover the middle of the coarse domain
        domain.refine(2);
    }

    Vector<MultiFab> density1(parms.nlevs);
    Vector<MultiFab> density2(parms.nlevs);
    for (int lev = 0; lev < parms.nlevs; lev++) {
        density1[lev].define(ba[lev], dm[lev], 1, 1);
        density1[lev].setVal(0.0);
        density2[lev].define(ba[lev], dm[lev], 1, 1);
        density2[lev].setVal(0.0);
    }

    MyParticleContainer myPC(geom, dm, ba, rr);
    myPC.SetVerbose(false);

    int num_particles = parms.nppc * parms.nx * parms.ny * parms.nz;
    amrex::Print() << "Total number of particles    : " << num_particles << '\n' << '\n';

    bool serialize = true;
    int iseed = 451;
    double mass = 10.0;

    MyParticleContainer::ParticleInitData pdata = {{mass}, {}, {}, {}};
    myPC.InitRandom(num_particles, iseed, pdata, serialize);

    //
    // Here we provide an example of one way to call ParticleToMesh
    //
    amrex::ParticleToMesh(myPC, GetVecOfPtrs(density1), 0, parms.nlevs-1,
        [=] AMREX_GPU_DEVICE (const MyParticleContainer::ParticleType& p,
                              amrex::Array4<amrex::Real> const& rho,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi) noexcept
        {
            ParticleInterpolator::Linear interp(p, plo, dxi);

            interp.ParticleToMesh(p, rho, 0, 0, 1,
                [=] AMREX_GPU_DEVICE (const MyParticleContainer::ParticleType& part, int comp)
                {
                    return part.rdata(comp);  // no weighting
                });
        });

    //
    // Here we provide an example of another way to call ParticleToMesh
    //
    int start_part_comp = 0;
    int start_mesh_comp = 0;
    int        num_comp = 1;

    amrex::ParticleToMesh(myPC,GetVecOfPtrs(density2),0,parms.nlevs-1,
                          TrilinearDeposition{start_part_comp,start_mesh_comp,num_comp});

    //
    // Now write the output from each into separate plotfiles for comparison
    //

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
        outputMF[lev] = &density1[lev];
        outputRR[lev] = IntVect(AMREX_D_DECL(2, 2, 2));
    }
    WriteMultiLevelPlotfile("plt00000_v1", output_levs, outputMF,
                            varnames, geom, 0.0, level_steps, outputRR);
    myPC.WritePlotFile("plt00000_v1", "particle0", particle_varnames);

    for (int lev = 0; lev < output_levs; ++lev) {
        outputMF[lev] = &density2[lev];
        outputRR[lev] = IntVect(AMREX_D_DECL(2, 2, 2));
    }
    WriteMultiLevelPlotfile("plt00000_v2", output_levs, outputMF,
                            varnames, geom, 0.0, level_steps, outputRR);
    myPC.WritePlotFile("plt00000_v2", "particle0", particle_varnames);
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
