#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"
#include "AMReX_PlotFileUtil.H"
#include <AMReX_ParticleMesh.H>
#include <AMReX_ParticleInterpolators.H>

using namespace amrex;

struct TestParams {
  int nx;
  int ny;
  int nz;
  int max_grid_size;
  int nppc;
  bool verbose;
};

void testParticleMesh (TestParams& parms)
{

  RealBox real_box;
  for (int n = 0; n < AMREX_SPACEDIM; n++) {
    real_box.setLo(n, 0.0);
    real_box.setHi(n, 1.0);
  }

  IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
  IntVect domain_hi(AMREX_D_DECL(parms.nx - 1, parms.ny - 1, parms.nz-1));
  const Box domain(domain_lo, domain_hi);

  // This sets the boundary conditions to be doubly or triply periodic
  int is_per[] = {AMREX_D_DECL(1,1,1)};
  Geometry geom(domain, &real_box, CoordSys::cartesian, is_per);

  BoxArray ba(domain);
  ba.maxSize(parms.max_grid_size);
  if (parms.verbose && ParallelDescriptor::IOProcessor()) {
    std::cout << "Number of boxes              : " << ba[0].size() << '\n' << '\n';
  }

  DistributionMapping dmap(ba);

  MultiFab partMF(ba, dmap, 1 + AMREX_SPACEDIM, 1);
  partMF.setVal(0.0);

  iMultiFab partiMF(ba, dmap, 1 + AMREX_SPACEDIM, 1);
  partiMF.setVal(0);

  typedef ParticleContainer<1 + 2*AMREX_SPACEDIM, 1> MyParticleContainer;
  MyParticleContainer myPC(geom, dmap, ba);
  myPC.SetVerbose(false);

  int num_particles = parms.nppc * parms.nx * parms.ny * parms.nz;
  if (ParallelDescriptor::IOProcessor())
    std::cout << "Total number of particles    : " << num_particles << '\n' << '\n';

  bool serialize = true;
  int iseed = 451;
  double mass = 10.0;

  MyParticleContainer::ParticleInitData pdata = {{mass, AMREX_D_DECL(1.0, 2.0, 3.0), AMREX_D_DECL(0.0, 0.0, 0.0)}, {},{},{}};
  myPC.InitRandom(num_particles, iseed, pdata, serialize);

  int nc = 1 + AMREX_SPACEDIM;
  const auto plo = geom.ProbLoArray();
  const auto dxi = geom.InvCellSizeArray();
  amrex::ParticleToMesh(myPC, partMF, 0,
                        [=] AMREX_GPU_DEVICE (const MyParticleContainer::ParticleTileType::ConstParticleTileDataType& ptd, int i,
                                              amrex::Array4<amrex::Real> const& rho)
      {
          auto p = ptd.m_aos[i];
          ParticleInterpolator::Linear interp(p, plo, dxi);

          interp.ParticleToMesh(p, rho, 0, 0, 1,
                      [=] AMREX_GPU_DEVICE (const MyParticleContainer::ParticleType& part, int comp)
                      {
                          return part.rdata(comp);  // no weighting
                      });

          interp.ParticleToMesh(p, rho, 1, 1, AMREX_SPACEDIM,
                      [=] AMREX_GPU_DEVICE (const MyParticleContainer::ParticleType& part, int comp)
                      {
                          return part.rdata(0) * p.rdata(comp);  // mass weight these comps
                      });
      });

  MultiFab acceleration(ba, dmap, AMREX_SPACEDIM, 1);
  acceleration.setVal(5.0);

  nc = AMREX_SPACEDIM;
  amrex::MeshToParticle(myPC, acceleration, 0,
      [=] AMREX_GPU_DEVICE (MyParticleContainer::ParticleType& p,
                            amrex::Array4<const amrex::Real> const& acc)
      {
          ParticleInterpolator::Linear interp(p, plo, dxi);

          interp.MeshToParticle(p, acc, 0, 1+AMREX_SPACEDIM, nc,
                  [=] AMREX_GPU_DEVICE (amrex::Array4<const amrex::Real> const& arr,
                                        int i, int j, int k, int comp)
                  {
                      return arr(i, j, k, comp);  // no weighting
                  },
                  [=] AMREX_GPU_DEVICE (MyParticleContainer::ParticleType& part,
                                        int comp, amrex::Real val)
                  {
                      part.rdata(comp) += ParticleReal(val);
                  });
      });

  // now also try the iMultiFab versions
  amrex::ParticleToMesh(myPC, partiMF, 0,
      [=] AMREX_GPU_DEVICE (const MyParticleContainer::SuperParticleType& p,
                            amrex::Array4<int> const& count)
      {
          ParticleInterpolator::Nearest interp(p, plo, dxi);

          interp.ParticleToMesh(p, count, 0, 0, 1,
              [=] AMREX_GPU_DEVICE (const MyParticleContainer::ParticleType& /*p*/, int /*comp*/) -> int
              {
                  return 1;  // just count the particles per cell
              });
      });

  amrex::MeshToParticle(myPC, partiMF, 0,
                        [=] AMREX_GPU_DEVICE (const MyParticleContainer::ParticleTileType::ParticleTileDataType& ptd, int ip,
                                              amrex::Array4<const int> const& count)
      {
          auto& p = ptd.m_aos[ip];
          ParticleInterpolator::Nearest interp(p, plo, dxi);

          interp.MeshToParticle(p, count, 0, 0, 1,
                  [=] AMREX_GPU_DEVICE (amrex::Array4<const int> const& arr,
                                        int i, int j, int k, int comp)
                  {
                      return arr(i, j, k, comp);  // no weighting
                  },
                  [=] AMREX_GPU_DEVICE (MyParticleContainer::ParticleType& part,
                                        int comp, int val)
                  {
                      part.idata(comp) = val;
                  });
      });

  WriteSingleLevelPlotfile("plot", partMF,
                           {"density", AMREX_D_DECL("vx", "vy", "vz")},
                           geom, 0.0, 0);

  myPC.WritePlotFile("plot", "particle0");
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
