#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"
#include "AMReX_PlotFileUtil.H"
#include <AMReX_ParticleMesh.H>

using namespace amrex;

struct TestParams {
  int nx;
  int ny;
  int nz;
  int max_grid_size;
  int nppc;
  bool verbose;
};

void testParticleMesh(TestParams& parms)
{

  RealBox real_box;
  for (int n = 0; n < BL_SPACEDIM; n++) {
    real_box.setLo(n, 0.0);
    real_box.setHi(n, 1.0);
  }

  IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
  IntVect domain_hi(AMREX_D_DECL(parms.nx - 1, parms.ny - 1, parms.nz-1));
  const Box domain(domain_lo, domain_hi);

  // This sets the boundary conditions to be doubly or triply periodic
  int is_per[BL_SPACEDIM];
  for (int i = 0; i < BL_SPACEDIM; i++) 
    is_per[i] = 1; 
  Geometry geom(domain, &real_box, CoordSys::cartesian, is_per);

  BoxArray ba(domain);
  ba.maxSize(parms.max_grid_size);
  if (parms.verbose && ParallelDescriptor::IOProcessor()) {
    std::cout << "Number of boxes              : " << ba[0].size() << '\n' << '\n';
  }

  DistributionMapping dmap(ba);

  MultiFab partMF(ba, dmap, 1 + BL_SPACEDIM, 1);
  partMF.setVal(0.0);

  typedef ParticleContainer<1 + 2*BL_SPACEDIM> MyParticleContainer;
  MyParticleContainer myPC(geom, dmap, ba);
  myPC.SetVerbose(false);

  int num_particles = parms.nppc * parms.nx * parms.ny * parms.nz;
  if (ParallelDescriptor::IOProcessor())
    std::cout << "Total number of particles    : " << num_particles << '\n' << '\n';

  bool serialize = true;
  int iseed = 451;
  Real mass = 10.0;

  MyParticleContainer::ParticleInitData pdata = {mass, AMREX_D_DECL(1.0, 2.0, 3.0), AMREX_D_DECL(0.0, 0.0, 0.0)};
  myPC.InitRandom(num_particles, iseed, pdata, serialize);

  int nc = 1 + BL_SPACEDIM;
  const auto plo = geom.ProbLoArray();
  const auto dxi = geom.InvCellSizeArray();
  amrex::ParticleToMesh(myPC, partMF, 0,
      [=] AMREX_GPU_DEVICE (const MyParticleContainer::ParticleType& p,
                            amrex::Array4<amrex::Real> const& rho)
      {
          amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
          amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
          amrex::Real lz = (p.pos(2) - plo[2]) * dxi[2] + 0.5;

          int i = std::floor(lx);
          int j = std::floor(ly);
          int k = std::floor(lz);

          amrex::Real xint = lx - i;
          amrex::Real yint = ly - j;
          amrex::Real zint = lz - k;

          amrex::Real sx[] = {1.-xint, xint};
          amrex::Real sy[] = {1.-yint, yint};
          amrex::Real sz[] = {1.-zint, zint};
    
          for (int kk = 0; kk <= 1; ++kk) { 
              for (int jj = 0; jj <= 1; ++jj) { 
                  for (int ii = 0; ii <= 1; ++ii) {
                      amrex::Gpu::Atomic::Add(&rho(i+ii-1, j+jj-1, k+kk-1, 0),
                                              sx[ii]*sy[jj]*sz[kk]*p.rdata(0));
                  }
              }
          }

          for (int comp=1; comp < nc; ++comp) { 
             for (int kk = 0; kk <= 1; ++kk) { 
                  for (int jj = 0; jj <= 1; ++jj) { 
                      for (int ii = 0; ii <= 1; ++ii) {
                          amrex::Gpu::Atomic::Add(&rho(i+ii-1, j+jj-1, k+kk-1, comp),
                                                  sx[ii]*sy[jj]*sz[kk]*p.rdata(0)*p.rdata(comp));
                      }
                  }
              }
          }
      });

  MultiFab acceleration(ba, dmap, BL_SPACEDIM, 1);
  acceleration.setVal(5.0);

  nc = BL_SPACEDIM;
  amrex::MeshToParticle(myPC, acceleration, 0,
      [=] AMREX_GPU_DEVICE (MyParticleContainer::ParticleType& p,
                            amrex::Array4<const amrex::Real> const& acc)
      {
          amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
          amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
          amrex::Real lz = (p.pos(2) - plo[2]) * dxi[2] + 0.5;

          int i = std::floor(lx);
          int j = std::floor(ly);
          int k = std::floor(lz);

          amrex::Real xint = lx - i;
          amrex::Real yint = ly - j;
          amrex::Real zint = lz - k;

          amrex::Real sx[] = {1.-xint, xint};
          amrex::Real sy[] = {1.-yint, yint};
          amrex::Real sz[] = {1.-zint, zint};

          for (int comp=0; comp < nc; ++comp) {                    
              for (int kk = 0; kk <= 1; ++kk) { 
                  for (int jj = 0; jj <= 1; ++jj) { 
                      for (int ii = 0; ii <= 1; ++ii) {
                          p.rdata(4+comp) += sx[ii]*sy[jj]*sz[kk]*acc(i+ii-1,j+jj-1,k+kk-1,comp);
                      }
                  }
              }
          }
      });
  
  WriteSingleLevelPlotfile("plot", partMF, 
                           {"density", "vx", "vy", "vz"},
                           geom, 0.0, 0);

  myPC.Checkpoint("plot", "particle0");
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
