#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Particles.H>
#include <AMReX_PlotFileUtil.H>

#include "MyParticleContainer.H"

using namespace amrex;

struct TestParams {
  int nx;
  int ny;
  int nz;
  int max_grid_size;
  int nppc;
  bool verbose;
};

// declare routines below
void solve_for_accel(const Array<MultiFab*>& rhs,
		     const Array<MultiFab*>& phi,
		     const Array<MultiFab*>& grad_phi,
                     const Array<Geometry>& geom,
		     int base_level, int finest_level, Real offset);

void field_solve(MultiFab& density, MultiFab& phi, MultiFab& E, const Geometry& geom) {

  BL_PROFILE("Field Solve.")
  
  MultiFab tmp(density.boxArray(), density.DistributionMap(), 1, 0);
  MultiFab::Copy(tmp, density, 0, 0, 1, 0);

  Array<MultiFab*> rhs_in(1);
  Array<MultiFab*> phi_in(1);
  Array<MultiFab*> gradphi_in(1);
  Array<Geometry> geom_in(1);

  Real offset = 209715200.0;

  rhs_in[0]     = &tmp;
  phi_in[0]     = &phi;
  gradphi_in[0] = &E;
  geom_in[0]    = geom;

  solve_for_accel(rhs_in, phi_in, gradphi_in, geom_in, 0, 0, offset);

}

void test_assign_density(TestParams& parms)
{

  RealBox real_box;
  for (int n = 0; n < BL_SPACEDIM; n++) {
    real_box.setLo(n, 0.0);
    real_box.setHi(n, 1.0);
  }

  IntVect domain_lo(0 , 0, 0); 
  IntVect domain_hi(parms.nx - 1, parms.ny - 1, parms.nz-1); 
  const Box domain(domain_lo, domain_hi);

  int coord = 0;
  int is_per[BL_SPACEDIM];
  for (int i = 0; i < BL_SPACEDIM; i++) 
    is_per[i] = 1; 
  Geometry geom(domain, &real_box, coord, is_per);

  BoxArray ba(domain);
  ba.maxSize(parms.max_grid_size);
  if (parms.verbose && ParallelDescriptor::IOProcessor()) {
    std::cout << "Number of boxes              : " << ba.size() << '\n' << '\n';
  }

  DistributionMapping dmap(ba);

  MultiFab density(ba, dmap, 1, 1);
  MultiFab E(ba, dmap, 3, 1);
  MultiFab phi(ba, dmap, 1, 1);

  density.setVal(0.0);
  phi.setVal(0.0);
  E.setVal(0.0);

  MyParticleContainer myPC(geom, dmap, ba);

  int num_particles = parms.nppc * parms.nx * parms.ny * parms.nz;
  if (ParallelDescriptor::IOProcessor())
    std::cout << "Total number of particles    : " << num_particles << '\n' << '\n';

  Real mass = 10.0;
  myPC.InitParticles(num_particles, mass);

  myPC.Deposit(density, E);

  std::cout << "Total particle mass is: " << myPC.sumParticleMass(0, 0) << std::endl;
  std::cout << "Total mesh mass is: " << density.sum(0) << std::endl;
  
  field_solve(density, phi, E, geom);

  WriteSingleLevelPlotfile("plt00000", density, {"density"}, geom, 0.0, 0);
  myPC.Checkpoint("plt00000", "particle0", true);
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
  
  test_assign_density(parms);
  
  amrex::Finalize();
}
