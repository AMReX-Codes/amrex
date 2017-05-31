#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Particles.H>
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

void test_assign_density(TestParams& parms)
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
    
    IntVect domain_lo(0 , 0, 0); 
    IntVect domain_hi(parms.nx - 1, parms.ny - 1, parms.nz-1); 
    const Box domain(domain_lo, domain_hi);

    // Define the refinement ratio
    Array<int> rr(nlevs-1);
    for (int lev = 1; lev < nlevs; lev++)
        rr[lev-1] = 2;

    // This says we are using Cartesian coordinates
    int coord = 0;

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per[i] = 1;

    // This defines a Geometry object which is useful for writing the plotfiles  
    Array<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, coord, is_per);
    for (int lev = 1; lev < nlevs; lev++) {
	geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
			 &real_box, coord, is_per);
    }

    Array<BoxArray> ba(nlevs);
    ba[0].define(domain);
    
    // Now we make the refined level be the center eighth of the domain
    if (nlevs > 1) {
        int n_fine = parms.nx*rr[0];
        IntVect refined_lo(n_fine/4,n_fine/4,n_fine/4); 
        IntVect refined_hi(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1);

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        ba[1].define(refined_patch);
    }
    
    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(parms.max_grid_size);
    }

    Array<DistributionMapping> dmap(nlevs);

    Array<std::unique_ptr<MultiFab> > partMF(nlevs);
    Array<std::unique_ptr<MultiFab> > density(nlevs);
    Array<std::unique_ptr<MultiFab> > acceleration(nlevs);
    for (int lev = 0; lev < nlevs; lev++) {
        dmap[lev] = DistributionMapping{ba[lev]};
        partMF[lev].reset(new MultiFab(ba[lev], dmap[lev], 1, 2));
        partMF[lev]->setVal(0.0, 2);
        density[lev].reset(new MultiFab(ba[lev], dmap[lev], 1, 0));
        density[lev]->setVal(0.0);
        acceleration[lev].reset(new MultiFab(ba[lev], dmap[lev], 3, 1));
        acceleration[lev]->setVal(5.0, 1);
    }

    typedef ParticleContainer<1> MyParticleContainer;
    MyParticleContainer myPC(geom, dmap, ba, rr);
    myPC.SetVerbose(false);

    int num_particles = parms.nppc * parms.nx * parms.ny * parms.nz;
    bool serialize = true;
    int iseed = 451;
    Real mass = 10.0;
    MyParticleContainer::ParticleInitData pdata = {mass};

    //    myPC.InitRandom(num_particles, iseed, pdata, serialize, fine_box);
    myPC.InitRandom(num_particles, iseed, pdata, serialize);

    myPC.AssignDensity(0, true, partMF, 0, 1, 1);
    //myPC.AssignDensityFort(0, partMF, 0, 1, nlevs-1);

    myPC.InterpolateFort(acceleration, 0, nlevs-1);

    for (int lev = 0; lev < nlevs; ++lev) {
        MultiFab::Copy(*density[lev], *partMF[lev], 0, 0, 1, 0);
    }

    Array<std::string> varnames;
    varnames.push_back("density");

    Array<std::string> particle_varnames;
    particle_varnames.push_back("mass");

    Array<int> level_steps;
    level_steps.push_back(0);
    level_steps.push_back(0);

    int output_levs = nlevs;

    Array<const MultiFab*> outputMF(output_levs);
    Array<IntVect> outputRR(output_levs);
    for (int lev = 0; lev < output_levs; ++lev) {
        outputMF[lev] = density[lev].get();
        outputRR[lev] = IntVect(2, 2, 2);
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
  
  test_assign_density(parms);
  
  amrex::Finalize();
}
