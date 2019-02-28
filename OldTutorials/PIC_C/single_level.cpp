#include <iostream>
#include <memory>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_BLFort.H>
#include <AMReX_MacBndry.H>
#include <AMReX_MultiFabUtil.H>

#include "AMReX_Particles.H"

using namespace amrex;

// declare routines below
void solve_for_accel(const Vector<MultiFab*>& rhs,
		     const Vector<MultiFab*>& phi,
		     const Vector<MultiFab*>& grad_phi,
                     const Vector<Geometry>& geom,
		     int base_level, int finest_level, Real offset);

int single_level(int nlevs, int nx, int ny, int nz, int max_grid_size, int nppc, bool verbose) 
{
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    Real strt_init, strt_assd, strt_solve, strt_mK;
    Real  end_init,  end_assd,  end_solve,  end_mK;

    // ********************************************************************************************
    // All of this defines the level 0 information -- size of box, type of boundary condition, etc.
    // ********************************************************************************************
    // This defines the physical size of the box.  Right now the box is [0,1] in each direction.
    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
       real_box.setLo(n,0.0);
       real_box.setHi(n,1.0);
    }

    // Define the lower and upper corner of a 3D domain
    IntVect domain_lo(0 , 0, 0); 
    IntVect domain_hi(nx-1,ny-1,nz-1); 
 
    // Build a box for the level 0 domain
    const Box domain(domain_lo, domain_hi);

    // This says we are using Cartesian coordinates
    int coord = 0;

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) is_per[i] = 1; 

    // This defines a Geometry object which is useful for writing the plotfiles  
    Vector<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, coord, is_per);

    // Build a BoxArray then initialize with the domain.
    Vector<BoxArray> ba(1);
    ba[0].define(domain);

    // Break the BoxArrays at both levels into max_grid_size^3 boxes
    ba[0].maxSize(max_grid_size);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
       std::cout << "Number of boxes              : " << ba[0].size() << '\n' << '\n';
    }

    // ********************************************************************************************
    // Set up the arrays for the solve
    // ********************************************************************************************

    // build a multifab for the rhs on the box array with 
    Vector<std::unique_ptr<MultiFab> > rhs(nlevs);
    Vector<std::unique_ptr<MultiFab> > phi(nlevs);
    Vector<std::unique_ptr<MultiFab> > grad_phi(nlevs);
    Vector<DistributionMapping> dmap(nlevs);

    int lev = 0;
    dmap[lev] = DistributionMapping{ba[lev]};
    rhs     [lev].reset(new MultiFab(ba[lev],dmap[lev],1          ,0));
    phi     [lev].reset(new MultiFab(ba[lev],dmap[lev],1          ,1));
    grad_phi[lev].reset(new MultiFab(ba[lev],dmap[lev],BL_SPACEDIM,1));

    rhs     [lev]->setVal(0.0);
    phi     [lev]->setVal(0.0);
    grad_phi[lev]->setVal(0.0);

    // Define a new particle container to hold my particles.
    // This holds a charge as well as three velocity components, three acceleration components  and three position components.
    typedef ParticleContainer<1+2*BL_SPACEDIM> MyParticleContainer;

    // We define the refinement ratio even though we are single level because
    //    we want to use the multilevel interface in the different calls.
    Vector<int> rr(nlevs-1);
    
    // Build a new particle container to hold my particles.
    std::unique_ptr<MyParticleContainer> MyPC(new MyParticleContainer(geom,dmap,ba,rr));

    MyPC->SetVerbose(0);

    // This allows us to write the gravitational acceleration into these components 
    int accel_comp = BL_SPACEDIM+1;
    Real dummy_dt  = 0.0;

    int num_particles = nppc * nx * ny * nz;
    int iseed = 10;
    Real mass  = 10.0;

    if (ParallelDescriptor::IOProcessor())
       std::cout << "Total number of particles    : " << num_particles << '\n' << '\n';

    // **************************************************************************
    // Do a single-level solve on level 0
    // **************************************************************************

    // Initialize "num_particles" number of particles, each with mass/charge "mass"
    bool serialize = false;

    strt_init = ParallelDescriptor::second();
    MyParticleContainer::ParticleInitData pdata = {mass, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    MyPC->InitRandom(num_particles,iseed,pdata,serialize);
    end_init = ParallelDescriptor::second() - strt_init;

    // Write out the positions, masses and accelerations of each particle.
    // if (verbose) MyPC->WriteAsciiFile("Particles_before");

    // **************************************************************************
    // Compute the total charge of all particles in order to compute the offset
    //     to make the Poisson equations solvable
    // **************************************************************************
    Real offset = 0.;
    if (geom[0].isAllPeriodic()) 
    {
        offset = MyPC->sumParticleMass(0,lev);
        offset /= geom[0].ProbSize();
    }

    // **************************************************************************

    // Define the density on level 0 from all particles at all levels
    int base_level   = 0;
    int finest_level = 0;

    Vector<std::unique_ptr<MultiFab> > PartMF(1);
    PartMF[0].reset(new MultiFab(ba[0],dmap[0],1,1));
    PartMF[0]->setVal(0.0);

    strt_assd = ParallelDescriptor::second();

    MyPC->AssignCellDensitySingleLevel(0, *PartMF[0], 0, 1, 0);

    end_assd = ParallelDescriptor::second() - strt_assd;

    MultiFab::Add(*rhs[0], *PartMF[0], 0, 0, 1, 0);
 
    // **************************************************************************
    // Define this to be solve at level 0 only
    // **************************************************************************

    base_level   = 0;
    finest_level = 0;

    strt_solve = ParallelDescriptor::second();

    // Use multigrid to solve Lap(phi) = rhs with periodic boundary conditions (set above)
    solve_for_accel(amrex::GetVecOfPtrs(rhs),
		    amrex::GetVecOfPtrs(phi),
		    amrex::GetVecOfPtrs(grad_phi),
		    geom,base_level,finest_level,offset);

    end_solve = ParallelDescriptor::second() - strt_solve;

    // Fill the particle data with the acceleration at the particle location
    // Note that we are calling moveKick with accel_comp > BL_SPACEDIM
    //      which means with dt = 0 we don't move the particle or set a velocity

    strt_mK = ParallelDescriptor::second();

    MyPC->moveKick(*grad_phi[0],nlevs-1,dummy_dt,1.0,1.0,accel_comp);

    end_mK = ParallelDescriptor::second() - strt_mK;

    // Write out the positions, masses and accelerations of each particle.
    // if (verbose) MyPC->WriteAsciiFile("Particles_after_level0_solve");

    ParallelDescriptor::ReduceRealMax(end_init ,IOProc);
    ParallelDescriptor::ReduceRealMax(end_assd,IOProc);
    ParallelDescriptor::ReduceRealMax(end_solve,IOProc);
    ParallelDescriptor::ReduceRealMax(end_mK   ,IOProc);
    if (verbose && ParallelDescriptor::IOProcessor())
    {
           std::cout << "Time in InitRandom   : " << end_init  << '\n';
           std::cout << "Time in AssignDensity: " << end_assd  << '\n';
           std::cout << "Time in Solve        : " << end_solve << '\n';
           std::cout << "Time in moveKick     : " << end_mK    << '\n';
    }
}
