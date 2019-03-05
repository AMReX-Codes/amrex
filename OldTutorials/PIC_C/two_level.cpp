#include <iostream>
#include <memory>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_BLFort.H>
#include <AMReX_MacBndry.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>

using namespace amrex;

// declare routines below
void solve_for_accel(const Vector<MultiFab*>& rhs,
		     const Vector<MultiFab*>& phi,
		     const Vector<MultiFab*>& grad_phi, 
                     const Vector<Geometry>& geom,
		     int base_level, int finest_level, Real offset);

void 
two_level(int nlevs, int nx, int ny, int nz, int max_grid_size, int nppc, bool verbose)
{
    BL_PROFILE("two_level");
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

    // This defines the physical size of the fine box which we will use to constrain where we put the particles.
    RealBox fine_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
       fine_box.setLo(n,0.4);
       fine_box.setHi(n,0.6);
    }

    // This defines the physical size of the fine box which we will use to constrain where we put the particles.
    RealBox left_corner;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
       left_corner.setLo(n,0.0);
       left_corner.setHi(n,1.0);
    }
    left_corner.setLo(0,0.7);
    left_corner.setHi(0,0.9);
    left_corner.setLo(1,0.7);
    left_corner.setHi(1,0.9);
    left_corner.setLo(2,0.7);
    left_corner.setHi(2,0.9);

    // Define the lower and upper corner of a 3D domain
    IntVect domain_lo(0 , 0, 0); 
    int n_cell = 64;
    IntVect domain_hi(n_cell-1,n_cell-1,n_cell-1); 
 
    // Build a box for the level 0 domain
    const Box domain(domain_lo, domain_hi);

    // Define the refinement ratio
    Vector<int> rr(nlevs-1);
    for (int lev = 1; lev < nlevs; lev++)
        rr[lev-1] = 2;

    // This says we are using Cartesian coordinates
    int coord = 0;

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) is_per[i] = 1; 

    // This defines a Geometry object which is useful for writing the plotfiles  
    Vector<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, coord, is_per);
    for (int lev = 1; lev < nlevs; lev++)
    {
	geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
			 &real_box, coord, is_per);
    }

    // ********************************************************************************************
    // This now defines the refinement information 
    // ********************************************************************************************

    // Build an array of BoxArrays,
    // then initialize the level 0 BoxArray with the domain.
    Vector<BoxArray> ba(nlevs);
    ba[0].define(domain);

    // Now we make the refined level be the center eighth of the domain
    if (nlevs > 1) 
    {
        int n_fine = n_cell*rr[0];
        IntVect refined_lo(n_fine/4,n_fine/4,n_fine/4); 
        IntVect refined_hi(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1);

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        ba[1].define(refined_patch);
    }

    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(max_grid_size);
    }

    // Break the BoxArrays at both levels into max_grid_size^3 boxes
    ba[0].maxSize(max_grid_size);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
       std::cout << "Number of boxes at level 0   : " << ba[0].size() << '\n' << '\n';
    }

    // ********************************************************************************************
    // Set up the arrays for the solve
    // ********************************************************************************************

    // build a multifab for the rhs on the box array with 
    Vector<std::unique_ptr<MultiFab> > rhs(nlevs); 
    Vector<std::unique_ptr<MultiFab> > phi(nlevs);
    Vector<std::unique_ptr<MultiFab> > grad_phi(nlevs);
    Vector<DistributionMapping> dmap(nlevs);

    for (int lev = 0; lev < nlevs; lev++)
    {
	dmap[lev] = DistributionMapping{ba[lev]};
	//                                                 # component # ghost cells
	rhs     [lev].reset(new MultiFab(ba[lev],dmap[lev],1          ,0));
	phi     [lev].reset(new MultiFab(ba[lev],dmap[lev],1          ,1));
	grad_phi[lev].reset(new MultiFab(ba[lev],dmap[lev],BL_SPACEDIM,1));

	rhs     [lev]->setVal(0.0);
	phi     [lev]->setVal(0.0);
	grad_phi[lev]->setVal(0.0);
    }

    // Define a new particle container to hold my particles.
    // This holds a charge as well as three velocity components, three acceleration components  and three position components.
    typedef ParticleContainer<1+2*BL_SPACEDIM> MyParticleContainer;
    
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
        std::cout << "Total number of particles = " << num_particles << std::endl;

    // Here we do a set of two experiments.  
    // 1) Do a single-level solve on level 0 as if there is no level 1, then 
    //    do a single-level solve on level 1 using the boundary conditions from level 0 from the previous step.
    // 2) Do a multi-level solve on levels 0 and 1.
    // We assume for all these tests that the particles within the area covered by the refined patch.

    // **************************************************************************
    // 1) Do a single-level solve on level 0, then do a solve on level 1 using the solution
    //    from level 0 for boundary condtions
    // **************************************************************************

    // Initialize "num_particles" number of particles, each with mass/charge "mass"
    bool serialize = false;
    MyParticleContainer::ParticleInitData pdata = {mass, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    MyPC->InitRandom(num_particles,iseed,pdata,serialize,fine_box);

    // Write out the positions, masses and accelerations of each particle.
    MyPC->WriteAsciiFile("Particles_before");

    // **************************************************************************
    // Compute the total charge of all particles in order to compute the offset
    //     to make the Poisson equations solvable
    // **************************************************************************

    Real offset = 0.;
    if (geom[0].isAllPeriodic()) 
    {
        for (int lev = 0; lev < nlevs; lev++)
            offset += MyPC->sumParticleMass(0,lev);
        offset /= geom[0].ProbSize();
    }

    // **************************************************************************

    // Define the density on level 0 from all particles at all levels
    int base_level   = 0;
    int finest_level = nlevs-1;

    Vector<std::unique_ptr<MultiFab> > PartMF;
    MyPC->AssignDensity(0, PartMF, base_level, 1, finest_level);
 
    // **************************************************************************
    // Define this to be solve at level 0 only
    // **************************************************************************
    base_level   = 0;
    finest_level = 0;

    // Use multigrid to solve Lap(phi) = rhs with periodic boundary conditions (set above)
    if (ParallelDescriptor::IOProcessor())
        std::cout << "Solving for phi at level 0 ... " << std::endl;
    solve_for_accel(amrex::GetVecOfPtrs(rhs),
		    amrex::GetVecOfPtrs(phi),
		    amrex::GetVecOfPtrs(grad_phi),
		    geom,base_level,finest_level,offset);
    if (ParallelDescriptor::IOProcessor())
        std::cout << "Solved  for phi at level 0 ... " << std::endl;
    
    // Fill the particle data with the acceleration at the particle location
    // Note that we are calling moveKick with accel_comp > BL_SPACEDIM
    //      which means with dt = 0 we don't move the particle or set a velocity
    MyPC->moveKick(*grad_phi[0],nlevs-1,dummy_dt,1.0,1.0,accel_comp);
    
    // Write out the positions, masses and accelerations of each particle.
    MyPC->WriteAsciiFile("Particles_after_level0_solve");
    
    if (nlevs > 1)
    {

    // **************************************************************************
    // Define this to be solve at level 1 only
    // **************************************************************************

    base_level   = 1;
    finest_level = 1;

    // Use multigrid to solve Lap(phi) = rhs with boundary conditions from level 0
    if (ParallelDescriptor::IOProcessor())
       std::cout << "Solving for phi at level 1 ... " << std::endl;
    solve_for_accel(amrex::GetVecOfPtrs(rhs),
		    amrex::GetVecOfPtrs(phi),
		    amrex::GetVecOfPtrs(grad_phi),
		    geom,base_level,finest_level,offset);
    if (ParallelDescriptor::IOProcessor())
       std::cout << "Solved  for phi at level 1 ... " << std::endl;

    // Fill the particle data with the acceleration at the particle location
    // Note that we are calling moveKick with accel_comp > BL_SPACEDIM
    //      which means with dt = 0 we don't move the particle or set a velocity
    for (int lev = 0; lev < nlevs; lev++)
        MyPC->moveKick(*grad_phi[lev],lev,dummy_dt,1.0,1.0,accel_comp);

    // Write out the positions, masses and accelerations of each particle.
    MyPC->WriteAsciiFile("Particles_after_level1_solve");

    // **************************************************************************
    // 2) Do a multi-level solve on levels 0 and 1.
    // **************************************************************************

    // Reset everything to 0
    for (int lev = 0; lev < nlevs; lev++)
    {
	phi[lev]->setVal(0.0);
	grad_phi[lev]->setVal(0.0);
    }

    // Define this to be solve at multi-level solve
    base_level   = 0;
    finest_level = 1;

    // Redistribute the particles since we include both levels now.
    MyPC->Redistribute();

    // Use the PIC approach to deposit the "mass" onto the grid
    MyPC->AssignDensity(0, rhs, base_level,1,finest_level);

    // Use multigrid to solve Lap(phi) = rhs with periodic boundary conditions (set above)
    if (ParallelDescriptor::IOProcessor())
       std::cout << "Solving for phi at levels 0 and 1 ... " << std::endl;
    solve_for_accel(amrex::GetVecOfPtrs(rhs),
		    amrex::GetVecOfPtrs(phi),
		    amrex::GetVecOfPtrs(grad_phi),
		    geom,base_level,finest_level,offset);
    if (ParallelDescriptor::IOProcessor())
       std::cout << "Solved  for phi at levels 0 and 1 ... " << std::endl;

    // Fill the particle data with the acceleration at the particle location
    // Note that we are calling moveKick with accel_comp > BL_SPACEDIM
    //      which means with dt = 0 we don't move the particle or set a velocity
    for (int lev = 0; lev < nlevs; lev++)
        MyPC->moveKick(*grad_phi[lev],lev,dummy_dt,1.0,1.0,accel_comp);

    // Write out the positions, masses and accelerations of each particle.
    MyPC->WriteAsciiFile("Particles_after_multilevel_solve");

    } // end if (nlevs > 1)

}
