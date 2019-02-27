#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_BLFort.H>
#include <AMReX_MacBndry.H>
#include <AMReX_MultiFabUtil.H>

#include "AMReX_AmrParticles.H"
#include "AMReX_Particles.H"

using namespace amrex;

// declare routines below
void solve_for_accel(const Vector<MultiFab*>& rhs,
		     const Vector<MultiFab*>& phi,
		     const Vector<MultiFab*>& grad_phi, 
                     const Vector<Geometry>& geom, int base_level, int finest_level, Real offset);

Real                getEfficiency  (const DistributionMapping& dm, const Vector<long>& cost);
DistributionMapping getCostCountDM (const Vector<long>& cost, const BoxArray& ba);
void                splitBoxes     (BoxArray& ba, Vector<long>& newcost, const Vector<long>& cost_in, int max_grid_size);

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
        
    int max_grid_size = 32;
    int min_grid_size = 4;

    // Define this as a single-level problem.
    int nlevs = 1;

    // Define this as a 2-level problem.
    // int nlevs = 2;

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
    left_corner.setLo(0,0.55);
    left_corner.setLo(1,0.55);
    left_corner.setLo(2,0.85);
    left_corner.setHi(0,0.95);
    left_corner.setHi(1,0.95);
    left_corner.setHi(2,0.95);

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

    // break the BoxArrays at both levels into 32^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(32);
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

	rhs[lev]->setVal(0.0);
	phi[lev]->setVal(0.0);
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

    int num_particles = 1000;
    int iseed = 10;
    Real mass  = 10.0;

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
    bool serialize = true;
    MyParticleContainer::ParticleInitData pdata = {mass, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    MyPC->InitRandom(num_particles,iseed,pdata,serialize,left_corner);

    // Write out the positions, masses and accelerations of each particle.
    MyPC->WriteAsciiFile("Particles_before");

    // **************************************************************************
    // Load Balance
    // **************************************************************************
    {
	const Real eff_target = 0.8;
	const int lev = 0;

	Vector<long> new_particle_cost = MyPC->NumberOfParticlesInGrid(lev);
	Real neweff = getEfficiency(dmap[0], new_particle_cost);

	if (neweff < eff_target) 
	{
	    Real oldeff;
	    Vector<long> old_particle_cost;
	    int heavy_grid_size = max_grid_size;  // for the most heavy boxes

	    do {
		oldeff = neweff;
		old_particle_cost = new_particle_cost;

		if (ParallelDescriptor::IOProcessor()) 
		{
		    std::cout << "*** " << std::endl;
		    std::cout << "*** Before remapping, # of boxes: " << new_particle_cost.size()
			      << ", efficiency: " << neweff << "\n";
		}

		BoxArray new_ba = MyPC->ParticleBoxArray(lev);
		// This returns new_particle_cost as an *estimate* of the new cost per grid, based just
		//      on dividing the cost proportionally as the grid is divided
		splitBoxes(new_ba, new_particle_cost, old_particle_cost, heavy_grid_size);
		heavy_grid_size /= 2;

		// We use this approximate cost to get a new DistrbutionMapping so we can go ahead 
		//      and move the particles
		DistributionMapping new_dm = getCostCountDM(new_particle_cost, new_ba);

		// We get an *estimate* of the new efficiency
		neweff = getEfficiency(new_dm, new_particle_cost);

		if (ParallelDescriptor::IOProcessor()) 
		{
		    std::cout << "*** If     remapping, # of boxes: " << new_particle_cost.size()
			      << ", approx. eff: " <<  neweff << "\n";
		}

		// Only if the new_ba and new_dm are expected to improve the efficiency, ...
		if (neweff > oldeff)  
		{
		    // Now we actually move the particles onto the new_ba with the new_dm
		    MyPC->SetParticleBoxArray(lev,new_ba);
		    MyPC->SetParticleDistributionMap(lev,new_dm);
		    MyPC->Redistribute();
		    
		    // This counts how many particles are *actually* in each grid of the 
		    //      ParticleContainer's new ParticleBoxArray
		    new_particle_cost = MyPC->NumberOfParticlesInGrid(lev);
		    
		    // Here we get the *actual* new efficiency
		    neweff = getEfficiency(new_dm, new_particle_cost);
		    
		    if (ParallelDescriptor::IOProcessor()) 
		    {
			std::cout << "*** After  remapping, # of boxes: " << new_particle_cost.size()
				  << ", actual  eff: " <<  neweff << "\n";
		    }
		}
	    } 
	    while (neweff < eff_target && neweff > oldeff && heavy_grid_size >= 2*min_grid_size);
	} 
	else {
	    if (ParallelDescriptor::IOProcessor()) 
	    {
		std::cout << "*** " << std::endl;
		std::cout << "*** No remapping required: # of boxes: " << ba[0].size() 
			  << ", efficiency: " <<  neweff << "\n";
		std::cout << "*** " << std::endl;
	    }
	} 
    }

    // **************************************************************************
    // Compute the total charge of all particles in order to compute the offset
    //     to make the Poisson equations solvable
    // **************************************************************************

    Real offset = 0.;
    if (geom[0].isAllPeriodic()) 
    {
        for (int lev = 0; lev < nlevs; lev++)
            offset = MyPC->sumParticleMass(0,lev);
        if (ParallelDescriptor::IOProcessor())
           std::cout << "Total charge of particles = " << offset << std::endl;
        offset /= geom[0].ProbSize();
    }

    // **************************************************************************

    // Define the density on level 0 from all particles at all levels
    int base_level   = 0;
    int finest_level = nlevs-1;

    Vector<std::unique_ptr<MultiFab> > PartMF(nlevs);
    PartMF[0].reset(new MultiFab(ba[0],dmap[0],1,1));
    PartMF[0]->setVal(0.0);

//  MyPC->AssignDensity(0, PartMF, false, base_level, 1, finest_level);
    MyPC->AssignCellDensitySingleLevel(0, *PartMF[0], 0, 1, 0);

    for (int lev = finest_level - 1 - base_level; lev >= 0; lev--)
        amrex::average_down(*PartMF[lev+1],*PartMF[lev],0,1,rr[lev]);

    for (int lev = 0; lev < nlevs; lev++)
        MultiFab::Add(*rhs[base_level+lev], *PartMF[lev], 0, 0, 1, 0);
 
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
    solve_for_accel(amrex::GetVecOfPtrs(rhs),
		    amrex::GetVecOfPtrs(phi),
		    amrex::GetVecOfPtrs(grad_phi),
		    geom,base_level,finest_level,offset);

    // Fill the particle data with the acceleration at the particle location
    // Note that we are calling moveKick with accel_comp > BL_SPACEDIM
    //      which means with dt = 0 we don't move the particle or set a velocity
    for (int lev = 0; lev < nlevs; lev++)
        MyPC->moveKick(*grad_phi[lev],lev,dummy_dt,1.0,1.0,accel_comp);

    // Write out the positions, masses and accelerations of each particle.
    MyPC->WriteAsciiFile("Particles_after_multilevel_solve");

    } // end if (nlevs > 1)

    }
    
    amrex::Finalize();
}

Real
getEfficiency(const DistributionMapping& dm, const Vector<long>& cost)
{
    Vector<long> cpr(ParallelDescriptor::NProcs(), 0);
    Real ctot=0;
    for (int i=0, N=cost.size(); i<N; i++) {
        ctot += cost[i];
        cpr[dm[i]] += cost[i];
    }
    long cmax = *std::max_element(cpr.begin(), cpr.end());
    Real cavg = ctot / ParallelDescriptor::NProcs();
    return cavg / cmax;
}

DistributionMapping
getCostCountDM (const Vector<long>& cost, const BoxArray& ba)
{
    DistributionMapping res;
    int nprocs = ParallelDescriptor::NProcs();
    const int factor = 1.5; // A process can get up to 'factor' times of the average number of boxes.
    int nmax = (cost.size()+nprocs-1) / nprocs * factor;
    Real eff;
    res.KnapSackProcessorMap(cost, nprocs, &eff, true, nmax);
    return res;
}
 
 
