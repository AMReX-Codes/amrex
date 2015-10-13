#include <iostream>

#include <BoxLib.H>
#include <MultiFab.H>
#include <MultiFabUtil.H>
#include <BLFort.H>
#include <MacBndry.H>
#include <MGT_Solver.H>
#include <mg_cpp_f.h>
#include <stencil_types.H>

#include "Particles.H"

// declare routines below
void solve_for_phi(MultiFab& rhs, MultiFab& grad_phi, const Geometry& geom);

int main(int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    // define the lower and upper corner of a 3D domain
    IntVect domain_lo(0 , 0, 0); 
    IntVect domain_hi(63,63,63); 
 
    // build a box for the domain
    Box domain(domain_lo, domain_hi);

    // build a box array from the 64^3 domain box
    BoxArray ba(domain);
    // break the box array into 32^3 boxes
    ba.maxSize(32);

    // build a multifab for the rhs on the box array with 1 component, 0 ghost cells
    MultiFab rhs(ba, 1, 0);  

    // build a multifab for grad_phi on the box array with BL_SPACEDIM components, 1 ghost cell
    MultiFab grad_phi(ba, BL_SPACEDIM, 1);  

    // Initialize phi to 0 everywhere.
    grad_phi.setVal(0.0);

    // This defines the physical size of the box.  Right now the box is [-1,1] in each direction.
    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
       real_box.setLo(n,0.0);
       real_box.setHi(n,1.0);
    }

    // This says we are using Cartesian coordinates
    int coord = 0;

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) is_per[i] = 1; 

    // This defines a Geometry object which is useful for writing the plotfiles  
    Geometry geom(domain,&real_box,coord,is_per);

    DistributionMapping dmap = rhs.DistributionMap();

    // Define a new particle container to hold my particles.
    // This holds a charge as well as three velocity components, three acceleration components  and three position components.
    typedef ParticleContainer<1+2*BL_SPACEDIM> MyParticleContainer;
    
    // Build a new particle container to hold my particles.
    MyParticleContainer* MyPC = new MyParticleContainer(geom,dmap,ba);

    int count = 1;
    int iseed = 10;
    Real mass  = 10.0;

    // Initialize "count" number of particles, each with mass/charge "mass"
    MyPC->InitRandom(count,iseed,mass);

    // Use the PIC approach to deposit the "mass" onto the grid
    MyPC->AssignDensitySingleLevel(rhs,0);

    // Write out the positions, masses and accelerations of each particle.
    MyPC->WriteAsciiFile("Particles_before");

    // Use multigrid to solve Lap(phi) = rhs with periodic boundary conditions (set above)
    solve_for_phi(rhs,grad_phi,geom);

    // Fill the particle data with the acceleration at the particle location
    int start_comp = BL_SPACEDIM+1;
    Real dummy_dt = 0.0;
    MyPC->moveKick(grad_phi,0,dummy_dt,1.0,1.0,start_comp);

    // Write out the positions, masses and accelerations of each particle.
    MyPC->WriteAsciiFile("Particles_after");

    BoxLib::Finalize();
}
