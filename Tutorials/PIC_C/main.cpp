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
void solve_for_accel(PArray<MultiFab>& rhs, PArray<MultiFab>& phi, PArray<MultiFab>& grad_phi, 
                     const PArray<Geometry>& geom, int base_level);

int main(int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    // Define this as a 2-level problem.
    int nlevs = 2;

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
    IntVect domain_hi(31,31,31); 
//  IntVect domain_hi( 7, 7, 7); 
 
    // Build a box for the level 0 domain
    Box domain(domain_lo, domain_hi);

    // This says we are using Cartesian coordinates
    int coord = 0;

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) is_per[i] = 1; 

    // Define the refinement ratio
    int rr[nlevs-1];
    rr[0] = 2;

    // This defines a Geometry object which is useful for writing the plotfiles  
    PArray<Geometry> geom;
    geom.resize(nlevs,PArrayManage);
    for (int n = 0; n < nlevs; n++)
    {
       geom.set(n, new Geometry(domain,&real_box,coord,is_per));
       if (n < nlevs-1) domain.refine(rr[n]);
    }

    // We will start each solve at level 0.
    int base_level = 0;

    // ********************************************************************************************
    // This now defines the refinement information 
    // ********************************************************************************************

    // For now we make the refined level be the center eighth of the domain
    IntVect refined_lo(15,15,15); 
    IntVect refined_hi(47,47,47); 

    // Build a box for the level 1 domain
    Box refined_patch(refined_lo, refined_hi);

    // Build an array of BoxArrays,
    // then initialize the level 0 BoxArray with the domain.
    PArray<BoxArray> ba;
    ba.resize(nlevs);
    ba.set(0,new BoxArray(domain));
    ba.set(1,new BoxArray(refined_patch));

    // break the BoxArrays at both levels into 32^3 boxes
    for (int n = 0; n < nlevs; n++)
        ba[n].maxSize(64);

    // build a multifab for the rhs on the box array with 1 component, 0 ghost cells
    PArray<MultiFab> rhs; 
    PArray<MultiFab> phi; 
    PArray<MultiFab> grad_phi;
    Array<DistributionMapping> dmap;

    rhs.resize(nlevs,PArrayManage);
    phi.resize(nlevs,PArrayManage);
    grad_phi.resize(nlevs,PArrayManage);
    dmap.resize(nlevs);

    for (int lev = 0; lev < nlevs; lev++)
    {
       rhs.set(lev,new MultiFab(ba[lev],1,0));
       phi.set(lev,new MultiFab(ba[lev],1,1));
       dmap[lev] = rhs[lev].DistributionMap();
    }

    // build a multifab for grad_phi on the box array with BL_SPACEDIM components, 1 ghost cell
    for (int n = 0; n < nlevs; n++)
    {
       grad_phi.set(n, new MultiFab(ba[n], BL_SPACEDIM, 1));
       grad_phi[n].setVal(0.0);
    }

    // Define a new particle container to hold my particles.
    // This holds a charge as well as three velocity components, three acceleration components  and three position components.
    typedef ParticleContainer<1+2*BL_SPACEDIM> MyParticleContainer;
    
    // Build a new particle container to hold my particles.
    MyParticleContainer* MyPC = new MyParticleContainer(geom,dmap,ba,rr);

    int count = 1;
    int iseed = 10;
    Real mass  = 10.0;

    // Initialize "count" number of particles, each with mass/charge "mass"
    MyPC->InitRandom(count,iseed,mass);
    MyPC->Redistribute();

    // Use the PIC approach to deposit the "mass" onto the grid
    MyPC->AssignDensitySingleLevel(rhs[0],0);

    // std::cout << "RHS " << rhs[0][0] << std::endl;

    // Write out the positions, masses and accelerations of each particle.
    MyPC->WriteAsciiFile("Particles_before");

    // Use multigrid to solve Lap(phi) = rhs with periodic boundary conditions (set above)
    solve_for_accel(rhs,phi,grad_phi,geom,base_level);

    // Fill the particle data with the acceleration at the particle location
    int start_comp = BL_SPACEDIM+1;
    Real dummy_dt = 0.0;
    MyPC->moveKick(grad_phi[0],0,dummy_dt,1.0,1.0,start_comp);

    // Write out the positions, masses and accelerations of each particle.
    MyPC->WriteAsciiFile("Particles_after");

    BoxLib::Finalize();
}
