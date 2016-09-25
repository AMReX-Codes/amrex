#include <iostream>

#include <BoxLib.H>
#include <MultiFab.H>
#include <MultiFabUtil.H>
#include <BLFort.H>
#include <MacBndry.H>
#include <MultiFabUtil.H>

#include "Particles.H"

void
single_level(int nlevs, int nx, int ny, int nz, int max_grid_size, int nppc, bool verbose) 
{
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    Real strt_init, strt_assb, strt_assd, strt_assc, strt_mK;
    Real  end_init,  end_assb,  end_assd,  end_assc,  end_mK;

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
    Array<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, coord, is_per);

    // Build a BoxArray then initialize with the domain.
    Array<BoxArray> ba(1);
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

    // Build a dummy MultiFab so we can use its DistributionMap
    MultiFab dummyMF(ba[0],1,0);
    Array<DistributionMapping> dmap(1);
    dmap[0] = dummyMF.DistributionMap();

    // Particle Attributes:
    //    1 -- w  (weight)
    //    2 -- px (x mom)
    //    3 -- py (y mom)
    //    4 -- pz (z mom)
    //    5 -- Ex 
    //    6 -- Ey 
    //    7 -- Ez 
    //    8 -- Bx 
    //    9 -- By 
    //   10 -- Bz 
    const int N_attributes = 10;

    // Define a new particle container to hold my particles.
    typedef ParticleContainer<N_attributes> MyParticleContainer;

    // We define the refinement ratio even though we are single level because
    //    we want to use the multilevel interface in the different calls.
    Array<int> rr(nlevs-1);
    
    // Build a new particle container to hold my particles.
    MyParticleContainer* MyPC = new MyParticleContainer(geom,dmap,ba,rr);

    MyPC->SetVerbose(0);

    // This allows us to write the gravitational acceleration into these components 

    int num_particles = nppc * nx * ny * nz;
    Real charge = 1.0;

    if (ParallelDescriptor::IOProcessor())
       std::cout << "Total number of particles    : " << num_particles << '\n' << '\n';

    // **************************************************************************
    // Do a single-level solve on level 0
    // **************************************************************************

    strt_init = ParallelDescriptor::second();

    // Randomly initialize "num_particles" number of particles, each with charge "charge"
    // bool serialize = false;
    // int iseed   = 10;
    // MyPC->InitRandom(num_particles,iseed,charge,serialize);

    // Initialize one particle at each cell center
    MultiFab dummy_mf(ba[0],1,0,Fab_allocate);
    MyPC->InitOnePerCell(0.5,0.5,0.5,charge,dummy_mf);

    end_init = ParallelDescriptor::second() - strt_init;

    // Write out the positions, masses and accelerations of each particle.
    if (verbose) MyPC->WriteAsciiFile("Particles_before");

    // **************************************************************************

    MultiFab ChargeMF;
    IntVect nodal(1,1,1);
    ChargeMF.define(ba[0],1,0,Fab_allocate,nodal);

    // **************************************************************************
    // First we test the PICSAR charge deposition
    // **************************************************************************

    strt_assd = ParallelDescriptor::second();

    // Initialize to zero
    ChargeMF.setVal(0.0);
    
    // Charge deposition
    MyPC->ChargeDeposition(ChargeMF,0); 

    end_assd = ParallelDescriptor::second() - strt_assd;

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "PICSAR:Min of ChargeMF " << ChargeMF.min(0,0) << std::endl;
        std::cout << "PICSAR:Max of ChargeMF " << ChargeMF.max(0,0) << std::endl;
        std::cout << "Time in PicsarChargeDeposition : " << end_assd << '\n';
        std::cout << " " << std::endl;
    }

    // **************************************************************************
    // Next we test the BoxLib charge deposition
    // **************************************************************************

    strt_assb = ParallelDescriptor::second();

    // Initialize to zero
    ChargeMF.setVal(0.0);
    
    // Charge deposition
    MyPC->AssignNodalDensitySingleLevel(ChargeMF,0,1); 

    end_assb = ParallelDescriptor::second() - strt_assb;

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "BoxLib:Min of ChargeMF " << ChargeMF.min(0,0) << std::endl;
        std::cout << "BoxLib:Max of ChargeMF " << ChargeMF.max(0,0) << std::endl;
        std::cout << "Time in BoxLibChargeDeposition : " << end_assb << '\n';
        std::cout << " " << std::endl;
    }

    // **************************************************************************
    // Define the different spatial centerings
    // **************************************************************************

    // These are the centers of the faces (used for B)
    IntVect x_face (1,0,0);
    IntVect y_face (0,1,0);
    IntVect z_face (0,0,1);

    // These are the centers of the edges (used for E and J)
    IntVect x_edge (0,1,1);
    IntVect y_edge (1,0,1);
    IntVect z_edge (1,1,0);

    // **************************************************************************
    // Now we test the PICSAR current deposition
    // **************************************************************************

    PArray<MultiFab> CurrentMF;
    CurrentMF.resize(BL_SPACEDIM,PArrayManage);

    CurrentMF.set(0,new MultiFab(ba[0],1,0,Fab_allocate,x_edge));
    CurrentMF.set(1,new MultiFab(ba[0],1,0,Fab_allocate,y_edge));
    CurrentMF.set(2,new MultiFab(ba[0],1,0,Fab_allocate,z_edge));

    CurrentMF[0].setVal(0.0);
    CurrentMF[1].setVal(0.0);
    CurrentMF[2].setVal(0.0);

    strt_assc = ParallelDescriptor::second();

    // Current deposition
    Real dummy_dt  = 1.0;
    MyPC->CurrentDeposition(CurrentMF,0,dummy_dt); 

    end_assc = ParallelDescriptor::second() - strt_assc;

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "Time in CurrentDeposition      : " << end_assc << '\n';
        std::cout << " " << std::endl;
    }

    // **************************************************************************
    // Create the B arrays -- these are on the centers of faces
    // **************************************************************************

    PArray<MultiFab> BfieldMF;
    BfieldMF.resize(BL_SPACEDIM,PArrayManage);

    BfieldMF.set(0,new MultiFab(ba[0],1,0,Fab_allocate,x_face));
    BfieldMF.set(1,new MultiFab(ba[0],1,0,Fab_allocate,y_face));
    BfieldMF.set(2,new MultiFab(ba[0],1,0,Fab_allocate,z_face));

    BfieldMF[0].setVal(0.0);
    BfieldMF[1].setVal(0.0);
    BfieldMF[2].setVal(0.0);

    // **************************************************************************
    // Create the E arrays -- these are on the centers of edges
    // **************************************************************************

    PArray<MultiFab> EfieldMF;
    EfieldMF.resize(BL_SPACEDIM,PArrayManage);

    EfieldMF.set(0,new MultiFab(ba[0],1,0,Fab_allocate,x_edge));
    EfieldMF.set(1,new MultiFab(ba[0],1,0,Fab_allocate,y_edge));
    EfieldMF.set(2,new MultiFab(ba[0],1,0,Fab_allocate,z_edge));

    EfieldMF[0].setVal(0.0);
    EfieldMF[1].setVal(0.0);
    EfieldMF[2].setVal(0.0);

    // **************************************************************************
    // Now we create the E and B arrays
    // **************************************************************************

    // Fill the particle data with the acceleration at the particle location
    // Note that we are calling moveKick with accel_comp > BL_SPACEDIM
    //      which means with dt = 0 we don't move the particle or set a velocity

    strt_mK = ParallelDescriptor::second();

    // MyPC->moveKick(grad_phi[0],nlevs-1,dummy_dt,1.0,1.0,accel_comp);

    end_mK = ParallelDescriptor::second() - strt_mK;

    // Write out the positions, masses and accelerations of each particle.
    // if (verbose) MyPC->WriteAsciiFile("Particles_after_level0_solve");

    delete MyPC;

    ParallelDescriptor::ReduceRealMax(end_init,IOProc);
    ParallelDescriptor::ReduceRealMax(end_assd,IOProc);
    ParallelDescriptor::ReduceRealMax(end_assc,IOProc);
    ParallelDescriptor::ReduceRealMax(end_mK  ,IOProc);
    if (verbose && ParallelDescriptor::IOProcessor())
    {
           std::cout << "Time in InitRandom             : " << end_init << '\n';
           std::cout << "Time in BoxLibChargeDeposition : " << end_assb << '\n';
           std::cout << "Time in PicsarChargeDeposition : " << end_assd << '\n';
           std::cout << "Time in CurrentDeposition      : " << end_assc << '\n';
           std::cout << "Time in moveKick               : " << end_mK   << '\n';
    }
}
