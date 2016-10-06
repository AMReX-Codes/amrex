#include <iostream>

#include <BoxLib.H>
#include <MultiFab.H>
#include <MultiFabUtil.H>
#include <BLFort.H>
#include <MultiFabUtil.H>

#include "ParticleContainer.H"

namespace {
    constexpr Real clight = 299792458.;
};

extern "C" {
  void pxrpush_em3d_evec_norder( Real* ex, Real* ey, Real* ez,
	Real* bx, Real* by, Real* bz,
	const Real* jx, const Real* jy, const Real* jz,
	const Real* mudt,
	const Real* dtsdx, const Real* dtsdy, const Real* dtsdz, 
	const long* nxlocal, const long* nylocal, const long* nzlocal,
	const long* nxorder, const long* nyorder, const long* nzorder,
	const long* nxguard, const long* nyguard, const long* nzguard, 
	const long* nxs, const long* nys, const long* nzs, const int* l_nodal );

  void pxrpush_em3d_bvec_norder( Real* ex, Real* ey, Real* ez,
	Real* bx, Real* by, Real* bz,
	const Real* dtsdx, const Real* dtsdy, const Real* dtsdz, 
	const long* nxlocal, const long* nylocal, const long* nzlocal,
	const long* nxorder, const long* nyorder, const long* nzorder,
	const long* nxguard, const long* nyguard, const long* nzguard, 
	const long* nxs, const long* nys, const long* nzs, const int* l_nodal );
}


void
single_level(int nlevs, int nx, int ny, int nz, int max_grid_size, int order, bool verbose) 
{
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    Real strt_init, strt_assb, strt_assd, strt_assc, strt_mK, strt_evec, strt_bvec;
    Real  end_init,  end_assb,  end_assd,  end_assc,  end_mK,  end_evec,  end_bvec;

    // ********************************************************************************************
    // All of this defines the level 0 information -- size of box, type of boundary condition, etc.
    // ********************************************************************************************

    BL_ASSERT(nlevs == 1);

    Array<Geometry> geom(1);
    Array<DistributionMapping> dmap(1);
    Array<BoxArray> ba(1);
    Real dx, dy, dz;
    Real dt;
    {
	// This defines the physical size of the box. 
	// Right now the box is [0,1] in each direction.
	RealBox real_box;
	for (int n = 0; n < BL_SPACEDIM; n++)
	{
	    real_box.setLo(n,0.0);
	    real_box.setHi(n,1.0);
	}
	dx = real_box.length(0)/nx;
	dy = real_box.length(1)/ny;
	dz = real_box.length(2)/nz;
	dt  = 0.95 * dx /clight ;
	
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
	geom[0].define(domain, &real_box, coord, is_per);
	
	// Build a BoxArray then initialize with the domain.
	ba[0].define(domain);
	
	// Break the BoxArrays at both levels into max_grid_size^3 boxes
	ba[0].maxSize(max_grid_size);

	// Build a dummy MultiFab so we can use its DistributionMap
	MultiFab dummyMF(ba[0],1,0);
	dmap[0] = dummyMF.DistributionMap();
    }

    if (verbose && ParallelDescriptor::IOProcessor())
    {
       std::cout << "Number of boxes              : " << ba[0].size() << '\n' << '\n';
    }

    // ********************************************************************************************
    // Set up the arrays for the solve
    // ********************************************************************************************

    // We define the refinement ratio even though we are single level because
    //    we want to use the multilevel interface in the different calls.
    Array<int> rr(nlevs-1);
    
    // Build a new particle container to hold my particles.
    MyParticleContainer MyPC(geom,dmap,ba,rr);

    MyPC.SetVerbose(0);

    // This allows us to write the gravitational acceleration into these components 

    int num_particles = nx * ny * nz;

    if (ParallelDescriptor::IOProcessor())
       std::cout << "Total number of particles    : " << num_particles << '\n' << '\n';

    // **************************************************************************
    // Do a single-level solve on level 0
    // **************************************************************************

    strt_init = ParallelDescriptor::second();

    // Initialize particles and distribute using the BoxArray and DistributionMap from dummy_mf
    {
	MultiFab dummy_mf(ba[0],1,0);
	MyPC.Init(dummy_mf);
    }

    end_init = ParallelDescriptor::second() - strt_init;

    // Write out the positions, masses and accelerations of each particle.
    if (verbose) MyPC.WriteAsciiFile("Particles_before");

    // **************************************************************************

    IntVect nodal(1,1,1);
    int ngrow = order-1;
    MultiFab ChargeMF(ba[0],1,ngrow,Fab_allocate,nodal);

    // **************************************************************************
    // First we test the PICSAR charge deposition
    // **************************************************************************

    strt_assd = ParallelDescriptor::second();

    // Initialize to zero
    ChargeMF.setVal(0.0);
    
    // Charge deposition
    int level = 0;
    MyPC.ChargeDeposition(ChargeMF,level,order); 

    end_assd = ParallelDescriptor::second() - strt_assd;

    if (verbose)
    {
	Real cmin = ChargeMF.min(0,0);
	Real cmax = ChargeMF.max(0.0);
	if (ParallelDescriptor::IOProcessor()) {
	    std::cout << "PICSAR:Min of ChargeMF " << cmin << std::endl;
	    std::cout << "PICSAR:Max of ChargeMF " << cmax << std::endl;
	    std::cout << " " << std::endl;
	}
    }

    // **************************************************************************
    // Next we test the BoxLib charge deposition
    // **************************************************************************

    strt_assb = ParallelDescriptor::second();

    // Initialize to zero
    ChargeMF.setVal(0.0);
    
    // Charge deposition
    MyPC.NodalDepositionSingleLevel(ChargeMF,0,1); 

    end_assb = ParallelDescriptor::second() - strt_assb;

    if (verbose)
    {
	Real cmin = ChargeMF.min(0,0);
	Real cmax = ChargeMF.max(0.0);
	if (ParallelDescriptor::IOProcessor()) {
	    std::cout << "BoxLib:Min of ChargeMF " << cmin << std::endl;
	    std::cout << "BoxLib:Max of ChargeMF " << cmax << std::endl;
	    std::cout << " " << std::endl;
	}
    }

    // **************************************************************************
    // Now we test the PICSAR current deposition - 
    //     these are on the centers of edge but we declare the arrays as nodal
    // **************************************************************************

    PArray<MultiFab> CurrentMF(BL_SPACEDIM,PArrayManage);

    CurrentMF.set(0,new MultiFab(ba[0],1,1,Fab_allocate,nodal));
    CurrentMF.set(1,new MultiFab(ba[0],1,1,Fab_allocate,nodal));
    CurrentMF.set(2,new MultiFab(ba[0],1,1,Fab_allocate,nodal));

    CurrentMF[0].setVal(0.0);
    CurrentMF[1].setVal(0.0);
    CurrentMF[2].setVal(0.0);

    strt_assc = ParallelDescriptor::second();

    // Current deposition
    MyPC.CurrentDeposition(CurrentMF,0,dt); 

    end_assc = ParallelDescriptor::second() - strt_assc;

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "Time in CurrentDeposition      : " << end_assc << '\n';
        std::cout << " " << std::endl;
    }

    // **************************************************************************
    // Create the B arrays -- 
    //     these are on the centers of faces but we declare the arrays as nodal
    // **************************************************************************

    PArray<MultiFab> BfieldMF(BL_SPACEDIM,PArrayManage);

    BfieldMF.set(0,new MultiFab(ba[0],1,1,Fab_allocate,nodal));
    BfieldMF.set(1,new MultiFab(ba[0],1,1,Fab_allocate,nodal));
    BfieldMF.set(2,new MultiFab(ba[0],1,1,Fab_allocate,nodal));

    BfieldMF[0].setVal(0.0);
    BfieldMF[1].setVal(0.0);
    BfieldMF[2].setVal(0.0);

    // **************************************************************************
    // Create the E arrays -- 
    //     these are on the centers of edges but we declare the arrays as nodal
    // **************************************************************************

    PArray<MultiFab> EfieldMF(BL_SPACEDIM,PArrayManage);

    EfieldMF.set(0,new MultiFab(ba[0],1,1,Fab_allocate,nodal));
    EfieldMF.set(1,new MultiFab(ba[0],1,1,Fab_allocate,nodal));
    EfieldMF.set(2,new MultiFab(ba[0],1,1,Fab_allocate,nodal));

    EfieldMF[0].setVal(0.0);
    EfieldMF[1].setVal(0.0);
    EfieldMF[2].setVal(0.0);
      
    Real mu0 = 1.2566370614359173e-06;
    Real mu_c2_dt = mu0*clight*clight * dt;
    /* Define coefficients of the stencil: for now, only the Yee grid */
    Real dtsdx[1], dtsdy[1], dtsdz[1], dtsdx_c2[1], dtsdy_c2[1],
    dtsdz_c2[1];
    dtsdx[0] = dt/dx;
    dtsdy[0] = dt/dy;
    dtsdz[0] = dt/dz;
    dtsdx_c2[0] = dt/dx * clight*clight;
    dtsdy_c2[0] = dt/dy * clight*clight;
    dtsdz_c2[0] = dt/dz * clight*clight;
    long norder = 2;
    long nguard = 1;
    int l_nodal = 0;

    strt_evec = ParallelDescriptor::second();

    for ( MFIter mfi(EfieldMF[0]); mfi.isValid(); ++mfi ) {
      const Box& bx = mfi.validbox();
      long nxlocal = bx.length(0)-1, nylocal = bx.length(1)-1, nzlocal = bx.length(2)-1; 
      
      pxrpush_em3d_evec_norder( EfieldMF[0][mfi].dataPtr(),
			        EfieldMF[1][mfi].dataPtr(),
  			        EfieldMF[2][mfi].dataPtr(),
  			        BfieldMF[0][mfi].dataPtr(),
  			        BfieldMF[1][mfi].dataPtr(),
			        BfieldMF[2][mfi].dataPtr(), 
			        CurrentMF[0][mfi].dataPtr(),
			        CurrentMF[1][mfi].dataPtr(),
			        CurrentMF[2][mfi].dataPtr(),
			        &mu_c2_dt, dtsdx_c2, dtsdy_c2, dtsdz_c2, 
			        &nxlocal, &nylocal, &nzlocal,
			        &norder, &norder, &norder,
			        &nguard, &nguard, &nguard,
			        &nguard, &nguard, &nguard,
			        &l_nodal );
    }

    end_evec = ParallelDescriptor::second() - strt_evec;

    strt_bvec = ParallelDescriptor::second();

    for ( MFIter mfi(EfieldMF[0]); mfi.isValid(); ++mfi ) {
      const Box& bx = mfi.validbox();
      long nxlocal = bx.length(0)-1, nylocal = bx.length(1)-1, nzlocal = bx.length(2)-1; 
      
      pxrpush_em3d_bvec_norder( EfieldMF[0][mfi].dataPtr(),
			        EfieldMF[1][mfi].dataPtr(),
  			        EfieldMF[2][mfi].dataPtr(),
  			        BfieldMF[0][mfi].dataPtr(),
  			        BfieldMF[1][mfi].dataPtr(),
			        BfieldMF[2][mfi].dataPtr(), 
			        dtsdx, dtsdy, dtsdz, 
			        &nxlocal, &nylocal, &nzlocal,
			        &norder, &norder, &norder,
			        &nguard, &nguard, &nguard,
			        &nguard, &nguard, &nguard,
			        &l_nodal );
    }

    end_bvec = ParallelDescriptor::second() - strt_bvec;
    
    // **************************************************************************
    // Now we scatter E and B onto particles
    // **************************************************************************

    long gather_order     = 1;
    long field_gathe_algo = 1;

    // Fill the particle data with E, B at the particle locations
    MyPC.FieldGather( EfieldMF[0], EfieldMF[1], EfieldMF[2],
                       BfieldMF[0], BfieldMF[1], BfieldMF[2],
                       gather_order,field_gathe_algo);

    // **************************************************************************
    // Now we scatter E and B onto particles
    // **************************************************************************

    // Note that we are calling moveKick with accel_comp > BL_SPACEDIM
    //      which means with dt = 0 we don't move the particle or set a velocity

    strt_mK = ParallelDescriptor::second();

    MyPC.ParticlePush(dt);

    end_mK = ParallelDescriptor::second() - strt_mK;

    // After the particles' positions have been updated, update which grid and process own them.
    MyPC.Redistribute();

    // Write out the positions, masses and accelerations of each particle.
    if (verbose) MyPC.WriteAsciiFile("Particles_after");

    ParallelDescriptor::ReduceRealMax(end_init,IOProc);
    ParallelDescriptor::ReduceRealMax(end_assb,IOProc);
    ParallelDescriptor::ReduceRealMax(end_assc,IOProc);
    ParallelDescriptor::ReduceRealMax(end_assd,IOProc);
    ParallelDescriptor::ReduceRealMax(end_evec,IOProc);
    ParallelDescriptor::ReduceRealMax(end_bvec,IOProc);
    ParallelDescriptor::ReduceRealMax(end_mK  ,IOProc);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
           std::cout << "Time in Init                   : " << end_init << '\n';
           std::cout << "Time in BoxLibChargeDeposition : " << end_assb << '\n';
           std::cout << "Time in PicsarChargeDeposition : " << end_assd << '\n';
           std::cout << "Time in CurrentDeposition      : " << end_assc << '\n';
           std::cout << "Time in E-field                : " << end_evec << '\n';
           std::cout << "Time in B-field                : " << end_bvec << '\n';
           std::cout << "Time in moveKick               : " << end_mK   << '\n';
    }
}
