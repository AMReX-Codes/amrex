
#include <random>

#include <BoxLib.H>
#include <ParmParse.H>
#include <Array.H>
#include <MultiFab.H>
#include <PlotFileUtil.H>

#include <WarpXConst.H>
#include <WarpX_f.H>

int main(int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    {
	long nox=1, noy=1, noz=1;
	{
	    ParmParse pp("interpolation");
	    pp.query("nox", nox);
	    pp.query("noy", noy);
	    pp.query("noz", noz);  
	    if (nox != noy || nox != noz) {
		BoxLib::Abort("warpx.nox, noy and noz must be equal");
	    }
	    if (nox < 1) {
		BoxLib::Abort("warpx.nox must >= 1");
	    }
	}

	std::mt19937 rand_eng(42);
	std::uniform_real_distribution<Real> rand_dis(0.0,100.0);

	long nx = 64, ny = 64, nz = 64;
	Real dx[3] = {1.0/nx, 1.0/ny, 1.0/nz};
	Real dt = 1.e-6;

	Box cc_domain{IntVect{D_DECL(0,0,0)},IntVect{D_DECL(nx-1,ny-1,nz-1)}};
	BoxArray grids{cc_domain};
	grids.maxSize(32);

	const int ng = nox;

	Array<std::unique_ptr<MultiFab> > current(3);
	Array<std::unique_ptr<MultiFab> > Efield(3);
	Array<std::unique_ptr<MultiFab> > Bfield(3);
	for (int i = 0; i < 3; ++i) {
	    const IntVect& nodalflag = IntVect::TheUnitVector();
	    current[i].reset(new MultiFab(grids,1,ng,Fab_allocate,nodalflag));
	    Efield [i].reset(new MultiFab(grids,1,ng,Fab_allocate,nodalflag));
	    Bfield [i].reset(new MultiFab(grids,1,ng,Fab_allocate,nodalflag));
	}

	for (MFIter mfi(*current[0]); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.fabbox();
	    FArrayBox& exfab = (*Efield [0])[mfi];
	    FArrayBox& eyfab = (*Efield [1])[mfi];
	    FArrayBox& ezfab = (*Efield [2])[mfi];
	    FArrayBox& bxfab = (*Bfield [0])[mfi];
	    FArrayBox& byfab = (*Bfield [1])[mfi];
	    FArrayBox& bzfab = (*Bfield [2])[mfi];
	    FArrayBox& jxfab = (*current[0])[mfi];
	    FArrayBox& jyfab = (*current[1])[mfi];
	    FArrayBox& jzfab = (*current[2])[mfi];
	    for (IntVect cell=bx.smallEnd(); cell <= bx.bigEnd(); bx.next(cell))
	    {
		exfab(cell) = rand_dis(rand_eng);
		eyfab(cell) = rand_dis(rand_eng);
		ezfab(cell) = rand_dis(rand_eng);
		bxfab(cell) = rand_dis(rand_eng);
		byfab(cell) = rand_dis(rand_eng);
		bzfab(cell) = rand_dis(rand_eng);
		jxfab(cell) = rand_dis(rand_eng);
		jyfab(cell) = rand_dis(rand_eng);
		jzfab(cell) = rand_dis(rand_eng);
	    }
	}


	{ // Evolve B
	    Real dtsdx[3];
#if (BL_SPACEDIM == 3)
	    dtsdx[0] = dt / dx[0];
	    dtsdx[1] = dt / dx[1];
	    dtsdx[2] = dt / dx[2];
#elif (BL_SPACEDIM == 2)
	    dtsdx[0] = dt / dx[0];
	    dtsdx[1] = std::numeric_limits<Real>::quiet_NaN();
	    dtsdx[2] = dt / dx[1];
#endif

	    long norder = 2;
	    long nstart = 0;
	    int l_nodal = false;
	    long nguard = Efield[0]->nGrow();

#if (BL_SPACEDIM == 3)
	    long nxguard = nguard;
	    long nyguard = nguard;
	    long nzguard = nguard; 
#elif (BL_SPACEDIM == 2)
	    long nxguard = nguard;
	    long nyguard = 0;
	    long nzguard = nguard; 
#endif

	    for ( MFIter mfi(*Bfield[0]); mfi.isValid(); ++mfi )
	    {
		const Box& bx = BoxLib::enclosedCells(mfi.validbox());
#if (BL_SPACEDIM == 3)
		long nx = bx.length(0);
		long ny = bx.length(1);
		long nz = bx.length(2); 
#elif (BL_SPACEDIM == 2)
		long nx = bx.length(0);
		long ny = 0;
		long nz = bx.length(1); 
#endif

		warpx_push_bvec( (*Efield[0])[mfi].dataPtr(),
				 (*Efield[1])[mfi].dataPtr(),
				 (*Efield[2])[mfi].dataPtr(),
				 (*Bfield[0])[mfi].dataPtr(),
				 (*Bfield[1])[mfi].dataPtr(),
				 (*Bfield[2])[mfi].dataPtr(), 
				 dtsdx, dtsdx+1, dtsdx+2,
				 &nx, &ny, &nz,
				 &norder, &norder, &norder,
				 &nxguard, &nyguard, &nzguard,
				 &nstart, &nstart, &nstart,
				 &l_nodal );
	    }	    
	}

	{ // evolve E
	    Real mu_c2_dt = (PhysConst::mu0*PhysConst::c*PhysConst::c) * dt;

	    Real dtsdx_c2[3];
#if (BL_SPACEDIM == 3)
	    dtsdx_c2[0] = (PhysConst::c*PhysConst::c) * dt / dx[0];
	    dtsdx_c2[1] = (PhysConst::c*PhysConst::c) * dt / dx[1];
	    dtsdx_c2[2] = (PhysConst::c*PhysConst::c) * dt / dx[2];
#else
	    dtsdx_c2[0] = (PhysConst::c*PhysConst::c) * dt / dx[0];
	    dtsdx_c2[1] = std::numeric_limits<Real>::quiet_NaN();
	    dtsdx_c2[2] = (PhysConst::c*PhysConst::c) * dt / dx[1];
#endif

	    long norder = 2;
	    long nstart = 0;
	    int l_nodal = false;	    
	    long nguard = Efield[0]->nGrow();

#if (BL_SPACEDIM == 3)
	    long nxguard = nguard;
	    long nyguard = nguard;
	    long nzguard = nguard; 
#elif (BL_SPACEDIM == 2)
	    long nxguard = nguard;
	    long nyguard = 0;
	    long nzguard = nguard; 
#endif

	    for ( MFIter mfi(*Efield[0]); mfi.isValid(); ++mfi )
	    {
		const Box & bx = BoxLib::enclosedCells(mfi.validbox());
#if (BL_SPACEDIM == 3)
		long nx = bx.length(0);
		long ny = bx.length(1);
		long nz = bx.length(2); 
#elif (BL_SPACEDIM == 2)
		long nx = bx.length(0);
		long ny = 0;
		long nz = bx.length(1); 
#endif

		warpx_push_evec( (*Efield[0])[mfi].dataPtr(),
				 (*Efield[1])[mfi].dataPtr(),
				 (*Efield[2])[mfi].dataPtr(),
				 (*Bfield[0])[mfi].dataPtr(),
				 (*Bfield[1])[mfi].dataPtr(),
				 (*Bfield[2])[mfi].dataPtr(), 
				 (*current[0])[mfi].dataPtr(),
				 (*current[1])[mfi].dataPtr(),
				 (*current[2])[mfi].dataPtr(),
				 &mu_c2_dt, dtsdx_c2, dtsdx_c2+1, dtsdx_c2+2,
				 &nx, &ny, &nz,
				 &norder, &norder, &norder,
				 &nxguard, &nyguard, &nzguard,
				 &nstart, &nstart, &nstart,
				 &l_nodal );
	    }
	}

	MultiFab plotmf(grids, 6, 0);
	const auto& src_typ = (*Efield[0]).boxArray().ixType();
	const auto& dst_typ = grids.ixType();
	for (MFIter mfi(plotmf); mfi.isValid(); ++mfi)
	{
	    FArrayBox& dfab = plotmf[mfi];
	    dfab.SetBoxType(src_typ);
	    dfab.copy((*Efield[0])[mfi], 0, 0, 1);
	    dfab.copy((*Efield[1])[mfi], 0, 1, 1);
	    dfab.copy((*Efield[2])[mfi], 0, 2, 1);
	    dfab.copy((*Bfield[0])[mfi], 0, 3, 1);
	    dfab.copy((*Bfield[1])[mfi], 0, 4, 1);
	    dfab.copy((*Bfield[2])[mfi], 0, 5, 1);
	    dfab.SetBoxType(dst_typ);
	}

	Real xyzmin[3] = {0.0,0.0,0.0};
	RealBox realbox{cc_domain, dx, xyzmin};
	int is_per[3] = {0,0,0};
	Geometry geom{cc_domain, &realbox, 0, is_per};
	std::string plotname{"plt00000"};
	Array<std::string> varnames{"Ex", "Ey", "Ez", "Bx", "By", "Bz"};
	BoxLib::WriteSingleLevelPlotfile(plotname, plotmf, varnames, geom, 0.0, 0);
    }

    BoxLib::Finalize();
}
