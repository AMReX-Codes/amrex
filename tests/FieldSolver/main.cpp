
#include <limits>
#include <random>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

#include <WarpXConst.H>
#include <WarpX_f.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
	long nox=1, noy=1, noz=1;
	{
	    ParmParse pp("interpolation");
	    pp.query("nox", nox);
	    pp.query("noy", noy);
	    pp.query("noz", noz);  
	    if (nox != noy || nox != noz) {
		amrex::Abort("warpx.nox, noy and noz must be equal");
	    }
	    if (nox < 1) {
		amrex::Abort("warpx.nox must >= 1");
	    }
	}

	std::mt19937 rand_eng(42);
	std::uniform_real_distribution<Real> rand_dis(0.0,1.0);

	long nx = 64, ny = 64, nz = 64;
	Real dx[3] = {1.0/nx, 1.0/ny, 1.0/nz};
	Real dt = 1.e-6;

	Box cc_domain{IntVect{D_DECL(0,0,0)},IntVect{D_DECL(nx-1,ny-1,nz-1)}};
	BoxArray grids{cc_domain};
	grids.maxSize(32);

	DistributionMapping dmap {grids};

	const int ng = nox;

#if (BL_SPACEDIM == 3)
        IntVect Bx_nodal_flag(1,0,0);
        IntVect By_nodal_flag(0,1,0);
        IntVect Bz_nodal_flag(0,0,1);
#elif (BL_SPACEDIM == 2)
        IntVect Bx_nodal_flag(1,0);  // x is the first dimension to AMReX
        IntVect By_nodal_flag(0,0);  // y is the missing dimension to 2D AMReX
        IntVect Bz_nodal_flag(0,1);  // z is the second dimension to 2D AMReX
#endif
        
#if (BL_SPACEDIM == 3)
        IntVect Ex_nodal_flag(0,1,1);
        IntVect Ey_nodal_flag(1,0,1);
        IntVect Ez_nodal_flag(1,1,0);
#elif (BL_SPACEDIM == 2)
        IntVect Ex_nodal_flag(0,1);  // x is the first dimension to AMReX
        IntVect Ey_nodal_flag(1,1);  // y is the missing dimension to 2D AMReX
        IntVect Ez_nodal_flag(1,0);  // z is the second dimension to 2D AMReX
#endif
        
#if (BL_SPACEDIM == 3)
        IntVect jx_nodal_flag(0,1,1);
        IntVect jy_nodal_flag(1,0,1);
        IntVect jz_nodal_flag(1,1,0);
#elif (BL_SPACEDIM == 2)
        IntVect jx_nodal_flag(0,1);  // x is the first dimension to AMReX
        IntVect jy_nodal_flag(1,1);  // y is the missing dimension to 2D AMReX
        IntVect jz_nodal_flag(1,0);  // z is the second dimension to 2D AMReX
#endif

	Vector<std::unique_ptr<MultiFab> > current(3);
	Vector<std::unique_ptr<MultiFab> > Efield(3);
	Vector<std::unique_ptr<MultiFab> > Bfield(3);
        // Create the MultiFabs for B
        Bfield[0].reset( new MultiFab(amrex::convert(grids,Bx_nodal_flag),dmap,1,ng));
        Bfield[1].reset( new MultiFab(amrex::convert(grids,By_nodal_flag),dmap,1,ng));
        Bfield[2].reset( new MultiFab(amrex::convert(grids,Bz_nodal_flag),dmap,1,ng));
        
        // Create the MultiFabs for E
        Efield[0].reset( new MultiFab(amrex::convert(grids,Ex_nodal_flag),dmap,1,ng));
        Efield[1].reset( new MultiFab(amrex::convert(grids,Ey_nodal_flag),dmap,1,ng));
        Efield[2].reset( new MultiFab(amrex::convert(grids,Ez_nodal_flag),dmap,1,ng));
        
        // Create the MultiFabs for the current
        current[0].reset( new MultiFab(amrex::convert(grids,jx_nodal_flag),dmap,1,ng));
        current[1].reset( new MultiFab(amrex::convert(grids,jy_nodal_flag),dmap,1,ng));
        current[2].reset( new MultiFab(amrex::convert(grids,jz_nodal_flag),dmap,1,ng));
        
        Bfield[0]->setVal(0.0);
        Bfield[1]->setVal(0.0);
        Bfield[2]->setVal(0.0);
        Efield[0]->setVal(0.0);
        Efield[1]->setVal(0.0);
        Efield[2]->setVal(0.0);
        current[0]->setVal(0.0);
        current[1]->setVal(0.0);
        current[2]->setVal(0.0);

	for (MFIter mfi(*current[0]); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.validbox();
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
		exfab(cell) = rand_dis(rand_eng)*1.e5;
		eyfab(cell) = rand_dis(rand_eng)*1.e5;
		ezfab(cell) = rand_dis(rand_eng)*1.e5;
		bxfab(cell) = rand_dis(rand_eng)*1.e-5;
		byfab(cell) = rand_dis(rand_eng)*1.e-5;
		bzfab(cell) = rand_dis(rand_eng)*1.e-5;
		jxfab(cell) = rand_dis(rand_eng)*1.e10;
		jyfab(cell) = rand_dis(rand_eng)*1.e10;
		jzfab(cell) = rand_dis(rand_eng)*1.e10;
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

	    const int norder = 2;

	    for ( MFIter mfi(*Bfield[0],true); mfi.isValid(); ++mfi )
	    {
                const Box& tbx  = mfi.tilebox(Bx_nodal_flag);
                const Box& tby  = mfi.tilebox(By_nodal_flag);
                const Box& tbz  = mfi.tilebox(Bz_nodal_flag);
                
                // Call picsar routine for each tile
                WRPX_PXR_PUSH_BVEC(
                    tbx.loVect(), tbx.hiVect(),
                    tby.loVect(), tby.hiVect(),
                    tbz.loVect(), tbz.hiVect(),
                    BL_TO_FORTRAN_3D((*Efield[0])[mfi]),
                    BL_TO_FORTRAN_3D((*Efield[1])[mfi]),
                    BL_TO_FORTRAN_3D((*Efield[2])[mfi]),
                    BL_TO_FORTRAN_3D((*Bfield[0])[mfi]),
                    BL_TO_FORTRAN_3D((*Bfield[1])[mfi]),
                    BL_TO_FORTRAN_3D((*Bfield[2])[mfi]),
                    &dtsdx[0], &dtsdx[1], &dtsdx[2],
                    &norder);
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

	    const int norder = 2;

	    for ( MFIter mfi(*Efield[0],true); mfi.isValid(); ++mfi )
	    {
                const Box& tex  = mfi.tilebox(Ex_nodal_flag);
                const Box& tey  = mfi.tilebox(Ey_nodal_flag);
                const Box& tez  = mfi.tilebox(Ez_nodal_flag);
                
                // Call picsar routine for each tile
                WRPX_PXR_PUSH_EVEC(
                    tex.loVect(), tex.hiVect(),
                    tey.loVect(), tey.hiVect(),
                    tez.loVect(), tez.hiVect(),
                    BL_TO_FORTRAN_3D((*Efield[0])[mfi]),
                    BL_TO_FORTRAN_3D((*Efield[1])[mfi]),
                    BL_TO_FORTRAN_3D((*Efield[2])[mfi]),
                    BL_TO_FORTRAN_3D((*Bfield[0])[mfi]),
                    BL_TO_FORTRAN_3D((*Bfield[1])[mfi]),
                    BL_TO_FORTRAN_3D((*Bfield[2])[mfi]),
                    BL_TO_FORTRAN_3D((*current[0])[mfi]),
                    BL_TO_FORTRAN_3D((*current[1])[mfi]),
                    BL_TO_FORTRAN_3D((*current[2])[mfi]),
                    &mu_c2_dt,
                    &dtsdx_c2[0], &dtsdx_c2[1], &dtsdx_c2[2],
                    &norder);
	    }
	}

	MultiFab plotmf(grids, dmap, 6, 0);
        amrex::average_edge_to_cellcenter(plotmf, 0, {Efield[0].get(),Efield[1].get(),Efield[2].get()});
        amrex::average_face_to_cellcenter(plotmf, 3, {Bfield[0].get(),Bfield[1].get(),Bfield[2].get()});

	Real xyzmin[3] = {0.0,0.0,0.0};
	RealBox realbox{cc_domain, dx, xyzmin};
	int is_per[3] = {0,0,0};
	Geometry geom{cc_domain, &realbox, 0, is_per};
	std::string plotname{"plt00000"};
	Vector<std::string> varnames{"Ex", "Ey", "Ez", "Bx", "By", "Bz"};
	amrex::WriteSingleLevelPlotfile(plotname, plotmf, varnames, geom, 0.0, 0);
    }

    amrex::Finalize();
}
