
#include <random>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BoxIterator.H>
#include <AMReX_PlotFileUtil.H>

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

	long field_gathering_algo = 1;

	{
	    ParmParse pp("algo");
	    pp.query("field_gathering", field_gathering_algo);
	}

	long nx = 64, ny = 64, nz = 64;
	long np = nx*ny*nz;

	Real xyzmin[3] = {0.5, 1.4, 0.3};
        int ixyzmin[3] = {0, 0, 0};

	Vector<Real> xp, yp, zp, Exp, Eyp, Ezp, Bxp, Byp, Bzp;
	Real dx[3] = {1.0/nx, 1.0/ny, 1.0/nz};

	std::mt19937 rand_eng(42);
	std::uniform_real_distribution<Real> rand_dis(0.0,dx[0]);

	for (int k = 0; k < nz; ++k) {
	    for (int j = 0; j < ny; ++j) {
		for (int i = 0; i < nx; ++i) {
		    Real x0 = xyzmin[0] + i*dx[0];
		    Real y0 = xyzmin[1] + j*dx[1];
		    Real z0 = xyzmin[2] + k*dx[2];
		    xp.push_back(x0 + rand_dis(rand_eng));
		    yp.push_back(y0 + rand_dis(rand_eng));
		    zp.push_back(z0 + rand_dis(rand_eng));
		}
	    }
	}

	Exp.resize(np,0.0);
	Eyp.resize(np,0.0);
	Ezp.resize(np,0.0);
	Bxp.resize(np,0.0);
	Byp.resize(np,0.0);
	Bzp.resize(np,0.0);

	const int ng = nox;
	Box domain_box {IntVect{D_DECL(0,0,0)}, IntVect{D_DECL(nx-1,ny-1,nz-1)}};

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

        BoxArray ba{domain_box};
        DistributionMapping dm{ba};

        MultiFab Bx(amrex::convert(ba,Bx_nodal_flag), dm, 1, ng);
        MultiFab By(amrex::convert(ba,By_nodal_flag), dm, 1, ng);
        MultiFab Bz(amrex::convert(ba,Bz_nodal_flag), dm, 1, ng);
        MultiFab Ex(amrex::convert(ba,Ex_nodal_flag), dm, 1, ng);
        MultiFab Ey(amrex::convert(ba,Ey_nodal_flag), dm, 1, ng);
        MultiFab Ez(amrex::convert(ba,Ez_nodal_flag), dm, 1, ng);

	FArrayBox& exfab = Ex[0];
	FArrayBox& eyfab = Ey[0];
	FArrayBox& ezfab = Ez[0];
	FArrayBox& bxfab = Bx[0];
	FArrayBox& byfab = By[0];
	FArrayBox& bzfab = Bz[0];

	std::uniform_real_distribution<Real> rand_dis2(0.0,100.0);

	for (BoxIterator bxi(exfab.box()); bxi.ok(); ++bxi)
	{
	    exfab(bxi()) = rand_dis2(rand_eng);
	}

	for (BoxIterator bxi(eyfab.box()); bxi.ok(); ++bxi)
	{
	    eyfab(bxi()) = rand_dis2(rand_eng);
	}

	for (BoxIterator bxi(ezfab.box()); bxi.ok(); ++bxi)
	{
	    ezfab(bxi()) = rand_dis2(rand_eng);
	}

	for (BoxIterator bxi(bxfab.box()); bxi.ok(); ++bxi)
	{
	    bxfab(bxi()) = rand_dis2(rand_eng);
	}

	for (BoxIterator bxi(byfab.box()); bxi.ok(); ++bxi)
	{
	    byfab(bxi()) = rand_dis2(rand_eng);
	}

	for (BoxIterator bxi(bzfab.box()); bxi.ok(); ++bxi)
	{
	    bzfab(bxi()) = rand_dis2(rand_eng);
	}

	long ngx = ng;
	long ngy = ng;
	long ngz = ng;

	const int ll4symtry          = false;
	const int l_lower_order_in_v = true;
	long lvect_fieldgathe = 64;
	warpx_geteb_energy_conserving(&np, xp.data(), yp.data(), zp.data(),
				      Exp.data(),Eyp.data(),Ezp.data(),
				      Bxp.data(),Byp.data(),Bzp.data(),
                                      ixyzmin,
                                      &xyzmin[0], &xyzmin[1], &xyzmin[2],
                                      &dx[0], &dx[1], &dx[2],
				      &nox, &noy, &noz, 
                                      BL_TO_FORTRAN_ANYD(exfab),
                                      BL_TO_FORTRAN_ANYD(eyfab),
                                      BL_TO_FORTRAN_ANYD(ezfab),
                                      BL_TO_FORTRAN_ANYD(bxfab),
                                      BL_TO_FORTRAN_ANYD(byfab),
                                      BL_TO_FORTRAN_ANYD(bzfab),
				      &ll4symtry, &l_lower_order_in_v,
				      &lvect_fieldgathe,
				      &field_gathering_algo);

	MultiFab plotmf(ba, dm, 6, 0);
	FArrayBox& plotfab = plotmf[0];
	plotfab.setVal(0.0);
	for (int k = 0; k < nz; ++k) {
	    for (int j = 0; j < ny; ++j) {
		for (int i = 0; i < nx; ++i) {
		    IntVect cell{i,j,k};
		    int ip = i + j*nx + k*nx*ny;
		    plotfab(cell,0) += Exp[ip];
		    plotfab(cell,1) += Eyp[ip];
		    plotfab(cell,2) += Ezp[ip];
		    plotfab(cell,3) += Bxp[ip];
		    plotfab(cell,4) += Byp[ip];
		    plotfab(cell,5) += Bzp[ip];
		}
	    }
	}

	RealBox realbox{domain_box, dx, xyzmin};
	int is_per[3] = {0,0,0};
	Geometry geom{domain_box, &realbox, 0, is_per};
	std::string plotname{"plt00000"};
	Vector<std::string> varnames{"Ex", "Ey", "Ez", "Bx", "By", "Bz"};
	amrex::WriteSingleLevelPlotfile(plotname, plotmf, varnames, geom, 0.0, 0);
    }

    amrex::Finalize();
}
