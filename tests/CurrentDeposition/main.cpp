
#include <random>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

#include <WarpX_f.H>
#include <WarpXConst.H>

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

	Real charge = -PhysConst::q_e;
	Real weight = 10.0;
	Real dt = 1.e-10;

	long current_deposition_algo = 3;
	{
	    ParmParse pp("algo");
	    pp.query("current_deposition", current_deposition_algo);
	}

	long nx = 64, ny = 64, nz = 64;
	long np = nx*ny*nz;

	Real xyzmin[3] = {0.5, 1.4, 0.3};

	Vector<Real> xp, yp, zp, uxp, uyp, uzp, giv, wp;
	Real dx[3] = {1.0/nx, 1.0/ny, 1.0/nz};

	std::mt19937 rand_eng(42);
	std::uniform_real_distribution<Real> rand_dis(0.0,1.0);

	for (int k = 0; k < nz; ++k) {
	    for (int j = 0; j < ny; ++j) {
		for (int i = 0; i < nx; ++i) {
		    Real x0 = xyzmin[0] + i*dx[0];
		    Real y0 = xyzmin[1] + j*dx[1];
		    Real z0 = xyzmin[2] + k*dx[2];
		    xp.push_back(x0 + dx[0]*rand_dis(rand_eng));
		    yp.push_back(y0 + dx[1]*rand_dis(rand_eng));
		    zp.push_back(z0 + dx[2]*rand_dis(rand_eng));
		    wp.push_back(weight);
		    Real vx,vy,vz,v2;
		    do {
			vx = rand_dis(rand_eng);
			vy = rand_dis(rand_eng);
			vz = rand_dis(rand_eng);
			v2 = vx*vx + vy*vy + vz*vz;
		    } while(v2 >= 0.999999);
		    Real gam = 1.0/sqrt(1.0-v2);
		    uxp.push_back(vx*gam);
		    uyp.push_back(vy*gam);
		    uzp.push_back(vz*gam);
		    giv.push_back(1.0/gam);
		}
	    }
	}

	const int ng = nox;
	Box domain_box {IntVect{D_DECL(0,0,0)}, IntVect{D_DECL(nx-1,ny-1,nz-1)}};

#if (BL_SPACEDIM == 3)
        IntVect jx_nodal_flag(0,1,1);
        IntVect jy_nodal_flag(1,0,1);
        IntVect jz_nodal_flag(1,1,0);
#elif (BL_SPACEDIM == 2)
        IntVect jx_nodal_flag(0,1);  // x is the first dimension to AMReX
        IntVect jy_nodal_flag(1,1);  // y is the missing dimension to 2D AMReX
        IntVect jz_nodal_flag(1,0);  // z is the second dimension to 2D AMReX
#endif

        BoxArray ba{domain_box};
        DistributionMapping dm{ba};

        MultiFab jx(amrex::convert(ba,jx_nodal_flag), dm, 1, ng);
        MultiFab jy(amrex::convert(ba,jy_nodal_flag), dm, 1, ng);
        MultiFab jz(amrex::convert(ba,jz_nodal_flag), dm, 1, ng);

	FArrayBox& jxfab = jx[0];
	FArrayBox& jyfab = jy[0];
	FArrayBox& jzfab = jz[0];
	jxfab.setVal(0.0);
	jyfab.setVal(0.0);
	jzfab.setVal(0.0);

	long ngx = ng;
	long ngy = ng;
	long ngz = ng;
	long lvect = 8;
	warpx_current_deposition(jxfab.dataPtr(), &ngx, jxfab.length(),
                                 jyfab.dataPtr(), &ngy, jyfab.length(),
                                 jzfab.dataPtr(), &ngz, jzfab.length(),
				 &np, xp.data(), yp.data(), zp.data(), 
				 uxp.data(), uyp.data(), uzp.data(), 
				 giv.data(), wp.data(), &charge, 
				 &xyzmin[0], &xyzmin[1], &xyzmin[2], 
				 &dt, &dx[0], &dx[1], &dx[2],
				 &nox, &noy,&noz,
				 &lvect, &current_deposition_algo);

	MultiFab plotmf(ba, dm, 3, 0);
        amrex::average_edge_to_cellcenter(plotmf, 0, {&jx, &jy, &jz});

	RealBox realbox{domain_box, dx, xyzmin};
	int is_per[3] = {0,0,0};
	Geometry geom{domain_box, &realbox, 0, is_per};
	std::string plotname{"plt00000"};
	Vector<std::string> varnames{"jx", "jy", "jz"};
	amrex::WriteSingleLevelPlotfile(plotname, plotmf, varnames, geom, 0.0, 0);
    }

    amrex::Finalize();
}
