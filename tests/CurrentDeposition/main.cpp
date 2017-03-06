
#include <random>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
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

	Array<Real> xp, yp, zp, uxp, uyp, uzp, giv, wp;
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
		    zp.push_back(z0 + dx[3]*rand_dis(rand_eng));
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
	Box domain_box {IntVect{D_DECL(0,0,0)}, IntVect{D_DECL(nx,ny,nz)}};
	Box grown_box = amrex::grow(domain_box, ng);

	long ngx = ng;
	long ngy = ng;
	long ngz = ng;

	FArrayBox jxfab(grown_box);
	FArrayBox jyfab(grown_box);
	FArrayBox jzfab(grown_box);
	jxfab.setVal(0.0);
	jyfab.setVal(0.0);
	jzfab.setVal(0.0);

	long lvect = 8;
	warpx_current_deposition(jxfab.dataPtr(), jyfab.dataPtr(), jzfab.dataPtr(),
				 &np, xp.data(), yp.data(), zp.data(), 
				 uxp.data(), uyp.data(), uzp.data(), 
				 giv.data(), wp.data(), &charge, 
				 &xyzmin[0], &xyzmin[1], &xyzmin[2], 
				 &dt, &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
				 &ngx, &ngy, &ngz, 
				 &nox, &noy,&noz,
				 &lvect, &current_deposition_algo);

	Box plotbox{IntVect{D_DECL(0,0,0)},IntVect{D_DECL(nx-1,ny-1,nz-1)}};
	BoxArray plotba {plotbox};
	DistributionMapping plotdm {plotba};
	MultiFab plotmf(plotba, plotdm, 3, 0);
	plotmf[0].copy(jxfab,0,0,1);
	plotmf[0].copy(jyfab,0,1,1);
	plotmf[0].copy(jzfab,0,2,1);

	RealBox realbox{plotbox, dx, xyzmin};
	int is_per[3] = {0,0,0};
	Geometry geom{plotbox, &realbox, 0, is_per};
	std::string plotname{"plt00000"};
	Array<std::string> varnames{"jx", "jy", "jz"};
	amrex::WriteSingleLevelPlotfile(plotname, plotmf, varnames, geom, 0.0, 0);
    }

    amrex::Finalize();
}
