
#include <limits>
#include <random>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <WarpXConst.H>
#include <WarpX_f.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
	long particle_pusher_algo = 0;
	{
	    ParmParse pp("algo");
	    pp.query("particle_pusher", particle_pusher_algo);
	}

	long nx = 64, ny = 64, nz = 64;
	long np = nx*ny*nz;

	Real xyzmin[3] = {0.5, 1.4, 0.3};

	Vector<Real> xp, yp, zp, uxp, uyp, uzp, giv, Exp, Eyp, Ezp, Bxp, Byp, Bzp;
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
		    Real vx,vy,vz,v2;
		    do {
			vx = rand_dis(rand_eng);
			vy = rand_dis(rand_eng);
			vz = rand_dis(rand_eng);
			v2 = vx*vx + vy*vy + vz*vz;
		    } while(v2 >= 0.999999);
		    Real gam = 1.0/sqrt(1.0-v2);
		    uxp.push_back(vx*PhysConst::c*gam);
		    uyp.push_back(vy*PhysConst::c*gam);
		    uzp.push_back(vz*PhysConst::c*gam);
		    giv.push_back(std::numeric_limits<Real>::quiet_NaN());
		}
	    }
	}

	Exp.resize(np,0.0);
	Eyp.resize(np,0.0);
	Ezp.resize(np,0.0);
	Bxp.resize(np,0.0);
	Byp.resize(np,0.0);
	Bzp.resize(np,0.0);
	for (int i = 0; i < np; ++i) {
	    Exp[i] = rand_dis(rand_eng)*1.0e5;
	    Eyp[i] = rand_dis(rand_eng)*1.0e5;
	    Ezp[i] = rand_dis(rand_eng)*1.0e5;
	    Bxp[i] = rand_dis(rand_eng)*1.0e-5;
	    Byp[i] = rand_dis(rand_eng)*1.0e-5;
	    Bzp[i] = rand_dis(rand_eng)*1.0e-5;
	}

	Real charge = -PhysConst::q_e;
	Real mass   =  PhysConst::m_e;
	Real dt     = 1.e-10;
	warpx_particle_pusher(&np, xp.data(), yp.data(), zp.data(),
			      uxp.data(), uyp.data(), uzp.data(), giv.data(),
			      Exp.dataPtr(), Eyp.dataPtr(), Ezp.dataPtr(),
			      Bxp.dataPtr(), Byp.dataPtr(), Bzp.dataPtr(),
			      &charge, &mass, &dt, 
			      &particle_pusher_algo);

	Box plotbox{IntVect{D_DECL(0,0,0)},IntVect{D_DECL(nx-1,ny-1,nz-1)}};
	BoxArray plotba {plotbox};
	DistributionMapping plotdm {plotba};
	MultiFab plotmf(plotba, plotdm, 7, 0);
	FArrayBox& plotfab = plotmf[0];
	plotfab.setVal(0.0);
	for (int k = 0; k < nz; ++k) {
	    for (int j = 0; j < ny; ++j) {
		for (int i = 0; i < nx; ++i) {
		    IntVect cell{i,j,k};
		    int ip = i + j*nx + k*nx*ny;
		    plotfab(cell,0) += xp[ip];
		    plotfab(cell,1) += yp[ip];
		    plotfab(cell,2) += zp[ip];
		    plotfab(cell,3) += uxp[ip];
		    plotfab(cell,4) += uyp[ip];
		    plotfab(cell,5) += uzp[ip];
		    plotfab(cell,6) += 1.0/giv[ip];
		}
	    }
	}

	RealBox realbox{plotbox, dx, xyzmin};
	int is_per[3] = {0,0,0};
	Geometry geom{plotbox, &realbox, 0, is_per};
	std::string plotname{"plt00000"};
	Vector<std::string> varnames{"x", "y", "z", "ux", "uy", "uz", "gamma"};
	amrex::WriteSingleLevelPlotfile(plotname, plotmf, varnames, geom, 0.0, 0);
    }

    amrex::Finalize();
}
