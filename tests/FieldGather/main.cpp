
#include <random>

#include <BoxLib.H>
#include <ParmParse.H>
#include <Array.H>
#include <MultiFab.H>

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

	long field_gathering_algo = 1;

	{
	    ParmParse pp("algo");
	    pp.query("field_gathering", field_gathering_algo);
	}

	long nx = 64, ny = 64, nz = 64;
	long np = nx*ny*nz;

	Real xyzmin[3] = {0.5, 1.4, 0.3};

	Array<Real> xp, yp, zp, Exp, Eyp, Ezp, Bxp, Byp, Bzp;
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
	Box domain_box {IntVect{D_DECL(0,0,0)}, IntVect{D_DECL(nx,ny,nz)}};
	Box grown_box = BoxLib::grow(domain_box, ng);

	long ngx = ng;
	long ngy = ng;
	long ngz = ng;

	FArrayBox exfab(grown_box);
	FArrayBox eyfab(grown_box);
	FArrayBox ezfab(grown_box);
	FArrayBox bxfab(grown_box);
	FArrayBox byfab(grown_box);
	FArrayBox bzfab(grown_box);

	for (IntVect cell=grown_box.smallEnd(); cell <= grown_box.bigEnd(); grown_box.next(cell))
	{
	    exfab(cell) = rand_dis(rand_eng);
	    eyfab(cell) = rand_dis(rand_eng);
	    ezfab(cell) = rand_dis(rand_eng);
	    bxfab(cell) = rand_dis(rand_eng);
	    byfab(cell) = rand_dis(rand_eng);
	    bzfab(cell) = rand_dis(rand_eng);
	}

	const int ll4symtry          = false;
	const int l_lower_order_in_v = true;
	warpx_geteb_energy_conserving(&np, xp.data(), yp.data(), zp.data(),
				      Exp.data(),Eyp.data(),Ezp.data(),
				      Bxp.data(),Byp.data(),Bzp.data(),
				      &xyzmin[0], &xyzmin[1], &xyzmin[2],
				      &dx[0], &dx[1], &dx[2],
				      &nx, &ny, &nz, &ngx, &ngy, &ngz, 
				      &nox, &noy, &noz, 
				      exfab.dataPtr(), eyfab.dataPtr(), ezfab.dataPtr(),
				      bxfab.dataPtr(), byfab.dataPtr(), bzfab.dataPtr(),
				      &ll4symtry, &l_lower_order_in_v,
				      &field_gathering_algo);
    }

    BoxLib::Finalize();
}
