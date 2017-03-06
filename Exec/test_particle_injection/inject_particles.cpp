
#include <random>
#include <algorithm>

#include <WarpX.H>
#include <WarpXWrappers.h>

using namespace amrex;

namespace {
    void make_particles (Array<Real>& x, Array<Real>& y, Array<Real>& z,
			 Array<Real>& vx, Array<Real>& vy, Array<Real>& vz,
			 Array<Real>& w);
}

void inject_particles ()
{
    Array<Real> x, y, z, vx, vy, vz, w;
    make_particles(x, y, z, vx, vy, vz, w);

    int nattr = 1;
    int uniqueparticles = 1;

    addNParticles(0, x.size(), x.data(), y.data(), z.data(), vx.data(), vy.data(), vz.data(),
		  nattr, w.data(), uniqueparticles);
}

namespace
{
    void make_particles (Array<Real>& x, Array<Real>& y, Array<Real>& z,
			 Array<Real>& vx, Array<Real>& vy, Array<Real>& vz,
			 Array<Real>& w)
    {
	static std::mt19937 rand_eng(42+ParallelDescriptor::MyProc());
	static std::uniform_real_distribution<Real> rand_real(0.0,1.0);
	static std::poisson_distribution<> rand_int(0.5);
	
	auto& warpx = WarpX::GetInstance();

	static bool first = true;
	static BoxArray ba;
	static Array<int> procmap;

	const Real weight = 1.57;

	if (first)
	{
	    first = false;

	    const Geometry& geom = warpx.Geom(0);
	    const Box& domain = geom.Domain();
	    
	    ba = BoxArray{domain};
	    const int nprocs = ParallelDescriptor::NProcs();
	    // chopping grids into at least nprocs pieces
	    IntVect chunk(domain.size());
	    while (ba.size() < nprocs)
	    {
		for (int j = BL_SPACEDIM-1; j >= 0 ; j--)
		{
		    chunk[j] /= 2;
		    
		    if (ba.size() < nprocs) {
			ba.maxSize(chunk);
		    }
		}
	    }
	
	    DistributionMapping dm {ba, nprocs};
	    procmap.assign(dm.ProcessorMap().begin(), dm.ProcessorMap().begin()+ba.size());
	}

	// To make things more interesting, let's permutate it everytime;.
	std::next_permutation(procmap.begin(), procmap.end());

	const int myproc = ParallelDescriptor::MyProc();

	const Geometry& geom = warpx.Geom(0);
	const Real* problo = geom.ProbLo();
	const Real* dx = geom.CellSize();

	for (int i = 0, n = ba.size(); i < n; ++i)
	{
	    if (procmap[i] == myproc)
	    {
		const Box& bx = ba[i];
		for (IntVect cell = bx.smallEnd(); cell <= bx.bigEnd(); bx.next(cell))
		{
		    Real x0 = problo[0] + cell[0]*dx[0];
		    Real y0 = problo[1] + cell[1]*dx[1];
		    Real z0 = problo[2] + cell[2]*dx[2];
		    // number of particles to add to this cell
		    int np = rand_int(rand_eng);
		    for (int ip = 0; ip < np; ++ip)
		    {
			x.push_back(x0 + dx[0]*rand_real(rand_eng));
			y.push_back(x0 + dx[1]*rand_real(rand_eng));
			z.push_back(x0 + dx[2]*rand_real(rand_eng));
			Real vxt,vyt,vzt,v2;
			do {
			    vxt = rand_real(rand_eng);
			    vyt = rand_real(rand_eng);
			    vzt = rand_real(rand_eng);
			    v2 = vxt*vxt + vyt*vyt + vzt*vzt;
			} while(v2 >= 0.999999);
			Real gam = 1.0/sqrt(1.0-v2);
			vx.push_back(vxt*gam);
			vy.push_back(vyt*gam);
			vz.push_back(vzt*gam);
			w.push_back(weight);
		    }
		}
	    }
	}
    }
}
