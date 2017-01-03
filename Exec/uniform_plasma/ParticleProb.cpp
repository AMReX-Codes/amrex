
//
// Each problem must have its own version of MyParticleContainer::InitData()
// to initialize the particle data on this level
//

#include <cmath>

#include <AMReX_BLProfiler.H>

#include <ParticleContainer.H>
#include <WarpXConst.H>
#include <random>

void
MyParticleContainer::InitData()
{
    BL_PROFILE("MyPC::InitData()");

    charge = -PhysConst::q_e;
    mass = PhysConst::m_e;

    m_particles.resize(m_gdb->finestLevel()+1);

    const int lev = 0;

    const Geometry& geom = m_gdb->Geom(lev);
    const Real* dx  = geom.CellSize();

    Real weight, ux, uy, uz, u_th;
    Real particle_shift_x, particle_shift_y, particle_shift_z;
    int n_part_per_cell;
    {
      ParmParse pp("uniform_plasma");
      n_part_per_cell = 1;
      pp.query("num_particles_per_cell", n_part_per_cell);
      weight = 1.e25;
      pp.query("n_e", weight);
      #if BL_SPACEDIM==3
      weight *= dx[0]*dx[1]*dx[2]/n_part_per_cell;
      #elif BL_SPACEDIM==2
      weight *= dx[0]*dx[1]/n_part_per_cell;
      #endif

      u_th = 0.;
      pp.query("u_th", u_th);
      u_th *= PhysConst::c;
    }

    const BoxArray& ba = m_gdb->ParticleBoxArray(lev);
    const DistributionMapping& dm = m_gdb->ParticleDistributionMap(lev);

    MultiFab dummy_mf(ba, 1, 0, dm, Fab_noallocate);

    for (MFIter mfi(dummy_mf,false); mfi.isValid(); ++mfi)
    {
	int gid = mfi.index();
        Box grid = ba[gid];
        RealBox grid_box { grid,dx,geom.ProbLo() };

#if (BL_SPACEDIM == 3)
	int nx = grid.length(0), ny = grid.length(1), nz = grid.length(2);
#elif (BL_SPACEDIM == 2)
	int nx = grid.length(0), ny = 1, nz = grid.length(1);
#endif

  // Initialize random generator for normal distribution
  std::default_random_engine generator;
  std::normal_distribution<Real> velocity_distribution(0.0, u_th);
  std::uniform_real_distribution<Real> position_distribution(0.0,1.0);

	for (int k = 0; k < nz; k++) {
	  for (int j = 0; j < ny; j++) {
	    for (int i = 0; i < nx; i++) {
	      for (int i_part=0; i_part<n_part_per_cell;i_part++) {

          // Randomly generate the speed according to a normal distribution
          ux = velocity_distribution(generator);
          uy = velocity_distribution(generator);
          uz = velocity_distribution(generator);

          // Randomly generate the positions (uniformly inside each cell)
          particle_shift_x = position_distribution(generator);
          particle_shift_y = position_distribution(generator);
          particle_shift_z = position_distribution(generator);
#if (BL_SPACEDIM == 3)
		Real x = grid_box.lo(0) + (i + particle_shift_x)*dx[0];
		Real y = grid_box.lo(1) + (j + particle_shift_y)*dx[1];
		Real z = grid_box.lo(2) + (k + particle_shift_z)*dx[2];
#elif (BL_SPACEDIM == 2)
		Real x = grid_box.lo(0) + (i + particle_shift_x)*dx[0];
		Real y = 0.0;
		Real z = grid_box.lo(1) + (k + particle_shift_z)*dx[1];
#endif

		ParticleType p;

		p.m_id  = ParticleBase::NextID();
		p.m_cpu = ParallelDescriptor::MyProc();
		p.m_lev = lev;
		p.m_grid = gid;

#if (BL_SPACEDIM == 3)
		p.m_pos[0] = x;
		p.m_pos[1] = y;
		p.m_pos[2] = z;
#elif (BL_SPACEDIM == 2)
		p.m_pos[0] = x;
		p.m_pos[1] = z;
#endif

		for (int i = 0; i < BL_SPACEDIM; i++) {
		  BL_ASSERT(p.m_pos[i] < grid_box.hi(i));
		}

		p.m_data[PIdx::w] = weight;

		for (int i = 1; i < PIdx::nattribs; i++) {
		  p.m_data[i] = 0;
		}

		p.m_data[PIdx::ux] = ux;
		p.m_data[PIdx::uy] = uy;
		p.m_data[PIdx::uz] = uz;

		if (!ParticleBase::Where(p,m_gdb)) // this will set m_cell
		{
		    BoxLib::Abort("invalid particle");
		}

		BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_gdb->finestLevel());
		//
		// Add it to the appropriate PBox at the appropriate level.
		//
		m_particles[p.m_lev][p.m_grid].push_back(p);
	      }
	    }
	  }
        }
    }

    //
    // We still need to redistribute in order to define each particle's cell, grid and level, but this
    //    shouldn't require any inter-node communication because the particles should already be in the right grid.
    //
    Redistribute(true);
}
