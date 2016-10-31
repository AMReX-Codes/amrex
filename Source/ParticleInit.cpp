
#include <cmath>

#include <BLProfiler.H>

#include <ParticleContainer.H>
#include <PICSAR_f.H>
#include <WarpXConst.H>

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

    Real weight, ux;
    Real particle_xmin, particle_xmax, particle_ymin, particle_ymax, particle_zmin, particle_zmax;
    int n_part_per_cell;
    {
      ParmParse pp("langmuirwave");
      n_part_per_cell = 1;
      pp.query("num_particles_per_cell", n_part_per_cell);
      weight = 1.e25;
      pp.query("n_e", weight);
      weight *= dx[0]*dx[1]*dx[2]/n_part_per_cell;

      pp.query("particle_xmin", particle_xmin);
      pp.query("particle_xmax", particle_xmax);
      pp.query("particle_ymin", particle_ymin);
      pp.query("particle_ymax", particle_ymax);
      pp.query("particle_zmin", particle_zmin);
      pp.query("particle_zmax", particle_zmax);   
 
      ux = 0.01;
      pp.query("ux", ux);

      Real gamma = 1.0 / std::sqrt(1.0 - ux*ux);
      ux *= PhysConst::c*gamma;
    }

    const BoxArray& ba = m_gdb->ParticleBoxArray(lev);
    const DistributionMapping& dm = m_gdb->ParticleDistributionMap(lev);

    MultiFab dummy_mf(ba, 1, 0, dm, Fab_noallocate);

    for (MFIter mfi(dummy_mf,false); mfi.isValid(); ++mfi)
    {
	int gid = mfi.index();
        Box grid = ba[gid];
        RealBox grid_box = RealBox(grid,dx,geom.ProbLo());

	int nx = grid.length(0), ny = grid.length(1), nz = grid.length(2); 

	for (int k = 0; k < nz; k++) {
	  for (int j = 0; j < ny; j++) {
	    for (int i = 0; i < nx; i++) {
	      for (int i_part=0; i_part<n_part_per_cell;i_part++) {
		Real particle_shift = (0.5+i_part)/n_part_per_cell;
		Real x = grid_box.lo(0) + (i + particle_shift)*dx[0];
		Real y = grid_box.lo(1) + (j + particle_shift)*dx[1];
		Real z = grid_box.lo(2) + (k + particle_shift)*dx[2];
   
		if (x >= particle_xmax || x < particle_xmin ||
		    y >= particle_ymax || y < particle_ymin ||
		    z >= particle_zmax || z < particle_zmin ) continue;
	      
		ParticleType p;

		p.m_id  = ParticleBase::NextID();
		p.m_cpu = ParallelDescriptor::MyProc();
		p.m_lev = lev;
		p.m_grid = gid; 

		p.m_pos[0] = x;
		p.m_pos[1] = y;
		p.m_pos[2] = z;
		
		for (int i = 0; i < BL_SPACEDIM; i++) {
		  BL_ASSERT(p.m_pos[i] < grid_box.hi(i));
		}
		
		p.m_data[PIdx::w] = weight;
	      
		for (int i = 1; i < PIdx::nattribs; i++) {
		  p.m_data[i] = 0;
		}
	      
		p.m_data[PIdx::ux] = ux;

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
