
#include <BLProfiler.H>

#include <ParticleContainer.H>
#include <PICSAR_f.H>

void
MyParticleContainer::Init(MultiFab& dummy_mf)
{
    BL_PROFILE("MyPC::Init()");

    charge = -q_e;
    mass = m_e;

    m_particles.reserve(15);
    m_particles.resize(m_gdb->finestLevel()+1);

    const Geometry& geom = m_gdb->Geom(0);
    const Real* dx  = geom.CellSize();

    Real weight, ux;
    {
      ParmParse pp("langmuirwave");
      weight = 1.e25;
      pp.query("n_e", weight);
      weight *= dx[0]*dx[1]*dx[2];

      ux = 0.01;
      pp.query("ux", ux);
      ux *= 2.9961479e8;
    }
    
    for (MFIter mfi(dummy_mf,false); mfi.isValid(); ++mfi)
    {
        Box grid = m_gdb->ParticleBoxArray(0)[mfi.index()];
        RealBox grid_box = RealBox(grid,dx,geom.ProbLo());

	int nx = grid.length(0), ny = grid.length(1), nz = grid.length(2); 

	for (int k = 0; k < nz; k++) {
	  for (int j = 0; j < ny; j++) {
	    for (int i = 0; i < nx; i++) {
	      Real x = grid_box.lo(0) + (i+0.5)*dx[0];
	      if (x > 0) continue;
	      
	      ParticleType p;

	      p.m_pos[0] = grid_box.lo(0) + (0.5 + i)*dx[0];
	      p.m_pos[1] = grid_box.lo(1) + (0.5 + j)*dx[1];
	      p.m_pos[2] = grid_box.lo(2) + (0.5 + k)*dx[2];
	      
	      for (int i = 0; i < BL_SPACEDIM; i++) {
		BL_ASSERT(p.m_pos[i] < grid_box.hi(i));
	      }
	      
	      p.m_data[PIdx::w] = weight;
	      
	      for (int i = 1; i < PIdx::nattribs; i++) {
		p.m_data[i] = 0;
	      }
	      
	      p.m_data[PIdx::ux] = ux;

	      p.m_id  = ParticleBase::NextID();
	      p.m_cpu = ParallelDescriptor::MyProc();

	      if (!ParticleBase::Where(p,m_gdb))
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

    //
    // We still need to redistribute in order to define each particle's cell, grid and level, but this 
    //    shouldn't require any inter-node communication because the particles should already be in the right grid.
    //
    Redistribute(true);
}
