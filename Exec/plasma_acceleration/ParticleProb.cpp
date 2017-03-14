
//
// Each problem must have its own version of PhysicalParticleContainer::InitData()
// to initialize the particle data.  It must also initialize charge and mass.
//

#include <cmath>

#include <AMReX_BLProfiler.H>

#include <ParticleContainer.H>
#include <WarpXConst.H>

using namespace amrex;

void
PhysicalParticleContainer::InitData()
{
    BL_PROFILE("SPC::InitData()");

    // species_id 0 : the beam
    // species_id 1 : the electrons of the plasma
    // Note: the ions of the plasma are implicitly motionless, and so are not part of the simulation
    if (species_id == 0 or species_id == 1) {
	charge = -PhysConst::q_e;
	mass = PhysConst::m_e;
    } else {
	amrex::Abort("PhysicalParticleContainer::InitData(): species_id must be 0 or 1");
    }

    const int lev = 0;

    const Geometry& geom = Geom(lev);
    const Real* dx  = geom.CellSize();

    Real weight, gamma, uz;
    Real particle_xmin, particle_xmax, particle_ymin, particle_ymax, particle_zmin, particle_zmax;
    int n_part_per_cell;
    {
      ParmParse pp("pwfa");

      // Calculate the particle weight
      n_part_per_cell = 1;
      pp.query("num_particles_per_cell", n_part_per_cell);
      weight = 1.e22;
      if (species_id == 0) {
	pp.query("beam_density", weight);
      } else if (species_id == 1) {
	pp.query("plasma_density", weight);
      }
      #if BL_SPACEDIM==3
      weight *= dx[0]*dx[1]*dx[2]/n_part_per_cell;
      #elif BL_SPACEDIM==2
      weight *= dx[0]*dx[1]/n_part_per_cell;
      #endif

      // Calculate the limits between which the particles will be initialized
      particle_xmin = particle_ymin = particle_zmin = -2.e-5;
      particle_xmax = particle_ymax = particle_zmax =  2.e-5;
      if (species_id == 0) {
	pp.query("beam_xmin", particle_xmin);
	pp.query("beam_xmax", particle_xmax);
	pp.query("beam_ymin", particle_ymin);
	pp.query("beam_ymax", particle_ymax);
	pp.query("beam_zmin", particle_zmin);
	pp.query("beam_zmax", particle_zmax);
      } else if (species_id == 1) {
	pp.query("plasma_xmin", particle_xmin);
	pp.query("plasma_xmax", particle_xmax);
	pp.query("plasma_ymin", particle_ymin);
	pp.query("plasma_ymax", particle_ymax);
	pp.query("plasma_zmin", particle_zmin);
	pp.query("plasma_zmax", particle_zmax);
      }
      uz = 0.;
      if (species_id == 0) { // relativistic beam
	gamma = 1.;
	pp.query("beam_gamma", gamma);
	uz = std::sqrt( gamma*gamma - 1 );
      }     
      uz *= PhysConst::c;
    }

    std::array<Real,PIdx::nattribs> attribs;
    attribs.fill(0.0);
    attribs[PIdx::w ] = weight;
    attribs[PIdx::uz] = uz;

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box = mfi.tilebox();
        RealBox tile_real_box { tile_box, dx, geom.ProbLo() };

        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();

        const auto& boxlo = tile_box.smallEnd();
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
        {
            for (int i_part=0; i_part<n_part_per_cell;i_part++)
            {
                Real particle_shift = (0.5+i_part)/n_part_per_cell;
#if (BL_SPACEDIM == 3)
                Real x = tile_real_box.lo(0) + (iv[0]-boxlo[0] + particle_shift)*dx[0];
                Real y = tile_real_box.lo(1) + (iv[1]-boxlo[1] + particle_shift)*dx[1];
                Real z = tile_real_box.lo(2) + (iv[2]-boxlo[2] + particle_shift)*dx[2];
#elif (BL_SPACEDIM == 2)
                Real x = tile_real_box.lo(0) + (iv[0]-boxlo[0] + particle_shift)*dx[0];
                Real y = 0.0;
                Real z = tile_real_box.lo(1) + (iv[1]-boxlo[1] + particle_shift)*dx[1];
#endif   
                
                if (x >= particle_xmax || x < particle_xmin ||
                    y >= particle_ymax || y < particle_ymin ||
                    z >= particle_zmax || z < particle_zmin ) 
                {
                    continue;
                }
                else
                {
                    AddOneParticle(lev, grid_id, tile_id, x, y, z, attribs);
                }
            } 
        }
    }
}
