
//
// Each problem must have its own version of PhysicalParticleContainer::InitData()
// to initialize the particle data.  It must also initialize charge and mass.
//

#include <cmath>

#include <AMReX_BLProfiler.H>

#include <ParticleContainer.H>
#include <WarpXConst.H>
#include <random>

using namespace amrex;

void
PhysicalParticleContainer::InitData()
{
    BL_PROFILE("PPC::InitData()");

    charge = -PhysConst::q_e;
    mass = PhysConst::m_e;

    const int lev = 0;

    const Geometry& geom = Geom(lev);
    const Real* dx  = geom.CellSize();

    Real weight, u_th;
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

    std::array<Real,PIdx::nattribs> attribs;
    attribs.fill(0.0);
    attribs[PIdx::w ] = weight;

    // Initialize random generator for normal distribution
    std::default_random_engine generator;
    std::normal_distribution<Real> velocity_distribution(0.0, u_th);
    std::uniform_real_distribution<Real> position_distribution(0.0,1.0);

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
                // Randomly generate the speed according to a normal distribution
                Real ux = velocity_distribution(generator);
                Real uy = velocity_distribution(generator);
                Real uz = velocity_distribution(generator);

                // Randomly generate the positions (uniformly inside each cell)
                Real particle_shift_x = position_distribution(generator);
                Real particle_shift_y = position_distribution(generator);
                Real particle_shift_z = position_distribution(generator);

#if (BL_SPACEDIM == 3)
                Real x = tile_real_box.lo(0) + (iv[0]-boxlo[0] + particle_shift_x)*dx[0];
                Real y = tile_real_box.lo(1) + (iv[1]-boxlo[1] + particle_shift_y)*dx[1];
                Real z = tile_real_box.lo(2) + (iv[2]-boxlo[2] + particle_shift_z)*dx[2];
#elif (BL_SPACEDIM == 2)
                Real x = tile_real_box.lo(0) + (iv[0]-boxlo[0] + particle_shift_x)*dx[0];
                Real y = 0.0;
                Real z = tile_real_box.lo(1) + (iv[1]-boxlo[1] + particle_shift_z)*dx[1];
#endif   

                attribs[PIdx::ux] = ux;
                attribs[PIdx::uy] = uy;
                attribs[PIdx::uz] = uz;

                AddOneParticle(lev, grid_id, tile_id, x, y, z, attribs);
            }
        }
    }
}
