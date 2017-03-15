
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
    BL_PROFILE("PhysicalParticleContainer::InitData()");

    charge = plasma_injector->getCharge();
    mass = plasma_injector->getMass();

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real* dx  = geom.CellSize();

    Real scale_fac;
    int n_part_per_cell = plasma_injector->numParticlesPerCell();

#if BL_SPACEDIM==3
    scale_fac = dx[0]*dx[1]*dx[2]/n_part_per_cell;
#elif BL_SPACEDIM==2
    scale_face = dx[0]*dx[1]/n_part_per_cell;
#endif

    std::array<Real,PIdx::nattribs> attribs;
    attribs.fill(0.0);
    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
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
                
                if (plasma_injector->insideBounds(x, y, z)) {
                    Real weight;
                    Real u[3];
                    weight = plasma_injector->getDensity(x, y, z) * scale_fac;
                    plasma_injector->getMomentum(u);
                    attribs[PIdx::w ] = weight;
                    attribs[PIdx::ux] = u[0];
                    attribs[PIdx::uy] = u[1];
                    attribs[PIdx::uz] = u[2];
                    AddOneParticle(lev, grid_id, tile_id, x, y, z, attribs);
                }
            }
        }
    }
}
