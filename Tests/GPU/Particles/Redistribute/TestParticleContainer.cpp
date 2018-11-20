#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>

#include "TestParticleContainer.H"
#include "Constants.H"

#include "test_F.H"

using namespace amrex;

namespace
{    
    void get_position_unit_cell(Real* r, const IntVect& nppc, int i_part)
    {
        int nx = nppc[0];
        int ny = nppc[1];
        int nz = nppc[2];
        
        int ix_part = i_part/(ny * nz);
        int iy_part = (i_part % (ny * nz)) % ny;
        int iz_part = (i_part % (ny * nz)) / ny;
        
        r[0] = (0.5+ix_part)/nx;
        r[1] = (0.5+iy_part)/ny;
        r[2] = (0.5+iz_part)/nz;
    }
    
    void get_gaussian_random_momentum(Real* u, Real u_mean, Real u_std) {
        Real ux_th = amrex::RandomNormal(0.0, u_std);
        Real uy_th = amrex::RandomNormal(0.0, u_std);
        Real uz_th = amrex::RandomNormal(0.0, u_std);
        
        u[0] = u_mean + ux_th;
        u[1] = u_mean + uy_th;
        u[2] = u_mean + uz_th;
    }
}

TestParticleContainer::
TestParticleContainer(const Geometry            & a_geom,
                      const DistributionMapping & a_dmap,
                      const BoxArray            & a_ba)
    : amrex::ParticleContainer<0, 0, PIdx::nattribs>(a_geom, a_dmap, a_ba)
{}

void
TestParticleContainer::
InitParticles(const IntVect& a_num_particles_per_cell,
              const Real     a_thermal_momentum_std,
              const Real     a_thermal_momentum_mean,
              const Real     a_density,
              const RealBox& a_bounds)
{
    BL_PROFILE("TestParticleContainer::InitParticles");
    BL_ASSERT(finestLevel() == 0);

    const int lev = 0;
    const Geometry& geom = Geom(lev);

    const Real* dx = geom.CellSize();

    const int num_ppc = AMREX_D_TERM( a_num_particles_per_cell[0],
                                     *a_num_particles_per_cell[1],
                                     *a_num_particles_per_cell[2]);
    const Real scale_fac = dx[0]*dx[1]*dx[2]/num_ppc;

    std::array<Real, PIdx::nattribs> attribs;
    attribs.fill(0.0);

    BuildRedistributeMask(0, 1);

    for(MFIter mfi(*redistribute_mask_ptr); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const Real* plo = geom.ProbLo();
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particles = GetParticles(lev)[std::make_pair(grid_id, tile_id)];
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
            for (int i_part=0; i_part<num_ppc;i_part++) {
                Real r[3];
                Real u[3];

                get_position_unit_cell(r, a_num_particles_per_cell, i_part);

                get_gaussian_random_momentum(u, a_thermal_momentum_mean, a_thermal_momentum_std);

                Real x = plo[0] + (iv[0] + r[0])*dx[0];
                Real y = plo[1] + (iv[1] + r[1])*dx[1];
                Real z = plo[2] + (iv[2] + r[2])*dx[2];
                
                if (x >= a_bounds.hi(0) || x < a_bounds.lo(0) ||
                    y >= a_bounds.hi(1) || y < a_bounds.lo(1) ||
                    z >= a_bounds.hi(2) || z < a_bounds.lo(2) ) continue;
                
                ParticleType p;
                p.pos(0) = x;
                p.pos(1) = y;
                p.pos(2) = z;
                
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                
                attribs[PIdx::vx] = u[0] * PhysConst::c;
                attribs[PIdx::vy] = u[1] * PhysConst::c;
                attribs[PIdx::vz] = u[2] * PhysConst::c;
                attribs[PIdx::w ] = a_density * scale_fac;
                
                particles.push_back(p);
                particles.push_back_real(attribs);
            }
        }
    }
}

void TestParticleContainer::
MoveParticles()
{
    BL_PROFILE("TestParticleContainer::MoveParticles");

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();

    for (TestParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        auto& ptile       = GetParticles(lev)[std::make_pair(grid_id, tile_id)];
        auto& particles   = ptile.GetArrayOfStructs();
        const int np      = ptile.numParticles();

        if (np == 0) continue;

        auto& vx = ptile.GetStructOfArrays().GetRealData(PIdx::vx);
        auto& vy = ptile.GetStructOfArrays().GetRealData(PIdx::vy);
        auto& vz = ptile.GetStructOfArrays().GetRealData(PIdx::vz);

	FTOC(move_particles)(np,
                             particles.data(),
                             vx.data(), vy.data(), vz.data(), dx);
    }
}

