#include "EMParticleContainer.H"
#include "Constants.H"

#include "em_pic_K.H"

using namespace amrex;

namespace
{    
    AMREX_GPU_HOST_DEVICE void get_position_unit_cell(Real* r, const IntVect& nppc, int i_part)
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

    AMREX_GPU_HOST_DEVICE void get_gaussian_random_momentum(Real* u, Real u_mean, Real u_std) {
        Real ux_th = amrex::RandomNormal(0.0, u_std);
        Real uy_th = amrex::RandomNormal(0.0, u_std);
        Real uz_th = amrex::RandomNormal(0.0, u_std);
        
        u[0] = u_mean + ux_th;
        u[1] = u_mean + uy_th;
        u[2] = u_mean + uz_th;
    }
}

EMParticleContainer::
EMParticleContainer(const Geometry            & a_geom,
                    const DistributionMapping & a_dmap,
                    const BoxArray            & a_ba,
                    const int                   a_species_id,
                    const Real                  a_charge,
                    const Real                  a_mass)
    : ParticleContainer<0, 0, PIdx::nattribs, 0>(a_geom, a_dmap, a_ba),
    m_species_id(a_species_id), m_charge(a_charge), m_mass(a_mass)
{}

void
EMParticleContainer::
InitParticles(const IntVect& a_num_particles_per_cell,
              const Real     a_thermal_momentum_std,
              const Real     a_thermal_momentum_mean,
              const Real     a_density,
              const RealBox& a_bounds,
              const int      a_problem)
{
    BL_PROFILE("EMParticleContainer::InitParticles");

    const int lev = 0;   
    const auto dx = Geom(lev).CellSizeArray();
    const auto plo = Geom(lev).ProbLoArray();
    
    const int num_ppc = AMREX_D_TERM( a_num_particles_per_cell[0],
                                      *a_num_particles_per_cell[1],
                                      *a_num_particles_per_cell[2]);
    const Real scale_fac = dx[0]*dx[1]*dx[2]/num_ppc;
    
    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        
        // count number of particles we will add
        int num_to_add = 0;
        Gpu::DeviceScalar<int> num_to_add_dev(num_to_add);
        int* p_num_to_add_dev = num_to_add_dev.dataPtr();
        
        amrex::ParallelFor(tile_box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for (int i_part=0; i_part<num_ppc;i_part++)
            {
                Real r[3];
                
                get_position_unit_cell(r, a_num_particles_per_cell, i_part);
                
                Real x = plo[0] + (i + r[0])*dx[0];
                Real y = plo[1] + (j + r[1])*dx[1];
                Real z = plo[2] + (k + r[2])*dx[2];
                
                if (x >= a_bounds.hi(0) || x < a_bounds.lo(0) ||
                    y >= a_bounds.hi(1) || y < a_bounds.lo(1) ||
                    z >= a_bounds.hi(2) || z < a_bounds.lo(2) ) continue;
                
                Cuda::Atomic::Add(p_num_to_add_dev, 1);
            }
        });

        num_to_add = num_to_add_dev.dataValue();

        amrex::CheckSeedArraySizeAndResize(tile_box.numPts());

        auto& particles = GetParticles(lev);
        auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];

        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + num_to_add;
        particle_tile.resize(new_size);
        
        ParticleType* pstruct = particle_tile.GetArrayOfStructs()().data();
        
        auto arrdata = particle_tile.GetStructOfArrays().realarray();
        
        int procID = ParallelDescriptor::MyProc();

        amrex::ParallelFor(tile_box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for (int i_part=0; i_part<num_ppc;i_part++)
            {
                Real r[3];
                Real u[3];
                
                get_position_unit_cell(r, a_num_particles_per_cell, i_part);
                
                Real x = plo[0] + (i + r[0])*dx[0];
                Real y = plo[1] + (j + r[1])*dx[1];
                Real z = plo[2] + (k + r[2])*dx[2];
                
                if (x >= a_bounds.hi(0) || x < a_bounds.lo(0) ||
                    y >= a_bounds.hi(1) || y < a_bounds.lo(1) ||
                    z >= a_bounds.hi(2) || z < a_bounds.lo(2) ) continue;

                if (a_problem == 0) {
                    get_gaussian_random_momentum(u, a_thermal_momentum_mean,
                                                 a_thermal_momentum_std);
                }
                else if (a_problem == 1 ) {
                    u[0] = 0.01;
                    u[1] = 0.0;
                    u[2] = 0.0;
                } else {
                    amrex::Abort("problem type not valid");
                }
                
                long pidx = num_ppc * tile_box.index(IntVect(i, j, k)) + i_part;

                ParticleType& p = pstruct[pidx];
                p.id()  = 0;
                p.cpu() = procID;
                p.pos(0) = x;
                p.pos(1) = y;
                p.pos(2) = z;
                
                arrdata[PIdx::ux  ][pidx] = u[0] * PhysConst::c;
                arrdata[PIdx::uy  ][pidx] = u[1] * PhysConst::c;
                arrdata[PIdx::uz  ][pidx] = u[2] * PhysConst::c;
                arrdata[PIdx::w   ][pidx] = a_density * scale_fac;
                arrdata[PIdx::Ex  ][pidx] = 0.0;
                arrdata[PIdx::Ey  ][pidx] = 0.0;
                arrdata[PIdx::Ez  ][pidx] = 0.0;
                arrdata[PIdx::Bx  ][pidx] = 0.0;
                arrdata[PIdx::By  ][pidx] = 0.0;
                arrdata[PIdx::Bz  ][pidx] = 0.0;
                arrdata[PIdx::ginv][pidx] = 0.0;
            }            
        });
    }
}

void EMParticleContainer::
PushAndDeposeParticles(const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                       const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                       MultiFab& jx, MultiFab& jy, MultiFab& jz, Real dt)
{
    BL_PROFILE("EMParticleContainer::PushAndDeposeParticles");

    const int lev = 0;

    const auto dxi = Geom(lev).InvCellSizeArray();
    const auto plo = Geom(lev).ProbLoArray();

    for (EMParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const int np  = pti.numParticles();

        ParticleType * pstruct = &(pti.GetArrayOfStructs()[0]);

        auto& attribs = pti.GetAttribs();
        Real *  wp  = attribs[PIdx::w].data();
        Real * uxp  = attribs[PIdx::ux].data();
        Real * uyp  = attribs[PIdx::uy].data();
        Real * uzp  = attribs[PIdx::uz].data();
/*
        Real       * AMREX_RESTRICT Exp  = attribs[PIdx::Ex].data();
        Real       * AMREX_RESTRICT Eyp  = attribs[PIdx::Ey].data();
        Real       * AMREX_RESTRICT Ezp  = attribs[PIdx::Ez].data();
        Real       * AMREX_RESTRICT Bxp  = attribs[PIdx::Bx].data();
        Real       * AMREX_RESTRICT Byp  = attribs[PIdx::By].data();
        Real       * AMREX_RESTRICT Bzp  = attribs[PIdx::Bz].data();
        Real       * AMREX_RESTRICT ginv = attribs[PIdx::ginv].data();
*/

        auto const Exarr = Ex.array(pti);
        auto const Eyarr = Ey.array(pti);
        auto const Ezarr = Ez.array(pti);
        auto const Bxarr = Bx.array(pti);
        auto const Byarr = By.array(pti);
        auto const Bzarr = Bz.array(pti);
        auto       jxarr = jx.array(pti);
        auto       jyarr = jy.array(pti);
        auto       jzarr = jz.array(pti);

        Real q = m_charge;
        Real m = m_mass;
        AMREX_FOR_1D ( np, i,
        {
            amrex::Real Exp;
            amrex::Real Eyp;
            amrex::Real Ezp;
            amrex::Real Bxp;
            amrex::Real Byp;
            amrex::Real Bzp;
            amrex::Real ginv;
            gather_fields(pstruct[i], Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                          Exarr, Eyarr, Ezarr, Bxarr, Byarr, Bzarr, plo, dxi);

            push_momentum_boris(uxp[i], uyp[i], uzp[i], ginv, Exp, Eyp, Ezp,
                                Bxp, Byp, Bzp, q, m, dt);

            push_position_boris(pstruct[i], uxp[i], uyp[i], uzp[i], ginv, dt);

            deposit_current(jxarr, jyarr, jzarr, pstruct[i], uxp[i], uyp[i], uzp[i],
                            ginv, wp[i], q, dt, plo, dxi);
        });
    }
}

void EMParticleContainer::
PushParticleMomenta(const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                    const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz, Real dt)
{   
    BL_PROFILE("EMParticleContainer::PushParticleMomenta");

    const int lev = 0;
    
    const auto dxi = Geom(lev).InvCellSizeArray();
    const auto plo = Geom(lev).ProbLoArray();

    for (EMParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const int np  = pti.numParticles();
        
        ParticleType const* AMREX_RESTRICT pstruct = &(pti.GetArrayOfStructs()[0]);

        auto& attribs = pti.GetAttribs();
        Real * AMREX_RESTRICT uxp  = attribs[PIdx::ux].data();
        Real * AMREX_RESTRICT uyp  = attribs[PIdx::uy].data();
        Real * AMREX_RESTRICT uzp  = attribs[PIdx::uz].data();
        /*
        Real * AMREX_RESTRICT Exp  = attribs[PIdx::Ex].data();
        Real * AMREX_RESTRICT Eyp  = attribs[PIdx::Ey].data();
        Real * AMREX_RESTRICT Ezp  = attribs[PIdx::Ez].data();
        Real * AMREX_RESTRICT Bxp  = attribs[PIdx::Bx].data();
        Real * AMREX_RESTRICT Byp  = attribs[PIdx::By].data();
        Real * AMREX_RESTRICT Bzp  = attribs[PIdx::Bz].data();
        Real * AMREX_RESTRICT ginv = attribs[PIdx::ginv].data();
        */

        auto const Exarr = Ex.array(pti);
        auto const Eyarr = Ey.array(pti);
        auto const Ezarr = Ez.array(pti);
        auto const Bxarr = Bx.array(pti);
        auto const Byarr = By.array(pti);
        auto const Bzarr = Bz.array(pti);

        Real q = m_charge;
        Real m = m_mass;
        AMREX_PARALLEL_FOR_1D ( np, i,
        {
            amrex::Real Exp;
            amrex::Real Eyp;
            amrex::Real Ezp;
            amrex::Real Bxp;
            amrex::Real Byp;
            amrex::Real Bzp;
            amrex::Real ginv;
            gather_fields(pstruct[i], Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                          Exarr, Eyarr, Ezarr, Bxarr, Byarr, Bzarr, plo, dxi);

            push_momentum_boris(uxp[i], uyp[i], uzp[i], ginv, Exp, Eyp, Ezp,
                                Bxp, Byp, Bzp, q, m, dt);
        });
    }
}

