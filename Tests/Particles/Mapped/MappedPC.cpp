#include "MappedPC.H"

#include <AMReX_TracerParticle_mod_K.H>

using namespace amrex;

static constexpr int NSR = 2*AMREX_SPACEDIM;
static constexpr int NSI = AMREX_SPACEDIM;
static constexpr int NAR = 0;
static constexpr int NAI = 0;


void
MappedPC::
InitParticles (MultiFab& a_xyz_loc)
{
    BL_PROFILE("MappedPC::InitParticles");

    const int lev = 0;

    for(MFIter mfi(a_xyz_loc); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = enclosedCells(mfi.tilebox());

        Gpu::HostVector<ParticleType> host_particles;
        std::array<Gpu::HostVector<ParticleReal>, NAR> host_real;
        std::array<Gpu::HostVector<int>, NAI> host_int;

        auto loc_arr = a_xyz_loc.array(mfi);

        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
        {
            int i = iv[0]; int j = iv[1];
#if (AMREX_SPACEDIM == 2)
            int k = 0;
#elif (AMREX_SPACEDIM == 3)
            int k = iv[2];
#endif
            if (iv[0] == 0) {

                // This is the physical location of the center of the cell
                Real x = 0.25*( loc_arr(i  ,j,k,0) + loc_arr(i  ,j+1,k,0)
                               +loc_arr(i+1,j,k,0) + loc_arr(i+1,j+1,k,0));
                Real y = 0.25*( loc_arr(i  ,j,k,1) + loc_arr(i  ,j+1,k,1)
                               +loc_arr(i+1,j,k,1) + loc_arr(i+1,j+1,k,1));
                amrex::Print() << "Particle at (x,y) OF " << iv << " " << x << " " << y << std::endl;

                ParticleType p;
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                p.pos(0) = x;
                p.pos(1) = y;

                p.rdata(MappedRealIdx::vx) =Real(0.0);
                p.rdata(MappedRealIdx::vy) = Real(0.0);

                p.idata(MappedIntIdx::i) = iv[0];  // particles carry their i-index
                p.idata(MappedIntIdx::j) = iv[1];  // particles carry their j-index

#if (AMREX_SPACEDIM == 3)
                Real z = Real(0.25)*(loc_arr(i  ,j,k,2) + loc_arr(i  ,j+1,k,2) +
                                     loc_arr(i+1,j,k,2) + loc_arr(i+1,j+1,k,2));
                p.pos(2) = z;
                p.rdata(MappedRealIdx::vz) = Real(0.);
                p.idata(MappedIntIdx::k) = iv[2];  // particles carry their k-index
#endif

                host_particles.push_back(p);
                for (int nr = 0; nr < NAR; ++nr)
                    host_real[nr].push_back(p.rdata(nr));
                for (int ni = 0; ni < NAI; ++ni)
                    host_int[ni].push_back(p.idata(ni));

           }
        }

            auto& particle_tile = DefineAndReturnParticleTile(lev, mfi.index(), mfi.LocalTileIndex());
            auto old_size = particle_tile.GetArrayOfStructs().size();
            auto new_size = old_size + host_particles.size();
            particle_tile.resize(new_size);

            Gpu::copyAsync(Gpu::hostToDevice,
                           host_particles.begin(),
                           host_particles.end(),
                           particle_tile.GetArrayOfStructs().begin() + old_size);

            auto& soa = particle_tile.GetStructOfArrays();
            for (int i = 0; i < NAR; ++i)
            {
                Gpu::copyAsync(Gpu::hostToDevice,
                               host_real[i].begin(),
                               host_real[i].end(),
                               soa.GetRealData(i).begin() + old_size);
            }

            for (int i = 0; i < NAI; ++i)
            {
                Gpu::copyAsync(Gpu::hostToDevice,
                               host_int[i].begin(),
                               host_int[i].end(),
                               soa.GetIntData(i).begin() + old_size);
            }

            Gpu::streamSynchronize();
    }
    RedistributeLocal();
}

/*
  /brief Uses midpoint method to advance particles using cell-centered velocity.
*/
void
MappedPC::AdvectWithUCC (MultiFab& vel_cc, int lev, Real dt, const MultiFab& a_xyz_loc)
{
    BL_PROFILE("MappedPC::AdvectWithCC()");
    AMREX_ASSERT(lev >= 0 && lev < GetParticles().size());

    auto probhi = this->ParticleContainerBase::Geom(0).ProbHi();
    auto problo = this->ParticleContainerBase::Geom(0).ProbLo();

    // Center of the annulus
    Real cx = 0.5 * (problo[0]+probhi[0]);
    Real cy = 0.5 * (problo[1]+probhi[1]);

    for (int ipass = 0; ipass < 2; ipass++)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& ptile = ParticlesAt(lev, pti);
            auto& aos  = ptile.GetArrayOfStructs();
            const int n = aos.numParticles();
            auto *p_pbox = aos().data();

            const auto loc_arr = a_xyz_loc.array(pti);
            const auto vel_cc_arr = vel_cc.const_array(pti);

            ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = p_pbox[i];

                if (p.id() <= 0) { return; }

                ParticleReal v[AMREX_SPACEDIM];
                cic_interpolate_mapped(p, vel_cc_arr, loc_arr, v);

                if (ipass == 0)
                {
                    Real r = std::sqrt((p.pos(0)-cx)*(p.pos(0)-cx)  +(p.pos(1)-cy)*(p.pos(1)-cy));
                    Real theta = atan((p.pos(1)-cy)/(p.pos(0)-cx));

                    amrex::Print() << "UPDATING FROM " << p.pos(0) << " " << p.pos(1) <<
                                       " WITH RADIUS   " << r <<
                                       " WITH  THETA   " << theta*180/3.1415926 <<
                                       " AND VEL " << v[0] << " " << v[1] << std::endl;

                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.rdata(dim) = p.pos(dim);
                        p.pos(dim) += static_cast<ParticleReal>(ParticleReal(0.5)*dt*v[dim]);
                    }
                    r = std::sqrt((p.pos(0)-cx)*(p.pos(0)-cx)  +(p.pos(1)-cy)*(p.pos(1)-cy));
                    amrex::Print() << "          TO  " << p.pos(0) << " " << p.pos(1) <<
                                       " WITH RADIUS   " << r <<
                                       " AND VEL " << v[0] << " " << v[1] << std::endl;
                }
                else
                {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.pos(dim) = p.rdata(dim) + static_cast<ParticleReal>(dt*v[dim]);
                        p.rdata(dim) = v[dim];
                    }

                    // also update z-coordinate here
                    IntVect iv(
                       AMREX_D_DECL(p.idata(0),
                                    p.idata(1),
                                    p.idata(2)));

                    auto xlo = loc_arr(AMREX_D_DECL(iv[0]  , iv[1]  , iv[2]),0);
                    auto xhi = loc_arr(AMREX_D_DECL(iv[0]+1, iv[1]  , iv[2]),0);
                    auto ylo = loc_arr(AMREX_D_DECL(iv[0]  , iv[1]  , iv[2]),1);
                    auto yhi = loc_arr(AMREX_D_DECL(iv[0]  , iv[1]+1, iv[2]),1);

                    if (p.pos(0) > xhi) { // need to be careful here
                        p.idata(0) += 1;
                    } else if (p.pos(0) <= xlo) {
                        p.idata(0) -= 1;
                    }
                    if (p.pos(1) > yhi) { // need to be careful here
                        p.idata(1) += 1;
                    } else if (p.pos(1) <= ylo) {
                        p.idata(1) -= 1;
                    }

#if (AMREX_SPACEDIM == 3)
                    auto zlo = loc_arr(iv[0], iv[1], iv[2],2);
                    auto zhi = loc_arr(iv[0], iv[1], iv[2]+1,2);
                    if (p.pos(2) > zhi) { // need to be careful here
                        p.idata(2) += 1;
                    } else if (p.pos(2) <= zlo) {
                        p.idata(2) -= 1;
                    }
#endif
                }
            });
        } // ParIter
    } // ipass

    Redistribute();
}

/*
  /brief Uses midpoint method to advance particles using cell-centered velocity.
*/
void
MappedPC::AdvectWithUND (MultiFab& vel_nd, int lev, Real dt, const MultiFab& a_xyz_loc)
{
    BL_PROFILE("MappedPC::AdvectWithND()");
    AMREX_ASSERT(lev >= 0 && lev < GetParticles().size());

    const auto dxi = this->ParticleContainerBase::Geom(0).InvCellSizeArray();

    auto probhi = this->ParticleContainerBase::Geom(0).ProbHiArray();
    auto problo = this->ParticleContainerBase::Geom(0).ProbLoArray();

    // Center of the annulus
    Real cx = 0.5 * (problo[0]+probhi[0]);
    Real cy = 0.5 * (problo[1]+probhi[1]);

    for (int ipass = 0; ipass < 2; ipass++)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& ptile = ParticlesAt(lev, pti);
            auto& aos  = ptile.GetArrayOfStructs();
            const int n = aos.numParticles();
            auto *p_pbox = aos().data();

            const auto loc_arr = a_xyz_loc.array(pti);
            const auto vel_nd_arr = vel_nd.const_array(pti);

            ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = p_pbox[i];

                if (p.id() <= 0) { return; }

                ParticleReal v[AMREX_SPACEDIM];
                cic_interpolate_mapped_z(p, problo, dxi, vel_nd_arr, loc_arr, v);

                if (ipass == 0)
                {
                    Real r = std::sqrt((p.pos(0)-cx)*(p.pos(0)-cx)  +(p.pos(1)-cy)*(p.pos(1)-cy));
                    Real theta = atan((p.pos(1)-cy)/(p.pos(0)-cx));

                    amrex::Print() << "UPDATING FROM " << p.pos(0) << " " << p.pos(1) <<
                                       " WITH RADIUS   " << r <<
                                       " WITH  THETA   " << theta*180/3.1415926 <<
                                       " AND VEL " << v[0] << " " << v[1] << std::endl;

                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.rdata(dim) = p.pos(dim);
                        p.pos(dim) += static_cast<ParticleReal>(ParticleReal(0.5)*dt*v[dim]);
                    }
                    r = std::sqrt((p.pos(0)-cx)*(p.pos(0)-cx)  +(p.pos(1)-cy)*(p.pos(1)-cy));
                    amrex::Print() << "          TO  " << p.pos(0) << " " << p.pos(1) <<
                                       " WITH RADIUS   " << r <<
                                       " AND VEL " << v[0] << " " << v[1] << std::endl;
                }
                else
                {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.pos(dim) = p.rdata(dim) + static_cast<ParticleReal>(dt*v[dim]);
                        p.rdata(dim) = v[dim];
                    }

                    // also update z-coordinate here
                    IntVect iv(
                       AMREX_D_DECL(p.idata(0),
                                    p.idata(1),
                                    p.idata(2)));

                    auto xlo = loc_arr(AMREX_D_DECL(iv[0]  , iv[1]  , iv[2]),0);
                    auto xhi = loc_arr(AMREX_D_DECL(iv[0]+1, iv[1]  , iv[2]),0);
                    auto ylo = loc_arr(AMREX_D_DECL(iv[0]  , iv[1]  , iv[2]),1);
                    auto yhi = loc_arr(AMREX_D_DECL(iv[0]  , iv[1]+1, iv[2]),1);

                    if (p.pos(0) > xhi) { // need to be careful here
                        p.idata(0) += 1;
                    } else if (p.pos(0) <= xlo) {
                        p.idata(0) -= 1;
                    }
                    if (p.pos(1) > yhi) { // need to be careful here
                        p.idata(1) += 1;
                    } else if (p.pos(1) <= ylo) {
                        p.idata(1) -= 1;
                    }

#if (AMREX_SPACEDIM == 3)
                    auto zlo = loc_arr(iv[0], iv[1], iv[2],2);
                    auto zhi = loc_arr(iv[0], iv[1], iv[2]+1,2);
                    if (p.pos(2) > zhi) { // need to be careful here
                        p.idata(2) += 1;
                    } else if (p.pos(2) <= zlo) {
                        p.idata(2) -= 1;
                    }
#endif
                }
            });
        } // ParIter
    } // ipass

    Redistribute();
}
