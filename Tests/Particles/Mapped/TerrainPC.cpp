#include "TerrainPC.H"

#include <AMReX_TracerParticle_mod_K.H>

using namespace amrex;

static constexpr int NSR = 2*AMREX_SPACEDIM;
static constexpr int NSI = 1;
static constexpr int NAR = 0;
static constexpr int NAI = 0;


void
TerrainPC::
InitParticles (MultiFab& a_z_loc)
{
    BL_PROFILE("TerrainPC::InitParticles");

    const int lev = 0;

    //auto problo = this->ParticleContainerBase::Geom(0).ProbLo();

    const auto dx     = Geom(lev).CellSizeArray();
    const auto domain = Geom(lev).Domain();

    auto dom_lo = lbound(domain);
    auto dom_hi = ubound(domain);

    for(MFIter mfi(a_z_loc); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = enclosedCells(mfi.tilebox());

        Gpu::HostVector<ParticleType> host_particles;
        std::array<Gpu::HostVector<ParticleReal>, NAR> host_real;
        std::array<Gpu::HostVector<int>, NAI> host_int;

        auto height_arr = a_z_loc.array(mfi);

        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
        {
            int i = iv[0]; int j = iv[1];
#if (AMREX_SPACEDIM == 2)
            int k = 0;
            if (iv[0] == 0 && iv[1] >= dom_lo.y && iv[1] <= dom_hi.y) {
#elif (AMREX_SPACEDIM == 3)
            int k = iv[2];
            if (iv[0] == 0 && iv[1] == 0 && iv[2] >= dom_lo.z && iv[2] <= dom_hi.z) {
#endif

                // This is the physical location of the center of the cell
                Real x = (i+0.5) * dx[0];

                ParticleType p;
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();

                p.rdata(TerrainRealIdx::vx) = Real(0.0);
                p.rdata(TerrainRealIdx::vy) = Real(0.0);

                p.pos(0) = x;

#if (AMREX_SPACEDIM == 2)
                // Put the particle at the cell center for now
                Real y = 0.25*(height_arr(i,j  ,k) + height_arr(i+1,j  ,k) +
                               height_arr(i,j+1,k) + height_arr(i+1,j+1,k) );
                p.pos(1) = y;;
                p.idata(0) = j;  // particles carry their j-index
                amrex::Print() << "MAKING PARTICLE WITH IDATA POS " << p.idata(0) << " " << p.pos(0) << " " << p.pos(1) << std::endl;
#elif (AMREX_SPACEDIM == 3)

                Real y = (j+0.5) * dx[1];
                p.pos(1) = y;

                // Put the particle at the cell center for now
                Real z = 0.125*( height_arr(i  ,j,k  ) + height_arr(i  ,j+1,k  ) +
                                 height_arr(i+1,j,k  ) + height_arr(i+1,j+1,k  ) +
                                 height_arr(i  ,j,k+1) + height_arr(i  ,j+1,k+1) +
                                 height_arr(i+1,j,k+1) + height_arr(i+1,j+1,k+1) );
                p.pos(2) = z;

                p.idata(0) = k;  // particles carry their k-index

                p.rdata(TerrainRealIdx::vz) = Real(0.);
                amrex::Print() << "MAKING PARTICLE WITH IDATA POS " << p.idata(0) << " " << p.pos(0) << " " << p.pos(2) << std::endl;
#endif

                host_particles.push_back(p);
                for (int nr = 0; nr < NAR; ++nr) {
                    host_real[nr].push_back(p.rdata(nr));
                }
                for (int ni = 0; ni < NAI; ++ni) {
                    host_int[ni].push_back(p.idata(ni));
                }

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
  /brief Uses midpoint method to advance particles using umac.
*/
void
TerrainPC::AdvectWithUmac (MultiFab* umac, int lev, Real dt, const MultiFab& a_z_loc)
{
    BL_PROFILE("TerrainPC::AdvectWithUmac()");
    AMREX_ASSERT(lev >= 0 && lev < GetParticles().size());

    AMREX_D_TERM(AMREX_ASSERT(umac[0].nGrow() >= 1);,
                 AMREX_ASSERT(umac[1].nGrow() >= 1);,
                 AMREX_ASSERT(umac[2].nGrow() >= 1););

    const Geometry& geom = m_gdb->Geom(lev);

    const auto dxi = geom.InvCellSizeArray();

    auto problo = this->ParticleContainerBase::Geom(0).ProbLoArray();

    AMREX_ALWAYS_ASSERT(OnSameGrids(lev, umac[0]));

    Vector<MultiFab*> umac_pointer(AMREX_SPACEDIM);
    {
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            umac_pointer[i] = &umac[i];
        }
    }

    for (int ipass = 0; ipass < 2; ipass++)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti)
        {
            int grid    = pti.index();
            auto& ptile = ParticlesAt(lev, pti);
            auto& aos  = ptile.GetArrayOfStructs();
            const int n = aos.numParticles();
            auto *p_pbox = aos().data();
            const FArrayBox* fab[AMREX_SPACEDIM] = { AMREX_D_DECL(&((*umac_pointer[0])[grid]),
                                                                  &((*umac_pointer[1])[grid]),
                                                                  &((*umac_pointer[2])[grid])) };

            const auto height_arr = a_z_loc.array(pti);

            //array of these pointers to pass to the GPU
            amrex::GpuArray<amrex::Array4<const Real>, AMREX_SPACEDIM>
                const umacarr {{AMREX_D_DECL( (*fab[0]).array(),
                                              (*fab[1]).array(),
                                              (*fab[2]).array() )}};

            ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = p_pbox[i];

                if (p.id() <= 0) { return; }

                ParticleReal v[AMREX_SPACEDIM];
                // amrex::Print() << "CALLING MAC INTERP WITH IDATA " << p.pos(0) << " " <<
                //                    p.pos(AMREX_SPACEDIM-1) << " " << p.idata(0) << std::endl;
                mac_interpolate_mapped_z(p, problo, dxi, umacarr, height_arr, v);

                if (ipass == 0)
                {
                    amrex::Print() << "FROM  " << p.pos(0) << " " << p.pos(AMREX_SPACEDIM-1) << std::endl;
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.rdata(dim) = p.pos(dim);
                        p.pos(dim) += static_cast<ParticleReal>(ParticleReal(0.5)*dt*v[dim]);
                    }
                }
                else
                {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.pos(dim) = p.rdata(dim) + static_cast<ParticleReal>(dt*v[dim]);
                        p.rdata(dim) = v[dim];
                    }
                    amrex::Print() << " TO  " << p.pos(0) << " " << p.pos(AMREX_SPACEDIM-1) << std::endl;

#if (AMREX_SPACEDIM == 2)
                    IntVect iv( int(amrex::Math::floor((p.pos(0)-problo[0])*dxi[0])), p.idata(0) );
                    iv[0] += geom.Domain().smallEnd()[0];

                    auto ylo = Real(0.5) * (height_arr(iv[0],iv[1]  ,0) + height_arr(iv[0]+1,iv[1]  ,0));
                    auto yhi = Real(0.5) * (height_arr(iv[0],iv[1]+1,0) + height_arr(iv[0]+1,iv[1]+1,0));
                    amrex::Print() << "OLD IDATA " << p.idata(0) << std::endl;
                    if (p.pos(1) > yhi) { // need to be careful here
                        p.idata(0) += 1;
                    } else if (p.pos(1) <= ylo) {
                        p.idata(0) -= 1;
                    }
                    amrex::Print() << "NEW IDATA " << p.idata(0) << std::endl;

#elif (AMREX_SPACEDIM == 3)
                    IntVect iv( int(amrex::Math::floor((p.pos(0)-problo[0])*dxi[0])),
                                int(amrex::Math::floor((p.pos(1)-problo[1])*dxi[1])),
                                p.idata(0) );
                    iv[0] += geom.Domain().smallEnd()[0];
                    iv[1] += geom.Domain().smallEnd()[1];

                    auto zlo = Real(0.25) * (height_arr(iv[0],iv[1]  ,iv[2]  ) + height_arr(iv[0]+1,iv[1]  ,iv[2]  ) +
                                             height_arr(iv[0],iv[1]+1,iv[2]  ) + height_arr(iv[0]+1,iv[1]+1,iv[2]  ));
                    auto zhi = Real(0.25) * (height_arr(iv[0],iv[1]  ,iv[2]+1) + height_arr(iv[0]+1,iv[1]  ,iv[2]+1) +
                                             height_arr(iv[0],iv[1]+1,iv[2]+1) + height_arr(iv[0]+1,iv[1]+1,iv[2]+1));
                    if (p.pos(2) > zhi) { // need to be careful here
                        p.idata(0) += 1;
                    } else if (p.pos(2) <= zlo) {
                        p.idata(0) -= 1;
                    }
#endif
                }
            });
        } // pti
    } // ipass

    Redistribute();
}

/*
  /brief Uses midpoint method to advance particles using cell-centered velocity.
*/
void
TerrainPC::AdvectWithUCC (MultiFab& vel_cc, int lev, Real dt, const MultiFab& a_z_loc)
{
    BL_PROFILE("TerrainPC::AdvectWithCC()");
    AMREX_ASSERT(lev >= 0 && lev < GetParticles().size());

    const Geometry& geom = m_gdb->Geom(lev);

    const auto dxi = this->ParticleContainerBase::Geom(lev).InvCellSizeArray();

    auto problo = this->ParticleContainerBase::Geom(0).ProbLoArray();

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

            const auto height_arr = a_z_loc.array(pti);
            const auto vel_cc_arr = vel_cc.const_array(pti);

            ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = p_pbox[i];

                if (p.id() <= 0) { return; }

                ParticleReal v[AMREX_SPACEDIM];
                cic_interpolate_mapped_z(p, problo, dxi, vel_cc_arr, height_arr, v);

                if (ipass == 0)
                {
                    amrex::Print() << "        FROM  " << p.pos(0) << " " << p.pos(1) << std::endl;
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.rdata(dim) = p.pos(dim);
                        p.pos(dim) += static_cast<ParticleReal>(ParticleReal(0.5)*dt*v[dim]);
                    }
                    amrex::Print() << "          TO  " << p.pos(0) << " " << p.pos(1) << std::endl;
                }
                else
                {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.pos(dim) = p.rdata(dim) + static_cast<ParticleReal>(dt*v[dim]);
                        p.rdata(dim) = v[dim];
                    }

#if (AMREX_SPACEDIM == 2)
                    IntVect iv( int(amrex::Math::floor((p.pos(0)-problo[0])*dxi[0])), p.idata(0) );
                    iv[0] += geom.Domain().smallEnd()[0];

                    auto ylo = Real(0.5) * (height_arr(iv[0],iv[1]  ,0) + height_arr(iv[0]+1,iv[1]  ,0));
                    auto yhi = Real(0.5) * (height_arr(iv[0],iv[1]+1,0) + height_arr(iv[0]+1,iv[1]+1,0));
                    if (p.pos(1) > yhi) { // need to be careful here
                        p.idata(0) += 1;
                    } else if (p.pos(1) <= ylo) {
                        p.idata(0) -= 1;
                    }

#elif (AMREX_SPACEDIM == 3)
                    IntVect iv( int(amrex::Math::floor((p.pos(0)-problo[0])*dxi[0])),
                                int(amrex::Math::floor((p.pos(1)-problo[1])*dxi[1])),
                                p.idata(0) );
                    iv[0] += geom.Domain().smallEnd()[0];
                    iv[1] += geom.Domain().smallEnd()[1];

                    auto zlo = Real(0.25) * (height_arr(iv[0],iv[1]  ,iv[2]  ) + height_arr(iv[0]+1,iv[1]  ,iv[2]  ) +
                                             height_arr(iv[0],iv[1]+1,iv[2]  ) + height_arr(iv[0]+1,iv[1]+1,iv[2]  ));
                    auto zhi = Real(0.25) * (height_arr(iv[0],iv[1]  ,iv[2]+1) + height_arr(iv[0]+1,iv[1]  ,iv[2]+1) +
                                             height_arr(iv[0],iv[1]+1,iv[2]+1) + height_arr(iv[0]+1,iv[1]+1,iv[2]+1));
                    if (p.pos(2) > zhi) { // need to be careful here
                        p.idata(0) += 1;
                    } else if (p.pos(2) <= zlo) {
                        p.idata(0) -= 1;
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
TerrainPC::AdvectWithUND (MultiFab& vel_nd, int lev, Real dt, const MultiFab& a_z_loc)
{
    BL_PROFILE("TerrainPC::AdvectWithCC()");
    AMREX_ASSERT(lev >= 0 && lev < GetParticles().size());

    const Geometry& geom = m_gdb->Geom(lev);

    const auto dxi = this->ParticleContainerBase::Geom(0).InvCellSizeArray();

    auto problo = this->ParticleContainerBase::Geom(0).ProbLoArray();

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

            const auto height_arr = a_z_loc.array(pti);
            const auto vel_nd_arr = vel_nd.const_array(pti);

            ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = p_pbox[i];

                if (p.id() <= 0) { return; }

                ParticleReal v[AMREX_SPACEDIM];
                cic_interpolate_mapped_z(p, problo, dxi, vel_nd_arr, height_arr, v);

                if (ipass == 0)
                {
                    amrex::Print() << "        FROM  " << p.pos(0) << " " << p.pos(AMREX_SPACEDIM-1) << std::endl;
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.rdata(dim) = p.pos(dim);
                        p.pos(dim) += static_cast<ParticleReal>(ParticleReal(0.5)*dt*v[dim]);
                    }
                }
                else
                {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.pos(dim) = p.rdata(dim) + static_cast<ParticleReal>(dt*v[dim]);
                        p.rdata(dim) = v[dim];
                    }
                    amrex::Print() << "          TO  " << p.pos(0) << " " << p.pos(AMREX_SPACEDIM-1) << std::endl;


#if (AMREX_SPACEDIM == 2)
                    IntVect iv( int(amrex::Math::floor((p.pos(0)-problo[0])*dxi[0])), p.idata(0) );
                    iv[0] += geom.Domain().smallEnd()[0];

                    auto ylo = Real(0.5) * (height_arr(iv[0],iv[1]  ,0) + height_arr(iv[0]+1,iv[1]  ,0));
                    auto yhi = Real(0.5) * (height_arr(iv[0],iv[1]+1,0) + height_arr(iv[0]+1,iv[1]+1,0));
                    if (p.pos(1) > yhi) { // need to be careful here
                        p.idata(0) += 1;
                    } else if (p.pos(1) <= ylo) {
                        p.idata(0) -= 1;
                    }

#elif (AMREX_SPACEDIM == 3)
                    IntVect iv( int(amrex::Math::floor((p.pos(0)-problo[0])*dxi[0])),
                                int(amrex::Math::floor((p.pos(1)-problo[1])*dxi[1])),
                                p.idata(0) );
                    auto zlo = Real(0.25) * (height_arr(iv[0],iv[1]  ,iv[2]  ) + height_arr(iv[0]+1,iv[1]  ,iv[2]  ) +
                                             height_arr(iv[0],iv[1]+1,iv[2]  ) + height_arr(iv[0]+1,iv[1]+1,iv[2]  ));
                    auto zhi = Real(0.25) * (height_arr(iv[0],iv[1]  ,iv[2]+1) + height_arr(iv[0]+1,iv[1]  ,iv[2]+1) +
                                             height_arr(iv[0],iv[1]+1,iv[2]+1) + height_arr(iv[0]+1,iv[1]+1,iv[2]+1));
                    if (p.pos(2) > zhi) { // need to be careful here
                        p.idata(0) += 1;
                    } else if (p.pos(2) <= zlo) {
                        p.idata(0) -= 1;
                    }
#endif
                }
            });
        } // ParIter
    } // ipass

    Redistribute();
}
