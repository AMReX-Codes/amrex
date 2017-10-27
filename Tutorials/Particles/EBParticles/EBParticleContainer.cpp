#include <math.h>

#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>

#include "EBParticleContainer.H"
#include "ebparticles_F.H"

using namespace amrex;

EBParticleContainer::
EBParticleContainer(const Geometry            & geom,
                              const DistributionMapping & dmap,
                              const BoxArray            & ba)
    : ParticleContainer<2*BL_SPACEDIM>(geom, dmap, ba)
{}

void EBParticleContainer::InitParticles(const Real& radius, 
                                                  const RealVect& center) {

    BL_PROFILE("EBParticleContainer::InitParticles");

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real* dx  = geom.CellSize();
    
    std::mt19937 mt(0451);
    std::uniform_real_distribution<double> dist(0.0, 2*M_PI);

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const Box& tile_box = mfi.tilebox();
        const RealBox tile_real_box { tile_box, dx, geom.ProbLo() };

        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id, tile_id)];
        
        const auto& boxlo = tile_box.smallEnd();
        ParticleType p;
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
            
            p.id() = ParticleType::NextID();
            p.cpu() = ParallelDescriptor::MyProc();
            
            p.pos(0) = tile_real_box.lo(0) + (iv[0]- boxlo[0] + 0.5)*dx[0];
            p.pos(1) = tile_real_box.lo(1) + (iv[1]- boxlo[1] + 0.5)*dx[1];
#if (BL_SPACEDIM == 3)
            p.pos(2) = tile_real_box.lo(2) + (iv[2]- boxlo[2] + 0.5)*dx[2];
#endif

            Real theta = dist(mt);
            Real speed = 20.0;

            p.rdata(0) = speed*cos(theta);
            p.rdata(1) = speed*sin(theta);
#if (BL_SPACEDIM == 3)
            p.rdata(2) = dist(mt);
#endif
            
            p.rdata(BL_SPACEDIM)   = 0;
            p.rdata(BL_SPACEDIM+1) = 0;
#if (BL_SPACEDIM == 3)
            p.rdata(BL_SPACEDIM+2) = 0;
#endif
            
            Real d = std::sqrt(AMREX_D_TERM((  p.pos(0) - center[0]) * (p.pos(0) - center[0]),
                                            + (p.pos(1) - center[1]) * (p.pos(1) - center[1]),
                                            + (p.pos(2) - center[2]) * (p.pos(2) - center[2])));
            
            if (d <= (radius - 1e-4)) {
                particle_tile.push_back(p);
            }
        }
    }
}

void EBParticleContainer::moveParticles(const Real dt) {

    BL_PROFILE("EBParticleContainer::moveParticles");

    const int lev = 0;
    const RealBox& prob_domain = Geom(lev).ProbDomain();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();
        amrex_move_particles(particles.data(), &Np, &dt,
                             prob_domain.lo(), prob_domain.hi());
    }
}

void EBParticleContainer::bounceWalls(const MultiFab& dummy,
                                      const MultiFab* volfrac,
                                      const MultiCutFab* bndrycent,
                                      std::array<const MultiCutFab*, AMREX_SPACEDIM>& areafrac,
                                      const Real dt) 
{
    
    BL_PROFILE("EBParticleContainer::bounceWalls");
    
    const int lev = 0;
    const Geometry& gm          = Geom(lev);
    const Real*     plo         = gm.ProbLo();
    const Real*     dx          = gm.CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MyConstParIter pti(*this, lev); pti.isValid(); ++pti) {
        const Box& bx = pti.tilebox();
        const AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        const auto& sfab = dynamic_cast <EBFArrayBox const&>((dummy)[pti]);
        const auto& flag = sfab.getEBCellFlagFab();
        
        if (flag.getType(bx) == FabType::regular) {
            continue;
        } else {
            amrex_bounce_walls(particles.data(), &Np, plo, dx, &dt,
                               flag.dataPtr(), flag.loVect(), flag.hiVect(),
                               (*bndrycent)[pti].dataPtr(), 
                               (*bndrycent)[pti].loVect(), (*bndrycent)[pti].hiVect(),
                               (*areafrac[0])[pti].dataPtr(), 
                               (*areafrac[0])[pti].loVect(), (*areafrac[0])[pti].hiVect(),
                               (*areafrac[1])[pti].dataPtr(), 
                               (*areafrac[1])[pti].loVect(), (*areafrac[1])[pti].hiVect());
        }
    }
}

void EBParticleContainer::writeParticles(int n) {

    BL_PROFILE("EBParticleContainer::writeParticles");

    const std::string& pltfile = amrex::Concatenate("particles", n, 5);
    WriteAsciiFile(pltfile);
}
