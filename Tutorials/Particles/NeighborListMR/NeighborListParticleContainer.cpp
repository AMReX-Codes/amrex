#include "NeighborListParticleContainer.H"

#include "neighbor_list_F.H"

using namespace amrex;

constexpr Real NeighborListParticleContainer::min_r;
constexpr Real NeighborListParticleContainer::cutoff;

std::ostream& amrex::operator <<(std::ostream& stream, const Tag& tag)
{
    stream << "(" << tag.lev << ", " << tag.grid << ", " << tag.tile << ")";
    return stream;
}

NeighborListParticleContainer::
NeighborListParticleContainer(const Geometry            & geom,
                              const DistributionMapping & dmap,
                              const BoxArray            & ba,
                              int                         ncells)
    : NeighborParticleContainer<2*BL_SPACEDIM, 0> 
    (geom, dmap, ba, ncells)
{}

NeighborListParticleContainer::
NeighborListParticleContainer(const Vector<Geometry>            & geom,
                              const Vector<DistributionMapping> & dmap,
                              const Vector<BoxArray>            & ba,
                              const Vector<int>                 & rr,
                              int                               ncells)
    : NeighborParticleContainer<2*BL_SPACEDIM, 0> 
    (geom, dmap, ba, rr, ncells)
{}

void NeighborListParticleContainer::InitParticles()
{
    BL_PROFILE("NeighborListParticleContainer::InitParticles");
    
    const int lev = 0;  // initialize the particles on level 0
    const Geometry& geom = Geom(lev);
    const Real* dx  = geom.CellSize();
    
    std::mt19937 mt(0451);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

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
            p.rdata(0) = dist(mt);
            p.rdata(1) = dist(mt);
#if (BL_SPACEDIM == 3)
            p.rdata(2) = dist(mt);
#endif
            
            p.rdata(BL_SPACEDIM)   = 0;
            p.rdata(BL_SPACEDIM+1) = 0;
#if (BL_SPACEDIM == 3)
            p.rdata(BL_SPACEDIM+2) = 0;
#endif
            
            particle_tile.push_back(p);
        }
    }

    Redistribute();
}

void NeighborListParticleContainer::computeForces()
{
    BL_PROFILE("NeighborListParticleContainer::computeForces");
    
    int num_levs = finestLevel() + 1;
    for (int lev = 0; lev < num_levs; ++lev) {
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            AoS& particles = pti.GetArrayOfStructs();
            int Np = particles.size();
            PairIndex index(pti.index(), pti.LocalTileIndex());
            int Nn = neighbors[lev][index].size();
            amrex_compute_forces(particles.data(), &Np, 
                                 neighbors[lev][index].dataPtr(), &Nn, 
                                 &cutoff, &min_r);
        }
    }
}

void NeighborListParticleContainer::computeForcesNL()
{
    BL_PROFILE("NeighborListParticleContainer::computeForcesNL");
    
    buildNeighborList(CheckPair);
    
    int num_levs = finestLevel() + 1;
    for (int lev = 0; lev < num_levs; ++lev) {
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MyParIter pti(*this, lev, MFItInfo().SetDynamic(false)); pti.isValid(); ++pti) {
            PairIndex index(pti.index(), pti.LocalTileIndex());
            AoS& particles = pti.GetArrayOfStructs();
            int Np = particles.size();
            int Nn = neighbors[lev][index].size();
            int size = neighbor_list[lev][index].size();
            amrex_compute_forces_nl(particles.data(), &Np, 
                                    neighbors[lev][index].dataPtr(), &Nn,
                                    neighbor_list[lev][index].dataPtr(), &size, 
                                    &cutoff, &min_r);
        }
    }
}

void NeighborListParticleContainer::moveParticles(const Real dt)
{
    BL_PROFILE("NeighborListParticleContainer::moveParticles");
    
    int num_levs = finestLevel() + 1;
    for (int lev = 0; lev < num_levs; ++lev) {
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
}

void NeighborListParticleContainer::writeParticles(int n)
{
    BL_PROFILE("NeighborListParticleContainer::writeParticles");
    const std::string& pltfile = amrex::Concatenate("particles", n, 5);
    WriteAsciiFile(pltfile);
}

Vector<Tag>
NeighborListParticleContainer::getNeighborGrids(const ParticleType& p,
                                                const Tag src,
                                                const int nGrow)
{
    return getNeighborGrids(p, src, IntVect(AMREX_D_DECL(nGrow, nGrow, nGrow)));
}

Vector< Tag >
NeighborListParticleContainer::getNeighborGrids(const ParticleType& p,
                                                const Tag src,
                                                const IntVect& nGrow)
{
    Vector< Tag > grids;
    std::vector< std::pair<int, Box> > isects;
    Box tbx;
    IntVect ref_fac = IntVect(AMREX_D_DECL(1,1,1));
    const int num_levels = finestLevel() + 1;
    for (int lev = 0; lev < num_levels; ++lev)
    {
        if (lev > 0) ref_fac *= this->GetParGDB()->refRatio(lev-1);
        const BoxArray& ba = this->ParticleBoxArray(lev);
        const IntVect& iv = this->Index(p, lev);
        bool first_only = false;
        Box pbox = Box(iv, iv);
        ba.intersections(pbox, isects, first_only, ref_fac*nGrow);
        for (const auto& isec : isects)
        {
            for (IntVect cell = pbox.smallEnd(); cell <= pbox.bigEnd(); pbox.next(cell))
            {
                if (isec.second.contains(cell))
                {
                    int tile = getTileIndex(cell, isec.second,
                                            this->do_tiling, this->tile_size, tbx);
                    auto nbor = Tag(lev, isec.first, tile);
                    if (src != nbor) grids.push_back(nbor);
                }
            }
        }
    }
    return grids;
}

void
NeighborListParticleContainer::getNeighborParticles(const int nGrow)
{
    BL_PROFILE("getNeighborParticles()");

    const int num_levels = finestLevel() + 1;
    for (int lev = 0; lev < num_levels; ++lev)
    {
        for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            AoS& particles = pti.GetArrayOfStructs();
            auto grid = pti.index();
            auto tile = pti.LocalTileIndex();
            int np = particles.size();
            for (int i = 0; i < np; ++i)
            {
                auto src = Tag(lev, grid, tile);
                Vector<Tag> nbors = getNeighborGrids(particles[i], src, nGrow);
                if (nbors.size() == 0) continue;
                
                amrex::Print() << "Particle " << particles[i].id() << " is on " << src
                               << " and goes to: ";
                for (const auto& nbor : nbors)
                {
                    BL_ASSERT(nbor != src);
                    amrex::Print() << nbor << " ";
                }
                amrex::Print() << "\n";
            }
        }
    }
}
