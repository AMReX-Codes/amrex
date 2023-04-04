#ifndef AMREX_PARITER_H_
#define AMREX_PARITER_H_
#include <AMReX_Config.H>

#include <AMReX_MFIter.H>
#include <AMReX_Gpu.H>

namespace amrex
{

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
class ParticleContainer;

template <bool is_const, int NStructReal, int NStructInt=0, int NArrayReal=0, int NArrayInt=0,
          template<class> class Allocator=DefaultAllocator>
class ParIterBase
    : public MFIter
{
private:

    using PCType = ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>;
    using ContainerRef    = typename std::conditional<is_const, PCType const&, PCType&>::type;
    using ParticleTileRef = typename std::conditional
        <is_const, typename PCType::ParticleTileType const&, typename PCType::ParticleTileType &>::type;
    using ParticleTilePtr = typename std::conditional
        <is_const, typename PCType::ParticleTileType const*, typename PCType::ParticleTileType *>::type;
    using AoSRef          = typename std::conditional
        <is_const, typename PCType::AoS const&, typename PCType::AoS&>::type;
    using SoARef          = typename std::conditional
        <is_const, typename PCType::SoA const&, typename PCType::SoA&>::type;

public:

    using ContainerType    = ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>;
    using ParticleTileType = typename ContainerType::ParticleTileType;
    using AoS              = typename ContainerType::AoS;
    using SoA              = typename ContainerType::SoA;
    using ParticleType     = typename ContainerType::ParticleType;
    using RealVector       = typename SoA::RealVector;
    using IntVector        = typename SoA::IntVector;
    using ParticleVector   = typename ContainerType::ParticleVector;

    ParIterBase (ContainerRef pc, int level);

    ParIterBase (ContainerRef pc, int level, MFItInfo& info);

#ifdef AMREX_USE_OMP
    void operator++ ()
    {
        if (dynamic) {
#pragma omp atomic capture
            m_pariter_index = nextDynamicIndex++;
        } else {
            ++m_pariter_index;
        }
        currentIndex = m_valid_index[m_pariter_index];
    }
#else
    void operator++ ()
    {
        ++m_pariter_index;
        currentIndex = m_valid_index[m_pariter_index];
#ifdef AMREX_USE_GPU
        Gpu::Device::setStreamIndex(currentIndex);
#endif
    }
#endif

    [[nodiscard]] ParticleTileRef GetParticleTile () const { return *m_particle_tiles[m_pariter_index]; }

    [[nodiscard]] AoSRef GetArrayOfStructs () const { return GetParticleTile().GetArrayOfStructs(); }

    [[nodiscard]] SoARef GetStructOfArrays () const { return GetParticleTile().GetStructOfArrays(); }

    [[nodiscard]] int numParticles () const { return GetArrayOfStructs().numParticles(); }

    [[nodiscard]] int numRealParticles () const { return GetArrayOfStructs().numRealParticles(); }

    [[nodiscard]] int numNeighborParticles () const { return GetArrayOfStructs().numNeighborParticles(); }

    [[nodiscard]] int GetLevel () const { return m_level; }

    [[nodiscard]] std::pair<int, int> GetPairIndex () const { return std::make_pair(this->index(), this->LocalTileIndex()); }

    [[nodiscard]] const Geometry& Geom (int lev) const { return m_pc.Geom(lev); }

protected:

    int m_level;
    int m_pariter_index;
    Vector<int> m_valid_index;
    Vector<ParticleTilePtr> m_particle_tiles;
    ContainerRef m_pc;
};

template <int NStructReal, int NStructInt=0, int NArrayReal=0, int NArrayInt=0,
          template<class> class Allocator=DefaultAllocator>
class ParIter
    : public ParIterBase<false,NStructReal,NStructInt, NArrayReal, NArrayInt, Allocator>
{
public:

    using ContainerType    = ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt,
                                               Allocator>;
    using ParticleTileType = typename ContainerType::ParticleTileType;
    using AoS              = typename ContainerType::AoS;
    using SoA              = typename ContainerType::SoA;
    using ParticleType     = typename ContainerType::ParticleType;
    using RealVector       = typename SoA::RealVector;
    using IntVector        = typename SoA::IntVector;

    ParIter (ContainerType& pc, int level)
        : ParIterBase<false, NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>(pc,level)
        {}

    ParIter (ContainerType& pc, int level, MFItInfo& info)
        : ParIterBase<false, NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>(pc,level,info)
        {}
};

template <int NStructReal, int NStructInt=0, int NArrayReal=0, int NArrayInt=0,
          template<class> class Allocator=DefaultAllocator>
class ParConstIter
    : public ParIterBase<true,NStructReal,NStructInt, NArrayReal, NArrayInt, Allocator>
{
public:

    using ContainerType    = ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt,
                                               Allocator>;
    using ParticleTileType = typename ContainerType::ParticleTileType;
    using AoS              = typename ContainerType::AoS;
    using SoA              = typename ContainerType::SoA;
    using ParticleType     = typename ContainerType::ParticleType;
    using RealVector       = typename SoA::RealVector;
    using IntVector        = typename SoA::IntVector;

    ParConstIter (ContainerType const& pc, int level)
        : ParIterBase<true,NStructReal,NStructInt,NArrayReal,NArrayInt,Allocator>(pc,level)
        {}

    ParConstIter (ContainerType const& pc, int level, MFItInfo& info)
        : ParIterBase<true,NStructReal,NStructInt,NArrayReal,NArrayInt,Allocator>(pc,level,info)
        {}
};

template <bool is_const, int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
ParIterBase<is_const, NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>::ParIterBase
  (ContainerRef pc, int level, MFItInfo& info)
    :
      MFIter(*pc.m_dummy_mf[level], pc.do_tiling ? info.EnableTiling(pc.tile_size) : info),
      m_level(level),
      m_pariter_index(0),
      m_pc(pc)
{
    auto& particles = pc.GetParticles(level);

    int start = dynamic ? 0 : beginIndex;
    for (int i = start; i < endIndex; ++i)
    {
        int grid = (*index_map)[i];
        int tile = local_tile_index_map ? (*local_tile_index_map)[i] : 0;
        auto key = std::make_pair(grid,tile);
        auto f = particles.find(key);
        if (f != particles.end() && f->second.numParticles() > 0)
        {
            m_valid_index.push_back(i);
            m_particle_tiles.push_back(&(f->second));
        }
    }

    if (m_valid_index.empty())
    {
        endIndex = beginIndex;
    }
    else
    {
        currentIndex = beginIndex = m_valid_index.front();
        if (dynamic) {
#ifdef AMREX_USE_OMP
            int ind = omp_get_thread_num();
            m_pariter_index += ind;
            if (ind < m_valid_index.size()) {
                currentIndex = beginIndex = m_valid_index[ind];
            } else {
                currentIndex = endIndex;
            }
            for (int i = 0; i < omp_get_num_threads(); ++i) {
                m_valid_index.push_back(endIndex);
            }
#endif
        }
        m_valid_index.push_back(endIndex);
    }
}

template <bool is_const, int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
ParIterBase<is_const, NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>::ParIterBase
  (ContainerRef pc, int level)
    :
    MFIter(*pc.m_dummy_mf[level],
           pc.do_tiling ? pc.tile_size : IntVect::TheZeroVector()),
    m_level(level),
    m_pariter_index(0),
    m_pc(pc)
{
    auto& particles = pc.GetParticles(level);

    for (int i = beginIndex; i < endIndex; ++i)
    {
        int grid = (*index_map)[i];
        int tile = local_tile_index_map ? (*local_tile_index_map)[i] : 0;
        auto key = std::make_pair(grid,tile);
        auto f = particles.find(key);
        if (f != particles.end() && f->second.numParticles() > 0)
        {
            m_valid_index.push_back(i);
            m_particle_tiles.push_back(&(f->second));
        }
    }

    if (m_valid_index.empty())
    {
        endIndex = beginIndex;
    }
    else
    {
        currentIndex = beginIndex = m_valid_index.front();
        m_valid_index.push_back(endIndex);
    }
}

}

#endif
