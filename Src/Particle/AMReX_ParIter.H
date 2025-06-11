#ifndef AMREX_PARITER_H_
#define AMREX_PARITER_H_
#include <AMReX_Config.H>

#include <AMReX_MFIter.H>
#include <AMReX_Gpu.H>

namespace amrex
{

template <typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
class ParticleContainer_impl;

template <int T_NReal, int T_NInt>
struct Particle;

template <int NArrayReal, int NArrayInt>
struct SoAParticle;

struct DefaultAssignor;

// for backwards compatibility
template <int T_NStructReal, int T_NStructInt=0, int T_NArrayReal=0, int T_NArrayInt=0,
          template<class> class Allocator=DefaultAllocator, class CellAssignor=DefaultAssignor>
using ParticleContainer = ParticleContainer_impl<Particle<T_NStructReal, T_NStructInt>, T_NArrayReal, T_NArrayInt, Allocator, CellAssignor>;

template <bool is_const, typename T_ParticleType, int NArrayReal=0, int NArrayInt=0,
          template<class> class Allocator=DefaultAllocator, class CellAssignor=DefaultAssignor>
class ParIterBase_impl
    : public MFIter
{

private:

    using PCType = ParticleContainer_impl<T_ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>;
    using ContainerRef    = std::conditional_t<is_const, PCType const&, PCType&>;
    using ContainerPtr    = std::conditional_t<is_const, PCType const*, PCType*>;
    using ParticleTileRef = std::conditional_t
        <is_const, typename PCType::ParticleTileType const&, typename PCType::ParticleTileType &>;
    using ParticleTilePtr = std::conditional_t
        <is_const, typename PCType::ParticleTileType const*, typename PCType::ParticleTileType *>;
    using AoSRef          = std::conditional_t
        <is_const, typename PCType::AoS const&, typename PCType::AoS&>;
    using SoARef          = std::conditional_t
        <is_const, typename PCType::SoA const&, typename PCType::SoA&>;

public:

    using ContainerType    = ParticleContainer_impl<T_ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>;
    using ParticleTileType = typename ContainerType::ParticleTileType;
    using AoS              = typename ContainerType::AoS;
    using SoA              = typename ContainerType::SoA;
    using ParticleType     = T_ParticleType;
    using RealVector       = typename SoA::RealVector;
    using IntVector        = typename SoA::IntVector;
    using ParticleVector   = typename ContainerType::ParticleVector;
    static constexpr int NStructReal = ParticleType::NReal;
    static constexpr int NStructInt = ParticleType::NInt;

    ParIterBase_impl (ContainerRef pc, int level);

    ParIterBase_impl (ContainerRef pc, int level, MFItInfo& info);

    void operator++ ()
    {
        m_particle_current_tile = nullptr;
        auto& particles = m_pc->GetParticles(m_level);
        while (true) {
            MFIter::operator++();
            if (isValid()) {
                auto key = std::make_pair(index(), LocalTileIndex());
                auto f = particles.find(key);
                if (f != particles.end() && f->second.numParticles() > 0) {
                    m_particle_current_tile = &(f->second);
                    break; // break if isValid() and non-empty particle tile found
                }
            } else {
                break; // break if isValid() returns false
            }
        }
    }

    [[nodiscard]] ParticleTileRef GetParticleTile () const { return *m_particle_current_tile; }

    [[nodiscard]] AoSRef GetArrayOfStructs () const { return GetParticleTile().GetArrayOfStructs(); }

    [[nodiscard]] SoARef GetStructOfArrays () const { return GetParticleTile().GetStructOfArrays(); }

    [[nodiscard]] int numParticles () const { return GetParticleTile().numParticles(); }

    [[nodiscard]] int numRealParticles () const { return GetParticleTile().numRealParticles(); }

    [[nodiscard]] int numNeighborParticles () const { return GetParticleTile().numNeighborParticles(); }

    [[nodiscard]] int GetLevel () const { return m_level; }

    [[nodiscard]] std::pair<int, int> GetPairIndex () const { return std::make_pair(this->index(), this->LocalTileIndex()); }

    [[nodiscard]] const Geometry& Geom (int lev) const { return m_pc->Geom(lev); }

protected:

    int m_level;
    ParticleTilePtr m_particle_current_tile = nullptr;
    ContainerPtr m_pc;
};

template <typename T_ParticleType, int NArrayReal=0, int NArrayInt=0,
          template<class> class Allocator=DefaultAllocator, class T_CellAssignor=DefaultAssignor>
class ParIter_impl
    : public ParIterBase_impl<false, T_ParticleType, NArrayReal, NArrayInt, Allocator, T_CellAssignor>
{
public:

    using ParticleType=T_ParticleType;
    using CellAssignor=T_CellAssignor;
    static constexpr int NStructReal = ParticleType::NReal;
    static constexpr int NStructInt = ParticleType::NInt;

    using ContainerType    = ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>;
    using ConstParticleType = typename ContainerType::ConstParticleType;
    using ParticleTileType = typename ContainerType::ParticleTileType;
    using AoS              = typename ContainerType::AoS;
    using SoA              = typename ContainerType::SoA;
    using RealVector       = typename SoA::RealVector;
    using IntVector        = typename SoA::IntVector;

    ParIter_impl (ContainerType& pc, int level)
        : ParIterBase_impl<false, T_ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>(pc,level)
        {}

    ParIter_impl (ContainerType& pc, int level, MFItInfo& info)
        : ParIterBase_impl<false, T_ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>(pc,level,info)
        {}
};

template <typename T_ParticleType, int NArrayReal=0, int NArrayInt=0,
          template<class> class Allocator=DefaultAllocator, class T_CellAssignor=DefaultAssignor>
class ParConstIter_impl
    : public ParIterBase_impl<true,T_ParticleType, NArrayReal, NArrayInt, Allocator, T_CellAssignor>
{
public:

    using ParticleType     = T_ParticleType;
    using CellAssignor     = T_CellAssignor;
    using ContainerType    = ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>;
    using ParticleTileType = typename ContainerType::ParticleTileType;
    using AoS              = typename ContainerType::AoS;
    using SoA              = typename ContainerType::SoA;
    using RealVector       = typename SoA::RealVector;
    using IntVector        = typename SoA::IntVector;

    ParConstIter_impl (ContainerType const& pc, int level)
        : ParIterBase_impl<true,ParticleType,NArrayReal,NArrayInt,Allocator,CellAssignor>(pc,level)
        {}

    ParConstIter_impl (ContainerType const& pc, int level, MFItInfo& info)
        : ParIterBase_impl<true,ParticleType,NArrayReal,NArrayInt,Allocator, CellAssignor>(pc,level,info)
        {}
};

template <bool is_const, typename ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
ParIterBase_impl<is_const, ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::ParIterBase_impl
  (ContainerRef pc, int level, MFItInfo& info)
    :
      MFIter(*pc.m_dummy_mf[level], pc.do_tiling ? info.EnableTiling(pc.tile_size) : info),
      m_level(level),
      m_pc(&pc)
{
    if (isValid()) {
        auto& particles = m_pc->GetParticles(m_level);
        while (true) {
            auto key = std::make_pair(index(), LocalTileIndex());
            auto f = particles.find(key);
            if (f != particles.end() && f->second.numParticles() > 0) {
                m_particle_current_tile = &(f->second);
                break; // break if isValid() and non-empty particle tile found
            }
            MFIter::operator++();
            if (!isValid()) { break; } // break if isValid() returns false
        }
    }
}

template <bool is_const, typename T_ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
ParIterBase_impl<is_const, T_ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::ParIterBase_impl
  (ContainerRef pc, int level)
    :
    MFIter(*pc.m_dummy_mf[level],
           pc.do_tiling ? pc.tile_size : IntVect::TheZeroVector()),
    m_level(level),
    m_pc(&pc)
{
    if (isValid()) {
        auto& particles = m_pc->GetParticles(m_level);
        while (true) {
            auto key = std::make_pair(index(), LocalTileIndex());
            auto f = particles.find(key);
            if (f != particles.end() && f->second.numParticles() > 0) {
                m_particle_current_tile = &(f->second);
                break;
            }
            MFIter::operator++();
            if (!isValid()) { break; }
        }
    }
}

template <bool is_const, int T_NStructReal, int T_NStructInt, int T_NArrayReal=0, int T_NArrayInt=0,
          template<class> class Allocator=DefaultAllocator, class CellAssignor=DefaultAssignor>
using ParIterBase = ParIterBase_impl<is_const, Particle<T_NStructReal, T_NStructInt>, T_NArrayReal, T_NArrayInt, Allocator, CellAssignor>;

template <bool is_const, int T_NArrayReal=0, int T_NArrayInt=0,
          template<class> class Allocator=DefaultAllocator, class CellAssignor=DefaultAssignor>
using ParIterBaseSoA = ParIterBase_impl<is_const,SoAParticle<T_NArrayReal, T_NArrayInt>, T_NArrayReal, T_NArrayInt, Allocator, CellAssignor>;

template <int T_NStructReal, int T_NStructInt=0, int T_NArrayReal=0, int T_NArrayInt=0,
          template<class> class Allocator=DefaultAllocator, class CellAssignor=DefaultAssignor>
using ParConstIter = ParConstIter_impl<Particle<T_NStructReal, T_NStructInt>, T_NArrayReal, T_NArrayInt, Allocator, CellAssignor>;

template <int T_NArrayReal, int T_NArrayInt, template<class> class Allocator=DefaultAllocator, class CellAssignor=DefaultAssignor>
using ParConstIterSoA = ParConstIter_impl<SoAParticle<T_NArrayReal, T_NArrayInt>, T_NArrayReal, T_NArrayInt, Allocator, CellAssignor>;

template <int T_NStructReal, int T_NStructInt=0, int T_NArrayReal=0, int T_NArrayInt=0,
          template<class> class Allocator=DefaultAllocator, class CellAssignor=DefaultAssignor>
using ParIter = ParIter_impl<Particle<T_NStructReal, T_NStructInt>, T_NArrayReal, T_NArrayInt, Allocator, CellAssignor>;

template <int T_NArrayReal, int T_NArrayInt, template<class> class Allocator=DefaultAllocator, class CellAssignor=DefaultAssignor>
using ParIterSoA = ParIter_impl<SoAParticle<T_NArrayReal, T_NArrayInt>, T_NArrayReal, T_NArrayInt, Allocator, CellAssignor>;

}

#endif
