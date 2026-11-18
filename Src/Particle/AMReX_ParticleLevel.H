#ifndef AMREX_PARTICLELEVEL_H_
#define AMREX_PARTICLELEVEL_H_
#include <AMReX_Config.H>

#include <AMReX_Arena.H>

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace amrex {

/**
 * \brief What a particle tile needs to know about its particle container: the
 *        number and names of its Struct-of-Arrays components, and the memory arena
 *        it allocates from.
 *
 * A ParticleContainer owns one of these and shares it with its ParticleLevels, so
 * that the tiles they create are set up like the container. The container keeps it
 * up to date when components are added or the arena is changed.
 */
struct ParticleTileConfig
{
    //! number of Real and int SoA components added at runtime
    int num_runtime_real = 0;
    int num_runtime_int = 0;
    //! names of all the Real and int SoA components, compile-time ones first
    std::vector<std::string> soa_rdata_names;
    std::vector<std::string> soa_idata_names;
    //! arena of containers with a PolymorphicArenaAllocator (nullptr: not set yet)
    Arena* arena = nullptr;
};

/**
 * \brief The particle tiles of one AMR level, keyed by (grid index, local tile index).
 *
 * This is a std::map, and can be used like one. In addition, a tile created with
 * operator[] has all the components of its ParticleContainer, including the ones
 * added at runtime. For a ParticleContainer with a PolymorphicArenaAllocator, the
 * tile also uses the container's arena, if it has been set.
 *
 * Other ways of adding a tile, such as insert or emplace, add it as it is given.
 */
template <class Tile>
class ParticleLevelMap
    : public std::map<std::pair<int, int>, Tile>
{
    using Base = std::map<std::pair<int, int>, Tile>;

public:
    using key_type = typename Base::key_type;

    ParticleLevelMap () = default;

    explicit ParticleLevelMap (std::shared_ptr<ParticleTileConfig> a_config) noexcept
        : m_config(std::move(a_config))
    {}

    //! Return the tile at key, creating it with the container's components if needed.
    Tile& operator[] (key_type const& key)
    {
        auto [it, inserted] = Base::try_emplace(key);
        if (inserted && m_config) { defineTile(it->second); }
        return it->second;
    }

    void swap (ParticleLevelMap& other) noexcept
    {
        Base::swap(other);
        std::swap(m_config, other.m_config);
    }

    //! Set what the tiles created from now on are set up with.
    void setConfig (std::shared_ptr<ParticleTileConfig> a_config) noexcept
    {
        m_config = std::move(a_config);
    }

    [[nodiscard]] std::shared_ptr<ParticleTileConfig> const& config () const noexcept
    {
        return m_config;
    }

private:

    void defineTile (Tile& tile) const
    {
        auto& c = *m_config;
        if constexpr (Tile::ParticleType::is_rtsoa_particle) {
            // a ParticleTileRT needs its arena when it is defined
            if (c.arena != nullptr) {
                tile.define(c.num_runtime_real, c.num_runtime_int,
                            &c.soa_rdata_names, &c.soa_idata_names, c.arena);
            }
        } else if constexpr (Tile::has_polymorphic_allocator) {
            if (c.arena != nullptr) {
                tile.define(c.num_runtime_real, c.num_runtime_int,
                            &c.soa_rdata_names, &c.soa_idata_names, c.arena);
            } else {
                // add the components now; the arena is set when the tile is defined
                tile.defineComponents(c.num_runtime_real, c.num_runtime_int,
                                      &c.soa_rdata_names, &c.soa_idata_names);
            }
        } else {
            tile.define(c.num_runtime_real, c.num_runtime_int,
                        &c.soa_rdata_names, &c.soa_idata_names);
        }
    }

    std::shared_ptr<ParticleTileConfig> m_config;
};

}

#endif
