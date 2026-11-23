#ifndef AMREX_PARTICLELEVEL_H_
#define AMREX_PARTICLELEVEL_H_
#include <AMReX_Config.H>

#include <AMReX_Arena.H>

#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
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
 * \brief The allocator of a ParticleLevel. It allocates memory like std::allocator,
 *        and sets up the particle tiles that the map creates without being given one.
 *
 * std::map creates every element through its allocator. When a tile is created from
 * nothing, as std::map::operator[] and try_emplace(key) do, it is given all the
 * components in the ParticleTileConfig the allocator holds, including the ones added
 * at runtime. Tiles that are copied or moved into the map are left as they are.
 *
 * A copy of a ParticleLevel, and a ParticleLevel it is assigned or swapped to, share
 * its configuration.
 */
template <class T, class Tile>
class ParticleLevelAllocator
{
public:
    using value_type = T;
    using propagate_on_container_copy_assignment = std::true_type;
    using propagate_on_container_move_assignment = std::true_type;
    using propagate_on_container_swap = std::true_type;
    //! the memory of every instance is interchangeable; only the config differs
    using is_always_equal = std::true_type;

    template <class U>
    struct rebind { using other = ParticleLevelAllocator<U, Tile>; };

    ParticleLevelAllocator () noexcept = default;

    explicit ParticleLevelAllocator (std::shared_ptr<ParticleTileConfig> a_config) noexcept
        : m_config(std::move(a_config))
    {}

    // Copy the config instead of moving it, so that a map that is moved from can
    // still create tiles.
    ParticleLevelAllocator (ParticleLevelAllocator const&) noexcept = default;
    ParticleLevelAllocator& operator= (ParticleLevelAllocator const&) noexcept = default;
    ~ParticleLevelAllocator () = default;

    template <class U>
    ParticleLevelAllocator (ParticleLevelAllocator<U, Tile> const& other) noexcept // NOLINT(google-explicit-constructor)
        : m_config(other.config())
    {}

    [[nodiscard]] T* allocate (std::size_t n) { return std::allocator<T>{}.allocate(n); }

    void deallocate (T* p, std::size_t n) noexcept { std::allocator<T>{}.deallocate(p, n); }

    template <class U, class... Args>
    void construct (U* p, Args&&... args)
    {
        ::new (static_cast<void*>(p)) U(std::forward<Args>(args)...);
        if constexpr (createsTile<U, Args...>()) {
            if (m_config) {
                try {
                    defineTile(p->second);
                } catch (...) {
                    p->~U();
                    throw;
                }
            }
        }
    }

    [[nodiscard]] std::shared_ptr<ParticleTileConfig> const& config () const noexcept
    {
        return m_config;
    }

    template <class U>
    friend bool operator== (ParticleLevelAllocator const&, ParticleLevelAllocator<U, Tile> const&) noexcept
    {
        return true;
    }

private:

    //! Whether constructing a U from Args creates a tile from nothing, as
    //! value_type(piecewise_construct, forward_as_tuple(key), tuple<>()) does.
    template <class U, class... Args>
    static constexpr bool createsTile ()
    {
        if constexpr (std::is_same_v<U, std::pair<const std::pair<int, int>, Tile>> &&
                      sizeof...(Args) == 3) {
            using ArgTuple = std::tuple<std::remove_cvref_t<Args>...>;
            return std::is_same_v<std::tuple_element_t<0, ArgTuple>, std::piecewise_construct_t> &&
                   std::is_same_v<std::tuple_element_t<2, ArgTuple>, std::tuple<>>;
        } else {
            return false;
        }
    }

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

/**
 * \brief The particle tiles of one AMR level, keyed by (grid index, local tile index).
 *
 * This is a std::map with a ParticleLevelAllocator. A tile it creates without being
 * given one, with operator[] or try_emplace(key), has all the components of its
 * ParticleContainer, including the ones added at runtime. For a ParticleContainer
 * with a PolymorphicArenaAllocator, the tile also uses the container's arena, if it
 * has been set. Tiles added with insert, or with emplace from a tile, are added as
 * they are given.
 *
 * Code that names this type should use ParticleContainer::ParticleLevel (or auto),
 * not std::map<std::pair<int,int>, Tile>, which is a different type.
 */
template <class Tile>
using ParticleLevelMap = std::map<std::pair<int, int>, Tile, std::less<std::pair<int, int>>,
                                  ParticleLevelAllocator<std::pair<const std::pair<int, int>, Tile>, Tile>>;

/**
 * \brief Make the tiles that pmap creates from now on be set up from a_config. The
 *        tiles pmap already has are kept as they are.
 */
template <class Tile>
void setParticleLevelConfig (ParticleLevelMap<Tile>& pmap,
                             std::shared_ptr<ParticleTileConfig> const& a_config)
{
    using Map = ParticleLevelMap<Tile>;
    if (pmap.get_allocator().config() != a_config) {
        // The allocators are equal, so this moves the nodes, not the tiles.
        pmap = Map(std::move(pmap), typename Map::allocator_type(a_config));
    }
}

}

#endif
