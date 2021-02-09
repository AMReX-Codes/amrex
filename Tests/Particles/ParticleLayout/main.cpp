#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_PODVector.H>
#include <AMReX_Tuple.H>
#include <AMReX_GpuContainers.H>

#include <functional>
#include <tuple>

/**
   /brief Particle type is defined as a GpuTuple.
   x, y, z, id, and cpu methods are required,
   but can be implemented freely, i.e. using
   a cell + offset representation of position.
 */
template <typename... T>
struct Particle : public amrex::GpuTuple<T...>
{
    using amrex::GpuTuple<T...>::GpuTuple;
    auto& x () { return amrex::get<0>(*this); }
    auto& y () { return amrex::get<1>(*this); }
    auto& z () { return amrex::get<2>(*this); }
    auto& id () { return amrex::get<3>(*this); }
    auto& cpu () { return amrex::get<4>(*this); }
};

/**
   A tag that defines the data layout policy used by
   particle tiles.
 */
enum class DataLayout
{
    AoS = 0,
    SoA
};

template <typename ContainerType>
class ParticleIterator;

template <template <typename...> class ContainerType,
          typename ParticleType,
          DataLayout DataLayoutTag>
struct DataLayoutPolicy;

template <typename ParticleType, DataLayout DataLayoutTag>
struct DataLayoutPolicyRaw;

template <typename ParticleType, DataLayout DataLayoutTag>
struct ParticleTileRaw;

template <typename T>
struct ref_wrapper : public std::reference_wrapper<T>
{
    operator T& () const noexcept { return this->get(); }
    ref_wrapper (T& a_other) : std::reference_wrapper<T>(a_other) {}
    void operator = (T&& a_other) {this->get()=a_other;}
};

/**
   Implementation of the AoS policy. Pretty much a
   straightforward wrapper around ContainterType<ParticleType>
 */
template <template <typename...> class ContainerType,
          template<typename...> class ParticleType,
          typename... Types>
struct DataLayoutPolicy<ContainerType, ParticleType<Types...>, DataLayout::AoS>
{
    using container_type = ContainerType<ParticleType<Types...>>;
    using value_type = ParticleType<Types...>&;
    using pointer_type = ParticleType<Types...>*;

    constexpr static pointer_type get_data_pointers (container_type& a_container)
    {
        return pointer_type(static_cast<ParticleType<Types...>*>(&a_container[0]));
    }

    constexpr static void resize (container_type& a_container, std::size_t a_size)
    {
        a_container.resize(a_size);
    }

    template <typename ValueType>
    constexpr static void push_back (container_type& a_container, ValueType&& a_value)
    {
        a_container.push_back(a_value);
    }

    static constexpr std::size_t size (container_type& a_container)
    {
        return a_container.size();
    }
};

/**
   A non-owning verion of AoS policy for passing to the GPU.
 */
template <template<typename...> class ParticleType, typename... Types>
struct DataLayoutPolicyRaw<ParticleType<Types...>, DataLayout::AoS>
{
    using pointer_type = ParticleType<Types...>*;
    using value_type = ParticleType<Types...>&;

    constexpr static value_type get (pointer_type a_ptr, std::size_t a_index)
    {
        return value_type(*static_cast<ParticleType<Types...>*>(&a_ptr[a_index]));
    }
};

/**
   Implementation of the SoA policy. The underlying data structure
   is a Tuple<ContainerType<ParticleType>>. Note that unlike the AoS,
   this container works with a "ref_wrap"ed version of the particle data,
   so we can modify the particle data in the tile.
 */
template <template <typename...> class ContainerType,
          template<typename...> class ParticleType,
          typename... Types>
struct DataLayoutPolicy<ContainerType, ParticleType<Types...>, DataLayout::SoA> {
    using container_type = std::tuple<ContainerType<Types>...>;
    using value_type = ParticleType<ref_wrapper<Types>...>;
    using pointer_type = amrex::GpuTuple<Types*...>;

    constexpr static pointer_type get_data_pointers (container_type& a_container)
    {
        return get_data_pointers_impl(a_container, amrex::makeIndexSequence<sizeof...(Types)>());
    }

    constexpr static void resize (container_type& a_container, std::size_t a_size)
    {
        resize_impl(a_container, a_size, amrex::makeIndexSequence<sizeof...(Types)>());
    }

    template <typename ValueType>
    constexpr static void push_back (container_type& a_container, ValueType&& a_value)
    {
        push_back_impl(a_container, std::forward<ValueType>(a_value),
                       amrex::makeIndexSequence<sizeof...(Types)>());
    }

    static constexpr std::size_t size (container_type& a_container)
    {
        return std::get<0>(a_container).size();
    }

private:

    template <std::size_t... Is>
    constexpr static auto get_data_pointers_impl (container_type& a_container,
                                                  amrex::IndexSequence<Is...>)
    {
        return pointer_type{static_cast<Types*>(&std::get<Is>(a_container)[0])... };
    }

    template <std::size_t... Is>
    constexpr static void resize_impl (container_type& a_container, std::size_t a_size,
                                       amrex::IndexSequence<Is...>)
    {
        using expander = int[];
        (void) expander{ 0, (std::get<Is>(a_container).resize(a_size), 0)... };
    }

    template <typename ValueType, std::size_t... Is>
    constexpr static void push_back_impl(container_type& a_container, ValueType&& a_value,
                                         amrex::IndexSequence<Is...>)
    {
        using expander = int[];
        (void) expander{ 0, (std::get<Is>(a_container).push_back(
                                 std::get<Is>(std::forward<ValueType>(a_value))), 0)... };
    }
};

/**
   A non-owning verion of SoA policy for passing to the GPU.
 */
template <template<typename...> class ParticleType, typename... Types>
struct DataLayoutPolicyRaw<ParticleType<Types...>, DataLayout::SoA> {
    using pointer_type = amrex::GpuTuple<Types*...>;
    using value_type = ParticleType<ref_wrapper<Types>...>;

    constexpr static value_type get (pointer_type& a_tuple, std::size_t a_index)
    {
        return get_impl(a_tuple, a_index, amrex::makeIndexSequence<sizeof...(Types)>());
    }

private:

    template <std::size_t... Is>
    constexpr static auto get_impl (pointer_type& a_tuple, std::size_t a_index,
                                    amrex::IndexSequence<Is...>)
    {
        return value_type{ref_wrapper<Types>(amrex::get<Is>(a_tuple)[a_index])... };
    }
};

/**
   Tile implementation, it basically just forwards to the policy's methods.
 */
template <template <typename ValueType> class ContainerType,
          typename ParticleType,
          DataLayout DataLayoutTag>
struct ParticleTile
{
    using policy_type = DataLayoutPolicy<ContainerType, ParticleType, DataLayoutTag>;
    using raw_type = ParticleTileRaw<ParticleType, DataLayoutTag>;
    using value_type = typename policy_type::value_type;
    using container_type = typename policy_type::container_type;

    auto static constexpr data_layout = DataLayoutTag;

    ParticleTile ()
    {
        resize(0);
    }

    ParticleTile (size_t a_size)
    {
        resize(a_size);
    }

    template <typename ValueType>
    void push_back (ValueType&& val)
    {
        policy_type::push_back(m_data, std::forward<ValueType>(val));
    }

    std::size_t size ()
    {
        return policy_type::size(m_data);
    }

    void resize (size_t a_size)
    {
        policy_type::resize(m_data, a_size);
    }

    raw_type get_particle_data ()
    {
        return raw_type(size(), policy_type::get_data_pointers(m_data));
    }

private:

    container_type m_data;
};

/**
   A Version of ParticleTile that contains only raw pointers, so that it can
   be copied by value into GPU kernels.
 */
template <typename ParticleType, DataLayout DataLayoutTag>
struct ParticleTileRaw
{
    using policy_type  = DataLayoutPolicyRaw<ParticleType, DataLayoutTag>;
    using iterator     = ParticleIterator<ParticleTileRaw<ParticleType, DataLayoutTag>>;
    using value_type   = typename policy_type::value_type;
    using pointer_type = typename policy_type::pointer_type;

    auto static constexpr data_layout = DataLayoutTag;

    ParticleTileRaw (std::size_t a_size, pointer_type a_data)
        : m_size(a_size), m_data(a_data)
    {}

    value_type operator[] (std::size_t a_index) { return policy_type::get(m_data, a_index); }

    std::size_t size () { return m_size; }

    iterator begin() { return iterator(this, 0); }
    iterator end()   { return iterator(this, m_size); }

private:

    std::size_t m_size;
    pointer_type m_data;
};

/**
   Iterator over the particles in a tile.

 */
template <typename ParticleTileType>
class ParticleIterator
{
private:
    using ptile_type = ParticleTileType;

public:
    using policy_type = typename ptile_type::policy_type;
    using value_type  = typename policy_type::value_type;

    template<typename U>
    ParticleIterator (U* a_ptile, std::size_t a_index = 0):
        m_ptile(a_ptile),
        m_index(a_index)
    {}

    ParticleIterator& operator= (ParticleIterator const& a_other )
    {
        m_index = a_other.m_index;
    }

    friend bool operator!= (ParticleIterator const& lhs, ParticleIterator const& rhs)
    {
        return lhs.m_index != rhs.m_index;
    }

    friend bool operator== (ParticleIterator const& lhs, ParticleIterator const& rhs)
    {
        return !operator!=(lhs, rhs);
    }

    template <typename T>
    void operator+= (T a_offset) { m_index += a_offset; }

    template <typename T>
    void operator-= (T a_offset) { m_index -= a_offset; }

    void operator++ () { return operator +=(1); }
    void operator-- () { return operator -=(1); }

    value_type operator* () { return (*m_ptile)[m_index]; }

private:

    ptile_type* m_ptile = nullptr;
    std::size_t m_index = std::numeric_limits<std::size_t>::infinity();
};

template <DataLayout DataLayoutTag>
void testLayout ()
{
    using ParticleType = Particle<double, double, double, int, int>;
    using ParticleTileType = ParticleTile<amrex::Gpu::DeviceVector, ParticleType, DataLayoutTag>;

    ParticleTileType ptile;
    ptile.resize(1000);

    auto particles = ptile.get_particle_data();

    amrex::ParallelFor( particles.size(), [=] AMREX_GPU_DEVICE (int i) noexcept
    {
	particles[i].x() = 7.0;
        particles[i].y() = 8.0;
	particles[i].z() = 9.0;
    });

    amrex::ParallelFor( particles.size(), [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        auto&& p = particles[i];
        p.x() -= 6.0;
        p.y() -= 7.0;
        p.z() -= 8.0;
    });

    amrex::ParallelFor( particles.size(), [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        auto&& p = particles[i];
        AMREX_ALWAYS_ASSERT(p.x() == p.y() == p.z() == 1.0);
    });
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    amrex::Print() << "Running data layout test \n";
    testLayout<DataLayout::AoS>();
    testLayout<DataLayout::SoA>();

    amrex::Finalize();
}
