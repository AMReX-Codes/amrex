#ifndef AMREX_PARTICLEARRAY_H_
#define AMREX_PARTICLEARRAY_H_
#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_Tuple.H>

#include <functional>
#include <tuple>
#include <type_traits>

namespace amrex
{
/**
   A tag that defines the data layout policy used by
   particle tiles.
 */
enum class DataLayout
{
    AoS = 0,
    SoA
};

/**
   Forward declarations
 */
template <template <typename...> class ContainerType,
          typename ParticleType,
          DataLayout DataLayoutTag>
struct DataLayoutPolicy;

template <typename ParticleType, DataLayout DataLayoutTag>
struct DataLayoutPolicyRaw;

template <typename ParticleType, DataLayout DataLayoutTag>
struct ParticleArrayAccessor;

/**
   Because std::reference_wrapper is not __host__ __device__
 */
template <class T>
class ref_wrapper {
public:
    typedef T type;

    AMREX_GPU_HOST_DEVICE
    ref_wrapper(T& ref) noexcept : _ptr(&ref) {}
    ref_wrapper(T&&) = delete;

    ~ref_wrapper () = default;
    ref_wrapper (const ref_wrapper&) noexcept = default;
    ref_wrapper (ref_wrapper&&) noexcept = default;

    AMREX_GPU_HOST_DEVICE
    ref_wrapper& operator= (T&& a_other) { this->get()=a_other; return *this; }

    ref_wrapper& operator= (const ref_wrapper&) noexcept = default;
    ref_wrapper& operator= (ref_wrapper&&) noexcept = default;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    operator T& () const noexcept { return *_ptr; }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    T& get() const noexcept { return *_ptr; }

private:
    T* _ptr;
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
    using raw_type = ParticleType<Types...>*;
    using value_type = ParticleType<Types...>&;

    static constexpr raw_type get_raw_data (container_type& a_container)
    {
        return raw_type(static_cast<ParticleType<Types...>*>(&a_container[0]));
    }

    static constexpr void resize (container_type& a_container, std::size_t a_size)
    {
        a_container.resize(a_size);
    }

    template <typename ValueType>
    static constexpr void push_back (container_type& a_container, ValueType&& a_value)
    {
        a_container.push_back(a_value);
    }

    static constexpr std::size_t size (container_type& a_container)
    {
        return a_container.size();
    }
};

/**
   A non-owning version of AoS policy for passing to the GPU.
 */
template <template<typename...> class ParticleType, typename... Types>
struct DataLayoutPolicyRaw<ParticleType<Types...>, DataLayout::AoS>
{
    using raw_type = ParticleType<Types...>*;
    using value_type = ParticleType<Types...>&;

    AMREX_GPU_HOST_DEVICE
    static constexpr value_type get (raw_type a_ptr, std::size_t a_index)
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
    using raw_type = amrex::GpuTuple<Types*...>;
    using value_type = ParticleType<ref_wrapper<Types>...>;

    static constexpr raw_type get_raw_data (container_type& a_container)
    {
        return get_raw_data_impl(a_container, std::make_index_sequence<sizeof...(Types)>());
    }

    static constexpr void resize (container_type& a_container, std::size_t a_size)
    {
        resize_impl(a_container, a_size, std::make_index_sequence<sizeof...(Types)>());
    }

    template <typename ValueType>
    static constexpr void push_back (container_type& a_container, ValueType&& a_value)
    {
        push_back_impl(a_container, std::forward<ValueType>(a_value),
                       std::make_index_sequence<sizeof...(Types)>());
    }

    static constexpr std::size_t size (container_type& a_container)
    {
        return std::get<0>(a_container).size();
    }

private:

    template <std::size_t... Is>
    static constexpr auto get_raw_data_impl (container_type& a_container,
                                             std::index_sequence<Is...>)
    {
        return raw_type{static_cast<Types*>(&std::get<Is>(a_container)[0])... };
    }

    template <std::size_t... Is>
    static constexpr void resize_impl (container_type& a_container, std::size_t a_size,
                                       std::index_sequence<Is...>)
    {
        using expander = int[];
        (void) expander{ 0, (std::get<Is>(a_container).resize(a_size), 0)... };
    }

    template <typename ValueType, std::size_t... Is>
    static constexpr void push_back_impl(container_type& a_container, ValueType&& a_value,
                                         std::index_sequence<Is...>)
    {
        using expander = int[];
        (void) expander{ 0, (std::get<Is>(a_container).push_back(
                                 std::get<Is>(std::forward<ValueType>(a_value))), 0)... };
    }
};

/**
   A non-owning version of SoA policy for passing to the GPU.
 */
template <template<typename...> class ParticleType, typename... Types>
struct DataLayoutPolicyRaw<ParticleType<Types...>, DataLayout::SoA> {
    using raw_type = amrex::GpuTuple<Types*...>;
    using value_type = ParticleType<ref_wrapper<Types>...>;

    AMREX_GPU_HOST_DEVICE
    static constexpr value_type get (raw_type const& a_tuple, std::size_t a_index)
    {
        return get_impl(a_tuple, a_index, std::make_index_sequence<sizeof...(Types)>());
    }

private:

    template <std::size_t... Is>
    AMREX_GPU_HOST_DEVICE
    static constexpr auto get_impl (raw_type const& a_tuple, std::size_t a_index,
                                    std::index_sequence<Is...>)
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
struct ParticleArray
{
    using policy_type = DataLayoutPolicy<ContainerType, ParticleType, DataLayoutTag>;
    using accessor_type = ParticleArrayAccessor<ParticleType, DataLayoutTag>;
    using value_type = typename policy_type::value_type;
    using container_type = typename policy_type::container_type;

    static constexpr auto data_layout = DataLayoutTag;

    ParticleArray () { resize(0); }

    ParticleArray (size_t a_size) { resize(a_size); }

    template <typename ValueType>
    void push_back (ValueType&& val)
    {
        policy_type::push_back(m_data, std::forward<ValueType>(val));
    }

    std::size_t size () { return policy_type::size(m_data); }

    void resize (size_t a_size) { policy_type::resize(m_data, a_size); }

    accessor_type get_particle_data ()
    {
        return accessor_type(size(), policy_type::get_raw_data(m_data));
    }

private:

    container_type m_data;
};

/**
   A Version of ParticleArray that contains only raw pointers, so that it can
   be copied by value into GPU kernels.
 */
template <typename ParticleType, DataLayout DataLayoutTag>
struct ParticleArrayAccessor
{
    // ParticleType: Particle<double, double, double, int, int>
    // DataLayoutTag: amrex::DataLayout::SoA
    // value_type: Particle<amrex::ref_wrapper<double>, amrex::ref_wrapper<double>, amrex::ref_wrapper<double>, amrex::ref_wrapper<int>, amrex::ref_wrapper<int> >

    using policy_type  = DataLayoutPolicyRaw<ParticleType, DataLayoutTag>;
    using value_type   = typename policy_type::value_type;
    using raw_type     = typename policy_type::raw_type;

    static constexpr auto  data_layout = DataLayoutTag;

    ParticleArrayAccessor (std::size_t a_size, raw_type a_data)
        : m_size(a_size), m_data(a_data)
    {}

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    value_type operator[] (std::size_t a_index) const { return policy_type::get(m_data, a_index); }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    std::size_t size () const { return m_size; }

private:
    std::size_t m_size;
    raw_type m_data;
};
}

#endif
