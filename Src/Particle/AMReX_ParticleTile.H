#ifndef AMREX_PARTICLETILE_H_
#define AMREX_PARTICLETILE_H_

#include <AMReX_Particle.H>
#include <AMReX_ArrayOfStructs.H>
#include <AMReX_StructOfArrays.H>
#include <AMReX_Vector.H>
#include <AMReX_IndexSequence.H>

#include <tuple>
#include <array>
#include <cassert>

#ifdef AMREX_USE_CUDA
#include <thrust/tuple.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/device_ptr.h>
#endif

namespace amrex {

namespace detail {
    
    template <size_t I, typename T>
    struct tuple_n{
        template< typename...Args> using type = typename tuple_n<I-1, T>::template type<T, Args...>;
    };
    
    template <typename T>
    struct tuple_n<0, T> {
        template<typename...Args> using type = std::tuple<Args...>;
    };
    
    template <size_t I, typename T> using tuple_of = typename tuple_n<I,T>::template type<>;
    
    template<size_t N, typename T, std::size_t... Is>
    tuple_of<N, T*> make_pointer_tuple(std::array<T, N>& arr, IndexSequence<Is...>)
    {
        return std::make_tuple(&arr[Is]...);
    }
    
    template<size_t N, typename T>
    tuple_of<N, T*> make_pointer_tuple(std::array<T, N>& arr)
    {
        return make_pointer_tuple(arr, makeIndexSequence<N>());
    }    
    
#ifdef AMREX_USE_CUDA
    template<typename... Args, std::size_t... Is>
    auto to_thrust_tuple(const std::tuple<Args...>& t, IndexSequence<Is...>)
        -> thrust::tuple<typename thrust::device_ptr<typename std::remove_pointer<Args>::type::value_type>...>
    {
        return thrust::make_tuple(thrust::device_ptr<typename std::remove_pointer<Args>::type::value_type>(std::get<Is>(t)->data())...);
    }
#endif // AMREX_USE_CUDA
}; // namespace detail

#ifdef AMREX_USE_CUDA
template<typename... Args>
auto to_thrust_tuple(const std::tuple<Args...>& t)
    -> thrust::tuple<typename thrust::device_ptr<typename std::remove_pointer<Args>::type::value_type>...>
{
    return detail::to_thrust_tuple(t, makeIndexSequence<sizeof...(Args)>());
}
#endif // AMREX_USE_CUDA

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
struct ParticleTile
{
    using Real = double;
    using ParticleType = Particle<NStructReal, NStructInt>;

    using AoS = ArrayOfStructs<NStructReal, NStructInt>;
    using ParticleVector = typename AoS::ParticleVector;

    using SoA = StructOfArrays<NArrayReal, NArrayInt>;
    using RealVector = typename SoA::RealVector;
    using IntVector = typename SoA::IntVector;

#ifdef AMREX_USE_CUDA
    using struct_tuple_t = detail::tuple_of<1, ParticleVector*>;
    using real_tuple_t = detail::tuple_of<NArrayReal, RealVector*>;
    using int_tuple_t  = detail::tuple_of<NArrayInt,  IntVector* >;
    using tuple_t = decltype(std::tuple_cat(struct_tuple_t(), real_tuple_t(), int_tuple_t()));
    using thrust_zip_iterator_t = decltype(thrust::make_zip_iterator(to_thrust_tuple(tuple_t())));
#endif

    AoS&       GetArrayOfStructs ()       { return m_aos_tile; }
    const AoS& GetArrayOfStructs () const { return m_aos_tile; }

    SoA&       GetStructOfArrays ()       { return m_soa_tile; }
    const SoA& GetStructOfArrays () const { return m_soa_tile; }

    bool empty () const { return m_aos_tile.empty(); }
    
    std::size_t size () const { return m_aos_tile.size(); }

    int numParticles () const { return m_aos_tile.numParticles(); }

    ///
    /// Add one particle to this tile.
    ///
    void push_back (const ParticleType& p) { m_aos_tile().push_back(p); }

    ///
    /// Add a Real value to the struct-of-arrays at index comp.
    /// This sets the data for one particle.
    ///
    void push_back_real (int comp, Real v) { 
        m_soa_tile.GetRealData(comp).push_back(v);
    }

    ///
    /// Add Real values to the struct-of-arrays, for all comps at once.
    /// This sets the data for one particle.
    ///
    void push_back_real (const std::array<Real, NArrayReal>& v) { 
        for (int i = 0; i < NArrayReal; ++i) {
            m_soa_tile.GetRealData(i).push_back(v[i]);
        }
    }

    ///
    /// Add a range of Real values to the struct-of-arrays for the given comp.
    /// This sets the data for several particles at once.
    ///
    void push_back_real (int comp, const Real* beg, const Real* end) {
        auto it = m_soa_tile.GetRealData(comp).end();
        m_soa_tile.GetRealData(comp).insert(it, beg, end);
    }

    ///
    /// Add npar copies of the Real value v to the struct-of-arrays for the given comp.
    /// This sets the data for several particles at once.
    ///
    void push_back_real (int comp, std::size_t npar, Real v) {
        auto new_size = m_soa_tile.GetRealData(comp).size() + npar;
        m_soa_tile.GetRealData(comp).resize(new_size, v);
    }

    ///
    /// Add an int value to the struct-of-arrays at index comp.
    /// This sets the data for one particle.
    ///
    void push_back_int (int comp, int v) { 
        m_soa_tile.GetIntData(comp).push_back(v);
    }
    
    ///
    /// Add int values to the struct-of-arrays, for all comps at once.
    /// This sets the data for one particle.
    ///
    void push_back_int (const std::array<int, NArrayInt>& v) { 
        for (int i = 0; i < NArrayInt; ++i) {
            m_soa_tile.GetIntData(i).push_back(v[i]);
        }
    }
    
    ///
    /// Add a range of int values to the struct-of-arrays for the given comp.
    /// This sets the data for several particles at once.
    ///
    void push_back_int (int comp, const int* beg, const int* end) {
        auto it = m_soa_tile.GetIntData(comp).end();
        m_soa_tile.GetIntData(comp).insert(it, beg, end);
    }
    
    ///
    /// Add npar copies of the int value v to the struct-of-arrays for the given comp.
    /// This sets the data for several particles at once.
    ///
    void push_back_int (int comp, std::size_t npar, int v) {
        auto new_size = m_soa_tile.GetIntData(comp).size() + npar;
        m_soa_tile.GetIntData(comp).resize(new_size, v);
    }

#ifdef AMREX_USE_CUDA
    tuple_t data ()
    {
        auto struct_tuple = std::make_tuple(&m_aos_tile());
        auto real_tuple = detail::make_pointer_tuple(m_soa_tile.GetRealData());
        auto int_tuple = detail::make_pointer_tuple(m_soa_tile.GetIntData());
        return std::tuple_cat(struct_tuple, real_tuple, int_tuple);
    }
    
    thrust_zip_iterator_t begin ()
    {
        return thrust::make_zip_iterator(to_thrust_tuple(data()));
    }
    
    thrust_zip_iterator_t end ()
    {
        auto it = begin();
        thrust::advance(it, size());
        return it;
    }
#endif // AMREX_USE_CUDA

private:

    AoS m_aos_tile;
    SoA m_soa_tile;
};

} // namespace amrex;

#endif // AMREX_PARTICLETILE_H_
