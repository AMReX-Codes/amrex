#ifndef AMREX_STRUCTOFARRAYS_H_
#define AMREX_STRUCTOFARRAYS_H_

#include <array>
#include <tuple>

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_CudaAllocators.H>
#include <AMReX_IndexSequence.H>

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
};

#ifdef AMREX_USE_CUDA
template<typename... Args>
auto to_thrust_tuple(const std::tuple<Args...>& t)
    -> thrust::tuple<typename thrust::device_ptr<typename std::remove_pointer<Args>::type::value_type>...>
{
    return detail::to_thrust_tuple(t, makeIndexSequence<sizeof...(Args)>());
}
#endif // AMREX_USE_CUDA

template <int NReal, int NInt>
struct StructOfArrays {
  
#ifdef AMREX_USE_CUDA
    using RealVector = Vector<Real, CudaManagedAllocator<Real>>;
    using  IntVector = Vector<int , CudaManagedAllocator<int >>;
#else
    using RealVector = Vector<Real>;
    using  IntVector = Vector<int >;
#endif

#ifdef AMREX_USE_CUDA
    using real_tuple_t = detail::tuple_of<NReal, RealVector*>;
    using int_tuple_t  = detail::tuple_of<NInt,  IntVector* >;
    using tuple_t = decltype(std::tuple_cat(real_tuple_t(), int_tuple_t()));
    using thrust_zip_iterator_t = decltype(thrust::make_zip_iterator(to_thrust_tuple(tuple_t())));
#endif

    std::array<RealVector, NReal>& GetRealData () { return m_rdata; }
    std::array< IntVector,  NInt>& GetIntData  () { return m_idata; }

    const std::array<RealVector, NReal>& GetRealData () const { return m_rdata; }
    const std::array< IntVector,  NInt>& GetIntData  () const { return m_idata; }

    RealVector& GetRealData (const int index) {
        AMREX_ASSERT(index < NReal);
        return m_rdata[index];
    }
    
    const RealVector& GetRealData (const int index) const {
        AMREX_ASSERT(index < NReal);
        return m_rdata[index];
    }

    IntVector& GetIntData (const int index) {
        AMREX_ASSERT(index < NInt);
        return m_idata[index];
    }
    
    const IntVector& GetIntData (const int index) const {
        AMREX_ASSERT(index < NInt);
        return m_idata[index];
    }

    int numParticles () const { 
        if (NReal > 0)
            return m_rdata[0].size();
        else if (NInt > 0) 
            return m_idata[0].size();
        else
            return 0;
    }

    void resize(size_t count)
    {
        for (int i = 0; i < NReal; ++i) m_rdata[i].resize(count);
        for (int i = 0; i < NInt;  ++i) m_idata[i].resize(count); 
    }

#ifdef AMREX_USE_CUDA    
    tuple_t data()
    {
        auto real_tuple = detail::make_pointer_tuple(m_rdata);
        auto int_tuple = detail::make_pointer_tuple(m_idata);
        return std::tuple_cat(real_tuple, int_tuple);
    }
    
    thrust_zip_iterator_t begin()
    {
        return thrust::make_zip_iterator(to_thrust_tuple(data()));
    }
    
    thrust_zip_iterator_t end()
    {
        auto it = begin();
        thrust::advance(it, numParticles());
        return it;
    }
#endif

private:
    std::array<RealVector, NReal> m_rdata;
    std::array< IntVector,  NInt> m_idata;    
};

} // namespace amrex

#endif // AMREX_STRUCTOFARRAYS_H_
