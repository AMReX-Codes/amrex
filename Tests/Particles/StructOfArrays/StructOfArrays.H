#include <tuple>
#include <iostream>
#include <array>
#include <vector>
#include <cassert>

#include <thrust/tuple.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/device_ptr.h>

#include <CudaManagedAllocator.H>

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
    
    template<int... Is>
    struct seq { };
    
    template<int N, int... Is>
    struct gen_seq : gen_seq<N - 1, N - 1, Is...> { };
    
    template<int... Is>
    struct gen_seq<0, Is...> : seq<Is...> { };
    
    template<size_t N, typename T, int... Is>
    tuple_of<N, T*> make_pointer_tuple(std::array<T, N>& arr, seq<Is...>)
    {
        return std::make_tuple(&arr[Is]...);
    }
    
    template<size_t N, typename T>
    tuple_of<N, T*> make_pointer_tuple(std::array<T, N>& arr)
    {
        return make_pointer_tuple(arr, gen_seq<N>());
    }    
    
    template<typename... Args, int... Is>
    auto to_thrust_tuple(const std::tuple<Args...>& t, seq<Is...>)
        -> thrust::tuple<typename thrust::device_ptr<typename std::remove_pointer<Args>::type::value_type>...>
    {
        return thrust::make_tuple(thrust::device_ptr<typename std::remove_pointer<Args>::type::value_type>(std::get<Is>(t)->data())...);
    }

    template<typename TupleType, typename IteratorType, typename FunctionType>
    void for_each(TupleType&&, IteratorType&&, IteratorType&&, FunctionType
                  , std::integral_constant<size_t, std::tuple_size<typename std::remove_reference<TupleType>::type >::value>) {}
    
    template<std::size_t I, typename TupleType, typename IteratorType, typename FunctionType
             , typename = typename std::enable_if<I!=std::tuple_size<typename std::remove_reference<TupleType>::type>::value>::type >
    void for_each(TupleType&& t, IteratorType&& first, IteratorType&& last,
                  FunctionType f, std::integral_constant<size_t, I>)
    {
        f(std::get<I>(t), thrust::get<I>(first), thrust::get<I>(last));
        for_each(std::forward<TupleType>(t), 
                 std::forward<IteratorType>(first), 
                 std::forward<IteratorType>(last), 
                 f, std::integral_constant<size_t, I + 1>());
    }    
};

template<typename TupleType, typename IteratorType, typename FunctionType>
void for_each(TupleType&& t, IteratorType&& first, IteratorType&& last, FunctionType f)
{
    detail::for_each(std::forward<TupleType>(t), 
                     std::forward<IteratorType>(first), 
                     std::forward<IteratorType>(last), 
                     f, std::integral_constant<size_t, 0>());
}

template<typename... Args>
auto to_thrust_tuple(const std::tuple<Args...>& t)
    -> thrust::tuple<typename thrust::device_ptr<typename std::remove_pointer<Args>::type::value_type>...>
{
    return detail::to_thrust_tuple(t, detail::gen_seq<sizeof...(Args)>());
}

template <int NReal, int NInt>
struct StructOfArrays {
    
    using RealVector = std::vector<double, CudaManagedAllocator<double> >;
    using IntVector  = std::vector<int,    CudaManagedAllocator<int   > >;
    
    using real_tuple_t = detail::tuple_of<NReal, RealVector*>;
    using int_tuple_t  = detail::tuple_of<NInt,  IntVector* >;
    using tuple_t = decltype(std::tuple_cat(real_tuple_t(), int_tuple_t()));
    using thrust_zip_iterator_t = decltype(thrust::make_zip_iterator(to_thrust_tuple(tuple_t())));

    std::array<RealVector, NReal>& GetRealData () { return m_rdata; }
    std::array<IntVector,  NInt>&  GetIntData  () { return m_idata; }
    
    const std::array<RealVector, NReal>& GetRealData () const { return m_rdata; }
    const std::array<IntVector,  NInt >& GetIntData  () const { return m_idata; }
    
    RealVector& GetRealData (const int index)
    {
        assert(index < NReal);
        return m_rdata[index];
    }
    
    const RealVector& GetRealData (const int index) const 
    {
        assert(index < NReal);
        return m_rdata[index];
    }
    
    IntVector& GetIntData (const int index) 
    {
        assert(index < NInt);
        return m_idata[index];
    }
    
    const IntVector& GetIntData (const int index) const 
    {
        assert(index < NInt);
        return m_idata[index];
    }
    
    int size () const
    {
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
        auto it = thrust::make_zip_iterator(to_thrust_tuple(data()));
        thrust::advance(it, size());
        return it;
    }
    
    struct vector_insert_functor
    {
        size_t offset;
        
        vector_insert_functor(size_t pos_index) : offset(pos_index) {} 
        
        template<typename DstIt, typename SrcIt>
        void operator () (DstIt const &v, SrcIt const& first, SrcIt const& last)
        {
            v->insert(v->begin() + offset, &first, &last);
        }
    };
    
    void insert(thrust_zip_iterator_t pos, thrust_zip_iterator_t first, thrust_zip_iterator_t last) 
    {
        const size_t pos_index = thrust::distance(begin(), pos);
        for_each(data(), first[0], last[0], vector_insert_functor(pos_index));
    }
    
private:
    
    std::array<RealVector, NReal> m_rdata;
    std::array<IntVector,  NInt>  m_idata;
};
