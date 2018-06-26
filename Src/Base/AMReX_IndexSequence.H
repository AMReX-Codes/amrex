#ifndef AMREX_INDEX_SEQUENCE_H_
#define AMREX_INDEX_SEQUENCE_H_

#include <utility>

namespace amrex {

#if __cplusplus >= 201402L

    template <std::size_t... ISs>
    using IndexSequence = std::integer_sequence<std::size_t, ISs...>;

    template <std::size_t N>
    using makeIndexSequence = std::make_index_sequence<N>;

#else

    template <std::size_t... ISs> struct IndexSequence {};

    template <std::size_t N, std::size_t... ISs>
    struct makeIndexSequenceHelper
    {
        using type = typename makeIndexSequenceHelper<N-1, N-1, ISs...>::type;
    };

    template <std::size_t... ISs>
    struct makeIndexSequenceHelper<0, ISs...>
    {
        using type = IndexSequence<ISs...>;
    };

    template <std::size_t N>
    using makeIndexSequence = typename makeIndexSequenceHelper<N>::type;

#endif

}

#endif
