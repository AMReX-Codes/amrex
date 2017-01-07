
#ifndef BL_TUPLE_H
#define BL_TUPLE_H

#include <array>

namespace amrex {
    template <class T, std::size_t N>
    using Tuple = std::array<T,N>;
}

#endif /*BL_TUPLE_H*/
