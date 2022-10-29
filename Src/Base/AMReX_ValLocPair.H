#ifndef AMREX_VALLOCPAIR_H_
#define AMREX_VALLOCPAIR_H_

#include <limits>

namespace amrex {

template <typename TV, typename TI>
struct ValLocPair
{
    TV value;
    TI index;

    static constexpr ValLocPair<TV,TI> max () {
        return ValLocPair<TV,TI>{std::numeric_limits<TV>::max(), TI()};
    }

    static constexpr ValLocPair<TV,TI> lowest () {
        return ValLocPair<TV,TI>{std::numeric_limits<TV>::lowest(), TI()};
    }

    friend constexpr bool operator< (ValLocPair<TV,TI> const& a, ValLocPair<TV,TI> const& b)
    {
        return a.value < b.value;
    }

    friend constexpr bool operator> (ValLocPair<TV,TI> const& a, ValLocPair<TV,TI> const& b)
    {
        return a.value > b.value;
    }
};

}

#endif
