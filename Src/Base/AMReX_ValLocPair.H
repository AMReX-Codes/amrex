#ifndef AMREX_VALLOCPAIR_H_
#define AMREX_VALLOCPAIR_H_

/**
 * \file AMReX_ValLocPair.H
 * \brief Value-location pair type used in reductions.
 */

#include <limits>

namespace amrex {

/**
 * \class amrex::ValLocPair
 * \ingroup amrex_utilities
 * \brief A pair holding a value and an associated index, used in
 *        value-location reductions (e.g., min/max with location).
 *
 * \tparam TV Value type.
 * \tparam TI Index (location) type.
 */
template <typename TV, typename TI>
struct ValLocPair
{
    TV value; ///< The value.
    TI index; ///< The associated index (location).

    //! Return a mutable reference to the value (first element).
    constexpr TV&       first  ()       { return value; }
    //! Return a const reference to the value (first element).
    constexpr TV const& first  () const { return value; }
    //! Return a mutable reference to the index (second element).
    constexpr TI&       second ()       { return index; }
    //! Return a const reference to the index (second element).
    constexpr TI const& second () const { return index; }

    //! Return a ValLocPair with the maximum representable value and a default index.
    static constexpr ValLocPair<TV,TI> max () {
        return ValLocPair<TV,TI>{std::numeric_limits<TV>::max(), TI()};
    }

    //! Return a ValLocPair with the lowest representable value and a default index.
    static constexpr ValLocPair<TV,TI> lowest () {
        return ValLocPair<TV,TI>{std::numeric_limits<TV>::lowest(), TI()};
    }

    //! Return true if the value of \p a is less than the value of \p b.
    friend constexpr bool operator< (ValLocPair<TV,TI> const& a, ValLocPair<TV,TI> const& b)
    {
        return a.value < b.value;
    }

    //! Return true if the value of \p a is greater than the value of \p b.
    friend constexpr bool operator> (ValLocPair<TV,TI> const& a, ValLocPair<TV,TI> const& b)
    {
        return a.value > b.value;
    }
};

/// \ingroup amrex_utilities
//! Alias for ValLocPair using key/value terminology.
template <typename K, typename V> using KeyValuePair = ValLocPair<K,V>;

}

#endif
