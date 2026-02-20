#ifndef BL_MAKETYPE_H_
#define BL_MAKETYPE_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_MakeType.H
 * \brief Enumeration controlling alias vs. deep-copy construction.
 */

namespace amrex
{
    /// \ingroup amrex_utilities
    //! Controls whether an operation creates an alias or a deep copy.
    enum MakeType {make_alias, make_deep_copy};
}

#endif

