#ifndef AMREX_EXCEPTION_H_
#define AMREX_EXCEPTION_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_Exception.H
 * \brief AMReX exception type.
 */

#include <stdexcept>

namespace amrex
{
    /// \ingroup amrex_utilities
    /// Exception type thrown by AMReX on fatal errors.
    using RuntimeError = std::runtime_error;
}

#endif
