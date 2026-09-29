#ifndef AMREX_EXCEPTION_H_
#define AMREX_EXCEPTION_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_Exception.H
 * \brief AMReX exception types.
 */

#include <new>
#include <stdexcept>
#include <string>

namespace amrex
{
    /// \ingroup amrex_utilities
    /// Exception type used when AMReX is configured to throw on fatal errors.
    using RuntimeError = std::runtime_error;

    /// \ingroup amrex_utilities
    /// Thrown when an Arena cannot allocate the requested memory.
    /// `what()` describes the failed request and the current memory usage.
    class OutOfMemoryError
        : public std::bad_alloc
    {
    public:
        explicit OutOfMemoryError (std::string const& msg) : m_msg(msg) {}
        [[nodiscard]] const char* what () const noexcept override { return m_msg.what(); }
    private:
        std::runtime_error m_msg; // noexcept copy
    };
}

#endif
