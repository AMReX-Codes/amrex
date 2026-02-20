#ifndef AMREX_STRING_H_
#define AMREX_STRING_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_String.H
 * \brief String utility functions.
 */

#include <string>
#include <vector>

namespace amrex {

    /// \ingroup amrex_utilities
    //! Converts all characters of the string into lower case based on std::locale
    std::string toLower (std::string s);

    /// \ingroup amrex_utilities
    //! Converts all characters of the string into uppercase based on std::locale
    std::string toUpper (std::string s);

    /// \ingroup amrex_utilities
    //! Trim leading and trailing characters in the optional `space` argument.
    std::string trim (std::string s, std::string const& space = " \t");

    /// \ingroup amrex_utilities
    //! Returns rootNNNN where NNNN == num.
    std::string Concatenate (const std::string& root,
                             int                num,
                             int                mindigits = 5);

    /// \ingroup amrex_utilities
    //! Split a string using given tokens in `sep`.
    std::vector<std::string> split (std::string const& s,
                                    std::string const& sep = " \t");

    /// \ingroup amrex_utilities
    //! Join a vector of strings with given char `sep` as delimiter.
    std::string join (std::vector<std::string> const& sv, char sep);

    /// \ingroup amrex_utilities
    //! Join a vector of strings without delimiter.
    std::string join (std::vector<std::string> const& sv);
}

#endif
