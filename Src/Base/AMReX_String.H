#ifndef AMREX_STRING_H_
#define AMREX_STRING_H_
#include <AMReX_Config.H>

#include <string>
#include <vector>

namespace amrex {

    //! Converts all characters of the string into lower case based on std::locale
    std::string toLower (std::string s);

    //! Converts all characters of the string into uppercase based on std::locale
    std::string toUpper (std::string s);

    //! Trim leading and trailing characters in the optional `space`
    //! argument.
    std::string trim (std::string s, std::string const& space = " \t");

    //! Returns rootNNNN where NNNN == num.
    std::string Concatenate (const std::string& root,
                             int                num,
                             int                mindigits = 5);

    //! Split a string using given tokens in `sep`.
    std::vector<std::string> split (std::string const& s,
                                    std::string const& sep = " \t");
}

#endif
