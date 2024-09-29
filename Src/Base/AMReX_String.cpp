#include <AMReX_String.H>
#include <AMReX_BLassert.H>

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <sstream>

namespace amrex {

std::string toLower (std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return s;
}

std::string toUpper (std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    return s;
}

std::string trim(std::string s, std::string const& space)
{
    const auto sbegin = s.find_first_not_of(space);
    if (sbegin == std::string::npos) { return std::string{}; }
    const auto send = s.find_last_not_of(space);
    s = s.substr(sbegin, send-sbegin+1);
    return s;
}

std::string Concatenate (const std::string& root, int num, int mindigits)
{
    BL_ASSERT(mindigits >= 0);
    std::stringstream result;
    result << root << std::setfill('0') << std::setw(mindigits) << num;
    return result.str();
}

std::vector<std::string> split (std::string const& s, std::string const& sep)
{
    std::vector<std::string> result;
    std::size_t pos_begin, pos_end = 0;
    while ((pos_begin = s.find_first_not_of(sep,pos_end)) != std::string::npos) {
        pos_end = s.find_first_of(sep,pos_begin);
        result.push_back(s.substr(pos_begin,pos_end-pos_begin));
        if (pos_end == std::string::npos) { break; }
    }
    return result;
}

}
