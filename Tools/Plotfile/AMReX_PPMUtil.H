#ifndef AMREX_PPM_UTIL_H_
#define AMREX_PPM_UTIL_H_

#include <string>
#include <AMReX_Array.H>

namespace amrex {

int loadPalette (const std::string& filename,
                 Array<unsigned char,256>& r, Array<unsigned char,256>& g,
                 Array<unsigned char,256>& b, Array<unsigned char,256>& a);

void storePPM (const std::string& filename,
               unsigned char const* data, int width, int height,
               Array<unsigned char,256> const& r,
               Array<unsigned char,256> const& g,
               Array<unsigned char,256> const& b);

}

#endif
