#ifndef AMREX_MORTON_H_
#define AMREX_MORTON_H_

#include <AMReX_REAL.H>
#include <AMReX_Dim3.H>
#include <AMReX_Array.H>
#include <AMReX_GpuQualifiers.H>

#include <cstdint>

namespace amrex {
namespace Morton {

/**
 * \brief
 *  On input, the lowest 10 bits of x are filled, the next 22 are 0
 *  On output, x has been "stretched" to 30 bits, with each input bit
 *  followed by interleaving bits set to 0. In 2D, there is one
 *  interleaving bit, and in 3D, there are 2.
 *
 * \param x unsigned int holding the input to be split
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::uint32_t splitBits10 (std::uint32_t x) noexcept {
#if (AMREX_SPACEDIM == 3)
    // x            : 0000,0000,0000,0000,0000,00a9,8765,4321
    x = (x | (x << 16)) & 0x030000FF;
    // x << 16      : 0000,00a9,8765,4321,0000,0000,0000,0000
    // x | (x << 16): 0000,00a9,8765,4321,0000,00a9,8765,4321
    // 0x030000FF   : 0000,0011,0000,0000,0000,0000,1111,1111
    // x            : 0000,00a9,0000,0000,0000,0000,8765,4321
    x = (x | (x <<  8)) & 0x0300F00F;
    // x << 8       : 0000,0000,0000,0000,8765,4321,0000,0000
    // x | (x << 8) : 0000,00a9,0000,0000,8765,4321,8765,4321
    // 0x0300F00F   : 0000,0011,0000,0000,1111,0000,0000,1111
    // x            : 0000,00a9,0000,0000,8765,0000,0000,4321
    x = (x | (x <<  4)) & 0x030C30C3;
    // x << 4       : 00a9,0000,0000,8765,0000,0000,4321,0000
    // x | (x << 4) : 00a9,00a9,0000,8765,8765,0000,4321,4321
    // 0x030C30C3   : 0000,0011,0000,1100,0011,0000,1100,0011
    // x            : 0000,00a9,0000,8700,0065,0000,4300,0021
    x = (x | (x <<  2)) & 0x09249249;
    // x << 2       : 0000,a900,0087,0000,6500,0043,0000,2100
    // x | (x << 2) : 0000,a9a9,0087,8700,6565,0043,4300,2121
    // 0x09249249   : 0000,1001,0010,0100,1001,0010,0100,1001
    // x            : 0000,a009,0080,0700,6005,0040,0300,2001
    return x;
#elif (AMREX_SPACEDIM == 2)
    // x           : 0000,0000,0000,0000,gfed,cba9,8765,4321
    x = (x | (x << 8)) & 0x00FF00FF;
    // x << 8      : 0000,0000,gfed,cba9,8765,4321,0000,0000
    // x | (x << 8): 0000,0000,gfed,cba9,????,????,8765,4321
    // 0x00FF00FF  : 0000,0000,1111,1111,0000,0000,1111,1111
    // x           : 0000,0000,gfed,cba9,0000,0000,8765,4321
    x = (x | (x << 4)) & 0x0F0F0F0F;
    // x << 4      : 0000,gfed,cba9,0000,0000,8765,4321,0000
    // x | (x << 4): 0000,gfed,????,cba9,0000,8765,????,4321
    // 0x0F0F0F0F  : 0000,1111,0000,1111,0000,1111,0000,1111
    // x           : 0000,gfed,0000,cba9,0000,8765,0000,4321
    x = (x | (x << 2)) & 0x33333333;
    // x << 2      : 00gf,ed00,00cb,a900,0087,6500,0043,2100
    // x | (x << 2): 00gf,??ed,00cb,??a9,0087,??65,0043,??21
    // 0x33333333  : 0011,0011,0011,0011,0011,0011,0011,0011
    // x           : 00gf,00ed,00cb,00a9,0087,0065,0043,0021
    x = (x | (x << 1)) & 0x55555555;
    // x << 1      : 0gf0,0ed0,0cb0,0a90,0870,0650,0430,0210
    // x | (x << 1): 0g?f,0e?d,0c?b,0a?9,08?7,06?5,04?3,02?1
    // 0x55555555  : 0101,0101,0101,0101,0101,0101,0101,0101
    // x           : 0g0f,0e0d,0c0b,0a09,0807,0605,0403,0201
    return x;
#endif
}

/**
 * \brief
 *  Put the first 10 bits of x into retval with the 22 highest bits set to 0.
 *  It is assumed that xmin <= x < xmax.
 *
 * \param x the coordinate to convert.
 * \param xmin the low end of the domain.
 * \param xmax the high end of the domain.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::uint32_t packBits10 (amrex::Real x, amrex::Real xmin, amrex::Real xmax) noexcept {
    constexpr std::uint32_t code_offset = (1 << 10);
    return (std::uint32_t)(((x - xmin)/(xmax - xmin)) * code_offset);
}

/**
 * \brief
 * Given a 3D real-space coordinate, returns a 30 bit Morton code, stored in an
 * unsigned 32 bit int. This works by converting the coordinate to a 10-bit
 * nonnegative integer in each direction and then interleaving the bits.
 * It is assumed that plo <= p < phi has already been enforced.
 *
 * \param p the coordinate to convert.
 * \param plo the low end of the domain.
 * \param phi the high end of the domain.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::uint32_t get30BitCode (const XDim3& p,
                            const GpuArray<amrex::Real,AMREX_SPACEDIM>& plo,
                            const GpuArray<amrex::Real,AMREX_SPACEDIM>& phi) noexcept {
    std::uint32_t a = packBits10(p.x, plo[0], phi[0]);
    std::uint32_t b = packBits10(p.y, plo[1], phi[1]);
    std::uint32_t c = packBits10(p.z, plo[2], phi[2]);
    return splitBits10(a) | (splitBits10(b) << 1) | (splitBits10(c) << 2);
}

}
}
#endif
