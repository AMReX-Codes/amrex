
#ifndef BL_ARRAY_H
#define BL_ARRAY_H

#include <array>

namespace amrex {

    template <class T, std::size_t N>
    using Array = std::array<T,N>;

}

#ifndef FORTRAN_BOXLIB

#include <AMReX_SPACE.H>

namespace amrex
{
    template <class T>
    std::array<T*,AMREX_SPACEDIM> GetArrOfPtrs (std::array<T,AMREX_SPACEDIM>& a)
    {
        return {AMREX_D_DECL(&a[0], &a[1], &a[2])};
    }

    template <class T>
    std::array<T const*,AMREX_SPACEDIM> GetArrOfConstPtrs (const std::array<T,AMREX_SPACEDIM>& a)
    {
        return {AMREX_D_DECL(&a[0], &a[1], &a[2])};
    }
}

#endif

#endif /*BL_ARRAY_H*/
