#ifndef AMREX_DIM3_H_
#define AMREX_DIM3_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_Dim3.H
 * \brief Simple structs for 3-component integer and real-valued tuples.
 */

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

#include <iostream>

namespace amrex {

/**
 * \ingroup amrex_index_space
 * \brief A simple struct holding 3 int values for a 3D index.
 *
 * Used to represent a 3D integer index, independent of AMREX_SPACEDIM.
 */
struct Dim3 { int x; int y; int z; };

/**
 * \ingroup amrex_real_space
 * \brief A simple struct holding 3 Real values for a 3D coordinate.
 *
 * Used to represent a 3D real-valued coordinate, independent of AMREX_SPACEDIM.
 */
struct XDim3 { Real x; Real y; Real z; };

/**
 * \brief Return the dot product of two XDim3 vectors.
 * \param a First operand.
 * \param b Second operand.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real dot_product (XDim3 const& a, XDim3 const& b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

/**
 * \brief Return the cross product of two XDim3 vectors.
 * \param a First operand.
 * \param b Second operand.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
XDim3 cross_product (XDim3 const& a, XDim3 const& b)
{
    return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}

/**
 * \brief Component-wise addition of two XDim3 vectors.
 * \param a First operand.
 * \param b Second operand.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
XDim3 operator+ (XDim3 const& a, XDim3 const& b)
{
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

/**
 * \brief Component-wise subtraction of two XDim3 vectors.
 * \param a Minuend vector.
 * \param b Subtrahend vector.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
XDim3 operator- (XDim3 const& a, XDim3 const& b)
{
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

/**
 * \brief Write a Dim3 or XDim3 to an output stream in \c (x,y,z) format.
 * \param os Stream receiving formatted output.
 * \param d  Coordinate tuple to print.
 */
template <typename T,
          std::enable_if_t <std::is_same_v<T,  Dim3> ||
                                   std::is_same_v<T, XDim3>>* = nullptr>
std::ostream& operator<< (std::ostream& os, const T& d)
{
    os << '(' << d.x << ',' << d.y << ',' << d.z << ')';
    return os;
}

}

#endif
