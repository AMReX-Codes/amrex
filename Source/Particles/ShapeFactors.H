/* Copyright 2019 Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef SHAPEFACTORS_H_
#define SHAPEFACTORS_H_

/**
 *  Compute shape factor and return index of leftmost cell where
 *  particle writes.
 *  Specialized templates are defined below for orders 0 to 3.
 */
template <int depos_order>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor(amrex::Real* const sx, amrex::Real xint)
{
    return 0;
};

/**
 *  Compute shape factor and return index of leftmost cell where
 *  particle writes.
 *  Specialization for order 0
 */
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor <0> (amrex::Real* const sx, amrex::Real xmid){
    const auto j = static_cast<int>(xmid+0.5);
    sx[0] = 1.0;
    return j;
}

/**
 *  Compute shape factor and return index of leftmost cell where
 *  particle writes.
 *  Specialization for order 1
 */
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor <1> (amrex::Real* const sx, amrex::Real xmid){
    const auto j = static_cast<int>(xmid);
    const amrex::Real xint = xmid-j;
    sx[0] = 1.0 - xint;
    sx[1] = xint;
    return j;
}

/**
 *  Compute shape factor and return index of leftmost cell where
 *  particle writes.
 *  Specialization for order 2
 */
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor <2> (amrex::Real* const sx, amrex::Real xmid){
    const auto j = static_cast<int>(xmid+0.5);
    const amrex::Real xint = xmid-j;
    sx[0] = 0.5*(0.5-xint)*(0.5-xint);
    sx[1] = 0.75-xint*xint;
    sx[2] = 0.5*(0.5+xint)*(0.5+xint);
    // index of the leftmost cell where particle deposits
    return j-1;
}

/**
 *  Compute shape factor and return index of leftmost cell where
 *  particle writes.
 *  Specialization for order 3
 */
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor <3> (amrex::Real* const sx, amrex::Real xmid){
    const auto j = static_cast<int>(xmid);
    const amrex::Real xint = xmid-j;
    sx[0] = 1.0/6.0*(1.0-xint)*(1.0-xint)*(1.0-xint);
    sx[1] = 2.0/3.0-xint*xint*(1-xint/2.0);
    sx[2] = 2.0/3.0-(1-xint)*(1-xint)*(1.0-0.5*(1-xint));
    sx[3] = 1.0/6.0*xint*xint*xint;
    // index of the leftmost cell where particle deposits
    return j-1;
}

/**
 *  Compute shifted shape factor and return index of leftmost cell where
 *  particle writes, for Esirkepov algorithm.
 *  Specialized templates are defined below for orders 1, 2 and 3.
 */
template <int depos_order>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shifted_shape_factor (amrex::Real* const sx,
                                  const amrex::Real x_old,
                                  const int i_new);

/**
 *  Compute shifted shape factor and return index of leftmost cell where
 *  particle writes, for Esirkepov algorithm.
 *  Specialization for order 1
 */
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shifted_shape_factor <1> (amrex::Real* const sx,
                                      const amrex::Real x_old,
                                      const int i_new){
    const auto i = static_cast<int>(x_old);
    const int i_shift = i - i_new;
    const amrex::Real xint = x_old - i;
    sx[1+i_shift] = 1.0 - xint;
    sx[2+i_shift] = xint;
    return i;
}

/**
 *  Compute shifted shape factor and return index of leftmost cell where
 *  particle writes, for Esirkepov algorithm.
 *  Specialization for order 2
 */
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shifted_shape_factor <2> (amrex::Real* const sx,
                                      const amrex::Real x_old,
                                      const int i_new){
    const auto i = static_cast<int>(x_old+0.5);
    const int i_shift = i - (i_new + 1);
    const amrex::Real xint = x_old - i;
    sx[1+i_shift] = 0.5*(0.5-xint)*(0.5-xint);
    sx[2+i_shift] = 0.75-xint*xint;
    sx[3+i_shift] = 0.5*(0.5+xint)*(0.5+xint);
    // index of the leftmost cell where particle deposits
    return i-1;
}

/**
 *  Compute shifted shape factor and return index of leftmost cell where
 *  particle writes, for Esirkepov algorithm.
 *  Specialization for order 3
 */
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shifted_shape_factor <3> (amrex::Real* const sx,
                                      const amrex::Real x_old,
                                      const int i_new){
    const auto i = static_cast<int>(x_old);
    const int i_shift = i - (i_new + 1);
    const amrex::Real xint = x_old - i;
    sx[1+i_shift] = 1.0/6.0*(1.0-xint)*(1.0-xint)*(1.0-xint);
    sx[2+i_shift] = 2.0/3.0-xint*xint*(1-xint/2.0);
    sx[3+i_shift] = 2.0/3.0-(1-xint)*(1-xint)*(1.0-0.5*(1-xint));
    sx[4+i_shift] = 1.0/6.0*xint*xint*xint;
    // index of the leftmost cell where particle deposits
    return i-1;
}

#endif // SHAPEFACTORS_H_
