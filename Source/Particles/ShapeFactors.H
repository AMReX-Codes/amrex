#ifndef SHAPEFACTORS_H_
#define SHAPEFACTORS_H_

using namespace amrex;

// Compute shape factor and return index of leftmost cell where 
// particle writes.
// Specialized templates are defined below for orders 1, 2 and 3.
template <int depos_order>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor(Real* const sx, Real xint) {return 0;};

// Compute shape factor for order 0.
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor <0> (Real* const sx, Real xmid){
    int j = (int) (xmid+0.5);
    sx[0] = 1.0;
    return j;
}

// Compute shape factor for order 1.
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor <1> (Real* const sx, Real xmid){
    int j = (int) xmid;
    Real xint = xmid-j;
    sx[0] = 1.0 - xint;
    sx[1] = xint;
    return j;
}

// Compute shape factor for order 2.
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor <2> (Real* const sx, Real xmid){
    int j = (int) (xmid+0.5);
    Real xint = xmid-j;
    sx[0] = 0.5*(0.5-xint)*(0.5-xint);
    sx[1] = 0.75-xint*xint;
    sx[2] = 0.5*(0.5+xint)*(0.5+xint);
    // index of the leftmost cell where particle deposits
    return j-1;
}

// Compute shape factor for order 3.
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor <3> (Real* const sx, Real xmid){
    int j = (int) xmid;
    Real xint = xmid-j;
    sx[0] = 1.0/6.0*(1.0-xint)*(1.0-xint)*(1.0-xint);
    sx[1] = 2.0/3.0-xint*xint*(1-xint/2.0);
    sx[2] = 2.0/3.0-(1-xint)*(1-xint)*(1.0-0.5*(1-xint));
    sx[3] = 1.0/6.0*xint*xint*xint;
    // index of the leftmost cell where particle deposits
    return j-1;
}

#endif // SHAPEFACTORS_H_
