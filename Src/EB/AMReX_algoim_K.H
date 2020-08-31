/*
 * This file contains a modified version of Algoim developed by
 * R. Saye, SIAM J. Sci. Comput., Vol. 37, No. 2, pp. A993-A1019,
 * http://dx.doi.org/10.1137/140966290, https://algoim.github.io/.
 *
 * Algoim Copyright (c) 2018, The Regents of the University of
 * California, through Lawrence Berkeley National Laboratory (subject
 * to receipt of any required approvals from the U.S. Dept. of
 * Energy). All rights reserved.
 */

#ifndef AMREX_ALGOIM_K_H_
#define AMREX_ALGOIM_K_H_

#include <AMReX_algoim.H>
#include <AMReX_Array.H>
#include <limits>
#include <cfloat>

namespace amrex { namespace algoim {

struct EBPlane
{
    GpuArray<Real,3> cent{};
    GpuArray<Real,3> norm{};

    EBPlane () = default;

    AMREX_GPU_HOST_DEVICE
    constexpr
    EBPlane (GpuArray<Real,3> const& c, GpuArray<Real,3> const n) noexcept
        : cent(c), norm(n)
        {}

    AMREX_GPU_HOST_DEVICE
    constexpr
    EBPlane (Real cx, Real cy, Real cz, Real nx, Real ny, Real nz) noexcept
        : cent{cx,cy,cz}, norm{nx,ny,nz}
        {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real operator() (Real x, Real y, Real z) const noexcept
    {
        return (x-cent[0])*norm[0] + (y-cent[1])*norm[1] + (z-cent[2])*norm[2];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real operator() (GpuArray<Real,3> const& p) const noexcept
    {
        return (p[0]-cent[0])*norm[0] + (p[1]-cent[1])*norm[1]
            +  (p[2]-cent[2])*norm[2];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real grad (int d) const noexcept
    {
        return norm[d];
    }
};

struct QuadratureRule
{
    struct Node
    {
        Real x, y, z, w;
    };
    int nnodes = 0;
    GpuArray<Node,256> nodes;

    template<typename F>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real
    operator() (const F& f) const noexcept
    {
        Real r = 0;
        for (int i = 0; i < nnodes; ++i) {
            r += f(nodes[i].x,nodes[i].y,nodes[i].z) * nodes[i].w;
        }
        return r;
    }

    template<typename F>
    AMREX_GPU_HOST AMREX_FORCE_INLINE
    Real
    eval (const F& f) const noexcept
    {
        Real r = 0;
        for (int i = 0; i < nnodes; ++i) {
            r += f(nodes[i].x,nodes[i].y,nodes[i].z) * nodes[i].w;
        }
        return r;
    }

    // evalIntegrand records quadrature nodes when it is invoked by
    // ImplicitIntegral
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void
    evalIntegrand (const GpuArray<Real,3>& x_, Real w_) noexcept
    {
        nodes[nnodes].x = x_[0];
        nodes[nnodes].y = x_[1];
        nodes[nnodes].z = x_[2];
        nodes[nnodes].w = w_;
        ++nnodes;
    }
};

namespace detail
{
    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void swap (T& a, T& b) noexcept {
        T c(a); a=b; b=c;
    }

    // Determines the sign conditions for restricting a (possibly
    // already restricted) level set function, i.e., sgn_L and sgn_U
    // in [R. Saye, High-Order Quadrature Methods for Implicitly
    // Defined Surfaces and Volumes in Hyperrectangles, SIAM
    // J. Sci. Comput., Vol. 37, No. 2, pp. A993-A1019,
    // http://dx.doi.org/10.1137/140966290].
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void
    determineSigns (bool positiveAbove, int sign, int& bottomSign, int& topSign) noexcept
    {
        // When evaluating a volume integral:
        //   If 0 integrating over the positive part:
        //      if positive above the height function: bottom = +/-, top = +
        //      if positive below the height function: bottom = +, top = +/-
        //   If integrating over the negative part:
        //      if positive above the height function: bottom = -, top = +/-
        //      if positive below the height function: bottom = +/-, top = -
        //   If integrating over both the positive and negative part (i.e. unrestricted in sign), keep it alive
        if (sign == 1)
        {
            bottomSign = positiveAbove? 0 : 1;
            topSign = positiveAbove? 1 : 0;
        }
        else if (sign == -1)
        {
            bottomSign = positiveAbove? -1 : 0;
            topSign = positiveAbove? 0 : -1;
        }
        else {
            bottomSign = topSign = 0;
        }
    }
} // namespace detail

// PsiCode encodes sign information of restricted level set functions
// on particular sides of a hyperrectangle in a packed array of
// bits. The first N bits encodes side information, the N+1st bit is
// true iff the sign == 0, while the N+2nd bit stores the sign if sign
// != 0.
template<int N>
struct PsiCode
{
    unsigned char bits = 0;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    PsiCode () noexcept {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    PsiCode (const GpuArray<int,N>& sides, int sign) noexcept
    {
        static_assert(N <= 3, "algoim::PsiCode: N must be <= 3");
        for (int dim = 0; dim < N; ++dim) {
            if (sides[dim] == 1) {
                bits |= (1 << dim);
            }
        }
        if (sign == 0) {
            bits |= (1 << N);
        } else {
            bits &= ~(1 << N);
            if (sign == 1) {
                bits |= (1 << (N+1));
            }
        }
    }

    // Modify an existing code by restriction in a particular dimension.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    PsiCode (const PsiCode& i, int dim, int side, int sign) noexcept
        : bits(i.bits)
    {
        if (side == 1) {
            bits |= (1 << dim);
        }
        if (sign == 0) {
            bits |= (1 << N);
        } else {
            bits &= ~(1 << N);
            if (sign == 1) {
                bits |= (1 << (N+1));
            }
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int side (int dim) const noexcept
    {
        return bits & (1 << dim);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int sign () const noexcept
    {
        return (bits & (1 << N)) ? 0 : ((bits & (1 << (N+1))) ? 1 : -1);
    }
};

template <int N>
struct BoundingBox
{
    static constexpr Real min (int /*dim*/) noexcept { return -0.5; }

    static constexpr Real max (int /*dim*/) noexcept { return  0.5; }

    static constexpr Real extent (int /*dim*/) noexcept { return 1.0; }

    static constexpr GpuArray<Real,N> midpoint () noexcept {
        return GpuArray<Real,N>{};
    }

    static constexpr Real midpoint (int /*dim*/) noexcept { return 0.0; }

    constexpr Real operator() (int side, int /*dim*/) const noexcept {
        return (side == 0) ? -0.5 : 0.5;
    }
};

struct Interval
{
    Real alpha;
};

// M-dimensional integral of an N-dimensional function restricted to
// given implicitly defined domains
template<int M, int N, typename Phi, typename F>
struct ImplicitIntegral
{
    const Phi& phi;
    F& f;
    const GpuArray<bool,N> free;
    GpuArray<PsiCode<N>,1 << (N - 1)> psi;
    int psiCount;
    const BoundingBox<N> xrange;
    int e0;
    GpuArray<Interval,N> xint;
    static constexpr int p = 4;

    // Prune the given set of functions by checking for the existence
    // of the interface. If a function is uniformly positive or
    // negative and is consistent with specified sign, it can be
    // removed. If a function is uniformly positive or negative but
    // inconsistent with specified sign, the domain of integration is
    // empty.
    AMREX_GPU_HOST_DEVICE
    bool prune () noexcept
    {
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
        static constexpr Real eps = std::numeric_limits<Real>::epsilon();
#else
        static constexpr Real eps = DBL_EPSILON;
#endif
        static constexpr Real almostone = 1.0-10.*eps;
        for (int i = 0; i < psiCount; )
        {
            GpuArray<Real,N> mid = xrange.midpoint();
            Real dphi_max = 0.0;
            for (int dim = 0; dim < N; ++dim) {
                if (free[dim]) {
                    dphi_max += amrex::Math::abs(phi.grad(dim));
                } else {
                    mid[dim] = xrange(psi[i].side(dim),dim);
                }
            }
            dphi_max *= 0.5*almostone;
            const Real phi_0 = phi(mid);
            bool uniform_sign = (phi_0 > dphi_max) or (phi_0 < -dphi_max);
            if (uniform_sign)
            {
                if ((phi_0 >= 0.0 and psi[i].sign() >= 0) or
                    (phi_0 <= 0.0 and psi[i].sign() <= 0) )
                {
                    --psiCount;
                    detail::swap(psi[i], psi[psiCount]);
                }
                else {
                    return false;
                }
            }
            else {
                ++i;
            }
        }
        return true;
    }

    // Gaussian quadrature for when the domain of integration is
    // determined to be the entire M-dimensional cube.
    AMREX_GPU_HOST_DEVICE
    void tensorProductIntegral () noexcept
    {
        constexpr Real gauss_x[4]={0.069431844202973712388026755553595247452,
                                   0.33000947820757186759866712044837765640,
                                   0.66999052179242813240133287955162234360,
                                   0.93056815579702628761197324444640475255};
        constexpr Real gauss_w[4]={0.173927422568726928686531974610999703618,
                                   0.326072577431273071313468025389000296382,
                                   0.326072577431273071313468025389000296382,
                                   0.173927422568726928686531974610999703618};

        int nloops = 1;
        for (int i = 0; i < M; ++i) {
            nloops *= p;
        }
        GpuArray<int,M> i{}; // i is initialized to zero.
        for (int iloop = 0; iloop < nloops; ++iloop)
        {
            GpuArray<Real,N> x{};
            Real w = 1.0;
            for (int dim = 0, k = 0; dim < N; ++dim) {
                if (free[dim]) {
                    x[dim] = xrange.min(dim) + xrange.extent(dim) * gauss_x[i[k]];
                    w *= xrange.extent(dim) * gauss_w[i[k]];
                    ++k;
                }
            }
            f.evalIntegrand(x, w);

            for (int m = M-1; m >= 0; --m) {
                ++(i[m]);
                if (i[m] < p) {
                    break;
                } else {
                    i[m] = 0;
                }
            }
        }
    }

    // Given x, valid in all free variables barring e0, root find in
    // the e0 direction on each of the implicit functions, and apply
    // Gaussian quadrature to each segment. Weights are multiplied
    // upon going back up the tree of recursive calls.
    AMREX_GPU_HOST_DEVICE
    void evalIntegrand (GpuArray<Real,N> x, Real w) const noexcept
    {
        constexpr Real gauss_x[4]={0.069431844202973712388026755553595247452,
                                   0.33000947820757186759866712044837765640,
                                   0.66999052179242813240133287955162234360,
                                   0.93056815579702628761197324444640475255};
        constexpr Real gauss_w[4]={0.173927422568726928686531974610999703618,
                                   0.326072577431273071313468025389000296382,
                                   0.326072577431273071313468025389000296382,
                                   0.173927422568726928686531974610999703618};

        // Each thread has its own storage for computing roots. These
        // are not corrupted by the recursive chain of evalIntegrand()
        // calls since each call is in a different templated
        // namespace. Moreover, if using OpenMP tasking, each task is
        // assumed tied and so one thread can only execute
        // "evalIntegrand" from start to end uninterrupted.
        GpuArray<Real,4> roots;

        // Partition [xmin(e0), xmax(e0)] by roots found
        Real x_min = xrange.min(e0);
        Real x_max = xrange.max(e0);
        roots[0] = x_min;
        int nroots = 1;
        if (phi.grad(e0) != Real(0.0)) {
            for (int i = 0; i < psiCount; ++i)
            {
                for (int dim = 0; dim < N; ++dim) {
                    if (!free[dim]) {
                        x[dim] = xrange(psi[i].side(dim),dim);
                    }
                }
                // x is now valid in all variables except e0
                x[e0] = x_min;
                Real phi_lo = phi(x);
                Real xroot = x[e0] - phi_lo / phi.grad(e0);
                if (xroot > x_min and xroot < x_max) {
                    roots[nroots++] = xroot;
                }
            }
            AMREX_ASSERT(nroots <= 3);
        }
        if (nroots == 3 and roots[1] > roots[2]) {
            detail::swap(roots[1],roots[2]);
        }
        roots[nroots++] = x_max;

        // In rare cases, degenerate segments can be found, filter out with a tolerance
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
        static constexpr Real eps = std::numeric_limits<Real>::epsilon();
#else
        static constexpr Real eps = DBL_EPSILON;
#endif
        constexpr Real tol = 50.0*eps;

        // Loop over segments of divided interval
        for (int i = 0; i < nroots - 1; ++i)
        {
            if (roots[i+1] - roots[i] < tol) continue;

            // Evaluate sign of phi within segment and check for consistency with psi
            bool okay = true;
            x[e0] = (roots[i] + roots[i+1])*0.5;
            for (int j = 0; j < psiCount && okay; ++j)
            {
                for (int dim = 0; dim < N; ++dim) {
                    if (!free[dim]) {
                        x[dim] = xrange(psi[j].side(dim), dim);
                    }
                }
                bool new_ok = (phi(x) > 0.0) ? (psi[j].sign() >= 0) : (psi[j].sign() <= 0);
                okay = okay and new_ok;
            }
            if (!okay) continue;

            for (int j = 0; j < p; ++j)
            {
                x[e0] = roots[i] + (roots[i+1] - roots[i]) * gauss_x[j];
                Real gw = (roots[i+1] - roots[i]) * gauss_w[j];
                f.evalIntegrand(x, w * gw);
            }
        }
    }

    // Main calling engine; parameters with underscores are copied
    // upon entry but modified internally in the ctor
    AMREX_GPU_HOST_DEVICE
    ImplicitIntegral (const Phi& phi_, F& f_, const GpuArray<bool,N>& free_,
                      const GpuArray<PsiCode<N>,1 << (N-1)>& psi_,
                      int psiCount_) noexcept
        : phi(phi_), f(f_), free(free_), psi(psi_), psiCount(psiCount_),
          xrange()
    {
        // For the one-dimensional base case, evaluate the
        // bottom-level integral.
        if (M == 1)
        {
            for (int dim = 0; dim < N; ++dim) {
                if (free[dim]) {
                    e0 = dim;
                }
            }
            evalIntegrand(GpuArray<Real,N>{}, 1.0);
            return;
        }

        // Establish interval bounds for prune() and remaining part of ctor.
        for (int dim = 0; dim < N; ++dim)
        {
            if (free[dim])
            {
                xint[dim].alpha = xrange.midpoint(dim);
            }
            else
            {
                xint[dim].alpha = 0.0;
            }
        }

        // Prune list of psi functions: if prune procedure returns
        // false, then the domain of integration is empty.
        if (!prune()) {
            return;
        }

        // If all psi functions were pruned, then the volumetric
        // integral domain is the entire hyperrectangle.
        if (psiCount == 0)
        {
            tensorProductIntegral();
            return;
        }

        // Among all monotone height function directions, choose the
        // one that makes the associated height function look as flat
        // as possible.  This is a modification to the criterion
        // presented in [R. Saye, High-Order Quadrature Methods for
        // Implicitly Defined Surfaces and Volumes in Hyperrectangles,
        // SIAM J. Sci. Comput., Vol. 37, No. 2, pp. A993-A1019,
        // http://dx.doi.org/10.1137/140966290].
        e0 = -1;
        {
            Real gmax = -1.;
            for (int dim = 0; dim < N; ++dim) {
                if (free[dim] && amrex::Math::abs(phi.grad(dim)) > gmax) {
                    gmax = amrex::Math::abs(phi.grad(dim));
                    e0 = dim;
                }
            }
        }

        // Check compatibility with all implicit functions whilst
        // simultaneously constructing new implicit functions.
        GpuArray<PsiCode<N>,1 << (N-1)> newPsi;
        int newPsiCount = 0;
        for (int i = 0; i < psiCount; ++i)
        {
            // Evaluate gradient in an interval
            for (int dim = 0; dim < N; ++dim) {
                if (!free[dim]) {
                    xint[dim].alpha = xrange(psi[i].side(dim), dim);
                }
            }

            // w.z. We assume grad does not change.
            int bottomSign, topSign;
            detail::determineSigns(phi.grad(e0) > 0.0, psi[i].sign(),
                                   bottomSign, topSign);
            // w.z. There are sides 0 and 1.
            newPsi[newPsiCount++] = PsiCode<N>(psi[i], e0, 0, bottomSign);
            newPsi[newPsiCount++] = PsiCode<N>(psi[i], e0, 1, topSign);
        }

        // Dimension reduction call
        GpuArray<bool,N> new_free = free;
        new_free[e0] = false;
        ImplicitIntegral<M-1,N,Phi,ImplicitIntegral<M,N,Phi,F> >
            (phi, *this, new_free, newPsi, newPsiCount);
    }
};

// Partial specialisation on M=0 as a dummy base case for the compiler
template<int N, typename Phi, typename F>
struct ImplicitIntegral<0,N,Phi,F>
{
    AMREX_GPU_HOST_DEVICE
    ImplicitIntegral (const Phi&, F&, const GpuArray<bool,N>&,
                      const GpuArray<PsiCode<N>,1 << (N-1)>&,
                      int) noexcept {}
};

AMREX_GPU_HOST_DEVICE inline
QuadratureRule
quadGen (EBPlane const& phi) noexcept
{
    QuadratureRule q;
    GpuArray<bool,3> free{true,true,true};
    GpuArray<PsiCode<3>,4> psi;
    psi[0] = PsiCode<3>({0,0,0}, -1);
    ImplicitIntegral<3,3,EBPlane,QuadratureRule>(phi, q, free, psi, 1);
    return q;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void set_regular (int i, int j, int k, Array4<Real> const& intg) noexcept
{
    constexpr Real twelfth = 1./12.;
    constexpr Real   offth = 1./144.;
    intg(i,j,k,i_S_x    ) = 0.0;
    intg(i,j,k,i_S_y    ) = 0.0;
    intg(i,j,k,i_S_z    ) = 0.0;
    intg(i,j,k,i_S_x2   ) = twelfth;
    intg(i,j,k,i_S_y2   ) = twelfth;
    intg(i,j,k,i_S_z2   ) = twelfth;
    intg(i,j,k,i_S_x_y  ) = 0.0;
    intg(i,j,k,i_S_x_z  ) = 0.0;
    intg(i,j,k,i_S_y_z  ) = 0.0;
    intg(i,j,k,i_S_x2_y ) = 0.0;
    intg(i,j,k,i_S_x2_z ) = 0.0;
    intg(i,j,k,i_S_x_y2 ) = 0.0;
    intg(i,j,k,i_S_y2_z ) = 0.0;
    intg(i,j,k,i_S_x_z2 ) = 0.0;
    intg(i,j,k,i_S_y_z2 ) = 0.0;
    intg(i,j,k,i_S_x2_y2) = offth;
    intg(i,j,k,i_S_x2_z2) = offth;
    intg(i,j,k,i_S_y2_z2) = offth;
    intg(i,j,k,i_S_xyz  ) = 0.0;
}

}}

#endif
