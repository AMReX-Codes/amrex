#ifndef AMREX_EB2_GEOMETRYSHOP_H_
#define AMREX_EB2_GEOMETRYSHOP_H_
#include <AMReX_Config.H>

#include <AMReX_EB2_IF_Base.H>
#include <AMReX_EB2_Graph.H>
#include <AMReX_Geometry.H>
#include <AMReX_BaseFab.H>
#include <AMReX_Print.H>
#include <AMReX_Array.H>
#include <memory>
#include <type_traits>
#include <cmath>
#include <limits>
#include <utility>

namespace amrex::EB2 {

template <class F>
requires (IsGPUable<F>::value)
AMREX_GPU_HOST_DEVICE
Real
IF_f (F const& f, GpuArray<Real,AMREX_SPACEDIM> const& p) noexcept
{
    return f(AMREX_D_DECL(p[0],p[1],p[2]));
}

template <class F>
requires (!IsGPUable<F>::value)
AMREX_GPU_HOST_DEVICE
Real
IF_f (F const& f, GpuArray<Real,AMREX_SPACEDIM> const& p)
{
    AMREX_IF_ON_DEVICE((
        amrex::ignore_unused(f,p);
        amrex::Error("EB2::GeometryShop: how did this happen?");
        return 0.0;
    ))
    AMREX_IF_ON_HOST((return f({AMREX_D_DECL(p[0],p[1],p[2])});))
}

template <class F>
AMREX_GPU_HOST_DEVICE
Real
BrentRootFinder (GpuArray<Real,AMREX_SPACEDIM> const& lo,
                 GpuArray<Real,AMREX_SPACEDIM> const& hi,
                 int rangedir, F const& f)
{
    // The tolerance must not exceed a fraction of the edge length.
    const Real h = std::abs(hi[rangedir] - lo[rangedir]);
#ifdef AMREX_USE_FLOAT
    const Real tol = amrex::min(1.e-4_rt, 1.e-2_rt*h);
    const Real EPS = 1.0e-6_rt;
#else
    const Real tol = amrex::min(1.e-12, 1.e-10*h);
    const Real EPS = 3.0e-15;
#endif
    const int MAXITER = 100;

    Real p, q, r, s;

    GpuArray<Real,AMREX_SPACEDIM> aPt = lo;
    GpuArray<Real,AMREX_SPACEDIM> bPt = hi;

    Real fa = IF_f(f, aPt);
    Real fb = IF_f(f, bPt);
    Real c = bPt[rangedir];
    Real fc = fb;

    if (fb*fa > 0.0_rt) {
//        amrex::AllPrint() << "fa " << fa << " fb " << fb << "\n";
        amrex::Error("BrentRootFinder. Root must be bracketed, but instead the supplied end points have the same sign.");
        return 0.0_rt;
    } else if (fa == 0.0_rt) {
        return aPt[rangedir];
    } else if (fb == 0.0_rt) {
        return bPt[rangedir];
    }

    Real d = 0.0_rt, e = 0.0_rt;
    int i;
    for (i = 0; i < MAXITER; ++i)
    {
        if (fb*fc > 0)
        {
            //  Rename a, b, c and adjust bounding interval d
            c = aPt[rangedir];
            fc  = fa;
            d = bPt[rangedir] - aPt[rangedir];
            e = d;
        }

        if (std::abs(fc) < std::abs(fb))
        {
            aPt[rangedir] = bPt[rangedir];
            bPt[rangedir] = c;
            c = aPt[rangedir];
            fa  = fb;
            fb  = fc;
            fc  = fa;
        }

        //  Convergence check
        Real tol1  = 2.0_rt * EPS * std::abs(bPt[rangedir]) + 0.5_rt * tol;
        Real xm    = 0.5_rt * (c - bPt[rangedir]);

        if (std::abs(xm) <= tol1 || fb == 0.0_rt)
        {
            break;
        }

        if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb))
        {
            //  Attempt inverse quadratic interpolation
            s = fb / fa;
            if (aPt[rangedir] == c)
            {
                p = 2.0_rt * xm * s;
                q = 1.0_rt - s;
            }
            else
            {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0_rt * xm * q * (q-r) - (bPt[rangedir]-aPt[rangedir]) * (r-1.0_rt));
                q = (q-1.0_rt) * (r-1.0_rt) * (s-1.0_rt);
            }

            //  Check whether in bounds
            if (p > 0) { q = -q; }

            p = std::abs(p);

            if (2.0_rt * p < amrex::min(3.0_rt*xm*q-std::abs(tol1*q), 1.0_rt*std::abs(e*q)))
            {
                //  Accept interpolation
                e = d;
                d = p / q;
            }
            else
            {
                //  Interpolation failed, use bisection
                d = xm;
                e = d;
            }
        }
        else
        {
            //  Bounds decreasing too slowly, use bisection
            d = xm;
            e = d;
        }

        //  Move last best guess to a
        aPt[rangedir] = bPt[rangedir];
        fa  = fb;

        //  Evaluate new trial root
        if (std::abs(d) > tol1)
        {
            bPt[rangedir] = bPt[rangedir] + d;
        }
        else
        {
            if (xm < 0) { bPt[rangedir] = bPt[rangedir] - tol1; }
            else        { bPt[rangedir] = bPt[rangedir] + tol1; }
        }

        fb = IF_f(f, bPt);
    }

    if (i >= MAXITER)
    {
        amrex::Error("BrentRootFinder: exceeding maximum iterations.");
    }

    return bPt[rangedir];
}

// Inverse squared distance to a neighboring node, by the bits of its
// nonzero index offsets, for IF_roundoff.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuArray<Real,8> IF_invdist2 (GpuArray<Real,AMREX_SPACEDIM> const& dx)
{
    Real const d2x = dx[0]*dx[0];
    Real const d2y = dx[1]*dx[1];
    Real const d2z = (AMREX_SPACEDIM == 3) ? dx[AMREX_SPACEDIM-1]*dx[AMREX_SPACEDIM-1] : 0.0_rt;
    GpuArray<Real,8> invd2;
    invd2[0] = 0.0_rt;
    for (int m = 1; m < 8; ++m) {
        Real const den = Real(m&1)*d2x + Real((m>>1)&1)*d2y + Real((m>>2)&1)*d2z;
        invd2[m] = (den > 0.0_rt) ? Real(1.0)/den : 0.0_rt; // z entries unused in 2D
    }
    return invd2;
}

// Size of the coordinate operands anywhere in the domain, for IF_roundoff.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real IF_size (GpuArray<Real,AMREX_SPACEDIM> const& problo,
              GpuArray<Real,AMREX_SPACEDIM> const& probhi)
{
    Real size = 0.0_rt;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        size += std::abs(problo[d]) + (probhi[d]-problo[d]);
    }
    return size;
}

// Upper bound of IF_roundoff over a fab whose values lie in [vmin,vmax]
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real IF_tolmax (Real vmin, Real vmax, GpuArray<Real,8> const& invd2, Real size)
{
    Real const dv = vmax - vmin;
    Real invd2max = invd2[1];
    for (int m = 2; m < 8; ++m) { invd2max = amrex::max(invd2max, invd2[m]); }
    return Real(16.0) * std::numeric_limits<Real>::epsilon() * std::sqrt(dv*dv*invd2max) * size;
}

// Tolerance below which f at node iv is roundoff: 16 epsilon times the
// largest slope of f to a neighboring node, diagonals included so that a
// kink of the surface at the node does not hide the slope, times size (see
// IF_size).  size is uniform over the domain so that a surface at the same
// distance from a whole row of nodes is classified the same along the row.
// Neighbor values are taken from fab where it has them.  Nodes outside
// [blo,bhi] take the neighbor values of the clamped node.
template <class F>
AMREX_GPU_HOST_DEVICE
Real
IF_roundoff (F const& f, IntVect const& iv, Real fv, Array4<Real const> const& fab,
             GpuArray<Real,AMREX_SPACEDIM> const& problo,
             GpuArray<Real,AMREX_SPACEDIM> const& dx, GpuArray<Real,8> const& invd2,
             Real size, IntVect const& blo, IntVect const& bhi)
{
    IntVect ivc = iv;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { ivc[d] = amrex::Clamp(iv[d], blo[d], bhi[d]); }
    constexpr int kmax = (AMREX_SPACEDIM == 3) ? 1 : 0;
    int const ic = ivc[0];
    int const jc = ivc[1];
    int const kc = (AMREX_SPACEDIM == 3) ? ivc[AMREX_SPACEDIM-1] : 0;
    Real slope2 = 0.0_rt;
    IntVect const ivlo = ivc - 1;
    IntVect const ivhi = ivc + 1;
    if (ivlo.allGE(blo) && ivhi.allLE(bhi) && fab.contains(ivlo) && fab.contains(ivhi)) {
        // Same result as the loop below, grouped by invd2.
        Real const* q = fab.ptr(ic,jc,kc);
        Long const js = fab.template get_stride<1>();
        auto sq = [=] (Real a) { return (a-fv)*(a-fv); };
        slope2 = amrex::max(sq(q[-1]), sq(q[1])) * invd2[1];
        slope2 = amrex::max(slope2, amrex::max(sq(q[-js]), sq(q[js])) * invd2[2]);
        slope2 = amrex::max(slope2, amrex::max(sq(q[-1-js]), sq(q[1-js]),
                                               sq(q[-1+js]), sq(q[1+js])) * invd2[3]);
#if (AMREX_SPACEDIM == 3)
        Long const ks = fab.template get_stride<2>();
        slope2 = amrex::max(slope2, amrex::max(sq(q[-ks]), sq(q[ks])) * invd2[4]);
        slope2 = amrex::max(slope2, amrex::max(sq(q[-1-ks]), sq(q[1-ks]),
                                               sq(q[-1+ks]), sq(q[1+ks])) * invd2[5]);
        slope2 = amrex::max(slope2, amrex::max(sq(q[-js-ks]), sq(q[js-ks]),
                                               sq(q[-js+ks]), sq(q[js+ks])) * invd2[6]);
        slope2 = amrex::max(slope2, amrex::max(sq(q[-1-js-ks]), sq(q[1-js-ks]),
                                               sq(q[-1+js-ks]), sq(q[1+js-ks]),
                                               sq(q[-1-js+ks]), sq(q[1-js+ks]),
                                               sq(q[-1+js+ks]), sq(q[1+js+ks])) * invd2[7]);
#endif
    } else {
        for (int kk = -kmax; kk <= kmax; ++kk) {
        for (int jj = -1; jj <= 1; ++jj) {
        for (int ii = -1; ii <= 1; ++ii) {
            int const m = int(ii != 0) | (int(jj != 0) << 1) | (int(kk != 0) << 2);
            if (m != 0) {
                IntVect ivn = ivc + IntVect(AMREX_D_DECL(ii,jj,kk));
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    ivn[d] = amrex::Clamp(ivn[d], blo[d], bhi[d]);
                }
                Real fn;
                if (fab.contains(ivn)) {
                    fn = fab(ivn);
                } else {
                    GpuArray<Real,AMREX_SPACEDIM> p;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        p[d] = node_coord(problo[d], ivn[d], dx[d]);
                    }
                    fn = IF_f(f, p);
                }
                Real const df = fn - fv;
                slope2 = amrex::max(slope2, df*df*invd2[m]);
            }
        }}}
    }
    return Real(16.0) * std::numeric_limits<Real>::epsilon() * std::sqrt(slope2) * size;
}

// Intercept on the edge from node ivz, which is on the surface (zero level
// set), to the fluid node one cell away in direction rangedir.  This is
// where f leaves the surface, which is not ivz if the edge runs along it.
template <class F>
AMREX_GPU_HOST_DEVICE
Real
ZeroEndIntercept (F const& f, IntVect const& ivz, int rangedir, int side,
                  Array4<Real const> const& lst,
                  GpuArray<Real,AMREX_SPACEDIM> const& problo,
                  GpuArray<Real,AMREX_SPACEDIM> const& dx, Real size,
                  IntVect const& blo, IntVect const& bhi)
{
    GpuArray<Real,AMREX_SPACEDIM> p;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        p[d] = node_coord(problo[d], amrex::Clamp(ivz[d],blo[d],bhi[d]), dx[d]);
    }
    Real const tol = IF_roundoff(f, ivz, IF_f(f, p), lst, problo, dx, IF_invdist2(dx), size,
                                 blo, bhi);
    IntVect ivf = ivz;
    ivf[rangedir] += side;
    Real a = p[rangedir];
    Real b = node_coord(problo[rangedir], amrex::Clamp(ivf[rangedir],blo[rangedir],bhi[rangedir]),
                        dx[rangedir]);
    Real const xtol = Real(4.0) * std::numeric_limits<Real>::epsilon()
        * amrex::max(std::abs(a), std::abs(b), std::abs(problo[rangedir]));
    Real const x0 = a;
    p[rangedir] = b;
    Real const fb = IF_f(f, p);
    Real const h = std::abs(b - a);
    // Within roundoff of the node, the surface just passes through it.
    // Roundoff here is where f drops below -tol at the slope along the edge.
    Real const node_tol = (fb != 0.0_rt) ? Real(2.0)*tol*h/std::abs(fb) : 0.0_rt;
    // The usual case: f is already below -tol at node_tol, where the
    // bisection would end.
    if (node_tol >= h) { return x0; }
    p[rangedir] = x0 + Real(side)*node_tol;
    if (IF_f(f, p) < -tol) { return x0; }
    for (int iter = 0; iter < 200 && std::abs(b-a) > xtol; ++iter) {
        p[rangedir] = Real(0.5)*(a+b);
        if (IF_f(f, p) >= -tol) {
            a = p[rangedir];
        } else {
            b = p[rangedir];
        }
    }
    return (std::abs(a - x0) <= node_tol) ? x0 : a;
}

/**
 * \brief Samples an implicit function on boxes and fills EB data.
 *
 * \tparam F Functor providing the implicit function value.
 * \tparam R Optional resource handle kept alive as long as the shop exists.
 */
template <class F, class R = int>
class GeometryShop
{
public:

    static constexpr int in_fluid = -1;
    static constexpr int on_boundary = 0;
    static constexpr int in_body = 1;
    //
    static constexpr int allregular = -1;
    static constexpr int mixedcells = 0;
    static constexpr int allcovered = 1;

    using FunctionType = F;

    /**
     * \brief Construct a GeometryShop that owns an implicit function \p f.
     *
     * \param f Implicit function whose sign distinguishes inside/outside.
     */
    explicit GeometryShop (F f)
        : m_f(std::move(f)), m_resource()
        {}

    /**
     * \brief Construct a GeometryShop with an implicit function and resource handle.
     *
     * \param f Implicit function being sampled.
     * \param r Auxiliary resource (e.g., Parser data) kept alive while \p f is in use.
     */
    GeometryShop (F f, R r)
        : m_f(std::move(f)), m_resource(std::move(r))
        {}

    //! Immutable access to the underlying implicit function.
    [[nodiscard]] F const& GetImpFunc () const& { return m_f; }
    F&& GetImpFunc () && { return std::move(m_f); }

    // Nodes are classified as in fillFab, with the same neighbors and
    // clamping for the roundoff tolerance, so that a box is regular or
    // covered only if fillFab would make no cut cells in it.  Both signs
    // make the box mixed.  If all nodes are fluid, only those within
    // IF_tolmax can be body, which needs the roundoff stencil.
    [[nodiscard]] int getBoxType_Cpu (const Box& bx, Geometry const& geom,
                                      Box const& bounding_box) const noexcept
    {
        const auto& problo = geom.ProbLoArray();
        const auto& dx = geom.CellSizeArray();
        const Box gbx = amrex::grow(bx, 1);
        const auto& blo = amrex::lbound(bounding_box);
        const auto& bhi = amrex::ubound(bounding_box);
        BaseFab<Real> ffab(gbx, 1, The_Cpu_Arena());
        const auto& fv = ffab.array();
        const auto& fc = ffab.const_array();
        const auto& lo = amrex::lbound(bx);
        const auto& hi = amrex::ubound(bx);
        int npos = 0, nneg = 0;
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            Real const v = m_f(RealArray{AMREX_D_DECL(node_coord(problo[0], amrex::Clamp(i,blo.x,bhi.x), dx[0]),
                                                      node_coord(problo[1], amrex::Clamp(j,blo.y,bhi.y), dx[1]),
                                                      node_coord(problo[2], amrex::Clamp(k,blo.z,bhi.z), dx[2]))});
            fv(i,j,k) = v;
            if (v >= 0.0_rt) { ++npos; } else { ++nneg; }
            if (npos > 0 && nneg > 0) { return mixedcells; }
        }}}
        if (nneg == 0) { return allcovered; }

        Real vmin = std::numeric_limits<Real>::max();
        Real vmax = std::numeric_limits<Real>::lowest();
        amrex::LoopOnCpu(gbx, [&] (int i, int j, int k) noexcept
        {
            amrex::ignore_unused(j,k);
            if (!bx.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                fv(i,j,k) = m_f(RealArray{AMREX_D_DECL(node_coord(problo[0], amrex::Clamp(i,blo.x,bhi.x), dx[0]),
                                                       node_coord(problo[1], amrex::Clamp(j,blo.y,bhi.y), dx[1]),
                                                       node_coord(problo[2], amrex::Clamp(k,blo.z,bhi.z), dx[2]))});
            }
            vmin = amrex::min(vmin, fc(i,j,k));
            vmax = amrex::max(vmax, fc(i,j,k));
        });
        IntVect const glo = bounding_box.smallEnd();
        IntVect const ghi = bounding_box.bigEnd();
        auto const invd2 = IF_invdist2(dx);
        auto const size = IF_size(problo, geom.ProbHiArray());
        Real const tolmax = IF_tolmax(vmin, vmax, invd2, size);
        if (vmax < -tolmax) { return allregular; }
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            IntVect const iv(AMREX_D_DECL(i,j,k));
            Real const v = fc(i,j,k);
            if (std::abs(v) <= tolmax &&
                std::abs(v) <= IF_roundoff(m_f, iv, v, fc, problo, dx, invd2, size, glo, ghi)) {
                return mixedcells;
            }
        }}}
        return allregular;
    }

    /**
     * \brief Classify \p bx as regular, covered, or mixed (GPU-capable case).
     *        Nodes on the surface within roundoff count as body, as in fillFab().
     *
     * \param bx   Box to inspect.
     * \param geom Geometry providing physical coordinates.
     * \param run_on Whether to run on CPU or GPU.
     * \param bounding_box Bounding box limiting the sampling region.
     */
    template <class U=F>
    requires (IsGPUable<U>::value)
    [[nodiscard]] int getBoxType (const Box& bx, const Geometry& geom, RunOn run_on,
                                  Box const& bounding_box) const noexcept
    {
        if (run_on == RunOn::Gpu && Gpu::inLaunchRegion())
        {
            const auto& problo = geom.ProbLoArray();
            const auto& dx = geom.CellSizeArray();
            auto const& f = m_f;
            // Node classification as in getBoxType_Cpu.
            const Box gbx = amrex::grow(bx, 1);
            const auto& blo = amrex::lbound(bounding_box);
            const auto& bhi = amrex::ubound(bounding_box);
            BaseFab<Real> ffab(gbx, 1, The_Async_Arena());
            auto const& fv = ffab.array();
            auto const& fc = ffab.const_array();
            ReduceOps<ReduceOpSum,ReduceOpSum,ReduceOpMin,ReduceOpMax> reduce_op;
            ReduceData<int,int,Real,Real> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
            reduce_op.eval(gbx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                amrex::ignore_unused(j,k);
                Real const v = f(AMREX_D_DECL(node_coord(problo[0], amrex::Clamp(i,blo.x,bhi.x), dx[0]),
                                              node_coord(problo[1], amrex::Clamp(j,blo.y,bhi.y), dx[1]),
                                              node_coord(problo[2], amrex::Clamp(k,blo.z,bhi.z), dx[2])));
                fv(i,j,k) = v;
                bool const inbx = bx.contains(IntVect(AMREX_D_DECL(i,j,k)));
                return {int(inbx && v >= 0.0_rt), int(inbx && v < 0.0_rt), v, v};
            });
            ReduceTuple hv = reduce_data.value(reduce_op);
            int const npos = amrex::get<0>(hv);
            int const nneg = amrex::get<1>(hv);
            if (npos > 0 && nneg > 0) { return mixedcells; }
            if (nneg == 0) { return allcovered; }

            IntVect const glo = bounding_box.smallEnd();
            IntVect const ghi = bounding_box.bigEnd();
            auto const invd2 = IF_invdist2(dx);
            auto const size = IF_size(problo, geom.ProbHiArray());
            Real const tolmax = IF_tolmax(amrex::get<2>(hv), amrex::get<3>(hv), invd2, size);
            if (amrex::get<3>(hv) < -tolmax) { return allregular; }
            ReduceOps<ReduceOpSum> reduce_op2;
            ReduceData<int> reduce_data2(reduce_op2);
            using ReduceTuple2 = typename decltype(reduce_data2)::Type;
            reduce_op2.eval(bx, reduce_data2,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple2
            {
                amrex::ignore_unused(j,k);
                IntVect const iv(AMREX_D_DECL(i,j,k));
                Real const v = fc(iv);
                bool const body = std::abs(v) <= tolmax &&
                    std::abs(v) <= IF_roundoff(f, iv, v, fc, problo, dx, invd2, size, glo, ghi);
                return {int(body)};
            });
            int const nbody = amrex::get<0>(reduce_data2.value(reduce_op2));
            return (nbody > 0) ? mixedcells : allregular;
        }
        else
        {
            return getBoxType_Cpu(bx, geom, bounding_box);
        }
    }

    //! CPU-only fallback for getBoxType().
    template <class U=F>
    requires (!IsGPUable<U>::value)
    [[nodiscard]] int getBoxType (const Box& bx, const Geometry& geom, RunOn,
                                  Box const& bounding_box) const noexcept
    {
        return getBoxType_Cpu(bx, geom, bounding_box);
    }

    template <class U=F>
    requires (IsGPUable<U>::value)
    static constexpr bool isGPUable () noexcept { return true; }

    template <class U=F>
    requires (!IsGPUable<U>::value)
    static constexpr bool isGPUable () noexcept { return false; }

    /**
     * \brief Fill a level-set FAB by sampling the implicit function on each point in the FAB.
     *
     * \param levelset Destination BaseFab.
     * \param geom     Geometry providing location metadata.
     * \param run_on   Execution target.
     * \param bounding_box Bounding box limiting the sampling region.
     */
    template <class U=F>
    requires (IsGPUable<U>::value)
    void fillFab (BaseFab<Real>& levelset, const Geometry& geom, RunOn run_on,
                  Box const& bounding_box) const noexcept
    {
        const auto& problo = geom.ProbLoArray();
        const auto& dx = geom.CellSizeArray();
        const Box& bx = levelset.box();
        // One more node layer for the neighbors in the roundoff tolerance.
        const Box gbx = amrex::grow(bx, 1);
        // Nodes outside the bounding box take the values of the nearest node in it.
        const Box ibx = bx & bounding_box;
        const Box ibx2 = gbx & bounding_box;
        IntVect const bloiv = bounding_box.smallEnd();
        IntVect const bhiiv = bounding_box.bigEnd();
        const auto& f = m_f;
        bool const on_gpu = (run_on == RunOn::Gpu) && Gpu::inLaunchRegion();
        Arena* const arena = on_gpu ? The_Async_Arena() : The_Cpu_Arena();
        BaseFab<Real> ffab(gbx, 1, arena);
        auto const& fv = ffab.array();
        auto const& fc = ffab.const_array();
        // Range of the sampled values, for IF_tolmax, without a host sync.
        auto* const pmm = static_cast<Real*>(arena->alloc(2*sizeof(Real)));
#ifdef AMREX_USE_GPU
        if (on_gpu) {
            amrex::single_task([=] AMREX_GPU_DEVICE () noexcept
            {
                pmm[0] = std::numeric_limits<Real>::max();
                pmm[1] = std::numeric_limits<Real>::lowest();
            });
            amrex::ParallelFor(Gpu::KernelInfo().setReduction(true), ibx2,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, Gpu::Handler const& h) noexcept
            {
                amrex::ignore_unused(j,k);
                Real const v = f(AMREX_D_DECL(node_coord(problo[0], i, dx[0]),
                                              node_coord(problo[1], j, dx[1]),
                                              node_coord(problo[2], k, dx[2])));
                fv(i,j,k) = v;
                Gpu::deviceReduceMin(pmm, v, h);
                Gpu::deviceReduceMax(pmm+1, v, h);
            });
        } else
#endif
        {
            Real vmin = std::numeric_limits<Real>::max();
            Real vmax = std::numeric_limits<Real>::lowest();
            amrex::LoopOnCpu(ibx2, [&] (int i, int j, int k) noexcept
            {
                amrex::ignore_unused(j,k);
                Real const v = f(AMREX_D_DECL(node_coord(problo[0], i, dx[0]),
                                              node_coord(problo[1], j, dx[1]),
                                              node_coord(problo[2], k, dx[2])));
                fv(i,j,k) = v;
                vmin = amrex::min(vmin, v);
                vmax = amrex::max(vmax, v);
            });
            pmm[0] = vmin;
            pmm[1] = vmax;
        }
        if (ibx2 != gbx) {
            AMREX_HOST_DEVICE_FOR_3D_FLAG(run_on, gbx, i, j, k,
            {
                amrex::ignore_unused(k);
                IntVect const iv(AMREX_D_DECL(i,j,k));
                if (!ibx2.contains(iv)) {
                    IntVect ivc = iv;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        ivc[d] = amrex::Clamp(iv[d], bloiv[d], bhiiv[d]);
                    }
                    fv(iv) = fv(ivc);
                }
            });
        }

        // Nodes on the surface within roundoff get a zero level set.
        auto const& s = levelset.array();
        auto const invd2 = IF_invdist2(dx);
        auto const size = IF_size(problo, geom.ProbHiArray());
        AMREX_HOST_DEVICE_FOR_3D_FLAG(run_on, ibx, i, j, k,
        {
            IntVect const iv(AMREX_D_DECL(i,j,k));
            Real const v = fc(i,j,k);
            Real const tolmax = IF_tolmax(pmm[0], pmm[1], invd2, size);
            s(i,j,k) = (v != 0.0_rt && std::abs(v) <= tolmax &&
                        std::abs(v) <= IF_roundoff(f, iv, v, fc, problo, dx, invd2, size,
                                                   bloiv, bhiiv))
                ? 0.0_rt : v;
        });
        arena->free(pmm);
        if (ibx != bx) {
            AMREX_HOST_DEVICE_FOR_3D_FLAG(run_on, bx, i, j, k,
            {
                amrex::ignore_unused(k);
                IntVect const iv(AMREX_D_DECL(i,j,k));
                if (!ibx.contains(iv)) {
                    IntVect ivc = iv;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        ivc[d] = amrex::Clamp(iv[d], bloiv[d], bhiiv[d]);
                    }
                    s(iv) = s(ivc);
                }
            });
        }
    }

    //! CPU-only overload of fillFab().
    template <class U=F>
    requires (!IsGPUable<U>::value)
    void fillFab (BaseFab<Real>& levelset, const Geometry& geom, RunOn,
                  Box const& bounding_box) const noexcept
    {
#ifdef AMREX_USE_GPU
        if (!levelset.arena()->isHostAccessible()) {
            const Box& bx = levelset.box();
            BaseFab<Real> h_levelset(bx, levelset.nComp(), The_Pinned_Arena());
            fillFab_Cpu(h_levelset, geom, bounding_box);
            Gpu::htod_memcpy_async(levelset.dataPtr(), h_levelset.dataPtr(),
                                   levelset.nBytes(bx, levelset.nComp()));
            Gpu::streamSynchronize();
        }
        else
#endif
        {
            fillFab_Cpu(levelset, geom, bounding_box);
        }
    }

    void fillFab_Cpu (BaseFab<Real>& levelset, const Geometry& geom,
                      Box const& bounding_box) const noexcept
    {
        const auto& problo = geom.ProbLoArray();
        const auto& dx = geom.CellSizeArray();
        const Box& bx = levelset.box();
        // See the GPU version.
        const Box gbx = amrex::grow(bx, 1);
        const Box ibx = bx & bounding_box;
        const Box ibx2 = gbx & bounding_box;
        IntVect const bloiv = bounding_box.smallEnd();
        IntVect const bhiiv = bounding_box.bigEnd();

        BaseFab<Real> ffab(gbx, 1, The_Cpu_Arena());
        const auto& fv = ffab.array();
        const auto& fc = ffab.const_array();
        Real vmin = std::numeric_limits<Real>::max();
        Real vmax = std::numeric_limits<Real>::lowest();
        amrex::LoopOnCpu(ibx2, [&] (int i, int j, int k) noexcept
        {
            amrex::ignore_unused(j,k);
            Real const v = m_f(RealArray{AMREX_D_DECL(node_coord(problo[0], i, dx[0]),
                                                      node_coord(problo[1], j, dx[1]),
                                                      node_coord(problo[2], k, dx[2]))});
            fv(i,j,k) = v;
            vmin = amrex::min(vmin, v);
            vmax = amrex::max(vmax, v);
        });
        if (ibx2 != gbx) {
            amrex::LoopOnCpu(gbx, [&] (int i, int j, int k) noexcept
            {
                amrex::ignore_unused(k);
                IntVect const iv(AMREX_D_DECL(i,j,k));
                if (!ibx2.contains(iv)) {
                    IntVect ivc = iv;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        ivc[d] = amrex::Clamp(iv[d], bloiv[d], bhiiv[d]);
                    }
                    fv(iv) = fv(ivc);
                }
            });
        }

        // Nodes on the surface within roundoff get a zero level set.
        auto const& s = levelset.array();
        auto const invd2 = IF_invdist2(dx);
        auto const size = IF_size(problo, geom.ProbHiArray());
        Real const tolmax = IF_tolmax(vmin, vmax, invd2, size);
        amrex::LoopOnCpu(ibx, [&] (int i, int j, int k) noexcept
        {
            IntVect const iv(AMREX_D_DECL(i,j,k));
            Real const v = fc(i,j,k);
            s(i,j,k) = (v != 0.0_rt && std::abs(v) <= tolmax &&
                        std::abs(v) <= IF_roundoff(m_f, iv, v, fc, problo, dx, invd2, size,
                                                   bloiv, bhiiv))
                ? 0.0_rt : v;
        });
        if (ibx != bx) {
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
            {
                amrex::ignore_unused(k);
                IntVect const iv(AMREX_D_DECL(i,j,k));
                if (!ibx.contains(iv)) {
                    IntVect ivc = iv;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        ivc[d] = amrex::Clamp(iv[d], bloiv[d], bhiiv[d]);
                    }
                    s(iv) = s(ivc);
                }
            });
        }
    }


    /**
     * \brief Compute EB intercepts for irregular edges (GPU-capable case).
     *
     * \param inter_arr Array of intercept FABs per direction.
     * \param type_arr  Edge-type arrays marking regular/covered/irregular edges.
     * \param lst       Level set; edges from a zero node are cut where the
     *                  surface leaves them.
     * \param geom      Geometry providing metrics.
     * \param run_on    Execution target.
     * \param bounding_box Bounding box limiting evaluation.
     */
    template <class U=F>
    requires (IsGPUable<U>::value)
    void getIntercept (Array<Array4<Real>,AMREX_SPACEDIM> const& inter_arr,
                       Array<Array4<Type_t const>,AMREX_SPACEDIM> const& type_arr,
                       Array4<Real const> const& lst, Geometry const& geom, RunOn run_on,
                       Box const& bounding_box) const noexcept
    {
        auto const& dx = geom.CellSizeArray();
        auto const& problo = geom.ProbLoArray();
        const auto& blo = amrex::lbound(bounding_box);
        const auto& bhi = amrex::ubound(bounding_box);
        IntVect const bloiv = bounding_box.smallEnd();
        IntVect const bhiiv = bounding_box.bigEnd();
        Real const size = IF_size(problo, geom.ProbHiArray());
        auto const& f = m_f;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            Array4<Real> const& inter = inter_arr[idim];
            Array4<Type_t const> const& type = type_arr[idim];
            const Box bx{inter};
            AMREX_HOST_DEVICE_FOR_3D_FLAG(run_on, bx, i, j, k,
            {
                if (type(i,j,k) == Type::irregular) {
                    IntVect ivlo(AMREX_D_DECL(i,j,k));
                    IntVect ivhi(AMREX_D_DECL(i,j,k));
                    ivhi[idim] += 1;
                    if (lst(ivlo) == 0.0_rt) {
                        inter(i,j,k) = ZeroEndIntercept(f, ivlo, idim, 1, lst, problo, dx, size,
                                                        bloiv, bhiiv);
                    } else if (lst(ivhi) == 0.0_rt) {
                        inter(i,j,k) = ZeroEndIntercept(f, ivhi, idim, -1, lst, problo, dx, size,
                                                        bloiv, bhiiv);
                    } else {
                        inter(i,j,k) = BrentRootFinder
                            ({AMREX_D_DECL(node_coord(problo[0], amrex::Clamp(ivlo[0],blo.x,bhi.x), dx[0]),
                                           node_coord(problo[1], amrex::Clamp(ivlo[1],blo.y,bhi.y), dx[1]),
                                           node_coord(problo[2], amrex::Clamp(ivlo[2],blo.z,bhi.z), dx[2]))},
                             {AMREX_D_DECL(node_coord(problo[0], amrex::Clamp(ivhi[0],blo.x,bhi.x), dx[0]),
                                           node_coord(problo[1], amrex::Clamp(ivhi[1],blo.y,bhi.y), dx[1]),
                                           node_coord(problo[2], amrex::Clamp(ivhi[2],blo.z,bhi.z), dx[2]))},
                                idim, f);
                    }
                } else {
                    inter(i,j,k) = no_intercept;
                }
            });
        }
    }

    //! CPU-only overload of getIntercept().
    template <class U=F>
    requires (!IsGPUable<U>::value)
    void getIntercept (Array<Array4<Real>,AMREX_SPACEDIM> const& inter_arr,
                       Array<Array4<Type_t const>,AMREX_SPACEDIM> const& type_arr,
                       Array4<Real const> const& lst, Geometry const& geom, RunOn,
                       Box const& bounding_box) const noexcept
    {
#ifdef AMREX_USE_GPU
        BaseFab<Real> h_lst_fab;
        Array4<Real const> h_lst = lst;
        if (!The_Arena()->isHostAccessible()) {
            const Box lbx{lst};
            h_lst_fab.resize(lbx, 1, The_Pinned_Arena());
            Gpu::dtoh_memcpy_async(h_lst_fab.dataPtr(), lst.dataPtr(), h_lst_fab.nBytes(lbx, 1));
            Gpu::streamSynchronize();
            h_lst = h_lst_fab.const_array();
        }
#endif
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            Array4<Real> const& inter = inter_arr[idim];
            Array4<Type_t const> const& type = type_arr[idim];
            // When GShop is not GPUable, the intercept is either in managed memory or
            // in device memory and needs to be copied to host and setup there.
#ifdef AMREX_USE_GPU
            if (!The_Arena()->isHostAccessible()) {
                const Box bx{inter};
                BaseFab<Real>  h_inter_fab(bx, 1, The_Pinned_Arena());
                BaseFab<Type_t> h_type_fab(bx, 1, The_Pinned_Arena());
                Gpu::dtoh_memcpy_async(h_type_fab.dataPtr(), type.dataPtr(),
                                       h_type_fab.nBytes(bx, 1));
                Gpu::streamSynchronize();
                const auto& h_inter = h_inter_fab.array();
                const auto& h_type  = h_type_fab.const_array();
                getIntercept_Cpu(h_inter, h_type, h_lst, geom, bounding_box, idim);
                Gpu::htod_memcpy_async(inter.dataPtr(), h_inter.dataPtr(),
                                       h_inter_fab.nBytes(bx, 1));
                Gpu::streamSynchronize();
            }
            else
#endif
            {
                getIntercept_Cpu(inter, type, lst, geom, bounding_box, idim);
            }

        }
    }

    void getIntercept_Cpu (Array4<Real> const& inter,
                           Array4<Type_t const> const& type,
                           Array4<Real const> const& lst,
                           Geometry const& geom,
                           Box const& bounding_box,
                           const int idim) const noexcept
    {
        auto const& dx = geom.CellSizeArray();
        auto const& problo = geom.ProbLoArray();
        const auto& blo = amrex::lbound(bounding_box);
        const auto& bhi = amrex::ubound(bounding_box);
        IntVect const bloiv = bounding_box.smallEnd();
        IntVect const bhiiv = bounding_box.bigEnd();
        Real const size = IF_size(problo, geom.ProbHiArray());
        const Box bx{inter};
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            if (type(i,j,k) == Type::irregular) {
                IntVect ivlo(AMREX_D_DECL(i,j,k));
                IntVect ivhi(AMREX_D_DECL(i,j,k));
                ivhi[idim] += 1;
                if (lst(ivlo) == 0.0_rt) {
                    inter(i,j,k) = ZeroEndIntercept(m_f, ivlo, idim, 1, lst, problo, dx, size,
                                                    bloiv, bhiiv);
                } else if (lst(ivhi) == 0.0_rt) {
                    inter(i,j,k) = ZeroEndIntercept(m_f, ivhi, idim, -1, lst, problo, dx, size,
                                                    bloiv, bhiiv);
                } else {
                    inter(i,j,k) = BrentRootFinder
                        ({AMREX_D_DECL(node_coord(problo[0], amrex::Clamp(ivlo[0],blo.x,bhi.x), dx[0]),
                                       node_coord(problo[1], amrex::Clamp(ivlo[1],blo.y,bhi.y), dx[1]),
                                       node_coord(problo[2], amrex::Clamp(ivlo[2],blo.z,bhi.z), dx[2]))},
                         {AMREX_D_DECL(node_coord(problo[0], amrex::Clamp(ivhi[0],blo.x,bhi.x), dx[0]),
                                       node_coord(problo[1], amrex::Clamp(ivhi[1],blo.y,bhi.y), dx[1]),
                                       node_coord(problo[2], amrex::Clamp(ivhi[2],blo.z,bhi.z), dx[2]))},
                         idim, m_f);
                }
            } else {
                inter(i,j,k) = no_intercept;
            }
        });
    }

    /**
     * \brief Update intercepts when level-set data has been adjusted (e.g., fixing small cells).
     *
     * \param inter_arr Array of intercept FABs per direction.
     * \param type_arr  Edge-type arrays marking regular/covered/irregular edges.
     * \param lst       Level-set data used to determine zero crossings.
     * \param geom      Geometry providing metrics.
     */
    void updateIntercept (Array<Array4<Real>,AMREX_SPACEDIM> const& inter_arr,
                          Array<Array4<Type_t const>,AMREX_SPACEDIM> const& type_arr,
                          Array4<Real const> const& lst, Geometry const& geom) const noexcept
    {
        auto const& dx = geom.CellSizeArray();
        auto const& problo = geom.ProbLoArray();
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            Array4<Real> const& inter = inter_arr[idim];
            Array4<Type_t const> const& type = type_arr[idim];
            const Box bx{inter};
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
            {
                if (type(i,j,k) == Type::irregular) {
                    if (inter(i,j,k) == no_intercept) {
                        // Only edges without an intercept are updated: they
                        // became irregular because small cell fixing set lst
                        // to zero, or lst changed in FillBoundary at periodic
                        // boundaries.  A node that was zero from the start
                        // keeps the intercept from getIntercept.  Other
                        // edges are reset, so an edge that is irregular
                        // again gets a new intercept.
                        IntVect ivlo(AMREX_D_DECL(i,j,k));
                        IntVect ivhi = ivlo;
                        ivhi[idim] += 1;
                        // An irregular edge has one end below zero and one at or
                        // above.  The intercept is at the body end.
                        bool const at_lo = lst(ivlo) >= Real(0.0);
                        int const inode = at_lo ? ivlo[idim] : ivhi[idim];
                        inter(i,j,k) = node_coord(problo[idim], inode, dx[idim]);
                    }
                } else {
                    inter(i,j,k) = no_intercept;
                }
            });
        }
    }

private:

    F m_f;
    R m_resource;  // We use this to hold the ownership of resource for F if needed,
                   // because F needs to be a simply type suitable for GPU.
};

/**
 * \brief Helper that constructs a GeometryShop owning the supplied implicit function.
 *
 * \param f Implicit function callable.
 */
template <class F>
GeometryShop<std::decay_t<F>>
makeShop (F&& f)
{
    return GeometryShop<std::decay_t<F>>(std::forward<F>(f));
}

/**
 * \brief Helper that constructs a GeometryShop for an implicit function plus an extra
 *        resource handle.
 *
 * \param f Implicit function callable.
 * \param r Extra resource (e.g., Parser data) kept alive for \p f.
 */
template <class F, class R>
GeometryShop<std::decay_t<F>, std::decay_t<R>>
makeShop (F&& f, R&& r)
{
    return GeometryShop<std::decay_t<F>,
                        std::decay_t<R>>
        (std::forward<F>(f), std::forward<R>(r));
}

}

#endif
