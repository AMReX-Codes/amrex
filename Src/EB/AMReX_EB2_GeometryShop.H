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

// Tolerance below which f at node iv is roundoff: 16 epsilon times the
// change of f to the neighboring nodes times the size of the operands of
// the node position.  Nodes are placed as in fillFab.  Neighbor values are
// taken from fab where it has them.
template <class F>
AMREX_GPU_HOST_DEVICE
Real
IF_roundoff (F const& f, IntVect const& iv, Real fv, Array4<Real const> const& fab,
             GpuArray<Real,AMREX_SPACEDIM> const& problo,
             GpuArray<Real,AMREX_SPACEDIM> const& dx,
             IntVect const& blo, IntVect const& bhi)
{
    // Outside the bounding box, nodes take the values at the clamped node,
    // so they must get its tolerance too.
    IntVect ivc = iv;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { ivc[d] = amrex::Clamp(iv[d], blo[d], bhi[d]); }
    Real tol = 0.0_rt;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        Real df = 0.0_rt;
        for (int s = -1; s <= 1; s += 2) {
            IntVect ivn = ivc;
            ivn[d] += s;
            Real fn;
            if (fab.contains(ivn)) {
                fn = fab(ivn);
            } else {
                GpuArray<Real,AMREX_SPACEDIM> p;
                for (int dd = 0; dd < AMREX_SPACEDIM; ++dd) {
                    p[dd] = problo[dd] + Real(amrex::Clamp(ivn[dd],blo[dd],bhi[dd]))*dx[dd];
                }
                fn = IF_f(f, p);
            }
            df = amrex::max(df, std::abs(fn - fv));
        }
        tol += df / dx[d] * (std::abs(problo[d]) + std::abs(Real(ivc[d])*dx[d]));
    }
    return Real(16.0) * std::numeric_limits<Real>::epsilon() * tol;
}

// Intercept on the edge from node ivz, which is on the surface (zero level
// set), to the fluid node one cell away in direction rangedir.  This is
// where f leaves the surface, which is not ivz if the edge runs along it.
template <class F>
AMREX_GPU_HOST_DEVICE
Real
ZeroEndIntercept (F const& f, IntVect const& ivz, int rangedir, int side,
                  GpuArray<Real,AMREX_SPACEDIM> const& problo,
                  GpuArray<Real,AMREX_SPACEDIM> const& dx,
                  IntVect const& blo, IntVect const& bhi)
{
    GpuArray<Real,AMREX_SPACEDIM> p;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        p[d] = problo[d] + Real(amrex::Clamp(ivz[d],blo[d],bhi[d]))*dx[d];
    }
    Real const tol = IF_roundoff(f, ivz, IF_f(f, p), Array4<Real const>{},
                                 problo, dx, blo, bhi);
    IntVect ivf = ivz;
    ivf[rangedir] += side;
    Real a = p[rangedir];
    Real b = problo[rangedir] + Real(amrex::Clamp(ivf[rangedir],blo[rangedir],bhi[rangedir]))
        * dx[rangedir];
    Real const xtol = Real(4.0) * std::numeric_limits<Real>::epsilon()
        * amrex::max(std::abs(a), std::abs(b), std::abs(problo[rangedir]));
    Real const x0 = a;
    p[rangedir] = b;
    Real const fb = IF_f(f, p);
    Real const h = std::abs(b - a);
    for (int iter = 0; iter < 200 && std::abs(b-a) > xtol; ++iter) {
        p[rangedir] = Real(0.5)*(a+b);
        if (IF_f(f, p) >= -tol) {
            a = p[rangedir];
        } else {
            b = p[rangedir];
        }
    }
    // Within roundoff of the node, the surface just passes through it.
    // Roundoff here is where f drops below -tol at the slope along the edge.
    Real const node_tol = (fb != 0.0_rt) ? Real(2.0)*tol*h/std::abs(fb) : 0.0_rt;
    return (std::abs(a - x0) <= node_tol) ? x0 : a;
}

//! Level set of a node on the surface that is placed on the fluid side
inline constexpr Real on_surface_fluid = -std::numeric_limits<Real>::min();

// A node on the surface (zero level set) whose diagonal neighbor on a face is
// body while the other two corners are fluid makes the face look cut twice,
// although the body only touches the face along a line.  Such a node is put
// on the fluid side, unless that makes the opposite pattern on another face.
// Of two diagonal nodes of a face whose other corners are body, only one
// moves, or the face would be cut twice.
inline void ResolveZeroSaddles (BaseFab<Real>& levelset, RunOn run_on)
{
    const Box& bx = levelset.box();
    bool const on_gpu = (run_on == RunOn::Gpu) && Gpu::inLaunchRegion();
    BaseFab<Real> tmp(bx, 1, on_gpu ? The_Async_Arena() : The_Cpu_Arena());
    auto const& t = tmp.array();
    auto const& tc = tmp.const_array();
    auto const& a = levelset.array();
    AMREX_HOST_DEVICE_FOR_3D_FLAG(run_on, bx, i, j, k,
    {
        t(i,j,k) = a(i,j,k);
    });
    AMREX_HOST_DEVICE_FOR_3D_FLAG(run_on, bx, i, j, k,
    {
        if (tc(i,j,k) == 0.0_rt) {
            IntVect const iv(AMREX_D_DECL(i,j,k));
            bool body_saddle = false;
            bool fluid_saddle = false;
            for (int d = 0; d < 3; ++d) { // face normal
                int const u = (d+1)%3;
                int const v = (d+2)%3;
                if (u >= AMREX_SPACEDIM || v >= AMREX_SPACEDIM) { continue; }
                for (int su = -1; su <= 1; su += 2) {
                    for (int sv = -1; sv <= 1; sv += 2) {
                        IntVect ib = iv;
                        IntVect ie = iv;
                        ib[u] += su;
                        ie[v] += sv;
                        IntVect ic = ib;
                        ic[v] += sv;
                        if (tc.contains(ib) && tc.contains(ie) && tc.contains(ic)) {
                            bool const cb = tc(ic) >= 0.0_rt;
                            bool const bb = tc(ib) >= 0.0_rt;
                            bool const eb = tc(ie) >= 0.0_rt;
                            body_saddle = body_saddle || (cb && !bb && !eb);
                            fluid_saddle = fluid_saddle || (!cb && bb && eb);
                        }
                    }
                }
            }
            if (body_saddle && !fluid_saddle) { a(i,j,k) = on_surface_fluid; }
        }
    });
    // a now has the candidates; the decisions go to t.
    AMREX_HOST_DEVICE_FOR_3D_FLAG(run_on, bx, i, j, k,
    {
        if (a(i,j,k) == on_surface_fluid) {
            IntVect const iv(AMREX_D_DECL(i,j,k));
            bool move = true;
            for (int d = 0; d < 3; ++d) {
                int const u = (d+1)%3;
                int const v = (d+2)%3;
                if (u >= AMREX_SPACEDIM || v >= AMREX_SPACEDIM) { continue; }
                for (int su = -1; su <= 1; su += 2) {
                    for (int sv = -1; sv <= 1; sv += 2) {
                        IntVect ib = iv;
                        IntVect ie = iv;
                        ib[u] += su;
                        ie[v] += sv;
                        IntVect ic = ib;
                        ic[v] += sv;
                        if (a.contains(ib) && a.contains(ie) && a.contains(ic)) {
                            bool const bb = a(ib) >= 0.0_rt || a(ib) == on_surface_fluid;
                            bool const eb = a(ie) >= 0.0_rt || a(ie) == on_surface_fluid;
                            // The one with the lower index moves.
                            if (a(ic) == on_surface_fluid && bb && eb &&
                                ic[amrex::max(u,v)] < iv[amrex::max(u,v)]) {
                                move = false;
                            }
                        }
                    }
                }
            }
            t(i,j,k) = move ? on_surface_fluid : 0.0_rt;
        }
    });
    AMREX_HOST_DEVICE_FOR_3D_FLAG(run_on, bx, i, j, k,
    {
        if (a(i,j,k) == on_surface_fluid) { a(i,j,k) = t(i,j,k); }
    });
    if (on_gpu) { Gpu::streamSynchronize(); }
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

    [[nodiscard]] int getBoxType_Cpu (const Box& bx, Geometry const& geom) const noexcept
    {
        const Real* problo = geom.ProbLo();
        const Real* dx = geom.CellSize();
        const auto& len3 = bx.length3d();
        const int* blo = bx.loVect();
        int nbody = 0, nzero = 0, nfluid = 0;
        for         (int k = 0; k < len3[2]; ++k) {
            for     (int j = 0; j < len3[1]; ++j) {
                for (int i = 0; i < len3[0]; ++i) {
                    RealArray xyz {AMREX_D_DECL(problo[0]+(i+blo[0])*dx[0],
                                                problo[1]+(j+blo[1])*dx[1],
                                                problo[2]+(k+blo[2])*dx[2])};
                    Real v = m_f(xyz);
                    if (v == 0.0_rt) {
                        ++nzero;
                    } else if (v > 0.0_rt) {
                        ++nbody;
                    } else {
                        ++nfluid;
                    }
                    if (nbody > 0 && nfluid > 0) { return mixedcells; }
                }
            }
        }
        amrex::ignore_unused(nzero);

        if (nbody == 0) {
            return allregular;
        } else if (nfluid == 0) {
            return allcovered;
        } else {
            return mixedcells;
        }
    }

    /**
     * \brief Classify \p bx as regular, covered, or mixed (GPU-capable case).
     *
     * \param bx   Box to inspect.
     * \param geom Geometry providing physical coordinates.
     * \param run_on Whether to run on CPU or GPU.
     */
    template <class U=F>
    requires (IsGPUable<U>::value)
    [[nodiscard]] int getBoxType (const Box& bx, const Geometry& geom, RunOn run_on) const noexcept
    {
        if (run_on == RunOn::Gpu && Gpu::inLaunchRegion())
        {
            const auto& problo = geom.ProbLoArray();
            const auto& dx = geom.CellSizeArray();
            auto const& f = m_f;
            ReduceOps<ReduceOpSum,ReduceOpSum,ReduceOpSum> reduce_op;
            ReduceData<int,int,int> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                amrex::ignore_unused(j,k);
                Real v = f(AMREX_D_DECL(problo[0]+Real(i)*dx[0],
                                        problo[1]+Real(j)*dx[1],
                                        problo[2]+Real(k)*dx[2]));
                int nbody = 0, nzero = 0, nfluid = 0;
                if (v == 0.0_rt) {
                    ++nzero;
                } else if (v > 0.0_rt) {
                    ++nbody;
                } else {
                    ++nfluid;
                }
                return {nbody,nzero,nfluid};
            });
            ReduceTuple hv = reduce_data.value(reduce_op);
            int nbody = amrex::get<0>(hv);
            // int nzero = amrex::get<1>(hv);
            int nfluid = amrex::get<2>(hv);

            if (nbody == 0) {
                return allregular;
            } else if (nfluid == 0) {
                return allcovered;
            } else {
                return mixedcells;
            }
        }
        else
        {
            return getBoxType_Cpu(bx, geom);
        }
    }

    //! CPU-only fallback for getBoxType().
    template <class U=F>
    requires (!IsGPUable<U>::value)
    [[nodiscard]] int getBoxType (const Box& bx, const Geometry& geom, RunOn) const noexcept
    {
        return getBoxType_Cpu(bx, geom);
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
        const auto& a = levelset.array();
        const auto& blo = amrex::lbound(bounding_box);
        const auto& bhi = amrex::ubound(bounding_box);
        const auto& f = m_f;
        AMREX_HOST_DEVICE_FOR_3D_FLAG(run_on, bx, i, j, k,
        {
            a(i,j,k) = f(AMREX_D_DECL(problo[0]+Real(amrex::Clamp(i,blo.x,bhi.x))*dx[0],
                                      problo[1]+Real(amrex::Clamp(j,blo.y,bhi.y))*dx[1],
                                      problo[2]+Real(amrex::Clamp(k,blo.z,bhi.z))*dx[2]));
        });

        // Nodes on the surface within roundoff get a zero level set.
        bool const on_gpu = (run_on == RunOn::Gpu) && Gpu::inLaunchRegion();
        BaseFab<Real> tmp(bx, 1, on_gpu ? The_Async_Arena() : The_Cpu_Arena());
        auto const& t = tmp.array();
        auto const& tc = tmp.const_array();
        AMREX_HOST_DEVICE_FOR_3D_FLAG(run_on, bx, i, j, k,
        {
            t(i,j,k) = a(i,j,k);
        });
        IntVect const bloiv = bounding_box.smallEnd();
        IntVect const bhiiv = bounding_box.bigEnd();
        AMREX_HOST_DEVICE_FOR_3D_FLAG(run_on, bx, i, j, k,
        {
            IntVect const iv(AMREX_D_DECL(i,j,k));
            Real const v = tc(iv);
            if (v != 0.0_rt &&
                std::abs(v) <= IF_roundoff(f, iv, v, tc, problo, dx, bloiv, bhiiv)) {
                a(i,j,k) = 0.0_rt;
            }
        });
        if (on_gpu) { Gpu::streamSynchronize(); }
        ResolveZeroSaddles(levelset, run_on);
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
        const auto& blo = amrex::lbound(bounding_box);
        const auto& bhi = amrex::ubound(bounding_box);

        const auto& a = levelset.array();
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            a(i,j,k) = m_f(RealArray{AMREX_D_DECL(problo[0]+Real(amrex::Clamp(i,blo.x,bhi.x))*dx[0],
                                                  problo[1]+Real(amrex::Clamp(j,blo.y,bhi.y))*dx[1],
                                                  problo[2]+Real(amrex::Clamp(k,blo.z,bhi.z))*dx[2])});
        });

        // Nodes on the surface within roundoff get a zero level set.
        BaseFab<Real> tmp(bx, 1, The_Cpu_Arena());
        tmp.copy<RunOn::Host>(levelset);
        auto const& tc = tmp.const_array();
        IntVect const bloiv = bounding_box.smallEnd();
        IntVect const bhiiv = bounding_box.bigEnd();
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            IntVect const iv(AMREX_D_DECL(i,j,k));
            Real const v = tc(iv);
            if (v != 0.0_rt &&
                std::abs(v) <= IF_roundoff(m_f, iv, v, tc, problo, dx, bloiv, bhiiv)) {
                a(i,j,k) = 0.0_rt;
            }
        });
        ResolveZeroSaddles(levelset, RunOn::Cpu);
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
                    if (lst(ivlo) == on_surface_fluid) {
                        inter(i,j,k) = problo[idim]
                            + Real(amrex::Clamp(ivlo[idim],bloiv[idim],bhiiv[idim]))*dx[idim];
                    } else if (lst(ivhi) == on_surface_fluid) {
                        inter(i,j,k) = problo[idim]
                            + Real(amrex::Clamp(ivhi[idim],bloiv[idim],bhiiv[idim]))*dx[idim];
                    } else if (lst(ivlo) == 0.0_rt) {
                        inter(i,j,k) = ZeroEndIntercept(f, ivlo, idim, 1, problo, dx,
                                                        bloiv, bhiiv);
                    } else if (lst(ivhi) == 0.0_rt) {
                        inter(i,j,k) = ZeroEndIntercept(f, ivhi, idim, -1, problo, dx,
                                                        bloiv, bhiiv);
                    } else {
                        inter(i,j,k) = BrentRootFinder
                            ({AMREX_D_DECL(problo[0]+Real(amrex::Clamp(ivlo[0],blo.x,bhi.x))*dx[0],
                                           problo[1]+Real(amrex::Clamp(ivlo[1],blo.y,bhi.y))*dx[1],
                                           problo[2]+Real(amrex::Clamp(ivlo[2],blo.z,bhi.z))*dx[2])},
                             {AMREX_D_DECL(problo[0]+Real(amrex::Clamp(ivhi[0],blo.x,bhi.x))*dx[0],
                                           problo[1]+Real(amrex::Clamp(ivhi[1],blo.y,bhi.y))*dx[1],
                                           problo[2]+Real(amrex::Clamp(ivhi[2],blo.z,bhi.z))*dx[2])},
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
        const Box bx{inter};
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            if (type(i,j,k) == Type::irregular) {
                IntVect ivlo(AMREX_D_DECL(i,j,k));
                IntVect ivhi(AMREX_D_DECL(i,j,k));
                ivhi[idim] += 1;
                if (lst(ivlo) == on_surface_fluid) {
                    inter(i,j,k) = problo[idim]
                        + Real(amrex::Clamp(ivlo[idim],bloiv[idim],bhiiv[idim]))*dx[idim];
                } else if (lst(ivhi) == on_surface_fluid) {
                    inter(i,j,k) = problo[idim]
                        + Real(amrex::Clamp(ivhi[idim],bloiv[idim],bhiiv[idim]))*dx[idim];
                } else if (lst(ivlo) == 0.0_rt) {
                    inter(i,j,k) = ZeroEndIntercept(m_f, ivlo, idim, 1, problo, dx,
                                                    bloiv, bhiiv);
                } else if (lst(ivhi) == 0.0_rt) {
                    inter(i,j,k) = ZeroEndIntercept(m_f, ivhi, idim, -1, problo, dx,
                                                    bloiv, bhiiv);
                } else {
                    inter(i,j,k) = BrentRootFinder
                        ({AMREX_D_DECL(problo[0]+Real(amrex::Clamp(ivlo[0],blo.x,bhi.x))*dx[0],
                                       problo[1]+Real(amrex::Clamp(ivlo[1],blo.y,bhi.y))*dx[1],
                                       problo[2]+Real(amrex::Clamp(ivlo[2],blo.z,bhi.z))*dx[2])},
                         {AMREX_D_DECL(problo[0]+Real(amrex::Clamp(ivhi[0],blo.x,bhi.x))*dx[0],
                                       problo[1]+Real(amrex::Clamp(ivhi[1],blo.y,bhi.y))*dx[1],
                                       problo[2]+Real(amrex::Clamp(ivhi[2],blo.z,bhi.z))*dx[2])},
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
                    bool no_inter = (inter(i,j,k) == no_intercept);
                    if (idim == 0) {
                        if (no_inter && lst(i,j,k) >= Real(0.0))
                        {
                            // Only edges without an intercept are updated: they
                            // became irregular because small cell fixing set lst
                            // to zero, or lst changed in FillBoundary at periodic
                            // boundaries.  A node that was zero from the start
                            // keeps the intercept from getIntercept.  Other
                            // edges are reset, so an edge that is irregular
                            // again gets a new intercept.
                            inter(i,j,k) = problo[0] + Real(i)*dx[0];
                        }
                        else if (no_inter && lst(i+1,j,k) >= Real(0.0))
                        {
                            inter(i,j,k) = problo[0] + Real(i+1)*dx[0];
                        }
                    } else if (idim == 1) {
                        if (no_inter && lst(i,j,k) >= Real(0.0))
                        {
                            inter(i,j,k) = problo[1] + Real(j)*dx[1];
                        }
                        else if (no_inter && lst(i,j+1,k) >= Real(0.0))
                        {
                            inter(i,j,k) = problo[1] + Real(j+1)*dx[1];
                        }
                    } else {
                        if (no_inter && lst(i,j,k) >= Real(0.0))
                        {
                            inter(i,j,k) = problo[2] + Real(k)*dx[2];
                        }
                        else if (no_inter && lst(i,j,k+1) >= Real(0.0))
                        {
                            inter(i,j,k) = problo[2] + Real(k+1)*dx[2];
                        }
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
