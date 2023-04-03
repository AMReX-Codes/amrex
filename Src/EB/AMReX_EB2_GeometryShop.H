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
#include <utility>

namespace amrex::EB2 {

template <class F, typename std::enable_if<IsGPUable<F>::value>::type* FOO = nullptr>
AMREX_GPU_HOST_DEVICE
Real
IF_f (F const& f, GpuArray<Real,AMREX_SPACEDIM> const& p) noexcept
{
    return f(AMREX_D_DECL(p[0],p[1],p[2]));
}

template <class F, typename std::enable_if<!IsGPUable<F>::value>::type* BAR = nullptr>
AMREX_GPU_HOST_DEVICE
Real
IF_f (F const& f, GpuArray<Real,AMREX_SPACEDIM> const& p) noexcept
{
#if AMREX_DEVICE_COMPILE
    amrex::ignore_unused(f,p);
    amrex::Error("EB2::GeometryShop: how did this happen?");
    return 0.0;
#else
    return f({AMREX_D_DECL(p[0],p[1],p[2])});
#endif
}

template <class F>
AMREX_GPU_HOST_DEVICE
Real
BrentRootFinder (GpuArray<Real,AMREX_SPACEDIM> const& lo,
                 GpuArray<Real,AMREX_SPACEDIM> const& hi,
                 int rangedir, F const& f) noexcept
{
#ifdef AMREX_USE_FLOAT
    const Real tol = 1.e-4_rt;
    const Real EPS = 1.0e-6_rt;
#else
    const Real tol = 1.e-12;
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
            if (p > 0) q = -q;

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
            if (xm < 0) bPt[rangedir] = bPt[rangedir] - tol1;
            else        bPt[rangedir] = bPt[rangedir] + tol1;
        }

        fb = IF_f(f, bPt);
    }

    if (i >= MAXITER)
    {
        amrex::Error("BrentRootFinder: exceeding maximum iterations.");
    }

    return bPt[rangedir];
}

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

    explicit GeometryShop (F f)
        : m_f(std::move(f)), m_resource()
        {}

    GeometryShop (F f, R r)
        : m_f(std::move(f)), m_resource(std::move(r))
        {}

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
                    if (nbody > 0 && nfluid > 0) return mixedcells;
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

    template <class U=F, typename std::enable_if<IsGPUable<U>::value>::type* FOO = nullptr >
    [[nodiscard]] int getBoxType (const Box& bx, const Geometry& geom, RunOn run_on) const noexcept
    {
        if (run_on == RunOn::Gpu && Gpu::inLaunchRegion())
        {
            const auto& problo = geom.ProbLoArray();
            const auto& dx = geom.CellSizeArray();
            auto f = m_f;
            ReduceOps<ReduceOpSum,ReduceOpSum,ReduceOpSum> reduce_op;
            ReduceData<int,int,int> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                amrex::ignore_unused(j,k);
                Real v = f(AMREX_D_DECL(problo[0]+i*dx[0],
                                        problo[1]+j*dx[1],
                                        problo[2]+k*dx[2]));
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

    template <class U=F, typename std::enable_if<!IsGPUable<U>::value>::type* BAR = nullptr >
    [[nodiscard]] int getBoxType (const Box& bx, const Geometry& geom, RunOn) const noexcept
    {
        return getBoxType_Cpu(bx, geom);
    }

    template <class U=F, typename std::enable_if<IsGPUable<U>::value>::type* FOO = nullptr >
    static constexpr bool isGPUable () noexcept { return true; }

    template <class U=F, typename std::enable_if<!IsGPUable<U>::value>::type* BAR = nullptr >
    static constexpr bool isGPUable () noexcept { return false; }

    template <class U=F, typename std::enable_if<IsGPUable<U>::value>::type* FOO = nullptr >
    void fillFab (BaseFab<Real>& levelset, const Geometry& geom, RunOn run_on,
                  Box const& bounding_box) const noexcept
    {
        const auto problo = geom.ProbLoArray();
        const auto dx = geom.CellSizeArray();
        const Box& bx = levelset.box();
        const auto& a = levelset.array();
        const auto blo = amrex::lbound(bounding_box);
        const auto bhi = amrex::ubound(bounding_box);
        auto f = m_f;
        AMREX_HOST_DEVICE_FOR_3D_FLAG(run_on, bx, i, j, k,
        {
            a(i,j,k) = f(AMREX_D_DECL(problo[0]+amrex::Clamp(i,blo.x,bhi.x)*dx[0],
                                      problo[1]+amrex::Clamp(j,blo.y,bhi.y)*dx[1],
                                      problo[2]+amrex::Clamp(k,blo.z,bhi.z)*dx[2]));
        });
    }

    template <class U=F, typename std::enable_if<!IsGPUable<U>::value>::type* BAR = nullptr >
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
        const auto problo = geom.ProbLoArray();
        const auto dx = geom.CellSizeArray();
        const Box& bx = levelset.box();
        const auto blo = amrex::lbound(bounding_box);
        const auto bhi = amrex::ubound(bounding_box);

        const auto& a = levelset.array();
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            a(i,j,k) = m_f(RealArray{AMREX_D_DECL(problo[0]+amrex::Clamp(i,blo.x,bhi.x)*dx[0],
                                                  problo[1]+amrex::Clamp(j,blo.y,bhi.y)*dx[1],
                                                  problo[2]+amrex::Clamp(k,blo.z,bhi.z)*dx[2])});
        });
    }


    template <class U=F, typename std::enable_if<IsGPUable<U>::value>::type* FOO = nullptr >
    void getIntercept (Array<Array4<Real>,AMREX_SPACEDIM> const& inter_arr,
                       Array<Array4<Type_t const>,AMREX_SPACEDIM> const& type_arr,
                       Array4<Real const> const&, Geometry const& geom, RunOn run_on,
                       Box const& bounding_box) const noexcept
    {
        auto const& dx = geom.CellSizeArray();
        auto const& problo = geom.ProbLoArray();
        const auto blo = amrex::lbound(bounding_box);
        const auto bhi = amrex::ubound(bounding_box);
        auto f = m_f;
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
                    inter(i,j,k) = BrentRootFinder
                        ({AMREX_D_DECL(problo[0]+amrex::Clamp(ivlo[0],blo.x,bhi.x)*dx[0],
                                       problo[1]+amrex::Clamp(ivlo[1],blo.y,bhi.y)*dx[1],
                                       problo[2]+amrex::Clamp(ivlo[2],blo.z,bhi.z)*dx[2])},
                         {AMREX_D_DECL(problo[0]+amrex::Clamp(ivhi[0],blo.x,bhi.x)*dx[0],
                                       problo[1]+amrex::Clamp(ivhi[1],blo.y,bhi.y)*dx[1],
                                       problo[2]+amrex::Clamp(ivhi[2],blo.z,bhi.z)*dx[2])},
                            idim, f);
                } else {
                    inter(i,j,k) = std::numeric_limits<Real>::quiet_NaN();
                }
            });
        }
    }

    template <class U=F, typename std::enable_if<!IsGPUable<U>::value>::type* BAR = nullptr >
    void getIntercept (Array<Array4<Real>,AMREX_SPACEDIM> const& inter_arr,
                       Array<Array4<Type_t const>,AMREX_SPACEDIM> const& type_arr,
                       Array4<Real const> const&, Geometry const& geom, RunOn,
                       Box const& bounding_box) const noexcept
    {
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
                getIntercept_Cpu(h_inter, h_type, geom, bounding_box, idim);
                Gpu::htod_memcpy_async(inter.dataPtr(), h_inter.dataPtr(),
                                       h_inter_fab.nBytes(bx, 1));
                Gpu::streamSynchronize();
            }
            else
#endif
            {
                getIntercept_Cpu(inter, type, geom, bounding_box, idim);
            }

        }
    }

    void getIntercept_Cpu (Array4<Real> const& inter,
                           Array4<Type_t const> const& type,
                           Geometry const& geom,
                           Box const& bounding_box,
                           const int idim) const noexcept
    {
        auto const& dx = geom.CellSizeArray();
        auto const& problo = geom.ProbLoArray();
        const auto blo = amrex::lbound(bounding_box);
        const auto bhi = amrex::ubound(bounding_box);
        const Box bx{inter};
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            if (type(i,j,k) == Type::irregular) {
                IntVect ivlo(AMREX_D_DECL(i,j,k));
                IntVect ivhi(AMREX_D_DECL(i,j,k));
                ivhi[idim] += 1;
                inter(i,j,k) = BrentRootFinder
                    ({AMREX_D_DECL(problo[0]+amrex::Clamp(ivlo[0],blo.x,bhi.x)*dx[0],
                                   problo[1]+amrex::Clamp(ivlo[1],blo.y,bhi.y)*dx[1],
                                   problo[2]+amrex::Clamp(ivlo[2],blo.z,bhi.z)*dx[2])},
                     {AMREX_D_DECL(problo[0]+amrex::Clamp(ivhi[0],blo.x,bhi.x)*dx[0],
                                   problo[1]+amrex::Clamp(ivhi[1],blo.y,bhi.y)*dx[1],
                                   problo[2]+amrex::Clamp(ivhi[2],blo.z,bhi.z)*dx[2])},
                     idim, m_f);
            } else {
                inter(i,j,k) = std::numeric_limits<Real>::quiet_NaN();
            }
        });
    }



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
                    bool is_nan = amrex::isnan(inter(i,j,k));
                    if (idim == 0) {
                        if (lst(i,j,k) == Real(0.0) ||
                            (lst(i,j,k) > Real(0.0) && is_nan))
                        {
                            // interp might still be quiet_nan because lst that
                            // was set to zero has been changed by FillBoundary
                            // at periodic bounadries.
                            inter(i,j,k) = problo[0] + i*dx[0];
                        }
                        else if (lst(i+1,j,k) == Real(0.0) ||
                                 (lst(i+1,j,k) > Real(0.0) && is_nan))
                        {
                            inter(i,j,k) = problo[0] + (i+1)*dx[0];
                        }
                    } else if (idim == 1) {
                        if (lst(i,j,k) == Real(0.0) ||
                            (lst(i,j,k) > Real(0.0) && is_nan))
                        {
                            inter(i,j,k) = problo[1] + j*dx[1];
                        }
                        else if (lst(i,j+1,k) == Real(0.0) ||
                                 (lst(i,j+1,k) > Real(0.0) && is_nan))
                        {
                            inter(i,j,k) = problo[1] + (j+1)*dx[1];
                        }
                    } else {
                        if (lst(i,j,k) == Real(0.0) ||
                            (lst(i,j,k) > Real(0.0) && is_nan))
                        {
                            inter(i,j,k) = problo[2] + k*dx[2];
                        }
                        else if (lst(i,j,k+1) == Real(0.0) ||
                                 (lst(i,j,k+1) > Real(0.0) && is_nan))
                        {
                            inter(i,j,k) = problo[2] + (k+1)*dx[2];
                        }
                    }
                }
            });
        }
    }

private:

    F m_f;
    R m_resource;  // We use this to hold the ownership of resource for F if needed,
                   // because F needs to be a simply type suitable for GPU.
};

template <class F>
GeometryShop<typename std::decay<F>::type>
makeShop (F&& f)
{
    return GeometryShop<typename std::decay<F>::type>(std::forward<F>(f));
}

template <class F, class R>
GeometryShop<typename std::decay<F>::type, typename std::decay<R>::type>
makeShop (F&& f, R&& r)
{
    return GeometryShop<typename std::decay<F>::type,
                        typename std::decay<R>::type>
        (std::forward<F>(f), std::forward<R>(r));
}

}

#endif
