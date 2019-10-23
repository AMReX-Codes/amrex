#ifndef AMREX_EB2_GEOMETRYSHOP_H_
#define AMREX_EB2_GEOMETRYSHOP_H_

#include <AMReX_EB2_IF_Base.H>
#include <AMReX_EB2_Graph.H>
#include <AMReX_Geometry.H>
#include <AMReX_BaseFab.H>
#include <AMReX_Print.H>
#include <AMReX_Array.H>
#include <memory>
#include <type_traits>
#include <cmath>

namespace amrex { namespace EB2 {

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
#ifdef AMREX_DEVICE_COMPILE
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
    const Real tol = 1.e-12;
    const int MAXITER = 100;
    const Real EPS = 3.0e-15;

    Real p, q, r, s;

    GpuArray<Real,AMREX_SPACEDIM> aPt = lo;
    GpuArray<Real,AMREX_SPACEDIM> bPt = hi;

    Real fa = IF_f(f, aPt);
    Real fb = IF_f(f, bPt);
    Real c = bPt[rangedir];
    Real fc = fb;

    if (fb*fa > 0.0) {
//        amrex::AllPrint() << "fa " << fa << " fb " << fb << "\n";
        amrex::Error("BrentRootFinder. Root must be bracketed, but instead the supplied end points have the same sign.");
        return 0.0;
    } else if (fa == 0.0) {
        return aPt[rangedir];
    } else if (fb == 0.0) {
        return bPt[rangedir];
    }

    Real d = 0.0, e = 0.0;
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
        Real tol1  = 2.0 * EPS * std::abs(bPt[rangedir]) + 0.5 * tol;
        Real xm    = 0.5 * (c - bPt[rangedir]);

        if (std::abs(xm) <= tol1 || fb == 0.0)
        {
            break;
        }

        if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb))
        {
            //  Attempt inverse quadratic interpolation
            s = fb / fa;
            if (aPt[rangedir] == c)
            {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            }
            else
            {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * xm * q * (q-r) - (bPt[rangedir]-aPt[rangedir]) * (r-1.0));
                q = (q-1.0) * (r-1.0) * (s-1.0);
            }

            //  Check whether in bounds
            if (p > 0) q = -q;

            p = std::abs(p);

            if (2.0 * p < amrex::min(3.0*xm*q-std::abs(tol1*q), 1.0*std::abs(e*q)))
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

template <class F>
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

    explicit GeometryShop (F const& f)
        : m_f(f)
        {}

    explicit GeometryShop (F && f)
        : m_f(std::move(f))
        {}

    ~GeometryShop() {}

    GeometryShop (GeometryShop<F> const& rhs)
        : m_f(rhs.m_f)
        {}

    GeometryShop (GeometryShop<F> && rhs)
        : m_f(std::move(rhs.m_f))
        {}

    GeometryShop<F>& operator= (GeometryShop<F> const& rhs) = delete;
    GeometryShop<F>& operator= (GeometryShop<F> && rhs) = delete;

    F const& GetImpFunc () const& { return m_f; }
    F&& GetImpFunc () && { return std::move(m_f); }

    int getBoxType_Cpu (const Box& bx, Geometry const& geom) const noexcept
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
                    if (v == 0.0) {
                        ++nzero;
                    } else if (v > 0.0) {
                        ++nbody;
                    } else {
                        ++nfluid;
                    }
                    if (nbody > 0 && nfluid > 0) return mixedcells;
                }
            }
        }

        if (nbody == 0) {
            return allregular;
        } else if (nfluid == 0) {
            return allcovered;
        } else {
            return mixedcells;
        }
    }

    template <class U=F, typename std::enable_if<IsGPUable<U>::value>::type* FOO = nullptr >
    int getBoxType (const Box& bx, const Geometry& geom, RunOn run_on) const noexcept
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
                Real v = f(AMREX_D_DECL(problo[0]+i*dx[0],
                                        problo[1]+j*dx[1],
                                        problo[2]+k*dx[2]));
                int nbody = 0, nzero = 0, nfluid = 0;
                if (v == 0.0) {
                    ++nzero;
                } else if (v > 0.0) {
                    ++nbody;
                } else {
                    ++nfluid;
                }
                return {nbody,nzero,nfluid};
            });
            ReduceTuple hv = reduce_data.value();
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
    int getBoxType (const Box& bx, const Geometry& geom, RunOn) const noexcept
    {
        return getBoxType_Cpu(bx, geom);
    }

    template <class U=F, typename std::enable_if<IsGPUable<U>::value>::type* FOO = nullptr >
    static constexpr bool isGPUable () noexcept { return true; }

    template <class U=F, typename std::enable_if<!IsGPUable<U>::value>::type* BAR = nullptr >
    static constexpr bool isGPUable () noexcept { return false; }

    template <class U=F, typename std::enable_if<IsGPUable<U>::value>::type* FOO = nullptr >
    void fillFab (BaseFab<Real>& levelset, const Geometry& geom, RunOn run_on) const noexcept
    {
        const auto problo = geom.ProbLoArray();
        const auto dx = geom.CellSizeArray();
        const Box& bx = levelset.box();
        const auto& a = levelset.array();
        auto f = m_f;
        bool run_on_gpu = (run_on == RunOn::Gpu && Gpu::inLaunchRegion());
        AMREX_HOST_DEVICE_FOR_3D_FLAG(run_on_gpu, bx, i, j, k,
        {
            a(i,j,k) = f(AMREX_D_DECL(problo[0]+i*dx[0],
                                      problo[1]+j*dx[1],
                                      problo[2]+k*dx[2]));
        });
    }

    template <class U=F, typename std::enable_if<!IsGPUable<U>::value>::type* BAR = nullptr >
    void fillFab (BaseFab<Real>& levelset, const Geometry& geom, RunOn) const noexcept
    {
        const auto problo = geom.ProbLoArray();
        const auto dx = geom.CellSizeArray();
        const Box& bx = levelset.box();
        const auto& a = levelset.array();
        amrex::Loop(bx, [=] (int i, int j, int k) noexcept
        {
            a(i,j,k) = m_f(RealArray{AMREX_D_DECL(problo[0]+i*dx[0],
                                                  problo[1]+j*dx[1],
                                                  problo[2]+k*dx[2])});
        });
    }

    template <class U=F, typename std::enable_if<IsGPUable<U>::value>::type* FOO = nullptr >
    void getIntercept (Array<BaseFab<Real>,AMREX_SPACEDIM>& inter_fab,
                       Array<BaseFab<Type_t>,AMREX_SPACEDIM> const& type_fab,
                       Geometry const& geom, RunOn run_on) const noexcept
    {
        auto const& dx = geom.CellSizeArray();
        auto const& problo = geom.ProbLoArray();
        auto f = m_f;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            Array4<Real> const& inter = inter_fab[idim].array();
            Array4<Type_t const> const& type = type_fab[idim].array();
            const Box& bx = inter_fab[idim].box();
            int dir = (AMREX_SPACEDIM==2) ? 1-idim : idim;
            bool run_on_gpu = (run_on == RunOn::Gpu && Gpu::inLaunchRegion());
            AMREX_HOST_DEVICE_FOR_3D_FLAG(run_on_gpu, bx, i, j, k,
            {
                if (type(i,j,k) == Type::irregular) {
                    IntVect ivlo(AMREX_D_DECL(i,j,k));
                    IntVect ivhi(AMREX_D_DECL(i,j,k));
                    ivhi[dir] += 1;
                    inter(i,j,k) = BrentRootFinder
                        ({AMREX_D_DECL(problo[0]+ivlo[0]*dx[0],
                                       problo[1]+ivlo[1]*dx[1],
                                       problo[2]+ivlo[2]*dx[2])},
                         {AMREX_D_DECL(problo[0]+ivhi[0]*dx[0],
                                       problo[1]+ivhi[1]*dx[1],
                                       problo[2]+ivhi[2]*dx[2])},
                            dir, f);
                } else {
                    inter(i,j,k) = std::numeric_limits<Real>::quiet_NaN();
                }
            });
        }
    }

    template <class U=F, typename std::enable_if<!IsGPUable<U>::value>::type* BAR = nullptr >
    void getIntercept (Array<BaseFab<Real>,AMREX_SPACEDIM>& inter_fab,
                       Array<BaseFab<Type_t>,AMREX_SPACEDIM> const& type_fab,
                       Geometry const& geom, RunOn) const noexcept
    {
        auto const& dx = geom.CellSizeArray();
        auto const& problo = geom.ProbLoArray();
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            Array4<Real> const& inter = inter_fab[idim].array();
            Array4<Type_t const> const& type = type_fab[idim].array();
            const Box& bx = inter_fab[idim].box();
            int dir = (AMREX_SPACEDIM==2) ? 1-idim : idim;
            amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
            {
                if (type(i,j,k) == Type::irregular) {
                    IntVect ivlo(AMREX_D_DECL(i,j,k));
                    IntVect ivhi(AMREX_D_DECL(i,j,k));
                    ivhi[dir] += 1;
                    inter(i,j,k) = BrentRootFinder
                        ({AMREX_D_DECL(problo[0]+ivlo[0]*dx[0],
                                       problo[1]+ivlo[1]*dx[1],
                                       problo[2]+ivlo[2]*dx[2])},
                         {AMREX_D_DECL(problo[0]+ivhi[0]*dx[0],
                                       problo[1]+ivhi[1]*dx[1],
                                       problo[2]+ivhi[2]*dx[2])},
                          dir, m_f);
                } else {
                    inter(i,j,k) = std::numeric_limits<Real>::quiet_NaN();
                }
            });
        }
    }

private:

    F m_f;

};

template <class F>
GeometryShop<typename std::decay<F>::type>
makeShop (F&& f)
{
    return GeometryShop<typename std::decay<F>::type>(std::forward<F>(f));
}

}}

#endif
