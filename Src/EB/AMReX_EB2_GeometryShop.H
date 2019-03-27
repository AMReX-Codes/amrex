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

    GeometryShop (F const& f)
        : m_f(f)
        {}

    GeometryShop (F && f)
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

    int getBoxType (const Box& bx, const Geometry& geom) const;

    template <class U=F, typename std::enable_if<IsGPUable<U>::value>::type* FOO = nullptr >
    void fillFab (BaseFab<Real>& levelset, const Geometry& geom) const noexcept
    {
        const auto problo = geom.ProbLoArray();
        const auto dx = geom.CellSizeArray();
        const Box& bx = levelset.box();
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);
        const auto& a = levelset.array();
        
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    a(i,j,k) = m_f(AMREX_D_DECL(problo[0]+i*dx[0],
                                                problo[1]+j*dx[1],
                                                problo[2]+k*dx[2]));
                }
            }
        }
    }

    template <class U=F, typename std::enable_if<!IsGPUable<U>::value>::type* = nullptr >
    void fillFab (BaseFab<Real>& levelset, const Geometry& geom) const noexcept
    {
        const Real* AMREX_RESTRICT problo = geom.ProbLo();
        const Real* AMREX_RESTRICT dx = geom.CellSize();
        
        const Box& bx = levelset.box();
        const auto len = amrex::length(bx);
        const auto lo  = amrex::lbound(bx);
        const auto dp  = levelset.view(lo);
        
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    RealArray xyz {AMREX_D_DECL(problo[0]+(i+lo.x)*dx[0],
                                                problo[1]+(j+lo.y)*dx[1],
                                                problo[2]+(k+lo.z)*dx[2])};
                    dp(i,j,k,0) = m_f(xyz);
                }
            }
        }
    }

    Real getIntercept (int edgedir, const IntVect& iv, const Geometry& geom) const;

    Real BrentRootFinder (const RealArray& lo, const RealArray& hi, int rangedir) const;

private:

    F m_f;

};

template <class F>
GeometryShop<typename std::decay<F>::type>
makeShop (F&& f)
{
    return GeometryShop<typename std::decay<F>::type>(std::forward<F>(f));
}

template <class F>
int
GeometryShop<F>::getBoxType (const Box& bx, const Geometry& geom) const
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

template <class F>
Real
GeometryShop<F>::getIntercept (int edgedir, const IntVect& iv, const Geometry& geom) const
{
    const Real* problo = geom.ProbLo();
    const Real* dx = geom.CellSize();

#if (AMREX_SPACEDIM == 1)

    amrex::Abort("GeometryShop: 1d not supported");
    return 0;

#elif (AMREX_SPACEDIM == 2)

    IntVect ivlo = iv;
    IntVect ivhi = iv;  ivhi[edgedir] += 1;
    Real r = this->BrentRootFinder(RealArray{problo[0]+ivlo[0]*dx[0],
                                             problo[1]+ivlo[1]*dx[1]},
                                   RealArray{problo[0]+ivhi[0]*dx[0],
                                             problo[1]+ivhi[1]*dx[1]},
                                   edgedir);
    return r;

#else

    IntVect ivlo = iv;
    IntVect ivhi = iv;  ivhi[edgedir] += 1;
    Real r = this->BrentRootFinder(RealArray{problo[0]+ivlo[0]*dx[0],
                                             problo[1]+ivlo[1]*dx[1],
                                             problo[2]+ivlo[2]*dx[2]},
                                   RealArray{problo[0]+ivhi[0]*dx[0],
                                             problo[1]+ivhi[1]*dx[1],
                                             problo[2]+ivhi[2]*dx[2]},
                                   edgedir);
    return r;

#endif
}

template <class F>
Real
GeometryShop<F>::BrentRootFinder (const RealArray& lo, const RealArray& hi, int rangedir) const
{
    const Real tol = 1.e-12;
    const int MAXITER = 100;
    const Real EPS = 3.0e-15;

    Real p, q, r, s;

    RealArray aPt = lo;
    RealArray bPt = hi;

    Real fa = m_f(aPt);
    Real fb = m_f(bPt);
    Real c = bPt[rangedir];
    Real fc = fb;

    if (fb*fa > 0.0) {
//        amrex::AllPrint() << "fa " << fa << " fb " << fb << "\n";
        amrex::Error("GeometryShop::BrentRootFinder. Root must be bracketed, but instead the supplied end points have the same sign.");
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
            
            if (2.0 * p < std::min(3.0*xm*q-std::abs(tol1*q), 1.0*std::abs(e*q)))
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
        
        fb = m_f(bPt);
    }

    if (i >= MAXITER)
    {
        amrex::Error("BrentRootFinder: exceeding maximum iterations.");
    }

    return bPt[rangedir];
}

}}

#endif
