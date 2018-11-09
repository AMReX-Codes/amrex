#ifndef AMREX_PhysBCFunct_H_
#define AMREX_PhysBCFunct_H_

#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

namespace amrex {

extern "C"
{
    typedef void (*BndryFuncDefault)(Real* data, AMREX_ARLIM_P(lo), AMREX_ARLIM_P(hi),
                                     const int* dom_lo, const int* dom_hi,
                                     const Real* dx, const Real* grd_lo,
                                     const Real* time, const int* bc);
    typedef void (*BndryFunc3DDefault)(Real* data, const int* lo, const int* hi,
                                       const int* dom_lo, const int* dom_hi,
                                       const Real* dx, const Real* grd_lo,
                                       const Real* time, const int* bc);

    // Note that i, j, k, nc, and time are passed by value.
    // Also note that prob_lo is the lo of the whole domain.
    typedef void (*UserFillCell) (const int i, const int j, const int k,
                                  Real* data, const int* dlo, const int* dhi, const int nc,
                                  const int* dom_lo, const int* dom_hi,
                                  const Real* dx, const Real* prob_lo,
                                  const Real time, const int* bc,
                                  const int orig_comp);

    // Note that nc and time are passed by value.
    // Also note that prob_lo is the lo of the whole domain.
    typedef void (*UserFillBox) (Real* data, const int* dlo, const int* dhi, const int nc,
                                 const int* dom_lo, const int* dom_hi,
                                 const Real* dx, const Real* prob_lo,
                                 const Real time, const int* bc,
                                 const int orig_comp);
}

// This version calls function working on array
class BndryFuncArray
{
public:
    BndryFuncArray () {}
    BndryFuncArray (BndryFuncDefault inFunc) : m_func(inFunc) {}
    BndryFuncArray (BndryFunc3DDefault inFunc) : m_func3D(inFunc) {}

    void operator() (Box const& bx, FArrayBox& data,
                     const int dcomp, const int numcomp,
                     Geometry const& geom, const Real time,
                     const Vector<BCRec>& bcr, const int bcomp,
                     const int orig_comp);

    bool RunOnGPU () const { return m_run_on_gpu; }
    void setRunOnGPU (bool b) { m_run_on_gpu = b; }

protected:
    BndryFuncDefault   m_func   = nullptr;
    BndryFunc3DDefault m_func3D = nullptr;
    bool m_run_on_gpu = false;
};

// This gpu version calls function working on FArrayBox
class GpuBndryFuncFab
{
public:
    GpuBndryFuncFab () = default;
    GpuBndryFuncFab (UserFillCell a_f) : f_user(a_f) {}

    void operator() (Box const& bx, FArrayBox& data,
                     const int dcomp, const int numcomp,
                     Geometry const& geom, const Real time,
                     const Vector<BCRec>& bcr, const int bcomp,
                     const int orig_comp);

protected:
    UserFillCell f_user = nullptr;
};

// This cpu version calls function working on FArrayBox
class CpuBndryFuncFab
{
public:
    CpuBndryFuncFab () = default;
    CpuBndryFuncFab (UserFillBox a_f) : f_user(a_f) {}

    void operator() (Box const& bx, FArrayBox& data,
                     const int dcomp, const int numcomp,
                     Geometry const& geom, const Real time,
                     const Vector<BCRec>& bcr, const int bcomp,
                     const int orig_comp);

protected:
    UserFillBox f_user = nullptr;
};

class PhysBCFunctBase {
public:
    // note this is pure virtual and must be implemented in a derived class
    // (see PhysBCFunct below for an example)
    virtual void FillBoundary (MultiFab& mf, int dcomp, int ncomp, Real time, int bccomp) = 0;
    virtual ~PhysBCFunctBase () = 0;
};

inline PhysBCFunctBase::~PhysBCFunctBase () {}

template <class F>
class PhysBCFunct
    : public PhysBCFunctBase
{
public:
    PhysBCFunct () {}

    PhysBCFunct (const Geometry& geom, const Vector<BCRec>& bcr, F const& f)
        : m_geom(geom), m_bcr(bcr), m_f(f)
        {}

    virtual ~PhysBCFunct () {}

    void define (const Geometry& geom, const Vector<BCRec>& bcr, F const& f) {
        m_geom = geom; m_bcr = bcr; m_f = f;
    }

    virtual void FillBoundary (MultiFab& mf, int icomp, int ncomp, Real time, int bccomp) override
    {
        if (m_geom.isAllPeriodic()) return;

        BL_PROFILE("PhysBCFunct::FillBoundary");

        const Box&     domain      = m_geom.Domain();
        const int*     dlo         = domain.loVect();
        const int*     dhi         = domain.hiVect();
        const Real*    dx          = m_geom.CellSize();
        const RealBox& prob_domain = m_geom.ProbDomain();
        const Real*    problo      = prob_domain.lo();
        
        // create a grown domain box containing valid + periodic cells
        Box gdomain = amrex::convert(domain, mf.boxArray().ixType());
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            if (m_geom.isPeriodic(i)) {
                gdomain.grow(i, mf.nGrow());
            }
        }
        
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            Vector<BCRec> bcrs(ncomp);
            for (MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                FArrayBox& dest = mf[mfi];
                const Box& bx = mfi.fabbox();
                
                // if there are cells not in the valid + periodic grown box
                // we need to fill them here
                if (!gdomain.contains(bx))
                {
                    // Based on BCRec for the domain, we need to make BCRec for this Box
                    amrex::setBC(bx, domain, bccomp, 0, ncomp, m_bcr, bcrs);

                    // Note that we pass 0 as starting component of bcrs.
                    m_f(bx, dest, icomp, ncomp, m_geom, time, m_bcr, 0, bccomp);
                }
            }
        }
    }

private:
    Geometry      m_geom;
    Vector<BCRec> m_bcr;
    F             m_f;
};

}

#endif
