#ifndef AMREX_PhysBCFunct_H_
#define AMREX_PhysBCFunct_H_

#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_FilCC_C.H>
#include <AMReX_TypeTraits.H>

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
}

typedef void (*UserFillBox) (Box const& bx, FArrayBox& dest,
                             const int dcomp, const int numcomp,
                             GeometryData const& geom, const Real time,
                             const BCRec* bcr, const int bcomp,
                             const int orig_comp);

//! This version calls function working on array
class BndryFuncArray
{
public:
    BndryFuncArray () noexcept {}
    BndryFuncArray (BndryFuncDefault inFunc) noexcept : m_func(inFunc) {}
    BndryFuncArray (BndryFunc3DDefault inFunc) noexcept : m_func3D(inFunc) {}

    void operator() (Box const& bx, FArrayBox& dest,
                     const int dcomp, const int numcomp,
                     Geometry const& geom, const Real time,
                     const Vector<BCRec>& bcr, const int bcomp,
                     const int orig_comp);

    bool RunOnGPU () const noexcept { return m_run_on_gpu; }
    void setRunOnGPU (bool b) noexcept { m_run_on_gpu = b; }

protected:
    BndryFuncDefault   m_func   = nullptr;
    BndryFunc3DDefault m_func3D = nullptr;
    bool m_run_on_gpu = false;
};

/**
* In this gpu version, F is provided by the user.  It needs to have a
* __device__ operator() that can work on a cell for boundaries not
* handled by amrex_fab_filcc.
*/
template <class F>
class GpuBndryFuncFab
{
public:
    GpuBndryFuncFab () = default;
    GpuBndryFuncFab (F a_f) : m_user_f(a_f) {}

    void operator() (Box const& bx, FArrayBox& dest,
                     const int dcomp, const int numcomp,
                     Geometry const& geom, const Real time,
                     const Vector<BCRec>& bcr, const int bcomp,
                     const int orig_comp);
protected:
    F m_user_f;
};

//! This cpu version calls function working on FArrayBox
class CpuBndryFuncFab
{
public:
    CpuBndryFuncFab () = default;
    CpuBndryFuncFab (UserFillBox a_f) : f_user(a_f) {}

    void operator() (Box const& bx, FArrayBox& dest,
                     const int dcomp, const int numcomp,
                     Geometry const& geom, const Real time,
                     const Vector<BCRec>& bcr, const int bcomp,
                     const int orig_comp);

protected:
    UserFillBox f_user = nullptr;
};

class PhysBCFunctBase {
public:
    //! note this is pure virtual and must be implemented in a derived class
    //! (see PhysBCFunct below for an example)
    //!
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

        //! create a grown domain box containing valid + periodic cells
        Box gdomain = amrex::convert(domain, mf.boxArray().ixType());
        const IntVect& ngrow = mf.nGrowVect();
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            if (m_geom.isPeriodic(i)) {
                gdomain.grow(i, ngrow[i]);
            }
        }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            Vector<BCRec> bcrs(ncomp);
            for (MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                FArrayBox* dest = mf.fabPtr(mfi);
                const Box& bx = mfi.fabbox();

                //! if there are cells not in the valid + periodic grown box
                //! we need to fill them here
                //!
                if (!gdomain.contains(bx))
                {
                    //! Based on BCRec for the domain, we need to make BCRec for this Box
                    amrex::setBC(bx, domain, bccomp, 0, ncomp, m_bcr, bcrs);

                    //! Note that we pass 0 as starting component of bcrs.
                    m_f(bx, *dest, icomp, ncomp, m_geom, time, bcrs, 0, bccomp);
                }
            }
        }
    }

private:
    Geometry      m_geom;
    Vector<BCRec> m_bcr;
    F             m_f;
};

template <class F>
void
GpuBndryFuncFab<F>::operator() (Box const& bx, FArrayBox& dest,
                                const int dcomp, const int numcomp,
                                Geometry const& geom, const Real time,
                                const Vector<BCRec>& bcr, const int bcomp,
                                const int orig_comp)
{
    FArrayBox* fab = &dest;
    const auto geomdata = geom.data();

    AsyncArray<BCRec> bcr_aa(bcr.data()+bcomp, numcomp);
    BCRec* bcr_p = bcr_aa.data();

    const IntVect& len = bx.length();

    Box gdomain = geom.Domain();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (Geometry::isPeriodic(idim)) {
            gdomain.grow(idim,len[idim]);
        }
    }

    if (gdomain.contains(bx)) return;

    const auto f_user = m_user_f;

    // filcc on the faces first
    {
        Array<Box,2*AMREX_SPACEDIM> dom_face_boxes
            = { AMREX_D_DECL(amrex::adjCellLo(gdomain, 0, len[0]),
                             amrex::adjCellLo(gdomain, 1, len[1]),
                             amrex::adjCellLo(gdomain, 2, len[2])),
                AMREX_D_DECL(amrex::adjCellHi(gdomain, 0, len[0]),
                             amrex::adjCellHi(gdomain, 1, len[1]),
                             amrex::adjCellHi(gdomain, 2, len[2])) };

        Vector<Box> face_boxes;
        for (const Box& b : dom_face_boxes) {
            Box tmp = b & bx;
            if (tmp.ok()) face_boxes.push_back(tmp);
        }
        const int n_face_boxes = face_boxes.size();
        if (n_face_boxes > 0) {
            AsyncArray<Box> face_boxes_aa(face_boxes.data(), n_face_boxes);
            Box* boxes_p = face_boxes_aa.data();

            long ncells = 0;
            for (const auto& b : face_boxes) {
                ncells += b.numPts();
            }

            AMREX_LAUNCH_DEVICE_LAMBDA (ncells, icell,
            {
                const auto& cell = getCell(boxes_p, n_face_boxes, icell);
                filcc_cell(cell, *fab, dcomp, numcomp, geomdata, time,
                           bcr_p, 0, orig_comp);
                f_user(cell, *fab, dcomp, numcomp, geomdata, time,
                       bcr_p, 0, orig_comp);
            });
        }
    }

#if (AMREX_SPACEDIM >= 2)
    // filcc on the edges
    {
#if (AMREX_SPACEDIM == 2)
        Array<Box,4> dom_edge_boxes
            = { amrex::adjCellLo(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),
                amrex::adjCellLo(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),
                amrex::adjCellHi(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),
                amrex::adjCellHi(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]) };
#else
        Array<Box,12> dom_edge_boxes
            = { amrex::adjCellLo(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),
                amrex::adjCellLo(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),
                amrex::adjCellHi(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),
                amrex::adjCellHi(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),
                //
                amrex::adjCellLo(amrex::adjCellLo(gdomain,0,len[0]),2,len[2]),
                amrex::adjCellLo(amrex::adjCellHi(gdomain,0,len[0]),2,len[2]),
                amrex::adjCellHi(amrex::adjCellLo(gdomain,0,len[0]),2,len[2]),
                amrex::adjCellHi(amrex::adjCellHi(gdomain,0,len[0]),2,len[2]),
                //
                amrex::adjCellLo(amrex::adjCellLo(gdomain,1,len[1]),2,len[2]),
                amrex::adjCellLo(amrex::adjCellHi(gdomain,1,len[1]),2,len[2]),
                amrex::adjCellHi(amrex::adjCellLo(gdomain,1,len[1]),2,len[2]),
                amrex::adjCellHi(amrex::adjCellHi(gdomain,1,len[1]),2,len[2]) };
#endif

        Vector<Box> edge_boxes;
        for (const Box& b : dom_edge_boxes) {
            Box tmp = b & bx;
            if (tmp.ok()) edge_boxes.push_back(tmp);
        }
        const int n_edge_boxes = edge_boxes.size();
        if (n_edge_boxes > 0) {
            AsyncArray<Box> edge_boxes_aa(edge_boxes.data(), n_edge_boxes);
            Box* boxes_p = edge_boxes_aa.data();

            long ncells = 0;
            for (const auto& b : edge_boxes) {
                ncells += b.numPts();
            }

            AMREX_LAUNCH_DEVICE_LAMBDA (ncells, icell,
            {
                const auto& cell = getCell(boxes_p, n_edge_boxes, icell);
                filcc_cell(cell, *fab, dcomp, numcomp, geomdata, time,
                           bcr_p, 0, orig_comp);
                f_user(cell, *fab, dcomp, numcomp, geomdata, time,
                       bcr_p, 0, orig_comp);
            });
        }
    }
#endif

#if (AMREX_SPACEDIM == 3)
    // filcc on corners
    {
        Array<Box,8> dom_corner_boxes
            = { amrex::adjCellLo(amrex::adjCellLo(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),2,len[2]),
                amrex::adjCellLo(amrex::adjCellLo(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),2,len[2]),
                amrex::adjCellLo(amrex::adjCellHi(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),2,len[2]),
                amrex::adjCellLo(amrex::adjCellHi(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),2,len[2]),
                amrex::adjCellHi(amrex::adjCellLo(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),2,len[2]),
                amrex::adjCellHi(amrex::adjCellLo(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),2,len[2]),
                amrex::adjCellHi(amrex::adjCellHi(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),2,len[2]),
                amrex::adjCellHi(amrex::adjCellHi(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),2,len[2]) };

        Vector<Box> corner_boxes;
        for (const Box& b : dom_corner_boxes) {
            Box tmp = b & bx;
            if (tmp.ok()) corner_boxes.push_back(tmp);
        }
        const int n_corner_boxes = corner_boxes.size();
        if (n_corner_boxes > 0) {
            AsyncArray<Box> corner_boxes_aa(corner_boxes.data(), n_corner_boxes);
            Box* boxes_p = corner_boxes_aa.data();

            long ncells = 0;
            for (const auto& b : corner_boxes) {
                ncells += b.numPts();
            }

            AMREX_LAUNCH_DEVICE_LAMBDA (ncells, icell,
            {
                const auto& cell = getCell(boxes_p, n_corner_boxes, icell);
                filcc_cell(cell, *fab, dcomp, numcomp, geomdata, time,
                           bcr_p, 0, orig_comp);
                f_user(cell, *fab, dcomp, numcomp, geomdata, time,
                       bcr_p, 0, orig_comp);
            });
        }
    }
#endif
}

}

#endif
