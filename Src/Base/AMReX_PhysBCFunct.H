#ifndef AMREX_PhysBCFunct_H_
#define AMREX_PhysBCFunct_H_

#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_FilCC_C.H>
#include <AMReX_FilND_C.H>
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

typedef void (*UserFillBox) (Box const& bx, Array4<Real> const& dest,
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
//    void setRunOnGPU (bool b) noexcept { m_run_on_gpu = b; }

protected:
    BndryFuncDefault   m_func   = nullptr;
    BndryFunc3DDefault m_func3D = nullptr;
    bool m_run_on_gpu = false;
};

/**
* In this gpu version, F is provided by the user.  It needs to have a
* __device__ operator() that can work on a cell/node for boundaries not
* handled by amrex::fab_filcc/fab_filnd.
*/
template <class F>
class GpuBndryFuncFab
{
public:
    GpuBndryFuncFab () = default;
    GpuBndryFuncFab (F const& a_f) : m_user_f(a_f) {}
    GpuBndryFuncFab (F&& a_f) : m_user_f(std::move(a_f)) {}

    void operator() (Box const& bx, FArrayBox& dest,
                     const int dcomp, const int numcomp,
                     Geometry const& geom, const Real time,
                     const Vector<BCRec>& bcr, const int bcomp,
                     const int orig_comp);

    void ccdoit (Box const& bx, FArrayBox& dest,
                 const int dcomp, const int numcomp,
                 Geometry const& geom, const Real time,
                 const Vector<BCRec>& bcr, const int bcomp,
                 const int orig_comp);

    void nddoit (Box const& bx, FArrayBox& dest,
                 const int dcomp, const int numcomp,
                 Geometry const& geom, const Real time,
                 const Vector<BCRec>& bcr, const int bcomp,
                 const int orig_comp);
protected:
    F m_user_f;
};

#if !(defined(AMREX_USE_CUDA) && defined(AMREX_USE_GPU_PRAGMA) && defined(AMREX_GPU_PRAGMA_NO_HOST))
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
#endif

class PhysBCFunctNoOp
{
public:
    void operator() (MultiFab& /*mf*/, int /*dcomp*/, int /*ncomp*/, IntVect const& /*nghost*/,
                     Real /*time*/, int /*bccomp*/) {}
};


template <class F>
class PhysBCFunct
{
public:
    PhysBCFunct () {}

    PhysBCFunct (const Geometry& geom, const Vector<BCRec>& bcr, F const& f)
        : m_geom(geom), m_bcr(bcr), m_f(f)
        {}

    PhysBCFunct (const Geometry& geom, const Vector<BCRec>& bcr, F&& f)
        : m_geom(geom), m_bcr(bcr), m_f(std::move(f))
        {}

    ~PhysBCFunct () {}

    void define (const Geometry& geom, const Vector<BCRec>& bcr, F const& f) {
        m_geom = geom; m_bcr = bcr; m_f = f;
    }

    void define (const Geometry& geom, const Vector<BCRec>& bcr, F&& f) {
        m_geom = geom; m_bcr = bcr; m_f = std::move(f);
    }

    void operator() (MultiFab& mf, int icomp, int ncomp, IntVect const& nghost,
                     Real time, int bccomp)
    {
        if (m_geom.isAllPeriodic()) return;

        BL_PROFILE("PhysBCFunct::()");

        //! create a grown domain box containing valid + periodic cells
        const Box& domain = m_geom.Domain();
        Box gdomain = amrex::convert(domain, mf.boxArray().ixType());
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            if (m_geom.isPeriodic(i)) {
                gdomain.grow(i, nghost[i]);
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

                //! if there are cells not in the valid + periodic grown box
                //! we need to fill them here
                //!
                if (!gdomain.contains(bx))
                {
                    //! Based on BCRec for the domain, we need to make BCRec for this Box
                    amrex::setBC(bx, domain, bccomp, 0, ncomp, m_bcr, bcrs);

                    //! Note that we pass 0 as starting component of bcrs.
                    m_f(bx, dest, icomp, ncomp, m_geom, time, bcrs, 0, bccomp);
                }
            }
        }
    }

    // For backward compatibility
    void FillBoundary (MultiFab& mf, int dcomp, int ncomp, IntVect const& nghost,
                       Real time, int bccomp) {
        this->operator()(mf,dcomp,ncomp,nghost,time,bccomp);
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
    if (bx.ixType().cellCentered()) {
        ccdoit(bx,dest,dcomp,numcomp,geom,time,bcr,bcomp,orig_comp);
    } else if (bx.ixType().nodeCentered()) {
        nddoit(bx,dest,dcomp,numcomp,geom,time,bcr,bcomp,orig_comp);
    } else {
        amrex::Abort("GpuBndryFuncFab: mixed type not supported");
    }
}

template <class F>
void
GpuBndryFuncFab<F>::ccdoit (Box const& bx, FArrayBox& dest,
                            const int dcomp, const int numcomp,
                            Geometry const& geom, const Real time,
                            const Vector<BCRec>& bcr, const int bcomp,
                            const int orig_comp)
{
    const IntVect& len = bx.length();

    Box gdomain = geom.Domain();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (geom.isPeriodic(idim)) {
            gdomain.grow(idim,len[idim]);
        }
    }

    if (gdomain.contains(bx)) return;

    Array4<Real> const& fab = dest.array();
    const auto geomdata = geom.data();

#ifdef AMREX_USE_GPU

    AsyncArray<BCRec> bcr_aa(bcr.data()+bcomp, numcomp);
    BCRec* bcr_p = bcr_aa.data();

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
        if (n_face_boxes == 1) {
            amrex::ParallelFor(face_boxes[0],
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                AMREX_D_PICK(amrex::ignore_unused(j,k),amrex::ignore_unused(k),(void)0);
                IntVect const cell(AMREX_D_DECL(i,j,k));
                filcc_cell(cell, fab, dcomp, numcomp, geomdata, time,
                           bcr_p, 0, orig_comp);
                f_user(cell, fab, dcomp, numcomp, geomdata, time,
                       bcr_p, 0, orig_comp);
            });
        } else if (n_face_boxes > 1) {
            AsyncArray<Box> face_boxes_aa(face_boxes.data(), n_face_boxes);
            Box* boxes_p = face_boxes_aa.data();
            Long ncells = 0;
            for (const auto& b : face_boxes) {
                ncells += b.numPts();
            }
            amrex::ParallelFor(ncells,
            [=] AMREX_GPU_DEVICE (Long icell) noexcept
            {
                const auto& cell = getCell(boxes_p, n_face_boxes, icell);
                filcc_cell(cell, fab, dcomp, numcomp, geomdata, time,
                           bcr_p, 0, orig_comp);
                f_user(cell, fab, dcomp, numcomp, geomdata, time,
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
        if (n_edge_boxes == 1) {
            amrex::ParallelFor(edge_boxes[0],
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                AMREX_D_PICK(amrex::ignore_unused(j,k),amrex::ignore_unused(k),(void)0);
                IntVect const cell(AMREX_D_DECL(i,j,k));
                filcc_cell(cell, fab, dcomp, numcomp, geomdata, time,
                           bcr_p, 0, orig_comp);
                f_user(cell, fab, dcomp, numcomp, geomdata, time,
                       bcr_p, 0, orig_comp);
            });
        } else if (n_edge_boxes > 1) {
            AsyncArray<Box> edge_boxes_aa(edge_boxes.data(), n_edge_boxes);
            Box* boxes_p = edge_boxes_aa.data();
            Long ncells = 0;
            for (const auto& b : edge_boxes) {
                ncells += b.numPts();
            }
            amrex::ParallelFor(ncells,
            [=] AMREX_GPU_DEVICE (Long icell) noexcept
            {
                const auto& cell = getCell(boxes_p, n_edge_boxes, icell);
                filcc_cell(cell, fab, dcomp, numcomp, geomdata, time,
                           bcr_p, 0, orig_comp);
                f_user(cell, fab, dcomp, numcomp, geomdata, time,
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
        if (n_corner_boxes == 1) {
            amrex::ParallelFor(corner_boxes[0],
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                IntVect const cell(AMREX_D_DECL(i,j,k));
                filcc_cell(cell, fab, dcomp, numcomp, geomdata, time,
                           bcr_p, 0, orig_comp);
                f_user(cell, fab, dcomp, numcomp, geomdata, time,
                       bcr_p, 0, orig_comp);
            });
        } else if (n_corner_boxes > 1) {
            AsyncArray<Box> corner_boxes_aa(corner_boxes.data(), n_corner_boxes);
            Box* boxes_p = corner_boxes_aa.data();
            Long ncells = 0;
            for (const auto& b : corner_boxes) {
                ncells += b.numPts();
            }
            amrex::ParallelFor(ncells,
            [=] AMREX_GPU_DEVICE (Long icell) noexcept
            {
                const auto& cell = getCell(boxes_p, n_corner_boxes, icell);
                filcc_cell(cell, fab, dcomp, numcomp, geomdata, time,
                           bcr_p, 0, orig_comp);
                f_user(cell, fab, dcomp, numcomp, geomdata, time,
                       bcr_p, 0, orig_comp);
            });
        }
    }
#endif

#else

    BCRec const* bcr_p = bcr.data()+bcomp;

    const auto& f_user = m_user_f;

    // filcc on the faces first
    {
        Array<Box,2*AMREX_SPACEDIM> dom_face_boxes
            = {{ AMREX_D_DECL(amrex::adjCellLo(gdomain, 0, len[0]),
                              amrex::adjCellLo(gdomain, 1, len[1]),
                              amrex::adjCellLo(gdomain, 2, len[2])),
                 AMREX_D_DECL(amrex::adjCellHi(gdomain, 0, len[0]),
                              amrex::adjCellHi(gdomain, 1, len[1]),
                              amrex::adjCellHi(gdomain, 2, len[2])) }};

        for (const Box& b : dom_face_boxes) {
            Box tmp = b & bx;
            amrex::For(tmp, [=] (int i, int j, int k) noexcept
            {
                amrex::ignore_unused(j,k);
                IntVect const cell(AMREX_D_DECL(i,j,k));
                filcc_cell(cell, fab, dcomp, numcomp, geomdata, time,
                           bcr_p, 0, orig_comp);
                f_user(cell, fab, dcomp, numcomp, geomdata, time,
                       bcr_p, 0, orig_comp);
            });
        }
    }

#if (AMREX_SPACEDIM >= 2)
    // filcc on the edges
    {
#if (AMREX_SPACEDIM == 2)
        Array<Box,4> dom_edge_boxes
            = {{ amrex::adjCellLo(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),
                 amrex::adjCellLo(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),
                 amrex::adjCellHi(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),
                 amrex::adjCellHi(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]) }};
#else
        Array<Box,12> dom_edge_boxes
            = {{ amrex::adjCellLo(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),
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
                 amrex::adjCellHi(amrex::adjCellHi(gdomain,1,len[1]),2,len[2]) }};
#endif

        for (const Box& b : dom_edge_boxes) {
            Box tmp = b & bx;
            amrex::For(tmp, [=] (int i, int j, int k) noexcept
            {
                amrex::ignore_unused(j,k);
                IntVect const cell(AMREX_D_DECL(i,j,k));
                filcc_cell(cell, fab, dcomp, numcomp, geomdata, time,
                           bcr_p, 0, orig_comp);
                f_user(cell, fab, dcomp, numcomp, geomdata, time,
                       bcr_p, 0, orig_comp);
            });
        }
    }
#endif

#if (AMREX_SPACEDIM == 3)
    // filcc on corners
    {
        Array<Box,8> dom_corner_boxes
            = {{ amrex::adjCellLo(amrex::adjCellLo(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),2,len[2]),
                 amrex::adjCellLo(amrex::adjCellLo(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),2,len[2]),
                 amrex::adjCellLo(amrex::adjCellHi(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),2,len[2]),
                 amrex::adjCellLo(amrex::adjCellHi(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),2,len[2]),
                 amrex::adjCellHi(amrex::adjCellLo(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),2,len[2]),
                 amrex::adjCellHi(amrex::adjCellLo(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),2,len[2]),
                 amrex::adjCellHi(amrex::adjCellHi(amrex::adjCellLo(gdomain,0,len[0]),1,len[1]),2,len[2]),
                 amrex::adjCellHi(amrex::adjCellHi(amrex::adjCellHi(gdomain,0,len[0]),1,len[1]),2,len[2]) }};

        for (const Box& b : dom_corner_boxes) {
            Box tmp = b & bx;
            amrex::For(tmp, [=] (int i, int j, int k) noexcept
            {
                IntVect const cell(AMREX_D_DECL(i,j,k));
                filcc_cell(cell, fab, dcomp, numcomp, geomdata, time,
                           bcr_p, 0, orig_comp);
                f_user(cell, fab, dcomp, numcomp, geomdata, time,
                       bcr_p, 0, orig_comp);
            });
        }
    }
#endif

#endif
}

template <class F>
void
GpuBndryFuncFab<F>::nddoit (Box const& bx, FArrayBox& dest,
                            const int dcomp, const int numcomp,
                            Geometry const& geom, const Real time,
                            const Vector<BCRec>& bcr, const int bcomp,
                            const int orig_comp)
{
    const IntVect& len = bx.length();

    Box const& domain = amrex::convert(geom.Domain(),IntVect::TheNodeVector());;
    Box gdomain = domain;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (geom.isPeriodic(idim)) {
            gdomain.grow(idim,len[idim]);
        }
    }

    if (gdomain.contains(bx)) return;

    Array4<Real> const& fab = dest.array();
    const auto geomdata = geom.data();

    AsyncArray<BCRec> bcr_aa(bcr.data()+bcomp, numcomp);
    BCRec* bcr_p = bcr_aa.data();

    const auto f_user = m_user_f;

    // xlo
    if (!geom.isPeriodic(0) && bx.smallEnd(0) < domain.smallEnd(0)) {
        Box bndry = bx;
        int dxlo = domain.smallEnd(0);
        bndry.setBig(0,dxlo-1);
        amrex::ParallelFor(bndry, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for (int n = dcomp; n < dcomp+numcomp; ++n) {
                fab(i,j,k,n) = fab(dxlo,j,k,n);
            }
            f_user(IntVect(AMREX_D_DECL(i,j,k)), fab, dcomp, numcomp, geomdata, time,
                   bcr_p, 0, orig_comp);
        });
    }

    // xhi
    if (!geom.isPeriodic(0) && bx.bigEnd(0) > domain.bigEnd(0)) {
        Box bndry = bx;
        int dxhi = domain.bigEnd(0);
        bndry.setSmall(0,dxhi+1);
        amrex::ParallelFor(bndry, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for (int n = dcomp; n < dcomp+numcomp; ++n) {
                fab(i,j,k,n) = fab(dxhi,j,k,n);
            }
            f_user(IntVect(AMREX_D_DECL(i,j,k)), fab, dcomp, numcomp, geomdata, time,
                   bcr_p, 0, orig_comp);
        });
    }

#if (AMREX_SPACEDIM >= 2)
    // ylo
    if (!geom.isPeriodic(1) && bx.smallEnd(1) < domain.smallEnd(1)) {
        Box bndry = bx;
        int dylo = domain.smallEnd(1);
        bndry.setBig(1,dylo-1);
        amrex::ParallelFor(bndry, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for (int n = dcomp; n < dcomp+numcomp; ++n) {
                fab(i,j,k,n) = fab(i,dylo,k,n);
            }
            f_user(IntVect(AMREX_D_DECL(i,j,k)), fab, dcomp, numcomp, geomdata, time,
                   bcr_p, 0, orig_comp);
        });
    }

    // yhi
    if (!geom.isPeriodic(1) && bx.bigEnd(1) > domain.bigEnd(1)) {
        Box bndry = bx;
        int dyhi = domain.bigEnd(1);
        bndry.setSmall(1,dyhi+1);
        amrex::ParallelFor(bndry, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for (int n = dcomp; n < dcomp+numcomp; ++n) {
                fab(i,j,k,n) = fab(i,dyhi,k,n);
            }
            f_user(IntVect(AMREX_D_DECL(i,j,k)), fab, dcomp, numcomp, geomdata, time,
                   bcr_p, 0, orig_comp);
        });
    }
#endif

#if (AMREX_SPACEDIM == 3)
    // zlo
    if (!geom.isPeriodic(2) && bx.smallEnd(2) < domain.smallEnd(2)) {
        Box bndry = bx;
        int dzlo = domain.smallEnd(2);
        bndry.setBig(2,dzlo-1);
        amrex::ParallelFor(bndry, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for (int n = dcomp; n < dcomp+numcomp; ++n) {
                fab(i,j,k,n) = fab(i,j,dzlo,n);
            }
            f_user(IntVect(AMREX_D_DECL(i,j,k)), fab, dcomp, numcomp, geomdata, time,
                   bcr_p, 0, orig_comp);
        });
    }

    // zhi
    if (!geom.isPeriodic(2) && bx.bigEnd(2) > domain.bigEnd(2)) {
        Box bndry = bx;
        int dzhi = domain.bigEnd(2);
        bndry.setSmall(2,dzhi+1);
        amrex::ParallelFor(bndry, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for (int n = dcomp; n < dcomp+numcomp; ++n) {
                fab(i,j,k,n) = fab(i,j,dzhi,n);
            }
            f_user(IntVect(AMREX_D_DECL(i,j,k)), fab, dcomp, numcomp, geomdata, time,
                   bcr_p, 0, orig_comp);
        });
    }
#endif
}

}

#endif
