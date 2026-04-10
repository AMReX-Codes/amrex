#ifndef AMREX_MF_INTERPOLATER_H_
#define AMREX_MF_INTERPOLATER_H_
#include <AMReX_Config.H>

#include <AMReX_Extension.H>
#include <AMReX_InterpBase.H>

/**
 * \file AMReX_MFInterpolater.H
 * \brief Interpolaters that operate directly on MultiFabs (cell- or node-centered).
 */

namespace amrex {

class MultiFab;
class Geometry;

class MFInterpolater
    : public InterpBase
{
public:
    /**
     * \brief Interpolate coarse MultiFab data onto a fine MultiFab region.
     *
     * \param crsemf      Source coarse MultiFab.
     * \param ccomp       Starting component in the coarse MF.
     * \param finemf      Destination fine MultiFab.
     * \param fcomp       Starting component in the fine MF.
     * \param ncomp       Number of components to interpolate.
     * \param ng          Number of grow cells to fill on the fine MF.
     * \param cgeom       Coarse Geometry.
     * \param fgeom       Fine Geometry.
     * \param dest_domain Region of the fine domain that should be filled.
     * \param ratio       Refinement ratio.
     * \param bcs         Boundary-condition records per component.
     * \param bcscomp     Starting component offset inside the BC array.
     */
    virtual void interp (MultiFab const& crsemf, int ccomp, MultiFab& finemf, int fcomp, int ncomp,
                         IntVect const& ng, Geometry const& cgeom, Geometry const& fgeom,
                         Box const& dest_domain, IntVect const& ratio,
                         Vector<BCRec> const& bcs, int bcscomp) = 0;
};

/**
* \brief Piecewise constant interpolation on cell-centered data.
*/
class MFPCInterp final
    : public MFInterpolater
{
public:
    Box CoarseBox (Box const& fine, int ratio) override;
    Box CoarseBox (Box const& fine, IntVect const& ratio) override;

    void interp (MultiFab const& crsemf, int ccomp, MultiFab& finemf, int fcomp, int ncomp,
                         IntVect const& ng, Geometry const& cgeom, Geometry const& fgeom,
                         Box const& dest_domain, IntVect const& ratio,
                         Vector<BCRec> const& bcs, int bcscomp) override;
};

/**
* \brief Linear conservative interpolation on cell centered data
*
* Linear conservative interpolation on cell centered data, i.e, conservative
* interpolation with a limiting scheme that preserves the value of any
* linear combination of the fab components. If sum_ivar
* a(ic,jc,ivar)*fab(ic,jc,ivar) = 0, then sum_ivar
* a(ic,jc,ivar)*fab(if,jf,ivar) = 0 is satisfied in all fine cells if,jf
* covering coarse cell ic,jc.
*/
class MFCellConsLinInterp
    : public MFInterpolater
{
public:
    explicit MFCellConsLinInterp (bool do_linear_limiting_)
        : do_linear_limiting(do_linear_limiting_) {}

    Box CoarseBox (Box const& fine, int ratio) override;
    Box CoarseBox (Box const& fine, IntVect const& ratio) override;

    void interp (MultiFab const& crsemf, int ccomp, MultiFab& finemf, int fcomp, int ncomp,
                         IntVect const& ng, Geometry const& cgeom, Geometry const& fgeom,
                         Box const& dest_domain, IntVect const& ratio,
                         Vector<BCRec> const& bcs, int bcomp) override;
protected:
    bool do_linear_limiting = true;
};

/**
* \brief Linear conservative interpolation on cell centered data
*
* Linear conservative interpolation on cell centered data, i.e, conservative
* interpolation with a limiting scheme that preserves the value of any
* linear combination of the fab components. If sum_ivar
* a(ic,jc,ivar)*fab(ic,jc,ivar) = 0, then sum_ivar
* a(ic,jc,ivar)*fab(if,jf,ivar) = 0 is satisfied in all fine cells if,jf
* covering coarse cell ic,jc. This version *also* ensures no new min/max
* are created by limiting the slope vector.
*/
class MFCellConsLinMinmaxLimitInterp
    : public MFInterpolater
{
public:
    Box CoarseBox (Box const& fine, int ratio) override;
    Box CoarseBox (Box const& fine, IntVect const& ratio) override;

    void interp (MultiFab const& crsemf, int ccomp, MultiFab& finemf, int fcomp, int ncomp,
                         IntVect const& ng, Geometry const& cgeom, Geometry const& fgeom,
                         Box const& dest_domain, IntVect const& ratio,
                         Vector<BCRec> const& bcs, int bcomp) override;
};

/**
 * \brief [Bi|Tri]linear interpolation on cell centered data.
 */
class MFCellBilinear final
    : public MFInterpolater
{
public:
    Box CoarseBox (Box const& fine, int ratio) override;
    Box CoarseBox (Box const& fine, IntVect const& ratio) override;

    void interp (MultiFab const& crsemf, int ccomp, MultiFab& finemf, int fcomp, int ncomp,
                         IntVect const& ng, Geometry const& cgeom, Geometry const& fgeom,
                         Box const& dest_domain, IntVect const& ratio,
                         Vector<BCRec> const& bcs, int bcscomp) override;
};

/*
 * \brief [Bi|Tri] linear interpolation on nodal data
 */
class MFNodeBilinear final
    : public MFInterpolater
{
public:
    Box CoarseBox (Box const& fine, int ratio) override;
    Box CoarseBox (Box const& fine, IntVect const& ratio) override;

    void interp (MultiFab const& crsemf, int ccomp, MultiFab& finemf, int fcomp, int ncomp,
                         IntVect const& ng, Geometry const& cgeom, Geometry const& fgeom,
                         Box const& dest_domain, IntVect const& ratio,
                         Vector<BCRec> const& bcs, int bcscomp) override;
};

extern AMREX_EXPORT MFPCInterp          mf_pc_interp;
extern AMREX_EXPORT MFCellConsLinInterp mf_cell_cons_interp;
extern AMREX_EXPORT MFCellConsLinInterp mf_lincc_interp;
extern AMREX_EXPORT MFCellConsLinMinmaxLimitInterp mf_linear_slope_minmax_interp;
extern AMREX_EXPORT MFCellBilinear      mf_cell_bilinear_interp;
extern AMREX_EXPORT MFNodeBilinear      mf_node_bilinear_interp;

}

#endif
