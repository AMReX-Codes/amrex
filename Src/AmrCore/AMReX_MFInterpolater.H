#ifndef AMREX_MF_INTERPOLATER_H_
#define AMREX_MF_INTERPOLATER_H_
#include <AMReX_Config.H>

#include <AMReX_Extension.H>
#include <AMReX_InterpBase.H>

namespace amrex {

class MultiFab;
class Geometry;

class MFInterpolater
    : public InterpBase
{
public:
    virtual ~MFInterpolater () = default;

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
    virtual ~MFPCInterp () = default;

    virtual Box CoarseBox (Box const& fine, int ratio);
    virtual Box CoarseBox (Box const& fine, IntVect const& ratio);

    virtual void interp (MultiFab const& crsemf, int ccomp, MultiFab& finemf, int fcomp, int ncomp,
                         IntVect const& ng, Geometry const& cgeom, Geometry const& fgeom,
                         Box const& dest_domain, IntVect const& ratio,
                         Vector<BCRec> const& bcs, int bcscomp);
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

    virtual ~MFCellConsLinInterp () = default;

    virtual Box CoarseBox (Box const& fine, int ratio) override;
    virtual Box CoarseBox (Box const& fine, IntVect const& ratio) override;

    virtual void interp (MultiFab const& crsemf, int ccomp, MultiFab& finemf, int fcomp, int ncomp,
                         IntVect const& ng, Geometry const& cgeom, Geometry const& fgeom,
                         Box const& dest_domain, IntVect const& ratio,
                         Vector<BCRec> const& bcs, int bcscomp) override;
protected:
    bool do_linear_limiting = true;
};

/**
 * \brief [Bi|Tri]linear interpolation on cell centered data.
 */
class MFCellBilinear final
    : public MFInterpolater
{
public:
    virtual ~MFCellBilinear () = default;

    virtual Box CoarseBox (Box const& fine, int ratio);
    virtual Box CoarseBox (Box const& fine, IntVect const& ratio);

    virtual void interp (MultiFab const& crsemf, int ccomp, MultiFab& finemf, int fcomp, int ncomp,
                         IntVect const& ng, Geometry const& cgeom, Geometry const& fgeom,
                         Box const& dest_domain, IntVect const& ratio,
                         Vector<BCRec> const& bcs, int bcscomp);
};

/*
 * \brief [Bi|Tri] linear interpolation on nodal data
 */
class MFNodeBilinear final
    : public MFInterpolater
{
public:
    virtual ~MFNodeBilinear () = default;

    virtual Box CoarseBox (Box const& fine, int ratio);
    virtual Box CoarseBox (Box const& fine, IntVect const& ratio);

    virtual void interp (MultiFab const& crsemf, int ccomp, MultiFab& finemf, int fcomp, int ncomp,
                         IntVect const& ng, Geometry const& cgeom, Geometry const& fgeom,
                         Box const& dest_domain, IntVect const& ratio,
                         Vector<BCRec> const& bcs, int bcscomp);
};

extern AMREX_EXPORT MFPCInterp          mf_pc_interp;
extern AMREX_EXPORT MFCellConsLinInterp mf_cell_cons_interp;
extern AMREX_EXPORT MFCellConsLinInterp mf_lincc_interp;
extern AMREX_EXPORT MFCellBilinear      mf_cell_bilinear_interp;
extern AMREX_EXPORT MFNodeBilinear      mf_node_bilinear_interp;

}

#endif
