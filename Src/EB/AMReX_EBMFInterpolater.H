#ifndef AMREX_EBMFINTERPLATER_H_
#define AMREX_EBMFINTERPLATER_H_
#include <AMReX_Config.H>

#include <AMReX_MFInterpolater.H>

/**
 * \file AMReX_EBMFInterpolater.H
 *
 * EB-aware MultiFab interpolater for conservative cell data.
 */

namespace amrex {

//! Cell-conservative MultiFab interpolater that honors EB geometry.
class EBMFCellConsLinInterp final
    : public MFCellConsLinInterp
{
public:
    using MFCellConsLinInterp::MFCellConsLinInterp;

    //! Interpolate coarse MultiFab data \p crsemf into \p finemf with EB awareness.
    void interp (MultiFab const& crsemf, int ccomp, MultiFab& finemf, int fcomp, int ncomp,
                         IntVect const& ng, Geometry const& cgeom, Geometry const& fgeom,
                         Box const& dest_domain, IntVect const& ratio,
                         Vector<BCRec> const& bcs, int bcomp) override;
};

//! EB-aware MultiFab cell-conservative interpolater singleton.
extern AMREX_EXPORT EBMFCellConsLinInterp eb_mf_cell_cons_interp;
//! Alias for backward compatibility with older naming.
extern AMREX_EXPORT EBMFCellConsLinInterp eb_mf_lincc_interp;

}

#endif
