#ifndef AMREX_EBMFINTERPLATER_H_
#define AMREX_EBMFINTERPLATER_H_
#include <AMReX_Config.H>

#include <AMReX_MFInterpolater.H>

namespace amrex {

class EBMFCellConsLinInterp final
    : public MFCellConsLinInterp
{
public:
    using MFCellConsLinInterp::MFCellConsLinInterp;

    void interp (MultiFab const& crsemf, int ccomp, MultiFab& finemf, int fcomp, int ncomp,
                         IntVect const& ng, Geometry const& cgeom, Geometry const& fgeom,
                         Box const& dest_domain, IntVect const& ratio,
                         Vector<BCRec> const& bcs, int bcomp) override;
};

extern AMREX_EXPORT EBMFCellConsLinInterp eb_mf_cell_cons_interp;
extern AMREX_EXPORT EBMFCellConsLinInterp eb_mf_lincc_interp;

}

#endif
