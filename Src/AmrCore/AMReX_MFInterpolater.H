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

class MFCellConsLinInterp final
    : public MFInterpolater
{
public:
    explicit MFCellConsLinInterp (bool do_linear_limiting_)
        : do_linear_limiting(do_linear_limiting_) {}

    virtual ~MFCellConsLinInterp () = default;

    virtual Box CoarseBox (Box const& fine, int ratio);
    virtual Box CoarseBox (Box const& fine, IntVect const& ratio);

    virtual void interp (MultiFab const& crsemf, int ccomp, MultiFab& finemf, int fcomp, int ncomp,
                         IntVect const& ng, Geometry const& cgeom, Geometry const& fgeom,
                         Box const& dest_domain, IntVect const& ratio,
                         Vector<BCRec> const& bcs, int bcscomp);
protected:
    bool do_linear_limiting = true;
};

extern AMREX_EXPORT MFCellConsLinInterp mf_cell_cons_interp;
extern AMREX_EXPORT MFCellConsLinInterp mf_lincc_interp;

}

#endif
