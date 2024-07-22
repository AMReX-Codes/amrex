#ifndef AMREX_EBINTERPOLATER_H_
#define AMREX_EBINTERPOLATER_H_
#include <AMReX_Config.H>

#include <AMReX_Interpolater.H>

namespace amrex {

class EBCellConservativeLinear
    : public CellConservativeLinear
{
public:

    using CellConservativeLinear::CellConservativeLinear;

    void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec> const& bcr,
                         int              actual_comp,
                         int              actual_state,
                         RunOn            runon) override;
};

extern AMREX_EXPORT EBCellConservativeLinear  eb_lincc_interp;
extern AMREX_EXPORT EBCellConservativeLinear  eb_cell_cons_interp;

}

#endif
