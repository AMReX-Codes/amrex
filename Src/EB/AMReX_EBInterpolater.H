#ifndef AMREX_EBINTERPOLATER_H_
#define AMREX_EBINTERPOLATER_H_
#include <AMReX_Config.H>

#include <AMReX_Interpolater.H>

/**
 * \file AMReX_EBInterpolater.H
 *
 * EB-aware conservative cell interpolater used for AMR prolongation.
 */

namespace amrex {

//! CellConservativeLinear interpolater customized for EB volume fractions.
class EBCellConservativeLinear
    : public CellConservativeLinear
{
public:

    using CellConservativeLinear::CellConservativeLinear;

    /**
     * \brief Conservatively interpolate coarse data \p crse into fine FAB \p fine.
     *
     * Identical signature to CellConservativeLinear::interp but respects EB metadata.
     */
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
