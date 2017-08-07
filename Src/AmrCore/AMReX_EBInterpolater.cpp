
#include <AMReX_EBInterpolater.H>

namespace amrex {

EBCellConservativeLinear  eb_lincc_interp;
EBCellConservativeLinear  eb_cell_cons_interp(0);

EBCellConservativeLinear::EBCellConservativeLinear (bool do_linear_limiting_)
    : CellConservativeLinear(do_linear_limiting_)
{
}

EBCellConservativeLinear::~EBCellConservativeLinear ()
{

}

void
EBCellConservativeLinear::interp (const FArrayBox& crse,
                                  int              crse_comp,
                                  FArrayBox&       fine,
                                  int              fine_comp,
                                  int              ncomp,
                                  const Box&       fine_region,
                                  const IntVect&   ratio,
                                  const Geometry&  crse_geom,
                                  const Geometry&  fine_geom,
                                  Array<BCRec>&    bcr,
                                  int              actual_comp,
                                  int              actual_state)
{

}

}
