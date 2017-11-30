#ifndef AMREX_ML_CELL_LINOP_H_H
#define AMREX_ML_CELL_LINOP_H_H

#include <AMReX_MLLinOp.H>

namespace amrex {

class MLCellLinOp
    : public MLLinOp
{
public:

    friend class MLMG;
    friend class MLCGSolver;

    MLCellLinOp ();
    virtual ~MLCellLinOp ();

    MLCellLinOp (const MLCellLinOp&) = delete;
    MLCellLinOp (MLCellLinOp&&) = delete;
    MLCellLinOp& operator= (const MLCellLinOp&) = delete;
    MLCellLinOp& operator= (MLCellLinOp&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo());

};

}

#endif
