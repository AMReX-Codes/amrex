#ifndef AMREX_ML_CELL_ABECLAP_H_
#define AMREX_ML_CELL_ABECLAP_H_

#include <AMReX_MLCellLinOp.H>

namespace amrex {

class MLCellABecLap
    : public MLCellLinOp
{
public:

    friend class MLMG;
    friend class MLCGSolver;

    MLCellABecLap ();
    virtual ~MLCellABecLap ();

    MLCellABecLap (const MLCellABecLap&) = delete;
    MLCellABecLap (MLCellABecLap&&) = delete;
    MLCellABecLap& operator= (const MLCellABecLap&) = delete;
    MLCellABecLap& operator= (MLCellABecLap&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

    virtual bool needsUpdate () const override {
        return MLCellLinOp::needsUpdate();
    }
    virtual void update () override;

protected:

    virtual void prepareForSolve () override;

};

}

#endif
