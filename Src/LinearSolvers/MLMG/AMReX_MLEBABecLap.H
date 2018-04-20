#ifndef AMREX_MLEBABECLAP_H_
#define AMREX_MLEBABECLAP_H_

#include <AMReX_EBFabFactory.H>
#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

namespace amrex {

// (alpha * a - beta * (del dot b grad)) phi

class MLEBABecLap
    : public MLCellLinOp
{
public:

    MLEBABecLap () {}
    MLEBABecLap (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory);

    virtual ~MLEBABecLap ();

    MLEBABecLap (const MLEBABecLap&) = delete;
    MLEBABecLap (MLEBABecLap&&) = delete;
    MLEBABecLap& operator= (const MLEBABecLap&) = delete;
    MLEBABecLap& operator= (MLEBABecLap&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory);

protected:

private:

};

}

#endif
