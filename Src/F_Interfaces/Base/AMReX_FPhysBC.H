#ifndef AMREX_FPHYSBC_H_
#define AMREX_FPHYSBC_H_

#include <AMReX_PhysBCFunct.H>

namespace amrex
{

class FPhysBC
{
public:
    using fill_physbc_funptr_t = void(*)(MultiFab* mf, int scomp, int ncomp, Real time, const Geometry* geom);

    FPhysBC (fill_physbc_funptr_t fill, const Geometry* geom_) : fill_physbc(fill), geom(geom_) {}

    ~FPhysBC () {}
    void operator() (MultiFab& mf, int scomp, int ncomp, IntVect const& nghost,
                     Real time, int bccomp);

private:
    fill_physbc_funptr_t fill_physbc = nullptr;
    const Geometry* geom = nullptr;
};

}
#endif
