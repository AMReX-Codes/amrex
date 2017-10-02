#ifndef AmrAdvPhysBC_H_
#define AmrAdvPhysBC_H_

#include <AMReX_PhysBCFunct.H>
#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

using namespace amrex;

class AmrCoreAdvPhysBC
    : public PhysBCFunctBase
{
public:
    AmrCoreAdvPhysBC () {}
    virtual ~AmrCoreAdvPhysBC () {}
    AmrCoreAdvPhysBC (const Geometry& geom, const Array<BCRec>& bcr);
    void define (const Geometry& geom, const Array<BCRec>& bcr);
    virtual void FillBoundary (MultiFab& mf, int dcomp, int ncomp, Real time) override;
    using PhysBCFunctBase::FillBoundary;
private:
    Geometry     m_geom;
    Array<BCRec> m_bcr;
};

#endif
