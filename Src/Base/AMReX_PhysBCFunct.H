#ifndef _PhysBCFunct_H_
#define _PhysBCFunct_H_

#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

extern "C"
{
    typedef void (*BndryFuncDefault)(amrex::Real* data, ARLIM_P(lo), ARLIM_P(hi),
                                     const int* dom_lo, const int* dom_hi,
                                     const amrex::Real* dx, const amrex::Real* grd_lo,
                                     const amrex::Real* time, const int* bc);
    typedef void (*BndryFunc3DDefault)(amrex::Real* data, const int* lo, const int* hi,
                                       const int* dom_lo, const int* dom_hi,
                                       const amrex::Real* dx, const amrex::Real* grd_lo,
                                       const amrex::Real* time, const int* bc);
}

namespace amrex {

class BndryFunctBase
{
public:
    BndryFunctBase ();
    BndryFunctBase (BndryFuncDefault inFunc);
    BndryFunctBase (BndryFunc3DDefault inFunc);
    virtual ~BndryFunctBase ();
    /**
    * \brief Return a ptr to a clone of this object.
    * It is the responsibility of the caller to delete the result.
    */    
    virtual BndryFunctBase* clone () const;
    //! Fill boundary cells
    virtual void operator () (Real* data, const int* lo, const int* hi,
			      const int* dom_lo, const int* dom_hi,
			      const Real* dx, const Real* grd_lo,
			      const Real* time, const int* bc) const;
protected:
    BndryFuncDefault    m_func;
    BndryFunc3DDefault  m_func3D;
};

class PhysBCFunctBase {
public:
    virtual void FillBoundary (MultiFab& mf, int dcomp, int ncomp, Real time) = 0;
    virtual void FillBoundary (MultiFab& mf, Real time)
	{ FillBoundary(mf, 0, mf.nComp(), time); }
    virtual ~PhysBCFunctBase () = 0;
};

inline PhysBCFunctBase::~PhysBCFunctBase () {}

class PhysBCFunct
    : public PhysBCFunctBase
{
public:
    PhysBCFunct () {}
    virtual ~PhysBCFunct () {}
    PhysBCFunct (const Geometry& geom, const BCRec& bcr, const BndryFunctBase& func);
    void define (const Geometry& geom, const BCRec& bcr, const BndryFunctBase& func);
    virtual void FillBoundary (MultiFab& mf, int, int, Real time) override;
    using PhysBCFunctBase::FillBoundary;
private:
    Geometry m_geom;
    BCRec    m_bcr;
    std::unique_ptr<BndryFunctBase> m_bc_func;
};

}

#endif
