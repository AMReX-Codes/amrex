#ifndef _PhysBCFunct_H_
#define _PhysBCFunct_H_

#include <MultiFab.H>

class PhysBCFunctBase {
public:
    virtual void FillBoundary (MultiFab& mf, int dcomp, int ncomp, Real time) = 0;
    virtual ~PhysBCFunctBase () = 0;
};

inline PhysBCFunctBase::~PhysBCFunctBase () { }

#endif
