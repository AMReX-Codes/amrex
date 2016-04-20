#ifndef _PhysBCFunct_H_
#define _PhysBCFunct_H_

#include <MultiFab.H>

class PhysBCFunct {
public:
    virtual void doit (MultiFab& mf, int dcomp, int ncomp, Real time) = 0;
    virtual ~PhysBCFunct () = 0;
};

inline PhysBCFunct::~PhysBCFunct () { }

#endif
