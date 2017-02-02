#ifndef AMREX_FAMRCORE_H_
#define AMREX_FAMRCORE_H_

#include <AMReX_AmrCore.H>

namespace amrex {

class FAmrCore
    : public AmrCore
{
public:
    FAmrCore ();
    virtual ~FAmrCore ();

protected:

    virtual void ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow) override;
    
};

}

#endif

