#ifndef MY_AMR_H_
#define MY_AMR_H_

#include <AMReX_AmrMesh.H>

class MyAmr
    : public amrex::AmrMesh
{
public:
    using amrex::AmrMesh::AmrMesh;  // inherit AmrMesh constructors

    virtual ~MyAmr() {}

    //! Manually tag.  Note that tags is built on level lev grids coarsened by bf_lev.
    virtual void ManualTagsPlacement (int lev, amrex::TagBoxArray& tags,
                                      const amrex::Vector<amrex::IntVect>& bf_lev) override;

    // We could also override 
    //! Tag cells for refinement.  TagBoxArray tags is built on level lev grids.
    // virtual void ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow) override;
};

#endif
