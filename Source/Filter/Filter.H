#include <AMReX_MultiFab.H>
#include <AMReX_CudaContainers.H>

#ifndef WARPX_FILTER_H_
#define WARPX_FILTER_H_

class Filter
{
public:
    Filter () = default;

    virtual void ComputeStencils() = 0;
    void ApplyStencil(amrex::MultiFab& dstmf, 
                              const amrex::MultiFab& srcmf, int scomp=0, 
                              int dcomp=0, int ncomp=10000);

    // public for cuda
    void DoFilter(const amrex::Box& tbx,
                          amrex::Array4<amrex::Real const> const& tmp,
                          amrex::Array4<amrex::Real      > const& dst,
                          int scomp, int dcomp, int ncomp);

    amrex::IntVect stencil_length_each_dir;

protected:
    amrex::Dim3 slen;
    amrex::Gpu::ManagedVector<amrex::Real> stencil_x, stencil_y, stencil_z;

private:

};
#endif // #ifndef WARPX_FILTER_H_
