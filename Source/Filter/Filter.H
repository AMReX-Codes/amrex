#include <AMReX_MultiFab.H>
#include <AMReX_CudaContainers.H>

#ifndef WARPX_FILTER_H_
#define WARPX_FILTER_H_

class Filter
{
public:
    Filter () = default;

    // Apply stencil on MultiFab.
    // Guard cells are handled inside this function
    void ApplyStencil(amrex::MultiFab& dstmf, 
                              const amrex::MultiFab& srcmf, int scomp=0, 
                              int dcomp=0, int ncomp=10000);

    // Apply stencil on a FabArray.
    void ApplyStencil (amrex::FArrayBox& dstfab, 
                       const amrex::FArrayBox& srcfab, const amrex::Box& tbx, 
                       int scomp=0, int dcomp=0, int ncomp=10000);

    // public for cuda
    void DoFilter(const amrex::Box& tbx,
                          amrex::Array4<amrex::Real const> const& tmp,
                          amrex::Array4<amrex::Real      > const& dst,
                          int scomp, int dcomp, int ncomp);

    // In 2D, stencil_length_each_dir = {length(stencil_x), length(stencil_z)}
    amrex::IntVect stencil_length_each_dir;

protected:
    // Stencil along each direction.
    // in 2D, stencil_y is not initialized.
    amrex::Gpu::ManagedVector<amrex::Real> stencil_x, stencil_y, stencil_z;
    // Length of each stencil.
    // In 2D, slen = {length(stencil_x), length(stencil_z), 1}
    amrex::Dim3 slen;

private:

};
#endif // #ifndef WARPX_FILTER_H_
