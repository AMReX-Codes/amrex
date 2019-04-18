#include <AMReX_MultiFab.H>
#include <AMReX_CudaContainers.H>

#ifndef WARPX_FILTER_H_
#define WARPX_FILTER_H_

class BilinearFilter
{
public:
    BilinearFilter () = default;
    
    void ComputeStencils();
    void ApplyStencil(amrex::MultiFab& dstmf, const amrex::MultiFab& srcmf, int scomp=0, int dcomp=0, int ncomp=10000);

    amrex::IntVect npass_each_dir;
    amrex::IntVect stencil_length_each_dir;

    void Filter(const amrex::Box& tbx, amrex::FArrayBox const& tmpfab, amrex::FArrayBox &dstfab, int scomp, int dcomp, int ncomp);

private:

    amrex::Gpu::ManagedVector<amrex::Real> stencil_x, stencil_y, stencil_z;

    amrex::Dim3 slen;
};
#endif
