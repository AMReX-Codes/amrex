#ifndef WARPX_FILTER_H_
#define WARPX_FILTER_H_

class BilinearFilter
{
public:
    BilinearFilter () = default;
    
    void ComputeStencils();
    void ApplyStencil(amrex::MultiFab& dstmf, const amrex::MultiFab& srcmf, int scomp=0, int dcomp=0, int ncomp=10000);
    void Filter2d(const amrex::Box& tbx, amrex::FArrayBox &tmpfab, amrex::FArrayBox &dstfab);
    void Filter2d_slow(const amrex::Box& tbx, amrex::FArrayBox &tmpfab, amrex::FArrayBox &dstfab);
    void Filter3d(const amrex::Box& tbx, amrex::FArrayBox &tmpfab, amrex::FArrayBox &dstfab);
    void Filter3d_slow(const amrex::Box& tbx, amrex::FArrayBox &tmpfab, amrex::FArrayBox &dstfab);

    amrex::IntVect npass_each_dir;
    amrex::IntVect stencil_length_each_dir;

private:

    amrex::Vector<amrex::Real> stencil_x, stencil_y, stencil_z;

};
#endif
