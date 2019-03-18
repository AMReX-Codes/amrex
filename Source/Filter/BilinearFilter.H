#ifndef WARPX_FILTER_H_
#define WARPX_FILTER_H_

class BilinearFilter
{
public:
    BilinearFilter () = default;
    
    void ComputeStencils();
    // void ApplyStencils();
    //void ApplyStencils(amrex::MultiFab* dstmf, amrex::MultiFab* srcmf);
    void ApplyStencil (amrex::MultiFab& dstmf, const amrex::MultiFab& srcmf, int scomp = 0, int dcomp = 0, int ncomp = 10000);
    void filter_2d(const amrex::Box& tbx, const amrex::Box& gbx, amrex::FArrayBox &tmpfab, amrex::FArrayBox &dstfab, int ncomp);

    amrex::IntVect npass_each_dir;
    amrex::IntVect stencil_length_each_dir;

private:

    amrex::Vector<amrex::Real> stencil_x, stencil_y, stencil_z;

};
#endif
