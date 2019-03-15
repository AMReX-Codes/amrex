#ifndef WARPX_FILTER_H_
#define WARPX_FILTER_H_

class BilinearFilter
{
public:
    BilinearFilter () = default;
    
    void ComputeStencils();
    void ApplyStencils();
    amrex::IntVect npass_each_dir;

private:

    amrex::Vector<amrex::Real> stencil_x, stencil_y, stencil_z;

};
#endif
