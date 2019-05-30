#include <Filter.H>

#ifndef WARPX_BILIN_FILTER_H_
#define WARPX_BILIN_FILTER_H_

class BilinearFilter : public Filter
{
public:

    BilinearFilter() = default;
    
    void ComputeStencils();

    amrex::IntVect npass_each_dir;

};
#endif // #ifndef WARPX_BILIN_FILTER_H_
