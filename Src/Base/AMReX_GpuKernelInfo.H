#ifndef AMREX_GPU_KERNEL_INFO_H_
#define AMREX_GPU_KERNEL_INFO_H_

namespace amrex {
namespace Gpu {

class KernelInfo
{
public:
    KernelInfo& setFusable (bool flag) { fusable = flag; return *this; }
    bool isFusable () const { return fusable; }
private:
    bool fusable = false;
};

}}

#endif
