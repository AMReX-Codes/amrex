#ifndef AMREX_GPU_KERNEL_INFO_H_
#define AMREX_GPU_KERNEL_INFO_H_
#include <AMReX_Config.H>

namespace amrex {
namespace Gpu {

class KernelInfo
{
public:
    KernelInfo& setReduction (bool flag) { has_reduction = flag; return *this; }
    bool hasReduction () const { return has_reduction; }
private:
    bool has_reduction = false;
};

}}

#endif
