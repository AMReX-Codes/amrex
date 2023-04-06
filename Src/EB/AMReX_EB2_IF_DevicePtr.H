#ifndef AMREX_EB2_IF_DEVICEPTR_H_
#define AMREX_EB2_IF_DEVICEPTR_H_

namespace amrex::EB2 {

template <typename T>
class DevicePtrIF
    : public GPUable
{
public:

    DevicePtrIF(T* a_p) : p(a_p) {}

    AMREX_GPU_HOST_DEVICE
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        return (*p)(AMREX_D_DECL(x,y,z));
    }

    Real operator() (const RealArray& a) const noexcept
    {
        return this->operator() (AMREX_D_DECL(a[0], a[1], a[2]));
    }

protected:
    T* p = nullptr;
};

}

#endif
