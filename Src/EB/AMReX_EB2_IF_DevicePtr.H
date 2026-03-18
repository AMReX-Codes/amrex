#ifndef AMREX_EB2_IF_DEVICEPTR_H_
#define AMREX_EB2_IF_DEVICEPTR_H_

/**
 * \file AMReX_EB2_IF_DevicePtr.H
 *
 * Utility that wraps an implicit function referenced through a raw device pointer.
 */

namespace amrex::EB2 {

//! Wrapper that forwards operator() calls to an implicit function pointer.
template <typename T>
class DevicePtrIF
    : public GPUable
{
public:

    //! Store the pointer that owns the actual implicit-function operator.
    DevicePtrIF(T* a_p) : p(a_p) {}

    //! Device-friendly overload dispatching to the pointee.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        return (*p)(AMREX_D_DECL(x,y,z));
    }

    //! Array-based helper overload.
    [[nodiscard]] Real operator() (const RealArray& a) const noexcept
    {
        return this->operator() (AMREX_D_DECL(a[0], a[1], a[2]));
    }

protected:
    T* p = nullptr;
};

}

#endif
