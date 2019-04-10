
#ifndef AMREX_BOUNDCOND_H_
#define AMREX_BOUNDCOND_H_

#include <AMReX_BLassert.H>
#include <AMReX_GpuQualifiers.H>

namespace amrex {

/**
   \brief Maintain an identifier for boundary condition types.

   This is a placeholder for more extensive boundary condition
   implementations, which might include stencils, etc.  Presently,
   boundary conditions are specified via an integer identifier.
   This class maintains that integer.
*/

class BoundCond
{
public:
    //
    //! Initializes boundary condition type to default.
    BoundCond () noexcept : bctype(-1) {}
    //
    //! Initializes boundary condition type to specified value.
    BoundCond (int _bctype) noexcept : bctype(_bctype) {}
    //
    //! Return index of boundary condition type.
    AMREX_GPU_HOST_DEVICE
    operator int () const noexcept
    {
        BL_ASSERT(bctype != -1);
        return bctype;
    }
    //
    //! Set index of boundary condition type.
    BoundCond& operator= (const int _bctype) noexcept
    {
        bctype = _bctype;
        return *this;
    }
private:
    //
    //! Index of boundary condition type.
    int bctype;
};

}

#endif /*_BOUNDCOND_H_*/
