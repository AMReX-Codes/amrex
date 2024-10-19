#ifndef AMREX_FFT_HELPER_H_
#define AMREX_FFT_HELPER_H_
#include <AMReX_Config.H>

#include <AMReX_DistributionMapping.H>

namespace amrex::FFT
{

enum struct Direction { forward, backward, both };

struct Info
{
    //! Supported only in 3D. When batch_mode is true, FFT is performed on
    //! the first two dimensions only and the third dimension size is the
    //! batch size.
    bool batch_mode = false;

    Info& setBatchMode (bool x) { batch_mode = x; return *this; }
};

namespace detail
{
    DistributionMapping make_iota_distromap (Long n);
}

}

#endif
