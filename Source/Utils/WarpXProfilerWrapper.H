/* Copyright 2020 Axel Huebl, Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_PROFILERWRAPPER_H_
#define WARPX_PROFILERWRAPPER_H_

#include "AMReX_BLProfiler.H"
#include "AMReX_GpuDevice.H"

AMREX_FORCE_INLINE
void doDeviceSynchronize(int do_device_synchronize)
{
    if ( do_device_synchronize )
        amrex::Gpu::synchronize();
}

#define WARPX_PROFILE(fname) doDeviceSynchronize(WarpX::do_device_synchronize_before_profile); BL_PROFILE(fname)
#define WARPX_PROFILE_VAR(fname, vname) doDeviceSynchronize(WarpX::do_device_synchronize_before_profile); BL_PROFILE_VAR(fname, vname)
#define WARPX_PROFILE_VAR_NS(fname, vname) doDeviceSynchronize(WarpX::do_device_synchronize_before_profile); BL_PROFILE_VAR_NS(fname, vname)
#define WARPX_PROFILE_VAR_START(vname) doDeviceSynchronize(WarpX::do_device_synchronize_before_profile); BL_PROFILE_VAR_START(vname)
#define WARPX_PROFILE_VAR_STOP(vname) doDeviceSynchronize(WarpX::do_device_synchronize_before_profile); BL_PROFILE_VAR_STOP(vname)
#define WARPX_PROFILE_REGION(rname) doDeviceSynchronize(WarpX::do_device_synchronize_before_profile); BL_PROFILE_REGION(rname)

#endif // WARPX_PROFILERWRAPPER_H_
