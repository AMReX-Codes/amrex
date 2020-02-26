/* Copyright 2019 David Grote, Maxence Thevenet, Weiqun Zhang
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX_py.H"


extern "C" {

    WARPX_CALLBACK_PY_FUNC_0 warpx_py_afterinit = nullptr;
    WARPX_CALLBACK_PY_FUNC_0 warpx_py_beforeEsolve = nullptr;
    WARPX_CALLBACK_PY_FUNC_0 warpx_py_afterEsolve = nullptr;
    WARPX_CALLBACK_PY_FUNC_0 warpx_py_beforedeposition = nullptr;
    WARPX_CALLBACK_PY_FUNC_0 warpx_py_afterdeposition = nullptr;
    WARPX_CALLBACK_PY_FUNC_0 warpx_py_particlescraper = nullptr;
    WARPX_CALLBACK_PY_FUNC_0 warpx_py_particleloader = nullptr;
    WARPX_CALLBACK_PY_FUNC_0 warpx_py_beforestep = nullptr;
    WARPX_CALLBACK_PY_FUNC_0 warpx_py_afterstep = nullptr;
    WARPX_CALLBACK_PY_FUNC_0 warpx_py_afterrestart = nullptr;
    WARPX_CALLBACK_PY_FUNC_0 warpx_py_particleinjection = nullptr;
    WARPX_CALLBACK_PY_FUNC_0 warpx_py_appliedfields = nullptr;

}

