/* Copyright 2019 David Grote, Maxence Thevenet, Weiqun Zhang
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PY_H_
#define WARPX_PY_H_

#include "WarpXWrappers.h"


extern "C" {

    extern WARPX_CALLBACK_PY_FUNC_0 warpx_py_afterinit;
    extern WARPX_CALLBACK_PY_FUNC_0 warpx_py_beforeEsolve;
    extern WARPX_CALLBACK_PY_FUNC_0 warpx_py_afterEsolve;
    extern WARPX_CALLBACK_PY_FUNC_0 warpx_py_beforedeposition;
    extern WARPX_CALLBACK_PY_FUNC_0 warpx_py_afterdeposition;
    extern WARPX_CALLBACK_PY_FUNC_0 warpx_py_particlescraper;
    extern WARPX_CALLBACK_PY_FUNC_0 warpx_py_particleloader;
    extern WARPX_CALLBACK_PY_FUNC_0 warpx_py_beforestep;
    extern WARPX_CALLBACK_PY_FUNC_0 warpx_py_afterstep;
    extern WARPX_CALLBACK_PY_FUNC_0 warpx_py_afterrestart;
    extern WARPX_CALLBACK_PY_FUNC_0 warpx_py_particleinjection;
    extern WARPX_CALLBACK_PY_FUNC_0 warpx_py_appliedfields;

}

#endif
