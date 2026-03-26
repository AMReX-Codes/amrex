#ifndef AMREX_AMRCORE_FWD_H_
#define AMREX_AMRCORE_FWD_H_

/**
 * \file AMReX_AmrCoreFwd.H
 * \brief Forward declarations for AMR core types to help break include cycles.
 */

namespace amrex {

class AmrCore;
struct AmrInfo;
class AmrMesh;

class FluxRegister;
class Interpolater;
class MFInterpolater;

class TagBox;
class TagBoxArray;

class InterpFaceRegister;

}

#endif
