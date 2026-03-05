#ifndef AMREX_BOUNDARY_FWD_H_
#define AMREX_BOUNDARY_FWD_H_

/**
 * \file AMReX_BoundaryFwd.H
 * \brief Forward declarations for boundary-related classes to break cyclic dependencies.
 */

namespace amrex {

template <typename MF> class BndryDataT;
template <typename MF> class BndryRegisterT;
class BoundCond;
template <typename MF> class FabSetT;
template <typename MF> class InterpBndryDataT;
class Mask;
class MultiMask;
template <typename MF> class YAFluxRegisterT;
class EdgeFluxRegister;

}

#endif
