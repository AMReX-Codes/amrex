#ifndef AMREX_BC_TYPES_H_
#define AMREX_BC_TYPES_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_BC_TYPES.H
 * \brief Enumerations describing physical and mathematical boundary types.
 */

#if 0
---------------------------------------------------------------------
          CONVERSIONS FROM PHYSICAL TO MATHEMATICAL BCs
  Note that it is up to application codes to enforce conversions


PHYS BC        |   Velocity      |    Temperature   |    Scalars     |
----------------------------------------------------------------------
INTERIOR       |                 |                  |                |
   or          |    INT_DIR      |     INT_DIR      |    INT_DIR     |
PERIODIC       |                 |                  |                |
----------------------------------------------------------------------
               |                 |                  |                |
INFLOW         |    EXT_DIR      |     EXT_DIR      |    EXT_DIR     |
               |                 |                  |                |
----------------------------------------------------------------------
               |                 |                  |                |
OUTFLOW        |    FOEXTRAP     |     FOEXTRAP     |    FOEXTRAP    |
               |                 |                  |                |
----------------------------------------------------------------------
NO SLIP WALL   |                 |                  |                |
  with         |    EXT_DIR      |  REFLECT_EVEN    |    FOEXTRAP    |
ADIABATIC TEMP |    U = V = 0    |    dT/dn = 0     |                |
----------------------------------------------------------------------
NO SLIP WALL   |                 |                  |                |
  with         |    EXT_DIR      |     EXT_DIR      |    FOEXTRAP    |
FIXED TEMP     |    U = V = 0    |                  |                |
----------------------------------------------------------------------
SLIP WALL      |                 |                  |                |
  with         |  Un=0 EXT_DIR   |  REFLECT_EVEN    |    HOEXTRAP    |
ADIABATIC TEMP |  Ut   HOEXTRAP  |    dT/dn = 0     |                |
----------------------------------------------------------------------
SLIP WALL      |                 |                  |                |
  with         |  Un=0 EXT_DIR   |     EXT_DIR      |    HOEXTRAP    |
FIXED TEMP     |  Ut   HOEXTRAP  |                  |                |
----------------------------------------------------------------------
               |                 |                  |                |
SYMMETRY       |  Un REFLECT_ODD |  REFLECT_EVEN    | REFLECT_EVEN   |
               |  Ut REFLECT_EVEN|                  |                |
----------------------------------------------------------------------


INT_DIR      : data taken from other grids or interpolated
EXT_DIR      : data specified on EDGE (FACE)
EXT_DIR_CC   : data specified at cell center
HOEXTRAP     : higher order extrapolation to EDGE of bndry
HOEXTRAPCC   : linear extrapolation to CELL of bndry
FOEXTRAP     : first order extrapolation from last cell in interior
REFLECT_EVEN : F(-n) = F(n)   true reflection from interior cells
REFLECT_ODD  : F(-n) = -F(n)  true reflection from interior cells

#endif

#ifdef __cplusplus
namespace amrex {

/**
 * \namespace amrex::PhysBCType
 * \brief Physical boundary condition types assigned to domain faces.
 * \ingroup amrex_utilities
 */
namespace PhysBCType {
enum physicalBndryTypes : int {interior=0,inflow,outflow,symmetry,
                               slipwall,noslipwall,inflowoutflow};
}

/**
 * \namespace amrex::BCType
 * \brief Mathematical boundary condition types applied per face.
 *
 * These codes control how ghost cells are filled at domain boundaries.
 * See the conversion table at the top of this file for mappings from
 * physical to mathematical BCs.
 * \ingroup amrex_utilities
 */
namespace BCType {
enum mathematicalBndryTypes : int {
    bogus               = -666,  /**< Sentinel value for uninitialized BCs. */
    reflect_odd         = -1,    /**< F(-n) = -F(n); odd reflection. */
    int_dir             =  0,    /**< Interior / Dirichlet from interpolation. */
    reflect_even        =  1,    /**< F(-n) = F(n); even reflection. */
    foextrap            =  2,    /**< First-order extrapolation from last interior cell. */
    ext_dir             =  3,    /**< External Dirichlet value on face. */
    hoextrap            =  4,    /**< Higher-order extrapolation to face. */
    hoextrapcc          =  5,    /**< Linear extrapolation to cell center. */
    ext_dir_cc          =  6,    /**< External Dirichlet value at cell center. */
    direction_dependent =  7,    /**< Direction-dependent boundary treatment. */
    user_1              = 1001,  /**< First user-defined boundary type. */
    user_2              = 1002,  /**< Second user-defined boundary type. */
    user_3              = 1003   /**< Third user-defined boundary type. */
};
}

}
#endif

#endif
