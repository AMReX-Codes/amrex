#ifndef AMREX_BC_TYPES_H_
#define AMREX_BC_TYPES_H_
#include <AMReX_Config.H>

#if 0
---------------------------------------------------------------------
          CONVERSIONS FROM PHYSICAL TO MATHEMATICAL BCs


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
  with         |    EXT_DIR      |  REFLECT_EVEN    |    HOEXTRAP    |
ADIABATIC TEMP |    U = V = 0    |    dT/dn = 0     |                |
----------------------------------------------------------------------
NO SLIP WALL   |                 |                  |                |
  with         |    EXT_DIR      |     EXT_DIR      |    HOEXTRAP    |
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

namespace PhysBCType {
enum physicalBndryTypes : int {interior=0,inflow,outflow,symmetry,
                               slipwall,noslipwall};
}

namespace BCType {
enum mathematicalBndryTypes : int {
    bogus        = -666,
    reflect_odd  = -1,
    int_dir      =  0,
    reflect_even =  1,
    foextrap     =  2,
    ext_dir      =  3,
    hoextrap     =  4,
    hoextrapcc   =  5,
    ext_dir_cc   =  6,
    user_1       = 1001,
    user_2       = 1002,
    user_3       = 1003
};
}

}
#endif

#endif
