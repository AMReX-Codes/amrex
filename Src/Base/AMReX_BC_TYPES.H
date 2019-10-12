#ifndef AMREX_BC_TYPES_H_
#define AMREX_BC_TYPES_H_

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
EXT_DIR      : data specified on EDGE (FACE) of bndry for the linear solvers
               for application codes, the location of the EXT_DIR data will
               depend on how they do reconstruction, and may be edge or
               cell-centered.
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
    hoextrapcc   =  5
};
}

}
#endif

#ifndef AMREX_XSDK

#define BOGUS_BC -666

#define REFLECT_ODD -1
#define INT_DIR 0
#define REFLECT_EVEN 1
#define FOEXTRAP  2
#define EXT_DIR 3
#define HOEXTRAP 4
#define HOEXTRAPCC 5

#define Interior 0
#define Inflow 1
#define Outflow 2
#define Symmetry 3
#define SlipWall 4
#define NoSlipWall 5

#endif

#endif

