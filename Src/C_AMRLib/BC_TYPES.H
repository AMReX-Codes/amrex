#ifndef _BC_TYPES_H_
#define _BC_TYPES_H_ 

#if 0
---------------------------------------------------------------------
          CONVERSIONS FROM PHYSICAL TO MATHEMATICAL BCs


PHYS BC     |   Velocity      |    Temperature   |    Scalars     |
---------------------------------------------------------------------
INTERIOR    |                 |                  |                |
   or       |    INT_DIR      |     INT_DIR      |    INT_DIR     |
PERIODIC    |                 |                  |                |
---------------------------------------------------------------------
            |                 |                  |                |
INFLOW      |    EXT_DIR      |     EXT_DIR      |    EXT_DIR     |
            |                 |                  |                |
---------------------------------------------------------------------
            |                 |                  |                |
OUTFLOW     |    FOEXTRAP       |     FOEXTRAP       |    FOEXTRAP      |
            |                 |                  |                |
---------------------------------------------------------------------
NO SLIP WALL   |                 |                  |                |
  with         |    EXT_DIR      |  REFLECT_EVEN    |    HOEXTRAP    |
ADIABATIC TEMP |    U = V = 0    |    dT/dn = 0     |                |
---------------------------------------------------------------------
NO SLIP WALL |                 |                  |                |
  with       |    EXT_DIR      |     EXT_DIR      |    HOEXTRAP    |
FIXED TEMP   |    U = V = 0    |                  |                |
---------------------------------------------------------------------
SLIP WALL   |                 |                  |                |
  with      |  Un=0 EXT_DIR   |  REFLECT_EVEN    |    HOEXTRAP    |
ADIABATIC TEMP |  Ut   HOEXTRAP  |    dT/dn = 0     |                |
---------------------------------------------------------------------
SLIP WALL   |                 |                  |                |
  with      |  Un=0 EXT_DIR   |     EXT_DIR      |    HOEXTRAP    |
FIXED TEMP  |  Ut   HOEXTRAP  |                  |                |
---------------------------------------------------------------------
            |                 |                  |                |
SYMMETRY    |  Un REFLECT_ODD |  REFLECT_EVEN    | REFLECT_EVEN   |
            |  Ut REFLECT_EVEN|                  |                |
---------------------------------------------------------------------


INT_DIR      : data taken from other grids or interpolated
EXT_DIR      : data specified on EDGE (FACE) of bndry
HOEXTRAP     : higher order extrapolation to EDGE of bndry
FOEXTRAP       : first order extrapolation from last cell in interior
REFLECT_EVEN : F(-n) = F(n)   true reflection from interior cells
REFLECT_ODD  : F(-n) = -F(n)  true reflection from interior cells

#endif

#ifndef BL_LANG_FORT
enum physicalBndryTypes {Interior,Inflow,Outflow,Symmetry,
                         SlipWall,NoSlipWall};
#endif

#define BOGUS_BC -666

#define REFLECT_ODD -1
#define INT_DIR 0
#define REFLECT_EVEN 1
#define FOEXTRAP  2
#define EXT_DIR 3
#define HOEXTRAP 4

#define Interior 0
#define Inflow 1
#define Outflow 2
#define Symmetry 3
#define SlipWall 4
#define NoSlipWall 5

#endif

