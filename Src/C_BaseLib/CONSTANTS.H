#ifndef BL_CONSTANTS_H
#define BL_CONSTANTS_H

/*
 * $Id: CONSTANTS.H,v 1.10 2001-07-17 23:40:08 car Exp $
 */

#include <REAL.H>


#ifdef BL_LANG_FORT

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define bigreal 1.0e+30
#define zero    0.0
#define one     1.0
#define two     2.0
#define three   3.0
#define four    4.0
#define five    5.0
#define six     6.0
#define seven   7.0
#define eight   8.0
#define nine    9.0
#define ten     10.0
#define twelve  12.0
#define fifteen 15.0
#define sixteen 16.0
#define twenty  20.0
#define seventy 70.0
#define ninety  90.0
#define tenth   0.1
#define eighth  0.125
#define sixth   0.16666666666667
#define fifth   0.2
#define forth   0.25
#define fourth  0.25
#define third   0.33333333333333
#define half    0.5
#define two3rd  0.66666666666667
#define Pi      3.14159265358979

#else

#define bigreal 1.0d+30
#define zero    0.0d0
#define one     1.0d0
#define two     2.0d0
#define three   3.0d0
#define four    4.0d0
#define five    5.0d0
#define six     6.0d0
#define seven   7.0d0
#define eight   8.0d0
#define nine    9.0d0
#define ten     10.0d0
#define twelve  12.0d0
#define fifteen 15.0d0
#define sixteen 16.0d0
#define twenty  20.0d0
#define seventy 70.0d0
#define ninety  90.0d0
#define tenth   0.1d0
#define eighth  0.125d0
#define sixth   0.16666666666667d0
#define fifth   0.2d0
#define forth   0.25d0
#define fourth  0.25d0
#define third   0.33333333333333d0
#define half    0.5d0
#define two3rd  0.66666666666667d0
#define Pi      3.14159265358979d0

#endif /*BL_USE_DOUBLE*/

#endif /*BL_LANG_FORT*/

#endif /*BL_CONSTANTS_H*/
