/* Copyright 2019 Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_DTTYPE_H_
#define WARPX_DTTYPE_H_

enum struct DtType : int
{
    Full = 0,
    FirstHalf,
    SecondHalf
};

#endif // WARPX_DTTYPE_H_
