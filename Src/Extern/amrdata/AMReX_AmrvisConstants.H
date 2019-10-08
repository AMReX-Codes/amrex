
// ---------------------------------------------------------------
// AmrvisConstants.H
// ---------------------------------------------------------------
#ifndef AMREX_AMRVISCONSTANTS_H
#define AMREX_AMRVISCONSTANTS_H

#include <AMReX_REAL.H>
#include <cfloat>

#include <cstdio>

namespace amrex {

namespace Amrvis {

const int BUFSIZE = 8192;

const int LINELENGTH = 160;
const int MAXARGS = 50;

#if (BL_SPACEDIM == 1)
const int NPLANES = 1;
#endif
#if (BL_SPACEDIM == 2)
const int NPLANES = 1;
#endif
#if (BL_SPACEDIM == 3)
const int NPLANES = 3;
#endif

const int XDIR = 0;
const int YDIR = 1;
const int ZDIR = 2;

const int XY = 0;
const int XZ = 1;
const int YZ = 2;

const int XYPLANE = 0;
const int XZPLANE = 1;
const int YZPLANE = 2;

const int ZPLANE = 0;
const int YPLANE = 1;
const int XPLANE = 2;
const int ISOVIEW = 3;
const int PALETTEVIEW = 4;


const int WOFFSET = 4;

enum AnimDirection { ANIMPOSDIR = 0, ANIMNEGDIR };

// these are for PltApps wRangeRadioButton
enum MinMaxRangeTypeForButtons {
    BGLOBALMINMAX        = 0,
    BSUBREGIONMINMAX     = 1,
    BUSERMINMAX          = 2,
    BNUMBEROFMINMAX      = 3,
    BINVALIDMINMAX
};

enum MinMaxRangeType {
    GLOBALMINMAX        = 0,
    SUBREGIONMINMAX     = 1,
    USERMINMAX          = 2,
    FILEGLOBALMINMAX    = 3,
    FILESUBREGIONMINMAX = 4,
    FILEUSERMINMAX      = 5,
    NUMBEROFMINMAX      = 6,
    INVALIDMINMAX
};

enum ContourType {
  RASTERONLY     = 0,
  RASTERCONTOURS = 1,
  COLORCONTOURS  = 2,
  BWCONTOURS     = 3,
  VECTORS        = 4,
  NCONTOPTIONS   = 5,
  INVALIDCONTOURTYPE
};

enum FileType { INVALIDTYPE = 0, FAB, MULTIFAB, NEWPLT, PROFDATA, NUMDBPTYPES };

}  // end namespace Amrvis

}

#endif
// -------------------------------------------------------------------
// -------------------------------------------------------------------
