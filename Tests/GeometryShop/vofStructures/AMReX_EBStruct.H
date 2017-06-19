#ifndef _Node_H
#define _Node_H

#define COVERED_CELL  -1
#define REGULAR_CELL  -2
#define REGULAR_FACE  -2

#ifndef BL_LANG_FORT
//
// Structure to contain all necessary graph and geometry data for cut cells at an IntVect
//

typedef std::array<int,BL_SPACEDIM>  intDIM;

struct CutCell
{
    int Nnbr[BL_SPACEDIM][2];
    int nbr[BL_SPACEDIM][2][NCELLMAX];
    int faceID[BL_SPACEDIM][2][NCELLMAX];
    int parent;
};

struct CutFace
{
    int cellHi;
    int cellLo;
    int parent;
};

#endif

#endif
