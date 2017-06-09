#ifndef _Node_H
#define _Node_H

#define COVERED_CELL  -1
#define REGULAR_CELL  -2

#ifndef BL_LANG_FORT
//
// Structure to contain all necessary graph and geometry data for cut cells at an IntVect
//

typedef std::array<int,BL_SPACEDIM>  intDIM;

struct Node
{
    struct CutCell
    {
        int Nnbr[BL_SPACEDIM][2];
        int nbr[BL_SPACEDIM][2][NCELLMAX];
        int ebID;
    };
    int nCells;
    intDIM iv;
    CutCell cells[NCELLMAX];
};

#endif

#endif
