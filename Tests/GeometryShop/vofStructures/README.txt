Strategy for the vofStructuresA.cpp implementation:

typedef std::array<int,BL_SPACEDIM>  intDIM;
struct CNode
{
    struct CutCell
    {
        int Nnbr[BL_SPACEDIM][2];
        int nbr[BL_SPACEDIM][2][NCELLMAX];
        int faceID[BL_SPACEDIM][2][NCELLMAX];
        int ebCellID;
    };
    int nCells;
    intDIM iv;
    CutCell cells[NCELLMAX];
};


struct FNode
{
    struct CutFace
    {
        int cellHi;
        int cellLo;
        int ebFaceID;
    };
    int nFaces;
    intDIM iv;
    CutFace faces[NFACEMAX];
};

With NmvMAX=4, Node contains 224 integers for 3D.
Indexing on grid ID e.g., one has the following over a BoxArray:

    std::map<int,std::vector<Node> > graphCNodes;

in addition, one can make an iMultiFab mask that at each IntVect
that contains the value -2 (for regular fluid cells), -1 (for covered
cells), or the index, idx, of the Node in the graphNodes[gridIndex] vector
for cut cells.

If one can pass graphCNodes[index] to Fortran as gN, then one can know
all about neighbors at i,j,k.  E.g. for the right neighbor:

    m = mask(i,j,k)
    if (m==-1) ... cell is "covered" by the body
    if (m==-2) ... cell is "regular fluid"
    if (m>=0) then
       nXHI = gN(m)%Nnbr[0][0]
       ...if nXHI=0, no connected XHI neighbors
       do L=1,nXHI                    ! loop over right neighbors
          nbr_XHI = gN(m)%nbr[0][L]
          if (nbr_XHI >= 0) then
             ...Lth XHI neighbor is described by gN(nbr_XHI)
          else if (nbr_XHI==-2)
             ...then neighbor is regular, nXHI=1 and mask(i+1,j,k)=-2
          endif
       enddo
    endif

Note that often, this much info is not needed.  If all cut cells are
single-valued, gN(m)%Nnbr[X][Y] = 1, all state and geometrical data
can be stored in regular fabs.

Also, if m>=0, then gN(m).cells[L].ebID could be set to point to
the index in the array of Reals holding the multi-valued cut cell
values.  So, this gives a quasi-sparse-data based algorithm to store
cell-centered data effectively as u(i,j,k,L).

If this turns out to be too big, we could further squash it by
managing it as a single list of integers, but there would be a lot
more arithmetic to move around.

FACES: 
Added faceID to struct, where faceID is the box-wide index of the face
between this cut cell and the neighbor, with the same indexing as 
nbr.  The FNode structure includes cellLo and cellHI, the local indices
of the cells on either side of the face, and CutFace's faceID.

The thought is that ebFaceID will be the index of where in the data
containers to find data associated with that face, and ebCellID the 
corresponding location of data for that cell.  The hope is that we can
construct a cell-based divergence of face-based fluxes.



