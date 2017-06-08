Strategy for the vofStructuresA.cpp implementation:

typedef std::array<int,BL_SPACEDIM>  intDIM;
struct Node
{
    struct CutCellGraphNode
    {
        int Nnbr[BL_SPACEDIM][2];         // [DIM][lo/hi]
        int nbr[BL_SPACEDIM][2][NmvMAX];
        int ebID;
    };
    int nCells;
    intDIM iv;
    CutCellGraphNode cells[NmvMAX];
};

With NmvMAX=4, Node contains 87 integers (29 for n=2) for 2D.
Indexing on grid ID e.g., one has the following over a BoxArray:

    std::map<int,std::vector<Node> > graphNodes;

in addition, one can make an iMultiFab mask that at each IntVect
that contains the value -2 (for regular fluid cells), -1 (for covered
cells), or the index, idx, of the Node in the graphNodes[gridIndex] vector
for cut cells.

If one can pass graphNodes[index] to Fortran as gN, then one can know
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

If this turned out to be too big, we could further squash it by
managing it as a single list of integers, but there would be a lot
more arithmetic to move around.





