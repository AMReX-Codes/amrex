#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#include "amr_defs.H"

class amr_boundary_class;

class level_interface 
{
public:
    
#if (BL_SPACEDIM == 2)
#if 0
    enum 
    { 
	LL = 0x1, HL = 0x2, LH = 0x4, HH = 0x8,
	XL = 0x5, XH = 0xA, YL = 0x3, YH = 0xC,
	ALL = 0xF, LOW = 0x1 
    };
#else
    enum
    {
	ALL = 0xF, LOW = 0x01
    };
#endif
#else
    enum 
    {
	ALL = 0xFF, LOW = 0x01 
    };
#endif
    enum 
    {
	FACEDIM = BL_SPACEDIM-1
    };
    
    level_interface() 
	: status(-1) 
    { 
    }
    ~level_interface();
    // a pseudo constructor.
    void alloc(const BoxArray& Im, const Box& Domain, const amr_boundary_class* Bdy);
    // two pseudo copy operators/constructors.
    void copy(const level_interface& src);
    void alloc_coarsened(const BoxArray& Im, const amr_boundary_class* Bdy, const level_interface& src, const IntVect& rat);
    
    bool ok() const			    { return (status >= 0); }
    
    const Box& domain() const		    { return dom; }
    const BoxArray& interior_mesh() const   { return im; }
    const BoxArray& exterior_mesh() const   { return em; }
    int exterior_ref(int igrid) const	    { assert(igrid < -1); return grid_ref[-2-igrid]; }
    int direct_exterior_ref(int igrid)const { assert(igrid > -1); return grid_ref[   igrid]; }
    const Box& part_fine(int i) const	    { return pf[i]; }
    
    Array<int> geo_array(int idim, int i) const;

    int grid(int idim, int i, int j)  const;

    int fdim(int i) const		    { return fdm[i]; }

    int nboxes(int idim) const		    { return nbx[idim]; }
    const Box& box(int idim, int i) const   { return bx[idim][i]; }
    const Box& node_box(int idim, int i) const { return nodebx[idim][i]; }
    unsigned int geo(int idim, int i) const { return ge[idim][i]; }
    unsigned int aux(int idim, int i) const { return ax[idim][i]; }
    bool flag(int idim, int i) const	    { return flg[idim][i]; }
    static int ngrids(int idim);

private:
#if (BL_SPACEDIM == 2)
    enum
    { 
	N_FACE_GRIDS = 2,
	N_CORNER_GRIDS = 4
    };
#elif (BL_SPACEDIM == 3)
    enum 
    { 
	N_FACE_GRIDS = 2,
	N_EDGE_GRIDS = 4,
	N_CORNER_GRIDS = 8
    };
#endif
    int status;				// tells if some arrays are shared
    Box dom;				// BoxArray doesn't know domain
    BoxArray im;			// interior BoxArray
    BoxArray em;			// exterior BoxArray
    int* grid_ref;                      // one for each exterior box
    Box* pf;                            // expanded NODE boxes for interior grids
    int nbx[BL_SPACEDIM];               // numbers of interface boxes
    Box* bx[BL_SPACEDIM];               // array of interface boxes
    unsigned int* ge[BL_SPACEDIM];      // geometries of interface boxes
    int* ax[BL_SPACEDIM];               // if >= 0, gives the number of the
          // grid which must be cleared of bogus data along this edge, corner,
          // or face after computing on the part_fine section of that grid
    bool* flg[BL_SPACEDIM];             // flags dirichlet bdys
    int* fdm;                           // which is flat dim for this face
    Box* nodebx[BL_SPACEDIM];           // node-based versions of interface boxes
    int (*fgr)[N_FACE_GRIDS];           // grid neighbor data for faces
#if BL_SPACEDIM==3
    int (*egr)[N_EDGE_GRIDS];           // grid neighbor data for edges
#endif
    int (*cgr)[N_CORNER_GRIDS];         // grid neighbor data for corners
    
    void  add(List<Box>& bl, Box b, int startgrid);
    void xfer(const List<Box>& bl, int idim);
private:
    level_interface(const level_interface&);
    void operator=(const level_interface&);
};

inline int level_interface::grid(int idim, int i, int j)   const 
{ 
    assert(idim >= 0 && idim <= FACEDIM);
    assert(j >= 0 && i >= 0);
    int ret;
    switch(idim) 
    {
    case FACEDIM: assert(j < N_FACE_GRIDS && i < nboxes(FACEDIM)); ret = fgr[i][j]; break;
#if BL_SPACEDIM == 3
    case 1      : assert(j < N_EDGE_GRIDS && i < nboxes(1)); ret = egr[i][j]; break;
#endif
    case 0      : assert(j < N_CORNER_GRIDS && i < nboxes(0)); ret = cgr[i][j]; break;
    }
    return ret;
}

inline int level_interface::ngrids(int idim)
{
    assert(idim >= 0 && idim <= FACEDIM);
    int ret;
    switch(idim)
    {
    case FACEDIM: ret = N_FACE_GRIDS; break;
#if BL_SPACEDIM == 3
    case 1:       ret = N_EDGE_GRIDS; break;
#endif
    case 0:       ret = N_CORNER_GRIDS; break;
    }
    return ret;
}
#endif
