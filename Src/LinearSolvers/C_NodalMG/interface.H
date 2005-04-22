
#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#include "amr_defs.H"

class amr_boundary;

class level_interface
{
public:

#if (BL_SPACEDIM == 2)
    enum
    {
	ALL = 0xF, LOW = 0x01
    };
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

    level_interface() ;
    ~level_interface ();
    //
    // A pseudo constructor.
    //
    void alloc (const BoxArray&           Im,
		const Box&                Domain,
		const amr_boundary* Bdy);
    //
    // Two pseudo copy operators/constructors.
    //
    void copy (const level_interface& src);
    void alloc_coarsened (const BoxArray&           Im,
			  const amr_boundary* Bdy,
			  const level_interface&    src,
			  const IntVect&            rat);

    bool ok () const;

    const Box& domain () const;
    const BoxArray& interior_mesh () const;
    const BoxArray& exterior_mesh () const;
    int exterior_ref (int igrid) const;
    int direct_exterior_ref (int igrid) const;
    const Box& part_fine (int i) const;

    Array<int> geo_array (int idim, int i) const;

    int grid (int idim, int i, int j)  const;

    int fdim (int i) const;

    int nboxes (int idim) const;
    const Box& box (int idim, int i) const;
    const Box& node_box (int idim, int i) const;
    unsigned int geo (int idim, int i) const;
    int aux (int idim, int i) const;
    bool flag (int idim, int i) const;
    static int ngrids (int idim);

    mutable std::vector<int> m_fill_interface_1;
    mutable std::vector<int> m_fill_interface_2;
    mutable std::vector<int> m_fill_interface_3;
    mutable std::vector<int> m_fill_internal_borders_fn;
    mutable std::vector<int> m_fill_internal_borders_fc;

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

    void  add (std::list<Box>& bl, Box b, int startgrid);
    void xfer (const std::list<Box>& bl, int idim);
private:
    //
    // Not defined.
    //
    level_interface (const level_interface&);
    void operator= (const level_interface&);
};

extern const level_interface default_level_interface;

inline
level_interface::level_interface()
    :
    status(-1)
{}

inline
bool
level_interface::ok () const
{
    return (status >= 0);
}

inline
const Box&
level_interface::domain () const
{
    return dom;
}

inline
const BoxArray&
level_interface::interior_mesh () const
{
    return im;
}

inline
const BoxArray&
level_interface::exterior_mesh () const
{
    return em;
}

inline
int
level_interface::exterior_ref (int igrid) const
{
    BL_ASSERT(igrid < -1);
    return grid_ref[-2-igrid];
}

inline
int
level_interface::direct_exterior_ref (int igrid) const
{
    BL_ASSERT(igrid > -1); return grid_ref[   igrid];
}

inline
const Box&
level_interface::part_fine (int i) const
{
    return pf[i];
}

inline
int
level_interface::nboxes (int idim) const
{
    return nbx[idim];
}

inline
int
level_interface::grid (int idim, int i, int j) const
{
    BL_ASSERT(idim >= 0 && idim <= FACEDIM);
    BL_ASSERT(j >= 0 && i >= 0);
    int ret;
    switch(idim)
    {
    case FACEDIM:
	BL_ASSERT(j < N_FACE_GRIDS && i < nboxes(FACEDIM));
	ret = fgr[i][j];
	break;
#if BL_SPACEDIM == 3
    case 1      :
	BL_ASSERT(j < N_EDGE_GRIDS && i < nboxes(1));
	ret = egr[i][j];
	break;
#endif
    case 0      :
	BL_ASSERT(j < N_CORNER_GRIDS && i < nboxes(0));
	ret = cgr[i][j];
	break;
    }
    return ret;
}

inline
int
level_interface::ngrids(int idim)
{
    BL_ASSERT(idim >= 0 && idim <= FACEDIM);
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

inline
int
level_interface::fdim (int i) const
{
    return fdm[i];
}

inline
const Box&
level_interface::box (int idim, int i) const
{
    return bx[idim][i];
}

inline
const Box&
level_interface::node_box (int idim, int i) const
{
    return nodebx[idim][i];
}

inline
unsigned int
level_interface::geo (int idim, int i) const
{
    return ge[idim][i];
}

inline
int
level_interface::aux (int idim, int i) const
{
    return ax[idim][i];
}

inline
bool
level_interface::flag (int idim, int i) const
{
    return flg[idim][i];
}

#endif /*_INTERFACE_H_*/
