#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#include "amr_defs.H"

class amr_boundary_class;
//typedef const amr_boundary_class& amr_boundary;

extern const class level_interface null_level_interface;

class level_interface 
{
 public:

#if (BL_SPACEDIM == 2)
     enum { LL = 0x1, HL = 0x2, LH = 0x4, HH = 0x8,
	    XL = 0x5, XH = 0xA, YL = 0x3, YH = 0xC,
	    ALL = 0xF, LOW = 0x1 };
#elif (BL_SPACEDIM == 3)
     enum { LLL = 0x01, HLL = 0x02, LHL = 0x04, HHL = 0x08,
            LLH = 0x10, HLH = 0x20, LHH = 0x40, HHH = 0x80,
	    XL  = 0x55, XH  = 0xAA, YL  = 0x33, YH  = 0xCC,
            ZL  = 0x0F, ZH  = 0xF0, ALL = 0xFF, LOW = 0x01 };
#endif
     enum { FACEDIM  = BL_SPACEDIM - 1} ;
#if (BL_SPACEDIM == 2)
     enum { N_FACE_GRIDS = 2,
            N_EDGE_GRIDS = 2, 
            N_CORNER_GRIDS = 4 };
#else
     enum { N_FACE_GRIDS = 2,
            N_EDGE_GRIDS = 4,
            N_CORNER_GRIDS = 8 };
#endif

  level_interface() : status(-1) { }
  ~level_interface();

  void alloc(const BoxArray& Im, const Box& Domain, const amr_boundary_class& Bdy);
  level_interface(const BoxArray& Im, const Box& Domain, const amr_boundary_class& Bdy)
    : status(-1) {
      alloc(Im, Domain, Bdy);
  }
  void copy(const level_interface& src);
  void alloc_coarsened(const BoxArray& Im, const amr_boundary_class& Bdy,
		       const level_interface& src, const IntVect& rat);

  int ok()                          const { return (status >= 0); }
  int null()                        const { return (status <  0); }

  const Box& domain()               const { return dom; }
  const BoxArray& interior_mesh()   const { return im; }
  const BoxArray& exterior_mesh()   const { return em; }
  int exterior_ref(int igrid)       const { return grid_ref[-2-igrid]; }
  int direct_exterior_ref(int igrid)const { return grid_ref[igrid]; }
  const Box& part_fine(int i)       const { return pf[i]; }

  int nboxes(int idim)                 const { return nbx[idim]; }
  const Box& box(int idim, int i)      const { return bx[idim][i]; }
  const Box& node_box(int idim, int i) const { return nodebx[idim][i]; }
  unsigned geo(int idim, int i)        const { return ge[idim][i]; }
  unsigned aux(int idim, int i)        const { return ax[idim][i]; }
  unsigned flag(int idim, int i)       const { return flg[idim][i]; }

  void geo_array(int ga[], int idim, int i) const;

  int nfaces()                      const { return nbx[FACEDIM]; }
  const Box& face(int i)            const { return bx[FACEDIM][i]; }
  const Box& node_face(int i)       const { return nodebx[FACEDIM][i]; }
  unsigned fgeo(int i)              const { return ge[FACEDIM][i]; }
  int faux(int i)                   const { return ax[FACEDIM][i]; }
  unsigned fflag(int i)             const { return flg[FACEDIM][i]; }
  int fgrid(int i, int j)           const { return fgr[i][j]; }
  int fdim(int i)                   const { return fdm[i]; }

  int nedges()                      const { return nbx[1]; }
  const Box& edge(int i)            const { return bx[1][i]; }
  const Box& node_edge(int i)       const { return nodebx[1][i]; }
  unsigned egeo(int i)              const { return ge[1][i]; }
  int eaux(int i)                   const { return ax[1][i]; }
  unsigned eflag(int i)             const { return flg[1][i]; }
  int egrid(int i, int j)           const { return egr[i][j]; }

  int ncorners()                    const { return nbx[0]; }
  const Box& corner(int i)          const { return bx[0][i]; }
  unsigned cgeo(int i)              const { return ge[0][i]; }
  int caux(int i)                   const { return ax[0][i]; }
  unsigned cflag(int i)             const { return flg[0][i]; }
  int cgrid(int i, int j)           const { return cgr[i][j]; }

 protected:
  int status;                       // tells if some arrays are shared
  Box dom;                          // BoxArray doesn't know domain
  BoxArray im, em;                  // interior and exterior BoxArray
  int *grid_ref;                    // one for each exterior box
  Box *pf;                          // expanded NODE boxes for interior grids
  int nbx[BL_SPACEDIM];             // numbers of interface boxes
  Box *bx[BL_SPACEDIM];             // array of interface boxes
  unsigned *ge[BL_SPACEDIM];        // geometries of interface boxes
  int *ax[BL_SPACEDIM];             // if >= 0, gives the number of the
    // grid which must be cleared of bogus data along this edge, corner,
    // or face after computing on the part_fine section of that grid
  unsigned *flg[BL_SPACEDIM];       // flags dirichlet bdys
  int *fdm;                         // which is flat dim for this face
  Box *nodebx[BL_SPACEDIM];         // node-based versions of interface boxes
  int (*fgr)[2];                    // grid neighbor data for faces
  int (*egr)[N_EDGE_GRIDS];         // grid neighbor data for edges
  int (*cgr)[N_CORNER_GRIDS];       // grid neighbor data for corners

  void  add(List<Box>& bl, Box b, int startgrid = 0);
  void xfer(List<Box>& bl, int idim);

 public:
  static void ins(List<Box>& bl, const Box& b) {
    if (!bl.includes(b))
      bl.append(b);
  }
};

#endif
