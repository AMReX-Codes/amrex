#ifndef _AMR_GRAPH_H_
#define _AMR_GRAPH_H_

#ifndef _AMR_DEFS_H_
#  include "amr_defs.H"
#endif
#ifndef _CBASICS_H_
#  include "cbasics.H"
#endif

extern const Box      null_Box;
extern const BoxArray null_level_mesh;

#if (SPACEDIM == 1)
void real_coords(const IntVect& iv, const IntVect& type, int ratio, Real& x);
void real_coords(const Box& bx, int ratio, Real& x0, Real& x1);
#elif (SPACEDIM == 2)
void real_coords(const IntVect& iv, const IntVect& type, int ratio,
		 Real& x, Real& y);
void real_coords(const Box& bx, int ratio,
		 Real& x0, Real& x1, Real& y0, Real& y1);
#elif (SPACEDIM == 3)
void real_coords(const IntVect& iv, const IntVect& type, int ratio,
		 Real& x, Real& y, Real& z);
void real_coords(const Box& bx, int ratio,
		 Real& x0, Real& x1, Real& y0, Real& y1, Real& z0, Real& z1);
#endif

#if (SPACEDIM == 3)
void set_graphics_knobs(Real theta, Real phi);
#endif

// Set dimensions so the given box fills the image with the given margin.
void fit(Box b, int ratio = 1, Real margin = 0.05);

// Plot the given mesh object:

void plot(Box b, int ratio = 1);

void plot(const BoxArray& mesh, int ratio = 1);

void plot(const Array<BoxArray>& mesh, const Array<IntVect>& gen_ratio);

// The contour plotting routines take a pointer to a function like
// this one for setting the color of the ith contour from the range
// [0, n).  This default function chooses from a rainbow spectrum if
// n > 1, otherwise it leaves the current color unchanged.

inline void default_contour_color(int i, int n)
  { if (n > 1) colorscale((Real)i/(n-1)); }

// In the contour functions, if cflag is true the window will
// be cleared and the mesh plotted on it first.  If dflag is also
// true the window dimensions will be set to include the entire
// domain.  Passing false values for these flags permits user
// control of the plotting process.  (In the grid_real versions
// the domain is unknown, so instead of dflag there is an optional
// Box argument that specifies the domain.)

// Plot a single contour for the specified value:

void contour(const Fab& r, Box subbox, int ratio, Real value,
	     int cflag = 1, const Box& domain = null_Box,
	     void (*set_contour_color)(int,int) = default_contour_color,
	     const BoxArray& mask = null_level_mesh, int mrat = 1);

void contour(MultiFab& r, int ratio, Real value,
	     int cflag = 1, const Box& domain = null_Box,
	     void (*set_contour_color)(int,int) = default_contour_color,
	     const BoxArray& mask = null_level_mesh, int mrat = 1);

void contour(PArray<MultiFab>& r, const Array<IntVect>& gen_ratio, Real value,
	     int cflag = 1, const Box& domain0 = null_Box,
	     void (*set_contour_color)(int,int) = default_contour_color);

// Plot nval contours at the values specified by value[0...nval-1]:

void contour(const Fab& r, Box subbox, int ratio,
	     int nval, Real value[],
	     int cflag = 1, const Box& domain = null_Box,
	     void (*set_contour_color)(int,int) = default_contour_color,
	     const BoxArray& mask = null_level_mesh, int mrat = 1);

void contour(MultiFab& r, int ratio, int nval, Real value[],
	     int cflag = 1, const Box& domain = null_Box,
	     void (*set_contour_color)(int,int) = default_contour_color,
	     const BoxArray& mask = null_level_mesh, int mrat = 1);

void contour(PArray<MultiFab>& r, const Array<IntVect>& gen_ratio,
	     int nval, Real value[],
	     int cflag = 1, const Box& domain0 = null_Box,
	     void (*set_contour_color)(int,int) = default_contour_color);

// Use nval contours assigned according to a scheme selected through flag:
// flag =  0:  Equally spaced from -a to a.
//             If a == 0.0 then max_norm*(nval-1)/nval will be used.
// flag =  1:  Same except that if nval is odd the 0 contour is suppressed.
// flag =  2:  Spaced at points a+(b-a)*pow(i/(nval-1),p), i from [0,nval).
//             If a == b then alternate points from the pattern
//             of 2*nval+1 values from rmin to rmax will be used.

void contour(const Fab& r, Box subbox, int ratio,
	     int nval, int flag = 0,
	     Real a = 0.0, Real b = 0.0, Real p = 1.0,
	     int cflag = 1, const Box& domain = null_Box,
	     void (*set_contour_color)(int,int) = default_contour_color,
	     const BoxArray& mask = null_level_mesh, int mrat = 1);

void contour(MultiFab& r, int ratio, int nval, int flag = 0,
	     Real a = 0.0, Real b = 0.0, Real p = 1.0,
	     int cflag = 1, const Box& domain = null_Box,
	     void (*set_contour_color)(int,int) = default_contour_color,
	     const BoxArray& mask = null_level_mesh, int mrat = 1);

void contour(PArray<MultiFab>& r, const Array<IntVect>& gen_ratio,
	     int nval, int flag = 0,
	     Real a = 0.0, Real b = 0.0, Real p = 1.0,
	     int cflag = 1, const Box& domain0 = null_Box,
	     void (*set_contour_color)(int,int) = default_contour_color);

#if (SPACEDIM == 3)

// Slice routines: plot 2D slice out of 3D grid.  Only single-grid
// version currently exists.

// Set dimensions so the given box fills the image with the given margin.
void fit_slice(Box b, int idim, int ratio = 1, Real margin = 0.05);

// Plot the given mesh object:

void plot_slice(Box b, int idim, int ratio = 1);

void slice(const grid_real& r, Box subbox, int ratio,
           int idim, int islice,
           Real value,
           int cflag = 1, const Box& domain = null_Box,
           void (*set_contour_color)(int,int) = default_contour_color);

void slice(const grid_real& r, Box subbox, int ratio,
           int idim, int islice,
           int nval, Real value[],
           int cflag = 1, const Box& domain = null_Box,
           void (*set_contour_color)(int,int) = default_contour_color);

void slice(const grid_real& r, Box subbox, int ratio,
           int idim, int islice,
           int nval, int flag = 0,
           Real a = 0.0, Real b = 0.0, Real p = 1.0,
           int cflag = 1, const Box& domain = null_Box,
           void (*set_contour_color)(int,int) = default_contour_color);

#endif

#endif
