#ifndef _HypreABecLap_H_
#define _HypreABecLap_H_

#include <iomanip>
#include <iostream>

#include <AMReX_Tuple.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BndryData.H>
#include <AMReX_ParmParse.H>

#include <HYPRE_struct_ls.h>
#include "HYPRE_sstruct_ls.h"

// An AuxVar is an intermediate variable used in building up a
// linear system, defined as a linear combination of other
// variables which are either true independent degrees of freedom
// or other AuxVars.

class AuxVar {

 protected:

  class Connex {
   public:
    Connex() {
      other = NULL;
    }
    Connex(AuxVar* p, Real r) {
      val   = r;
      other = p;
    }
    Connex(int lev, const IntVect& iv, Real r) {
      val   = r;
      index = iv;
      level = lev;
      other = NULL;
    }
    bool same_target(const Connex& c) {
      return ((other != NULL)
              ? (other == c.other)
              : (c.other == NULL && level == c.level && index == c.index));
    }

    Real val;
    IntVect index;
    int level; //, grid; // may not need grid
    AuxVar *other;
  };

 public:

  AuxVar() : slave_flag(0) {
  }

  void push(AuxVar* p, Real r) {
    a.push_back(Connex(p,r));
  }

  void push(int lev, const IntVect& iv, Real r) {
    a.push_back(Connex(lev,iv,r));
  }

  void push_slave(AuxVar* p) {
    if (p->slave_flag == 0) { // don't count the same slave twice
      a.push_back(Connex(p,1.0));
      p->slave_flag = 1;
    }
  }

  bool empty() {
    return a.empty();
  }

  bool slave() {
    return slave_flag;
  }

  void collapse();

  int get_locations(Vector<int>& levels, Vector<IntVect>& cells);

  int get_coeffs(Vector<Real>& values);

  void clear();

 protected:

  std::list<Connex> a;
  int slave_flag;
};

// AuxVarBox is like a FAB, except that an AuxVar is stored at every
// grid point rather than a double.  It is a minimal container class
// with much less functionality than even a BaseFab.

// Could a BaseFab<AuxVar> be used instead? (AuxVar does not have arithmetic)

class AuxVarBox {

 public:

  AuxVarBox(const Box& bx) : domain(bx) {
    int numpts = domain.numPts();
    dptr = new AuxVar[numpts];
  }

  ~AuxVarBox() {
    delete[] dptr;
  }

  AuxVar& operator()(const IntVect& p) {
    BL_ASSERT(!(dptr == 0));
    BL_ASSERT(domain.contains(p));
    return dptr[domain.index(p)];
  }

  const Box& box() const {
    return domain;
  }

 protected:
  Box domain;
  AuxVar *dptr;
};

class BndryAuxVarBase {

 public:

  // Location type explanation ---
  //  INTERIOR: cells along insides of grid faces, overlaps marked "slave".
  //  EXTERIOR: cells along outsides of grid faces, overlaps (from other
  //            grids this time) marked "slave".
  //  GHOST: cells along outsides of grid faces, overlaps not marked.
  //  (Note that GHOST type is located at cell centers, but the overlap
  //  behavior more logically resembles that of face-centered data.  It
  //  is also analogous to ghost cells, in that ghost cells from
  //  different grids are distinct memory locations even if they have
  //  the same cell index.)

  enum Location { INTERIOR, EXTERIOR, GHOST };

  // This base class provides parallel iterator functionality,
  // since I'm too lazy to write true iterator classes:

  int firstLocal();
  int nextLocal(int i);
  bool isValid(int i) {
    return (i < distributionMap.ProcessorMap().size());
  }

 protected:

  DistributionMapping        distributionMap;
};

// BndryAuxVar is like a BndryRegister with AuxVars in it.  This too
// is a rather minimal container class.  It provides its own
// pseudo-iterator instead of building on something like FabSet.

class BndryAuxVar : public BndryAuxVarBase {

 public:

  BndryAuxVar(const BoxArray& _grids, Location loc);

  //~BndryAuxVar(); destructor not needed, PArrayManage is used

  PArray<AuxVarBox>& operator()(const Orientation ori) {
    return aux[ori];
  }

 protected:

  BoxArray            grids;
  PArray<AuxVarBox>   aux[2*BL_SPACEDIM];
};

// CrseBndryAuxVar is like a BndryAuxVar organized by coarse grid
// rather than by fine grid.  Entries only exist where the coarse grid
// is owned by the local processor.  Each entry has subarrays
// containing information for all fine grid edges falling within
// the local coarse grid.  Note that while cgrids defines the
// coarse grid locations and fgrids defines the fine grid locations,
// both of these BoxArrays must be at the same level of refinement!
// This class does not know the refinement ratio or problem domain.

// Mask arrays are also built, analogous to those in BndryData.

class CrseBndryAuxVar : public BndryAuxVarBase {

 public:

  CrseBndryAuxVar(const BoxArray& _cgrids,
                  const BoxArray& _fgrids, Location loc);

  CrseBndryAuxVar(const CrseBndryAuxVar& other, Location loc);

  CrseBndryAuxVar(const BoxArray& _cgrids,
                  const BoxArray& _fgrids,
                  const CrseBndryAuxVar& other, Location loc);

  void reinitialize_connections(Location loc); // use for new discretization

  void buildFaceData(IntVect& rat, int ncomp = 1);

  void rebuildFaceData(IntVect& rat, int ncomp = 1); // use when ncomp changes

  void loadFaceData(const Orientation ori,
                    MultiFab& src,
                    int srccomp,
                    int destcomp,
                    int numcomp = 1);

  //~CrseBndryAuxVar(); destructor not needed, PArrayManage is used

  Vector< PArray<AuxVarBox> >& operator()(const Orientation ori) {
    return aux[ori];
  }

  Vector< PArray<Mask> >& mask(const Orientation ori) {
    return msk[ori];
  }

  Vector< PArray<FArrayBox> >& faceData(const Orientation ori) {
    return face_data[ori];
  }

 protected:

  void initialize_slaves(Location loc);

  BoxArray                    cgrids, fgrids;
  Vector< PArray<AuxVarBox> >  aux[2*BL_SPACEDIM];
  Vector< PArray<Mask> >       msk[2*BL_SPACEDIM];
  Vector< Vector<int> >         fine_index[2*BL_SPACEDIM];
  Vector< PArray<FArrayBox> >  face_data[2*BL_SPACEDIM];
};


class HypreABecLap {

public:

  HypreABecLap(int _crse_level, int _fine_level, int _solver_flag);
  ~HypreABecLap();

  void addLevel(int level,
  		const Geometry& _geom,
  		const BoxArray& _grids,
  		IntVect         _crse_ratio);
  void setBndry(int level, BndryData& _bd);
  void buildMatrixStructure();

  void setScalars(Real sa, Real sb);
  void setACoeffs(int level, const MultiFab& alpha);
  void setBCoeffs(int level, const MultiFab& beta, int dir);
  void setRhs(int level, const MultiFab& rhs);
  void setInitGuess(int level, const MultiFab& guess);
  void solve(PArray<MultiFab>& soln, Real rel_tol, Real abs_tol, int max_iter);

protected:

  void loadMatrix();
  void setupSolver();
  void doIt();
  void getSolution(PArray<MultiFab>& soln);
  void clearSolver();

  int crse_level, fine_level, solver_flag;
  Vector<Geometry> geom;
  Vector<BoxArray> grids;
  Vector<IntVect> fine_ratio;
  PArray<BndryData> bd;

  int verbose;

  int is_periodic[BL_SPACEDIM];

  PArray<MultiFab> acoefs;
  PArray< Tuple<MultiFab, BL_SPACEDIM> > bcoefs;
  Real scalar_a, scalar_b;

  int maxiter;
  Real reltol, abstol;

  PArray<BndryAuxVar>     cintrp, ederiv;
  PArray<CrseBndryAuxVar> c_cintrp, c_ederiv, c_entry;

  HYPRE_SStructGrid hypreGrid;
  HYPRE_SStructStencil stencil;
  HYPRE_SStructGraph graph;
  HYPRE_SStructMatrix A;
  HYPRE_SStructVector b, x;
  HYPRE_Solver solver, precond;
  int ObjectType;

  int pfmg_rap_type, pfmg_relax_type, skip_relax;
  int num_pre_relax, num_post_relax;
  int print_level;

  bool x_loaded, b_loaded;

  int bho;

  // static utility functions follow:

  static int ispow2(int i)
  {
    return (i == 1) ? 1 : (((i <= 0) || (i & 1)) ? 0 : ispow2(i / 2));
  }
  
#if (BL_SPACEDIM == 1)
  static int vl[2];
  static int vh[2];
#endif

// The implicit assumption in the following is that vl[1], vh[1] are
// never changed from 0:

  static int* getV1(const IntVect& i) {
#if (BL_SPACEDIM == 1)
    vl[0] = i[0];
    return vl;
#else
    return (int*) i.getVect();
#endif
  }

  static int* getV2(const IntVect& i) {
#if (BL_SPACEDIM == 1)
    vh[0] = i[0];
    return vh;
#else
    return (int*) i.getVect();
#endif
  }

  static int* loV(const Box& b) {
#if (BL_SPACEDIM == 1)
    vl[0] = b.smallEnd(0);
    return vl;
#else
    return (int*) b.loVect();
#endif
  }
  
  static int* hiV(const Box& b) {
#if (BL_SPACEDIM == 1)
    vh[0] = b.bigEnd(0);
    return vh;
#else
    return (int*) b.hiVect();
#endif
  }

};
#endif

