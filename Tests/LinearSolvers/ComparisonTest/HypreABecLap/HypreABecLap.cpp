
#include <HypreABecLap.H>
#include <HypreABec_F.H>

#include "AMReX_LO_BCTYPES.H"

#include "_hypre_sstruct_mv.h"
#include "HYPRE_krylov.h"

using std::list;

#if (BL_SPACEDIM == 1)
int HypreMultiABec::vl[2] = { 0, 0 };
int HypreMultiABec::vh[2] = { 0, 0 };
#endif

void AuxVar::collapse()
{
  // "Flattens" the dependency list.  Any entry that points to another
  // AuxVar is replaced by a copy of the dependency list from that AuxVar.
  // The values in the inserted elements are multiplied by the value from
  // the entry being replaced.  AuxVar pointers that were in the dependency
  // list are removed from it, but the AuxVars they point to are not
  // themselves modified.  Entries pointing to the same cell are then
  // combined.

  // The loops in this function are tricky, since we're modifying
  // a list at the same time we're looping through it.  for statements
  // are used in somewhat unusual ways.  Erasing an entry invalidates
  // its iterator, so before we can do that we need to have another
  // iterator pointing to an adjacent entry.  While "past-the-end" is
  // a valid target for an iterator, we assume there is no corresponding
  // "before-the-beginning" position.

  // There is one property that can be useful in dealing with the hypre
  // semi-structured interface:  Call two AuxVar objects or two Connex
  // lists "similar" if they point to the same locations in the same order.
  // If two AuxVar objects are built whose Connex lists are similar up to
  // a point, and all the latter entries point only to locations already
  // referenced in the similar initial sections, then after being
  // collapsed the two AuxVars will be similar.  If one AuxVar is then
  // used for establishing the graph entries, and the other is used for
  // initializing the matrix coefficients, then it is ok for the second
  // AuxVar to include additional entries so long as they don't point
  // to any new locations not already established in the graph.

  for (list<Connex>::iterator it = a.begin(); it != a.end(); ) {
    AuxVar *o = it->other;
    if (o) {
      list<Connex>::iterator kt = it;
      a.insert(++kt, o->a.begin(), o->a.end());

      // multiply inserted entries by it->var here:
      for (list<Connex>::iterator jt = it; ++jt != kt; ) {
        jt->val *= it->val;
      }

      kt = it; // entry to be erased
      ++it;
      a.erase(kt);
    }
    else {
      ++it;
    }
  }

  // Combine entries that point to same cell:

  for (list<Connex>::iterator it = a.begin(); it != a.end(); ++it) {
    for (list<Connex>::iterator jt = it; ++jt != a.end(); ) {
      if (it->same_target(*jt)) {
        it->val += jt->val;
        list<Connex>::iterator kt = jt;
        --jt;
        a.erase(kt);
      }
    }
  }
}

void AuxVar::clear()
{
  a.clear();
  slave_flag = 0;
}

int AuxVar::get_locations(Vector<int>& levels, Vector<IntVect>& cells)
{
  if (slave()) {
    return 1; // failure
  }

  int n = a.size();
  levels.resize(n);
  cells.resize(n);

  int i = 0;
  for (list<Connex>::iterator it = a.begin(); it != a.end(); ++it) {
    if (it->other != NULL) {
      return 2; // failure
    }
    levels[i] = it->level;
    cells[i]  = it->index;
    i++;
  }

  return 0; // success
}

int AuxVar::get_coeffs(Vector<Real>& values)
{
  if (slave()) {
    return 1; // failure
  }

  int n = a.size();
  values.resize(n);

  int i = 0;
  for (list<Connex>::iterator it = a.begin(); it != a.end(); ++it) {
    if (it->other != NULL) {
      return 2; // failure
    }
    values[i] = it->val;
    i++;
  }

  return 0; // success
}

int BndryAuxVarBase::firstLocal()
{
  const int MyProc = ParallelDescriptor::MyProc();
  int i = 0;
  while (distributionMap[i] != MyProc) {
    i++;
  }
  return i;
}

int BndryAuxVarBase::nextLocal(int i)
{
  const int MyProc = ParallelDescriptor::MyProc();
  do {
    i++;
  } while (distributionMap[i] != MyProc);
  return i;
}

BndryAuxVar::BndryAuxVar(const BoxArray& _grids, Location loc)
  : grids(_grids)
{
  distributionMap.define(grids, ParallelDescriptor::NProcs());

  // For Location type EXTERIOR use CrseBndryAuxVar instead:
  BL_ASSERT(loc == INTERIOR || loc == GHOST);

  int inormal = (loc == INTERIOR) ? 0 : 1;

  for (OrientationIter oitr; oitr; ++oitr) {
    Orientation ori = oitr();
    aux[ori].resize(grids.size(), PArrayManage);
    int ishift = ori.isLow() ? 1-inormal : inormal-1;
    for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
      Box b = amrex::adjCell(grids[i], ori);
      aux[ori].set(i, new AuxVarBox(b.shift(ori.coordDir(), ishift)));
    }
  }

  // Make master-slave connections:

  if (loc == INTERIOR) {
    for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
#if 1
      // This is the current default implementation:
      for (OrientationIter omitr; omitr; ++omitr) {
        Orientation om = omitr();
        const Box& bm = aux[om][i].box();
        for (OrientationIter ositr; ositr; ++ositr) {
          Orientation os = ositr();
          if (os.coordDir() > om.coordDir()) {
            const Box& bs = aux[os][i].box();
            if (bm.intersects(bs)) {
              Box reg = (bm & bs);
              for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
                aux[om][i](v).push_slave(&aux[os][i](v));
              }
            }
          }
        }
      }
#elif 0
      // This was the original implementation, 2D only:
      Orientation oxlo(0, Orientation::low);
      Orientation oxhi(0, Orientation::high);
      Orientation oylo(1, Orientation::low);
      Orientation oyhi(1, Orientation::high);
      IntVect p = aux[oxlo][i].box().smallEnd();
      aux[oxlo][i](p).push_slave(&aux[oylo][i](p));
      p = aux[oxlo][i].box().bigEnd();
      aux[oxlo][i](p).push_slave(&aux[oyhi][i](p));
      p = aux[oxhi][i].box().smallEnd();
      aux[oxhi][i](p).push_slave(&aux[oylo][i](p));
      p = aux[oxhi][i].box().bigEnd();
      aux[oxhi][i](p).push_slave(&aux[oyhi][i](p));
#elif 0
      // This version is like the new default, except that
      // it loops through orientations in a different order.
      // Some master/slave pairs are therefore flipped, and in
      // the end the solvers return slightly different numbers.
      // (Results should be the same within the solver tolerance.)
      for (OrientationIter omitr; omitr; ++omitr) {
        Orientation om = omitr();
        const Box& bm = aux[om][i].box();
        for (OrientationIter ositr(om); ++ositr; ) {
          Orientation os = ositr();
          const Box& bs = aux[os][i].box();
          if (bm.intersects(bs)) {
            Box reg = (bm & bs);
            for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
              aux[om][i](v).push_slave(&aux[os][i](v));
            }
          }
        }
      }
#endif
    }
  }
}

CrseBndryAuxVar::CrseBndryAuxVar(const BoxArray& _cgrids,
                                 const BoxArray& _fgrids, Location loc)
  : cgrids(_cgrids), fgrids(_fgrids)
{
  distributionMap.define(cgrids, ParallelDescriptor::NProcs());

  // For Location type INTERIOR use BndryAuxVar instead:
  BL_ASSERT(loc == EXTERIOR || loc == GHOST);

  for (OrientationIter oitr; oitr; ++oitr) {
    Orientation ori = oitr();

    aux[ori].resize(cgrids.size());
    msk[ori].resize(cgrids.size());

    fine_index[ori].resize(cgrids.size());

    for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
      list<Box> bl;
      list<int> fi;
      for (int j = 0; j < fgrids.size(); j++) {
        Box face = amrex::adjCell(fgrids[j], ori);
        if (cgrids[i].intersects(face)) {
          bl.push_back(face & cgrids[i]);
          fi.push_back(j);
        }
      }

      // bl now has every fine grid face in it, even those entirely
      // covered by other fine grids.  We might optimize by removing
      // the covered ones here.

      int n = bl.size();
      aux[ori][i].resize(n, PArrayManage);
      msk[ori][i].resize(n, PArrayManage);

      fine_index[ori][i].resize(n);

      int j = 0;
      for (list<int>::iterator it = fi.begin(); it != fi.end(); ++it) {
        fine_index[ori][i][j++] = *it;
      }

      j = 0;
      for (list<Box>::iterator it = bl.begin(); it != bl.end(); ++it) {
        aux[ori][i].set(j, new AuxVarBox(*it));
        Box mask_box = *it;
        for (int dir = 0; dir < BL_SPACEDIM; dir++) {
          if (dir == ori.coordDir())
            continue;
          mask_box.grow(dir,1);
        }
        msk[ori][i].set(j, new Mask(mask_box));

        msk[ori][i][j].setVal(BndryData::outside_domain);

        // If we had the Geometry, we would have the problem domain
        // and would not have to loop over coarse grids below.  Also,
        // if we ever care to do this right for periodic domains we
        // will need the Geometry.  (See BndryData.cpp)

        for (int k = 0; k < cgrids.size(); k++) {
          if (cgrids[k].intersects(mask_box)) {
            msk[ori][i][j].setVal(BndryData::not_covered,
                                  (cgrids[k] & mask_box), 0);
          }
        }

        for (int k = 0; k < fgrids.size(); k++) {
          if (fgrids[k].intersects(mask_box)) {
            msk[ori][i][j].setVal(BndryData::covered,
                                  (fgrids[k] & mask_box), 0);
          }
        }

        j++;
      }
    }
  }

  initialize_slaves(loc);
}

CrseBndryAuxVar::CrseBndryAuxVar(const CrseBndryAuxVar& other, Location loc)
  : cgrids(other.cgrids), fgrids(other.fgrids)
{
  distributionMap = other.distributionMap;

  // For Location type INTERIOR use BndryAuxVar instead:
  BL_ASSERT(loc == EXTERIOR || loc == GHOST);

  for (OrientationIter oitr; oitr; ++oitr) {
    Orientation ori = oitr();

    aux[ori].resize(cgrids.size());
    msk[ori].resize(cgrids.size());

    fine_index[ori].resize(cgrids.size());

    for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
      int n = other.aux[ori][i].size();

      aux[ori][i].resize(n, PArrayManage);
      msk[ori][i].resize(n, PArrayManage);

      fine_index[ori][i].resize(n);

      for (int j = 0; j < n; j++) {
        fine_index[ori][i][j] = other.fine_index[ori][i][j];

        aux[ori][i].set(j, new AuxVarBox(other.aux[ori][i][j].box()));

        msk[ori][i].set(j, new Mask(other.msk[ori][i][j].box()));
        msk[ori][i][j].copy(other.msk[ori][i][j]);
      }
    }
  }

  initialize_slaves(loc);
}

CrseBndryAuxVar::CrseBndryAuxVar(const BoxArray& _cgrids,
                                 const BoxArray& _fgrids,
                                 const CrseBndryAuxVar& other, Location loc)
  : cgrids(_cgrids), fgrids(_fgrids)
{
  distributionMap = other.distributionMap;

  // For Location type INTERIOR use BndryAuxVar instead:
  BL_ASSERT(loc == EXTERIOR || loc == GHOST);

  for (OrientationIter oitr; oitr; ++oitr) {
    Orientation ori = oitr();

    aux[ori].resize(cgrids.size());
    msk[ori].resize(cgrids.size());

    fine_index[ori].resize(cgrids.size());

    for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
      int n = other.aux[ori][i].size();

      aux[ori][i].resize(n, PArrayManage);
      msk[ori][i].resize(n, PArrayManage);

      fine_index[ori][i].resize(n);

      for (int j = 0; j < n; j++) {
        fine_index[ori][i][j] = other.fine_index[ori][i][j];

        Box face = amrex::adjCell(fgrids[fine_index[ori][i][j]], ori);
        aux[ori][i].set(j, new AuxVarBox(face & cgrids[i]));

        Box mask_box = aux[ori][i][j].box();
        for (int dir = 0; dir < BL_SPACEDIM; dir++) {
          if (dir == ori.coordDir())
            continue;
          mask_box.grow(dir,1);
        }
        msk[ori][i].set(j, new Mask(mask_box));

        msk[ori][i][j].setVal(BndryData::outside_domain);

        // If we had the Geometry, we would have the problem domain
        // and would not have to loop over coarse grids below.  Also,
        // if we ever care to do this right for periodic domains we
        // will need the Geometry.  (See BndryData.cpp)

        for (int k = 0; k < cgrids.size(); k++) {
          if (cgrids[k].intersects(mask_box)) {
            msk[ori][i][j].setVal(BndryData::not_covered,
                                  (cgrids[k] & mask_box), 0);
          }
        }

        for (int k = 0; k < fgrids.size(); k++) {
          if (fgrids[k].intersects(mask_box)) {
            msk[ori][i][j].setVal(BndryData::covered,
                                  (fgrids[k] & mask_box), 0);
          }
        }

      }
    }
  }

  initialize_slaves(loc);
}

void CrseBndryAuxVar::reinitialize_connections(Location loc)
{
  for (OrientationIter oitr; oitr; ++oitr) {
    Orientation ori = oitr();
    for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
      int n = aux[ori][i].size();
      for (int j = 0; j < n; j++) {
        const Box& reg = aux[ori][i][j].box();
        for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
          aux[ori][i][j](v).clear();
        }
      }
    }
  }

  initialize_slaves(loc);
}

void CrseBndryAuxVar::initialize_slaves(Location loc)
{
  // Make master-slave connections:

  if (loc == EXTERIOR) {
    for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
      for (OrientationIter omitr; omitr; ++omitr) {
        Orientation om = omitr();
        for (OrientationIter ositr(om); ++ositr; ) {
          Orientation os = ositr();

          for (int jm = 0; jm < aux[om][i].size(); jm++) {
            const Box& bm = aux[om][i][jm].box();
            for (int js = 0; js < aux[os][i].size(); js++) {
              const Box& bs = aux[os][i][js].box();

              if (bm.intersects(bs)) {
                Box reg = (bm & bs);
                for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
                  aux[om][i][jm](v).push_slave(&aux[os][i][js](v));
                }
              }

            }
          }

        }
      }
    }
  }
}

void CrseBndryAuxVar::buildFaceData(IntVect& rat, int ncomp)
{
  for (OrientationIter oitr; oitr; ++oitr) {
    Orientation ori = oitr();

    face_data[ori].resize(cgrids.size());

    for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
      int n = aux[ori][i].size();
      face_data[ori][i].resize(n, PArrayManage);
      for (int j = 0; j < n; j++) {
        Box face_box = aux[ori][i][j].box();
        // face_box is known to be "adjacent cell", convert to face:
        face_box.shiftHalf(ori.coordDir(), (ori.isLow() ? 1 : -1));
        face_box.refine(rat);
        face_data[ori][i].set(j, new FArrayBox(face_box, ncomp));
      }
    }
  }
}

void CrseBndryAuxVar::rebuildFaceData(IntVect& rat, int ncomp)
{
  for (OrientationIter oitr; oitr; ++oitr) {
    Orientation ori = oitr();

    for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
      int n = face_data[ori][i].size();
      for (int j = 0; j < n; j++) {
        const Box& face_box = face_data[ori][i][j].box();
        face_data[ori][i][j].resize(face_box, ncomp);
      }
    }
  }
}

void CrseBndryAuxVar::loadFaceData(const Orientation ori,
                                   MultiFab& src,
                                   int srccomp,
                                   int destcomp,
                                   int numcomp)
{
  MultiFabCopyDescriptor mfcd;

  MultiFabId mfid = mfcd.RegisterMultiFab(&src);

  Vector< Vector<FillBoxId> > fbid(face_data[ori].size());
  for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
    int n = face_data[ori][i].size();
    fbid[i].resize(n);
    for (int j = 0; j < n; j++) {
      fbid[i][j] = mfcd.AddBox(mfid, face_data[ori][i][j].box(), NULL,
                               fine_index[ori][i][j],
                               srccomp, destcomp, numcomp);
    }
  }

  mfcd.CollectData();

  for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
    int n = face_data[ori][i].size();
    for (int j = 0; j < n; j++) {
      mfcd.FillFab(mfid, fbid[i][j], face_data[ori][i][j]);
    }
  }
}


HypreABecLap::HypreABecLap(int _crse_level, int _fine_level,
                               int _solver_flag)
  : crse_level(_crse_level), 
    fine_level(_fine_level),
    solver_flag(_solver_flag),
    geom(fine_level+1),
    grids(fine_level+1),
    fine_ratio(fine_level+1),
    bd(fine_level+1),
    acoefs(fine_level+1),
    bcoefs(fine_level+1),
    cintrp(fine_level+1),
    ederiv(fine_level+1),
    c_cintrp(fine_level+1),
    c_ederiv(fine_level+1),
    c_entry(fine_level+1),
    hypreGrid(NULL), stencil(NULL), graph(NULL), 
    A(NULL), b(NULL), x(NULL),
    solver(NULL), precond(NULL),
    x_loaded(false), b_loaded(false)
{
  bho = 0;

  ParmParse pp("hypre");

  verbose = 0;
  pp.query("verbose", verbose);

  solver_flag = 104;
//   pp.query("solver_flag", solver_flag);

  if (solver_flag == 104) {
    ObjectType = HYPRE_PARCSR;
  }
  else {
    amrex::Error("No such solver in HypreABecLap");
  }

  int nparts = fine_level - crse_level + 1;

#if (BL_SPACEDIM == 1)

  // Hypre doesn't support 1D directly, so we use 2D Hypre with
  // the second dimension collapsed.
  // (SMG reduces to cyclic reduction in this case, so it's an exact solve.)
  // (PFMG will not work.)

  HYPRE_SStructGridCreate(MPI_COMM_WORLD, 2, nparts, &hypreGrid);

#else

  HYPRE_SStructGridCreate(MPI_COMM_WORLD, BL_SPACEDIM, nparts, &hypreGrid);

#endif
}


HypreABecLap::~HypreABecLap()
{
  for (int level = crse_level; level <= fine_level; level++) {
    if (level > crse_level) {
      delete cintrp.remove(level);
      delete ederiv.remove(level);
      delete c_cintrp.remove(level);
      delete c_ederiv.remove(level);
      delete c_entry.remove(level);
    }
    delete acoefs.remove(level);
    delete bcoefs.remove(level);

    delete bd.remove(level);
  }

  HYPRE_SStructVectorDestroy(b);
  HYPRE_SStructVectorDestroy(x);
  
  HYPRE_SStructMatrixDestroy(A);

  HYPRE_SStructGraphDestroy(graph);
  HYPRE_SStructStencilDestroy(stencil);
  HYPRE_SStructGridDestroy(hypreGrid);
}


void HypreABecLap::addLevel(int             level,
			    const Geometry& _geom,
			    const BoxArray& _grids,
			    IntVect         _fine_ratio)
{
  int part = level - crse_level;

  geom[level]  = _geom;
  grids[level] = _grids;
  fine_ratio[level] = _fine_ratio;

#if (BL_SPACEDIM == 1)

  if (geom[level].isAnyPeriodic()) {
    BL_ASSERT(geom[level].isPeriodic(0));
    BL_ASSERT(geom[level].Domain().smallEnd(0) == 0);

    int is_periodic[2];
    is_periodic[0] = geom[level].period(0);
    is_periodic[1] = 0;
    BL_ASSERT(ispow2(is_periodic[0]));

    HYPRE_SStructGridSetPeriodic(hypreGrid, part, is_periodic);
  }

#else

  if (geom[level].isAnyPeriodic()) {
    int is_periodic[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) {
      is_periodic[i] = 0;
      if (geom[level].isPeriodic(i)) {
	is_periodic[i] = geom[level].period(i);
	BL_ASSERT(ispow2(is_periodic[i]));
	BL_ASSERT(geom[level].Domain().smallEnd(i) == 0);
      }
    }
    HYPRE_SStructGridSetPeriodic(hypreGrid, part, is_periodic);
  }

#endif

  int num_procs = ParallelDescriptor::NProcs();
  int myid      = ParallelDescriptor::MyProc();
  DistributionMapping distributionMap(grids[level], num_procs);

  for (int i = 0; i < grids[level].size(); i++) {
    if (distributionMap[i] == myid) {
      HYPRE_SStructGridSetExtents(hypreGrid, part,
				  loV(grids[level][i]),
				  hiV(grids[level][i]));
    }
  }

  // All variables are cell-centered
  HYPRE_SStructVariable vars[1] = {HYPRE_SSTRUCT_VARIABLE_CELL};
  HYPRE_SStructGridSetVariables(hypreGrid, part, 1, vars);
}


void HypreABecLap::setBndry(int level, BndryData& _bd)
{
  bd.clear(level);
  bd.set(level, &_bd);
}


static void 
TransverseInterpolant(AuxVarBox& cintrp, const Mask& msk,
                      const Box& reg, const Box& creg,
                      AMREX_D_DECL(const IntVect& rat, const IntVect& vj1, const IntVect& vk1),
                      AMREX_D_DECL(const IntVect& ve,  const IntVect& vjr, const IntVect& vkr),
                      AMREX_D_DECL(int idir, int jdir, int kdir),
                      int clevel)
{
  for (IntVect vc = creg.smallEnd(); vc <= creg.bigEnd(); creg.next(vc)) {
    IntVect vf = rat * vc;
    vf[idir] = reg.smallEnd(idir); // same as bigEnd(idir)
    Box face(vf, vf + ve);
    if (msk(vf) == BndryData::not_covered) {
#if (0)
      // force piecewise constant interpolation for debugging:
      for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
        cintrp(v).push(clevel, vc,     1.0);
      }
#elif (BL_SPACEDIM == 1)
      cintrp(vf).push(clevel, vc, 1.0);
#elif (BL_SPACEDIM == 2)
      if (msk(vf-vj1) != BndryData::not_covered &&
          msk(vf+vjr) == BndryData::not_covered) {
        // low direction not available, use linear interp upwards:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real xx = (v[jdir] - vf[jdir] - 0.5 * ve[jdir]) / rat[jdir];
          cintrp(v).push(clevel, vc+vj1, xx);
          cintrp(v).push(clevel, vc,     1.0 - xx);
        }
      }
      else if (msk(vf-vj1) == BndryData::not_covered &&
               msk(vf+vjr) == BndryData::not_covered) {
        // use piecewise quadratic interpolation whenever possible:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real xx = (v[jdir] - vf[jdir] - 0.5 * ve[jdir]) / rat[jdir];
          cintrp(v).push(clevel, vc+vj1, 0.5*xx*(xx+1));
          cintrp(v).push(clevel, vc,     1.0-xx*xx);
          cintrp(v).push(clevel, vc-vj1, 0.5*xx*(xx-1));
        }
      }
      else if (msk(vf-vj1) == BndryData::not_covered &&
               msk(vf+vjr) != BndryData::not_covered) {
        // high direction not available, use linear interp downwards:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real xx = (v[jdir] - vf[jdir] - 0.5 * ve[jdir]) / rat[jdir];
          cintrp(v).push(clevel, vc,     1.0 + xx);
          cintrp(v).push(clevel, vc-vj1, -xx);
        }
      }
      else {
        // neither direction available, drop back to piecewise const:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          cintrp(v).push(clevel, vc,     1.0);
        }
        //amrex::Error("Case not implemented");
      }
#elif (BL_SPACEDIM == 3)

      // First do the jdir direction, including piecewise-constant term:

      if (msk(vf-vj1) != BndryData::not_covered &&
          msk(vf+vjr) == BndryData::not_covered) {
        // low direction not available, use linear interp upwards:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real xx = (v[jdir] - vf[jdir] - 0.5 * ve[jdir]) / rat[jdir];
          cintrp(v).push(clevel, vc+vj1, xx);
          cintrp(v).push(clevel, vc,     1.0 - xx);
        }
      }
      else if (msk(vf-vj1) == BndryData::not_covered &&
               msk(vf+vjr) == BndryData::not_covered) {
        // use piecewise quadratic interpolation whenever possible:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real xx = (v[jdir] - vf[jdir] - 0.5 * ve[jdir]) / rat[jdir];
          cintrp(v).push(clevel, vc+vj1, 0.5*xx*(xx+1));
          cintrp(v).push(clevel, vc,     1.0-xx*xx);
          cintrp(v).push(clevel, vc-vj1, 0.5*xx*(xx-1));
        }
      }
      else if (msk(vf-vj1) == BndryData::not_covered &&
               msk(vf+vjr) != BndryData::not_covered) {
        // high direction not available, use linear interp downwards:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real xx = (v[jdir] - vf[jdir] - 0.5 * ve[jdir]) / rat[jdir];
          cintrp(v).push(clevel, vc,     1.0 + xx);
          cintrp(v).push(clevel, vc-vj1, -xx);
        }
      }
      else {
        // neither direction available, drop back to piecewise const:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          cintrp(v).push(clevel, vc,     1.0);
        }
      }

      // Then add on contributions from the kdir direction:

      if (msk(vf-vk1) != BndryData::not_covered &&
          msk(vf+vkr) == BndryData::not_covered) {
        // low direction not available, use linear interp upwards:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real yy = (v[kdir] - vf[kdir] - 0.5 * ve[kdir]) / rat[kdir];
          cintrp(v).push(clevel, vc+vk1, yy);
          cintrp(v).push(clevel, vc,    -yy);
        }
      }
      else if (msk(vf-vk1) == BndryData::not_covered &&
               msk(vf+vkr) == BndryData::not_covered) {
        // use piecewise quadratic interpolation whenever possible:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real yy = (v[kdir] - vf[kdir] - 0.5 * ve[kdir]) / rat[kdir];
          cintrp(v).push(clevel, vc+vk1, 0.5*yy*(yy+1));
          cintrp(v).push(clevel, vc,     -yy*yy);
          cintrp(v).push(clevel, vc-vk1, 0.5*yy*(yy-1));
        }
      }
      else if (msk(vf-vk1) == BndryData::not_covered &&
               msk(vf+vkr) != BndryData::not_covered) {
        // high direction not available, use linear interp downwards:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real yy = (v[kdir] - vf[kdir] - 0.5 * ve[kdir]) / rat[kdir];
          cintrp(v).push(clevel, vc,      yy);
          cintrp(v).push(clevel, vc-vk1, -yy);
        }
      }
      else {
        // neither direction available, drop back to piecewise const:
        // (do nothing, no need to explicitly add a zero contribution here)
      }

      // Finally add in cross derivative terms:

      if (msk(vf-vj1-vk1) == BndryData::not_covered &&
          msk(vf-vj1+vkr) == BndryData::not_covered &&
          msk(vf+vjr-vk1) == BndryData::not_covered &&
          msk(vf+vjr+vkr) == BndryData::not_covered) {
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real xx = (v[jdir] - vf[jdir] - 0.5 * ve[jdir]) / rat[jdir];
          Real yy = (v[kdir] - vf[kdir] - 0.5 * ve[kdir]) / rat[kdir];
          cintrp(v).push(clevel, vc-vj1-vk1,  0.25*xx*yy);
          cintrp(v).push(clevel, vc-vj1+vk1, -0.25*xx*yy);
          cintrp(v).push(clevel, vc+vj1-vk1, -0.25*xx*yy);
          cintrp(v).push(clevel, vc+vj1+vk1,  0.25*xx*yy);
        }
      }
#endif
    }
  }
}

static void 
NormalDerivative(AuxVarBox& ederiv, AuxVarBox& cintrp,
                 const Mask& msk, const Box& reg,
                 const IntVect& vin, Real h, int r, int bho, int flevel)
{
  if (bho == 1) {
    Real efacb = 8.0 / (h * (1 + r) * (3 + r));
    Real efac1 = (r - 3) / (h * (1 + r));
    Real efac2 = (1 - r) / (h * (3 + r));
    IntVect vi2 = 2 * vin;
    for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
      if (msk(v) == BndryData::not_covered) {
        ederiv(v).push(&cintrp(v),    efacb);
        ederiv(v).push(flevel, v-vin, efac1);
        ederiv(v).push(flevel, v-vi2, efac2);
      }
    }
  }
  else {
    Real efac = 2.0 / (h * (1 + r)); // normal derivative factor
    for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
      if (msk(v) == BndryData::not_covered) {
        ederiv(v).push(&cintrp(v),     efac);
        ederiv(v).push(flevel, v-vin, -efac);
      }
    }
  }
}


void HypreABecLap::buildMatrixStructure()
{
  // Build MultiFabs first, so that distribution maps are available
  // when building the graph below:

  for (int level = crse_level; level <= fine_level; level++) {
    int ncomp=1;
    int ngrow=0;
    acoefs.set(level, new MultiFab(grids[level], ncomp, ngrow));
    acoefs[level].setVal(0.0);

    bcoefs.set(level, new Tuple<MultiFab, BL_SPACEDIM>);
 
    for (int i = 0; i < BL_SPACEDIM; i++) {
      BoxArray edge_boxes(grids[level]);
      edge_boxes.surroundingNodes(i);
      bcoefs[level][i].define(edge_boxes, ncomp, ngrow, Fab_allocate);
      bcoefs[level][i].setVal(0.0);
    }
  }

  // This can be done now that addLevel has been called for each level:

  HYPRE_SStructGridAssemble(hypreGrid);

  // Setup stencils:

#if (BL_SPACEDIM == 1)
  // fake 1D as a 2D problem:
  int offsets[3][2] = {{-1,  0},      // 0
		       { 1,  0},      // 1
                       { 0,  0}};     // 2
#elif (BL_SPACEDIM == 2)
  int offsets[5][2] = {{-1,  0},      // 0
		       { 0, -1},      // 1
		       { 1,  0},      // 2
		       { 0,  1},      // 3
		       { 0,  0}};     // 4
#elif (BL_SPACEDIM == 3)
  int offsets[7][3] = {{-1,  0,  0},  // 0
		       { 0, -1,  0},  // 1
		       { 0,  0, -1},  // 2
		       { 1,  0,  0},  // 3
		       { 0,  1,  0},  // 4
		       { 0,  0,  1},  // 5
		       { 0,  0,  0}}; // 6
#endif

  BL_ASSERT(stencil == NULL);
#if (BL_SPACEDIM == 1)
  HYPRE_SStructStencilCreate(2, 3, &stencil);
#else
  HYPRE_SStructStencilCreate(BL_SPACEDIM, 2 * BL_SPACEDIM + 1, &stencil);
#endif

  for (int i = 0; i < 2 * BL_SPACEDIM + 1; i++) {
    HYPRE_SStructStencilSetEntry(stencil, i, offsets[i], 0);
  }

  BL_ASSERT(graph == NULL);
  HYPRE_SStructGraphCreate(MPI_COMM_WORLD, hypreGrid, &graph);
  HYPRE_SStructGraphSetObjectType(graph, ObjectType);

  for (int level = crse_level; level <= fine_level; level++) {
    int part = level - crse_level;
    HYPRE_SStructGraphSetStencil(graph, part, 0, stencil);
  }

  // Add non-stencil entries to the graph here:

  for (int level = crse_level + 1; level <= fine_level; level++) {
    int part = level - crse_level;
    cintrp.set(level, new BndryAuxVar(grids[level], BndryAuxVar::GHOST));
    ederiv.set(level, new BndryAuxVar(grids[level], BndryAuxVar::GHOST));
    BndryAuxVar entry(grids[level], BndryAuxVar::INTERIOR);
    IntVect rat = fine_ratio[level-1];
    
    for (OrientationIter oitr; oitr; ++oitr) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      IntVect vin = amrex::BASISV(idir);
      vin = (ori.isLow() ? -vin : vin); // outward normal unit vector
      Real h = geom[level].CellSize(idir); // normal fine grid spacing
      IntVect ve; // default constructor initializes to zero
#if (BL_SPACEDIM >= 2)
      int jdir = (idir + 1) % BL_SPACEDIM;
      IntVect vj1 = amrex::BASISV(jdir); // tangential unit vector
      IntVect vjr = rat * vj1;
      ve += (vjr - vj1);
#endif
#if (BL_SPACEDIM == 3)
      int kdir = (idir + 2) % 3;
      IntVect vk1 = amrex::BASISV(kdir);
      IntVect vkr = rat * vk1;
      ve += (vkr - vk1);
#endif
      for (int i = cintrp[level].firstLocal(); cintrp[level].isValid(i);
           i = cintrp[level].nextLocal(i)) {
        Box reg = amrex::adjCell(grids[level][i], ori);
        Box creg = amrex::coarsen(reg, rat); // coarse adjacent cells
        const Mask &msk = *(bd[level].bndryMasks(i)[ori]);

        TransverseInterpolant(cintrp[level](ori)[i], msk, reg, creg,
                              AMREX_D_DECL(rat,  vj1,  vk1),
                              AMREX_D_DECL(ve,   vjr,  vkr),
                              AMREX_D_DECL(idir, jdir, kdir),
                              level-1);

        NormalDerivative(ederiv[level](ori)[i],
                         cintrp[level](ori)[i],
                         msk, reg, vin, h, rat[idir], bho, level);

        // cintrp and ederiv are now complete.  entry will be done in
        // draft form to establish the graph connections, but must be
        // done again later when the edge coefficients are available.

        reg.shift(-vin); // fine interior cells
        for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
          if (msk(v+vin) == BndryData::not_covered) {
            // value not important, since coefficient not known.
            entry(ori)[i](v).push(&ederiv[level](ori)[i](v+vin), 1.0);
          }
        }
      }
    }

    for (OrientationIter oitr; oitr; ++oitr) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      IntVect vin = amrex::BASISV(idir);
      vin = (ori.isLow() ? -vin : vin); // outward normal unit vector
      for (int i = entry.firstLocal(); entry.isValid(i);
           i = entry.nextLocal(i)) {
        Box reg = amrex::adjCell(grids[level][i], ori);
        reg.shift(-vin); // fine interior cells
        for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
#if (0 && !defined(NDEBUG))
          if (msk(v+vin) == BndryData::not_covered &&
              entry(ori)[i](v).slave()) {
            cout << v << " is slave in orientation " << ori
                 << " on processor " << ParallelDescriptor::MyProc()
                 << endl;
          }
#endif
          // Even if this entry is covered, it could have a
          // not_covered slave:
          if (!entry(ori)[i](v).empty() &&
              !entry(ori)[i](v).slave()) {
            entry(ori)[i](v).collapse();
            Vector<int> levels;
            Vector<IntVect> cells;
            int retval = entry(ori)[i](v).get_locations(levels, cells);
            BL_ASSERT(retval == 0);
            for (int j = 0; j < levels.size(); j++) {
              // eliminate stencil-like connections:
              int not_stencil = 1;
              if (levels[j] == level) {
                IntVect d = cells[j] - v;
                for (int k = 0; k < 2 * BL_SPACEDIM + 1; k++) {
                  if (d == IntVect(offsets[k])) {
                    not_stencil = 0;
                  }
                }
              }
              if (not_stencil) {
                int to_part = levels[j] - crse_level;
                HYPRE_SStructGraphAddEntries(graph,
                                             part,    getV1(v),        0,
                                             to_part, getV2(cells[j]), 0);
              }
            }
          }
        }
      }
    }

    // Now add the graph entries seen by the coarse cells adjacent
    // to the coarse-fine interface.  These are averages of ederiv
    // over the fine faces making up each coarse face.  Since we
    // have to do this from the processor owning the coarse grid,
    // we recompute information using CrseBndryAuxVar

    const BoxArray& f_fgrids(grids[level]);
    const BoxArray& c_cgrids(grids[level-1]);
    BoxArray f_cgrids(c_cgrids);
    f_cgrids.refine(rat);
    BoxArray c_fgrids(f_fgrids);
    c_fgrids.coarsen(rat);

    c_cintrp.set(level, new CrseBndryAuxVar(f_cgrids, f_fgrids,
                                            BndryAuxVar::GHOST));
    //c_ederiv.set(level, new CrseBndryAuxVar(f_cgrids, f_fgrids,
    //                                        BndryAuxVar::GHOST));
    //c_entry.set( level, new CrseBndryAuxVar(c_cgrids, c_fgrids,
    //                                        BndryAuxVar::EXTERIOR));

    c_ederiv.set(level, new CrseBndryAuxVar(c_cintrp[level],
                                            BndryAuxVar::GHOST));
    c_entry.set( level, new CrseBndryAuxVar(c_cgrids, c_fgrids,
                                            c_cintrp[level],
                                            BndryAuxVar::EXTERIOR));
    
    for (OrientationIter oitr; oitr; ++oitr) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      IntVect vin = amrex::BASISV(idir);
      vin = (ori.isLow() ? -vin : vin); // outward normal unit vector
      Real h = geom[level].CellSize(idir); // normal fine grid spacing
      IntVect ve; // default constructor initializes to zero
#if (BL_SPACEDIM >= 2)
      int jdir = (idir + 1) % BL_SPACEDIM;
      IntVect vj1 = amrex::BASISV(jdir); // tangential unit vector
      IntVect vjr = rat * vj1;
      ve += (vjr - vj1);
#endif
#if (BL_SPACEDIM == 3)
      int kdir = (idir + 2) % 3;
      IntVect vk1 = amrex::BASISV(kdir);
      IntVect vkr = rat * vk1;
      ve += (vkr - vk1);
#endif
      for (int i = c_cintrp[level].firstLocal(); c_cintrp[level].isValid(i);
           i = c_cintrp[level].nextLocal(i)) {
        for (int j = 0; j < c_cintrp[level](ori)[i].size(); j++) {
          const Box& reg = c_cintrp[level](ori)[i][j].box(); // adjacent cells
          const Box& creg = c_entry[level](ori)[i][j].box(); // adjacent cells
          const Mask &msk = c_cintrp[level].mask(ori)[i][j]; // fine mask

          TransverseInterpolant(c_cintrp[level](ori)[i][j], msk, reg, creg,
                                AMREX_D_DECL(rat,  vj1,  vk1),
                                AMREX_D_DECL(ve,   vjr,  vkr),
                                AMREX_D_DECL(idir, jdir, kdir),
                                level-1);

          NormalDerivative(c_ederiv[level](ori)[i][j],
                           c_cintrp[level](ori)[i][j],
                           msk, reg, vin, h, rat[idir], bho, level);

          // c_cintrp and c_ederiv are now complete.  c_entry will be done
          // in draft form to establish the graph connections, but must be
          // done again later when the edge coefficients are available.

          for (IntVect vc = creg.smallEnd(); vc <= creg.bigEnd(); creg.next(vc)) {
            IntVect vf = rat * vc;
            vf[idir] = reg.smallEnd(idir); // same as bigEnd(idir)
            Box face(vf, vf + ve);
            if (msk(vf) == BndryData::not_covered) {
              for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
                // value not important, since coefficient not known.
                c_entry[level](ori)[i][j](vc)
                  .push(&c_ederiv[level](ori)[i][j](v), 1.0);
              }
            }
          }
        }
      }
    }

    for (OrientationIter oitr; oitr; ++oitr) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      for (int i = c_cintrp[level].firstLocal(); c_cintrp[level].isValid(i);
           i = c_cintrp[level].nextLocal(i)) {
        for (int j = 0; j < c_cintrp[level](ori)[i].size(); j++) {
          const Box& reg = c_cintrp[level](ori)[i][j].box(); // adjacent cells
          const Box& creg = c_entry[level](ori)[i][j].box(); // adjacent cells
          const Mask &msk = c_cintrp[level].mask(ori)[i][j]; // fine mask
          for (IntVect vc = creg.smallEnd(); vc <= creg.bigEnd(); creg.next(vc)) {
            IntVect vf = rat * vc;
            vf[idir] = reg.smallEnd(idir); // same as bigEnd(idir)
            // Unlike fine entry, it should not be possible for this
            // entry to be covered but have a not_covered slave:
            if (msk(vf) == BndryData::not_covered &&
                !c_entry[level](ori)[i][j](vc).slave()) {
              c_entry[level](ori)[i][j](vc).collapse();
              Vector<int> levels;
              Vector<IntVect> cells;
              int retval = c_entry[level](ori)[i][j](vc)
                .get_locations(levels, cells);
              BL_ASSERT(retval == 0);
              for (int jj = 0; jj < levels.size(); jj++) {
                // eliminate stencil-like connections:
                int not_stencil = 1;
                if (levels[jj] == level-1) {
                  IntVect d = cells[jj] - vc;
                  for (int k = 0; k < 2 * BL_SPACEDIM + 1; k++) {
                    if (d == IntVect(offsets[k])) {
                      not_stencil = 0;
                    }
                  }
                }
                if (not_stencil) {
                  int to_part = levels[jj] - crse_level;
                  HYPRE_SStructGraphAddEntries(graph,
                                               part-1,  getV1(vc),        0,
                                               to_part, getV2(cells[jj]), 0);
                }
              }
            }
          }
        }
      }
    }
  }

  HYPRE_SStructGraphAssemble(graph);

  BL_ASSERT(A == NULL);
  HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, graph, &A);
  HYPRE_SStructMatrixSetObjectType(A, ObjectType);
  //HYPRE_StructMatrixSetSymmetric(A, 1);
  //HYPRE_StructMatrixSetNumGhost(A, A_num_ghost);
  HYPRE_SStructMatrixInitialize(A);

  BL_ASSERT(b == NULL);
  HYPRE_SStructVectorCreate(MPI_COMM_WORLD, hypreGrid, &b);
  HYPRE_SStructVectorSetObjectType(b, ObjectType);

  BL_ASSERT(x == NULL);
  HYPRE_SStructVectorCreate(MPI_COMM_WORLD, hypreGrid, &x);
  HYPRE_SStructVectorSetObjectType(x, ObjectType);

  HYPRE_SStructVectorInitialize(b);
  HYPRE_SStructVectorInitialize(x);

  // According to Rob, the following is necessary in some cases before
  // we call the solver setup routine (and we may do that before loading
  // the vectors with data):

  HYPRE_SStructVectorAssemble(b);
  HYPRE_SStructVectorAssemble(x);
}


void HypreABecLap::setScalars(Real sa, Real sb)
{
  scalar_a = sa;
  scalar_b = sb;
}


void HypreABecLap::setACoeffs(int level, const MultiFab& alpha)
{
  BL_ASSERT( alpha.ok() );
  BL_ASSERT( alpha.boxArray() == acoefs[level].boxArray() );
  for (MFIter ai(alpha); ai.isValid(); ++ai) {
    acoefs[level][ai].copy(alpha[ai]);
  }
}

void HypreABecLap::setBCoeffs(int level, const MultiFab &b, int dir)
{
  BL_ASSERT( b.ok() );
  BL_ASSERT( b.boxArray() == bcoefs[level][dir].boxArray() );
  for (MFIter bi(b); bi.isValid(); ++bi) {
    bcoefs[level][dir][bi].copy(b[bi]);
  }
}

void HypreABecLap::setRhs(int level, const MultiFab& rhs)
{
  int part = level - crse_level;

  for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &reg = grids[level][i];

    FArrayBox *f;
    f = new FArrayBox(reg);
    f->copy(rhs[i]);

    Real* vec = f->dataPtr();

    const Vector< Vector<BoundCond> > & bcs_i = bd[level].bndryConds(i);
    const BndryData::RealTuple      & bcl_i = bd[level].bndryLocs(i);
    const BndryData::MaskTuple      & msk_i = bd[level].bndryMasks(i);

    // add b.c.'s to rhs
    const Box& domain = bd[level].getDomain();
    for (OrientationIter oitr; oitr; oitr++) {
      int cdir(oitr());
      int idim = oitr().coordDir();
      const BoundCond &bct = bcs_i[cdir][0];
      int bctype = bct;
      const Real      &bcl = bcl_i[cdir];
      const FArrayBox &bcv = bd[level].bndryValues(oitr())[i];
      const Mask      &msk = *msk_i[cdir];
      const int* blo = bcoefs[level][idim][mfi].loVect();
      const int* bhi = bcoefs[level][idim][mfi].hiVect();
      const int* mlo = msk.loVect();
      const int* mhi = msk.hiVect();
      const int* bvlo = bcv.loVect();
      const int* bvhi = bcv.hiVect();

      if (reg[oitr()] == domain[oitr()] || level == crse_level) {
	
	// Treat an exposed grid edge here as a boundary condition
	// for the linear solver:

	if (reg[oitr()] == domain[oitr()]) {
	  FORT_HPBVEC3(vec, bcoefs[level][idim][i].dataPtr(), ARLIM(blo), ARLIM(bhi),
		       reg.loVect(), reg.hiVect(), scalar_b, geom[level].CellSize(), 
		       cdir, bctype, bcl, 
		       msk.dataPtr(), ARLIM(mlo), ARLIM(mhi),
		       bcv.dataPtr(), ARLIM(bvlo), ARLIM(bvhi));
	}
	else {
	  FORT_HPBVEC (vec, bcoefs[level][idim][i].dataPtr(), ARLIM(blo), ARLIM(bhi),
		       reg.loVect(), reg.hiVect(), scalar_b, geom[level].CellSize(), 
		       cdir, bctype, bcl, 
		       msk.dataPtr(), ARLIM(mlo), ARLIM(mhi),
		       bcv.dataPtr(), ARLIM(bvlo), ARLIM(bvhi));	  
	}
      }
      // There is no else here, since we would then be at an
      // interior edge and would not need to add anything to vec.
    }

    HYPRE_SStructVectorSetBoxValues(b, part, loV(reg), hiV(reg), 0, vec);

    delete f; 
  }

  b_loaded = true;
}


void HypreABecLap::setInitGuess(int level, const MultiFab& guess)
{
  int part = level - crse_level;

  for (MFIter mfi(guess); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &reg = grids[level][i];

    FArrayBox *f;
    f = new FArrayBox(reg);
    f->copy(guess[mfi], 0, 0, 1);
    Real* vec = f->dataPtr();

    HYPRE_SStructVectorSetBoxValues(x, part, loV(reg), hiV(reg), 0, vec);

    delete f;   
  }

  x_loaded = true;
}


void HypreABecLap::solve(PArray<MultiFab>& soln, Real _reltol, Real _abstol, int _maxiter)
{
  if (!x_loaded || !b_loaded) {
    amrex::Error("Must setRhs and setInitGuess before calling solve");
  }
  
  HYPRE_SStructVectorAssemble(b);
  HYPRE_SStructVectorAssemble(x);

  loadMatrix();
  HYPRE_SStructMatrixAssemble(A);

  abstol = _abstol;
  reltol = _reltol;
  maxiter = _maxiter;

  setupSolver();

  doIt();

  getSolution(soln);

  // print out more information, such as number of iteration, residual, etc.
  // if (verbose > xxx )

  clearSolver();
}

void HypreABecLap::loadMatrix()
{
#if (BL_SPACEDIM == 1)
  // fake 1D as a 2D problem:
  int offsets[3][2] = {{-1,  0},      // 0
		       { 1,  0},      // 1
                       { 0,  0}};     // 2
#elif (BL_SPACEDIM == 2)
  int offsets[5][2] = {{-1,  0},      // 0
		       { 0, -1},      // 1
		       { 1,  0},      // 2
		       { 0,  1},      // 3
		       { 0,  0}};     // 4
#elif (BL_SPACEDIM == 3)
  int offsets[7][3] = {{-1,  0,  0},  // 0
		       { 0, -1,  0},  // 1
		       { 0,  0, -1},  // 2
		       { 1,  0,  0},  // 3
		       { 0,  1,  0},  // 4
		       { 0,  0,  1},  // 5
		       { 0,  0,  0}}; // 6
#endif

  const int size = 2 * BL_SPACEDIM + 1;

  int stencil_indices[size];

  for (int i = 0; i < size; i++) {
    stencil_indices[i] = i;
  }

  Real *mat;
  for (int level = crse_level; level <= fine_level; level++) {
    int part = level - crse_level;

    for (MFIter mfi(acoefs[level]); mfi.isValid(); ++mfi) {
      int i = mfi.index();
      const Box &reg = grids[level][i];

      int volume = reg.numPts();
      mat = hypre_CTAlloc(double, size*volume);

      // build matrix interior

      const int* alo = acoefs[level][mfi].loVect();
      const int* ahi = acoefs[level][mfi].hiVect();
      FORT_HPACOEF(mat, acoefs[level][mfi].dataPtr(), ARLIM(alo), ARLIM(ahi),
		   reg.loVect(), reg.hiVect(), scalar_a);
      
      for (int idim = 0; idim < BL_SPACEDIM; idim++) {
	const int* blo = bcoefs[level][idim][mfi].loVect();
	const int* bhi = bcoefs[level][idim][mfi].hiVect();
	FORT_HPBCOEF(mat, bcoefs[level][idim][mfi].dataPtr(), ARLIM(blo), ARLIM(bhi),
		     reg.loVect(), reg.hiVect(), scalar_b, geom[level].CellSize(), idim);
      }
      
      // add b.c.'s to matrix diagonal, and
      // zero out offdiag values at domain boundaries
      
      const Vector< Vector<BoundCond> > & bcs_i = bd[level].bndryConds(i);
      const BndryData::RealTuple      & bcl_i = bd[level].bndryLocs(i);
      const BndryData::MaskTuple      & msk_i = bd[level].bndryMasks(i);

      const Box& domain = bd[level].getDomain();
      for (OrientationIter oitr; oitr; oitr++) {
	int cdir(oitr());
	int idim = oitr().coordDir();
	const BoundCond &bct = bcs_i[cdir][0];
	int bctype = bct;
	const Real      &bcl = bcl_i[cdir];
	const Mask      &msk = *msk_i[cdir];
	const int* blo = bcoefs[level][idim][mfi].loVect();
	const int* bhi = bcoefs[level][idim][mfi].hiVect();
	const int* mlo = msk.loVect();
	const int* mhi = msk.hiVect();
      
	if (reg[oitr()] == domain[oitr()] || level == crse_level) {

          // Treat an exposed grid edge here as a boundary condition
          // for the linear solver:

	  if (reg[oitr()] == domain[oitr()]) {
	    FORT_HPMAT3(mat, bcoefs[level][idim][mfi].dataPtr(), ARLIM(blo), ARLIM(bhi),
			reg.loVect(), reg.hiVect(), scalar_b, 
			geom[level].CellSize(), cdir, bctype, bcl, 
			msk.dataPtr(), ARLIM(mlo), ARLIM(mhi));
	  }
	  else {
	    FORT_HPMAT(mat, bcoefs[level][idim][mfi].dataPtr(), ARLIM(blo), ARLIM(bhi),
		       reg.loVect(), reg.hiVect(), scalar_b, 
		       geom[level].CellSize(), cdir, bctype, bcl, 
		       msk.dataPtr(), ARLIM(mlo), ARLIM(mhi));	    
	  }
	}
	else {

	  // An exposed grid edge here actually borders the next coarser
          // level in the current linear system.  Zero out the interior
          // stencil using Neumann BC:

	  int bctype_coarse = LO_NEUMANN;
	  FORT_HPMAT(mat, bcoefs[level][idim][mfi].dataPtr(), ARLIM(blo), ARLIM(bhi),
		     reg.loVect(), reg.hiVect(), scalar_b, 
		     geom[level].CellSize(), cdir, bctype_coarse, bcl, 
		     msk.dataPtr(), ARLIM(mlo), ARLIM(mhi));	    
	}	  
      }

      // initialize matrix
      HYPRE_SStructMatrixSetBoxValues(A, part, loV(reg), hiV(reg), 0,
				      size, stencil_indices, mat);
      
      hypre_TFree(mat);
    }

    // Add coarse-fine interface entries to the matrix here:
    if (level == crse_level) {
      continue;
    }

    // First we do the entries as seen by the fine cells adjacent
    // to the interface, working on the fine processor:

    BndryAuxVar entry(grids[level], BndryAuxVar::INTERIOR);
    for (OrientationIter oitr; oitr; ++oitr) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      IntVect vin = amrex::BASISV(idir), ves;
      vin = (ori.isLow() ? -vin : vin);    // outward normal unit vector
      ves = (ori.isLow() ?  ves : vin);    // edge shift vector
      Real h = geom[level].CellSize(idir); // normal fine grid spacing
      Real ffac = (-scalar_b / h);             // divergence factor
      for (int i = cintrp[level].firstLocal(); cintrp[level].isValid(i);
           i = cintrp[level].nextLocal(i)) {
        Box reg = amrex::adjCell(grids[level][i], ori);
        reg.shift(-vin); // fine interior cells
	const Mask &msk = *(bd[level].bndryMasks(i)[ori]);
        for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
          if (msk(v+vin) == BndryData::not_covered) {
            entry(ori)[i](v).push(&ederiv[level](ori)[i](v+vin),
                                  ffac * bcoefs[level][idir][i](v+ves));
          }
        }
      }
    }

    for (OrientationIter oitr; oitr; ++oitr) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      IntVect vin = amrex::BASISV(idir);
      vin = (ori.isLow() ? -vin : vin);    // outward normal unit vector
      for (int i = cintrp[level].firstLocal(); cintrp[level].isValid(i);
           i = cintrp[level].nextLocal(i)) {
        Box reg = amrex::adjCell(grids[level][i], ori);
        reg.shift(-vin); // fine interior cells
        for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
          if (!entry(ori)[i](v).empty() &&
              !entry(ori)[i](v).slave()) {
            entry(ori)[i](v).collapse();
            Vector<int> levels;
            Vector<IntVect> cells;
            int retval = entry(ori)[i](v).get_locations(levels, cells);
            BL_ASSERT(retval == 0);
            Vector<Real> values;
            retval = entry(ori)[i](v).get_coeffs(values);
            BL_ASSERT(retval == 0);
            int ientry = 2 * BL_SPACEDIM + 1;
            for (int j = 0; j < levels.size(); j++) {
              // identify stencil-like connections for separate treatment:
              int not_stencil = 1;
              if (levels[j] == level) {
                IntVect d = cells[j] - v;
                for (int k = 0; k < 2 * BL_SPACEDIM + 1; k++) {
                  if (d == IntVect(offsets[k])) {
                    not_stencil = 0;
                    HYPRE_SStructMatrixAddToValues(A, part, getV1(v), 0,
                                                   1, &k, &values[j]);
                  }
                }
              }
              if (not_stencil) {
                HYPRE_SStructMatrixSetValues(A, part, getV1(v), 0,
                                             1, &ientry, &values[j]);
                ientry++;
              }
            }
          }
        }
      }
    }

    // Now add the matrix values seen by the coarse cells adjacent
    // to the coarse-fine interface.  These are averages of ederiv
    // over the fine faces making up each coarse face.  We use
    // CrseBndryAuxVar now because we need to work on the coarse
    // processor.

    IntVect rat = fine_ratio[level-1];
    const BoxArray& f_fgrids(grids[level]);
    const BoxArray& c_cgrids(grids[level-1]);
    BoxArray f_cgrids(c_cgrids);
    f_cgrids.refine(rat);
    BoxArray c_fgrids(f_fgrids);
    c_fgrids.coarsen(rat);
    //CrseBndryAuxVar c_entry(c_cgrids, c_fgrids, BndryAuxVar::EXTERIOR);
    c_entry[level].reinitialize_connections(BndryAuxVar::EXTERIOR);
    c_entry[level].buildFaceData(rat);

    for (OrientationIter oitr; oitr; ++oitr) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      c_entry[level].loadFaceData(ori, bcoefs[level][idir], 0, 0, 1);
      IntVect vin = amrex::BASISV(idir), ves;
      vin = (ori.isLow() ? -vin : vin); // outward normal unit vector
      ves = (ori.isLow() ? -vin : ves); // edge shift vector (diff from above)
      Real hc = geom[level-1].CellSize(idir); // normal coarse grid spacing
      Real cfac = (scalar_b / hc);                // divergence factor
      //cfac = 0.0;
      Real zfac = (-cfac / hc);               // factor for covered cell
      //cfac = 0.0;
      IntVect ve; // default constructor initializes to zero
#if (BL_SPACEDIM >= 2)
      int jdir = (idir + 1) % BL_SPACEDIM;
      ve += (rat[jdir] - 1) * amrex::BASISV(jdir);
      cfac /= rat[jdir]; // will average over fine cells in tangential dir
#endif
#if (BL_SPACEDIM == 3)
      int kdir = (idir + 2) % 3;
      ve += (rat[kdir] - 1) * amrex::BASISV(kdir);
      cfac /= rat[kdir]; // will average over fine cells in tangential dir
#endif
      for (int i = c_entry[level].firstLocal(); c_entry[level].isValid(i);
           i = c_entry[level].nextLocal(i)) {
        // parallel loop is tied to coarse grids
        for (int j = 0; j < c_entry[level](ori)[i].size(); j++) {
          const Box& reg = c_cintrp[level](ori)[i][j].box(); // adjacent cells
          const Box& creg = c_entry[level](ori)[i][j].box(); // adjacent cells
          const Mask& msk = c_cintrp[level].mask(ori)[i][j]; // fine mask
	  const FArrayBox& fbcoefs = c_entry[level].faceData(ori)[i][j];
          for (IntVect vc = creg.smallEnd(); vc <= creg.bigEnd(); creg.next(vc)) {
            IntVect vf = rat * vc;
            vf[idir] = reg.smallEnd(idir); // same as bigEnd(idir)
            Box face(vf, vf + ve);
            if (msk(vf) == BndryData::not_covered) {
              // Zero out connection to covered coarse cell:
              c_entry[level](ori)[i][j](vc).push(-1, vc-vin, 0.0);
              c_entry[level](ori)[i][j](vc).push(level-1, vc,
                                zfac * bcoefs[level-1][idir][i](vc+ves));
              // Add fine fluxes over face of coarse cell:
              for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
                c_entry[level](ori)[i][j](vc)
                  .push(&c_ederiv[level](ori)[i][j](v), cfac * fbcoefs(v+ves));
              }
            }
          }
        }
      }
    }

    for (OrientationIter oitr; oitr; ++oitr) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      for (int i = c_entry[level].firstLocal(); c_entry[level].isValid(i);
           i = c_entry[level].nextLocal(i)) {
        // parallel loop is tied to coarse grids
        for (int j = 0; j < c_entry[level](ori)[i].size(); j++) {
          const Box& reg = c_cintrp[level](ori)[i][j].box(); // adjacent cells
          const Box& creg = c_entry[level](ori)[i][j].box(); // adjacent cells
          const Mask &msk = c_cintrp[level].mask(ori)[i][j]; // fine mask
          for (IntVect vc = creg.smallEnd(); vc <= creg.bigEnd(); creg.next(vc)) {
            IntVect vf = rat * vc;
            vf[idir] = reg.smallEnd(idir); // same as bigEnd(idir)
            if (msk(vf) == BndryData::not_covered &&
                !c_entry[level](ori)[i][j](vc).slave()) {
              c_entry[level](ori)[i][j](vc).collapse();
              Vector<int> levels;
              Vector<IntVect> cells;
              int retval = c_entry[level](ori)[i][j](vc)
                .get_locations(levels, cells);
              BL_ASSERT(retval == 0);
              Vector<Real> values;
              retval = c_entry[level](ori)[i][j](vc).get_coeffs(values);
              BL_ASSERT(retval == 0);
              int ientry = 2 * BL_SPACEDIM + 1;
              for (int jj = 0; jj < levels.size(); jj++) {
                // identify stencil-like connections for separate treatment:
                int not_stencil = 1;
                if (levels[jj] == -1) {
                  // connection to covered coarse cell to be zeroed out:
                  IntVect d = cells[jj] - vc;
                  for (int k = 0; k < 2 * BL_SPACEDIM + 1; k++) {
                    if (d == IntVect(offsets[k])) {
                      not_stencil = 0;
                      HYPRE_SStructMatrixSetValues(A, part-1, getV1(vc), 0,
                                                   1, &k, &values[jj]);
                      BL_ASSERT(values[jj] == 0.0);
                    }
                  }
                  BL_ASSERT(not_stencil == 0);
                }
                if (levels[jj] == level-1) {
                  // other coarse-level entry, may or may not be in stencil:
                  IntVect d = cells[jj] - vc;
                  for (int k = 0; k < 2 * BL_SPACEDIM + 1; k++) {
                    if (d == IntVect(offsets[k])) {
                      not_stencil = 0;
                      HYPRE_SStructMatrixAddToValues(A, part-1, getV1(vc), 0,
                                                     1, &k, &values[jj]);
                    }
                  }
                }
                if (not_stencil) {
                  HYPRE_SStructMatrixSetValues(A, part-1, getV1(vc), 0,
                                               1, &ientry, &values[jj]);
                  ientry++;
                }
              }
            }
          }
        }
      }
    }
  }
}


void HypreABecLap::setupSolver()
{
  BL_ASSERT(solver         == NULL);
  BL_ASSERT(precond        == NULL);

  if (solver_flag == 104) {

    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b;
    HYPRE_ParVector par_x;

    HYPRE_SStructMatrixGetObject(A, (void**) &par_A);
    HYPRE_SStructVectorGetObject(b, (void**) &par_b);
    HYPRE_SStructVectorGetObject(x, (void**) &par_x);

    ParmParse pp("hypre");
    int kdim = 5; pp.query("kdim", kdim);

    HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);
    HYPRE_ParCSRGMRESSetKDim(solver, kdim);
    HYPRE_ParCSRGMRESSetMaxIter(solver, maxiter);
    HYPRE_ParCSRGMRESSetTol(solver, reltol);
    //HYPRE_ParCSRGMRESSetLogging(solver, 1);

    HYPRE_BoomerAMGCreate(&precond);
    HYPRE_BoomerAMGSetMaxIter(precond, 1);
    HYPRE_BoomerAMGSetTol(precond, reltol);
    //HYPRE_BoomerAMGSetStrongThreshold(precond, 0.6);
    //HYPRE_BoomerAMGSetPrintLevel(precond, 2);
    //HYPRE_BoomerAMGSetTruncFactor(precond, 0.5);

    // original version
    // HYPRE_BoomerAMGSetCoarsenType(precond, 6);
    //
    // first version by Ulrike Meier
    HYPRE_BoomerAMGSetCoarsenType(precond, 10);
    HYPRE_BoomerAMGSetAggNumLevels(precond,2);
    HYPRE_BoomerAMGSetNumPaths(precond,2);
    //
    // second version by Ulrike Meier
    //    HYPRE_BoomerAMGSetCoarsenType(precond, 10);
    //    HYPRE_BoomerAMGSetAggNumLevels(precond,1);
    //

    //HYPRE_BoomerAMGSetLogging(precond, 2);
    HYPRE_ParCSRGMRESSetPrecond(solver,
                                (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSolve,
                                (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSetup,
                                precond);

    const Real run_strt = ParallelDescriptor::second();

    HYPRE_ParCSRGMRESSetup(solver, par_A, par_b, par_x);

    Real run_time = ParallelDescriptor::second() - run_strt;

    ParallelDescriptor::ReduceRealMax(run_time, ParallelDescriptor::IOProcessorNumber());
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "   Time on HYPRE_ParCSRGMRESSetup: " << run_time << std::endl;
    }
  }
  else {
    amrex::Error("No such solver in HypreABecLap");
  }
}


void HypreABecLap::clearSolver()
{
  if (solver_flag == 104) {
    HYPRE_ParCSRGMRESDestroy(solver);
    HYPRE_BoomerAMGDestroy(precond);    
  }
  else {
    amrex::Error("No such solver in HypreABecLap");
  }  
}

void HypreABecLap::doIt()
{
  if (abstol > 0.0) {
    Real bnorm;
    hypre_SStructInnerProd((hypre_SStructVector *) b,
                           (hypre_SStructVector *) b,
                           &bnorm);
    bnorm = sqrt(bnorm);
    
    Real volume = 0.0;
    for (int level = crse_level; level <= fine_level; level++) {
      for (int i = 0; i < grids[level].size(); i++) {
        volume += grids[level][i].numPts();
      }
    }
    
    Real reltol_new = (bnorm > 0.0
		       ? abstol / bnorm * sqrt(volume)
		       : reltol);
    
    if (reltol_new > reltol) {
      if (solver_flag == 104) {
	HYPRE_ParCSRGMRESSetTol(solver, reltol_new);
        HYPRE_BoomerAMGSetTol(precond, reltol_new);
      }
      else {
	amrex::Error("No such solver in HypreABecLap");
      }  
    }
  }

  const Real run_strt = ParallelDescriptor::second();

  if (solver_flag == 104) {
    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b;
    HYPRE_ParVector par_x;
    
    HYPRE_SStructMatrixGetObject(A, (void**) &par_A);
    HYPRE_SStructVectorGetObject(b, (void**) &par_b);
    HYPRE_SStructVectorGetObject(x, (void**) &par_x);
    
    HYPRE_ParCSRGMRESSolve(solver, par_A, par_b, par_x);
  }
  else {
    amrex::Error("No such solver in HypreABecLap");    
  }

  HYPRE_SStructVectorGather(x);

  Real run_time = ParallelDescriptor::second() - run_strt;

  ParallelDescriptor::ReduceRealMax(run_time, ParallelDescriptor::IOProcessorNumber());
  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "   Hypre Solve Time: " << run_time << std::endl;
  }
}


void HypreABecLap::getSolution(PArray<MultiFab>& soln)
{
  for (int level = crse_level; level <= fine_level; level++) {

    int part = level - crse_level;

    for (MFIter mfi(soln[level]); mfi.isValid(); ++mfi) {
      int i = mfi.index();
      const Box &reg = grids[level][i];

      FArrayBox *f;
      if (soln[level].nGrow() == 0) { // need a temporary if dest is the wrong size
	f = &soln[level][mfi];
      }
      else {
	f = new FArrayBox(reg);
      }

      Real* vec = f->dataPtr();

      HYPRE_SStructVectorGetBoxValues(x, part, loV(reg), hiV(reg), 0, vec);

      if (soln[level].nGrow() != 0) {
	soln[level][mfi].copy(*f, 0, 0, 1);
	delete f;
      }
    }
  }
}
