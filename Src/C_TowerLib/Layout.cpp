#include <MFTower.H>
#include <Layout.H>
#include <VisMF.H>

std::ostream& operator<< (std::ostream&  os, const Node& node)
{
  os << "iv:" << node.iv << ",l:" << node.level << ",t:"; 
  if (node.type==Node::INIT)
    os << "I";
  else if (node.type==Node::COVERED)
    os << "C";
  else
    os << "V";
  if (os.fail())
    BoxLib::Error("operator<<(std::ostream&,Node&) failed");
  return os;
}

BoxArray
Layout::GetBndryCells (const BoxArray& ba,
                       const IntVect&  ngrow,
                       const Geometry& geom)
{
  //
  // First get list of all ghost cells.
  //
  BoxList gcells, bcells;

  for (int i = 0; i < ba.size(); ++i) {
    gcells.join(BoxLib::boxDiff(BoxLib::grow(ba[i],ngrow),ba[i]));
  }

  //
  // Now strip out intersections with original BoxArray (or, local boxes)
  //
  for (BoxList::const_iterator it = gcells.begin(); it != gcells.end(); ++it)
  {
    std::vector< std::pair<int,Box> > isects = ba.intersections(*it);

    if (isects.empty())
      bcells.push_back(*it);
    else
    {
      //
      // Collect all the intersection pieces.
      //
      BoxList pieces;
      for (int i = 0; i < isects.size(); i++)
        pieces.push_back(isects[i].second);
      BoxList leftover = BoxLib::complementIn(*it,pieces);
      bcells.catenate(leftover);
    }
  }
  //
  // Now strip out overlaps.
  //
  gcells.clear();
  gcells = BoxLib::removeOverlap(bcells);
  bcells.clear();

  if (geom.isAnyPeriodic())
  {
    Array<IntVect> pshifts(27);

    const Box& domain = geom.Domain();

    for (BoxList::const_iterator it = gcells.begin(); it != gcells.end(); ++it)
    {
      if (!domain.contains(*it))
      {
        //
        // Add in periodic ghost cells shifted to valid region.
        //
        geom.periodicShift(domain, *it, pshifts);

        for (int i = 0; i < pshifts.size(); i++)
        {
          const Box& shftbox = *it + pshifts[i];

          const Box& ovlp = domain & shftbox;
          BoxList bl = BoxLib::complementIn(ovlp,BoxList(ba));
          bcells.catenate(bl);
        }
      }
    }

    gcells.catenate(bcells);
  }

  return BoxArray(gcells);
}

std::ostream& operator<< (std::ostream&  os, const Stencil& a)
{
  os << "(size=" << a.size() << ") ";
  for (Stencil::const_iterator it=a.begin(), End=a.end(); it!=End; ++it) {
    os << "(" << it->first << "):(" << it->second << ") ";
  }
  return os;
}


const Stencil operator*(Real val, const Stencil& lhs)
{
  Stencil ret(lhs);
  ret *= val;
  return ret;
}

const Stencil operator*(const Stencil& lhs, Real val)
{
  Stencil ret(lhs);
  ret *= val;
  return ret;
}

const Stencil operator+(const Stencil& lhs, const Stencil& rhs)
{
  Stencil ret(lhs);
  ret += rhs;
  return ret;
}

const Stencil operator-(const Stencil& lhs, const Stencil& rhs)
{
  Stencil ret(lhs);
  ret -= rhs;
  return ret;
}

Stencil&
Stencil::operator+=(const Stencil& rhs)
{
  for (const_iterator it=rhs.begin(), End=rhs.end(); it!=End; ++it) {
    const_iterator mit=find(it->first);
    if (mit==end()) {
      insert(std::pair<Node,Real>(it->first,it->second));
    } else {
      (*this)[it->first] += it->second;
    }
  }
  return *this;
}

Stencil&
Stencil::operator-=(const Stencil& rhs)
{
  for (const_iterator it=rhs.begin(), End=rhs.end(); it!=End; ++it) {
    const_iterator mit=find(it->first);
    if (mit==end()) {
      insert(std::pair<Node,Real>(it->first,-it->second));
    } else {
      (*this)[it->first] -= it->second;
    }
  }
  return *this;
}

Stencil&
Stencil::operator*=(Real val)
{
  for (const_iterator it=begin(), End=end(); it!=End; ++it) {
    (*this)[it->first] *= val;
  }
  return *this;
}

Layout::Layout(const Array<BoxArray>&  aba,
               PArray<Geometry>& ag,
               const Array<IntVect>&   ar)
  : nGrow(1), initialized(false)
{
  Build(aba,ag,ar);
}

Layout::~Layout()
{
  Clear();
}

void
Layout::Clear()
{
  nodes.clear();
  crseNodes.clear();
  nodeIds.clear();
  crseIds.clear();
  geomArray.clear();
  gridArray.clear();
  volume.clear();
  for (int d=0; d<area.size(); ++d) {
    area[d].clear();
  }
  area.clear();
  initialized = false;
}

void 
Write(const Layout::MultiIntFab& imf,
      const std::string& filename,
      int comp,
      int nGrow) 
{
  BL_ASSERT(comp<=imf.nComp() && comp>=0);
  BL_ASSERT(nGrow<=imf.nGrow());
  MultiFab mf(imf.boxArray(), 1, imf.nGrow());
  for (MFIter mfi(imf); mfi.isValid(); ++mfi)
  {
    const int* iptr = imf[mfi].dataPtr(comp);
    Real* rptr = mf[mfi].dataPtr();
    int numpts = imf[mfi].box().numPts();
    for (int i=0; i<numpts; ++i) {
      rptr[i] = Real(iptr[i]);
    }
  }
  VisMF::Write(mf,filename);
}

void
Layout::BuildMetrics ()
{
  //
  // Build volume and face area arrays.
  //
  volume.resize(nLevs,PArrayManage);
  area.resize(BL_SPACEDIM);
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
  {
    area[dir].clear();
    area[dir].resize(nLevs,PArrayManage);
  }

  for (int lev=0; lev<nLevs; ++lev)
  {
    const Geometry& geom = geomArray[lev];
    volume.set(lev, new MultiFab());

    geom.GetVolume(volume[lev],gridArray[lev],nGrow);
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      area[dir].set(lev, new MultiFab());
      geom.GetFaceArea(area[dir][lev],gridArray[lev],dir,nGrow);
    }
  }
}

const MultiFab&
Layout::Volume(int lev) const
{
  return volume[lev];
}

const MultiFab&
Layout::Area(int lev, int dir) const
{
  return area[dir][lev];
}

void 
Layout::Build(const Array<BoxArray>&  aba,
              PArray<Geometry>& ag,
              const Array<IntVect>&   ar)
{
  nLevs = aba.size();
  BL_ASSERT(ag.size()>=nLevs);
  BL_ASSERT(ar.size()>=nLevs-1);

  Clear();

  gridArray.resize(nLevs);
  geomArray.resize(nLevs,PArrayManage);
  if (nLevs>1) {
    refRatio.resize(nLevs-1);
  }
  for (int i=0; i<nLevs; ++i) {
    gridArray[i] =  aba[i];
    geomArray.set(i,ag.remove(i));
    if (i<nLevs-1) {
      refRatio[i] = ar[i];
    }
  }

  BuildMetrics();

  NodeFab fnodeFab;

  nodes.resize(nLevs,PArrayManage);
  nodeIds.resize(nLevs,PArrayManage);
  bndryCells.resize(nLevs);
  crseNodes.resize(nLevs,PArrayManage);

  int cnt = 0; // local node counter
  Array<Array<int> > grid_cnt(nLevs);
  for (int lev=0; lev<nLevs; ++lev)
  {
    nodes.set(lev,new MultiNodeFab(gridArray[lev],1,nGrow));
    MultiNodeFab& nodesLev = nodes[lev];
    nodeIds.set(lev,new MultiIntFab(gridArray[lev],1,nGrow));


    for (MFIter fai(nodesLev); fai.isValid(); ++fai)
    {
      NodeFab&  ifab = nodesLev[fai];
      const Box& box = ifab.box() & gridArray[lev][fai.index()];
      for (IntVect iv=box.smallEnd(); iv<=box.bigEnd(); box.next(iv))
        ifab(iv,0) = Node(iv,lev,Node::VALID);
    }
            
    if (lev > 0)
    {
      BoxArray cba = BoxArray(gridArray[lev]).coarsen(refRatio[lev-1]);
      BoxArray cgba = BoxArray(cba).grow(1);
      crseNodes.set(lev,new MultiNodeFab(cgba,1,0));
      for (MFIter mfi(crseNodes[lev]); mfi.isValid(); ++mfi)
      {
        NodeFab& nfab=crseNodes[lev][mfi];
        const Box& box = nfab.box() & geomArray[lev-1].Domain();
        for (IntVect iv=box.smallEnd(); iv<=box.bigEnd(); box.next(iv)) {
          nfab(iv,0) = Node(iv,lev-1,Node::VALID);
        }
      }

      for (MFIter mfi(crseNodes[lev]); mfi.isValid(); ++mfi)
      {
        NodeFab& nfab=crseNodes[lev][mfi];
        const Box& box = nfab.box() & geomArray[lev-1].Domain();
        std::vector< std::pair<int,Box> > isects = cba.intersections(box);
        for (int i=0; i<isects.size(); ++i) {
          const Box& ibox=isects[i].second;
          for (IntVect iv=ibox.smallEnd(), End=ibox.bigEnd(); iv<=End; ibox.next(iv)) {
            nfab(iv,0).type = Node::COVERED;
          }
        }
      }

      const Box& rangeBox = Box(IntVect::TheZeroVector(),
                               refRatio[lev-1] - IntVect::TheUnitVector());
      bndryCells[lev] = BoxLib::intersect( GetBndryCells(nodesLev.boxArray(),refRatio[lev-1],geomArray[lev]),
                                           geomArray[lev].Domain() );
      BoxArray bndC = BoxArray(bndryCells[lev]).coarsen(refRatio[lev-1]);

      for (MFIter mfi(nodesLev); mfi.isValid(); ++mfi)
      {
        NodeFab& nodeFab = nodesLev[mfi];
        const Box& box = nodeFab.box();
        Box gbox = Box(mfi.validbox()).grow(refRatio[lev-1]);
        std::vector< std::pair<int,Box> > isects = bndryCells[lev].intersections(gbox);                
        for (int i=0; i<isects.size(); ++i) {
          const Box& fbox = isects[i].second;
          Box ovlp = fbox & box;
          if (ovlp.ok()) {
            fnodeFab.resize(fbox,1);
            Box cbox = Box(fbox).coarsen(refRatio[lev-1]);
            for (IntVect civ = cbox.smallEnd(), End=cbox.bigEnd(); civ<=End; cbox.next(civ)) {
              const IntVect& baseIV = refRatio[lev-1] * civ;
              for (IntVect ivt = rangeBox.smallEnd(), End=rangeBox.bigEnd(); ivt<=End;rangeBox.next(ivt)) {
                fnodeFab(baseIV + ivt,0) = crseNodes[lev][mfi](civ,0);
              }
            }
            nodeFab.copy(fnodeFab);
          }
        }
      }
    }

    // Block out cells covered by finer grid
    if (lev < nLevs-1)
    {
      const BoxArray coarsenedFineBoxes =
        BoxArray(gridArray[lev+1]).coarsen(refRatio[lev]);
                
      for (MFIter fai(nodesLev); fai.isValid(); ++fai)
      {
        NodeFab& ifab = nodesLev[fai];
        const Box& box = ifab.box();
        std::vector< std::pair<int,Box> > isects = coarsenedFineBoxes.intersections(box);
        for (int i = 0; i < isects.size(); i++)
        {
          const Box& ovlp = isects[i].second;
          for (IntVect iv=ovlp.smallEnd(); iv<=ovlp.bigEnd(); ovlp.next(iv))
            ifab(iv,0) = Node(iv,lev,Node::COVERED);
        }
      }
    }

    // Set nodeIds
    grid_cnt[lev].resize(nodesLev.size(),0);
    nodeIds[lev].setVal(-1);
    for (MFIter fai(nodesLev); fai.isValid(); ++fai)
    {
      NodeFab& ifab = nodes[lev][fai];
      IntFab& nfab = nodeIds[lev][fai];
      const Box& box = fai.validbox() & gridArray[lev][fai.index()];

      int start_cnt = cnt;
      for (IntVect iv(box.smallEnd()); iv<=box.bigEnd(); box.next(iv))
      {
        if (ifab(iv,0).type == Node::VALID)
        {
          if (ifab(iv,0).level != lev) 
            std::cout << "bad level: " << ifab(iv,0) << std::endl;
          nfab(iv,0) = cnt++;
        }
      }
      grid_cnt[lev][fai.index()] = cnt - start_cnt;
    }
    ParallelDescriptor::ReduceIntSum(grid_cnt[lev].dataPtr(),grid_cnt[lev].size());
  }
  nNodes_local = cnt;

#if 0
  // Note: This is an attempt to make the node number independent of the cpu partitioning.
  //       It is not yet functional because it needs an accompanying local-to-global mapping
  //       so that global numbers get mapped correctly to local numbers under the PETSc
  //       covers.  I havet done that yet, partly because I havet yet figured out how to 
  //       make sure the matrix is mapped the same way.  We punt for now, but leave the 
  //       code here for future finishing....grep for VecSetValuesLocal elsewhere

  // Compute total number of nodes
  nNodes_global = 0;
  for (int lev=0; lev<nLevs; ++lev) {
    for (int i=0; i<grid_cnt[lev].size(); ++i) {
      nNodes_global += grid_cnt[lev][i];
    }
  }

  // Replace grid_cnt with offset to first cell in box
  int n_offset = nNodes_global;
  for (int lev=nLevs-1; lev>=0; --lev) {
    for (int i=grid_cnt[lev].size()-1; i>=0; --i) {
      n_offset -= grid_cnt[lev][i];
      grid_cnt[lev][i] = n_offset;
    }
  }

  // Adjust node numbers for cells I own
  for (int lev=0; lev<nLevs; ++lev)
  {
    for (MFIter mfi(nodeIds[lev]); mfi.isValid(); ++mfi)
    {
      IntFab& nfab = nodeIds[lev][mfi];
      const Box& box = mfi.validbox() & gridArray[lev][mfi.index()];
      bool first = true;
      int offset;
      for (IntVect iv(box.smallEnd()); iv<=box.bigEnd(); box.next(iv))
      {
        if (nfab(iv,0) >= 0)
        {
          if (first) {
            first = false;
            offset = grid_cnt[lev][mfi.index()] - nfab(iv,0);
          }
          nfab(iv,0) += offset;
          if (nfab(iv,0) >= nNodes_global) {
            BoxLib::Error("Bad counting");
          }

        }
      }
    }
  }
#else

  nNodes_global = nNodes_local;

#if BL_USE_MPI
  // Adjust node numbers to be globally unique
  int num_size = ParallelDescriptor::NProcs();
  Array<int> num_nodes_p(num_size);
  BL_COMM_PROFILE(BLProfiler::Allgather, sizeof(int), BLProfiler::BeforeCall(),
                  BLProfiler::NoTag());
  BL_MPI_REQUIRE(MPI_Allgather(&nNodes_local, 1, MPI_INT, num_nodes_p.dataPtr(), 1, MPI_INT, 
                               ParallelDescriptor::Communicator()));
  BL_COMM_PROFILE(BLProfiler::Allgather, sizeof(int), BLProfiler::AfterCall(),
                  BLProfiler::NoTag());

  int offset = 0;
  for (int i=0; i<ParallelDescriptor::MyProc(); ++i) {
    offset += num_nodes_p[i];
  }

  // Adjust node numbers 
  for (int lev=0; lev<nLevs; ++lev)
  {
    for (MFIter mfi(nodeIds[lev]); mfi.isValid(); ++mfi)
    {
      IntFab& nfab = nodeIds[lev][mfi];
      const Box& box = mfi.validbox() & gridArray[lev][mfi.index()];
      for (IntVect iv(box.smallEnd()); iv<=box.bigEnd(); box.next(iv))
      {
        if (nfab(iv,0) >= 0)
        {
          nfab(iv,0) += offset;
        }
      }
    }
  }

  // Compute total number of nodes
  nNodes_global = 0;
  for (int i=0; i<ParallelDescriptor::NProcs(); ++i) {
    nNodes_global += num_nodes_p[i];
  }
#endif

#endif

  // Now communicate node numbering to neighbors grow cells  
  if (nLevs>1) {
    crseIds.resize(nLevs-1,PArrayManage);
    for (int lev=0; lev<nLevs; ++lev)
    {
      if (lev>0) 
      {
        const IntVect& ref = refRatio[lev-1];
        const IntVect& refm = ref - IntVect::TheUnitVector();
        BoxArray bndC = BoxArray(bndryCells[lev]).coarsen(ref);
        
        crseIds.set(lev-1,new MultiIntFab(bndC,1,0,Fab_allocate)); crseIds[lev-1].setVal(-1);
        crseIds[lev-1].copy(nodeIds[lev-1]); // parallel copy

        // "refine" crseIds using piecewise constant
        MultiIntFab fineIds(bndryCells[lev],1,0); fineIds.setVal(-1);
        for (MFIter mfi(crseIds[lev-1]); mfi.isValid(); ++mfi)
        {
          const IntFab& cIds = crseIds[lev-1][mfi];
          const Box& cbox = cIds.box();
          IntFab& fIds = fineIds[mfi];
          for (IntVect iv = cbox.smallEnd(), End=cbox.bigEnd(); iv<=End; cbox.next(iv)) {
            const IntVect& ll = ref * iv;
            const Box& fbox(ll,ll+refm);
            BL_ASSERT(fIds.box().contains(fbox));
            fIds.setVal(cIds(iv,0),fbox,0,1);
          }
        }

        nodeIds[lev].FillBoundary(0,1,geomArray[lev].periodicity()); 

        MultiIntFab ng(BoxArray(nodeIds[lev].boxArray()).grow(nodeIds[lev].nGrow()),1,0);
        for (MFIter mfi(nodeIds[lev]); mfi.isValid(); ++mfi) {
          ng[mfi].copy(nodeIds[lev][mfi]); // get valid + f-f (+periodic f-f) but NOT bndryCells
        }
        ng.copy(fineIds); // Parallel copy to get c-f from bndryCells (only)

        for (MFIter mfi(nodeIds[lev]); mfi.isValid(); ++mfi)
        {
          nodeIds[lev][mfi].copy(ng[mfi]); // put it all back
        }
      }
    }
  }
  else {
    for (int lev=0; lev<nLevs; ++lev)
    {
	nodeIds[lev].FillBoundary(0,1,geomArray[lev].periodicity());
    }
  }
  initialized = true;
}

void
Layout::SetNodeIds(BaseFab<int>& idFab, int lev, int grid, const Box& box) const
{
  BL_ASSERT(initialized);
  BL_ASSERT(lev < nodes.size());
  const MultiIntFab& nodeIds_lev = nodeIds[lev];
  const DistributionMapping dm = nodeIds_lev.DistributionMap();
  BL_ASSERT(grid < nodeIds_lev.size());
  BL_ASSERT(dm[grid] == ParallelDescriptor::MyProc());

  const IntFab& nFab = nodeIds_lev[grid];
  Box ovlp = nFab.box() & idFab.box() & box;
  if (ovlp.ok()) {
    idFab.copy(nFab,ovlp,0,ovlp,0,1);
  }
}

bool
Layout::IsCompatible(const MFTower& mft) const
{
  int numLevs = mft.NumLevels();
  bool isok = numLevs <= nLevs;
  if (isok) {
    for (int lev=0; lev<numLevs; ++lev)
    {
      isok &= mft.GetLayout().GridArray()[lev]==gridArray[lev];
    }
  }
  return isok;
}

PetscErrorCode
Layout::MFTowerToVec(Vec&           V,
                     const MFTower& mft,
                     int            comp)
{
  BL_ASSERT(initialized);
  BL_ASSERT(IsCompatible(mft));

  int numLevs = mft.NumLevels();
  BL_ASSERT(numLevs <= nLevs);

  PetscErrorCode ierr;
  FArrayBox fab;
  IntFab ifab;
  for (int lev=0; lev<numLevs; ++lev)
  {
    BL_ASSERT(comp<mft[lev].nComp()+mft.BaseComp());
    for (MFIter mfi(nodeIds[lev]); mfi.isValid(); ++mfi)
    {
      // NOTE: We built the node numbering so that the following operation works correctly
      const Box& vbox = mfi.validbox();

      BoxArray ba(vbox);
      if (lev<numLevs-1) {
        BoxArray cfba = BoxArray(mft[lev+1].boxArray()).coarsen(refRatio[lev]);
        ba = BoxLib::complementIn(vbox,cfba);
      }

      for (int i=0; i<ba.size(); ++i) {
        const Box& box = ba[i];
        int ni = box.numPts();

        ifab.resize(box,1);
        ifab.copy(nodeIds[lev][mfi]);
        BL_ASSERT(ifab.min()>=0);

        const int* ix = ifab.dataPtr();
                
        fab.resize(box,1);
        fab.copy(mft[lev][mfi],comp+mft.BaseComp(),0,1);
        const Real* y = fab.dataPtr();
                
        ierr = VecSetValues(V,ni,ix,y,INSERT_VALUES); CHKPETSC(ierr);                
        //ierr = VecSetValuesLocal(V,ni,ix,y,INSERT_VALUES); CHKPETSC(ierr);                
      }
    }
  }
  ierr = VecAssemblyBegin(V); CHKPETSC(ierr);
  ierr = VecAssemblyEnd(V); CHKPETSC(ierr);
  return ierr;
}

PetscErrorCode
Layout::VecToMFTower(MFTower&   mft,
                     const Vec& V,
                     int        comp) const
{
  BL_ASSERT(initialized);
  BL_ASSERT(IsCompatible(mft));

  int numLevs = mft.NumLevels();
  BL_ASSERT(numLevs<=nLevs);

  PetscErrorCode ierr;

  FArrayBox fab;
  IntFab ifab;
  ierr = VecAssemblyBegin(V); CHKPETSC(ierr);
  ierr = VecAssemblyEnd(V); CHKPETSC(ierr);
  for (int lev=0; lev<numLevs; ++lev)
  {
    BL_ASSERT(comp<mft[lev].nComp()+mft.BaseComp());

    if (mft[lev].nGrow()) {
      mft[lev].setBndry(0,comp+mft.BaseComp(),1);
    }

    for (MFIter mfi(nodeIds[lev]); mfi.isValid(); ++mfi)
    {
      const Box& vbox = mfi.validbox();
      BoxArray ba(vbox);
      if (lev<nLevs-1) {
        BoxArray cfba = BoxArray(mft[lev+1].boxArray()).coarsen(refRatio[lev]);
        ba = BoxLib::complementIn(vbox,cfba);
      }

      for (int i=0; i<ba.size(); ++i) {
        const Box& box = ba[i];
        int ni = box.numPts();

        ifab.resize(box,1);
        ifab.copy(nodeIds[lev][mfi]);
        const int* ix = ifab.dataPtr();
                
        fab.resize(box,1);
        Real* y = fab.dataPtr();

        ierr = VecGetValues(V,ni,ix,y); CHKPETSC(ierr);
        mft[lev][mfi].copy(fab,0,comp+mft.BaseComp(),1);
      }
    }
  }
  return ierr;
}

std::ostream& operator<< (std::ostream&  os, const Layout::IntFab& ifab)
{
  os << "IntFab" << '\n';
  const Box& box = ifab.box();
  os << "Box: " << box << '\n';
  for (IntVect iv(box.smallEnd()); iv<=box.bigEnd(); box.next(iv))
  {
    os << iv << " ";
    for (int n=0; n<ifab.nComp(); ++n) {
      os << ifab(iv,n) << " ";
    }
    os << '\n';
  } 
  return os;
}

std::ostream& operator<< (std::ostream&  os, const Layout::NodeFab& nfab)
{
  os << "NodeFab" << '\n';
  const Box& box = nfab.box();
  os << "Box: " << box << '\n';
  for (IntVect iv(box.smallEnd()); iv<=box.bigEnd(); box.next(iv))
  {
    os << iv << " ";
    for (int n=0; n<nfab.nComp(); ++n) {
      os << nfab(iv,n) << " ";
    }
    os << '\n';
  } 
  return os;
}

std::ostream& operator<< (std::ostream&  os, const Layout::MultiNodeFab& mnf)
{
  os << "MultiNodeFab" << '\n';
  const BoxArray& ba = mnf.boxArray();
  os << "BoxArray: " << ba << '\n';
  const DistributionMapping& dm = mnf.DistributionMap();
  for (int i=0; i<mnf.size(); ++i) {
    if (dm[i]==ParallelDescriptor::MyProc()) {
      os << mnf[i] << std::endl;
    }
    ParallelDescriptor::Barrier();
  }
  return os;
}


