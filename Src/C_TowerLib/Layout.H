#ifndef _Layout_H_
#define _Layout_H_

#include <winstd.H>
#include <petsc.h>

#include <PArray.H>
#include <Geometry.H>

#include <map>

class MFTower;

#define CHKPETSC(n) CHKERRABORT(ParallelDescriptor::Communicator(),n);

struct Node
{
  enum typeEnum{INIT=0, COVERED=1, VALID=2};
  Node()
    : iv(IntVect(D_DECL(-1,-1,-1))), level(-1), type(Node::INIT) {}
  Node(const IntVect& idx, int lev, typeEnum typ = INIT)
    : iv(idx), level(lev), type(typ) {}
  inline bool operator< (const Node& rhs) const
    {
      if (level < rhs.level) return true;
      if ((level == rhs.level) && IntVect::Compare()(iv,rhs.iv)) return true;
      return false;
    }
  inline bool operator!= (const Node& rhs) const
    {
      return ((*this) < rhs || rhs < (*this));
    }
  IntVect iv;
  int level;
  typeEnum type;
};


class Stencil
  : public std::map<Node,Real>
{
public:
  Stencil() {}
  typedef std::map<Node,Real>::const_iterator const_iterator;
  typedef std::map<Node,Real>::iterator iterator;
  Stencil& operator*=(Real val);
  Stencil& operator+=(const Stencil& rhs);
  Stencil& operator-=(const Stencil& rhs);
};

struct IVSMap
  : std::map<IntVect,Stencil,IntVect::Compare>
{
  // Really a typedef, but compiler was unhappy in some cases
};

typedef IVSMap::const_iterator IVScit;

const Stencil operator*(Real val, const Stencil& lhs);
const Stencil operator*(const Stencil& lhs, Real val);
const Stencil operator+(const Stencil& lhs, const Stencil& rhs);
const Stencil operator-(const Stencil& lhs, const Stencil& rhs);

struct Layout
{
  typedef BaseFab<Node> NodeFab;
  typedef FabArray<NodeFab> MultiNodeFab; 
  typedef BaseFab<int> IntFab;
  typedef FabArray<IntFab> MultiIntFab; 

  Layout(const Array<BoxArray>& aba,
         PArray<Geometry>& ag,
         const Array<IntVect>&   ar);

  ~Layout();
  void Clear();
  void SetNodeIds(BaseFab<int>& idFab, int lev, int grid, const Box& box) const;
  int NumberOfGlobalNodeIds() const {return nNodes_global;}
  int NumberOfLocalNodeIds() const {return nNodes_local;}
  int NumLevels() const {return nLevs;}
  const Array<BoxArray>& GridArray() const {return gridArray;}
  const PArray<Geometry>& GeomArray() const {return geomArray;}
  const Array<IntVect>& RefRatio() const {return refRatio;}
  const DistributionMapping& DistributionMap(int lev) const {return nodeIds[lev].DistributionMap();}
  const PArray<Layout::MultiNodeFab>& Nodes() const {return nodes;}
  const PArray<Layout::MultiNodeFab>& CrseNodes() const {return crseNodes;}
  const PArray<Layout::MultiIntFab>& NodeIds() const {return nodeIds;}
  const Array<BoxArray>& BndryCells() const {return bndryCells;}
  static BoxArray GetBndryCells (const BoxArray& ba,
                                 const IntVect&  ngrow,
                                 const Geometry& geom);
  
  bool IsCompatible(const MFTower& mft) const;
  const MultiFab& Volume(int lev) const;
  const MultiFab& Area(int lev, int dir) const;

  PetscErrorCode MFTowerToVec(Vec&           V,
                              const MFTower& mft,
                              int            comp);
  PetscErrorCode VecToMFTower(MFTower&   mft,
                              const Vec& V,
                              int        comp) const;


protected:
  void Build(const Array<BoxArray>& aba,
             PArray<Geometry>& ag,
             const Array<IntVect>&  ar);
  void BuildMetrics();

  PArray<MultiNodeFab> nodes;
  PArray<MultiNodeFab> crseNodes;
  PArray<MultiIntFab> nodeIds;
  PArray<MultiIntFab> crseIds;
  Array<BoxArray> bndryCells;

  PArray<Geometry> geomArray;
  Array<BoxArray> gridArray;
  Array<IntVect> refRatio;
  int nGrow, nLevs;
  bool initialized;
  int nNodes_global, nNodes_local;

  PArray<MultiFab> volume;
  Array<PArray<MultiFab> > area;
};


std::ostream& operator<< (std::ostream&  os, const Stencil& a);
std::ostream& operator<< (std::ostream&  os, const Node& node);
std::ostream& operator<< (std::ostream&  os, const Layout::IntFab& ifab);
std::ostream& operator<< (std::ostream&  os, const Layout::NodeFab& nfab);
std::ostream& operator<< (std::ostream&  os, const Layout::MultiNodeFab& mnf);

#endif
