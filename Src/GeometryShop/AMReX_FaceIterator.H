#ifndef AMREX_FACEITERATOR_H_
#define AMREX_FACEITERATOR_H_

#include "AMReX_IntVect.H"
#include "AMReX_FaceIndex.H"
#include "AMReX_IntVectSet.H"

namespace amrex
{
  class EBGraph;

///
/**
   Ennumeration class to distinguish where FaceIterator will stop.
*/
  class FaceStop
  {
  public:
    ///
    /**
       SurroundingWithBoundary  means stop at all faces on the high and
       low sides of IVS cells. \                                          \
       SurroundingNoBoundary  means stop at all faces on the high and
       low sides of IVS cells, excluding faces on the domain boundary.\   \
       LoWithBoundary  means stop at all faces on the
       low side of IVS cells.\\
       LoNoBoundary  means stop at all faces on the
       low side of IVS cells, excluding faces on the domain boundary.\    \
       HiWithBoundary  means stop at all faces on the
       high side of IVS cells.\\
       HiNoBoundary  means stop at all faces on the
       high side of IVS cells, excluding faces on the domain boundary.\   \
    */
    enum WhichFaces
    {
      Invalid=-1,
      SurroundingWithBoundary=0,
      SurroundingNoBoundary=1,
      HiWithBoundary,
      LoWithBoundary,
      HiNoBoundary,
      LoNoBoundary,
      LoBoundaryOnly,
      HiBoundaryOnly,
      AllBoundaryOnly,
      NUMTYPES
    };
  };


/// Iterator over faces within an IntVectSet and an Ebgraph
/**
   Iterator over within an IntVectSet and an Ebgraph.
   The enumeration tells which faces on which to stop.
*/
  class FaceIterator
  {
  public:

    /// general construction
    FaceIterator(const IntVectSet& a_ivs, const EBGraph& a_ebgraph,
                 int a_direction, const FaceStop::WhichFaces& a_location);

    ///
    FaceIterator();

    /// Going with the default...
    // FaceIterator( const FaceIterator & that);

    /// Going with the default...
    // FaceIterator & operator=( const FaceIterator & that);

    ///
    ~FaceIterator();


    ///
    void define(const IntVectSet& a_ivs, const EBGraph& a_ebgraph,
                int a_direction, const FaceStop::WhichFaces& a_location);

    ///
    void reset();

    ///
    void operator++();

    ///
    const FaceIndex& operator() () const;

    ///
    bool ok() const;

    ///
    bool isDefined() const;

    ///
    const Vector<FaceIndex>& getVector() const;

  private:
    ///
    /**
       Does the real work.  this->define() looks in the cache first.
    */
    void doDefine(const IntVectSet& a_ivs, const EBGraph& a_ebgraph,
                  int a_direction, const FaceStop::WhichFaces& a_location);

    bool                   m_isDefined;
    Vector<FaceIndex> m_faces;
    int                    m_iface;
    int                    m_direction;


  };
}

#endif
