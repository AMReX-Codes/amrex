/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */


#include "AMReX_FaceIterator.H"
#include "AMReX_EBGraph.H"

#include <set>

namespace amrex
{


/********************************/
  const Array<FaceIndex>&
  FaceIterator::getVector() const
  {
    return m_faces;
  }
/********************************/
  FaceIterator::FaceIterator()
  {
    m_isDefined = false;
  }
/********************************/
  FaceIterator::~FaceIterator()
  {
  }
/********************************/
  FaceIterator::FaceIterator(const IntVectSet&           a_ivs,
                             const EBGraph&              a_ebgraph,
                             int                         a_direction,
                             const FaceStop::WhichFaces& a_location)
  {
    define(a_ivs, a_ebgraph, a_direction, a_location);
  }




/********************************/
  void
  FaceIterator::define(const IntVectSet&           a_ivs,
                       const EBGraph&              a_ebgraph,
                       int                         a_direction,
                       const FaceStop::WhichFaces& a_location)
  {

    doDefine(a_ivs, a_ebgraph, a_direction, a_location);
  }

/********************************/
  void
  FaceIterator::doDefine(const IntVectSet&           a_ivs,
                         const EBGraph&              a_ebgraph,
                         int                         a_direction,
                         const FaceStop::WhichFaces& a_location)
  {
    //can't do this because minbox is broken
    //  CH_assert((a_ebgraph.getRegion().contains(a_ivs.minBox())||
    //          a_ivs.isEmpty()));
    m_isDefined = true;
    m_direction = a_direction;

    std::set<FaceIndex, std::less<FaceIndex>>  resultSet;

    bool doLo = ((a_location == FaceStop::SurroundingNoBoundary) ||
                 (a_location == FaceStop::SurroundingWithBoundary) ||
                 (a_location == FaceStop::LoWithBoundary) ||
                 (a_location == FaceStop::LoBoundaryOnly) ||
                 (a_location == FaceStop::AllBoundaryOnly) ||
                 (a_location == FaceStop::LoNoBoundary));
    bool doHi = ((a_location == FaceStop::SurroundingNoBoundary) ||
                 (a_location == FaceStop::SurroundingWithBoundary) ||
                 (a_location == FaceStop::HiWithBoundary) ||
                 (a_location == FaceStop::AllBoundaryOnly) ||
                 (a_location == FaceStop::HiBoundaryOnly) ||
                 (a_location == FaceStop::HiNoBoundary));
    bool doBoundary = ((a_location == FaceStop::HiWithBoundary) ||
                       (a_location == FaceStop::SurroundingWithBoundary) ||
                       (a_location == FaceStop::LoWithBoundary) ||
                       (a_location == FaceStop::LoBoundaryOnly) ||
                       (a_location == FaceStop::HiBoundaryOnly) ||
                       (a_location == FaceStop::AllBoundaryOnly));
    bool doBoundaryOnly = ((a_location == FaceStop::AllBoundaryOnly) ||
                           (a_location == FaceStop::LoBoundaryOnly) ||
                           (a_location == FaceStop::HiBoundaryOnly));
    bool doInterior = !doBoundaryOnly;

    //if this fails, invalid location.
    BL_ASSERT(doLo || doHi || doBoundaryOnly);

    Side::LoHiSide sides[2] =
      {
        doLo ? Side::Lo : Side::Invalid,
        doHi ? Side::Hi : Side::Invalid
      };

    for (IVSIterator ivsit(a_ivs); ivsit.ok(); ++ivsit)
    {
      Array<VolIndex> vofs = a_ebgraph.getVoFs(ivsit());
      for (int ivof=0; ivof<vofs.size(); ++ivof)
      {
        for (int iside=0; iside<2; ++iside )
        {
          if ( sides[iside] == Side::Invalid )
          {
            continue;
          }
          Array<FaceIndex> faces(
            a_ebgraph.getFaces( vofs[ivof], a_direction, sides[iside] ) );

          for (int iface=0; iface<faces.size(); ++iface)
          {
            const FaceIndex& face =  faces[iface];
            if (   (  face.isBoundary()  && doBoundary)
                   || ((!face.isBoundary()) && doInterior) )
            {
              resultSet.insert( face );
            }
          }
        }
      }
    }

    m_faces.clear();
    m_faces.reserve( resultSet.size() );
    for ( std::set<FaceIndex>::const_iterator i = resultSet.begin();
          i != resultSet.end();
          ++i )
    {
      m_faces.push_back( *i );
    }
    reset();
  }

/********************************/
  void
  FaceIterator::reset()
  {
    assert(isDefined());
    m_iface = 0;
  }

/********************************/
  void
  FaceIterator::operator++()
  {
    assert(isDefined());
    m_iface++;
  }

/********************************/
  const FaceIndex&
  FaceIterator::operator() () const
  {
    assert(isDefined());
    assert(m_iface < m_faces.size());
    return m_faces[m_iface];
  }

/********************************/
  bool
  FaceIterator::ok() const
  {
    return (m_iface < m_faces.size());
  }

/********************************/
  bool
  FaceIterator::isDefined() const
  {
    return m_isDefined;
  }

}

