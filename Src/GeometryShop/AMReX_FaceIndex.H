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

#ifndef _FACEINDEX_H_
#define _FACEINDEX_H_

#include "AMReX_LoHiSide.H"
#include "AMReX_VolIndex.H"
#include "AMReX_BaseIndex.H"
#include "AMReX_SPMD.H"
using std::ostream;
namespace amrex
{
  ///
  /**
     FaceIndex is a very lightweight object used
     to distinguish faces.  It has two VolIndex (s)
     and a direction and a way to access them.
  */
  class FaceIndex : public BaseIndex
  {
  public:
    ///
    /**
       if both cellindexes of the vofs are >= 0,
       then the face is interior.  otherwise the face is
       a boundary face.
    */
    FaceIndex(const VolIndex& a_vof1,
              const VolIndex& a_vof2,
              const int& a_direction);
    ///
    /**
       shorthand constructor.
    */
    FaceIndex(const VolIndex& a_vof1,
              const VolIndex& a_vof2);

    ///
    /**
       if both cellindexes of the vofs are >= 0,
       then the face is interior.  otherwise the face is
       a boundary face.
    */
    void define(const VolIndex& a_vof1,
                const VolIndex& a_vof2,
                const int& a_direction);

    ///
    /**
       Shorthand define. First determine if they are neighbors
       then find the direction to pass to full define
    */
    void define(const VolIndex& a_vof1,
                const VolIndex& a_vof2);
    ///
    /**
     */
    FaceIndex();

    ///
    /**
     */
    FaceIndex(const FaceIndex& a_facein);

    ///
    /**
     */
    void define(const FaceIndex& a_facein);

    ///
    /**
     */
    ~FaceIndex();

    ///
    /**
     */
    FaceIndex& operator= (const FaceIndex& a_facein);

    ///
    /**
     */
    bool operator== (const FaceIndex& a_facein) const {
	return ((m_loIndex == a_facein.m_loIndex) &&
		(m_hiIndex == a_facein.m_hiIndex) &&
		(m_loiv == a_facein.m_loiv) &&
		(m_hiiv == a_facein.m_hiiv));
    }

    ///
    /**
     */
    bool operator!= (const FaceIndex& a_facein) const;

    ///
    /**
     */
    const int& direction() const;

    ///
    /**
     */
    const bool& isDefined() const;

    ///
    /**
     */
    const bool& isBoundary() const;

    ///
    /**
       return cell index of the vof on the a_sd side.
       returns -1 if that side of the face is outside
       the domain (boundary face case).
    */
    const int& cellIndex(const Side::LoHiSide& a_sd) const;

    ///
    /**
       return grid index of the vof on the a_sd side.
       may be inside or ouside domain.
    */
    const IntVect& gridIndex(const Side::LoHiSide& a_sd) const;

    ///
    /**
       manufactures the appropriate volindex.
    */
    VolIndex getVoF(const Side::LoHiSide& a_sd) const;

    ///
    /**
       a_vof needs to be one of the two vofs connected by this face
       Standing on a_vof, face opens toward the other vof
       returns 1 if face opens toward positive axis
       returns -1 if face opens toward negative axis
       returns 0 if a_vof is not one of the member vofs.
    */
    int faceSign(const VolIndex& a_vof) const;

    ///for stencil gymnastics
    void shift(const IntVect& a_iv)
    {
      m_loiv += a_iv;
      m_hiiv += a_iv;
    }

    ///
    /** Lexicographic */
    friend bool operator<( const FaceIndex& f1, const FaceIndex& f2 );

    friend ostream& operator<<( ostream& out, const FaceIndex& fi );

  private:

    ///the intvect of the low vof
    /**
     */
    IntVect m_loiv;

    ///the intvect of the high vof
    /**
     */
    IntVect m_hiiv;

    ///
    /**
       The index of the low vof.
       This = -1 if the face is set to be a boundary
       face on the low side of the domain.
    */
    int m_loIndex;

    ///
    /**
       The index of the high vof.
       This = -1 if the face is set to be a boundary
       face on the high side of the domain.
    */
    int m_hiIndex;

    ///direction of the face.
    int m_direction;

    ///
    /**
       true if one of the define functions has been called.
    */
    bool m_isDefined;

    ///
    /**
       true if the face is constructed
       using the boundary face constructor
    */
    bool m_isBoundary;

  };

  template < >
  int linearSize(const FaceIndex& findex);

//FaceIndex specialization of linearIn
  template < >
  void linearIn(FaceIndex& a_outputT, const void* const inBuf);

//FaceIndex specialization of linearOut
  template < >
  void linearOut(void* const a_outBuf, const FaceIndex& a_inputT);



/** Friend of class FaceIndex */
inline
bool operator<( const FaceIndex& f1, const FaceIndex& f2 )
{
  if ( f1.m_loiv<( f2.m_loiv ) )
    {
      return true;
    }
  else
  if ( f2.m_loiv<( f1.m_loiv ) )
    {
      return false;
    }
  else
  if ( f1.m_loIndex < f2.m_loIndex )
    {
      return true;
    }
  else
  if ( f2.m_loIndex < f1.m_loIndex )
    {
      return false;
    }
  else
  if ( f1.m_hiiv<( f2.m_hiiv ) )
    {
      return true;
    }
  else
  if ( f2.m_hiiv<( f1.m_hiiv ) )
    {
      return false;
    }
  else
  if ( f1.m_hiIndex < f2.m_hiIndex )
    {
      return true;
    }
  else
  if ( f2.m_hiIndex < f1.m_hiIndex )
    {
      return false;
    }
  else
  if ( f1.m_direction < f2.m_direction )
    {
      return true;
    }
  else
  if ( f2.m_direction < f1.m_direction )
    {
      return false;
    }
  else
  if ( f1.m_isBoundary && (!f2.m_isBoundary) )
    {
      return true;
    }
  else
  if ( (!f1.m_isBoundary) && f2.m_isBoundary )
    {
      return false;
    }

  return false;
}
}
#endif
