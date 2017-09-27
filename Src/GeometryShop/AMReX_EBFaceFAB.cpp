
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

#include "AMReX_EBFaceFAB.H"
#include "AMReX_BoxIterator.H"

namespace amrex
{
  /**********************/
  /**********************/
  EBFaceFAB::EBFaceFAB():BaseEBFaceFAB<Real>()
  {
  }
               
  /**********************/
  /**********************/
  EBFaceFAB::EBFaceFAB(const EBISBox& a_ebisBox,
                       const Box& a_region,
                       int a_iDir, int a_nComp)
    :BaseEBFaceFAB<Real>(a_ebisBox, a_region, a_iDir, a_nComp)
  {
  }
  void
  EBFaceFAB::define(const EBISBox&  a_ebisBox,
                    const Box& a_region,
                    int a_iDir, int a_nComp)
  {
    BaseEBFaceFAB<Real>::define(a_ebisBox, a_region, a_iDir, a_nComp);
  }
               
  /**********************/
  /**********************/
  EBFaceFAB::~EBFaceFAB()
  {
  }
               
  /**********************/
  /**********************/
  const FArrayBox&
  EBFaceFAB::getFArrayBox() const
  {
    BL_ASSERT(isDefined());
    return (const FArrayBox&)m_regFAB;
  }
  /**********************/
  /**********************/
  FArrayBox&
  EBFaceFAB::getFArrayBox()
  {
    BL_ASSERT(isDefined());
    return (FArrayBox&)m_regFAB;
  }
               
  /**********************/
  /**********************/
  EBFaceFAB&
  EBFaceFAB::operator+=(const EBFaceFAB& a_src)
  {
    BL_ASSERT(a_src.m_nComp == m_nComp);
    plus(a_src, 0, 0, m_nComp);
               
    return *this;
  }
               
  EBFaceFAB&
  EBFaceFAB::plus(const EBFaceFAB& a_src,
                  int  a_srccomp,
                  int  a_dstcomp,
                  int  a_numcomp,
                  Real a_scale)
  {
    int numcomp = a_numcomp;
    if(numcomp < 0)
    {
      numcomp = nComp();
    }

    BL_ASSERT(isDefined());
    BL_ASSERT(a_src.isDefined());
    //  BL_ASSERT(m_ebisBox == a_src.m_ebisBox);
    BL_ASSERT(a_srccomp + numcomp <= a_src.m_nComp);
    BL_ASSERT(a_dstcomp + numcomp <= m_nComp);
               
    Box locRegionFace = a_src.m_regionFace & m_regionFace;
    //Box locRegion = a_src.m_region & m_region;
    if (!locRegionFace.isEmpty())
    {
      for(BoxIterator bit(locRegionFace); bit.ok(); ++bit)
      {
        for(int icomp = 0; icomp < numcomp; icomp++)
        {
          int isrc = a_srccomp + icomp;
          int idst = a_dstcomp + icomp;
          m_regFAB(bit(), idst) += a_scale*a_src.m_regFAB(bit(), isrc);
        }
      }
               
      const Array<FaceIndex>& faces = m_irrFAB.getFaces();
      for (int iface = 0; iface < faces.size(); iface++)
      {
        const FaceIndex& face = faces[iface];
        if (locRegionFace.contains(face.gridIndex(Side::Hi)))
        {
          for (int icomp = 0; icomp < numcomp; ++icomp)
          {
            m_irrFAB(face, a_dstcomp+icomp) +=
              a_scale*a_src.m_irrFAB(face, a_srccomp+icomp);
          }
        }
      }
    }
    return *this;
  }
               
  /**********************/
  /**********************/
  EBFaceFAB&
  EBFaceFAB::operator-=(const EBFaceFAB& a_src)
  {
    BL_ASSERT(a_src.m_nComp == m_nComp);
               
    minus(a_src, 0, 0, m_nComp);
               
    return *this;
  }
               
  EBFaceFAB&
  EBFaceFAB::minus(const EBFaceFAB& a_src,
                   int  a_srccomp,
                   int  a_dstcomp,
                   int  a_numcomp,
                   Real a_scale)
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(a_src.isDefined());
    int numcomp = a_numcomp;
    if(numcomp < 0)
    {
      numcomp = nComp();
    }
    // (DFM 7/28/05) This assertion fails for multifluid cases
    // where this and src have been defined using different
    // factories; what we really need is a better implementation
    // of the operator== for EBISBox
    //BL_ASSERT(m_ebisBox == a_src.m_ebisBox);
    BL_ASSERT(a_srccomp + numcomp <= a_src.m_nComp);
    BL_ASSERT(a_dstcomp + numcomp <= m_nComp);
               
    Box locRegionFace = a_src.m_regionFace & m_regionFace;
    //Box locRegion = a_src.m_region & m_region;
    if (!locRegionFace.isEmpty())
    {
      for(BoxIterator bit(locRegionFace); bit.ok(); ++bit)
      {
        for(int icomp = 0; icomp < numcomp; icomp++)
        {
          int isrc = a_srccomp + icomp;
          int idst = a_dstcomp + icomp;
          m_regFAB(bit(), idst) -= a_scale*a_src.m_regFAB(bit(), isrc);
        }
      }
               
      const Array<FaceIndex>& faces = m_irrFAB.getFaces();
      for (int iface = 0; iface < faces.size(); iface++)
      {
        const FaceIndex& face = faces[iface];
        if (locRegionFace.contains(face.gridIndex(Side::Hi)))
        {
          for (int icomp = 0; icomp < numcomp; ++icomp)
          {
            m_irrFAB(face, a_dstcomp+icomp) -=
              a_scale*a_src.m_irrFAB(face, a_srccomp+icomp);
          }
        }
      }
    }
    return *this;
  }
               
  /**********************/
  /**********************/
  EBFaceFAB&
  EBFaceFAB::operator*=(const EBFaceFAB& a_src)
  {
    BL_ASSERT(a_src.m_nComp == m_nComp);
               
    mult(a_src, 0,0,m_nComp);
               
    return *this;
  }
               
               
               
  EBFaceFAB&
  EBFaceFAB::mult(const EBFaceFAB& a_src,
                  int a_srccomp,
                  int a_dstcomp,
                  int a_numcomp)
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(a_src.isDefined());
    BL_ASSERT(m_ebisBox == a_src.m_ebisBox);
    int numcomp = a_numcomp;
    BL_ASSERT(a_srccomp + numcomp <= a_src.m_nComp);
    BL_ASSERT(a_dstcomp + numcomp <= m_nComp);
    if(numcomp < 0)
    {
      numcomp = nComp();
    }
               
    Box locRegionFace = a_src.m_regionFace & m_regionFace;
    //Box locRegion = a_src.m_region & m_region;
    if (!locRegionFace.isEmpty())
    {
      for(BoxIterator bit(locRegionFace); bit.ok(); ++bit)
      {
        for(int icomp = 0; icomp < numcomp; icomp++)
        {
          int isrc = a_srccomp + icomp;
          int idst = a_dstcomp + icomp;
          m_regFAB(bit(), idst) *= a_src.m_regFAB(bit(), isrc);
        }
      }
               
      const Array<FaceIndex>& faces = m_irrFAB.getFaces();
      for (int iface = 0; iface < faces.size(); iface++)
      {
        const FaceIndex& face = faces[iface];
        if (locRegionFace.contains(face.gridIndex(Side::Hi)))
        {
          for (int icomp = 0; icomp < numcomp; ++icomp)
          {
            m_irrFAB(face, a_dstcomp+icomp) *=
              a_src.m_irrFAB(face, a_srccomp+icomp);
          }
        }
      }
    }
    return *this;
  }
               
  /**********************/
  EBFaceFAB&
  EBFaceFAB::operator/=(const EBFaceFAB& a_src)
  {
    BL_ASSERT(a_src.m_nComp == m_nComp);
               
    divide(a_src, 0, 0, m_nComp);
               
    return *this;
  }
  /**********************/
  EBFaceFAB&
  EBFaceFAB::divide(const EBFaceFAB& a_src,
                    int a_srccomp,
                    int a_dstcomp,
                    int a_numcomp)
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(a_src.isDefined());
    BL_ASSERT(m_ebisBox == a_src.m_ebisBox);
    int numcomp = a_numcomp;
    BL_ASSERT(a_srccomp + numcomp <= a_src.m_nComp);
    BL_ASSERT(a_dstcomp + numcomp <= m_nComp);
    if(numcomp < 0)
    {
      numcomp = nComp();
    }
               
    Box locRegionFace = a_src.m_regionFace & m_regionFace;
    //Box locRegion = a_src.m_region & m_region;
    if (!locRegionFace.isEmpty())
    {
      for(BoxIterator bit(locRegionFace); bit.ok(); ++bit)
      {
        for(int icomp = 0; icomp < numcomp; icomp++)
        {
          int isrc = a_srccomp + icomp;
          int idst = a_dstcomp + icomp;
          m_regFAB(bit(), idst) /= a_src.m_regFAB(bit(), isrc);
        }
      }

      const Array<FaceIndex>& faces = m_irrFAB.getFaces();
      for (int iface = 0; iface < faces.size(); iface++)
      {
        const FaceIndex& face = faces[iface];
        if (locRegionFace.contains(face.gridIndex(Side::Hi)))
        {
          for (int icomp = 0; icomp < numcomp; ++icomp)
          {
            m_irrFAB(face, a_dstcomp+icomp) /=
              a_src.m_irrFAB(face, a_srccomp+icomp);
          }
        }
      }
    }
    return *this;
  }
  /**********************/
  EBFaceFAB&
  EBFaceFAB::operator+=(const Real& a_src)
  {
    BL_ASSERT(isDefined());
    for(BoxIterator bit(m_region); bit.ok(); ++bit)
    {
      for(int icomp = 0; icomp < nComp(); icomp++)
      {
        m_regFAB(bit(), icomp) += a_src;
      }
    }
               
    const Array<FaceIndex>& faces = m_irrFAB.getFaces();
    for (int iface = 0; iface < faces.size(); iface++)
    {
      const FaceIndex& face = faces[iface];
      for (int icomp = 0; icomp < m_nComp; ++icomp)
      {
        m_irrFAB(face, icomp) += a_src;
      }
    }
    return *this;
  }
               
  /**********************/
  EBFaceFAB&
  EBFaceFAB::operator*=(const Real& a_src)
  {
    BL_ASSERT(isDefined());
               
    for(BoxIterator bit(m_region); bit.ok(); ++bit)
    {
      for(int icomp = 0; icomp < nComp(); icomp++)
      {
        m_regFAB(bit(), icomp) *= a_src;
      }
    }
               
    const Array<FaceIndex>& faces = m_irrFAB.getFaces();
    for (int iface = 0; iface < faces.size(); iface++)
    {
      const FaceIndex& face = faces[iface];
      for (int icomp = 0; icomp < m_nComp; ++icomp)
      {
        m_irrFAB(face, icomp) *= a_src;
      }
    }
    return *this;
  }
  /**********************/
  Real
  EBFaceFAB::max(int a_comp) const
  {
    BL_ASSERT(isDefined());
    Real val = -1.0e30;
               
    for(BoxIterator bit(m_region); bit.ok(); ++bit)
    {
      val = std::max(val, m_regFAB(bit(), a_comp));
    }
               
    // Find the max on irregular faces.
    const Array<FaceIndex>& faces = m_irrFAB.getFaces();
    for (int iface = 0; iface < faces.size(); iface++)
    {
      const FaceIndex& face = faces[iface];
      for (int icomp = 0; icomp < m_nComp; ++icomp)
      {
        val = std::max(val, m_irrFAB(face, icomp));
      }
    }
    return val;
  }
  //-----------------------------------------------------------------------
  Real
  EBFaceFAB::min(int a_comp) const
  {
    BL_ASSERT(isDefined());
    Real val = 1.0e30;
               
    for(BoxIterator bit(m_region); bit.ok(); ++bit)
    {
      val = std::min(val, m_regFAB(bit(), a_comp));
    }
               
    // Find the max on irregular faces.
    const Array<FaceIndex>& faces = m_irrFAB.getFaces();
    for (int iface = 0; iface < faces.size(); iface++)
    {
      const FaceIndex& face = faces[iface];
      for (int icomp = 0; icomp < m_nComp; ++icomp)
      {
        val = std::min(val, m_irrFAB(face, icomp));
      }
    }
    return val;
  }
  //-----------------------------------------------------------------------
               
}
