#ifndef _BASEEBFACEFABI_H_
#define _BASEEBFACEFABI_H_

#include "AMReX_BoxIterator.H"

namespace amrex
{

  /**********************/
  template <class T>
  void
  BaseEBFaceFAB<T>::setCoveredFaceVal(const T&    a_val,
                                      const int&  a_comp,
                                      const bool& a_doMulti)
  {
    BL_ASSERT(a_comp >= 0);
    BL_ASSERT(a_comp < m_nComp);

    if (m_ebisBox.isAllRegular())
    {
      return;
    }
    else if (m_ebisBox.isAllCovered())
    {
      m_regFAB.setVal(a_val, a_comp);
    }
    else
    {
      for (BoxIterator bit(m_region); bit.ok(); ++bit)
      {
        const IntVect& ivCell = bit();
        for (SideIterator sit; sit.ok(); ++sit)
        {
          const Vector< FaceIndex >& faces = m_ebisBox.getAllFaces(ivCell,m_iDir,sit());
          if ( (faces.size()==0) || (faces.size()>1) )
          {
            const int& hilo =  sign(sit());
            IntVect ivFace = ivCell;
            if (hilo == 1)
            {
              ivFace[m_iDir] += 1;
            }
            m_regFAB(ivFace,a_comp) = a_val;
          }
        }
      }
    }
  }

  /*************/
  template <class T>
  std::size_t 
  BaseEBFaceFAB<T>::
  nBytes (const Box& bx, int start_comp, int ncomps) const
  {
    BL_PROFILE("BaseEBCellFAB::size");
    size_t retval =  m_regFAB.nBytes(bx, start_comp, ncomps);
    retval += m_irrFAB.nBytes(bx, start_comp, ncomps);
    return retval;
  }
      
  /*************/
  template <class T>
  std::size_t 
  BaseEBFaceFAB<T>::
  copyToMem (const Box& srcbox,
             int        srccomp,
             int        numcomp,
             void*      dst) const
  {
    BL_PROFILE("BaseEBCellFAB::linearOut");
    unsigned char* buffer = (unsigned char*) dst;
    size_t retval = 0;
    retval = m_regFAB.copyToMem(srcbox, srccomp, numcomp, buffer);
    buffer += retval;
    retval += m_irrFAB.copyToMem(srcbox, srccomp, numcomp, buffer);
    return retval;
  }
  /*************/
  template <class T>
  std::size_t 
  BaseEBFaceFAB<T>::
  copyFromMem (const Box&  dstbox,
               int         dstcomp,
               int         numcomp,
               const void* src)
  {
    BL_PROFILE("BaseEBCellFAB::linearIn");
    unsigned char* buffer = (unsigned char*) src;
    size_t retval = 0;
    retval = m_regFAB.copyFromMem(dstbox, dstcomp, numcomp, buffer);
    buffer += retval;
    retval += m_irrFAB.copyFromMem(dstbox, dstcomp, numcomp, buffer);
    return retval;
  }

  /*************/
  template <class T> inline
  const EBISBox&
  BaseEBFaceFAB<T>::getEBISBox() const
  {
    return m_ebisBox;
  }
  /*************/
  template <class T>
  BaseEBFaceFAB<T>&
  BaseEBFaceFAB<T>::
  copy(const BaseEBFaceFAB<T>& a_src,
       const Box&              a_srcbox,
       int                     a_srccomp,
       const Box&              a_dstbox,
       int                     a_dstcomp,
       int                     a_numcomp)
  {
    Box faceBox = surroundingNodes(a_dstbox, m_iDir);
    
    m_regFAB.copy(a_src.m_regFAB, faceBox , a_srccomp, faceBox , a_dstcomp, a_numcomp);
    m_irrFAB.copy(a_src.m_irrFAB, a_srcbox, a_srccomp, a_dstbox, a_dstcomp, a_numcomp);

    return *this;
  }
  /**********************/
  /**********************/
  template <class T> inline
  BaseEBFaceFAB<T>::BaseEBFaceFAB()
  {
    setDefaultValues();
  }

  /**********************/
  /**********************/
  template <class T> inline
  BaseEBFaceFAB<T>::BaseEBFaceFAB(const EBISBox& a_ebisBox,
                                  const Box& a_region,
                                  int a_iDir, int a_nComp)
  {
    setDefaultValues();
    define(a_ebisBox, a_region, a_iDir, a_nComp);
  }

  /**********************/
  /**********************/
  template <class T> inline
  void
  BaseEBFaceFAB<T>::define(const EBISBox&  a_ebisBox,
                           const Box& a_region,
                           int a_iDir, int a_nComp)
  {
    clear();
    m_isDefined = true;
    BL_ASSERT(a_region.cellCentered());
    BL_ASSERT(!a_region.isEmpty());
    BL_ASSERT(a_nComp > 0);
    BL_ASSERT((a_iDir >= 0) && (a_iDir < SpaceDim));

    m_ebisBox = a_ebisBox;
    m_region =  a_region;
    // GHM don't mask region
    m_region &= a_ebisBox.getDomain();
    BL_ASSERT(!m_region.isEmpty());

    m_nComp = a_nComp;
    m_iDir = a_iDir;
    m_regionFace= surroundingNodes(a_region, a_iDir);
    Box grownRegion = a_region;
    grownRegion.grow(a_iDir, 1);

    Box checkRegion = grownRegion;
    checkRegion &= a_ebisBox.getDomain();
    //if this fails, you need an ebisbox with more ghost cells
    BL_ASSERT(a_ebisBox.getRegion().contains(checkRegion));

    m_regFAB.resize(m_regionFace, m_nComp);

    IntVectSet ivsmulti = m_ebisBox.getMultiCells(a_region);
    m_irrFAB.define(ivsmulti, m_ebisBox.getEBGraph(), m_iDir, m_nComp);
  }

  /**********************/
  /**********************/
  template <class T> inline
  BaseEBFaceFAB<T>::~BaseEBFaceFAB()
  {
    clear();
  }

  /**********************/
  /**********************/
  template <class T> inline
  void
  BaseEBFaceFAB<T>::clear()
  {
    m_irrFAB.clear();
    m_regFAB.clear();
    m_isDefined = false;
  }

  /**********************/
  /**********************/
  template <class T> inline
  void
  BaseEBFaceFAB<T>::setVal(const T& value)
  {
    m_irrFAB.setVal(value);
    m_regFAB.setVal(value);
  }

  /**********************/
  template <class T> inline
  void
  BaseEBFaceFAB<T>::setVal(int ivar, const T& value)
  {
    m_irrFAB.setVal(ivar,  value);
    m_regFAB.setVal(value, m_regFAB.box(), ivar, 1);
  }

  /**********************/
  /**********************/
  template <class T> inline
  bool
  BaseEBFaceFAB<T>::isDefined() const
  {
    return (m_isDefined);
  }

  /**********************/
  /**********************/
  template <class T> inline
  int
  BaseEBFaceFAB<T>::nComp() const
  {
    BL_ASSERT(isDefined());
    return m_nComp;
  }

  /**********************/
  /**********************/
  template <class T> inline
  const BaseIFFAB<T>&
  BaseEBFaceFAB<T>::getMultiValuedFAB() const
  {
    BL_ASSERT(isDefined());
    return m_irrFAB;
  }

  /**********************/
  /**********************/
  template <class T> inline
  BaseIFFAB<T>&
  BaseEBFaceFAB<T>::getMultiValuedFAB()
  {
    BL_ASSERT(isDefined());
    return m_irrFAB;
  }

  /**********************/
  /**********************/
  template <class T> inline
  const BaseFab<T>&
  BaseEBFaceFAB<T>::getSingleValuedFAB() const
  {
    BL_ASSERT(isDefined());
    return m_regFAB;
  }

  /**********************/
  /**********************/
  template <class T> inline
  BaseFab<T>&
  BaseEBFaceFAB<T>::getSingleValuedFAB()
  {
    BL_ASSERT(isDefined());
    return m_regFAB;
  }

  /**********************/
  /**********************/
  template <class T> inline
  int
  BaseEBFaceFAB<T>::direction() const
  {
    BL_ASSERT(isDefined());
    return m_iDir;
  }

  /**********************/
  /**********************/
  template <class T> inline
  const Box&
  BaseEBFaceFAB<T>::getRegion() const
  {
    BL_ASSERT(isDefined());
    return m_regionFace;
  }
  /**********************/
  /**********************/
  template <class T> inline
  const Box&
  BaseEBFaceFAB<T>::getCellRegion() const
  {
    BL_ASSERT(isDefined());
    return m_region;
  }

  /**********************/
  /**********************/
  template <class T> inline
  T&
  BaseEBFaceFAB<T>::operator()(const FaceIndex& a_facein, int  a_nCompLoc)
  {
    BL_ASSERT(isDefined());
    BL_ASSERT((a_nCompLoc >= 0)&&(a_nCompLoc < m_nComp));
    const IntVect& ivhi = a_facein.gridIndex(Side::Hi);
    BL_ASSERT(m_regionFace.contains(ivhi));

    const Vector<FaceIndex>& multiFaces = m_irrFAB.getFaces();
    T* returnval;

    bool isFaceHere = false;
    for (int iface = 0; iface< multiFaces.size(); iface++)
    {
      const FaceIndex& face = multiFaces[iface];
      if (face == a_facein)
      {
        isFaceHere = true;
      }
    }
    if (isFaceHere)
    {
      returnval = &m_irrFAB(a_facein, a_nCompLoc);
    }
    else
    {
      returnval = &m_regFAB(ivhi, a_nCompLoc);
    }
    return *returnval;
  }

  /**********************/
  /**********************/
  template <class T> inline
  const T&
  BaseEBFaceFAB<T>::operator() (const FaceIndex& a_facein, int  a_nCompLoc) const
  {
    BL_ASSERT(isDefined());
    BL_ASSERT((a_nCompLoc >= 0)&&(a_nCompLoc < m_nComp));
    const IntVect& ivhi = a_facein.gridIndex(Side::Hi);
    BL_ASSERT(m_regionFace.contains(ivhi));

    const Vector<FaceIndex>& multiFaces = m_irrFAB.getFaces();
    const T* returnval;

    bool isFaceHere = false;
    for (int iface = 0; iface< multiFaces.size(); iface++)
    {
      const FaceIndex& face = multiFaces[iface];
      if (face == a_facein)
      {
        isFaceHere = true;
      }
    }
    if (isFaceHere)
    {
      returnval = &m_irrFAB(a_facein, a_nCompLoc);
    }
    else
    {
      returnval = &m_regFAB(ivhi, a_nCompLoc);
    }
    return *returnval;
  }

  /**********************/
  /**********************/
  template <class T> inline
  void
  BaseEBFaceFAB<T>::setDefaultValues()
  {
    m_isDefined = false;
    m_nComp = -1;
    m_iDir = -1;
  }
}
#endif
