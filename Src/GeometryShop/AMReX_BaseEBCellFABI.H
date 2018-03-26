#ifndef _BASEEBCELLFABI_H_
#define _BASEEBCELLFABI_H_

#include "AMReX_BoxIterator.H"

namespace amrex
{
/**********************/
  template <class T>
  void
  BaseEBCellFAB<T>::setCoveredCellVal(const T&    a_val,
                                      const int&  a_comp)
  {
    BL_PROFILE("BaseEBCellFAB::setCoveredCellVal");
    BL_ASSERT(a_comp >= 0);
    BL_ASSERT(a_comp < m_regFAB.nComp());

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
      for (BoxIterator bit(getRegion()); bit.ok(); ++bit)
      {
        const IntVect& iv = bit();
        if (m_ebisBox.isCovered(iv))
        {
          m_regFAB(iv, a_comp) = a_val;
        }
      }
      const Vector<VolIndex>& vofs = m_irrFAB.getVoFs(); 
      for (int i=0;  i< vofs.size(); i++)
      {
        const IntVect& iv = vofs[i].gridIndex();
        m_regFAB(iv, a_comp) = a_val;
      }
    }
  }
  /*************/
  template <class T>
  BaseEBCellFAB<T>&
  BaseEBCellFAB<T>::
  copy(const BaseEBCellFAB<T>& a_src)
  {

    copy(a_src, a_src.m_region, 0, m_region, 0, nComp());
    return *this;
  }
      
  /*************/
  /*************/
  template <class T>
  BaseEBCellFAB<T>&
  BaseEBCellFAB<T>::
  copy(const BaseEBCellFAB<T>& a_src,
       const Box&              a_srcbox,
       int                     a_srccomp,
       const Box&              a_dstbox,
       int                     a_dstcomp,
       int                     a_numcomp)
  {
    m_regFAB.copy(a_src.m_regFAB, a_srcbox, a_srccomp, a_dstbox, a_dstcomp, a_numcomp);
    if (m_hasMultiCells)
    {
      m_irrFAB.copy(a_src.m_irrFAB, a_srcbox, a_srccomp, a_dstbox, a_dstcomp, a_numcomp);
    }
    return *this;
  }
      
      
  /*************/
  /*************/
  template <class T>
  Vector<VolIndex>
  BaseEBCellFAB<T>::getMultiCells() const
  {
    Vector<VolIndex> retval;
    BL_ASSERT(isDefined());
    if(m_irrFAB.hasVoFs())
    {
      retval =  m_irrFAB.getVoFs();
    }
    return retval;
  }
      
  /*************/
  /*************/
  template <class T>
  BaseEBCellFAB<T>::BaseEBCellFAB()
  {
    setDefaultValues();
  }
      
  /*************/
  /*************/
  template <class T>
  BaseEBCellFAB<T>::BaseEBCellFAB(const EBISBox& a_ebisBox,
                                  const Box& a_region, int a_nVar)
  {
    setDefaultValues();
    define(a_ebisBox, a_region, a_nVar);
  }
      
  /*************/
  /*************/
  template <class T>
  void
  BaseEBCellFAB<T>::define(const EBISBox& a_ebisBox,
                           const Box& a_region, int a_nVar)
  {
    BL_PROFILE("BaseEBCellFAB::define");
    clear();
    m_isDefined = true;
    BL_ASSERT(a_region.cellCentered());
    BL_ASSERT(a_nVar > 0);
    BL_ASSERT(!a_region.isEmpty());
      
    m_region = a_region & a_ebisBox.getRegion();
    m_ebisBox = a_ebisBox;
      
    //make reg fab the full size of the box.
    //this makes ebalias possible
    m_regFAB.resize(a_region, a_nVar);
    //add all the multi-valued
    //cells to m_multiValued
    IntVectSet multiCells = m_ebisBox.getMultiCells(getRegion());
      
    m_hasMultiCells = !multiCells.isEmpty();
      
    //define irregular fab to be over multi-valued cells.
    m_irrFAB.define(multiCells, m_ebisBox.getEBGraph(), a_nVar);
  }
      
  /*************/
  /*************/
  template <class T>
  BaseEBCellFAB<T>::
  ~BaseEBCellFAB()
  {
    clear();
  }
      
  /*************/
  /*************/
  template <class T>
  void
  BaseEBCellFAB<T>::clear()
  {
    //m_irrFAB.clear();
    m_regFAB.clear();
    // m_multiCells.makeEmpty();
      
    m_isDefined = false;
  }
  /*************/
  /*************/
  template <class T>
  void
  BaseEBCellFAB<T>::setVal(const T& value)
  {
    m_irrFAB.setVal(value);
    m_regFAB.setVal(value);
  }
      
  /*************/
  template <class T>
  void
  BaseEBCellFAB<T>::setVal(int ivar,const  T& value)
  {
    m_irrFAB.setVal(ivar, value);
    m_regFAB.setVal(value, m_regFAB.box(), ivar, 1);
  }
  /*************/
  template <class T>
  void
  BaseEBCellFAB<T>::setVal(const T&   a_value,
                           const Box& a_box,
                           int        a_nstart,
                           int        a_numcomp)
  {
    BL_ASSERT(isDefined());
    if (!a_box.isEmpty())
    {
      m_irrFAB.setVal(a_value,a_box,a_nstart,a_numcomp);
      m_regFAB.setVal(a_value,a_box,a_nstart,a_numcomp);
    }
  }
      
  /*************/
  /*************/
  template <class T>
  bool
  BaseEBCellFAB<T>::isDefined() const
  {
    return (m_isDefined);
  }
  /*************/
  /*************/
  template <class T>
  int
  BaseEBCellFAB<T>::nComp() const
  {
    BL_ASSERT(isDefined());
    return m_regFAB.nComp();
  }
  /*************/
  /*************/
  template <class T>
  const Box&
  BaseEBCellFAB<T>::getRegion() const
  {
    BL_ASSERT(isDefined());
    return m_region;
  }
  /*************/
  /*************/
  template <class T>
  const Box&
  BaseEBCellFAB<T>::box() const
  {
    return getRegion();
  }
      
  /*************/
  /*************/
  template <class T>
  const BaseIVFAB<T>&
  BaseEBCellFAB<T>::getMultiValuedFAB() const
  {
    BL_ASSERT(isDefined());
    return m_irrFAB;
  }
      
  //
  /*************/
  /*************/
  template <class T>
  BaseIVFAB<T>&
  BaseEBCellFAB<T>::getMultiValuedFAB()
  {
    BL_ASSERT(isDefined());
    return m_irrFAB;
  }
      
  /*************/
  /*************/
  template <class T>
  const BaseFab<T>&
  BaseEBCellFAB<T>::getSingleValuedFAB() const
  {
    BL_ASSERT(isDefined());
    return m_regFAB;
  }
      
  //
  /*************/
  /*************/
  template <class T>
  BaseFab<T>&
  BaseEBCellFAB<T>::getSingleValuedFAB()
  {
    BL_ASSERT(isDefined());
    return m_regFAB;
  }
  /*************/
      
  /*************
        Arg a_isKnownMultiValued should be set to:
           -1 if arg a_ndin is known to be in a single-valued cell;
           +1 if arg a_ndin is known to be in a multi-valued cell;
            0 if you're not sure (in which case the expensive IntVectSet::contains()
              function gets called).
  */
  template <class T>
  const T&
  BaseEBCellFAB<T>::operator()(const VolIndex& a_ndin,int  a_nVarLoc,
                               int a_isKnownMultiValued ) const
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(!m_ebisBox.isCovered(a_ndin.gridIndex()));
    //BL_ASSERT(getRegion().contains(a_ndin.gridIndex())); caught in m_regFAB
    BL_ASSERT((a_nVarLoc >= 0)&&(a_nVarLoc < nComp()));
    BL_ASSERT( (a_isKnownMultiValued == -1)
               ||    (a_isKnownMultiValued ==  0)
               ||    (a_isKnownMultiValued == +1) );
      
    const T* returnval = NULL;
      
    switch( a_isKnownMultiValued )
    {
    case -1:
      returnval = &(m_regFAB(a_ndin.gridIndex(), a_nVarLoc));
      break;        // m_regFAB is a BaseFab
    case +1:
      returnval = &(m_irrFAB(a_ndin,a_nVarLoc));
      break;        // m_irrFAB is a BaseIVFAB
    case 0:
    {
      if ( m_ebisBox.isMultiValued(a_ndin.gridIndex()))
      {
        BL_ASSERT(m_irrFAB.getIVS().contains(a_ndin.gridIndex()));
        returnval = &(m_irrFAB(a_ndin,a_nVarLoc));
      } 
      else
      {
        BL_ASSERT(m_regFAB.box().contains(a_ndin.gridIndex()));
        returnval = &(m_regFAB(a_ndin.gridIndex(), a_nVarLoc));
      }
    }
    }
      
    return *returnval;
      
  }
      
  /************
        See comment for const version of this function.
  */
  template <class T>
  T&
  BaseEBCellFAB<T>::operator()(const VolIndex& a_ndin,int  a_nVarLoc,
                               int a_isKnownMultiValued )
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(getRegion().contains(a_ndin.gridIndex()));
    BL_ASSERT((a_nVarLoc >= 0)&&(a_nVarLoc < nComp()));
    BL_ASSERT(!m_ebisBox.isCovered(a_ndin.gridIndex()));
    BL_ASSERT( (a_isKnownMultiValued == -1)
               ||    (a_isKnownMultiValued ==  0)
               ||    (a_isKnownMultiValued == +1) );
      
    T* returnval = NULL;
      
    if ( (a_isKnownMultiValued  == -1) || !m_hasMultiCells)
    {
      returnval = &(m_regFAB(a_ndin.gridIndex(), a_nVarLoc));
    }
    else if (a_isKnownMultiValued  == 1)
    {
      returnval = &(m_irrFAB(a_ndin,a_nVarLoc));
    }
    else if (a_isKnownMultiValued  == 0)
    {
      if ( m_ebisBox.isMultiValued(a_ndin.gridIndex()) )
      {
        returnval = &(m_irrFAB(a_ndin,a_nVarLoc));
      }
      else
      {
        returnval = &(m_regFAB(a_ndin.gridIndex(), a_nVarLoc));
      }
    }
    else
    {
      amrex::Error("bogus a_isKnownMultiValued input");
    }
      
    return *returnval;
  }
      
      
      
  /*************/
  /*************/
  template <class T>
  void
  BaseEBCellFAB<T>:: setDefaultValues()
  {
    m_isDefined = false;
       
  }
  /*************/
  /*************/
  template <class T>
  const EBISBox&
  BaseEBCellFAB<T>::getEBISBox() const
  {
    return m_ebisBox;
  }
      
  /******************/
  /// Below lies serialization land.  Enter at thy own risk.
  /// Management is not responsible for any gibbering madness resulting 
  /// from ignoring this warning.
  /******************/
  template <class T>
  std::size_t 
  BaseEBCellFAB<T>::
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
  BaseEBCellFAB<T>::
  nBytesFull() const
  {
    BL_PROFILE("BaseEBCellFAB::nbytesfull");
    size_t retval = 0;
    size_t incrval = 0;
    //first the meta data
    //ncomp
    retval += sizeof(int);
    //region and regfab region
    retval += m_region.linearSize();
    retval += m_regFAB.box().linearSize();

    //these make the ebisbox
    incrval =  m_ebisBox.getEBGraph().nBytesFull();
    retval +=  incrval;
    incrval =  m_ebisBox.getEBData( ).nBytesFull();
    retval +=  incrval;
    //from the above we can call define
    //now for the data
    retval +=  m_regFAB.nBytes(m_regFAB.box(), 0, nComp());
    retval +=  m_irrFAB.nBytes(      m_region, 0, nComp());

    return retval;
  }
      
  /*************/
  template <class T>
  std::size_t 
  BaseEBCellFAB<T>::
  copyToMemFull(void*      dst) const
  {
    BL_PROFILE("BaseEBCellFAB::copyToMemFull");
    size_t retval = 0;
    //first the meta data
    //ncomp
    int* intbuf = (int *) dst;
    *intbuf = m_regFAB.nComp();    
    intbuf++;
    retval += sizeof(int);
    size_t incrval  = 0;

    //back to byte stream 
    unsigned char* buffer = (unsigned char*) intbuf;
    //region and box for fab
    Box  regfabbox = m_regFAB.box();

    m_region.linearOut(buffer);
    incrval = m_region.linearSize();
    buffer += incrval;
    retval += incrval;

    regfabbox.linearOut(buffer);
    incrval = regfabbox.linearSize();
    buffer += incrval;
    retval += incrval;


    //these make the ebisbox
    incrval =  m_ebisBox.getEBGraph().copyToMemFull(buffer);
    buffer += incrval;
    retval += incrval;
    incrval =  m_ebisBox.getEBData( ).copyToMemFull(buffer);
    buffer += incrval;
    retval += incrval;
    //from the above we can call define
    //now for the data

    incrval =  m_regFAB.copyToMem(regfabbox, 0, nComp(), buffer);
    buffer += incrval;
    retval += incrval;
    incrval =  m_irrFAB.copyToMem( m_region, 0, nComp(), buffer);
    buffer += incrval;
    retval += incrval;
    return retval;
  }
  /*************/
  template <class T>
  std::size_t 
  BaseEBCellFAB<T>::
  copyFromMemFull(const void* src)
  {
    BL_PROFILE("BaseEBCellFAB::copyFromMemFull");

    size_t retval = 0;

    //first the meta data
    //ncomp
    int* intbuf = (int *) src;
    int ncomp = *intbuf;
    *intbuf = m_regFAB.nComp();    
    intbuf++;
    retval += sizeof(int);

    unsigned char* buffer = (unsigned char*) intbuf;
    Box region;
    Box regfabbox;
    size_t incrval  = 0;

    //back to byte stream 
    region.linearIn(buffer);
    incrval = m_region.linearSize();
    buffer += incrval;
    retval += incrval;

    regfabbox.linearIn(buffer);
    incrval = regfabbox.linearSize();
    buffer += incrval;
    retval += incrval;

    //these make the ebisbox
    EBGraph ebgraph;
    EBData  ebdata;
    incrval =  ebgraph.copyFromMemFull(buffer);
    buffer += incrval;
    retval += incrval;
    incrval =   ebdata.copyFromMemFull(buffer);
    buffer += incrval;
    retval += incrval;
    //from the above we can call define
    EBISBox ebisBox(ebgraph, ebdata);
    //need to send regfabbox because that is what defines the reg fab
    //the region is intersected with the domain
    this->define(ebisBox, regfabbox, ncomp);

    BL_ASSERT(m_regFAB.box() == regfabbox);
    BL_ASSERT(        region == getRegion());

    //now for the data 
    incrval =  m_regFAB.copyFromMem(regfabbox, 0, ncomp, buffer);
    buffer += incrval;
    retval += incrval;
    incrval =  m_irrFAB.copyFromMem(region   , 0, ncomp, buffer);
    buffer += incrval;
    retval += incrval;
    return retval;
  }
  /*************/
  template <class T>
  std::size_t 
  BaseEBCellFAB<T>::
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
  BaseEBCellFAB<T>::
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
}
      
#endif
