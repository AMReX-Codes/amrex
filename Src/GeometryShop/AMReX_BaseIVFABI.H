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

#ifndef _BASEIVFABI_H_
#define _BASEIVFABI_H_

#include "AMReX_IntVectSet.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_SPMD.H"
#include "AMReX_parstream.H"    
namespace amrex
{
    
  /******************/
  template <class T> inline
  BaseIVFAB<T>::BaseIVFAB():m_nComp(0),m_data(0)
  {
    m_verbose = false;
  }
  /******************/
  template <class T> inline
  BaseIVFAB<T>::~BaseIVFAB()
  {
     
  }
    
  /******************/
  template <class T> inline
  BaseIVFAB<T>::BaseIVFAB(const IntVectSet& a_ivsin,
                          const EBGraph&    a_ebgraph,
                          const int&        a_nvarin)
  {
    define(a_ivsin, a_ebgraph, a_nvarin);
  }
  /******************/
  template <class T> inline
  void
  BaseIVFAB<T>::define(const IntVectSet& a_ivsin,
                       const EBGraph&    a_ebGraph,
                       const int&        a_nvarin)
  {
     
    assert(a_nvarin > 0);
    
    m_verbose = false;
    m_nComp=a_nvarin;
    
    //figure out how long vector has to be
    m_ivs = a_ivsin;
    VoFIterator vofit(a_ivsin, a_ebGraph);
    m_vofs = vofit.getVector();

    std::sort(m_vofs.begin(), m_vofs.end());
    
    int nVoFs = m_vofs.size();
    
    //now allocate the vector set the fab to go into
    //the pool of data at the first component of the first
    //vof
    m_Memory.resize(nVoFs*a_nvarin);
    if(nVoFs > 0)
    {
      m_data = &(m_Memory[0]);
    }
    else
    {
      m_data = NULL;
    }
                  
    
  }
    
  template <class T> inline
  BaseIVFAB<T>&
  BaseIVFAB<T>::copy(const BaseIVFAB<T> & a_src,
                     const Box&           a_srcbox,
                     int                  a_srccomp,
                     const Box&           a_dstbox,
                     int                  a_dstcomp,
                     int                  a_numcomp)
  {
    
    for(unsigned int ivof=0; ivof<m_vofs.size(); ivof++)
    {
      const VolIndex& vof = m_vofs[ivof];
      if(a_srcbox.contains(vof.gridIndex()))
      {
        if(a_src.m_ivs.contains(vof.gridIndex()))
        {
	  int src_ivof = a_src.getVolIndexIndex(vof);
          for(int icomp = 0; icomp < a_numcomp; icomp++)
          {
            int isrc = a_srccomp + icomp;
            int idst = a_dstcomp + icomp;
            (*this)(ivof, idst) = a_src(src_ivof, isrc);
          }
        }
      }
    }
    
    return *this;
  }
  template <class T> inline
  BaseIVFAB<T>&
  BaseIVFAB<T>::plus(const BaseIVFAB<T> & a_src,
                     const Box&        srcbox,
                     const Box&        dstbox,
                     int               srccomp,
                     int               dstcomp,
                     int               numcomp)
  {
    amrex::Error("not implemented");
  }
    
      
  /********************/
  template <class T> inline
  T*
  BaseIVFAB<T>::getIndex(const VolIndex& a_vof, const int& a_comp) const
  {    
    BL_ASSERT((a_comp >= 0) && (a_comp < this->m_nComp));
    int ivof = this->getVolIndexIndex(a_vof);
    return this->getIndex(ivof, a_comp);
  }
    
  /********************/
  template <class T> inline
  void BaseIVFAB<T>::setVal(const T& a_val)
  {
    for(int i=0; i<m_vofs.size()*m_nComp; i++)
    {
      m_data[i]=a_val;
    }
  }
  /********************/
  template <class T> inline
  void BaseIVFAB<T>::setVal(int a_comp, const T& a_val)
  {
    for(int ivof = 0; ivof < m_vofs.size(); ivof++)
    {
      (*this)(m_vofs[ivof], a_comp) = a_val;
    }
  }
    
  /********************/
  template <class T> inline
  int BaseIVFAB<T>::numVoFs() const
  {
    return m_vofs.size();
  }
  /********************/
  template <class T> inline
  const Array<VolIndex>&
  BaseIVFAB<T>::getVoFs() const
  {
    return m_vofs;
  }
  /***********/
  template <class T>
  template <typename F>
  void BaseIVFAB<T>::forall(const BaseIVFAB<T>& a_src, 
                            const Box& a_box, 
                            int a_srccomp, 
                            int a_destcomp, 
                            int a_numcomp, 
                            bool sameRegBox, 
                            const F& func)
  {
    if (sameRegBox)
    {
      Real* l = dataPtr(a_destcomp);
      const Real* r = a_src.dataPtr(a_srccomp);
      int nvof = m_vofs.size();
      for (int i=0; i<a_numcomp*nvof; i++)
      {
        func(l[i], r[i]);
      }
    }
    else
    {
      const Array<VolIndex>& vofs = m_vofs;
      for(int i=0; i<vofs.size(); i++)
      {
        const VolIndex& vof=vofs[i];
        if(a_box.contains(vof.gridIndex()))
        {
          for (int icomp = 0; icomp < a_numcomp; ++icomp)
          {
            func(this->operator()(vof, a_destcomp+icomp), a_src(vof, a_srccomp+icomp));
          }
        }
      }
    }
  }
    
  /******************/
  /// Below lies serialization land.  Enter at thy own risk.
  /// Management is not responsible for any gibbering madness resulting 
  /// from ignoring this warning.
  /******************/
  template<class T> inline
  void
  BaseIVFAB<T>::
  getVoFSubset(Array<VolIndex>& a_vofsubset, const Box& a_region) const
  {
    a_vofsubset.resize(0);
    for(unsigned int ivof=0; ivof<m_vofs.size(); ivof++)
    {
      const VolIndex& vof = m_vofs[ivof];
      if(a_region.contains(vof.gridIndex()))
      {
        a_vofsubset.push_back(vof);
      }
    }
  }
  /******************/
  template <class T> inline
  std::size_t 
  BaseIVFAB<T>::
  nBytes (const Box& a_region, int start_comp, int ncomps) const  
  {
    
    size_t retval = 0;

    Array<VolIndex> vofsubset;
    getVoFSubset(vofsubset, a_region);

    retval += linearListSize(vofsubset);
    for(unsigned int ivof=0; ivof< vofsubset.size(); ivof++)
    {
      const VolIndex& vof = vofsubset[ivof];

      for(int icomp = start_comp; icomp < start_comp+ncomps; icomp++)
      {
        const T& value = (*this)(vof, icomp);
        retval += linearSize(value);
      }
    }
    return retval;
  }

  /********************/
  template <class T> inline
  std::size_t 
  BaseIVFAB<T>::
  copyToMem (const Box& a_region,
             int        srccomp,
             int        numcomp,
             void*      a_buf) const
  {
    unsigned char* buf = (unsigned char*)a_buf;
    size_t retval = 0;
    size_t incval = 0;

    Array<VolIndex> vofsubset;
    getVoFSubset(vofsubset, a_region);

    linearListOut(buf, vofsubset);
    incval = linearListSize(vofsubset);
    
    buf    += incval;
    retval += incval;

    for(unsigned int ivof=0; ivof< vofsubset.size(); ivof++)
    {
      const VolIndex& vof = vofsubset[ivof];
      int endcomp = srccomp + numcomp;
      for(int icomp = srccomp; icomp < endcomp; icomp++)
      {
        const T& value = (*this)(vof, icomp);
        linearOut(buf, value);
        incval  = linearSize(value);
        buf    += incval;
        retval += incval;
      }
    }


    return retval;
  }
    
  /********************/
  template <class T> inline
  std::size_t 
  BaseIVFAB<T>::
  copyFromMem (const Box&  a_region,
               int         dstcomp,
               int         numcomp,
               const void* a_buf)
  {
    unsigned char* buf = (unsigned char*)a_buf;
    size_t retval = 0;
    size_t incval = 0;
    
    Array<VolIndex> vofsubset;

    linearListIn(vofsubset, buf);

    incval = linearListSize(vofsubset);
    buf    += incval;
    retval += incval;

    for(unsigned int ivof=0; ivof< vofsubset.size(); ivof++)
    {
      const VolIndex& vof = vofsubset[ivof];
      int endcomp = dstcomp + numcomp;
      for(int icomp = dstcomp; icomp < endcomp; icomp++)
      {
        T& value = (*this)(vof, icomp);
        linearIn(value, buf);
        incval  = linearSize(value);
        buf    += incval;
        retval += incval;
      }
    }


    return retval;
  }
    
    
  /******************/
  template <class T> inline
  std::size_t 
  BaseIVFAB<T>::
  nBytesFull() const  
  {
    
    size_t retval = 0;
    //this one is for ncomp
    retval = sizeof(int);

    //set
    retval += linearSize(m_ivs);
    //vofs
    retval += linearListSize(m_vofs);
    //data
    for(unsigned int ivof=0; ivof<m_vofs.size(); ivof++)
    {
      const VolIndex& vof = m_vofs[ivof];
      for(int icomp = 0; icomp < m_nComp; icomp++)
      {
        const T& value = (*this)(vof, icomp);
        retval +=linearSize(value);
      }
    }
    return retval;
  }
  /********************/
  template <class T> inline
  std::size_t 
  BaseIVFAB<T>::
  copyToMemFull(void*      a_buf) const
  {
    unsigned char* buf = (unsigned char*)a_buf;
    
    size_t retval = 0;
    size_t incval = 0;
    
    //this one is for ncomp
    int* intbuf = (int*)a_buf;
    *intbuf = m_nComp;
    incval = sizeof(int);
    retval += incval;
    buf    += incval;
    
    //set
    linearOut(buf, m_ivs);
    incval = linearSize(m_ivs);
    retval += incval;
    buf    += incval;
    
    //vofs
    linearListOut(buf, m_vofs);
    incval = linearListSize(m_vofs);
    retval += incval;
    buf    += incval;
    
    //data
    for(unsigned int ivof=0; ivof<m_vofs.size(); ivof++)
    {
      const VolIndex& vof = m_vofs[ivof];
      for(int icomp = 0; icomp < m_nComp; icomp++)
      {
        const T& value = (*this)(vof, icomp);
        linearOut(buf, value);
        incval = linearSize(value);
        retval += incval;
        buf    += incval;
      }
    }
    return retval;
  }
    
  /********************/
  template <class T> inline
  std::size_t 
  BaseIVFAB<T>::
  copyFromMemFull(const void* a_buf)
  {
    unsigned char* buf = (unsigned char*)a_buf;
    
    size_t retval = 0;
    size_t incval = 0;
    
    //this one is for ncomp
    int* intbuf = (int*)a_buf;
    m_nComp = *intbuf;
    incval = sizeof(int);
    retval += incval;
    buf    += incval;
    
    //set
    linearIn(m_ivs, buf);
    incval = linearSize(m_ivs);
    retval += incval;
    buf    += incval;
    
    //vofs
    linearListIn(m_vofs, buf);
    incval = linearListSize(m_vofs);
    retval += incval;
    buf    += incval;


    //data
    //define the memory and data holder
    int nVoFs = m_vofs.size();
    //now allocate the Array set the fab to go into
    //the pool of data at the first component of the first
    //vof
    m_Memory.resize(nVoFs*m_nComp);
    
    if(nVoFs > 0)
    {
      m_data = &(m_Memory[0]);
    }
    else
    {
      m_data = NULL;
    }
        
    for(unsigned int ivof=0; ivof<m_vofs.size(); ivof++)
    {
      const VolIndex& vof = m_vofs[ivof];
      for(int icomp = 0; icomp < m_nComp; icomp++)
      {
        T value;
        linearIn(value, buf);
        incval = linearSize(value);

        (*this)(vof, icomp) = value;
        retval += incval;
        buf    += incval;
      }
    }
    return retval;
  }
}
    
#endif
    

