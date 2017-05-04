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

namespace amrex
{

/******************/
  template <class T> inline
  BaseIVFAB<T>::BaseIVFAB():m_nComp(0),m_data(0)
  {
 
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

    m_nComp=a_nvarin;
    m_ebgraph = a_ebGraph;

    if (!a_ivsin.isEmpty())
    {
      int nVoFs = 0;
      //figure out how long vector has to be
      IVSIterator ivsit(a_ivsin);
      for (ivsit.begin(); ivsit.ok(); ++ivsit)
      {
        nVoFs += a_ebGraph.numVoFs(ivsit());
      }
      //now allocate the vector set the fab to go into
      //the pool of data at the first component of the first
      //vof
      if (nVoFs > 0)
      {
        if(!m_vofs)
        {
          m_vofs=shared_ptr<std::vector<VolIndex> >(new std::vector<VolIndex>(nVoFs));
          m_Memory= shared_ptr<std::vector<T> >(new std::vector<T>(nVoFs*a_nvarin));
          m_data = &(m_Memory->operator[](0));
              
          VoFIterator vofit(a_ivsin, a_ebGraph);
          *m_vofs = vofit.getVector();
        }
      }
      else
      {
        m_vofs.reset();
      }
    }
  }

  template <class T> inline
  BaseIVFAB<T>&
  BaseIVFAB<T>::copy(const BaseIVFAB<T> & a_src,
                     const Box&        srcbox,
                     int               srccomp,
                     const Box&        destbox,
                     int               destcomp,
                     int               numcomp)
  {
//    int istop = 0;
//    IntVect ivdeb(D_DECL(17,3,0));
//    if(srcbox.contains(ivdeb))
//    {
//      istop =1;
//    }
    T* to=NULL;
    if(m_vofs)
    {
      to = getIndex((*m_vofs)[0],destcomp);

      for(unsigned int i=0; i<m_vofs->size(); i++, to++)
      {
        const VolIndex& vof = (*m_vofs)[i];
        if(srcbox.contains(vof.gridIndex()))
        {
          const T* from= a_src.getIndex(vof,srccomp);

          for(unsigned int c=0; c<numcomp; c++)
          {
            to[c*m_vofs->size()]=from[c*a_src.m_vofs->size()];
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
                     const Box&        destbox,
                     int               srccomp,
                     int               destcomp,
                     int               numcomp)
  {
    T* to=NULL;
    if(m_vofs)
    {
      to = getIndex((*m_vofs)[0],destcomp);

      for(unsigned int i=0; i<m_vofs->size(); i++, to++)
      {
        const VolIndex& vof = (*m_vofs)[i];
        if(srcbox.contains(vof.gridIndex()))
        {
          const T* from= a_src.getIndex(vof,srccomp);

          for(unsigned int c=0; c<numcomp; c++)
          {
            to[c*m_vofs->size()]+=from[c*a_src.m_vofs->size()];
          }
        }
      }
    }
    return *this;
  }

  
  ///
/******************/
  template <class T> inline
  std::size_t 
  BaseIVFAB<T>::
  nBytes (const Box& a_region, int start_comp, int ncomps) const  
  {
    int count = 0;
    T tmp;
    //create set of cells in fab that are also in the input region
    if(m_vofs)
    {
      for (unsigned int i=0; i<m_vofs->size(); i++)
      {
        if (a_region.contains((*m_vofs)[i].gridIndex())) count++;
      }
    }
    if (count > 0)
    {
      return sizeof(int) + count*linearSize((*m_vofs)[0]) + count*ncomps*linearSize(tmp);
    }
    return sizeof(int);
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
    char* buffer = (char*)a_buf;
    buffer += sizeof(int);
    int count = 0;
    if (m_vofs)
    {
      const T* ptr = &(this->m_data[0]);
      for (unsigned int i=0; i<m_vofs->size(); i++, ptr++)
      {
        const VolIndex& v = (*m_vofs)[i];
        if (a_region.contains(v.gridIndex()))
        {
          count++;
          linearOut(buffer, v);
          buffer+= linearSize(v);
          int endcomp = srccomp + numcomp -1;
          for (int c= srccomp; c<=endcomp; c++)
          {
            linearOut(buffer, *(ptr+c*(m_vofs->size())));
            buffer += linearSize(*ptr);
          }
        }
      }
    }
    int* b = (int*)a_buf;
    *b = count;

    size_t retval = nBytes(a_region, srccomp, numcomp);
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
    int* b = (int*)a_buf;
    int count = *b;
    char* buffer = (char*)a_buf;
    buffer += sizeof(int);
    for (int i=0; i<count; i++)
    {
      VolIndex v;
      linearIn(v, buffer);
      buffer += linearSize(v);
      int endcomp = dstcomp + numcomp -1;
      for (int c= dstcomp; c<= endcomp; c++)
      {
        T* ptr = getIndex(v, c);
        linearIn(*ptr, buffer);
        buffer+= linearSize(*ptr);
      }
    }
    size_t retval = nBytes(a_region, dstcomp, numcomp);
    return retval;
  }
/********************/
  template <class T> inline
  T*
  BaseIVFAB<T>::getIndex(const VolIndex& a_vof, const int& a_comp) const
  {

    assert((a_comp >= 0) && (a_comp < this->m_nComp));

    T* dataPtr =  (T*)&(this->m_data[0]);
    int ivof = -1;
    bool found = false;
    for (unsigned int i=0; i<m_vofs->size(); ++i)
    {
      if (a_vof == (*m_vofs)[i])
      {
        found = true;
        ivof = i;
      }
      if(found) break;
    }
    if(!found)
    {
      amrex::Error("attempt to access data from vof that is not in BaseIVFAB");
    }
    dataPtr += ivof;
    dataPtr += a_comp*m_vofs->size();
    return dataPtr;
  }

/********************/
  template <class T> inline
  void BaseIVFAB<T>::setVal(const T& a_val)
  {
    if(m_vofs)
      for(int i=0; i<m_vofs->size()*m_nComp; i++) m_data[i]=a_val;
  }
/********************/
  template <class T> inline
  void BaseIVFAB<T>::setVal(int a_comp, const T& a_val)
  {
    if(m_vofs)
      for(int i=a_comp*m_nComp; i<(a_comp+1)*m_nComp; i++) m_data[i]=a_val;
  }

/********************/
  template <class T> inline
  void BaseIVFAB<T>::setVal(const T& a_val, const Box& a_region, int a_startcomp, int a_ncomp)
  {
    if(m_vofs)
      forall(*this, a_region, a_startcomp, a_startcomp, a_ncomp, false, [a_val](T& d, const T& s){d=a_val;});
  }



/********************/
  template <class T> inline
  int BaseIVFAB<T>::numVoFs() const
  {
    if(m_vofs)
      return m_vofs->size();
    return 0;
  }
/********************/
  template <class T> inline
  std::vector<VolIndex>
  BaseIVFAB<T>::getVoFs() const
  {
    //if this fails, nothing to get
    if(m_vofs)
    {
      return *m_vofs;
    }
    else
    {
      return std::vector<VolIndex>();
    }
  }
/***********/
  template <class T>
  template <typename F>
  void BaseIVFAB<T>::forall(const BaseIVFAB<T>& a_src, const Box& a_box, int a_srccomp, int a_destcomp, 
                            int a_numcomp, bool sameRegBox, const F& func)
  {
    if(!m_vofs) return;
    if (sameRegBox)
    {
      Real* l = dataPtr(a_destcomp);
      const Real* r = a_src.dataPtr(a_srccomp);
      int nvof = m_vofs->size();
      for (int i=0; i<a_numcomp*nvof; i++)
      {
        func(l[i], r[i]);
      }
    }
    else
    {
      const std::vector<VolIndex>& vofs = *m_vofs;
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

}

#endif

