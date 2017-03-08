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

#ifndef _BASEIFFABI_H_
#define _BASEIFFABI_H_
#include "AMReX_FaceIterator.H"

namespace amrex
{
  template <class T>
  bool BaseIFFAB<T>::s_verbose  = false;

/******************/
  template <class T> inline
  BaseIFFAB<T>::BaseIFFAB()
  {
    setDefaultValues();
  }
/******************/
  template <class T> inline
  BaseIFFAB<T>::~BaseIFFAB()
  {
    clear();
  }
/******************/
  template <class T> inline
  BaseIFFAB<T>::
  BaseIFFAB(const Box&     a_box,
            const EBGraph& a_ebgraph,
            const int&     a_direction,
            const int&     a_nvarin)
  {
    setDefaultValues();
    define(a_box, a_ebgraph, a_direction, a_nvarin);
  }
/******************/
  template <class T> inline
  void
  BaseIFFAB<T>::
  define(const Box&        a_box,
         const EBGraph&    a_ebgraph,
         const int&        a_direction,
         const int&        a_nvarin)
  {
    clear();
    m_isDefined = true;
    BL_assert(a_nvarin > 0);
    BL_assert((a_direction >= 0) && (a_direction < SpaceDim));

    m_direction = a_direction;
    m_ebgraph = a_ebgraph;
    m_nComp = a_nvarin;
    m_box = a_box;

    m_faces = m_ebgraph.getMultiValuedFaces(a_direction, a_box);
    m_nFaces = m_faces.size();

    m_data = new T[m_nComp*m_nFaces];
  }
/******************/
  template <class T> inline
  void
  BaseIFFAB<T>::
  setVal(const T& a_value)
  {
    BL_assert(isDefined());
    for (int ivec = 0; ivec < m_nFaces*m_nComp; ivec++)
    {
      m_data[ivec] = a_value;
    }
  }
/******************/
  template <class T> inline
  void
  BaseIFFAB<T>::
  setVal(int ivar, const T& a_value)
  {
    BL_assert(isDefined());
    BL_assert(ivar >= 0);
    BL_assert(ivar < m_nComp);

    for (int ivec = ivar; ivec < m_nFaces*m_nComp; ivec += m_nComp)
    {
      m_data[ivec] = a_value;
    }
  }
/******************/
  template <class T> inline
  BaseIFFAB<T>& 
  BaseIFFAB<T>::copy(const BaseIFFAB<T>& a_src,
                     const Box&          a_srcbox,
                     int                 a_srccomp,
                     const Box&          a_destbox,
                     int                 a_destcomp,
                     int                 a_numcomp)
  {
    BL_assert(isDefined());
    BL_assert(a_src.isDefined());
    BL_assert(srccomp >= 0);
    BL_assert(destcomp >= 0);
    BL_assert(numcomp   <= m_nComp);

    if ((!m_box.isEmpty()) && (!a_src.m_box.isEmpty()))
    {
      BL_assert(a_fromBox == a_toBox);
      Box intBox = a_toBox;
      Box boxIntersect = m_box;
      boxIntersect &=  a_src.m_box;
      boxIntersect &=  intBox;

      for (int iface = 0; iface < m_faces.size(); iface++)
      {
        const FaceIndex& face = m_faces[iface];
        const IntVect& ivHi = face.gridIndex(Side::Hi);
        const IntVect& ivLo = face.gridIndex(Side::Lo);
        if (boxIntersect.contains(ivHi) || boxIntersect.contains(ivLo))
        {
          for (int icomp = 0; icomp < a_numcomp; icomp++)
          {
            int isrccomp = a_srccomp  + icomp;
            int idstcomp = a_destcomp + icomp;
            (*this)(face, idstcomp) = a_src(face, isrccomp);
          }
        }
      }
    }
    return *this;
  }
/******************/
  template <class T> inline
  BaseIFFAB<T>& 
  BaseIFFAB<T>::plus(const BaseIFFAB<T>& a_src,
                     const Box&          a_srcbox,
                     int                 a_srccomp,
                     const Box&          a_destbox,
                     int                 a_destcomp,
                     int                 a_numcomp)
  {
    BL_assert(isDefined());
    BL_assert(a_src.isDefined());
    BL_assert(srccomp >= 0);
    BL_assert(destcomp >= 0);
    BL_assert(numcomp   <= m_nComp);

    if ((!m_box.isEmpty()) && (!a_src.m_box.isEmpty()))
    {
      BL_assert(a_fromBox == a_toBox);
      Box intBox = a_toBox;
      Box boxIntersect = m_box;
      boxIntersect &=  a_src.m_box;
      boxIntersect &=  intBox;

      for (int iface = 0; iface < m_faces.size(); iface++)
      {
        const FaceIndex& face = m_faces[iface];
        const IntVect& ivHi = face.gridIndex(Side::Hi);
        const IntVect& ivLo = face.gridIndex(Side::Lo);
        if (boxIntersect.contains(ivHi) || boxIntersect.contains(ivLo))
        {
          for (int icomp = 0; icomp < a_numcomp; icomp++)
          {
            int isrccomp = a_srccomp  + icomp;
            int idstcomp = a_destcomp + icomp;
            (*this)(face, idstcomp) += a_src(face, isrccomp);
          }
        }
      }
    }
    return *this;
  }
/*********/
  template <class T> inline
  T*
  BaseIFFAB<T>::
  getIndex(const FaceIndex& a_face, const int& a_comp) const
  {
    BL_assert(m_isDefined);
    BL_assert((a_comp >= 0) && (a_comp < this->m_nComp));

    T* dataPtr =  this->m_data;
    bool found = false;
    T* retval=NULL;
    for (int iface = 0; iface < m_faces.size(); ++iface)
    {
      if (a_face == m_faces[iface])
      {
        found = true;
        retval=dataPtr;
      }
      if (found) break;
      dataPtr++;
    }
    if (!found)
    {
      MayDay::Error("BaseIFFAB:index not found");
    }
    retval += a_comp*this->m_nFaces;
    return retval;

  }
/*********/
  template <class T> inline
  void
  BaseIFFAB<T>::
  clear()
  {
    if (m_data != NULL)
    {
      delete[] m_data;
      m_data = NULL;
    }
    m_faces.resize(0);
    setDefaultValues();
  }
/*************************/
  template <class T> inline
  T&
  BaseIFFAB<T>::operator() (const FaceIndex& a_ndin,
                            const int& a_comp)
  {
    BL_assert(isDefined());
    BL_assert(a_ndin.direction() == m_direction);

    T* dataPtr =   getIndex(a_ndin, a_comp);
    return(*dataPtr);
  }
/**************************/
  template <class T> inline
  const T&
  BaseIFFAB<T>::operator() (const FaceIndex& a_ndin,
                            const int& a_comp) const
  {
    BL_assert(isDefined());
    BL_assert(a_ndin.direction() == m_direction);

    T* dataPtr =   getIndex(a_ndin, a_comp);
    return(*dataPtr);
  }
/******************/
  template <class T> inline
  T*
  BaseIFFAB<T>::dataPtr(const int& a_comp)
  {
    BL_assert(a_comp >= 0);
    BL_assert(a_comp <= m_nComp);
    return m_data + a_comp*m_nFaces;
  }
/******************/
  template <class T> inline
  const T*
  BaseIFFAB<T>::dataPtr(const int& a_comp) const
  {
    BL_assert(a_comp >= 0);
    BL_assert(a_comp <= m_nComp);
    return m_data + a_comp*m_nFaces;
  }
/******************/
  template <class T> inline
  void
  BaseIFFAB<T>::setDefaultValues()
  {
    m_isDefined = false;
    m_data = NULL;
    m_nFaces = 0;
    m_nComp = 0;
    m_direction = -1;
  }
/******************/
  template <class T> inline
  void
  BaseIFFAB<T>::
  getFaceSubset(std::vector<FaceIndex>& a_subFaces, const Box& a_subBox) const
  {
    a_subFaces.resize(0);
    for (int iface = 0; iface < m_faces.size(); iface++)
    {
      if ((a_subBox.contains(m_faces[iface].gridIndex(Side::Lo))) ||
          (a_subBox.contains(m_faces[iface].gridIndex(Side::Hi))))
      {
        a_subFaces.push_back(m_faces[iface]);
      }
    }
  }
/*********************/
  template <class T> inline
  std::size_t nBytes (const Box& a_region, int start_comp, int ncomps) const;
  {
    BL_assert(isDefined());
    //create set of cells in fab that are also in the input region
    std::vector<FaceIndex> faces;
    getFaceSubset(faces, a_region);

    //account for list of faces
    std::size_t facesize = linearListSize(faces);

    //add for each data pt
    std::size_t datasize = 0;
    for (int iface = 0; iface < faces.size(); iface++)
    {
      int endcomp = start_comp + ncomps - 1;
      for (int icomp = start_comp; icomp <= endcomp; icomp++)
      {
        const T& dataPt = (*this)(faces[iface], icomp);
        int pointsize = linearSize(dataPt);
        datasize += pointsize;
      }
    }

    std::size_t retval = facesize + datasize;

    return retval;
  }
/********************/
  template <class T> inline
  std::size_t 
  BaseIFFAB<T>::
  copyToMem (const Box& a_region,
             int        srccomp,
             int        numcomp,
             void*      a_buf) const
  {
    BL_assert(isDefined());
    std::vector<FaceIndex> faces;
    getFaceSubset(faces, a_region);

    //output the faces.
    unsigned char* charbuffer = (unsigned char *) a_buf;
    linearListOut(charbuffer, faces);
    charbuffer += linearListSize(faces);

    //output the data
    const BaseIFFAB<T>& thisFAB = *this;
    for (int iface = 0; iface < faces.size(); iface++)
    {
      const FaceIndex& face = faces[iface];
      int endcomp = srccomp + numcomp - 1;
      for (int icomp = srccomp; icomp <= endcomp; icomp++)
      {
        //output the data into the buffer
        const T& dataPt =  thisFAB(face, icomp);
        linearOut(charbuffer, dataPt);
        //increment the buffer offset
        charbuffer += linearSize(dataPt);
      }
    }
    size_t retval = nBytes(a_region, srccomp, numcomp);
    return retval;
  }
/********************/
  template <class T> inline
  std::size_t
  BaseIFFAB<T>::
  copyFromMem (const Box&  dstbox,
               int         dstcomp,
               int         numcomp,
               const void* a_buf)
  {
    BL_assert(isDefined());

    unsigned char* charbuffer = (unsigned char *) a_buf;
    std::vector<FaceIndex> faces;
    linearListIn(faces, charbuffer);
    charbuffer += linearListSize(faces);

    for (int iface = 0; iface < faces.size(); iface++)
    {
      const FaceIndex& face = faces[iface];
      int endcomp = dstcomp + numcomp - 1;
      for (int icomp = dstcomp; icomp <= endcomp; icomp++)
      {
        //input the data
        T& dataPt = (*this)(face, icomp);
        CH_XD::linearIn(dataPt, charbuffer) ;
        //increment the buffer offset
        charbuffer += linearSize(dataPt);
      }
    }
    size_t retval = nBytes(a_region, dstcomp, numcomp);
    return retval;
  }
/********************/
/******************/

}
#endif
