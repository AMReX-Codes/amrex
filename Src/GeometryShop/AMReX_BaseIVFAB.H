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


#ifndef _BASEIVFAB_H_
#define _BASEIVFAB_H_

#include <cmath>
#include <cstdlib>
#include <memory>
#include "AMReX_IntVectSet.H"
#include "AMReX_VolIndex.H"
#include "AMReX_BaseFab.H"
#include "AMReX_EBGraph.H"
namespace amrex
{
///
/**
   BaseIVFAB is a templated
   data holder defined at the VoFs of an irregular domain.

   Implemented as just a raw vector of vofs and data, more optimized
   for smaller memory footprint and faster linearIn/linearOut.  will
   be more brutal for vof-by-vof indexing.
   bvs
*/
  template <class T>
  class BaseIVFAB 
  {
  public:
    typedef T value_type;
    ///
    /**
       Null constructor.  
    */
    BaseIVFAB();

    ///
    /**
       Defining constructor.  Specifies the irregular domain
       and the number of data components per VoF.  The irregular domain
       must lie completely within the EBGraph.  The
       contents are uninitialized.  Calls full define function.
    */
    BaseIVFAB( const IntVectSet& a_region,
               const EBGraph&    a_ebgraph,
               const int&        a_nvarin);

    ///
    bool hasVoFs() const
      {
        return (m_vofs.size() != 0);
      }

    ///
    ~BaseIVFAB();

    ///
    /**
       Full define function.  Specifies the irregular domain and the
       number of data components per VoF.  The irregular domain
       must lie completely within the EBGraph.   The contents
       are uninitialized.  If it has previously been defined, the old definition
       data is overwritten and lost.
    */
    void define(const IntVectSet& a_region,
                        const EBGraph&    a_ebgraph,
                        const int&        a_nvarin);

    ///default false.   can make it true for T where you  they are fixed size
    static bool preAllocatable () { return false; }
    static bool isCopyOMPSafe () { return true; }

    ///
    BaseIVFAB<T>& copy(const BaseIVFAB<T> & src,
                       const Box&        srcbox,
                       int               srccomp,
                       const Box&        destbox,
                       int               destcomp,
                       int               numcomp);



    BaseIVFAB<T>& plus(const BaseIVFAB<T>& src,
                       const Box&        srcbox,
                       const Box&        destbox,
                       int               srccomp,
                       int               destcomp,
                       int               numcomp=1);

    ///stuff to get copy and fillbndry to work
    std::size_t nBytes (const Box& bx, int start_comp, int ncomps) const;

    ///stuff to get copy and fillbndry to work
    std::size_t copyToMem (const Box& srcbox,
                           int        srccomp,
                           int        numcomp,
                           void*      dst) const;

    ///stuff to get copy and fillbndry to work
    std::size_t copyFromMem (const Box&  dstbox,
                             int         dstcomp,
                             int         numcomp,
                             const void* src);

    ///includes all meta data
    std::size_t nBytesFull() const;


    ///includes all meta data
    std::size_t copyToMemFull(void* src) const;

    ///includes all meta data
    std::size_t copyFromMemFull(const void* src);


    ///
    const Vector<VolIndex>&  getVoFs() const;

    ///slow!
    IntVectSet getIVS() const
      {
        IntVectSet retval;
        for (unsigned int i=0; i<m_vofs.size(); i++)
        {
          retval |= m_vofs[i].gridIndex();
        }
        return retval;

      }

    //for aggStencil
    const T* dataPtr(int a_dataType, int a_ivar) const
      {
        return dataPtr(a_ivar);
      }
    T* dataPtr(int a_dataType, int a_ivar) 
      {
        return dataPtr(a_ivar);
      }
    
    inline T* dataPtr(int a_ivar)
      {
        return m_data+a_ivar*(m_vofs.size());
      }
  
    const T* dataPtr(int a_ivar) const
      {
        return m_data + (a_ivar*m_vofs.size());
      }
  
    T* dataPtr () { return m_Memory.data(); }

    std::size_t size () const { return m_Memory.size(); }

    void setVal(const T& a_val);

    void setVal(int a_comp, const T& a_val);

    void setVal(const T& a_val, const Box& a_box, int a_startcomp, int a_ncomp);


    ///get index into vector
    /**
       needs to be public so that bulk stencils can be constructed
    */
    T* getIndex(const VolIndex& a_vof,const int& a_comp) const;

    int numDataTypes() const
      {
        return 1;
      }

    /// for AggStencil
    int   dataType(const BaseIndex& a_baseInd) const
      {
        return 0;
      }
    ///for AggStencil
    long offset(const BaseIndex& a_vof, const int& a_ivar) const
      {
        const VolIndex* vofptr = dynamic_cast< const VolIndex* >(&a_vof);
        if (vofptr == NULL) amrex::Error("cast failed:BaseIVFAB only takes vofs for indexing");

        int ivof = -1;
        bool found = false;
        const VolIndex& vof =  *vofptr;
        for(int ivec = 0;ivec < m_vofs.size(); ivec++)
        {
          if(vof == m_vofs[ivec])
          {
            ivof = ivec;
            found = true;
            break;
          }
        }
        if(!found)
        {
          amrex::Error("baseivfab::offset: vof not found in set");
        }
        ivof += a_ivar*m_vofs.size();
        return ivof;
      }

    ///
    T& operator() (const VolIndex& a_vof, int a_comp)
      {  return *(getIndex(a_vof, a_comp));}

    ///
    const T& operator() (const VolIndex& a_vof, int a_comp) const
      {
        return *(getIndex(a_vof, a_comp));
      }

    ///
    T& operator() (int a_ivof, int a_comp)
      {  return *(getIndex(a_ivof, a_comp));}

    ///
    const T& operator() (int a_ivof, int a_comp) const
      {
        return *(getIndex(a_ivof, a_comp));
      }

    int numVoFs() const;

    int nComp() const {return m_nComp;}

    mutable bool m_verbose;
    template<typename F>
    void forall(const BaseIVFAB<T>& a_src, const Box& a_box, int a_srccomp, int a_dstcomp, 
                int a_numcomp, bool sameRegBox, const F& func);

    void getVoFSubset(Vector<VolIndex>& a_vofs,
                      const Box& a_region) const;
  private:

  protected:
    int m_nComp = 0;

    IntVectSet m_ivs;
    Vector<VolIndex>  m_vofs;
    Vector<T>  m_Memory;
    T* m_data   = nullptr;

    // unfortunately, we return T*, not const T*
    T* getIndex (int ivof, int icomp) const {
	return m_data + ivof + icomp*m_vofs.size();
    }

    int getVolIndexIndex (const VolIndex& a_vof) const {
	auto low = std::lower_bound(m_vofs.begin(), m_vofs.end(), a_vof);
	if (low == m_vofs.end() || a_vof != *low) {
	    amrex::Error("attempt to access data from vof that is not in BaseIVFAB");
	}
	return low - m_vofs.begin();
    }
  };
}
#include "AMReX_BaseIVFABI.H"

#endif
