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
#include <vector>
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
        if(m_vofs) return true;
        return false;
      }

    ///
    virtual ~BaseIVFAB();

    ///
    /**
       Full define function.  Specifies the irregular domain and the
       number of data components per VoF.  The irregular domain
       must lie completely within the EBGraph.   The contents
       are uninitialized.  If it has previously been defined, the old definition
       data is overwritten and lost.
    */
    virtual void define(const IntVectSet& a_region,
                        const EBGraph&    a_ebgraph,
                        const int&        a_nvarin);

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

    std::size_t nBytes (const Box& bx, int start_comp, int ncomps) const;

    std::size_t copyToMem (const Box& srcbox,
                           int        srccomp,
                           int        numcomp,
                           void*      dst) const;


    std::size_t copyFromMem (const Box&  dstbox,
                             int         dstcomp,
                             int         numcomp,
                             const void* src);

    std::vector<VolIndex>  getVoFs() const;

    ///slow!
    IntVectSet getIVS() const
      {
        IntVectSet retval;
        for (unsigned int i=0; i<m_vofs->size(); i++)
        {
          retval |= (*m_vofs)[i].gridIndex();
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
        if(m_vofs) return m_data+a_ivar*m_vofs->size();
        return m_data;
      }
  
    const T* dataPtr(int a_ivar) const
      {
        if(m_vofs) return m_data+a_ivar*m_vofs->size();
        return m_data;
      }
  
    void setVal(const T& a_val);

    void setVal(int a_comp, const T& a_val);

    void setVal(const T& a_val, const Box& a_box, int a_startcomp, int a_ncomp);


    ///get index into vector
    /**
       needs to be public so that bulk stencils can be constructed
    */
    virtual T* getIndex(const VolIndex& a_vof,const int& a_comp) const;

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
        const VolIndex* vof = dynamic_cast< const VolIndex* >(&a_vof);
        if (vof == NULL) amrex::Error("cast failed:BaseIVFAB only takes vofs for indexing");
        unsigned int i = 0;
        for(;i<m_vofs->size(); i++)
          if(*vof == (*m_vofs)[i]) break;
        i+=a_ivar*m_vofs->size();
        return i;
      }

    ///
    T& operator() (const VolIndex& a_vof, int a_comp)
      {  return *(getIndex(a_vof, a_comp));}

    ///
    const T& operator() (const VolIndex& a_vof, int a_comp) const
      {
        return *(getIndex(a_vof, a_comp));
      }

    template <typename F>
    void forall(const BaseIVFAB& a_src, const Box& a_box, int a_srccomp, int a_destcomp, int a_numcomp, bool sameRegBox, const F& func);

    int numVoFs() const;

    int nComp() const {return m_nComp;}
  private:

  protected:
    EBGraph m_ebgraph;
    shared_ptr<std::vector<VolIndex> > m_vofs;
    shared_ptr<std::vector<T> > m_Memory;
    int m_nComp = 0;
    T* m_data   = nullptr;

  };
}
#include "AMReX_BaseIVFABI.H"

#endif
