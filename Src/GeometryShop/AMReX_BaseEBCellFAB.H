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



#ifndef _BASEEBCELLFAB_H_
#define _BASEEBCELLFAB_H_

#include <cmath>
#include <cstdlib>
#include "AMReX_EBISBox.H"
#include "AMReX_BaseIVFAB.H"
#include "AMReX_IntVectSet.H"
#include "AMReX_BaseFab.H"
#include "AMReX_BaseIVFAB.H"

namespace amrex
{
  ///  Array defined at the VolIndexs of an Box in an EBIS.
  /**
     An BaseEBCellFAB is a templated holder for cell-centered
     data over the intersection of a box and an EbisBox.
     At every uncovered VolIndex in this intersection,
     the BaseEBCellFAB contains a specified number of data values.
     At singly valued cells,
     the data is stored internally in a BaseFab.  At multiply-valued
     cells, the data is stored internally in an  BaseIVFAB.
     BaseEBCellFAB provides indexing by VolIndex.
     This class does not provide a copy constructor or assignment operator
     for all the usual  reasons.
  */
  template <class T>
  class BaseEBCellFAB
  {
  public:

    /// for AggStencil
    int numDataTypes() const
      {
        return 2;
      }

    /// for AggStencil
    int   dataType(const BaseIndex& a_baseInd) const
      {
        const VolIndex* vofPtr = dynamic_cast<const VolIndex *>(&a_baseInd);
        if (vofPtr == NULL) amrex::Error("Trying to index into an EBCellFAB with something that is not a VolIndex");
        int retval = 0;
        if (m_ebisBox.isMultiValued(vofPtr->gridIndex()))
        {
          retval = 1;
        }
        return retval;
      }

    ///for AggStencil
    long  offset(const BaseIndex& a_baseInd, const int a_ivar) const
      {
        const VolIndex* vofPtr = dynamic_cast<const VolIndex *>(&a_baseInd);
        if (vofPtr == NULL) amrex::Error("Trying to index into an EBCellFAB with something that is not a VolIndex");
        long retval = 0;
        if (m_ebisBox.isMultiValued(vofPtr->gridIndex()))
        {
          retval = m_irrFAB.offset(*vofPtr, a_ivar);
        }
        else
        {
          retval = m_regFAB.offset(vofPtr->gridIndex(), a_ivar);
        }
        return retval;
      }

    ///for AggStencil
    T* dataPtr(int a_dataType, int a_ivar)
      {
        if (a_dataType == 0)
        {
          return m_regFAB.dataPtr(a_ivar);
        }
        else if (a_dataType == 1)
        {
          return m_irrFAB.dataPtr(a_ivar);
        }
        else
        {
          amrex::Error("bogus datatype input into dataPtr");
        }
        return NULL; //just avoiding compiler warnings with this
      }
    ///for AggStencil
    const T* dataPtr(int a_dataType, int a_ivar) const
      {
        if (a_dataType == 0)
        {
          return m_regFAB.dataPtr(a_ivar);
        }
        else if (a_dataType == 1)
        {
          return m_irrFAB.dataPtr(a_ivar);
        }
        else
        {
          amrex::Error("bogus datatype input into dataPtr");
        }
        return NULL; //just avoiding compiler warnings with this
      }

    ///
    BaseEBCellFAB();

    ///
    BaseEBCellFAB(const EBISBox& a_ebisBox, const Box& a_region, int a_nVar);

    ///
    virtual void  define(const EBISBox& a_ebisBox, const Box& a_region, int a_nVar);


    ///
    virtual ~BaseEBCellFAB();

    ///
    void clear();

    ///
    void  setVal(const T& value);

    ///
    void  setVal(int ivar,const T& value);

    ///
    void  setVal(const T& value,
                 const Box& a_box,
                 int        a_nstart,
                 int        a_numcomp);

    ///
    bool isDefined() const;

    ///
    int nComp() const ;

    ///
    const Box& getRegion() const ;

    ///
    const Box& box() const ;

    /// has to be by value because sometimes the Array does not exist.
    Array<VolIndex> getMultiCells() const ;

    ///
    const BaseIVFAB<T>& getMultiValuedFAB() const;

    ///
    BaseIVFAB<T>& getMultiValuedFAB() ;

    ///
    const BaseFab<T>& getSingleValuedFAB() const;

    ///
    BaseFab<T>& getSingleValuedFAB();

    /**
       Arg a_isKnownMultiValued should be set to:
       -1 if arg a_ndin is known to be in a single-valued cell;
       +1 if arg a_ndin is known to be in a multi-valued cell;
       0 if you're not sure (in which case the expensive IntVectSet::contains()
       function gets called.

       Note: BVS says third arg is a wart that makes this operator look less like
       a logical array access.
    */
    const T& operator() (const VolIndex& a_ndin, int a_nVarLoc,
                         int a_isKnownMultiValued=0) const;

    /**
       Arg a_isKnownMultiValued should be set to:
       -1 if arg a_ndin is known to be in a single-valued cell;
       +1 if arg a_ndin is known to be in a multi-valued cell;
       0 if you're not sure (in which case the expensive IntVectSet::contains()
       function gets called.
    */
    T& operator() (const VolIndex& a_ndin,int  a_nVarLoc,
                   int a_isKnownMultiValued=0);

    ///
    BaseEBCellFAB<T>& copy(const BaseEBCellFAB<T>& source,
                           const Box&        srcbox,
                           int               srccomp,
                           const Box&        destbox,
                           int               destcomp,
                           int               numcomp);


    BaseEBCellFAB<T>& plus(const BaseEBCellFAB<T>& src,
                           const Box&        srcbox,
                           const Box&        destbox,
                           int               srccomp,
                           int               destcomp,
                           int               numcomp=1)
      {
        amrex::Error("this is not really a thing");
        return *this;
      }
    ///
    BaseEBCellFAB<T>& copy(const BaseEBCellFAB<T>& source);

    ///this is for copy and fillboundary when meta data is already filled
    std::size_t nBytes (const Box& bx, int start_comp, int ncomps) const;

    ///this includes all meta data
    std::size_t nBytesFull() const;

    ///this is for copy and fillboundary when meta data is already filled
    std::size_t copyToMem (const Box& srcbox,
                           int        srccomp,
                           int        numcomp,
                           void*      dst) const;

    ///default false.   can make it true for T where you  they are fixed size
    static bool preAllocatable () { return false; }
    static bool isCopyOMPSafe () { return true; }

    ///this includes all meta data
    std::size_t copyToMemFull(void*      dst) const;


    ///this is for copy and fillboundary when meta data is already filled
    std::size_t copyFromMem (const Box&  dstbox,
                             int         dstcomp,
                             int         numcomp,
                             const void* src);

    ///this includes all meta data
    std::size_t copyFromMemFull(const void* src);


    ///invalid but necessary for fabarray to compile in certain cases
    BaseEBCellFAB(const Box& a_region,
                  int        n = 1,
                  bool       alloc = true,
                  bool       shared = false)
      {
        //this strictly for checkpointing 
        m_region = a_region;
      }

    //same as above. 
    void resize (const Box& b,
                 int        N = 1)
      {
        m_region = b;
      }

    ///
    const EBISBox& getEBISBox() const;

    ///
    virtual void
    setCoveredCellVal(const T&    a_val,
                      const int&  a_comp);


  protected:

    ///data at multi-valued cells
    BaseIVFAB<T> m_irrFAB;

    ///data at single-valued cells
    BaseFab<T> m_regFAB;

    EBISBox m_ebisBox;

    Box m_region;

    bool m_hasMultiCells;

    ///has full define function been called?
    bool m_isDefined;

  protected:
    void
    setDefaultValues();

  private:
    //disallowed for perfomance reasons
    void operator= (const BaseEBCellFAB<T>& ebcin)
      {
        amrex::Error("operator= for baseebcellfab undefined");
      }
    BaseEBCellFAB (const BaseEBCellFAB<T>& ebcin)
      {
        amrex::Error("copy constructor for baseebcellfab undefined");
      }
  };

}

#include "AMReX_BaseEBCellFABI.H"

#endif
