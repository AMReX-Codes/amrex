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

#include "AMReX_VolIndex.H"
#include <cassert>
#include <iostream>



using std::ostream;
namespace amrex
{
  int VolIndex::linearSize() const
  {
    return (AMREX_SPACEDIM + 1)*sizeof(int);
  }

  void VolIndex::linearOut(void* const a_outBuf) const
  {
    assert(m_isDefined);
    int* buf = (int*)a_outBuf;
    AMREX_D_TERM(buf[0]=m_iv[0],; buf[1]=m_iv[1],; buf[2]=m_iv[2]);
    buf[AMREX_SPACEDIM]=m_cellIndex;

  }

  void VolIndex::linearIn(const void* const inBuf)
  {
    int* buf = (int*)inBuf;
    AMREX_D_TERM(m_iv[0]=buf[0],; m_iv[1]=buf[1],; m_iv[2]=buf[2]);
    m_cellIndex = buf[AMREX_SPACEDIM];
    m_isDefined = true;
  }
  /*****************************************/
  
  const IntVect&
  VolIndex::gridIndex() const
  {
    return m_iv;
  }

  /*****************************************/
  
  int
  VolIndex::cellIndex() const
  {
    return m_cellIndex;
  }
  /*****************************************/
  void
  VolIndex::define(const IntVect& a_ix,const int& a_vofID)
  {
    m_isDefined = true;
    m_iv = a_ix;
    m_cellIndex = a_vofID;
  }
  /*****************************************/
  
  VolIndex::VolIndex(const IntVect& a_ix,const int& a_vofID)
    :BaseIndex(),m_iv(a_ix), m_cellIndex(a_vofID), m_isDefined(true)
  {
    ;
  }
  /*****************************************/
  bool
  VolIndex::isDefined() const
  {
    return m_isDefined;
  }

  /*****************************************/
  
  VolIndex::VolIndex()
  {
    //set index to bogus number to make undefined
    //ones catchable.
    m_cellIndex= -1;
    m_isDefined = false;

  }

  bool
  VolIndex::operator!=(const VolIndex& rhs) const
  {
    return !(*this == rhs);
  }

  /*****************************************/
  void
  VolIndex::define(const VolIndex& a_vofin)
  {
    //must be allowed to propogate undefinednitude
    //because of vector<vol>
    //  CH_assert(a_vofin.isDefined());

    m_isDefined = a_vofin.m_isDefined;
    m_iv = a_vofin.m_iv;
    m_cellIndex = a_vofin.m_cellIndex;
  }

  VolIndex::~VolIndex()
  {
  }

/*****************************************/
  ostream&
  operator<< (ostream&       os,
              const VolIndex& p)
  {
    IntVect iv = p.gridIndex();
    int vofind = p.cellIndex();
    os <<  std::string("((");
    for(int idir = 0; idir < SpaceDim; idir ++)
    {
      os << iv[idir];
      if(idir != SpaceDim-1)
      {
        os << std::string(",");
      }
    }
    os << std::string(")(") << vofind << std::string("))") ;
    if (os.fail())
      amrex::Error("operator<<(ostream&,VolIndex&) failed");
    return os;
  }
}

