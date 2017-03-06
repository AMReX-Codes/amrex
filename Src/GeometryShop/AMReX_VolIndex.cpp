/*
 *      .o.       ooo        ooooo ooooooooo.             ooooooo  ooooo 
 *     .888.      `88.       .888' `888   `Y88.            `8888    d8'  
 *    .8"888.      888b     d'888   888   .d88'  .ooooo.     Y888..8P    
 *   .8' `888.     8 Y88. .P  888   888ooo88P'  d88' `88b     `8888'     
 *  .88ooo8888.    8  `888'   888   888`88b.    888ooo888    .8PY888.    
 * .8'     `888.   8    Y     888   888  `88b.  888    .o   d8'  `888b   
 *o88o     o8888o o8o        o888o o888o  o888o `Y8bod8P' o888o  o88888o 
 *
 */

#include "AMReX_VolIndex.H"
#include <cassert>



using std::ostream;
namespace amrex
{
  int VolIndex::linearSize() const
  {
    return (BL_SPACEDIM + 1)*sizeof(int);
  }

  void VolIndex::linearOut(void* const a_outBuf) const
  {
    assert(m_isDefined);
    int* buf = (int*)a_outBuf;
    D_TERM(buf[0]=m_iv[0],; buf[1]=m_iv[1],; buf[2]=m_iv[2]);
    buf[BL_SPACEDIM]=m_cellIndex;

  }

  void VolIndex::linearIn(const void* const inBuf)
  {
    int* buf = (int*)inBuf;
    D_TERM(m_iv[0]=buf[0],; m_iv[1]=buf[1],; m_iv[2]=buf[2]);
    m_cellIndex = buf[BL_SPACEDIM];
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
  VolIndex::operator== (const VolIndex& a_vofin) const
  {
    return ((m_iv == a_vofin.m_iv)&&
            (m_cellIndex == a_vofin.m_cellIndex));
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
}

