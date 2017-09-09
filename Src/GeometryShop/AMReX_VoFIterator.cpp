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

#include "AMReX_VoFIterator.H"
#include "AMReX_EBGraph.H"
#include <cassert>

namespace amrex
{

  /********************************/
  const Array<VolIndex> &
  VoFIterator::getVector() const
  {
    assert(m_isDefined);
    return m_vols;
  }
  /********************************/
  VoFIterator::VoFIterator()
  {
    m_isDefined = false;
  }
  /********************************/
  VoFIterator::~VoFIterator()
  {
  }
  /********************************/
  VoFIterator::VoFIterator(const IntVectSet& a_ivs,
                           const EBGraph   & a_ebgraph)
  {
    define(a_ivs, a_ebgraph);
  }
  /********************************/
  void
  VoFIterator::define(const IntVectSet& a_ivs,
                      const EBGraph   & a_ebgraph)
  {
    m_isDefined = true;
    m_vols.resize(0);
    for (IVSIterator ivsit(a_ivs); ivsit.ok(); ++ivsit)
    {
      Array<VolIndex> vols = (a_ebgraph.getVoFs(ivsit()));
      for(int ivol = 0; ivol < vols.size(); ivol++)
      {
        m_vols.push_back(vols[ivol]);
      }
    }
    reset();
  }
  /********************************/
  VoFIterator::VoFIterator(const IntVectSet& a_ivs,
                           const EBGraphImplem   & a_ebgraph)
  {
    m_isDefined = true;
    m_vols.resize(0);
    for (IVSIterator ivsit(a_ivs); ivsit.ok(); ++ivsit)
    {
      Array<VolIndex> vols = (a_ebgraph.getVoFs(ivsit()));
      for(int ivol = 0; ivol < vols.size(); ivol++)
      {
        m_vols.push_back(vols[ivol]);
      }
    }
    reset();
  }
  /********************************/
  void
  VoFIterator::reset()
  {
    assert(isDefined());
    m_ivol = 0;
  }
  /********************************/
  void
  VoFIterator::operator++()
  {
    assert(isDefined());
    m_ivol++;
  }
  /********************************/
  const VolIndex&
  VoFIterator::operator() () const
  {
    assert(isDefined());
    assert(m_ivol < m_vols.size());
    return m_vols[m_ivol];
  }
  /********************************/
  bool
  VoFIterator::ok() const
  {
    return (m_ivol < m_vols.size());
  }
  /********************************/
  bool
  VoFIterator::isDefined() const
  {
    return m_isDefined;
  }
}
