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

#ifndef _LayoutData__H_
#define _LayoutData__H_

#include "AMReX_FabArrayBase.H"
#include "AMReX_BoxArray.H"
#include "AMReX_IntVect.H"
#include "AMReX_MFIter.H"

namespace amrex
{
  /// a one-thingy-per-box distributed object
  template <class T>
  class LayoutData: public FabArrayBase
  {
  public:

    ///
    LayoutData()
      {
      }

    LayoutData(const BoxArray            & a_grids, 
               const DistributionMapping & a_dm)
      {
        define(a_grids, a_dm);
      }

    void define(const BoxArray           & a_grids,
                const DistributionMapping & a_dm)
      {
        FabArrayBase::define(a_grids, a_dm, 1, 0);

        m_data.resize(local_size());
      }

    virtual ~LayoutData()
      {
      }

    T& operator[](const MFIter& a_mfi)
      {
        BL_ASSERT(a_mfi.LocalIndex() < m_data.size());
        BL_ASSERT(DistributionMap() == a_mfi.DistributionMap());
        return (m_data[a_mfi.LocalIndex()]);
      }


    const T& operator[](const MFIter& a_mfi) const
      {
        BL_ASSERT(a_mfi.LocalIndex() < m_data.size());
        BL_ASSERT(DistributionMap() == a_mfi.DistributionMap());
        return (m_data[a_mfi.LocalIndex()]);
      }

    T& operator[](const int & a_localIndex)
      {
        BL_ASSERT(a_localIndex < m_data.size());
        return (m_data[a_localIndex]);
      }


    const T& operator[](const int & a_localIndex) const
      {
        BL_ASSERT(a_localIndex < m_data.size());
        return (m_data[a_localIndex]);
      }

  private:
    Array<T> m_data;
        
  };
}
#endif
