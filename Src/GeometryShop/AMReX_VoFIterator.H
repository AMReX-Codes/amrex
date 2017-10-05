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

#ifndef _VOFITERATOR_H_
#define _VOFITERATOR_H_

#include "AMReX_IntVect.H"
#include "AMReX_VolIndex.H"
#include "AMReX_IntVectSet.H"
#include <AMReX_Array.H>
#include <AMReX_Vector.H>

namespace amrex
{
  class EBGraph;
  class EBGraphImplem;

/// Iterator for all vofs within an IntVectSet and an Ebgraph
/**
   Iterator for all vofs within an IntVectSet and an Ebgraph
*/
  class VoFIterator
  {
  public:
    ///
    VoFIterator(const IntVectSet& a_ivs,
                const EBGraph   & a_ebgraph);

    VoFIterator(const IntVectSet   & a_ivs,
                const EBGraphImplem& a_ebgraph);

    ///
    VoFIterator();

    ///
    ~VoFIterator();

    ///
    void define(const IntVectSet& a_ivs,
                const EBGraph   & a_ebgraph);

    ///
    void reset();

    ///
    void operator++();

    ///
    const VolIndex& operator() () const;

    ///
    bool ok() const;

    ///
    bool isDefined() const;

    const Array<VolIndex>& getVector() const;

    int size() const
      {
        return m_vols.size();
      }
  private:
    bool m_isDefined;
    Array<VolIndex> m_vols;
    int m_ivol;
  };
}

#endif
