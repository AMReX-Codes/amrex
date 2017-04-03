
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

#ifndef BL_BASEIVFABFACTORY_H_
#define BL_BASEIVFABFACTORY_H_

#include <AMReX_FabArray.H>
#include <AMReX_FabFactory.H>
#include <AMReX_Box.H>
#include <AMReX_BaseIVFAB.H>
#include <AMReX_EBISLayout.H>

namespace amrex
{

  ///factory to create BaseIVFABs
  template <class T>
  class BaseIVFactory: public FabFactory<BaseIVFAB<T> >
  {
  public:
    BaseIVFactory()
      {
      }

    ///
    virtual ~BaseIVFactory () {}

    ///
    virtual BaseIVFAB<T>* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const
      {
        EBGraph& graph = (*m_graphs)[box_index];
        IntVectSet ivs = graph.getIrregCells(box);
        return new  BaseIVFAB<T>(ivs, graph, ncomps);
      }

    ///
    BaseIVFactory(shared_ptr<FabArray<EBGraph> > a_graphs)
      {
        m_graphs = a_graphs;
      }


    BaseIVFactory(const EBISLayout& a_ebisl)
      {
        m_graphs = a_ebisl.getAllGraphs();
      }

    ///
    void
    setGraphs(shared_ptr<FabArray<EBGraph> > a_graphs)
      {
        m_graphs = a_graphs;
      }

  private:
    shared_ptr<FabArray<EBGraph> > m_graphs;
  };

}

#endif
