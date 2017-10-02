
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
#include <AMReX_LayoutData.H>
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
    virtual BaseIVFAB<T>* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const override
      {
        EBGraph& graph       = (*m_graphs)[box_index];
        IntVectSet ivs;
        if(m_sets)
        {
          ivs   = (*m_sets)[box_index];
        }
        else
        {
          ivs = graph.getIrregCells(graph.getRegion());
        }

        BaseIVFAB<T>* retval = new  BaseIVFAB<T>(ivs, graph, ncomps);
        return retval;
      }

    virtual BaseIVFactory<T>* clone () const override {
        return new BaseIVFactory<T>(*this);
    }

    ///
    BaseIVFactory(shared_ptr<FabArray<EBGraph> > a_graphs)
      {
        m_graphs = a_graphs;
      }

    /// sets the irregulsr cells as the IntVectSet for the BaseIVFABs
    BaseIVFactory(const EBISLayout& a_ebisl)
      {
        m_graphs = a_ebisl.getAllGraphs();
      }

    ///
    BaseIVFactory(const EBISLayout                          & a_ebisl, 
                  const shared_ptr<LayoutData<IntVectSet> > & a_sets)
      {
        m_graphs = a_ebisl.getAllGraphs();
        m_sets   = a_sets;
      }
    ///
    void
    setGraphs(shared_ptr<FabArray<EBGraph> > a_graphs)
      {
        m_graphs = a_graphs;
      }

  protected:
    shared_ptr<FabArray<EBGraph> > m_graphs;
    shared_ptr<LayoutData<IntVectSet> > m_sets;
  };

}

#endif
