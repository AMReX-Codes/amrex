#ifndef BL_EBDataFACTORY_H_
#define BL_EBDataFACTORY_H_

#include <AMReX_FabArray.H>
#include <AMReX_FabFactory.H>
#include <AMReX_Box.H>
#include <AMReX_EBData.H>
#include <AMReX_EBGraph.H>

namespace amrex
{

  ///factory to create EBDatas
  class EBDataFactory : public FabFactory<EBData>
  {
  public:
    EBDataFactory()
      {
      }

    ///
    virtual ~EBDataFactory () {}

    ///
    virtual EBData* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const override
      {
        EBGraph& graph = (*m_graphs)[box_index];
        EBData* retval =  new  EBData();
        retval->define(graph, box, m_dx, m_hasMoments);
        return retval;
      }

    virtual EBDataFactory* clone () const override 
      {
        return new EBDataFactory(*this);
      }

    ///
    EBDataFactory(const shared_ptr<FabArray<EBGraph> >& a_graphs, bool a_hasMoments, Real a_dx)
      {
        m_graphs = a_graphs;
        m_hasMoments = a_hasMoments;
        m_dx = a_dx;
      }

    ///
    void
    setGraphs(const shared_ptr<FabArray<EBGraph> >& a_graphs)
      {
        m_graphs = a_graphs;
      }

  private:
    shared_ptr<FabArray<EBGraph> > m_graphs;
    Real m_dx;
    bool m_hasMoments;
  };

}

#endif
