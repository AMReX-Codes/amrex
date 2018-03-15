#ifndef BL_BASEIFFABFACTORY_H_
#define BL_BASEIFFABFACTORY_H_

#include <AMReX_FabFactory.H>
#include <AMReX_Box.H>
#include <AMReX_BaseIFFAB.H>

namespace amrex
{
  using std::shared_ptr;
  ///factory to create BaseIFFABs
  template <class T>
  class BaseIFFactory: public FabFactory<BaseIFFAB<T> >
  {
  public:
    ///
    BaseIFFactory()
      {
      }

    ///
    virtual ~BaseIFFactory () {}

    ///
    virtual BaseIFFAB<T>* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const
      {
        EBGraph& graph = (*m_graphs)[box_index];
        IntVectSet ivs = graph.getIrregCells(box);
        return new  BaseIFFAB<T>(ivs, graph, m_direction, ncomps);
      }

    virtual BaseIFFactory<T>* clone () const override {
        return new BaseIFFactory<T>(*this);
    }

    ///
    BaseIFFactory(shared_ptr<FabArray<EBGraph> > a_graphs, int a_direction)
      {
        m_graphs = a_graphs;
        m_direction = a_direction;
      }

    ///
    void
    setGraphs(shared_ptr<FabArray<EBGraph> > a_graphs, int a_direction)
      {
        m_graphs    = a_graphs;
        m_direction = a_direction;
      }

  private:
    int                            m_direction;
    shared_ptr<FabArray<EBGraph> > m_graphs;
  };

}

#endif
