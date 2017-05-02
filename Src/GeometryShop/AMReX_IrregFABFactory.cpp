#include "AMReX_IrregFABFactory.H"

namespace amrex
{

  IrregFAB*
  IrregFABFactory::
  create (const Box& box, int ncomps, const FabInfo& info, int box_index) const
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

    return new  IrregFAB(ivs, graph, ncomps);
  }
}
