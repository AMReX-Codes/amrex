#ifndef _IRREGFABFACTORY_H___
#define _IRREGFABFACTORY_H___

#include <cmath>
#include <cstdlib>
#include "AMReX_EBISLayout.H"
#include "AMReX_IrregFAB.H"
#include "AMReX_LayoutData.H"
#include "AMReX_FabFactory.H"

namespace amrex
{
/// Factory class to produce IrregFABs.
/**
   Factory class to produce IrregFABs.
   This is needed for LevelData construction.
*/
  class IrregFABFactory
    : public FabFactory<IrregFAB>
  {
  public:
    /// factory function.
    /**
       Creates a new baseIrregFAB object
       and returns a pointer to it.  Responsiblitly
       for calling operator 'delete' on this pointer is passed to the user.
    */
    virtual IrregFAB* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const override;

    virtual IrregFABFactory* clone () const override {
        return new IrregFABFactory(*this);
    }

    ///
    IrregFABFactory(const EBISLayout                          & a_ebisl, 
                    const shared_ptr<LayoutData<IntVectSet> > & a_sets)
      {
        m_graphs = a_ebisl.getAllGraphs();
        m_sets   = a_sets;
        m_useSets = true;
      }
    
    ///
    IrregFABFactory(const EBISLayout& a_ebisl)
      {
        m_graphs = a_ebisl.getAllGraphs();
        m_useSets = false;
      }

    ///
    virtual ~IrregFABFactory()
      { }

  protected:
    bool m_useSets;
    shared_ptr<FabArray<EBGraph> > m_graphs;
    shared_ptr<LayoutData<IntVectSet> > m_sets;

  private:
    //weak construction should be avoided.
    IrregFABFactory();
  };

}
#endif
