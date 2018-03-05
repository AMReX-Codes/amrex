#ifndef _EBCELLFACTORY_H_
#define _EBCELLFACTORY_H_

#include <cmath>
#include <cstdlib>
#include "AMReX_EBISLayout.H"
#include "AMReX_EBCellFAB.H"


/// Factory class to produce EBCellFABs.
/**
   Factory class to produce EBCellFABs.
   This is needed for LevelData construction.
*/
namespace amrex
{
  class EBCellFactory
    : public FabFactory<EBCellFAB>
  {
  public:
    /// factory function.
    /**
       Creates a new baseebcellfab object
       and returns a pointer to it.  Responsiblitly
       for calling operator 'delete' on this pointer is passed to the user.
    */
    virtual EBCellFAB* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const override
      {
        EBISBox ebisBox = m_ebisl[box_index];
        return new EBCellFAB(ebisBox, box, ncomps);
      }

    virtual EBCellFactory* clone () const override {
        return new EBCellFactory(*this);
    }

    EBCellFactory(const EBCellFactory& a_inputs) = default;

    ///
    /**
       create the factory with an ispace.  calls full define function
    */
    EBCellFactory(const EBISLayout& a_ebisl)
      {
        m_ebisl = a_ebisl;
      }

    ///
    virtual ~EBCellFactory()
      {
      }

  private:
    EBISLayout m_ebisl;

    //disallowed constructors and all that

    EBCellFactory();

    void operator=(const EBCellFactory& a_inputs);


  };
}
#endif
