#ifndef _BASEEBCELLFACTORY_H___
#define _BASEEBCELLFACTORY_H___

#include <cmath>
#include <cstdlib>
#include "AMReX_EBISLayout.H"
#include "AMReX_BaseEBCellFAB.H"


namespace amrex
{

/// Factory class to produce BaseEBCellFABs.
/**
   Factory class to produce BaseEBCellFABs.
   This is needed for LevelData construction.
*/
  template <class T>
  class BaseEBCellFactory
    : public FabFactory<BaseEBCellFAB<T> >
  {
  public:
    /// factory function.
    /**
       Creates a new baseebcellfab object
       and returns a pointer to it.  Responsiblitly
       for calling operator 'delete' on this pointer is passed to the user.
    */
    virtual BaseEBCellFAB<T>* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const override
      {
        EBISBox ebisBox = m_ebisl[box_index];
        return new BaseEBCellFAB<T>(ebisBox, box, ncomps);
      }

    virtual BaseEBCellFactory<T>* clone () const override {
        return new BaseEBCellFactory<T>(*this);
    }

    BaseEBCellFactory(const BaseEBCellFactory& a_inputs) = default;

    ///
    /**
       create the factory with an ispace.  calls full define function
    */
    BaseEBCellFactory(const EBISLayout& a_ebisl)
      {
        m_ebisl = a_ebisl;
      }

    ///
    virtual ~BaseEBCellFactory()
      {
      }

  private:
    EBISLayout m_ebisl;

    BaseEBCellFactory();

    void operator=(const BaseEBCellFactory& a_inputs);

  };

}

#endif
