

/*
 *      .o.       ooo        ooooo ooooooooo.             ooooooo  ooooo 
 *     .888.      `88.       .888' `888   `Y88.            `8888    d8'  
 *    .8"888.      888b     d'888   888   .d88'  .ooooo.     Y888..8P    
 *   .8' `888.     8 Y88. .P  888   888ooo88P'  d88' `88b     `8888'     
 *  .88ooo8888.    8  `888'   888   888`88b.    888ooo888    .8PY888.    
 * .8'     `888.   8    Y     888   888  `88b.  888    .o   d8'  `888b   
 *o88o     o8888o o8o        o888o o888o  o888o `Y8bod8P' o888o  o88888o 
 *
 */



#ifndef _EBCELLFAB_H_
#define _EBCELLFAB_H_

#include <cmath>
#include <cstdlib>
#include "AMReX_FArrayBox.H"
#include "AMReX_Box.H"
#include "AMReX_BaseIndex.H"
#include "AMReX_VolIndex.H"

namespace amrex
{
  ///
  /**
     EBCellFAB is a BaseEBCellFAB<Real>
     which includes arithmetic functions.
  */
  class EBCellFAB: public FArrayBox
  {
  public:

    ///
    EBCellFAB(): FArrayBox()
    {
    }

    ///
    /**
       Box going into this needs to be cell-centered.
       Data will exist over the surrounding nodes of the box.
    */
    EBCellFAB(const Box& a_region, int a_nVar):
      FArrayBox(a_region, a_nVar)
    {
    }
    

    //for aggstencil --this version of EBCellFAB is just a wrapper around FAB so
    // only one datatype
    int numDataTypes() const
    {
      return 1;
    }

    /// for AggStencil
    int   dataType(const BaseIndex& a_baseInd) const
    {
      return 0;
    }


    ///for AggStencil
    Real* dataPtr(int a_dataType, int a_ivar)
    {
      return FArrayBox::dataPtr(a_ivar);
    }

    ///for AggStencil
    const Real* dataPtr(int a_dataType, int a_ivar) const
    {
      return FArrayBox::dataPtr(a_ivar);
    }


    const Real* dataPtr(int a_ivar) const
    {
      return FArrayBox::dataPtr(a_ivar);
    }


    Real* dataPtr(int a_ivar)
    {
      return FArrayBox::dataPtr(a_ivar);
    }

    virtual long offset(const BaseIndex& a_baseInd, const int a_ivar) const;

    ///
    virtual ~EBCellFAB()
    {
    }


    const Real& operator()(const VolIndex& a_ndin, int a_nvarLoc) const;


    Real& operator()(const VolIndex& a_ndin, int a_nvarLoc) ;


    


  private:
    //disallowed for all the usual reasons
    EBCellFAB(const EBCellFAB& ebcin);

    void operator=(const EBCellFAB& fabin);


  };


}
#endif
