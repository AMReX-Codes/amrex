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


#ifndef _VOLINDEX_H_
#define _VOLINDEX_H_

#include "AMReX_REAL.H"
#include "AMReX_IntVect.H"
#include "AMReX_BaseIndex.H"

namespace amrex
{
  /// Volume of Fluid Index
  /**
     VolIndex is a very lightweight object
     to be used to distinguish between volumes of
     fluid in a cell.  We use VolIndex for
     array indexing.   Non-trivial constructors
     are all private because only EBISBox
     is supposed to make valid volume indicies.
  */
  class VolIndex: public BaseIndex
  {
  public:
    ///
    VolIndex();

    virtual ~VolIndex();

    ///
    //VolIndex& operator= (const VolIndex& a_vofin);

    ///
    inline
    bool operator== (const VolIndex& a_vofin) const;

    ///
    inline
    const IntVect& gridIndex() const;

    ///
    inline
    int cellIndex() const;

    ///
    inline
    bool isDefined() const;

    ///
    // VolIndex(const VolIndex& a_vofin);

    ///
    void define(const VolIndex& a_vofin);

    ///
    inline
    VolIndex(const IntVect& a_ix,const  int& a_vofID);

    ///
    inline
    void define(const IntVect& a_ix,const  int& a_vofID);

    inline
    bool operator<(const VolIndex& rhs) const;

    inline
    bool operator!=(const VolIndex& rhs) const;

    ///for stencil gymnastics
    void shift(const IntVect& a_iv)
    {
      m_iv += a_iv;
    }

  private:

    //cell of the VolIndex
    IntVect m_iv;

    //which VoF in the cell (0 to nVoF-1)
    int m_cellIndex;

    bool m_isDefined;
  };


}

#endif
