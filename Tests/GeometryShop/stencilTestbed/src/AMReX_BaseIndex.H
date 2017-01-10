

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

#ifndef _BASEINDEX_H_
#define _BASEINDEX_H_

#include "AMReX_REAL.H"
#include "AMReX_IntVect.H"

namespace amrex
{
  /// index for other indicies to inherit
  /**
     Aggstencil needed the concept of a unified
     index.   This is it.  Not much to it.
  */
  class BaseIndex
  {
  public:

    BaseIndex()
    {
      ;
    }

    virtual ~BaseIndex()
    {
      ;
    }
  };
}
#endif
