

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


#include "AMReX_EBCellFAB.H"


namespace amrex
{
  long
  EBCellFAB::
  offset(const BaseIndex& a_baseInd, const int a_ivar) const
  {
    const VolIndex* vofPtr = dynamic_cast<const VolIndex *>(&a_baseInd);
    if (vofPtr == NULL) Abort("Trying to index into an EBCellFAB with something that is not a VolIndex");
    long retval = 0;

    const IntVect& iv = vofPtr->gridIndex();
    IntVect ivDiff = iv - domain.smallEnd();
    IntVect ivSize = domain.size();
    retval = ivDiff[0];

#if BL_SPACEDIM > 1
    retval += ivDiff[1]*ivSize[0] ;
#endif
#if BL_SPACEDIM > 2
    retval += ivDiff[2]*ivSize[0]*ivSize[1];
#endif

    retval += numpts * a_ivar;

    return retval;
  }
  

  const Real& 
  EBCellFAB::
  operator() (const VolIndex& a_ndin, int a_nvarLoc) const
  {
    FArrayBox& fabcast = (FArrayBox&)(*this);
    return fabcast(a_ndin.gridIndex(),  a_nvarLoc);
  }


  Real& 
  EBCellFAB::
  operator() (const VolIndex& a_ndin, int a_nvarLoc) 
  {
    FArrayBox& fabcast = (FArrayBox&)(*this);
    return fabcast(a_ndin.gridIndex(),  a_nvarLoc);
  }


}

