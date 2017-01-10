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

#include "AMReX_TestbedUtil.H"

namespace amrex
{
  void 
  TestbedUtil::
  applyStencilPointwise(EBCellFAB                  & a_dst,
                        const EBCellFAB            & a_src,
                        const BaseFab<VoFStencil>  & a_stencil,
                        const BaseFab<int>         & a_regIrregCovered,
                        const Vector<IrregNode>    & a_nodes,
                        const Box                  & a_domain,
                        const Real                 & a_dx)
  { 
    BL_Abort("not implemented ");
  }


  ///apply a stencil using fortran on regular cells, pointwise irregular
  void 
  TestbedUtil::
  applyStencilFortranPlusPointwise(EBCellFAB                  & a_dst,
                                   const EBCellFAB            & a_src,
                                   const BaseFab<VoFStencil>  & a_stencil,
                                   const BaseFab<int>         & a_regIrregCovered,
                                   const Vector<IrregNode>    & a_nodes,
                                   const Box                  & a_domain,
                                   const Real                 & a_dx)
  { 
    BL_Abort("not implemented ");
  }


  ///apply a stencil using aggstencil everywhere
  void 
  TestbedUtil::
  applyStencilAllAggSten(EBCellFAB                  & a_dst,
                         const EBCellFAB            & a_src,
                         const BaseFab<VoFStencil>  & a_stencil,
                         const BaseFab<int>         & a_regIrregCovered,
                         const Vector<IrregNode>    & a_nodes,
                         const Box                  & a_domain,
                         const Real                 & a_dx)
  { 
    BL_Abort("not implemented ");
  }


  ///apply a stencil using aggstencil on irregular cells, fortran otherwise
  void 
  TestbedUtil::
  applyStencilFortranPlusAggSten(EBCellFAB                  & a_dst,
                                 const EBCellFAB            & a_src,
                                 const BaseFab<VoFStencil>  & a_stencil,
                                 const BaseFab<int>         & a_regIrregCovered,
                                 const Vector<IrregNode>    & a_nodes,
                                 const Box                  & a_domain,
                                 const Real                 & a_dx)
  { 
    BL_Abort("not implemented ");
  }

  //get the face stencil that goes from face centered fluxes  to centroid fluxes
  FaceStencil 
  TestbedUtil::
  getInterpStencil(const FaceIndex     & a_face,
                   const IrregNode     & a_node, 
                   const int           & a_arcIndex,
                   const BaseFab<int>  & a_regIrregCovered,
                   const Box           & a_domain,
                   const Real          & a_dx)
  { 
    BL_Abort("not implemented ");
    FaceStencil retval;
    return retval;
  }
}                                        

