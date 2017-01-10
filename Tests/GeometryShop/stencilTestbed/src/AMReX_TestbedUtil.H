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


#ifndef _TestbedUtil_H_
#define _TestbedUtil_H_

#include "AMReX_REAL.H"
#include "AMReX_IntVect.H"
#include "AMReX_Stencils.H"
#include "AMReX_BaseIndex.H"
#include "AMReX_IrregNode.H"
#include "AMReX_BaseFab.H"

namespace amrex
{
  /// collection of useful functions
  class TestbedUtil
  {
  public:
    ///
    TestbedUtil(){;}
    ~TestbedUtil(){;}

    ///apply a stencil in the slowest possible way
    static void applyStencilPointwise(EBCellFAB                       & a_dst,
                                      const EBCellFAB                 & a_src,
                                      const BaseFab<VoFStencil>       & a_stencil,
                                      const BaseFab<int>              & a_regIrregCovered,
                                      const std::vector<IrregNode>    & a_nodes,
                                      const Box                       & a_domain,
                                      const Real                      & a_dx);


    ///apply a stencil using fortran on regular cells, pointwise irregular
    static void applyStencilFortranPlusPointwise(EBCellFAB                       & a_dst,
                                                 const EBCellFAB                 & a_src,
                                                 const BaseFab<VoFStencil>       & a_stencil,
                                                 const BaseFab<int>              & a_regIrregCovered,
                                                 const std::vector<IrregNode>    & a_nodes,
                                                 const Box                       & a_domain,
                                                 const Real                      & a_dx);


    ///apply a stencil using aggstencil everywhere
    static void applyStencilAllAggSten(EBCellFAB                       & a_dst,
                                       const EBCellFAB                 & a_src,
                                       const BaseFab<VoFStencil>       & a_stencil,
                                       const BaseFab<int>              & a_regIrregCovered,
                                       const std::vector<IrregNode>    & a_nodes,
                                       const Box                       & a_domain,
                                       const Real                      & a_dx);


    ///apply a stencil using aggstencil on irregular cells, fortran otherwise
    static void applyStencilFortranPlusAggSten(EBCellFAB                       & a_dst,
                                               const EBCellFAB                 & a_src,
                                               const BaseFab<VoFStencil>       & a_stencil,
                                               const BaseFab<int>              & a_regIrregCovered,
                                               const std::vector<IrregNode>    & a_nodes,
                                               const Box                       & a_domain,
                                               const Real                      & a_dx);

    //get the face stencil that goes from face centered fluxes  to centroid fluxes
    static FaceStencil getInterpStencil(const FaceIndex     & a_face,
                                        const IrregNode     & a_node, 
                                        const int           & a_arcIndex,
                                        const BaseFab<int>  & a_regIrregCovered,
                                        const Box           & a_domain,
                                        const Real          & a_dx);
                                        




  };
}

#endif
