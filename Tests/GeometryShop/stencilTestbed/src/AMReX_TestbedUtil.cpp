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
#include "AMReX_BoxIterator.H"
#include "AMReX_Stencils.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_AggStencil.H"
#include "lapl_nd_F.H"
namespace amrex
{
  void 
  TestbedUtil::
  applyStencilPointwise(EBCellFAB                       & a_dst,
                        const EBCellFAB                 & a_src,
                        const BaseFab<VoFStencil>       & a_stencil,
                        const BaseFab<int>              & a_regIrregCovered,
                        const std::vector<IrregNode>    & a_nodes,
                        const Box                       & a_domain,
                        const Real                      & a_dx)
  { 
    for(BoxIterator boxit(a_domain); boxit.ok(); ++boxit)
    {
      const IntVect& iv = boxit();
      //skip covered cells because they are not part of the solution domain
      if(a_regIrregCovered(iv, 0) >= 0)
      {
        const VoFStencil& sten = a_stencil(iv, 0);
        Real dstVal = applyVoFStencil(sten, a_src);
        VolIndex dstVoF(iv,0);
        int ivar = 0;
        a_dst(dstVoF, ivar) = dstVal;
      }
    }
  }


  ///apply a stencil using aggstencil everywhere
  void 
  TestbedUtil::
  applyStencilAllAggSten(EBCellFAB                       & a_dst,
                         const EBCellFAB                 & a_src,
                         const BaseFab<VoFStencil>       & a_stencil,
                         const BaseFab<int>              & a_regIrregCovered,
                         const std::vector<IrregNode>    & a_nodes,
                         const Box                       & a_domain,
                         const Real                      & a_dx)
  { 
    Vector<std::shared_ptr<BaseIndex  > > dstVoFs;
    Vector<std::shared_ptr<BaseStencil> > vofStencils;

    BL_PROFILE_VAR("all aggsten prep",aap);
    for(BoxIterator boxit(a_domain); boxit.ok(); ++boxit)
    {
      const IntVect& iv = boxit();
      //exclude covered cells
      if(a_regIrregCovered(iv, 0) >= 0)
      {
        std::shared_ptr<BaseIndex>    vofptr(new VolIndex(iv, 0));
        std::shared_ptr<BaseStencil> stenptr(new VoFStencil(a_stencil(iv, 0)));
        dstVoFs.push_back(vofptr);
        vofStencils.push_back(stenptr);
      }
    }
    BL_PROFILE_VAR_STOP(aap);
	
    //define and apply have internal timers 
    BL_PROFILE_VAR("all aggsten ctr",aac);
    AggStencil<EBCellFAB, EBCellFAB> sten(dstVoFs, vofStencils, a_src, a_dst);
    BL_PROFILE_VAR_STOP(aac);

    BL_PROFILE_VAR("all aggsten apply",aaa);
    int isrc = 0; int idst = 0; int inco = a_dst.nComp(); bool incrOnly = false;
    sten.apply(a_dst, a_src, isrc, idst, inco, incrOnly);
    BL_PROFILE_VAR_STOP(aaa);

  }
  ///apply a stencil using fortran on regular cells, pointwise irregular
  void 
  TestbedUtil::
  applyStencilFortranPlusPointwise(EBCellFAB                       & a_dst,
                                   const EBCellFAB                 & a_src,
                                   const BaseFab<VoFStencil>       & a_stencil,
                                   const BaseFab<int>              & a_regIrregCovered,
                                   const std::vector<IrregNode>    & a_nodes,
                                   const Box                       & a_domain,
                                   const Real                      & a_dx)
  { 
    BL_PROFILE_VAR("fortran plus pointwise reg",fppr);
    const BaseFab<Real>& srcReg = a_src.getSingleValuedFAB();
    BaseFab<Real>      & dstReg = a_dst.getSingleValuedFAB();
    fort_lapl_simple(BL_TO_FORTRAN_N_3D(dstReg, 0), 
                     BL_TO_FORTRAN_N_3D(srcReg, 0), 
                     ARLIM_3D(a_domain.loVect()),
                     ARLIM_3D(a_domain.hiVect()), &a_dx);
    BL_PROFILE_VAR_STOP(fppr);


    BL_PROFILE_VAR("fortran plus pointwise eb",fppe);
    for(BoxIterator boxit(a_domain); boxit.ok(); ++boxit)
    {
      const IntVect& iv = boxit();
      //only applying at irregular cells
      if(a_regIrregCovered(iv, 0) == 0)
      {
        const VoFStencil& sten = a_stencil(iv, 0);
        Real dstVal = applyVoFStencil(sten, a_src);
        VolIndex dstVoF(iv,0);
        int ivar = 0;
        a_dst(dstVoF, ivar) = dstVal;
      }
    }
    BL_PROFILE_VAR_STOP(fppe);
  }




  ///apply a stencil using aggstencil on irregular cells, fortran otherwise
  void 
  TestbedUtil::
  applyStencilFortranPlusAggSten(EBCellFAB                       & a_dst,
                                 const EBCellFAB                 & a_src,
                                 const BaseFab<VoFStencil>       & a_stencil,
                                 const BaseFab<int>              & a_regIrregCovered,
                                 const std::vector<IrregNode>    & a_nodes,
                                 const Box                       & a_domain,
                                 const Real                      & a_dx)

  { 
    BL_PROFILE_VAR("fortran plus aggsten reg",fpar);
    const BaseFab<Real>& srcReg = a_src.getSingleValuedFAB();
    BaseFab<Real>      & dstReg = a_dst.getSingleValuedFAB();
    fort_lapl_simple(BL_TO_FORTRAN_N_3D(dstReg, 0), 
                     BL_TO_FORTRAN_N_3D(srcReg, 0), 
                     ARLIM_3D(a_domain.loVect()),
                     ARLIM_3D(a_domain.hiVect()), &a_dx);
    BL_PROFILE_VAR_STOP(fpar);

    Vector<std::shared_ptr<BaseIndex  > > dstVoFs;
    Vector<std::shared_ptr<BaseStencil> > vofStencils;

    BL_PROFILE_VAR("fortran plus aggsten agg prep",fpaap);
    for(BoxIterator boxit(a_domain); boxit.ok(); ++boxit)
    {
      const IntVect& iv = boxit();
      //only using aggsten for irregular cells
      if(a_regIrregCovered(iv, 0) == 0)
      {
        std::shared_ptr<BaseIndex>    vofptr(new VolIndex(iv, 0));
        std::shared_ptr<BaseStencil> stenptr(new VoFStencil(a_stencil(iv, 0)));
        dstVoFs.push_back(vofptr);
        vofStencils.push_back(stenptr);
      }
    }
    BL_PROFILE_VAR_STOP(fpaap);

    BL_PROFILE_VAR("fortran plus aggsten agg build",fpaab);
    AggStencil<EBCellFAB, EBCellFAB> sten(dstVoFs, vofStencils, a_src, a_dst);
    BL_PROFILE_VAR_STOP(fpaab);
    
    BL_PROFILE_VAR("fortran plus aggsten agg apply",fpaaa);
    int isrc = 0; int idst = 0; int inco = a_dst.nComp(); bool incrOnly = false;
    sten.apply(a_dst, a_src, isrc, idst, inco, incrOnly);
    BL_PROFILE_VAR_STOP(fpaaa);
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

    FaceStencil retval;
    int faceDir = a_face.direction();
    std::vector<RealVect> centroids = a_node.m_faceCentroid[a_arcIndex];
    //here i am ruling out both multivalued and pointing into a covered cell
    if(centroids.size() != 1)
    {
      Abort("face not found");
    }
    const RealVect& centroid = centroids[0];
    //need the centroids in the non-face direction (and the directions themselves)
    int   tanDirs[SpaceDim-1];
    Real tanCents[SpaceDim-1];
    int index = 0;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      if(faceDir != idir)
      {
        tanDirs[index] = idir;
        tanCents[index] = centroid[idir];
        index++;
      }
    }
    //all this stuff ignnores the connectivity walk and just uses box calculus
    //this will fail in some cases (where the EBArith method will not)
    if(SpaceDim == 2)
    {
      IntVect ivshift;
      if(tanCents[0] > 0)
      {
        ivshift = BASISV(tanDirs[0]);
      }
      else
      {
        ivshift = -BASISV(tanDirs[0]);
      }
      IntVect ivlo = a_face.gridIndex(Side::Lo) + ivshift;
      IntVect ivhi = a_face.gridIndex(Side::Hi) + ivshift;
      if((  a_domain.contains(ivlo) && (a_regIrregCovered(ivlo, 0) >= 0))
         &&(a_domain.contains(ivhi) && (a_regIrregCovered(ivhi, 0) >= 0)))
      {
        FaceIndex otherFace(VolIndex(ivlo, 0), VolIndex(ivhi, 0));
        Real dist = std::abs(tanCents[0]);
        Real thisWeight = 1.0 - dist;
        Real otherWeight =  dist;
        retval.add(a_face, thisWeight);
        retval.add(otherFace, otherWeight);
      }
      else
      {
        //the other face is not there
        retval.add(a_face, 1.0);
      }

    }
    else if(SpaceDim == 3)
    {
      FaceIndex otherface = a_face;
      IntVect ivshift[2];
      for(int ishift = 0; ishift < 2; ishift++)
      {
        if(tanCents[ishift] > 0)
        {
          ivshift[ishift] = BASISV(tanDirs[ishift]);
        }
        else
        {
          ivshift[ishift] = -BASISV(tanDirs[ishift]);
        }
      }
      IntVect ivlo00 = a_face.gridIndex(Side::Lo);
      IntVect ivhi00 = a_face.gridIndex(Side::Hi);
      IntVect ivlo10 = ivlo00 + ivshift[0];
      IntVect ivhi10 = ivhi00 + ivshift[0];
      IntVect ivlo01 = ivlo00 + ivshift[1];
      IntVect ivhi01 = ivhi00 + ivshift[1];
      IntVect ivlo11 = ivlo00 + ivshift[0]+ ivshift[1];
      IntVect ivhi11 = ivhi00 + ivshift[0]+ ivshift[1];
      if(
        (  a_domain.contains(ivlo01) && (a_regIrregCovered(ivlo01, 0) >= 0)) 
        &&(a_domain.contains(ivhi01) && (a_regIrregCovered(ivhi01, 0) >= 0)) &&
        (  a_domain.contains(ivlo10) && (a_regIrregCovered(ivlo10, 0) >= 0)) 
        &&(a_domain.contains(ivhi10) && (a_regIrregCovered(ivhi10, 0) >= 0)) &&
        (  a_domain.contains(ivlo11) && (a_regIrregCovered(ivlo11, 0) >= 0)) 
        &&(a_domain.contains(ivhi11) && (a_regIrregCovered(ivhi11, 0) >= 0))
        )
      {
        //if we have all the faces, do the second order thing.
        FaceIndex face00(VolIndex(ivlo00, 0), VolIndex(ivhi00,0));
        FaceIndex face01(VolIndex(ivlo01, 0), VolIndex(ivhi01,0));
        FaceIndex face10(VolIndex(ivlo10, 0), VolIndex(ivhi10,0));
        FaceIndex face11(VolIndex(ivlo11, 0), VolIndex(ivhi11,0));

        Real xbar = std::abs(tanCents[0]);
        Real ybar = std::abs(tanCents[1]);
        Real f00coef = 1.0 - xbar - ybar + xbar*ybar;
        Real f10coef = xbar - xbar*ybar;
        Real f01coef = ybar - xbar*ybar;
        Real f11coef = xbar*ybar;

        retval.add(face00, f00coef);
        retval.add(face01, f01coef);
        retval.add(face10, f10coef);
        retval.add(face11, f11coef);
      }
      else
      {
        //at least one of the other faces is not there
        retval.add(a_face, 1.0);
      }
    }
    else
    {
      Abort("bogus spacedim");
    }
    
    return retval;
  }
}                                        

