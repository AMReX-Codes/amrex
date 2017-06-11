/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */

#include "AMReX_EBData.H"
#include "AMReX_VoFIterator.H"
#include "AMReX_FaceIterator.H"
#include "AMReX_EBISBox.H"
#include "AMReX_PolyGeom.H"
#include "AMReX_EBDataVarMacros.H"

namespace amrex
{
  void null_deleter(EBDataImplem *)
 {}
/************************/
  void 
  EBDataImplem::
  addFullIrregularVoFs(const IntVectSet     &    a_vofsToChange,
                       const Box            &    a_region)
  {
    IntVectSet subset = a_vofsToChange;
    subset &= a_region;
    if (!subset.isEmpty())
    {
      //calculate set by adding in new intvects
      const IntVectSet& ivsOld = m_volData.getIVS();
      IntVectSet ivsNew = m_graph.getIrregCells(m_region);
      ivsNew |= a_vofsToChange;

      int nfaccomp = F_FACENUMBER;
      int nvolcomp = V_VOLNUMBER;

      //for copying purposes
      Box minBox = m_region;

      //save old data into temporarys
      BaseIVFAB<Real>  oldVolData(ivsOld, m_graph, nvolcomp);
      oldVolData.copy(m_volData, minBox, 0, minBox, 0, nvolcomp);

      BaseIFFAB<Real> oldFaceData[SpaceDim];
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        oldFaceData[idir].define(ivsOld, m_graph, idir, nfaccomp);
        oldFaceData[idir].copy(m_faceData[idir], minBox, 0, minBox, 0, nfaccomp);
      }

      //redefine member data with new IVS and copy old data back
      //into it.
      m_volData.define(ivsNew, m_graph, 1);
      m_volData.copy(oldVolData, minBox, 0, minBox, 0, nvolcomp);
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        m_faceData[idir].define(ivsNew, m_graph, idir, 1);
        m_faceData[idir].copy(oldFaceData[idir], minBox, 0, minBox, 0, nfaccomp);
      }

      //now put correct data into new vofs.  the volume fraction of a formally
      //regular cell is always unity
      for (VoFIterator vofit(subset, m_graph); vofit.ok(); ++vofit)
      {
        m_volData(vofit(), V_VOLFRAC) = 1;
        m_volData(vofit(), V_BNDAREA) = 0;
        for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_volData(vofit(), V_VOLCENTROIDX + idir) = 0.;
          m_volData(vofit(), V_BNDCENTROIDX + idir) = 0.;
          m_volData(vofit(), V_NORMALX      + idir) = 0.;
        }
        //no real boundary here but putting something valid for the normal
        m_volData(vofit(), V_NORMALX) = 1.;

      }

      //there are rare cases that arise from coarsening where the area fractions will not be one.
      //we can check this by seeing the number of faces on the side
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        FaceIterator faceit(subset, m_graph, idir, FaceStop::SurroundingWithBoundary);
        for (faceit.reset(); faceit.ok(); ++faceit)
        {

          Real areaFrac = 1;
          RealVect faceCentroid = RealVect::Zero;


          if (!faceit().isBoundary())
          {
            IntVect ivhi = faceit().gridIndex(Side::Hi);
            IntVect ivlo = faceit().gridIndex(Side::Lo);
            std::vector<FaceIndex> allFaces = m_graph.getAllFaces(ivhi, idir, Side::Lo);
            if (allFaces.size() > 1)
            {
              //now we have the wacky case where two full, single-valued faces were coarsened
              //to make a multi-valued face on the coarse level.
              //also know as the bjorn geometry.  We need to infer the centroids from the
              //volume centroids of the vofs involved.   The value of the centroid will be
              //0.5 so i just need to get the sign.
              //doing the best i can here.   there might be cases where this gets the
              //centroid wrong but i cannot think of them
              VolIndex otherVoF;
              if (subset.contains(ivlo) && (!subset.contains(ivhi)))
              {
                otherVoF = faceit().getVoF(Side::Hi);
              }
              else if (subset.contains(ivhi) && (!subset.contains(ivlo)))
              {
                otherVoF = faceit().getVoF(Side::Lo);
              }
              else
              {
                //you really should only be changing one of the vofs if there is a multivalued
                //face between them
                amrex::Error("vofsToChange contains internal mutlivalued face");
              }
              faceCentroid[idir] = 0.0;
              for (int tandir = 0; tandir < SpaceDim; tandir++)
              {
                if(tandir != idir)
                {
                  Real tancentroid = m_volData(otherVoF, V_VOLCENTROIDX + tandir);
                  faceCentroid[tandir] = tancentroid;
                }
              }
            }
           }
          m_faceData[idir](faceit(), F_AREAFRAC) = areaFrac;
          for (int idir = 0; idir < SpaceDim; idir++)
          {
            m_faceData[idir](faceit(), F_FACECENTROIDX+idir) = faceCentroid[idir];
          }
        }
      }
    }
  }
/************************/
/************************/
  EBDataImplem::
  EBDataImplem()
  {
    m_isDefined = false;
  }
/************************/
  EBDataImplem::
  ~EBDataImplem()
  {
  }
/************************/
  void 
  EBDataImplem::
  define(const Box& box, int comps)
  {
  }
/************************/
  EBDataImplem::
  EBDataImplem(const Box& a_box, int a_comps)
  {
  }
/************************/
  EBDataImplem& 
  EBDataImplem::
  copy(const EBDataImplem&   a_src,
       const Box&            a_srcbox,
       int                   a_srccomp,
       const Box&            a_dstbox,
       int                   a_dstcomp,
       int                   a_numcomp)
  {
    assert(m_isDefined);
    assert(a_src.m_isDefined);

    m_volData.copy(a_src.m_volData, a_srcbox, 0, a_dstbox, 0,  V_VOLNUMBER);
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_faceData[idir].copy(a_src.m_faceData[idir], a_srcbox, 0, a_dstbox, 0,  F_FACENUMBER);
    }
    computeNormalsAndBoundaryAreas(m_graph, m_region);
    return *this;
  }
/************************/
  void 
  EBDataImplem::
  define(const EBGraph&           a_graph,
         const Box&               a_region)
  {
    m_graph = a_graph;
    m_region = a_region & a_graph.getDomain();
    IntVectSet ivsIrreg = a_graph.getIrregCells(a_region);
    m_isDefined = true;
    m_volData.define(ivsIrreg, a_graph, V_VOLNUMBER);
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_faceData[idir].define(ivsIrreg, a_graph, idir, F_FACENUMBER);
    }
  }
/************************/
  void
  EBDataImplem::
  define(const EBGraph&           a_graph,
         const std::vector<IrregNode>& a_irregGraph,
         const Box&               a_validBox,
         const Box&               a_region)

  {
    BL_PROFILE("EBDataImpem::define");
    m_region = a_region & a_graph.getDomain();
    define( a_graph, a_region);
    if (a_graph.hasIrregular())
    {
      for (int inode = 0; inode < a_irregGraph.size(); inode++)
      {
        const IrregNode& node = a_irregGraph[inode];
        const IntVect& iv = node.m_cell;
        if (a_validBox.contains(iv))
        {
          const int&  cellInd = node.m_cellIndex;
          VolIndex vof(iv, cellInd);

          m_volData(vof, V_VOLFRAC) = node.m_volFrac;
          for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
          {
            m_volData(vof, V_VOLCENTROIDX+faceDir) =   node.m_volCentroid[faceDir];
            m_volData(vof, V_BNDCENTROIDX+faceDir) = node.m_bndryCentroid[faceDir];
          }
          for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
          {
            for (SideIterator sit; sit.ok(); ++sit)
            {
              std::vector<FaceIndex> faces = a_graph.getFaces(vof, faceDir, sit());
              int nodeind = node.index(faceDir, sit());
              std::vector<Real> areaFracs         = node.m_areaFrac[nodeind];
              std::vector<RealVect> faceCentroids = node.m_faceCentroid[nodeind];
              for (int iface = 0; iface < faces.size(); iface++)
              {
                const Real&     areaFracNode     = areaFracs[iface];
                const RealVect& faceCentroidNode = faceCentroids[iface];
                const FaceIndex& face = faces[iface];

                m_faceData[faceDir](face,F_AREAFRAC)  = areaFracNode;
                for(int idir = 0; idir < SpaceDim; idir++)
                {
                  m_faceData[faceDir](face,F_FACECENTROIDX+idir)= faceCentroidNode[idir];
                }
              }
            }
          }
        }
      }
    }
    computeNormalsAndBoundaryAreas(a_graph, a_validBox);
  }
/*******************************/
  const Real& EBDataImplem::volFrac(const VolIndex& a_vof) const
  {
    return m_volData(a_vof, V_VOLFRAC);
  }

/*******************************/
  const Real& EBDataImplem::bndryArea(const VolIndex& a_vof) const
  {
    return m_volData(a_vof, V_BNDAREA);
  }

/*******************************/
  RealVect EBDataImplem::normal(const VolIndex& a_vof) const
  {
    RealVect retval;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      retval[idir] = m_volData(a_vof, V_NORMALX+idir);
    }
    return retval;

  }
/*******************************/
  RealVect EBDataImplem::centroid(const VolIndex& a_vof) const
  {
    RealVect retval;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      retval[idir] = m_volData(a_vof, V_VOLCENTROIDX+idir);
    }
    return retval;
  }
/*******************************/
  RealVect EBDataImplem::bndryCentroid(const VolIndex& a_vof) const
  {
    RealVect retval;

    for(int idir = 0; idir < SpaceDim; idir++)
    {
      retval[idir] = m_volData(a_vof, V_BNDCENTROIDX+idir);
    }
    return retval;
  }
/*******************************/
  RealVect EBDataImplem::centroid(const FaceIndex& a_face) const
  {
    int faceDir = a_face.direction();
    RealVect retval;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      retval[idir] = m_faceData[faceDir](a_face, F_FACECENTROIDX + idir);
    }
    return retval;
  }
/*******************************/
  const Real& EBDataImplem::areaFrac(const FaceIndex& a_face) const
  {
    int faceDir = a_face.direction();
    return m_faceData[faceDir](a_face, F_AREAFRAC);
  }
/*******************************/
  void 
  EBDataImplem::
  computeNormalsAndBoundaryAreas(const EBGraph& a_graph,
                                 const Box&     a_region)
  {
    EBISBox ebisBox;
    EBData tempData(this);
    ebisBox.define(a_graph, tempData);

    if (a_graph.hasIrregular())
    {
      IntVectSet ivsIrreg = a_graph.getIrregCells(a_region);
      for (VoFIterator vofit(ivsIrreg, a_graph); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real bndryArea  =  PolyGeom::bndryArea(vof, ebisBox);
        RealVect normal =  PolyGeom::normal(   vof, ebisBox, bndryArea);

        //copy results to volData
        m_volData(vof,V_BNDAREA) = bndryArea;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          m_volData(vof,V_NORMALX+idir) = normal[idir];
        }
      }
    }
  }
/*******************************/
  void
  EBDataImplem::
  coarsenBoundaryAreaAndNormal(Real&                    a_bndryAreaCoar,
                               RealVect&                a_normalCoar,
                               const std::vector<Real>&      a_bndryAreaFine,
                               const std::vector<RealVect>&  a_normalFine)
  {
    BL_PROFILE("EBDataImplem::coarsenBoundaryAreaAndNormal");

    Real faceCoarsenFactor = AMREX_D_TERM(1.0, * 0.5, *0.5);
    // Real voluCoarsenFactor = AMREX_D_TERM(0.5, * 0.5, *0.5);

    //A^B_xi, C = sum(A^B_xi, F/2^D-1)
    RealVect bndryAreaVec= RealVect::Zero;
    for (int ifine = 0; ifine < a_normalFine.size(); ifine++)
    {
      RealVect normalFine = a_normalFine[ifine];
      Real bndryAreaFine =  a_bndryAreaFine[ifine];
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        bndryAreaVec[idir] += normalFine[idir]*bndryAreaFine;
      }
    }
    bndryAreaVec *= faceCoarsenFactor;
    //A^B,C = ||A^B_xi,C||
    a_bndryAreaCoar = 0.;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_bndryAreaCoar += bndryAreaVec[idir]*bndryAreaVec[idir];
    }
    a_bndryAreaCoar = sqrt(a_bndryAreaCoar);

    //n_xi, C = A^B_xi,C/AB,c
    if (a_bndryAreaCoar > 0.)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        a_normalCoar[idir] = bndryAreaVec[idir]/a_bndryAreaCoar;
      }
    }
    else
    {
      a_normalCoar = RealVect::Zero;
    }
  }


  void EBDataImplem::
  coarsenVoFs(const EBDataImplem&  a_fineEBDataImplem,
              const EBGraph&       a_fineGraph,
              const EBGraph&       a_coarGraph,
              const Box&           a_validRegion)
  {
    BL_PROFILE("EBDataImplem::coarsenVoFs");
    //unlike before, define needs to be called first
    assert(m_isDefined);

    if (a_coarGraph.hasIrregular())
    {
      IntVectSet ivsIrreg = a_coarGraph.getIrregCells(a_validRegion);

      for (VoFIterator vofit(ivsIrreg, a_coarGraph); vofit.ok(); ++vofit)
      {
        BL_PROFILE("EBDataImplem::coarsenVoFs_VoFIterator");
        const VolIndex& vofCoar = vofit();
        std::vector<VolIndex> vofsFine = a_coarGraph.refine(vofCoar);
        int nFine = vofsFine.size();
        std::vector<Real> bndryAreaFine(nFine);
        std::vector<Real> volFracFine(nFine);
        std::vector<int>  phase(nFine);
        std::vector<RealVect> bndryCentroidFine(nFine);
        std::vector<RealVect> volCentroidFine(nFine);
        std::vector<RealVect> normalFine(nFine);

        for (int ifine = 0; ifine < nFine; ifine++)
        {
          BL_PROFILE("EBDataImplem::coarsenVoFs_fine");
          const VolIndex& vofFine =vofsFine[ifine];

          if (a_fineGraph.isIrregular(vofFine.gridIndex()))
          {
            bndryAreaFine[ifine] = a_fineEBDataImplem.bndryArea(vofFine);

            volFracFine[ifine]       = a_fineEBDataImplem.volFrac(vofFine);
            bndryCentroidFine[ifine] = a_fineEBDataImplem.bndryCentroid(vofFine);
            volCentroidFine[ifine]   = a_fineEBDataImplem.centroid(vofFine);
            normalFine[ifine]        = a_fineEBDataImplem.normal(vofFine);
          }
          else
          {
            assert(a_fineGraph.isRegular(vofFine.gridIndex()));
            bndryAreaFine[ifine] = 0.0;
            volFracFine[ifine] = 1.0;
            bndryCentroidFine[ifine] = RealVect::Zero;
            volCentroidFine[ifine] = RealVect::Zero;
            normalFine[ifine] = RealVect::Zero;
          }
        }

        Real volFracCoar, bndryAreaCoar;
        RealVect volCentroidCoar, bndryCentroidCoar, normalCoar;

        coarsenVolFracAndCentroid(volFracCoar,
                                  volCentroidCoar,
                                  volFracFine,
                                  volCentroidFine,
                                  vofsFine,
                                  vofCoar);

        coarsenBoundaryAreaAndNormal(bndryAreaCoar,
                                     normalCoar,
                                     bndryAreaFine,
                                     normalFine);

        coarsenBndryCentroid(bndryCentroidCoar,
                             bndryCentroidFine,
                             bndryAreaFine,
                             vofsFine,
                             vofCoar);

        m_volData(vofCoar, V_VOLFRAC)     = volFracCoar;
        m_volData(vofCoar, V_BNDAREA)     = bndryAreaCoar;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          m_volData(vofCoar, V_VOLCENTROIDX+idir) =   volCentroidCoar[idir];
          m_volData(vofCoar, V_BNDCENTROIDX+idir) = bndryCentroidCoar[idir];
          m_volData(vofCoar, V_NORMALX     +idir) =        normalCoar[idir];
        }
      }
    }
  }
/*******************************/
  void EBDataImplem::
  coarsenFaces(const EBDataImplem& a_fineEBDataImplem,
               const EBGraph&      a_fineGraph,
               const EBGraph&      a_coarGraph,
               const Box&          a_validRegion)
  {
    BL_PROFILE("EBDataImplem::coarsenFaces");

    //unlike before, the define function has to be called first.
    assert(m_isDefined);

    IntVectSet ivsIrreg = a_coarGraph.getIrregCells(a_validRegion);
    Box fineRegion = a_fineGraph.getRegion();
    if (a_coarGraph.hasIrregular())
    {
      for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
      {
        BL_PROFILE("EBDataImplem::coarsenFaces_faceDir");

        FaceIterator faceit(ivsIrreg, a_coarGraph, faceDir,
                            FaceStop::SurroundingWithBoundary);

        for (faceit.reset(); faceit.ok(); ++faceit)
        {
          BL_PROFILE("EBDataImplem::coarsenFaces_FaceIterator");

          const FaceIndex&  faceCoar  = faceit();
          std::vector<FaceIndex> facesFine = a_coarGraph.refine(faceCoar, a_fineGraph);

          std::vector<Real>     areaFracsFine(facesFine.size());
          std::vector<RealVect> centroidsFine(facesFine.size());
          for (int ifine = 0; ifine < facesFine.size(); ifine++)
          {
            BL_PROFILE("EBDataImplem::coarsenFaces_fine");

            const FaceIndex& faceFine = facesFine[ifine];
            IntVect loiv = faceFine.gridIndex(Side::Lo);
            IntVect hiiv = faceFine.gridIndex(Side::Hi);
            if ((fineRegion.contains(loiv) && a_fineGraph.isIrregular(loiv)) ||
                (fineRegion.contains(hiiv) && a_fineGraph.isIrregular(hiiv)))
            {
              areaFracsFine[ifine] = a_fineEBDataImplem.areaFrac(faceFine);
              centroidsFine[ifine] = a_fineEBDataImplem.centroid(faceFine);
            }
            else
            {
              areaFracsFine[ifine] = 1.0;
              centroidsFine[ifine] = RealVect::Zero;
            }
          }
          Real areaFracCoar;
          RealVect centroidCoar;
          coarsenAreaFrac(areaFracCoar, areaFracsFine);

          coarsenFaceCentroid(centroidCoar, centroidsFine,
                              areaFracsFine, facesFine,
                              faceCoar);

          m_faceData[faceDir](faceCoar, F_AREAFRAC)       = areaFracCoar;
          for(int idir = 0; idir < SpaceDim; idir++)
          {
            m_faceData[faceDir](faceCoar, F_FACECENTROIDX+idir) = centroidCoar[idir];
          }
        } //end loop over faces
      } //end loop over face directions
    }
  }
/*******************************/
  void EBDataImplem::
  coarsenFaceCentroid(RealVect&                a_centroidCoar,
                      const std::vector<RealVect>&  a_centroidsFine,
                      const std::vector<Real>&      a_areaFracFine,
                      const std::vector<FaceIndex>& a_facesFine,
                      const FaceIndex&         a_faceCoar)
  {
    BL_PROFILE("EBDataImplem::coarsenFaceCentroid");

    Real totalFaceArea = 0.0;
    for (int ifineFace = 0; ifineFace < a_facesFine.size(); ifineFace++)
    {
      totalFaceArea += a_areaFracFine[ifineFace];
    }

    a_centroidCoar = RealVect::Zero;
    if (totalFaceArea > 0.0)
    {
      const IntVect& coarIV = a_faceCoar.gridIndex(Side::Lo);
      for (int ifineFace = 0; ifineFace < a_facesFine.size(); ifineFace++)
      {
        const FaceIndex& fineFace        = a_facesFine[ifineFace];
        IntVect fineIV = fineFace.gridIndex(Side::Lo);

        RealVect centroidFine = a_centroidsFine[ifineFace];
        centroidFine = fineToCoarseTransform(centroidFine,
                                             coarIV, fineIV);

        const Real&   areaFracFine= a_areaFracFine[ifineFace];
        for (int idir = 0; idir < SpaceDim; idir++)
        {
          a_centroidCoar[idir] +=
            centroidFine[idir]*(areaFracFine/totalFaceArea);
        }
      }
    }
  }
/*******************************/
  void EBDataImplem::
  coarsenAreaFrac(Real& a_areaFracCoar,
                  const std::vector<Real>& a_areaFracFine)
  {
    BL_PROFILE("EBDataImplem::coarsenAreaFrac");
    //this is the factor by which the area of a fine
    //face is smaller than the area of a coarse face.
    Real faceCoarsenFactor = 1.0;
    for (int dir = 0; dir < SpaceDim-1; ++dir)
    {
      faceCoarsenFactor *= 0.5;
    }
    a_areaFracCoar = 0.0;
    for (int ifine = 0; ifine < a_areaFracFine.size(); ifine++)
    {
      a_areaFracCoar += a_areaFracFine[ifine];
    }
    a_areaFracCoar *= faceCoarsenFactor;
  }
/*******************************/
  void
  EBDataImplem::
  coarsenVolFracAndCentroid(Real&                   a_volFracCoar,
                            RealVect&               a_volCentroidCoar,
                            const std::vector<Real>&     a_volFracFine,
                            const std::vector<RealVect>& a_volCentroidFine,
                            const std::vector<VolIndex>& a_fineVoFs,
                            const VolIndex&         a_coarVoF)
  {
    BL_PROFILE("EBDataImplem::coarsenVolFracAndCentroid");

    Real volCoarsenFactor = 1.0;
    for (int dir = 0; dir < SpaceDim; ++dir)
    {
      volCoarsenFactor *= 0.5;
    }

    Real totalVol = 0;
    for (int ifine = 0; ifine < a_fineVoFs.size(); ifine++)
    {
      totalVol  += a_volFracFine[ifine];
    }
    a_volFracCoar = volCoarsenFactor*totalVol;

    a_volCentroidCoar = RealVect::Zero;
    for (int ifine = 0; ifine < a_fineVoFs.size(); ifine++)
    {
      const VolIndex& fineVoF= a_fineVoFs[ifine];
      const IntVect& coarIV = a_coarVoF.gridIndex();
      const IntVect& fineIV = fineVoF.gridIndex();
      RealVect volCentroidFine = a_volCentroidFine[ifine];
      Real volFracFine=  a_volFracFine[ifine];
      volCentroidFine = fineToCoarseTransform(volCentroidFine,
                                              coarIV, fineIV);

      if (totalVol > 0.0)
      {
        for (int idir = 0; idir < SpaceDim; idir++)
        {
          a_volCentroidCoar[idir] +=
            volCentroidFine[idir]*(volFracFine/totalVol);
        }
      }
    }
  }
/*******************************/
  RealVect
  EBDataImplem::
  fineToCoarseTransform(const RealVect& a_finePoint,
                        const IntVect&  a_coarCell,
                        const IntVect&  a_fineCell)
  {
    RealVect retval;
    //assuming nref = 2. make dxf = 1
    Real dxc = 2.0;
    RealVect fineCellLoc, coarCellLoc;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      fineCellLoc[idir] = Real(a_fineCell[idir]) + 0.5;
      coarCellLoc[idir] = dxc*(Real(a_coarCell[idir]) + 0.5);
    }
    retval = a_finePoint+ fineCellLoc;
    retval -= coarCellLoc;
    //to put it into a space where dxc = 1, the normalized answer
    retval /= dxc;
    return retval;
  }
/*******************************/
  void
  EBDataImplem::
  coarsenBndryCentroid(RealVect&               a_bndryCentroidCoar,
                       const std::vector<RealVect>& a_bndryCentroidFine,
                       const std::vector<Real>&     a_bndryAreaFine,
                       const std::vector<VolIndex>& a_fineVoFs,
                       const VolIndex&         a_coarVoF)
  {
    BL_PROFILE("EBDataImplem::coarsenBndryCentroid");

    Real totalArea = 0;
    for (int ifine = 0; ifine < a_fineVoFs.size(); ifine++)
    {
      totalArea += a_bndryAreaFine[ifine];
    }

    a_bndryCentroidCoar = RealVect::Zero;
    for (int ifine = 0; ifine < a_fineVoFs.size(); ifine++)
    {
      const VolIndex& fineVoF= a_fineVoFs[ifine];
      const IntVect& fineIV=  fineVoF.gridIndex();
      const IntVect& coarIV=a_coarVoF.gridIndex();

      Real bndryAreaFine = a_bndryAreaFine[ifine];

      RealVect bndryCentroidFine = a_bndryCentroidFine[ifine];

      bndryCentroidFine =
        fineToCoarseTransform(bndryCentroidFine,coarIV, fineIV);
      if (totalArea > 0.0)
      {
        for (int idir = 0; idir < SpaceDim; idir++)
        {
          a_bndryCentroidCoar[idir] +=
            bndryCentroidFine[idir]*(bndryAreaFine/totalArea);
        }
      }
    }
  }
/******************/
  /// Below lies serialization land.  Enter at thy own risk.
  /// Management is not responsible for any gibbering madness resulting 
  /// from ignoring this warning.
/*******************************/
  std::size_t 
  EBDataImplem::
  nBytes (const Box& bx, int a_srccomp, int a_ncomps) const
  {
    size_t retval = m_volData.nBytes(bx, 0, V_VOLNUMBER);
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      retval += m_faceData[idir].nBytes(bx, 0, F_FACENUMBER);
    }
    return retval;
  }
/*******************************/

  std::size_t 
  EBDataImplem::
  copyToMem (const Box& bx,
             int        srccomp,
             int        ncomps,
             void*      dst) const
  {
    size_t retval = 0;
    unsigned char* buf = (unsigned char*) dst;
    size_t incrval = m_volData.copyToMem(bx, 0, V_VOLNUMBER, buf);
    retval += incrval;
    buf    += incrval;

    for(int idir = 0; idir < SpaceDim; idir++)
    {
      incrval = m_faceData[idir].copyToMem(bx, 0, F_FACENUMBER, buf);
      retval += incrval;
      buf    += incrval;
    }
    return retval;
  }

/*******************************/

  std::size_t 
  EBDataImplem::
  copyFromMem (const Box&  bx,
               int         dstcomp,
               int         numcomp,
               const void* src)
  {
    size_t retval = 0;
    unsigned char* buf = (unsigned char*) src;
    size_t incrval = m_volData.copyFromMem(bx, 0, V_VOLNUMBER, buf);
    retval += incrval;
    buf    += incrval;

    for(int idir = 0; idir < SpaceDim; idir++)
    {
      incrval = m_faceData[idir].copyFromMem(bx, 0, F_FACENUMBER, buf);
      retval += incrval;
      buf    += incrval;
    }
    return retval;
  }
/*******************************/
  std::size_t 
  EBDataImplem::
  nBytesFull() const
  {
    //first the region
    size_t retval = m_region.linearSize();
    retval +=   m_graph.nBytesFull();
    retval += m_volData.nBytesFull();
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      retval += m_faceData[idir].nBytesFull();
    }
    return retval;
  }
/*******************************/

  std::size_t 
  EBDataImplem::
  copyToMemFull(void* dst) const
  {
    size_t retval  = 0;
    size_t incrval = 0;
    unsigned char* buf = (unsigned char*) dst;
    m_region.linearOut(buf);
    incrval = m_region.linearSize();
    buf += incrval;
    retval += incrval;

    incrval = m_graph.copyToMemFull(buf);
    retval += incrval;
    buf    += incrval;

    incrval = m_volData.copyToMemFull(buf);
    retval += incrval;
    buf    += incrval;

    for(int idir = 0; idir < SpaceDim; idir++)
    {
      incrval = m_faceData[idir].copyToMemFull(buf);
      retval += incrval;
      buf    += incrval;
    }
    return retval;
  }

/*******************************/

  std::size_t 
  EBDataImplem::
  copyFromMemFull(const void* src)
  {
    
    size_t retval  = 0;
    size_t incrval = 0;
    unsigned char* buf = (unsigned char*) src;
    m_region.linearIn(buf);
    incrval = m_region.linearSize();
    buf += incrval;
    retval += incrval;

    incrval = m_graph.copyFromMemFull(buf);
    retval += incrval;
    buf    += incrval;

    incrval = m_volData.copyFromMemFull(buf);
    retval += incrval;
    buf    += incrval;

    for(int idir = 0; idir < SpaceDim; idir++)
    {
      incrval = m_faceData[idir].copyFromMemFull(buf);
      retval += incrval;
      buf    += incrval;
    }

    m_isDefined = true;
    return retval;
  }

/*******************************/
}
