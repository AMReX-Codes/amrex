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
#include "AMReX_EBDataVarMacros.H"

namespace amrex
{
/************************/
  void 
  EBData::
  addFullIrregularVoFs(const IntVectSet& a_vofsToChange,
                       const EBGraph&    a_newGraph,
                       const BaseIVFAB<VolData>& a_grownData,
                       const EBGraph&    a_oldGraph)
  {
    m_implem->addFullIrregularVoFs(a_vofsToChange, a_newGraph, a_grownData, a_oldGraph);
  }
/************************/
  void
  EBData::
  addEmptyIrregularVoFs(const IntVectSet& a_vofsToChange,
                        const EBGraph&    a_newGraph)
  {
    m_implem->addEmptyIrregularVoFs(a_vofsToChange, a_newGraph);
  }
/************************/
  void 
  EBDataImplem::
  addFullIrregularVoFs(const IntVectSet     &    a_vofsToChange,
                       const EBGraph        &    a_newGraph,
                       const BaseIVFAB<Real>&    a_grownData,
                       const EBGraph        &    a_oldGraph)
  {
    if (!a_vofsToChange.isEmpty())
    {
      //calculate set by adding in new intvects
      const IntVectSet& ivsOld = m_volData.getIVS();
      IntVectSet ivsNew = ivsOld | a_vofsToChange;
      int nfaccomp = F_FACENUMBER;
      int nvolcomp = V_VOLNUMBER;

      //for copying purposes
      Box minBox = ivsOld.minBox();
      Interval interv(0,0);

      //save old data into temporarys
      BaseIVFAB<Real>  oldVolData(ivsOld, a_oldGraph, nvolcomp);
      oldVolData.copy(m_volData, minBox, 0, minBox, 0, nvolcomp);

      BaseIFFAB<Real> oldFaceData[SpaceDim];
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        oldFaceData[idir].define(ivsOld, a_oldGraph, idir, nfaccomp);
        oldFaceData[idir].copy(m_faceData[idir], minBox, 0, minBox, 0, nfaccomp);
      }

      //redefine member data with new IVS and copy old data back
      //into it.
      m_volData.define(ivsNew, a_newGraph, 1);
      m_volData.copy(oldVolData, minBox, 0, minBox, 0, nvolcomp);
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        m_faceData[idir].define(ivsNew, a_newGraph, idir, 1);
        m_faceData[idir].copy(oldFaceData[idir], minBox, 0, minBox, 0, nfaccomp);
      }

      //now put correct data into new vofs.  the volume fraction of a formally
      //regular cell is always unity
      for (VoFIterator vofit(a_vofsToChange, a_newGraph); vofit.ok(); ++vofit)
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
        FaceIterator faceit(a_vofsToChange, a_newGraph, idir, FaceStop::SurroundingWithBoundary);
        for (faceit.reset(); faceit.ok(); ++faceit)
        {

          Real areaFrac = 1;
          RealVect faceCentroid = RealVect::Zero;


          if (!faceit().isBoundary())
          {
            IntVect ivhi = faceit().gridIndex(Side::Hi);
            IntVect ivlo = faceit().gridIndex(Side::Lo);
            Vector<FaceIndex> allFaces = a_newGraph.getAllFaces(ivhi, idir, Side::Lo);
            if (allFaces.size() > 1)
            {
              //now we have the wacky case where two full, single-valued faces were coarsened
              //to make a multi-valued face on the coarse level.
              //also know as the bjorn geometry.  We need to infer the centroids from the
              //volume centroids of the vofs involved.   The value of the centroid will be
              //0.5 so i just need to get the sign.
              //doing the best i can here.   there might be cases where this gets the
              //centroid wrong but i cannot think of them
              m_areaFrac = 1.0/Real(allFaces.size());
              VolIndex otherVoF;
              if (a_vofsToChange.contains(ivlo) && (!a_vofsToChange.contains(ivhi)))
              {
                otherVoF = faceit().getVoF(Side::Hi);
              }
              else if (a_vofsToChange.contains(ivhi) && (!a_vofsToChange.contains(ivlo)))
              {
                otherVoF = faceit().getVoF(Side::Lo);
              }
              else
              {
                //you really should only be changing one of the vofs if there is a multivalued
                //face between them
                MayDay::Error("vofsToChange contains internal mutlivalued face");
              }
              for (int tandir = 0; tandir < SpaceDim; tandir++)
              {
                faceCentroid[tandir] = 1.0/Real(allFaces.size());
                Real tancentroid = a_grownData(otherVoF, V_VOLCENTROIDX + tandir);
                if (tancentroid  < 0)
                {
                  fullFace.m_faceCentroid[tandir] *= -1;
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
  void EBDataImplem::
  addEmptyIrregularVoFs(const IntVectSet& a_vofsToChange,
                        const EBGraph&    a_newGraph)
  {
    if (!a_vofsToChange.isEmpty())
    {
      //calculate set by adding in new intvects
      const IntVectSet& ivsOld = m_volData.getIVS();
      IntVectSet ivsNew = ivsOld | a_vofsToChange;

      //for copying purposes
      Box minBox = ivsOld.minBox();
      Interval interv(0,0);

      //save old data into temporarys
      BaseIVFAB<VolData>  oldVolData(ivsOld, a_newGraph, 1);
      oldVolData.copy(minBox, interv, minBox, m_volData, interv);
      BaseIFFAB<FaceData> oldFaceData[SpaceDim];
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        oldFaceData[idir].define(ivsOld, a_newGraph, idir, 1);
        oldFaceData[idir].copy(minBox, interv, minBox, m_faceData[idir], interv);
      }

      //redefine member data with new IVS and copy old data back
      //into it.
      m_volData.define(ivsNew, a_newGraph, 1);
      m_volData.copy(minBox, interv, minBox, oldVolData, interv);
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        m_faceData[idir].define(ivsNew, a_newGraph, idir, 1);
        m_faceData[idir].copy(minBox, interv, minBox, oldFaceData[idir], interv);
      }

      //now put correct data into new vofs.  the volume fraction of a formally
      //regular cell is always unity
      for (VoFIterator vofit(a_vofsToChange, a_newGraph); vofit.ok(); ++vofit)
      {
        m_volData(vofit(), V_VOLFRAC) = 0;
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
    }
  }
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
       int                   a_srccomp
       const Box&            a_dstbox,
       int                   a_dstcomp
       int                   a_numcomp)
  {
    assert(m_isDefined);
    assert(a_source.m_isDefined);

    m_volData.copy(a_src.m_volData, a_srcbox, 0, a_dstcomp, 0,  V_VOLNUMBER);
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_faceData[idir].copy(a_src.m_faceData[idir], a_srcbox, 0, a_dstcomp, 0,  F_FACENUMBER);
    }
  }
/************************/
  void 
  EBDataImplem::
  define(const EBGraph&           a_graph,
         const Box&               a_region)
  {
    m_graph = a_graph;
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
         const Vector<IrregNode>& a_irregGraph,
         const Box&               a_validBox,
         const Box&               a_region)

  {
    BL_PROFILE("EBDataImpem::define");
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

          VolData& vol = m_volData(vof, 0);
          vol.m_volFrac       = node.m_volFrac;
          vol.m_volCentroid   = node.m_volCentroid;
          vol.m_averageFace.m_volIndex      = vof;
          vol.m_averageFace.m_bndryCentroid = node.m_bndryCentroid;
          vol.m_averageFace.m_bndryPhase    = -1;  // EB

          for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
          {
            for (SideIterator sit; sit.ok(); ++sit)
            {
              Vector<FaceIndex> faces = a_graph.getFaces(vof, faceDir, sit());
              int nodeind = node.index(faceDir, sit());
              Vector<Real> areaFracs         = node.m_areaFrac[nodeind];
              Vector<RealVect> faceCentroids = node.m_faceCentroid[nodeind];
              for (int iface = 0; iface < faces.size(); iface++)
              {
                const Real&     areaFracNode     = areaFracs[iface];
                const RealVect& faceCentroidNode = faceCentroids[iface];
                const FaceIndex& face = faces[iface];

                m_faceData[faceDir](face,0).m_areaFrac     = areaFracNode;
                m_faceData[faceDir](face,0).m_faceCentroid = faceCentroidNode;
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
    return m_volData(a_vof, 0).m_volFrac;
  }

/*******************************/
  const Real& EBDataImplem::bndryArea(const VolIndex& a_vof, int face) const
  {
    const VolData& v =  m_volData(a_vof, 0);
    if (v.m_phaseFaces.size()>0)
      return v.m_phaseFaces[face].m_bndryArea;
    assert(face == 0);
    return v.m_averageFace.m_bndryArea;
  }

  const Real& EBDataImplem::bndryArea(const VolIndex& a_vof) const
  {
    static Real zero = 0;
    if (m_volData.getIVS().contains(a_vof.gridIndex()))
      return m_volData(a_vof, 0).m_averageFace.m_bndryArea;

    return zero;
  }

/*******************************/
  const RealVect& EBDataImplem::normal(const VolIndex& a_vof) const
  {
    return m_volData(a_vof, 0).m_averageFace.m_normal;

  }
  const RealVect& EBDataImplem::normal(const VolIndex& a_vof, int face) const
  {
    const VolData& v =  m_volData(a_vof, 0);
    if (v.m_phaseFaces.size()>0)
      return v.m_phaseFaces[face].m_normal;
    assert(face == 0);
    return v.m_averageFace.m_normal;
  }

/*******************************/
  const RealVect& EBDataImplem::centroid(const VolIndex& a_vof) const
  {
    return m_volData(a_vof, 0).m_volCentroid;
  }
/*******************************/
  const RealVect& EBDataImplem::bndryCentroid(const VolIndex& a_vof) const
  {
    return m_volData(a_vof, 0).m_averageFace.m_bndryCentroid;
  }
  const RealVect& EBDataImplem::bndryCentroid(const VolIndex& a_vof, int face) const
  {
    const VolData& v =  m_volData(a_vof, 0);
    if (v.m_phaseFaces.size()>0)
      return v.m_phaseFaces[face].m_bndryCentroid;
    assert(face == 0);
    return v.m_averageFace.m_bndryCentroid;

  }

  /// used by multi-fluid code
  int EBDataImplem::facePhase(const VolIndex& a_vof, int face) const
  {
    return m_volData(a_vof, 0).m_phaseFaces[face].m_bndryPhase;
  }

  /// used by multi-fluid code
  const VolIndex& EBDataImplem::faceIndex(const VolIndex& a_vof, int face) const
  {
    return m_volData(a_vof, 0).m_phaseFaces[face].m_volIndex;
  }

  /// used by multi-fluid code
  void EBDataImplem::setFacePhase(const VolIndex& a_vof, int face, int phase)
  {
    VolData& voldat = m_volData(a_vof, 0);
    voldat.m_phaseFaces[face].m_bndryPhase=phase;
  }

  /// used by multi-fluid code
  void EBDataImplem::setFaceIndex(const VolIndex& a_vof, int face, const VolIndex& index)
  {
    m_volData(a_vof, 0).m_phaseFaces[face].m_volIndex = index;
  }

///
  int EBDataImplem::numFacePhase(const VolIndex& a_vof) const
  {
    return  m_volData(a_vof, 0).m_phaseFaces.size();
  }

/*******************************/
  const RealVect& EBDataImplem::centroid(const FaceIndex& a_face) const
  {
    int faceDir = a_face.direction();
    return m_faceData[faceDir](a_face, 0).m_faceCentroid;
  }
/*******************************/
  const Real& EBDataImplem::areaFrac(const FaceIndex& a_face) const
  {
    int faceDir = a_face.direction();
    return m_faceData[faceDir](a_face, 0).m_areaFrac;
  }
/*******************************/
  void EBData::
  computeNormalsAndBoundaryAreas(const EBGraph& a_graph,
                                 const Box&     a_validRegion)
  {
    EBISBox ebisBox;
    ebisBox.define(a_graph, *this);
    BaseIVFAB<VolData>& volData = m_implem->getVolData();
    if (a_graph.hasIrregular())
    {
      IntVectSet ivsIrreg = a_graph.getIrregCells(a_validRegion);
      for (VoFIterator vofit(ivsIrreg, a_graph); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real bndryArea  =  PolyGeom::bndryArea(vof, ebisBox);
        RealVect normal =  PolyGeom::normal(   vof, ebisBox, bndryArea);

        //copy results to volData
        volData(vof,0).m_averageFace.m_bndryArea = bndryArea;
        volData(vof,0).m_averageFace.m_normal    = normal;
      }
    }

  }
/*******************************/
  void
  EBDataImplem::
  coarsenBoundaryAreaAndNormal(Real&                    a_bndryAreaCoar,
                               RealVect&                a_normalCoar,
                               const Vector<Real>&      a_bndryAreaFine,
                               const Vector<RealVect>&  a_normalFine)
  {
    BL_PROFILE("EBDataImplem::coarsenBoundaryAreaAndNormal");

    Real faceCoarsenFactor = D_TERM(1.0, * 0.5, *0.5);
    // Real voluCoarsenFactor = D_TERM(0.5, * 0.5, *0.5);

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
//  for (int ifine = 0; ifine < a_normalFine.size(); ifine++)
//    {
//      a_bndryAreaCoar += a_bndryAreaFine[ifine];
//    }
//  a_bndryAreaCoar *= faceCoarsenFactor;


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
      std::list<BoundaryData>  boundary;

      for (VoFIterator vofit(ivsIrreg, a_coarGraph); vofit.ok(); ++vofit)
      {
        BL_PROFILE("EBDataImplem::coarsenVoFs_VoFIterator");
        const VolIndex& vofCoar = vofit();
        Vector<VolIndex> vofsFine = a_coarGraph.refine(vofCoar);
        int nFine = vofsFine.size();
        Vector<Real> bndryAreaFine(nFine);
        Vector<Real> volFracFine(nFine);
        Vector<int>  phase(nFine);
        Vector<RealVect> bndryCentroidFine(nFine);
        Vector<RealVect> volCentroidFine(nFine);
        Vector<RealVect> normalFine(nFine);

        for (int ifine = 0; ifine < nFine; ifine++)
        {
          BL_PROFILE("EBDataImplem::coarsenVoFs_fine");
          const VolIndex& vofFine =vofsFine[ifine];

          if (a_fineGraph.isIrregular(vofFine.gridIndex()))
          {
            const VolData& vol = a_fineEBDataImplem.m_volData(vofFine,0);
            bndryAreaFine[ifine] = vol.m_averageFace.m_bndryArea;
            //a_fineEBDataImplem.bndryArea(vofFine);

            volFracFine[ifine] = vol.m_volFrac;
            //  a_fineEBDataImplem.volFrac(vofFine);
            bndryCentroidFine[ifine] = vol.m_averageFace.m_bndryCentroid;
            // a_fineEBDataImplem.bndryCentroid(vofFine);
            volCentroidFine[ifine] = vol.m_volCentroid;
            // a_fineEBDataImplem.centroid(vofFine);
            normalFine[ifine] = vol.m_averageFace.m_normal;
            //  a_fineEBDataImplem.normal(vofFine);
            if (vol.m_phaseFaces.size()>0)
            {
              for (int i=0; i<vol.m_phaseFaces.size(); i++)
              {
                boundary.push_back(vol.m_phaseFaces[i]);
              }
            }
            else
            {
              boundary.push_back(vol.m_averageFace);
            }

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

        VolData& vol = m_volData(vofCoar,0);

        vol.m_volFrac       = volFracCoar;
        vol.m_averageFace.m_bndryArea     = bndryAreaCoar;
        vol.m_averageFace.m_volIndex      = vofCoar;
        vol.m_averageFace.m_bndryCentroid = bndryCentroidCoar;
        vol.m_averageFace.m_normal        = normalCoar;
        vol.m_volCentroid   = volCentroidCoar;

        vol.m_phaseFaces.resize(boundary.size());
        for (int i=0; i<vol.m_phaseFaces.size(); ++i)
        {
          vol.m_phaseFaces[i] = boundary.front();
          boundary.pop_front();
          //vol.m_bndryPhase = phase[0];  // two-phase flow assumption
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
          Vector<FaceIndex> facesFine = a_coarGraph.refine(faceCoar, a_fineGraph);

          Vector<Real>     areaFracsFine(facesFine.size());
          Vector<RealVect> centroidsFine(facesFine.size());
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

          m_faceData[faceDir](faceCoar, 0).m_areaFrac     = areaFracCoar;
          m_faceData[faceDir](faceCoar, 0).m_faceCentroid = centroidCoar;
        } //end loop over faces
      } //end loop over face directions
    }
  }
/*******************************/
  void EBDataImplem::
  coarsenFaceCentroid(RealVect&                a_centroidCoar,
                      const Vector<RealVect>&  a_centroidsFine,
                      const Vector<Real>&      a_areaFracFine,
                      const Vector<FaceIndex>& a_facesFine,
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
                  const Vector<Real>& a_areaFracFine)
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
                            const Vector<Real>&     a_volFracFine,
                            const Vector<RealVect>& a_volCentroidFine,
                            const Vector<VolIndex>& a_fineVoFs,
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
                       const Vector<RealVect>& a_bndryCentroidFine,
                       const Vector<Real>&     a_bndryAreaFine,
                       const Vector<VolIndex>& a_fineVoFs,
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
/*******************************/
  std::size_t 
  EBDataImplem::
  nBytes (const Box& bx, int start_comp, int ncomps) const
  {
    amrex::Error("not implemented");
  }
/*******************************/

  std::size_t copyToMem (const Box& srcbox,
                         int        srccomp,
                         int        numcomp,
                         void*      dst) const
  {
    amrex::Error("not implemented");
  }

/*******************************/

  std::size_t copyFromMem (const Box&  dstbox,
                           int         dstcomp,
                           int         numcomp,
                           const void* src)
  {
    amrex::Error("not implemented");
  }

/*******************************/
