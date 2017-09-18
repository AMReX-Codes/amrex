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
#include "AMReX_parstream.H"
#include "AMReX_Print.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_MemPool.H"

#include <limits>

namespace amrex
{

  bool EBDataImplem::s_verbose = false;

  static const IntVect   ebd_debiv(D_DECL(15, 6, 0));
  static const VolIndex  ebd_debvof(ebd_debiv, 0);
  static const IntVect   ebd_debivlo(D_DECL(190,15,0));
  static const IntVect   ebd_debivhi(D_DECL(191,15,0));
  static const VolIndex  ebd_debvoflo(ebd_debivlo, 0);
  static const VolIndex  ebd_debvofhi(ebd_debivhi, 0);
  static const FaceIndex ebd_debface(ebd_debvoflo, ebd_debvofhi);

  /******/
  void checkFaceData (const BaseIFFAB<Real> a_faceData[SpaceDim], const Box& a_valid, const string& a_identifier)
  {
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      const BaseIFFAB<Real> & data  = a_faceData[idir];
      const Array<FaceIndex>& faces = data.getFaces();
      for(int iface = 0; iface < faces.size(); iface++)
      {
        if(faces[iface] == ebd_debface)
        {
          amrex::Print() << a_identifier << ", valid = " << a_valid << ", areaFrac(" << ebd_debface << ")=" << data(ebd_debface,0) << endl;
        }
      }
    }
  }
  /******/
  extern void null_deleter_ebdi (EBDataImplem* a_input)
  {
  }
  /******/
  void
  EBDataImplem::
  setCoveredAndRegular ()
  {
    //regular cell is always unity
    const vector<VolIndex>& vofs = m_volData.getVoFs();
    for(int ivof= 0; ivof < vofs.size(); ivof++)
    {
      const IntVect& iv = vofs[ivof].gridIndex();
      bool regular = m_graph.isRegular(iv);
      bool covered = m_graph.isCovered(iv);
      if(regular || covered)
      {
        if(regular)
        {
          m_volData(vofs[ivof], V_VOLFRAC) = 1;
        }
        else if(covered)
        {
          m_volData(vofs[ivof], V_VOLFRAC) = 0;
        }
        m_volData(vofs[ivof], V_BNDAREA) = 0;
        for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_volData(vofs[ivof], V_VOLCENTROIDX + idir) = 0.;
          m_volData(vofs[ivof], V_BNDCENTROIDX + idir) = 0.;
          m_volData(vofs[ivof], V_NORMALX      + idir) = 0.;
        }
        //no real boundary here but putting something valid for the normal
        m_volData(vofs[ivof], V_NORMALX) = 1.;
      }
    }
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      BaseIFFAB<Real>& faceData = m_faceData[idir];
      const vector<FaceIndex>& faces = faceData.getFaces();
      for(int iface = 0; iface < faces.size(); iface++)
      {
        const FaceIndex& face = faces[iface];
        const IntVect  & ivlo = face.gridIndex(Side::Lo);
        const IntVect  & ivhi = face.gridIndex(Side::Hi);
        bool regular, covered;
        if(!face.isBoundary())
        {
          regular = (m_graph.isRegular(ivlo) ||  m_graph.isRegular(ivhi));
          covered = (m_graph.isCovered(ivlo) ||  m_graph.isCovered(ivhi));
        }
        else
        {
          int cellLo = face.cellIndex(Side::Lo);
//          int cellHi = face.cellIndex(Side::Hi);
          IntVect ivCheck;
          if(cellLo >=0)
          {
            ivCheck = ivlo;
          }
          else
          {
            
            ivCheck = ivhi;
          }
          regular = m_graph.isRegular(ivCheck);
          covered = m_graph.isCovered(ivCheck);
        }
        if(regular || covered)
        {
          if(regular)
          {
            faceData(face, F_AREAFRAC) = 1.;
          }
          else if(covered)
          {
            faceData(face, F_AREAFRAC) = 0.;
          }
          for (int idir = 0; idir < SpaceDim; idir++)
          {
            faceData(face, F_FACECENTROIDX + idir) = 0.;
          }
        }
      }
    }
  }
 /************************/
  void 
  EBDataImplem::
  addFullIrregularVoFs (const IntVectSet     &    a_vofsToChange,
                        const Box            &    a_region)
  {
    IntVectSet subset = a_vofsToChange;
    subset &= a_region;
    if (!subset.isEmpty())
    {
      //calculate set by adding in new intvects
      const IntVectSet& ivsOld = m_volData.getIVS();
      const IntVectSet& ivsIrreg = m_graph.getIrregCells(m_region);
      IntVectSet ivscomp1 = ivsIrreg;
      ivscomp1 -= ivsOld;
      IntVectSet ivscomp2 = ivsOld;
      ivscomp2 -= ivsIrreg;
      IntVectSet ivsNew = ivsOld;
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
      m_volData.define(ivsNew, m_graph, V_VOLNUMBER);
      m_volData.copy(oldVolData, minBox, 0, minBox, 0, nvolcomp);
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        m_faceData[idir].define(ivsNew, m_graph, idir, F_FACENUMBER);
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
            const Array<FaceIndex>& allFaces = m_graph.getAllFaces(ivhi, idir, Side::Lo);
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
          //idir = face direction, jdir= centroid direction
          for (int jdir = 0; jdir < SpaceDim; jdir++)
          {
            m_faceData[idir](faceit(), F_FACECENTROIDX+jdir) = faceCentroid[idir];
          }
        }
      }
    }
  }
/************************/
/************************/
  EBDataImplem::
  EBDataImplem ()
  {
    m_isDefined = false;
    m_hasVolumeMask = false;

  }
/************************/
  EBDataImplem::
  ~EBDataImplem ()
  {
  }
/************************/
  void 
  EBDataImplem::
  define (const Box& box, int comps)
  {
  }
/************************/
  EBDataImplem::
  EBDataImplem (const Box& a_box, int a_comps)
  {
    m_isDefined = false;
    m_hasVolumeMask = false;
  }
/************************/
  EBDataImplem& 
  EBDataImplem::
  copy (const EBDataImplem&   a_src,
        const Box&            a_srcbox,
        int                   a_srccomp,
        const Box&            a_dstbox,
        int                   a_dstcomp,
        int                   a_numcomp)
  {
    assert(m_isDefined);
    assert(a_src.m_isDefined);
    //if(a_src.m_volData.getIVS().contains(ebd_debiv))
    //{
    //  pout() << "ebdata copy src = " << a_src.m_graph.getDomain() << ",region "  << a_src.m_region << ", a_region = " << a_dstbox <<  ", data(" << ebd_debvof << ",1) = " << a_src.m_volData(ebd_debvof, 1) << endl;
    //}
    //if(a_src.m_faceData[1].hasFace(ebd_debface))
    //{
    //  pout() << "ebdata copy src = " << a_src.m_graph.getDomain() << ",region "  << a_src.m_region << "\t, srcbox = " << a_srcbox <<  ", data(" << ebd_debface << ",0) = " << a_src.m_faceData[1](ebd_debface, 0) << endl;
    //}

    BL_PROFILE_VAR("EBData_copy_voldata", copy_voldata);
    m_volData.copy(a_src.m_volData, a_srcbox, 0, a_dstbox, 0,  V_VOLNUMBER);
    BL_PROFILE_VAR_STOP(copy_voldata);




    BL_PROFILE_VAR("EBData_copy_facedata", copy_facedata);
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box grownBox = a_srcbox;
      grownBox.grow(1);
      grownBox &= m_graph.getDomain();
      m_faceData[idir].copy(a_src.m_faceData[idir], grownBox, 0, grownBox, 0,  F_FACENUMBER);

    }
    BL_PROFILE_VAR_STOP(copy_facedata);
    
    if(m_hasVolumeMask && a_src.m_hasVolumeMask)
    {
      m_volMask.copy(a_src.m_volMask, a_srcbox, 0, a_dstbox, 0, m_volMask.nComp());
    }
    return *this;
  }
/************************/
  void 
  EBDataImplem::
  define (const EBGraph&           a_graph,
          const Box&               a_region)
  {
    m_graph = a_graph;
    m_region = a_region & a_graph.getDomain();
    const IntVectSet& ivsIrreg = a_graph.getIrregCells(a_region);
    m_isDefined = true;
    m_volData.define(ivsIrreg, a_graph, V_VOLNUMBER);
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      //this directional grow is to accomodate the fact that there 
      //might be an irregular cell lurking on the other side of region
      //and it will have faces that EBData will have to own.
      Box regionG1D = a_region;
      regionG1D.grow(idir, 1);
      regionG1D &= a_graph.getDomain();
      const IntVectSet& ivsIrregG1D = a_graph.getIrregCells(regionG1D);
      m_faceData[idir].define(ivsIrregG1D, a_graph, idir, F_FACENUMBER);
    }

//#if !defined(NDEBUG) || defined(BL_TESTING)
    //   init_snan();  // Currently we rely on this to indicate bad data.
    // wz. How could std::isnan fail on snan?  This is gcc 4.8.4.
    init_qnan(); // So we have to initialize it to qnan;
//#endif
    setCoveredAndRegular();
    defineVolumeMask();
  }
/************************/
  void
  EBDataImplem::
  define (const EBGraph&           a_graph,
          const Array<IrregNode>&  a_irregGraph,
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
              const Array<FaceIndex>& faces = a_graph.getFaces(vof, faceDir, sit());
              int nodeind = node.index(faceDir, sit());
              const Array<Real>& areaFracs         = node.m_areaFrac[nodeind];
              const Array<RealVect>& faceCentroids = node.m_faceCentroid[nodeind];
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

//begin debug
//    checkFaceData(m_faceData, a_validBox, string("right after irregnode init"));
//end debug
    computeNormalsAndBoundaryAreas(a_graph, a_validBox);
    //if(m_volData.getIVS().contains(ebd_debiv))
    //{
    //  pout() << "ebdata initial define::domain = " << m_graph.getDomain() << ",region "  << m_region << ", a_region = " << a_region <<  ", data(" << ebd_debvof << ",1) = " << m_volData(ebd_debvof, 1) << endl;
    //}
    //if(m_faceData[1].hasFace(ebd_debface))
    //{
    //  pout() << "ebdata initial define domain  = " << m_graph.getDomain() << ",region "  << m_region << ", a_region = " << a_region <<  ", data(" << ebd_debface << ",0) = " << m_faceData[1](ebd_debface, 0) << endl;
    //}

    defineVolumeMask();
    setVolumeMask();
  }


  void EBDataImplem::init_snan ()
  {
#ifdef BL_USE_DOUBLE
      std::size_t n = m_volData.size();
      if (n > 0) {
          amrex_array_init_snan(m_volData.dataPtr(), n);
      }
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
          n = m_faceData[idim].size();
          if (n > 0) {
              amrex_array_init_snan(m_faceData[idim].dataPtr(), n);
          }
      }
#else
      init_qnan();
#endif
  }

  void EBDataImplem::init_qnan ()
  {
      std::size_t n = m_volData.size();
      if (n > 0) {
          std::fill(m_volData.dataPtr(), m_volData.dataPtr() + n,
                    std::numeric_limits<Real>::quiet_NaN());
      }
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
          n = m_faceData[idim].size();
          if (n > 0) {
              std::fill(m_faceData[idim].dataPtr(), m_faceData[idim].dataPtr() + n,
                        std::numeric_limits<Real>::quiet_NaN());
          }
      }
  }

/*******************************/
  const Real& EBDataImplem::volFrac (const VolIndex& a_vof) const
  {
    return m_volData(a_vof, V_VOLFRAC);
  }

/*******************************/
  const Real& EBDataImplem::bndryArea (const VolIndex& a_vof) const
  {
    return m_volData(a_vof, V_BNDAREA);
  }

/*******************************/
  RealVect EBDataImplem::normal (const VolIndex& a_vof) const
  {
    RealVect retval;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      retval[idir] = m_volData(a_vof, V_NORMALX+idir);
    }
    return retval;

  }
/*******************************/
  RealVect EBDataImplem::centroid (const VolIndex& a_vof) const
  {
    RealVect retval;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      retval[idir] = m_volData(a_vof, V_VOLCENTROIDX+idir);
    }
    return retval;
  }
/*******************************/
  RealVect EBDataImplem::bndryCentroid (const VolIndex& a_vof) const
  {
    RealVect retval;

    for(int idir = 0; idir < SpaceDim; idir++)
    {
      retval[idir] = m_volData(a_vof, V_BNDCENTROIDX+idir);
    }
    return retval;
  }
/*******************************/
  RealVect EBDataImplem::centroid (const FaceIndex& a_face) const
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
  computeNormalsAndBoundaryAreas (const EBGraph& a_graph,
                                  const Box&     a_region)
  {
    EBISBox ebisBox;
    EBData tempData(this);
    ebisBox.define(a_graph, tempData);

    if (a_graph.hasIrregular())
    {
      const IntVectSet& ivsIrreg = a_graph.getIrregCells(a_region);
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
  coarsenBoundaryAreaAndNormal (Real&                    a_bndryAreaCoar,
                                RealVect&                a_normalCoar,
                                const Array<Real>&      a_bndryAreaFine,
                                const Array<RealVect>&  a_normalFine)
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
  coarsenVoFs (const EBDataImplem&  a_fineEBDataImplem,
               const EBGraph&       a_fineGraph,
               const EBGraph&       a_coarGraph,
               const Box&           a_validRegion)
  {
    BL_PROFILE("EBDataImplem::coarsenVoFs");
    //unlike before, define needs to be called first
    assert(m_isDefined);

    if (a_coarGraph.hasIrregular())
    {
      const IntVectSet& ivsIrreg = a_coarGraph.getIrregCells(a_validRegion);

      for (VoFIterator vofit(ivsIrreg, a_coarGraph); vofit.ok(); ++vofit)
      {
        BL_PROFILE("EBDataImplem::coarsenVoFs_VoFIterator");
        const VolIndex& vofCoar = vofit();
        const Array<VolIndex>& vofsFine = a_coarGraph.refine(vofCoar);
        int nFine = vofsFine.size();
        Array<Real> bndryAreaFine(nFine);
        Array<Real> volFracFine(nFine);
        Array<int>  phase(nFine);
        Array<RealVect> bndryCentroidFine(nFine);
        Array<RealVect> volCentroidFine(nFine);
        Array<RealVect> normalFine(nFine);

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
  coarsenFaces (const EBDataImplem& a_fineEBDataImplem,
                const EBGraph&      a_fineGraph,
                const EBGraph&      a_coarGraph,
                const Box&          a_validRegion)
  {
    BL_PROFILE("EBDataImplem::coarsenFaces");

    //unlike before, the define function has to be called first.
    assert(m_isDefined);

    const IntVectSet& ivsIrreg = a_coarGraph.getIrregCells(a_validRegion);
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
          const Array<FaceIndex>& facesFine = a_coarGraph.refine(faceCoar, a_fineGraph);

          Array<Real>     areaFracsFine(facesFine.size());
          Array<RealVect> centroidsFine(facesFine.size());
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
  coarsenFaceCentroid (RealVect&                a_centroidCoar,
                       const Array<RealVect>&  a_centroidsFine,
                       const Array<Real>&      a_areaFracFine,
                       const Array<FaceIndex>& a_facesFine,
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
  coarsenAreaFrac (Real& a_areaFracCoar,
                   const Array<Real>& a_areaFracFine)
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
  coarsenVolFracAndCentroid (Real&                   a_volFracCoar,
                             RealVect&               a_volCentroidCoar,
                             const Array<Real>&     a_volFracFine,
                             const Array<RealVect>& a_volCentroidFine,
                             const Array<VolIndex>& a_fineVoFs,
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
  fineToCoarseTransform (const RealVect& a_finePoint,
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
  coarsenBndryCentroid (RealVect&               a_bndryCentroidCoar,
                        const Array<RealVect>& a_bndryCentroidFine,
                        const Array<Real>&     a_bndryAreaFine,
                        const Array<VolIndex>& a_fineVoFs,
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
    //do we have a volume mask that intersects this box?
    retval += sizeof(int);
    if(m_graph.hasIrregular(bx))
    {
      retval += m_volMask.nBytes(bx, 0, V_VOLNUMBER);
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

    //if(m_volData.getIVS().contains(ebd_debiv))
    //{
    //  pout() << "ebdata copyto  mem::domain = " << m_graph.getDomain() << ",region "  << m_region << ", a_bx = " << bx <<  ", data(" << ebd_debvof << ",1) = " << m_volData(ebd_debvof, 1) << endl;
    //}
    //if(m_faceData[1].hasFace(ebd_debface))
    //{
    //  pout() << "ebdata copyto mem::domain  = " << m_graph.getDomain() << ",region "  << m_region << ", a_bx = " << bx <<  ", data(" << ebd_debface << ",0) = " << m_faceData[1](ebd_debface, 0) << endl;
    //}

    size_t incrval = m_volData.copyToMem(bx, 0, V_VOLNUMBER, buf);
    retval += incrval;
    buf    += incrval;

    for(int idir = 0; idir < SpaceDim; idir++)
    {

      incrval = m_faceData[idir].copyToMem(bx, 0, F_FACENUMBER, buf);

      retval += incrval;
      buf    += incrval;

    }

    //do we have a volume mask that intersects this box?
    int hasIrreg = 0;
    if(m_graph.hasIrregular(bx))
    {
      hasIrreg = 1;
    }
    int* intbuf = (int*) buf;
    *intbuf = hasIrreg;

    buf    += sizeof(int);
    retval += sizeof(int);
    if(hasIrreg == 1)
    {
      retval += m_volMask.copyToMem(bx, 0, V_VOLNUMBER, buf);
    }
    //pout() << "ebdata:: copytomem:   bx = " << bx  << ", retval = " << retval << endl;
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

    //do we have a volume mask that intersects this box?
    int hasIrreg = *buf;
    buf    += sizeof(int);
    retval += sizeof(int);
    if(hasIrreg == 1)
    {
      retval += m_volMask.copyFromMem(bx, 0, V_VOLNUMBER, buf);
    }
    //if(m_volData.getIVS().contains(ebd_debiv))
    //{
    //  pout() << "ebdata copyfrommem::domain = " << m_graph.getDomain() << ",region "  << m_region << ", a_bx = " << bx <<  ", data(" << ebd_debvof << ",1) = " << m_volData(ebd_debvof, 1) << endl;
    //}
    //pout() << "ebdata:: copyfrommem: bx = " << bx  << ", retval = " << retval << endl;
    //if(m_faceData[1].hasFace(ebd_debface))
    //{
    //  pout() << "ebdata copyfrommem::domain  = " << m_graph.getDomain() << ",region "  << m_region << ", a_bx = " << bx <<  ", data(" << ebd_debface << ",0) = " << m_faceData[1](ebd_debface, 0) << endl;
    //}
    return retval;
  }
/*******************************/
  std::size_t 
  EBDataImplem::
  nBytesFull () const
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
  copyToMemFull (void* dst) const
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
  copyFromMemFull (const void* src)
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

    defineVolumeMask();
    setVolumeMask();

    m_isDefined = true;
    return retval;
  }

/*******************************/
  void 
  EBDataImplem::
  defineVolumeMask ()
  {
    if(m_graph.hasIrregular(m_region))
    {
      m_hasVolumeMask = true;
      m_volMask.resize(m_region, V_VOLNUMBER);
    }
  }
/*******************************/
  void 
  EBDataImplem::
  setVolumeMask ()
  {
    if(m_hasVolumeMask)
    {
      //set everything to zero. then loop and set regular volume fractions to 1.
      m_volMask.setVal(0.0);
      for(BoxIterator bit(m_region); bit.ok(); ++bit)
      {
        if(m_graph.isRegular(bit()))
        {
          m_volMask(bit(), V_VOLFRAC) = 1.0;
        }
      }

      //  now loop through the volumes  and set stuff (if it is not nan)
      const Array<VolIndex>& vofs = m_volData.getVoFs();
      for(int ivof = 0; ivof < vofs.size(); ivof++)
      {
        const VolIndex& vof = vofs[ivof];
        for(int icomp = 0; icomp < V_VOLNUMBER; icomp++)
        {
          if(!std::isnan(m_volData(vof, icomp)))
          {
            m_volMask(vof.gridIndex(), icomp) = m_volData(vof, icomp);
          }
        }
      }
    }
  }
/*******************************/
}
