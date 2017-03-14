
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

#ifndef _EBDATA_H_
#define _EBDATA_H_

#include "AMReX_REAL.H"
#include "AMReX_RealVect.H"
#include "AMReX_EBGraph.H"
#include "AMReX_IrregNode.H"
#include "AMReX_BaseIVFAB.H"
#include "AMReX_BaseIFFAB.H"

namespace amrex
{
  class EBIndexSpace;

///
/**
   This class contains all the geometric information
   for an ebisbox.
*/
  class EBDataImplem
  {

  public:

    ///
    EBDataImplem();

    ///
    ~EBDataImplem();

    ///
    /**
       Noop.
    */
    void define(const Box& box, int comps);

    ///
    /**
       Noop
    */
    EBDataImplem(const Box& a_box, int a_comps);

    ///
    EBDataImplem& copy(const EBData&   a_source,
                       const Box&      a_regionFrom,
                       int             a_srccomp,
                       const Box&      a_regionto,
                       int             a_destcomp,
                       int             a_numcomp);



    EBDataImplem& plus(const EBData    & a_src,
                       const Box&        a_srcbox,
                       const Box&        a_destbox,
                       int               a_srccomp,
                       int               a_destcomp,
                       int               a_numcomp=1);
    ///
    /**
       Copy the information from a_source to the over the intersection
       of the box a_region, the box of the current EBDataImplem and
       the box of a_source.  The Interval arguments are ignored.
       This function is required by LevelData.
    */
    void copy(const Box&      a_regionFrom,
              const Interval& a_Cd,
              const Box&      a_regionto,
              const EBDataImplem&   a_source,
              const Interval& a_Cs);

    /// each data holder is defined over the irregular cells of the graph
    void
    defineVoFData(const EBGraph& a_graph, const Box& a_region);

    /// each data holder is defined over the irregular cells of the graph
    void
    defineFaceData(const EBGraph& a_graph, const Box& a_region);

    ///define the whole thing
    void
    define(const EBGraph&           a_graph,
           const Vector<IrregNode>& a_irregData,
           const Box&               a_validBox);

    ///
    void
    coarsenVoFs(const EBDataImplem&   a_fineEBDataImplem,
                const EBGraph&        a_fineGraph,
                const EBGraph&        a_coarGraph,
                const Box&           a_validRegion);

    ///
    void
    coarsenFaces(const EBDataImplem& a_fineEBDataImplem,
                 const EBGraph&      a_fineGraph,
                 const EBGraph&      a_coarGraph,
                 const Box&          a_validRegion);
    ///
    const Real& volFrac(const VolIndex& a_vof) const;

    ///
    const Real& areaFrac(const FaceIndex& a_face1) const;

    ///
    const RealVect& centroid(const FaceIndex& facein) const;

    ///
    const RealVect& centroid(const VolIndex& a_vof) const;

    ///
    const RealVect& bndryCentroid(const VolIndex& a_vof) const;

    ///
    const Real& bndryArea(const VolIndex& a_vof) const;

    ///
    const RealVect& normal(const VolIndex& a_vof) const;

    std::size_t nBytes (const Box& bx, int start_comp, int ncomps) const;

    std::size_t copyToMem (const Box& srcbox,
                           int        srccomp,
                           int        numcomp,
                           void*      dst) const;


    std::size_t copyFromMem (const Box&  dstbox,
                             int         dstcomp,
                             int         numcomp,
                             const void* src);

    /**
       This stuff required by LevelData in parallel:
    int size(const Box& R, const Interval& comps) const;

    ///
    void linearOut(void* buf, const Box& R, const Interval& comps) const;

    ///
    void linearIn(void* buf, const Box& R, const Interval& comps);
    */

    ///
    void addFullIrregularVoFs(const IntVectSet& a_vofsToChange,
                              const EBGraph&    a_newGhostGraph,
                              const BaseIVFAB<VolData>&     a_newGhostData,
                              const EBGraph&    a_oldGhostGraph);

    ///multifluid angels dancing on the heads of pins.
    void addEmptyIrregularVoFs(const IntVectSet& a_vofsToChange,
                               const EBGraph&    a_newGraph);

  private:

    ///
    BaseIVFAB<Real>   m_volData;

    ///
    BaseIFFAB<Real>   m_faceData[SpaceDim];

    ///
    bool m_isFaceDataDefined;

    ///
    bool m_isVoFDataDefined;

    void operator=(const EBDataImplem& ebiin)
      {;}

    EBDataImplem(const EBDataImplem& ebiin)
      {;}

    void
    coarsenFaceCentroid(RealVect&                a_centroidCoar,
                        const Vector<RealVect>&  a_centroidsFine,
                        const Vector<Real>&      a_areaFracFine,
                        const Vector<FaceIndex>& a_facesFine,
                        const FaceIndex&         a_faceCoar);
    void
    coarsenAreaFrac(Real& a_areaFracCoar,
                    const Vector<Real>& a_areaFracFine);

    void
    coarsenVolFracAndCentroid(Real&                   a_volFracCoar,
                              RealVect&               a_volCentroidCoar,
                              const Vector<Real>&     a_volFracFine,
                              const Vector<RealVect>& a_volCentroidFine,
                              const Vector<VolIndex>& a_fineVoFs,
                              const VolIndex&         a_coarVoF);

    void
    coarsenBoundaryAreaAndNormal(Real&                    a_bndryAreaCoar,
                                 RealVect&                a_normalCoar,
                                 const Vector<Real>&      a_bndryAreaFine,
                                 const Vector<RealVect>&  a_normalFine);

    RealVect
    fineToCoarseTransform(const RealVect& a_finePoint,
                          const IntVect&  a_coarCell,
                          const IntVect&  a_fineCell);

    void
    coarsenBndryCentroid(RealVect&               a_bndryCentroidCoar,
                         const Vector<RealVect>& a_bndryCentroidFine,
                         const Vector<Real>&     a_bndryAreaFine,
                         const Vector<VolIndex>& a_fineVoFs,
                         const VolIndex&         a_coarVoF);

    friend class EBIndexSpace;
    friend class EBISLevel;
  };

///
/**
   Ref-counted version of EBDataImplem.
*/
  class EBData
  {
  public:

    ///
    EBData();

    ///
    ~EBData();

    ///
    /**
       Noop.
    */
    void define(const Box& box, int comps);

    ///
    /**
       Noop
    */
    EBData(const Box& a_box, int a_comps);

    ///
    /**
       Copy the information from a_source to the over the intersection
       of the box a_region, the box of the current EBData and
       the box of a_source.  The Interval arguments are ignored.
       This function is required by LevelData.
    */
    EBData& copy(const EBData&   a_source,
                 const Box&      a_regionFrom,
                 int             a_srccomp,
                 const Box&      a_regionto,
                 int             a_destcomp,
                 int             a_numcomp)
      {
        m_implem->copy(a_source, a_regionFrom, a_srccomp, a_destcomp, a_numcomp);
        return *this;
      }



    EBData& plus(const EBData    & a_src,
                 const Box&        a_srcbox,
                 const Box&        a_destbox,
                 int               a_srccomp,
                 int               a_destcomp,
                 int               a_numcomp=1)
      {
        m_implem->plus(a_src, a_srcbox, a_destbox, a_srccomp, a_destcomp, a_numcomp);
        return *this;
      }

    /// each data holder is defined over the irregular cells of the graph
    void
    defineVoFData(const EBGraph& a_graph, const Box& a_region);

    /// each data holder is defined over the irregular cells of the graph
    void
    defineFaceData(const EBGraph& a_graph, const Box& a_region);

    ///define the whole thing
    void
    define(const EBGraph&           a_graph,
           const Vector<IrregNode>& a_irregData,
           const Box&               a_validBox);

    ///
    void
    coarsenVoFs(const EBData&   a_fineEBData,
                const EBGraph&  a_fineGraph,
                const EBGraph&  a_coarGraph,
                const Box&      a_validRegion);

    void
    coarsenFaces(const EBData&  a_fineEBData,
                 const EBGraph& a_fineGraph,
                 const EBGraph& a_coarGraph,
                 const Box&     a_validRegion);
    ///
    const Real& volFrac(const VolIndex& a_vof) const;

    ///
    const Real& areaFrac(const FaceIndex& a_face1) const;

    ///
    const RealVect& centroid(const FaceIndex& facein) const;

    ///
    const RealVect& centroid(const VolIndex& a_vof) const;

    ///
    const RealVect& bndryCentroid(const VolIndex& a_vof) const;
 
    ///
    const Real& bndryArea(const VolIndex& a_vof) const;
 
    ///
    const RealVect& normal(const VolIndex& a_vof) const;
 
    ///
    EBData(const EBData& a_ebiin);

    ///
    /**
       This is a pointer comparison.
    */
    bool operator==(const EBData& a_ebiin);

    ///
    EBData& operator=(const EBData& a_ebiin);

    std::size_t nBytes (const Box& bx, int start_comp, int ncomps) const
      {
        return m_implem->nBytes(bx, start_comp, ncomps);
      }

    std::size_t copyToMem (const Box& srcbox,
                           int        srccomp,
                           int        numcomp,
                           void*      dst) const
      {
        return m_implem->copyToMem(srcbox, srccomp, numcomp, dst);
      }


    std::size_t copyFromMem (const Box&  dstbox,
                             int         dstcomp,
                             int         numcomp,
                             const void* src)
      {
        return m_implem->copyFromMem(dstbox, dstcomp, numcomp, dst);
      }

    ///
    /**
       This stuff required by LevelData in parallel:
    int size(const Box& R, const Interval& comps) const;

    ///
    void linearOut(void* buf, const Box& R, const Interval& comps) const;

    ///
    void linearIn(void* buf, const Box& R, const Interval& comps);
    */

    ///
    void addFullIrregularVoFs(const IntVectSet& a_vofsToChange,
                              const EBGraph&    a_newGhostGraph,
                              const BaseIVFAB<VolData>& a_grownData,
                              const EBGraph&    a_oldGhostGraph);

    ///multifluid angels dancing on the heads of pins.
    void addEmptyIrregularVoFs(const IntVectSet& a_vofsToChange,
                               const EBGraph&    a_newGraph);

    ///
    void
    computeNormalsAndBoundaryAreas(const EBGraph& a_graph,
                                   const Box&     a_validRegion);

    BaseIVFAB<VolData>& getVolData()
      {
        return m_implem->getVolData();
      }

    const BaseIVFAB<VolData>& getVolData() const
      {
        return m_implem->getVolData();
      }

    ///
    static int preAllocatable()
      {
        return 2; // dyanmic allocatable.
      }
  private:

    ///
    std::shared_ptr<EBDataImplem> m_implem;

    friend class EBIndexSpace;
    friend class EBISLevel;
  };

/*******************************/
  inline const Real& EBData::volFrac(const VolIndex& a_vof) const
  {
    return m_implem->volFrac(a_vof);
  }
/*******************************/
  inline const Real& EBData::bndryArea(const VolIndex& a_vof) const
  {
    return m_implem->bndryArea(a_vof);
  }
/*******************************/
  inline const RealVect& EBData::normal(const VolIndex& a_vof) const
  {
    return m_implem->normal(a_vof);
  }
/*******************************/
  inline const RealVect& EBData::centroid(const VolIndex& a_vof) const
  {
    return m_implem->centroid(a_vof);
  }
/*******************************/
  inline const RealVect& EBData::bndryCentroid(const VolIndex& a_vof) const
  {
    return m_implem->bndryCentroid(a_vof);
  }
/*******************************/
  inline const RealVect& EBData::centroid(const FaceIndex& a_face) const
  {
    return m_implem->centroid(a_face);
  }
/*******************************/
  inline const Real& EBData::areaFrac(const FaceIndex& a_face) const
  {
    return m_implem->areaFrac(a_face);
  }
/*******************************/

}
#endif
