
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
#include "AMReX_FArrayBox.H"

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
    static bool s_verbose;
    ///sure hope this does not matter
    typedef int value_type;

    ///
    EBDataImplem();

    ///
    ~EBDataImplem();

    ///this is only defined if there are irregular cells
    //const FArrayBox& getVolumeMask()
    //  {
    //    return m_volMask;
    //  }

    ///
    /**
       Noop.
    */
    void define(const Box& box, int comps);

    //this is the define to use when you are just going to copy the data into it.  region should include ghost cells
    //this just defines the data holders
    void define(const EBGraph&           a_graph,
                const Box&               a_region);

    ///
    /**
       Noop
    */
    EBDataImplem(const Box& a_box, int a_comps);

    ///
    EBDataImplem& copy(const EBDataImplem&   a_source,
                       const Box&            a_regionFrom,
                       int                   a_srccomp,
                       const Box&            a_regionto,
                       int                   a_destcomp,
                       int                   a_numcomp);



    EBDataImplem& plus(const EBDataImplem & a_src,
                       const Box&           a_srcbox,
                       const Box&           a_destbox,
                       int                  a_srccomp,
                       int                  a_destcomp,
                       int                  a_numcomp=1)
      {
        amrex::Error("this operation makes no sense here");
        return *this;
      }



    ///define the whole thing. valid box= no ghost, region has ghost, defines data holders and fills them
    void
    define(const EBGraph&           a_graph,
           const Vector<IrregNode>& a_irregData,
           const Box&               a_validBox,
           const Box&               a_region);

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
    RealVect centroid(const FaceIndex& facein) const;

    ///
    RealVect centroid(const VolIndex& a_vof) const;

    ///
    RealVect bndryCentroid(const VolIndex& a_vof) const;

    ///
    const Real& bndryArea(const VolIndex& a_vof) const;

    ///
    static bool preAllocatable () { return false; }
    static bool isCopyOMPSafe () { return false; }

    ///
    RealVect normal(const VolIndex& a_vof) const;

    ///this is for copy and fillboundary where the metadata is already filled in
    std::size_t nBytes (const Box& bx, int start_comp, int ncomps) const;

    ///this is for copy and fillboundary where the metadata is already filled in
    std::size_t copyToMem (const Box& srcbox,
                           int        srccomp,
                           int        numcomp,
                           void*      dst) const;


    ///this is for copy and fillboundary where the metadata is already filled in
    std::size_t copyFromMem (const Box&  dstbox,
                             int         dstcomp,
                             int         numcomp,
                             const void* src);



    ///this includes all metadata
    std::size_t nBytesFull() const;

    ///this includes all metadata
    std::size_t copyToMemFull(void*      dst) const;


    ///this includes all metadata
    std::size_t copyFromMemFull(const void* src);


    ///
    void addFullIrregularVoFs(const IntVectSet& a_vofsToChange,
                              const Box       & a_region);

    ///
    void
    computeNormalsAndBoundaryAreas(const EBGraph& a_graph,
                                   const Box&     a_validRegion);

    BaseIVFAB<Real>& getVolData()
      {
        return m_volData;
      }

    const BaseIVFAB<Real>& getVolData() const
      {
        return m_volData;
      }

    BaseIFFAB<Real>& getFaceData(int a_idir)
      {
        return m_faceData[a_idir];
      }

    const BaseIFFAB<Real>& getFaceData(int a_idir) const
      {
        return m_faceData[a_idir];
      }
    const Box& getRegion() const
      {
        return m_region;
      }

    ///invalid but necessary for fabarray to compile in certain cases
    EBDataImplem(const Box& a_region,
                 int        n ,
                 bool       alloc,
                 bool       shared)
      {
        //this strictly for checkpointing 
        //m_hasVolumeMask = false;
        m_region = a_region;
      }
private:


    //internal.  initializes places where we know the answer
    void setCoveredAndRegular();

    //void defineVolumeMask();
    //void setVolumeMask();

    Box m_region;

    EBGraph m_graph;
    ///
    BaseIVFAB<Real>   m_volData;

    ///
    BaseIFFAB<Real>   m_faceData[SpaceDim];

    ///
    bool m_isDefined;

    /// only defined if there are cut cells somewhere in m_region
    //FArrayBox m_volMask;

    ///true if there are cut cells somewhere in m_region
    //bool m_hasVolumeMask;

    void operator=(const EBDataImplem& ebiin) = delete;
    EBDataImplem(const EBDataImplem& ebiin) = delete;

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

    void
    init_snan ();

    void
    init_qnan ();

    friend class EBIndexSpace;
    friend class EBISLevel;
  };

  extern void null_deleter_ebdi(EBDataImplem *);

///
/**
   Ref-counted version of EBDataImplem.
*/
  class EBData
  {
  public:

    ///sure hope this does not matter
    typedef int value_type;

    ///
    EBData() : m_implem( std::make_shared<EBDataImplem>() )
      {
      }

    ///invalid but necessary for fabarray to compile in certain cases
    EBData(const Box& a_region,
           int        n,
           bool       alloc,
           bool       shared)
        : m_implem( std::make_shared<EBDataImplem>(a_region, n, alloc, shared) )

      {
      }
    ///
    ~EBData()
      {
      }


    EBData (EBData&& rhs) noexcept = default;

    ///this is only defined if there are irregular cells
    //const FArrayBox& getVolumeMask()
    //  {
    //    return m_implem->getVolumeMask();
    //  }

    ///
    const Box& getRegion() const
      {
        return m_implem->getRegion();
      }

    ///this is just a way around smart pointers being too smart.
    ///this keeps the input pointer from being deleted.  Not for the faint of heart.
    EBData(EBDataImplem* a_implemPtr)
      {
        m_implem = std::shared_ptr<EBDataImplem>(a_implemPtr, &null_deleter_ebdi);
      }
    ///
    void define(const Box& a_box, int a_comps)
      {
        m_implem->define(a_box, a_comps);
      }

    void resize(const Box& a_box, int a_comps)
      {
        m_implem->define(a_box, a_comps);
      }

    ///
    static bool preAllocatable () { return false; }
    static bool isCopyOMPSafe () { return false; }

    ///
    /**
    */
    EBData(const Box& a_box, int a_comps)
        : m_implem( std::make_shared<EBDataImplem>(a_box, a_comps) )
      {
      }

    ///this includes all metadata
    std::size_t nBytesFull() const
      {
        return m_implem->nBytesFull();
      }


    ///this includes all metadata
    std::size_t copyToMemFull(void*      dst) const
      {
        return m_implem->copyToMemFull(dst);
      }


    ///this includes all metadata
    std::size_t copyFromMemFull(const void* src)
      {
        return m_implem->copyFromMemFull(src);
      }

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
                 const Box&      a_regionTo,
                 int             a_destcomp,
                 int             a_numcomp)
      {
        m_implem->copy(*a_source.m_implem, a_regionFrom, a_srccomp, a_regionTo, a_destcomp, a_numcomp);
        return *this;
      }



    EBData& plus(const EBData    & a_src,
                 const Box&        a_srcbox,
                 const Box&        a_destbox,
                 int               a_srccomp,
                 int               a_destcomp,
                 int               a_numcomp=1)
      {
        m_implem->plus(*a_src.m_implem, a_srcbox, a_destbox, a_srccomp, a_destcomp, a_numcomp);
        return *this;
      }

    ///define the whole thing from primitive objects--region should include ghost cells--
    //but it is ok for the valid box to be smaller
    void
    define(const EBGraph&           a_graph,
           const Vector<IrregNode>& a_irregData,
           const Box&               a_validBox,
           const Box&               a_region)
    
      {
        m_implem->define(a_graph, a_irregData, a_validBox, a_region);
      }

    //this is the define to use when you are just going to copy the data into it.  region should include ghost cells
    void define(const EBGraph&           a_graph,
                const Box&               a_region)
      {
        m_implem->define(a_graph, a_region);
      }

    ///
    void
    coarsenVoFs(const EBData&   a_fineEBData,
                const EBGraph&  a_fineGraph,
                const EBGraph&  a_coarGraph,
                const Box&      a_validRegion)
      {
        m_implem->coarsenVoFs(*a_fineEBData.m_implem, a_fineGraph, a_coarGraph, a_validRegion);
      }

    void
    coarsenFaces(const EBData&  a_fineEBData,
                 const EBGraph& a_fineGraph,
                 const EBGraph& a_coarGraph,
                 const Box&     a_validRegion)
      {
        m_implem->coarsenFaces(*a_fineEBData.m_implem, a_fineGraph, a_coarGraph, a_validRegion);
      }

    ///
    const Real& volFrac(const VolIndex& a_vof) const;

    ///
    const Real& areaFrac(const FaceIndex& a_face1) const;

    ///
    RealVect centroid(const FaceIndex& facein) const;

    ///
    RealVect centroid(const VolIndex& a_vof) const;

    ///
    RealVect bndryCentroid(const VolIndex& a_vof) const;
 
    ///
    const Real& bndryArea(const VolIndex& a_vof) const;
 
    ///
    RealVect normal(const VolIndex& a_vof) const;
 
    ///
    EBData(const EBData& a_ebiin)
      {
        m_implem = a_ebiin.m_implem;
      }

    ///
    /**
       This is a pointer comparison.
    */
    bool operator==(const EBData& a_ebiin)
      {
          return m_implem.get() == a_ebiin.m_implem.get();
      }

    ///
    EBData& operator=(const EBData& a_ebiin) = default;

    ///this is for copy or fillboundary where the meta data is already filled
    std::size_t nBytes (const Box& bx, int start_comp, int ncomps) const
      {
        return m_implem->nBytes(bx, start_comp, ncomps);
      }

    ///this is for copy or fillboundary where the meta data is already filled
    std::size_t copyToMem (const Box& srcbox,
                           int        srccomp,
                           int        numcomp,
                           void*      dst) const
      {
        return m_implem->copyToMem(srcbox, srccomp, numcomp, dst);
      }


    ///this is for copy or fillboundary where the meta data is already filled
    std::size_t copyFromMem (const Box&  dstbox,
                             int         dstcomp,
                             int         numcomp,
                             const void* src)
      {
        return m_implem->copyFromMem(dstbox, dstcomp, numcomp, src);
      }




    ///
    void addFullIrregularVoFs(const IntVectSet&      a_vofsToChange,
                              const Box       &      a_region)
      {
        m_implem->addFullIrregularVoFs(a_vofsToChange, a_region);
      }


    ///
    void
    computeNormalsAndBoundaryAreas(const EBGraph& a_graph,
                                   const Box&     a_validRegion)
      {
        return m_implem->computeNormalsAndBoundaryAreas(a_graph, a_validRegion);
      }

    BaseIVFAB<Real>& getVolData()
      {
        return m_implem->getVolData();
      }

    const BaseIVFAB<Real>& getVolData() const
      {
        return m_implem->getVolData();
      }


    BaseIFFAB<Real>& getFaceData(int a_idir)
      {
        return m_implem->getFaceData(a_idir);
      }

    const BaseIFFAB<Real>& getFaceData(int a_idir) const
      {
        return m_implem->getFaceData(a_idir);
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
  inline RealVect EBData::normal(const VolIndex& a_vof) const
  {
    return m_implem->normal(a_vof);
  }
/*******************************/
  inline RealVect EBData::centroid(const VolIndex& a_vof) const
  {
    return m_implem->centroid(a_vof);
  }
/*******************************/
  inline RealVect EBData::bndryCentroid(const VolIndex& a_vof) const
  {
    return m_implem->bndryCentroid(a_vof);
  }
/*******************************/
  inline  RealVect EBData::centroid(const FaceIndex& a_face) const
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
