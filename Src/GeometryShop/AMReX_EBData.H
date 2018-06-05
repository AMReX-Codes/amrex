#ifndef AMREX_EBDATA_H_
#define AMREX_EBDATA_H_

#include "AMReX_REAL.H"
#include "AMReX_RealVect.H"
#include "AMReX_EBGraph.H"
#include "AMReX_IrregNode.H"
#include "AMReX_BaseIVFAB.H"
#include "AMReX_BaseIFFAB.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_EB_TYPEDEFS.H"
#include "AMReX_IndexedMoments.H"

namespace amrex
{
  class EBIndexSpace;

  class EBIndexSpace;
  class BaseIF;

  template <class T>
  int linearSize(const T& inputT);

  template <class T>
  void linearIn(T& a_outputT, const void* const inBuf);

  template <class T>
  void linearOut(void* const a_outBuf, const T& inputT);

  ///
  /**
     Class that encapsulates the geometric information at an irregular face
  */
  class BoundaryData
  {
  public:
    ///
    BoundaryData()
      {
      }

    /// the moment at  the irregular face associated with the  monomial with the input exponents 
    /**
       Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
       returns integral_over_irregular_area((x^p y^q z^r) dA)
    **/
    IndMomSpaceDim m_EBMoments;


    /// the derivatives of the normal
    /**
       (d_p)(n_d)
    **/
    IndMomSpaceDim m_normalPartialDeriv[SpaceDim];

    /// the derivatives of the normal
    /**
       (d_p)(n_d)
    **/
    IndMomSpaceDim m_EBNormalMoments[SpaceDim];

    BoundaryData& operator=(const BoundaryData& a_in)
      {
        if(&a_in != this)
        {
          m_EBMoments = a_in.m_EBMoments;
          for(int idir = 0; idir < SpaceDim; idir++)
          {
            m_normalPartialDeriv[idir] = a_in.m_normalPartialDeriv[idir];
            m_EBNormalMoments   [idir] = a_in.m_EBNormalMoments   [idir];
          }
        }
        return *this;
      }
    static size_t linearSize();
  
    void linearIn( const char* const a_inbuf);
    void linearOut(const char* const a_inbuf) const;
  };

///
/**
 */
  class VolData
  {
  public:
    IndMomSpaceDim                    m_volumeMoments;
    BoundaryData                      m_averageFace;
  };

/// VolData specializations of linear* functions used by LevelData
//[NOTE: see also SPMD.H. <dbs>]
  template < >
  int linearSize(const VolData& vdata);
  template < >
  void linearIn(VolData& a_outputT, const void* const inBuf);
  template < >
  void linearOut(void* const a_outBuf, const VolData& a_inputT);



///
/**
 */
  class FaceData
  {
  public:
    IndMomSDMinOne  m_faceMoments;
  };

  template < >
  int linearSize(const FaceData& vdata);
  template < >
  void linearIn(FaceData& a_outputT, const void* const inBuf);
  template < >
  void linearOut(void* const a_outBuf, const FaceData& a_inputT);
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
                const Box&               a_region,
                const Real&              a_dx,
                bool a_hasMoments);

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
           const Box&               a_region,
           const Real&              a_dx,
           bool a_hasMoments);

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


    /// get the moment at  the VoF associated with the  monomial with the input exponents
    /**
       Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
       returns integral_over_VoF(x^p y^q z^r dV) for all p q
    **/
    IndMomSpaceDim getVolumeMoments(const VolIndex& a_vof) const;

    /// get the normal the irregular face associated with the  monomial with the input exponents 
    /**
       Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
       returns integral_over_VoF((x^p y^q z^r) dV) for p q r
    **/
    IndMomSpaceDim getEBMoments(const VolIndex& a_vof) const;

    /// get the moment at  the face associated with the  monomial with the input exponents
    /**
       Given face variables x, y,  p = mono(0), q = mono(1)),
       returns integral_over_face(x^p y^q dA) for all p q
    **/
    IndMomSDMinOne getFaceMoments(const FaceIndex& a_face) const;

    /// I am sick of guessing this name wrong
    IndMomSDMinOne getAreaMoments(const FaceIndex& a_face) const
      {
        return getFaceMoments(a_face);
      }

    /// 
    /**
       gets the partial derivs of the normal component.   
    **/
    IndMomSpaceDim getEBNormalPartialDerivs(const VolIndex& a_vof, int normalComponent) const;

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

    ///is this vof in my irregular (cut cell) data
    bool irregVoF(const VolIndex& a_vol) const;

    ///is this face in my irregular (cut cell) data
    bool irregFace(const FaceIndex& a_face) const;

private:


    void setVolumeMomentsToZero(const VolIndex& a_vof);

    void setAreaMomentsToZero(const FaceIndex& a_face);

    void setCoveredAndRegular();

    void shiftAndIncrement(IndMomSpaceDim& a_output, const IndMomSpaceDim& a_input,const RealVect& a_shiftRV);
    void shiftAndIncrement(IndMomSDMinOne& a_output, const IndMomSDMinOne& a_input,const RealVect& a_shiftRV, int faceDir);

    Box m_region;

    EBGraph m_graph;
    ///
    BaseIVFAB<Real>   m_volData;

    ///
    BaseIFFAB<Real>   m_faceData[SpaceDim];

    ///
    bool m_isDefined;


    //this stuff only defined if the geometryservice was able to provide them
    ///
    BaseIVFAB<VolData>   m_volMoments;

    ///
    BaseIFFAB<FaceData>  m_faceMoments[SpaceDim];

    ///
    IndMomSDMinOne m_regularAreaMoments;
 
    ///
    IndMomSpaceDim m_regularVolumeMoments;

    /// true if this object has higher order moments
    bool m_hasMoments;

    Real m_dx;
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
           const Box&               a_region,
           const Real&              a_dx,
           bool a_hasMoments)
    
      {
        m_implem->define(a_graph, a_irregData, a_validBox, a_region, a_dx, a_hasMoments);
      }

    //this is the define to use when you are just going to copy the data into it.  region should include ghost cells
    void define(const EBGraph&           a_graph,
                const Box&               a_region,
                const Real&              a_dx,
                bool a_hasMoments) 
      {
        m_implem->define(a_graph, a_region, a_dx, a_hasMoments);
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

    /// get the moment at  the face associated with the  monomial with the input exponents
    /**
       Given face variables x, y,  p = mono(0), q = mono(1)),
       returns integral_over_face(x^p y^q dA) for all p q
    **/
    IndMomSDMinOne getFaceMoments(const FaceIndex& a_face) const;

    /// get the moment at  the VoF associated with the  monomial with the input exponents
    /**
       Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
       returns integral_over_VoF(x^p y^q z^r dV) for all p q
    **/
    IndMomSpaceDim getVolumeMoments(const VolIndex& a_vof) const;

    /// get the normal*moment at  the irregular face associated with the  monomial with the input exponents 
    /**
       Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
       returns integral_over_VoF((x^p y^q z^r) dV) for p q r
    **/
    IndMomSpaceDim getEBMoments(const VolIndex& a_vof) const;

    /// get the partial derivatives at  the irregular face associated with the  monomial with the input exponents 
    /**
       Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
       returns partial_(x^p y^q z^r)(n_i)  for p q r
    **/
    IndMomSpaceDim getEBNormalPartialDerivs(const VolIndex& a_vof, int normalComponent) const;

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
  inline  IndMomSDMinOne EBData::getFaceMoments(const FaceIndex& a_face) const
  {
    return m_implem->getFaceMoments(a_face);
  }
/*******************************/
  inline  IndMomSpaceDim EBData::getVolumeMoments(const VolIndex& a_vof) const
  {
    return m_implem->getVolumeMoments(a_vof);
  }
/*******************************/
  inline  IndMomSpaceDim EBData::getEBMoments(const VolIndex& a_vof) const
  {
    return m_implem->getEBMoments(a_vof);
  }

/*******************************/
  inline  IndMomSpaceDim EBData::getEBNormalPartialDerivs(const VolIndex& a_vof, int a_normalComponent) const
  {
    return m_implem->getEBNormalPartialDerivs(a_vof, a_normalComponent);
  }
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
