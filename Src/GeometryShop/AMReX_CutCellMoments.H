
#ifndef _CUTCELLMOMENTS_H_
#define _CUTCELLMOMENTS_H_


#include <map>
using std::map;

#include "AMReX_Vector.H"
#include "AMReX_REAL.H"
#include "AMReX_IndexTM.H"
#include "AMReX_Notation.H"
#include "AMReX_IFData.H"

namespace amrex
{

  template <int dim> class CutCellMoments
  {
  public:
    typedef IndexTM<int,dim>  IvDim;
    typedef IndexTM<Real,dim> RvDim;

    typedef map<IvDim,Real >  PthMoment;

    typedef map<IndexTM<int,dim-1>,Real > PthMomentLesserDimension;

    typedef map<IndexTM<int,1>,Real >     OneDMoments;
    typedef map<int,IvDim>                LocPthMoment;
    typedef map<IvDim,int > PthMomentLoc;

    typedef map<Iv2,CutCellMoments<dim-1> > BdCutCellMoments;

    // Constructors
    CutCellMoments();
    CutCellMoments(const CutCellMoments<dim>& a_cutCellMoments);

    // This is used to build bd CutCellMoments
    CutCellMoments(const IFData<dim>& a_info);

    // Destructor
    ~CutCellMoments();

    // Returns a boundary element of one dim less
    const CutCellMoments<dim - 1> getBdCutCellMoments(const Iv2& a_bdId) const;

    Real changeMomentCoordinates(PthMoment               & a_refinedMomentMap,
                                 const IndexTM<int,dim>  & a_monomial,
                                 const IndexTM<Real,dim> & a_refinedCenterDelta);

    void changeMomentCoordinatesToCellCenter();

    void changeMomentCoordinatesToParentCenter();

    void initialize(CutCellMoments<dim> & a_refinedCutCell);

    void initializeMap(PthMoment & a_map1,
                       PthMoment & a_map2);

    void initializeMap(PthMomentLesserDimension & a_map1,
                       PthMomentLesserDimension & a_map2);

  
    // Get the value of the moment depending on the values of the booleans
    // m_allVerticesIn/Out these functions are used under refinement, when
    // moments on the faces are needed even on a covered/regular cell
    Real getBdMoment(const IvDim             & a_mono,
                     const IFData<dim+1>     & a_IFData,
                     const IndexTM<Real,dim> & a_refinedCenterDelta,
                     PthMoment                 a_fullCellMap = PthMoment());

    Real getBdEBMoment(const IvDim             & a_mono,
                       const IFData<dim+1>     & a_IFData,
                       const IndexTM<Real,dim> & a_refinedCenterDelta);

    void addBdMoments(CutCellMoments<dim>     & a_coarseCutCell,
                      const IFData<dim+1>     & a_IFData,
                      const int               & a_degreePmax,
                      const bool              & a_useConstraints,
                      const IndexTM<Real,dim> & a_refinedCenterDelta,
                      const IndexTM<int,dim>  & a_localHilo);

    // Integrates a monomial over a full cell
    Real fullCellQuadrature(const IndexTM<int,dim>      & a_mono,
                            const CoordinateSystem<dim> & a_coord);

    // Output methods that check a_mono is in the map
    Real getMoment(const IvDim   & a_mono,
                   const EBorVol & a_EBorVOL) const;

    // Methods for reading geom data that do sanity checks on results
    Real getVol(const EBorVol & a_EBorVol) const;

    RvDim getCentroid(const EBorVol & a_EBorVOL) const;

    Real getResidual(const int & a_iDegree,
                     const int & a_normJ) const;

    void setResidual(const Real& a_value,
                     const int & a_iDegree,
                     const int & a_normJ);

    Vector<Real> sliceResidual(const int & a_iDegree) const;

    bool isCovered() const;

    bool isRegular() const;

    // Output
    void print(ostream& out) const;

    void dump() const;

    // Operators
    void operator=(const CutCellMoments<dim>& a_cutCellMoments);

    //volume moments
    PthMoment             m_moments;

    //eb moments
    PthMoment             m_EBmoments;

    //lower dimensional cut cells
    BdCutCellMoments      m_bdCutCellMoments;

    //edge interesections,normals and derivatives of normals
    IFData<dim>           m_IFData;

    //indicates that a boundary CutCellMoment coincides with the interface
    bool                  m_bdCCOn;

    //residual from the least squares problem at the highest dimension
    Vector<Vector<Real> > m_residual;

    //number of active constraints
    int                   m_numActiveBounds;

    // records whether this cutCellMoment or any elements of BdCutCellMoments are using the zero vector for a normal
    bool                  m_badNormal;
  };

// One dimensional cutCellMoments
  template <> class CutCellMoments<1>
  {
  public:
    typedef map<IndexTM<int,1>,Real> OneDMoments;

    // Constructors
    CutCellMoments();
    CutCellMoments(const CutCellMoments<1> & a_cutCellMoments);

    CutCellMoments(const IFData<1>& a_info);

    // Destructor
    ~CutCellMoments();

    Real changeMomentCoordinates(OneDMoments           & a_refinedMap,
                                 const IndexTM<int,1>  & a_monomial,
                                 const IndexTM<Real,1> & a_refinedCenterDelta);

    void changeMomentCoordinatesToCellCenter();

    void changeMomentCoordinatesToParentCenter();

    void initialize(CutCellMoments<1> & a_refinedCutCell);

    void initializeMap(OneDMoments & a_map1,
                       OneDMoments & a_map2);

    Real getBdMoment(const IndexTM<int,1>  & a_mono,
                     const IFData<2>       & a_IFData,
                     const IndexTM<Real,1> & a_refinedCenterDelta,
                     OneDMoments             a_fullCellMap = OneDMoments());

    Real getBdEBMoment(const IndexTM<int,1>  & a_mono,
                       const IFData<2>       & a_IFData,
                       const IndexTM<Real,1> & a_refinedCenterDelta);

    void addBdMoments(CutCellMoments<1>     & a_coarseCutCell,
                      const IFData<2>       & a_IFData,
                      const int             & a_degreePmax,
                      const bool            & a_useConstraints,
                      const IndexTM<Real,1> & a_refinedCenterDelta,
                      const IndexTM<int,1>  & a_localHilo);

    // Output method that check a_mono is in the map
    Real getMoment(const IndexTM<int,1> & a_mono,
                   const EBorVol        & a_EBorVOL) const;

    Real getMoment(const IndexTM<int,1>& a_mono) const;

    // Methods for reading geom data that do sanity checks on results
    Real getVol(const EBorVol& a_EBorVol) const;

    IndexTM<Real,1> getCentroid(const EBorVol& a_EBorVol) const;

    bool isCovered() const;

    bool isRegular() const;

    // Output
    void print(ostream& out) const;

    void dump() const;

    // Operators
    void operator=(const CutCellMoments<1>& a_cutCellMoments);

    // Member data
    OneDMoments  m_moments;
    IFData<1>    m_IFData;
    bool         m_bdCCOn;
    int          m_numActiveBounds;
    bool         m_badNormal;

    // This is a dummy quantity: no EB in 1D
    OneDMoments  m_EBmoments;
  };

}
#include "AMReX_CutCellMomentsImplem.H"

#endif
