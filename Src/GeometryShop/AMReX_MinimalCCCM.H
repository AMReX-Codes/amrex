
#ifndef _MINIMALCCCM_H_
#define _MINIMALCCCM_H_

#include <map>
using std::map;

#include "AMReX_Vector.H"
#include "AMReX_REAL.H"
#include "AMReX_IndexTM.H"
#include "AMReX_Notation.H"
#include "AMReX_LSquares.H"
#include "AMReX_IFData.H"
#include "AMReX_CutCellMoments.H"
#include "AMReX_RefinementCriterion.H"

namespace amrex
{
  template <int dim> class LSProblem;

  template <int dim> class MinimalCCCM
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
    MinimalCCCM();
    MinimalCCCM(const MinimalCCCM<dim>& a_MinimalCCCM);

    // This is used to build bd MinimalCCCM
    MinimalCCCM(const IFData<dim>& a_info);

    // Destructor
    ~MinimalCCCM();

    void computeMoments(const int                & a_orderPmax,
                        const int                & a_degreePmax);

    void computeMomentsRecursively(const int                & a_orderPmax,
                                   const int                & a_degreePmax);

    Vector<Real> computeRhs(LSProblem<dim> & a_lsp,
                            const int      & a_order);


    // Output
    void print(ostream& out) const;

    void dump() const;

    // Operators
    void operator=(const MinimalCCCM<dim>& a_MinimalCCCM);

    Real factorial(const IvDim & a_multiIndex) const;

    // All the moments (and other data)
    CutCellMoments<dim> m_cutCellMoments;

    // Flag whether the boundary moments have been computed already
    bool m_boundaryMomentsComputed;
  };

// One dimensional MinimalCCCM
  template <> class MinimalCCCM<1>
  {
  public:
    typedef map<IndexTM<int,1>,Real> OneDMoments;

    // Constructors
    MinimalCCCM();
    MinimalCCCM(const MinimalCCCM<1> & a_MinimalCCCM);
    MinimalCCCM(const IFData<1>& a_info);

    // Destructor
    ~MinimalCCCM();

    void computeMoments(const int              & a_orderPmax,
                        const int              & a_degreePmax);


    void simpleComputeMoments(const Real & a_loPt,
                              const Real & a_hiPt,
                              const int  & a_degreePmax);

    void computeMomentsUsingBinomial(const Real & a_loPt,
                                     const Real & a_hiPt,
                                     const int  & a_loSign,
                                     const int  & a_hiSign,
                                     const int  & a_degreePmax);

    // Output
    void print(ostream& out) const;

    void dump() const;

    // Operators
    void operator=(const MinimalCCCM<1>& a_MinimalCCCM);

    // Member data
    CutCellMoments<1> m_cutCellMoments;
  };


}

#include "AMReX_MinimalCCCMImplem.H"

#endif
