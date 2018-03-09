

#ifndef _LSPROBLEM_H_
#define _LSPROBLEM_H_

#include <map>
using std::map;

#include "AMReX_Vector.H"
#include "AMReX_REAL.H"
#include "AMReX_LSquares.H"
#include "AMReX_IndexTM.H"
#include "AMReX_Notation.H"
#include "AMReX_IFData.H"


namespace amrex
{

  template <int dim> class CutCellMoments;

  template <int dim> class LSProblem
  {
    typedef IndexTM<int,dim>IvDim;
    typedef map<IvDim,int > PthMomentLoc;
    typedef map<int,IvDim> LocPthMoment;

  public:

    //destructor will free m_matrix and free matrixes in m_rhs
    ~LSProblem();

    //constructor for just accessing fill map function
    LSProblem(const int  & a_degreeP,
              const bool & a_useConstraints);

    //constructor
    LSProblem(const int                & a_orderAccuracy,
              const int                & a_degreeP,
              const bool               & a_useConstraints,
              const IndexTM<Real,dim>  & a_normal);

    //Solves over determined system matrix*unknown = rhs
    int invertNormalEq(const Vector<Real> & a_rhs,
                       Vector<Real>       & a_residual);


    //accessor functions
    int getDegree()
      {
        return m_degreeP;
      }

    int getOrderAccuracy()
      {
        return m_order;
      }

    void getMatrix(Real ** a_matrix);

    void getRhs(Vector<Real>& a_rhs);

    //getUnknowns will return the geometry info at the current level
    void getUnknowns(Vector<Real>& a_unknown);

    // Compute max and min values of monomial over the cell in local coordinates
    void monoMaxMin(Real                   & a_maxVal,
                    Real                   & a_minVal,
                    const IndexTM<int,dim> & a_mono,
                    const IFData<dim>      & a_IFData);

    void computeBounds(const IndexTM<Real,dim>   & a_dx,
                       const CutCellMoments<dim> & a_ccm);

    void print(ostream & a_out) const;

    void outputMatrix() const;
    void outputRhs() const;
    void outputUnknowns() const;
    void outputBounds() const;

    const LocPthMoment& getLocMonomialMapDegreeP() const
      {
        return (const LocPthMoment&)m_locMonoP;
      }

    const PthMomentLoc& getMonomialLocMapDegreeP() const
      {
        return (const PthMomentLoc&)m_monoLocP;
      }

    const PthMomentLoc& getMonomialLocMapDegreePLess1() const
      {
        return (const PthMomentLoc&)m_monoLocPLess1;
      }

    Real getUnknown(int loc)
      {
        return m_unknowns[loc];
      }

    int getNumberDegP()
      {
        return m_numP;
      }

    int getNumberDegPLess1()
      {
        return m_numPLess1;
      }

    int numActiveBounds() const
      {
        return m_numActiveBounds;
      }




    void setRhs(const Vector<Real>& a_rhs); //deprecated

    //constructs matrix of overdetermined system
    void setMatrix();

    void momentBounds(Real              & a_lobnd,
                      Real              & a_hibnd,
                      const IvDim       & a_mono,
                      const IFData<dim> & a_ifData);

    //\partial/\partial xidir (mono) = coeff*Dmono
    void differentiate(int        & a_coeff,
                       IvDim      & a_Dmono,
                       int        & a_idir,
                       const IvDim& a_mono);

    int nChooseR(int a_n,
                 int a_r);

    //makes a list of monomials of degreeP
    void fillMap(PthMomentLoc & a_monoLoc,
                 LocPthMoment & a_locMono,
                 const int    & a_degree);

    //computes number of monomials of degree a_monoDegree
    int numMonomials(const int & a_monoDegree);

    int factorial(const int & a_n,
                  const int & a_m=0);

    //array tools
    void allocArray(const int & a_rows,
                    const int & a_cols,
                    Real**    & a_A);

    void freeArray(const int & a_rows,
                   const int & a_cols,
                   Real**    & a_A);

    //order of accuracy
    int m_order;

    //degree of monomial
    int m_degreeP;

    //number of active constraints
    int m_numActiveBounds;

    //use constrained LS
    bool m_useConstraints;

    //normal
    IndexTM<Real,dim> m_normal;

    //monomials of degree P
    PthMomentLoc m_monoLocP;
    LocPthMoment m_locMonoP;

    //number of monomials of degreeP
    int m_numP;

    //monomials of degree P-1
    PthMomentLoc m_monoLocPLess1;
    LocPthMoment m_locMonoPLess1;

    //number of monomials of degreeP
    int m_numPLess1;


    //matrix
    Real**  m_matrix;
    Vector<Real> m_unknowns;
    Vector<Real> m_rhs;

    Vector<Real> m_lowerBound;
    Vector<Real> m_upperBound;
  };

  template<>class LSProblem<1>
  {
  public:

    //number of monomials of degree P
    int m_numP;

    //number of monomials of degree P-1
    int m_numPLess1;

    //monomials of degree P
    Real m_pMoments;

    //used in filling the map
    int m_ithMoment;

    //Destructor
    ~LSProblem();

    //copy constructor
    LSProblem(const LSProblem<1> & a_lsProblem);

    //empty Constructor
    LSProblem();

    //count num monomials
    int recursiveCount(const int & a_degreeP);

    void setNumMonomials();

    //print
    void print(ostream & a_out)const;

    //equals operator
    void operator=(const LSProblem & a_lSProblem);
  };

}

#include "AMReX_LSProblemImplem.H"

#endif
