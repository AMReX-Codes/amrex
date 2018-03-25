// We can assume that template class T has null construction.

#ifndef _AMREBMULTIGRID_H_
#define _AMREBMULTIGRID_H_

#include "AMReX_EBMultiGrid.H"
#include "AMReX_NoOpSolver.H"
#include "AMReX_REAL.H"
#include "AMReX_Box.H"
#include "AMReX_SPMD.H"
#include "AMReX_BLProfiler.H"

#include <iomanip>
using std::vector;
namespace amrex
{
///
/**
   Operator class for AMR EBMultiGrid
*/
  template <typename T>
  class AMRLevelOp : public MGLevelOp<T>
  {
  public:


    //! Constructor.
    AMRLevelOp()
      :MGLevelOp<T>()
      {
      }


    //! Destructor.
    virtual ~AMRLevelOp()
      {
      }

    virtual Real AMRNorm(const T&   a_coarResid,
                         const T&   a_fineResid,
                         const int& a_refRat,
                         const int& a_ord)
      {
        return this->norm(a_coarResid, 0);
      }


    ///
    /**
       a_residual = a_rhs - L(a_phiFine, a_phi, a_phiCoarse)
       The objects are undefined when there is nothing in that direction.
       (i.e. no finer means a_phiFine is undefined.)
      
    */
    virtual void AMRResidual(T& a_residual, const T& a_phiFine, const T& a_phi,
                             const T& a_phiCoarse, const T& a_rhs,
                             bool a_homogeneousDomBC,
                             AMRLevelOp<T>*  a_finerOp) = 0;

    ///
    /**
       a_residual = a_rhs - L^nf(a_phi, a_phiCoarse)
       assume no finer AMR level
    */
    virtual void AMRResidualNF(T& a_residual, const T& a_phi, const T& a_phiCoarse,
                               const T& a_rhs, bool a_homogeneousBC) = 0;

    ///
    /**
       a_residual = a_rhs - L(a_phiFine, a_phi)
       assume no coarser AMR level
    */
    virtual void AMRResidualNC(T& a_residual, const T& a_phiFine, const T& a_phi,
                               const T& a_rhs, bool a_homogeneousBC,
                               AMRLevelOp<T>* a_finerOp) = 0;
    ///
    /**
       Apply the AMR operator, including coarse-fine matching
    */
    virtual void AMROperator(T& a_LofPhi,
                             const T& a_phiFine, const T& a_phi,
                             const T& a_phiCoarse,
                             bool a_homogeneousDomBC,
                             AMRLevelOp<T>*  a_finerOp) = 0;

    ///
    /** a_resCoarse = I[h-2h]( a_residual - L(a_correction, a_coarseCorrection)) */
    virtual void AMRRestrict(T& a_resCoarse, const T& a_residual, const T& a_correction,
                             const T& a_coarseCorrection) = 0;

    ///
    /** a_correction += I[2h->h](a_coarseCorrection) */
    virtual void AMRProlong(T& a_correction, const T& a_coarseCorrection) = 0;

    ///
    /** a_residual = a_residual - L(a_correction, a_coarseCorrection) */
    virtual void AMRUpdateResidual(T& a_residual, const T& a_correction,
                                   const T& a_coarseCorrection) = 0;

    ///
    /**
     */
    virtual void createCoarsened(T        & a_lhs,
                                 const T  & a_rhs,
                                 const int& a_refRat) = 0;


    //a_rhs holds a coarsened-fine residual.   the zero and assign is designed to
    //fill the bit that is covered by a finer grid with zeros.
    virtual void zeroCovered(T& a_lhs, T& a_rhs)
      {
        this->setToZero(a_rhs);
        this->assign(a_lhs, a_rhs);
      }
    virtual Real localMaxNorm(const T& a_phi)
      {
        return this->norm(a_phi, 0);
      }

  };

///
/**
   Factory to create AMRLevelOps
*/
  template <class T>
  class AMRLevelOpFactory : public MGLevelOpFactory<T>
  {
  public:
    virtual ~AMRLevelOpFactory()
      {
      }

    ///
    /**
       return a new operator.  this is done with a new call.
       caller is responsible for deletion
    */
    virtual AMRLevelOp<T>* AMRnewOp(const Box& a_indexSpace)=0;

  };


///
/**
   Class to solve elliptic equations using the Martin and Cartwright algorithm.
*/
  template <class T>
  class AMREBMultiGrid
  {
  public:

    AMREBMultiGrid();

    virtual ~AMREBMultiGrid();

    ///
    /**
       Define the solver.
       a_coarseDomain is the index space on the coarsest AMR level.
       a_factory is the operator factory through which all special information is conveyed.
       a_bottomSolver is the solver to be used at the termination of EBMultiGrid coarsening.
       It is the client's responsibility to free up the dynamically-allocated memory.
       a_numLevels is the number of AMR levels.
       If enableLevelSolves= false, only works for lbase==0 (saves a bunch of definition
       time if you are only ever solving over the whole hierarchy).
    */
    virtual void define(const Box                    & a_coarseDomain,
                        AMRLevelOpFactory<T>         & a_factory,
                        shared_ptr<LinearSolver<T> > & a_bottomSolver,
                        const std::vector<int>       & a_refRat,
                        bool a_enableLevelSolves,
                        int a_numLevels);


    ///
    /**
       Solve L(phi) = rho from l_base to l_max.  To solve over all levels,
       l_base = 0 and l_max = max_level = numLevels-1.
    */
    virtual void solve( std::vector<T*>& a_phi, const std::vector<T*>& a_rhs,
                        int l_max, int l_base, bool a_zeroPhi=true,
                        bool forceHomogeneous = false );

    ///
    /** same as "solve" except user has taken the reponsibility of having previously
        called "init" so solver can allocate temporary holders.
    */
    virtual void solveNoInit(std::vector<T*>& a_phi, const std::vector<T*>& a_rhs,
                             int l_max, int l_base, bool a_zeroPhi=true,
                             bool forceHomogeneous = false);

    ///use if you want final residual
    virtual void solveNoInitResid(std::vector<T*>& a_phi, std::vector<T*>& a_finalResid,
                                  const std::vector<T*>& a_rhs,
                                  int l_max, int l_base, bool a_zeroPhi=true,
                                  bool forceHomogeneous = false);

    virtual void AMRVCycle(std::vector<T*>& a_correction,
                           std::vector<T*>& a_residual,
                           int l, int l_max, int l_base);

    void setMGCycle(int a_numMG);

    virtual void init(const std::vector<T*>& a_phi, const std::vector<T*>& a_rhs,
                      int l_max, int l_base);
    //init messes with EBMultiGrid depth.  this puts it back
    void revert(const std::vector<T*>& a_phi, const std::vector<T*>& a_rhs,
                int l_max, int l_base);

    int getRefToCoarser(int a_depth)
      {
        int retval = 2;
        if(a_depth > 0)
        {
          retval = m_refRat[a_depth-1];
        }
        return retval;
      }

    std::vector<int> m_refRat;
    Real m_eps, m_hang, m_normThresh;
    bool   m_solverParamsSet, m_enableLevelSolves;
    int m_imin, m_iterMax, m_iterMin, m_verbosity, m_exitStatus;
    int m_pre, m_post, m_bottom, m_numMG;
    /// max no. of coarsenings -- -1 (default) means coarsen as far as possible
    /** If using a value besides the default, need to set it _before_
        define function is called */
    int m_maxDepth;
    AMRLevelOp<T>& levelOp(int level);
    // default m_convergenceMetric = 0.:  initial residual will be set to
    // result of computeAMRResidual.
    // if m_convergenceMetric > 0., then initial residual will be set to
    // m_convergenceMetric.
    Real m_convergenceMetric;

    // used to give an additional cushion in the EPS used for bottom solves
    Real m_bottomSolverEpsCushion;

    ///
    /**
       resid = L(phi) - rhs
    */
    Real computeAMRResidual(std::vector<T*>&       a_resid,
                            std::vector<T*>&       a_phi,
                            const std::vector<T*>& a_rhs,
                            int               l_max,
                            int               l_base,
                            bool              a_homogeneousBC=false,
                            bool              a_computeNorm=true);

    /** just return the normed value of computeAMRResidual.  used for benchmarking */
    Real computeAMRResidual(std::vector<T*>&      a_phi,
                            const std::vector<T*>& a_rhs,
                            int l_max,
                            int l_min);

    ///
    /**
       lph = L(phi)
    */
    void computeAMROperator(std::vector<T*>&       a_lph,
                            std::vector<T*>&       a_phi,
                            int               l_max,
                            int               l_base,
                            bool              a_homogeneousBC=false);

    ///
    /**
       For changing coefficients.  Use at thy own peril.
    */
    std::vector< MGLevelOp<T>* > getAllOperators();
    std::vector< MGLevelOp<T>* > getOperatorsOp();
    std::vector< std::vector< MGLevelOp<T>* > > getOperatorsMG();
    std::vector< AMRLevelOp<T> * >& getAMROperators()
      {
        return m_op;
      }

    ///
    /**
       Set parameters of the solve.
       a_pre is the number of smoothings before averaging.
       a_post is the number of smoothings after averaging.
       a_bottom is the number of smoothings at the bottom level.
       a_numMG = 1 for vcycle, =2 for wcycle (use 1).
       a_itermax is the max number of v cycles.
       a_hang is the minimum amount of change per vcycle.
       a_eps is the solution tolerance.
       a_normThresh is how close to zero eps*resid is allowed to get.
    */
    void setSolverParameters(const int&   a_pre,
                             const int&   a_post,
                             const int&   a_bottom,
                             const int&   a_numMG,
                             const int&   a_iterMax,
                             const Real&  a_eps,
                             const Real&  a_hang,
                             const Real&  a_normThresh);

    /// set up bottom solver for internal MG solver
    /**
       This function is normally called by the solve(...) function.
       However, it must be called if solve will not be called (in particular,
       if only the V-cycle is being used)
    */
    void setBottomSolver(int l_max, int l_base);

    void setBottomSolverEpsCushion(Real a_bottomSolverEpsCushion);


  protected:

    void relax(T& phi, T& R, int depth, int nRelax = 2);

    void computeAMRResidualLevel(std::vector<T*>&       a_resid,
                                 std::vector<T*>&       a_phi,
                                 const std::vector<T*>& a_rhs,
                                 int l_max, int l_base, int ilev,
                                 bool a_homogeneousBC);

    std::vector<AMRLevelOp<T>*>          m_op;
    std::vector<EBMultiGrid<T> *>          m_mg;
    std::vector<T*>  m_correction;
    std::vector<T*>  m_residual;
    std::vector<T*>  m_resC;

    shared_ptr<LinearSolver<T> >   m_nosolve;

    shared_ptr<LinearSolver<T> > m_bottomSolver;

    std::vector<char> m_hasInitBeenCalled;

    void clear();

  private:

    // Forbidden copiers.
    AMREBMultiGrid(const AMREBMultiGrid<T>&);
    AMREBMultiGrid& operator=(const AMREBMultiGrid<T>&);
  };

//*******************************************************
// AMREBMultiGrid Implementation
//*******************************************************

//===================================================================


  template <class T>
  void
  AMREBMultiGrid<T>::setSolverParameters(const int&   a_pre,
                                         const int&   a_post,
                                         const int&   a_bottom,
                                         const int&   a_numMG,
                                         const int&   a_iterMax,
                                         const Real&  a_eps,
                                         const Real&  a_hang,
                                         const Real&  a_normThresh)
  {
    m_solverParamsSet = true;
    m_pre        =    a_pre;
    m_post       =    a_post;
    m_bottom     =    a_bottom;
    m_eps        =    a_eps;
    m_hang       =    a_hang;
    m_normThresh =    a_normThresh;
    m_iterMax    =    a_iterMax;
    for (int img = 0; img < m_mg.size(); img++)
    {
      m_mg[img]->m_pre    = a_pre;
      m_mg[img]->m_post   = a_post;
      m_mg[img]->m_bottom   = a_bottom;
    }
    setMGCycle(a_numMG);
    m_bottomSolverEpsCushion = 1.0;
  }
  template <class T>
  std::vector< MGLevelOp<T> * >
  AMREBMultiGrid<T>::getAllOperators()
  {
    std::vector< MGLevelOp<T>* > retval;
    for (int iop = 0;  iop < m_op.size(); iop++)
    {
      MGLevelOp<T>* operPtr = (MGLevelOp<T>*) m_op[iop];
      retval.push_back(operPtr);
    }

    for (int img = 0; img < m_mg.size(); img++)
    {
      std::vector< MGLevelOp<T>* > mgOps = m_mg[img]->getAllOperators();
      retval.append(mgOps);
    }
    return retval;
  }

  template <class T>
  std::vector< MGLevelOp<T> * >
  AMREBMultiGrid<T>::getOperatorsOp()
  {
    std::vector< MGLevelOp<T>* > retval;
    for (int iop = 0;  iop < m_op.size(); iop++)
    {
      MGLevelOp<T>* operPtr = (MGLevelOp<T>*) m_op[iop];
      retval.push_back(operPtr);
    }
    return retval;
  }

  template <class T>
  std::vector<std::vector< MGLevelOp<T> * > >
  AMREBMultiGrid<T>::getOperatorsMG()
  {
    std::vector< std::vector< MGLevelOp<T>* > > retval(m_mg.size());

    for (int img = 0; img < m_mg.size(); img++)
    {
      retval[img] = m_mg[img]->getAllOperators();
    }
    return retval;
  }

  template <class T>
  AMREBMultiGrid<T>::AMREBMultiGrid()
    :m_eps(1E-6),
     m_hang(1E-15),
     m_normThresh(1E-30),
     m_imin(5),
     m_iterMax(20),
     m_iterMin(-1),
     m_verbosity(3),
     m_pre(2),
     m_post(2),
     m_bottom(2),
     m_numMG(1),
     m_maxDepth(-1),
     m_convergenceMetric(0.),
     m_bottomSolverEpsCushion(1.0),
     m_bottomSolver(NULL)
  {
    m_solverParamsSet = false;

    //m_maxDepth = 4;
    //pout() << "       AMREBMultiGrid<T>::AMREBMultiGrid() = " << m_maxDepth << endl;
  }
  template <class T>
  void AMREBMultiGrid<T>::setMGCycle(int a_numMG)
  {
    for (int ilev = 0; ilev < m_op.size(); ilev++)
    {
      m_mg[ilev]->m_numMG = a_numMG;
      m_mg[ilev]->m_cycle = a_numMG;
    }
    m_numMG = a_numMG;
  }
  template <class T>
  void AMREBMultiGrid<T>::relax(T& a_correction, T& a_residual, int depth, int a_numSmooth)
  {
    BL_PROFILE("AMREBMultiGrid::relax");

    int refRat = getRefToCoarser(depth);
    if (refRat > 2)
    {
      // intermediate EBMultiGrid levels exist between this level and the
      // next coarser AMR level, so do a mini v-cyle to smooth on those

      int intermediateDepth = 0;
      int r = refRat;
      while (r>2)
      {
        r/=2;
        intermediateDepth++;
      }
      int tmp = m_mg[depth]->m_depth;
      m_mg[depth]->m_depth = intermediateDepth;
      //ok, there is an intermediate EBMultiGrid level that is
      // not an AMR level.  We use regular EBMultiGrid smoothing algorithm for this.
      // note that EBMultiGrid::cycle handles relaxation on this level
      // as well, so there's no need to call the LinearOp's relax function.
      m_mg[depth]->cycle(0, a_correction, a_residual);
      m_mg[depth]->m_depth = tmp;
    }
    else
    {
      // no intermediate EBMultiGrid levels between AMR levels, so
      // just call relax on this level and be done with it.
      m_op[depth]->relax(a_correction, a_residual, a_numSmooth);  //numSmoothDown
    }

  }
/************/
  template <class T>
  AMREBMultiGrid<T>::~AMREBMultiGrid()
  {
    BL_PROFILE("~AMREBMultiGrid");
    clear();
  }

  template <class T>
  AMRLevelOp<T>&  AMREBMultiGrid<T>::levelOp(int level)
  {
    return *(m_op[level]);
  }

/************/
  template <class T>
  Real AMREBMultiGrid<T>::computeAMRResidual(std::vector<T*>&       a_resid,
                                           std::vector<T*>&       a_phi,
                                           const std::vector<T*>& a_rhs,
                                           int               l_max,
                                           int               l_base,
                                           bool              a_homogeneousBC,
                                           bool              a_computeNorm)
  {
    BL_PROFILE("AMREBMultiGrid::computeAMRResidual");

    Real rnorm = 0;
    Real localNorm = 0;
    for (int ilev = l_base; ilev <= l_max; ilev++)
    {
      //always used at top level where bcs are inhomogeneous
      computeAMRResidualLevel(a_resid,
                              a_phi,
                              a_rhs,
                              l_max, l_base, ilev, a_homogeneousBC);
      if (a_computeNorm)
      {
        if (ilev == l_max)
        {
          localNorm =m_op[ilev]->localMaxNorm(*a_resid[ilev]);
        }
        else
        {
          m_op[ilev]->zeroCovered(*a_resid[ilev], *m_resC[ilev+1]);
          localNorm = m_op[ilev]->localMaxNorm(*a_resid[ilev]);
        }
        const int normType = 0;
        if (normType==2)
          localNorm = m_op[ilev]->norm(*a_resid[ilev],normType);
        rnorm = std::max(localNorm, rnorm);
      }
    }
#ifdef BL_USE_MPI
    if (a_computeNorm)
    {
      BL_PROFILE("MPI_Allreduce");
      Real recv;
      int result = MPI_Allreduce(&rnorm, &recv, 1, MPI_CH_REAL,
                                 MPI_MAX, MPI_COMM_WORLD);
      if (result != MPI_SUCCESS)
      {
        //bark!!!
        amrex::Error("sorry, but I had a communcation error on norm");
      }
      rnorm = recv;
    }
#endif

    return rnorm; // if a_computeNorm is false, then this just returns zero.
  }

  template <class T>
  Real AMREBMultiGrid<T>::computeAMRResidual(std::vector<T*>&       a_phi,
                                           const std::vector<T*>&       a_rhs,
                                           int l_max,
                                           int l_min)
  {
    return computeAMRResidual(m_residual,
                              a_phi,
                              a_rhs,
                              l_max,
                              l_min,
                              false,
                              true);
  }

/************/
// old inneficient way of doing this
/*
  template <class T>
  void AMREBMultiGrid<T>::computeAMROperator(std::vector<T*>&       a_lph,
  std::vector<T*>&       a_phi,
  int               l_max,
  int               l_base,
  bool              a_homogeneousBC)
  {
  
  for (int ilev = l_base; ilev <= l_max; ilev++)
  {
  m_op[ilev]->create(*m_residual[ilev], *a_lph[ilev]);
  m_op[ilev]->setToZero(*m_residual[ilev]);
  }
  computeAMRResidual(a_lph, a_phi, m_residual, l_max, l_base, a_homogeneousBC, false);
  // fixing the bug of not negating the result --- Qinghai Zhang 12/10/2009
  // In loving memory of the two weeks I lost over this bug!!!
  for (int ilev = l_base; ilev <= l_max; ilev++)
  {
  m_op[ilev]->scale(*a_lph[ilev], -1.);
  }
  }
*/
  template <class T>
  void AMREBMultiGrid<T>::computeAMROperator(std::vector<T*>&       a_lphi,
                                           std::vector<T*>&       a_phi,
                                           int               l_max,
                                           int               l_base,
                                           bool              a_homogeneousBC)
  {
    BL_PROFILE("AMREBMultiGrid<T>::computeAMROperator");
 
    for(int ilev=l_base; ilev<=l_max; ilev++)
    { 
      if (l_max != l_base)
      {
        if (ilev == l_max)
        {
          m_op[l_max]->AMROperatorNF(*(a_lphi[l_max]),*(a_phi[l_max]),
                                     *(a_phi[l_max-1]), a_homogeneousBC);
        }
        else if (ilev == l_base)
        {
          if (l_base == 0)
          {
            m_op[l_base]->AMROperatorNC(*(a_lphi[l_base]), *(a_phi[l_base+1]),
                                        *(a_phi[l_base]),
                                        a_homogeneousBC, m_op[l_base+1]);
          }
          else
          {
            m_op[l_base]->AMROperator(*a_lphi[l_base], *a_phi[l_base+1], *a_phi[l_base],
                                      *a_phi[l_base-1],
                                      a_homogeneousBC, m_op[l_base+1]);
          }
        }
        else
        {
          m_op[ilev]->AMROperator(*a_lphi[ilev], *a_phi[ilev+1], *a_phi[ilev],
                                  *a_phi[ilev-1], 
                                  a_homogeneousBC, m_op[ilev+1]);
        }
      }
      else
      {
        BL_ASSERT(ilev == l_base);
        if (l_base == 0)
        {
          m_op[l_max]->applyOp(*a_lphi[l_max], *a_phi[l_max], a_homogeneousBC);
        }
        else
        {
          m_op[l_max]->AMROperatorNF(*(a_lphi[l_max]), *(a_phi[l_max]),
                                     *(a_phi[l_max-1]), 
                                     a_homogeneousBC);
        }
      }
    }
  } 

/************/
  template <class T>
  void AMREBMultiGrid<T>::computeAMRResidualLevel(std::vector<T*>&       a_resid,
                                                std::vector<T*>&       a_phi,
                                                const std::vector<T*>& a_rhs,
                                                int l_max, int l_base, int ilev,
                                                bool a_homogeneousBC)
  {
    BL_PROFILE("AMREBMultiGrid<T>::computeAMRResidualLevel");
    BL_ASSERT(m_hasInitBeenCalled[ilev]=='t');

    //m_op[ilev]->setToZero(*(a_resid[l_max]));
    if (l_max != l_base)
    {
      if (ilev == l_max)
      {
        m_op[l_max]->AMRResidualNF(*(a_resid[l_max]), *(a_phi[l_max]),
                                   *(a_phi[l_max-1]), *(a_rhs[l_max]),
                                   a_homogeneousBC);
      }
      else if (ilev == l_base)
      {
        if (l_base == 0)
        {
          m_op[l_base]->AMRResidualNC(*(a_resid[l_base]), *(a_phi[l_base+1]),
                                      *(a_phi[l_base]),  *(a_rhs[l_base]),
                                      a_homogeneousBC, m_op[l_base+1]);
        }
        else
        {
          m_op[l_base]->AMRResidual(*a_resid[l_base], *a_phi[l_base+1], *a_phi[l_base],
                                    *a_phi[l_base-1], *a_rhs[l_base],
                                    a_homogeneousBC, m_op[l_base+1]);
        }
      }
      else
      {
        m_op[ilev]->AMRResidual(*a_resid[ilev], *a_phi[ilev+1], *a_phi[ilev],
                                *a_phi[ilev-1], *a_rhs[ilev],
                                a_homogeneousBC, m_op[ilev+1]);
      }
    }
    else
    {
      BL_ASSERT(ilev == l_base);
      if (l_base == 0)
      {
        m_op[l_max]->residual(*a_resid[l_max], *a_phi[l_max], *a_rhs[l_max],a_homogeneousBC);
      }
      else
      {
        m_op[l_max]->AMRResidualNF(*(a_resid[l_max]), *(a_phi[l_max]),
                                   *(a_phi[l_max-1]), *(a_rhs[l_max]),
                                   a_homogeneousBC);
      }
    }

  }

  template<class T>
  void AMREBMultiGrid<T>::solve(std::vector<T*>& a_phi, const std::vector<T*>& a_rhs,
                              int l_max, int l_base, bool a_zeroPhi,
                              bool a_forceHomogeneous)
  {
    BL_PROFILE("AMREBMultiGrid::solve");
    init(a_phi, a_rhs, l_max, l_base);
    solveNoInit(a_phi, a_rhs, l_max, l_base, a_zeroPhi, a_forceHomogeneous);
    //puts EBMultiGrid depth back where it started so the solver can be reused
    revert(a_phi, a_rhs, l_max, l_base);
  }

  template<class T>
  void AMREBMultiGrid<T>::solveNoInit(std::vector<T*>& a_phi,
                                    const std::vector<T*>& a_rhs,
                                    int l_max, int l_base, bool a_zeroPhi,
                                    bool a_forceHomogeneous)
  {
    BL_PROFILE("AMREBMultiGrid::solveNoInit");
    std::vector<T*> uberResidual(a_rhs.size());
    int lowlim = l_base;
    if (l_base > 0)  // we need an ubercorrection one level lower than l_base
      lowlim--;
    // for (int ilev = l_base; ilev <= l_max; ilev++)
    for (int ilev = lowlim; ilev <= l_max; ilev++)
    {
      uberResidual[ilev] = new T();
      if (ilev >= l_base)
      {
        m_op[ilev]->create(*uberResidual[ilev], *a_rhs[ilev]);
      }
    }
    solveNoInitResid(a_phi, uberResidual, a_rhs, l_max, l_base, a_zeroPhi, a_forceHomogeneous);
    for (int i = lowlim; i <= l_max; i++)
    {
      m_op[i]->clear(*uberResidual[i]);
      delete uberResidual[i];
    }
  }

  template<class T>
  void AMREBMultiGrid<T>::solveNoInitResid(std::vector<T*>& a_phi,  std::vector<T*>& uberResidual,
                                         const std::vector<T*>& a_rhs,
                                         int l_max, int l_base, bool a_zeroPhi,
                                         bool a_forceHomogeneous)
  {
    BL_PROFILE("AMREBMultiGrid::solveNo-InitResid");
    setBottomSolver(l_max, l_base);

    BL_ASSERT(l_base <= l_max);
    BL_ASSERT(a_rhs.size() == a_phi.size());

    //these correspond to the residual and correction
    //that live in AMRSolver
    std::vector<T*> uberCorrection(a_rhs.size());

    int lowlim = l_base;
    if (l_base > 0)  // we need an ubercorrection one level lower than l_base
      lowlim--;
    // for (int ilev = l_base; ilev <= l_max; ilev++)

    for (int ilev = lowlim; ilev <= l_max; ilev++)
    { BL_PROFILE("uberCorrection");
      uberCorrection[ilev] = new T();
      m_op[ilev]->create(*uberCorrection[ilev], *a_phi[ilev]);
      if (ilev >= l_base)
      {
        m_op[ilev]->create(*uberResidual[ilev], *a_rhs[ilev]);
        m_op[ilev]->setToZero(*(uberResidual[ilev]));
      }
      m_op[ilev]->setToZero(*(uberCorrection[ilev]));
    }

    if (a_zeroPhi)
      for (int ilev = l_base; ilev <=l_max; ++ilev)
      {
        m_op[ilev]->setToZero(*(a_phi[ilev]));
      }
    //compute initial residual and initialize internal residual to it

    Real initial_rnorm = 0;
    {
      BL_PROFILE("Initial AMR Residual");
      initial_rnorm=computeAMRResidual(uberResidual, a_phi, a_rhs, l_max, l_base, a_forceHomogeneous, true);
    }

    if (m_convergenceMetric != 0.)
    {
      initial_rnorm = m_convergenceMetric;
    }

    Real rnorm = initial_rnorm;
    Real norm_last = 2*initial_rnorm;

    /// set bottom solver convergence norm and solver tolerance
    m_bottomSolver->setConvergenceMetrics(initial_rnorm, m_bottomSolverEpsCushion*m_eps);

    int iter=0;
    if (m_verbosity >= 2) ////
    {
      pout() << "    AMREBMultiGrid:: iteration = " << iter << ", residual norm = " << rnorm << std::endl;
    }

    bool goNorm = rnorm > m_normThresh;                        //iterate if norm is not small enough
    bool goRedu = rnorm > m_eps*initial_rnorm;                 //iterate if initial norm is not reduced enough
    bool goIter = iter < m_iterMax;                            //iterate if iter < max iteration count
    bool goHang = iter < m_imin || rnorm <(1-m_hang)*norm_last;//iterate if we didn't hang
    bool goMin = iter < m_iterMin ; // iterate if iter < min
    while (goMin || (goIter && goRedu && goHang && goNorm))
    {

      norm_last = rnorm;

      //this generates a correction from the current residual

      AMRVCycle(uberCorrection, uberResidual, l_max, l_max, l_base);


      //increment phi by correction and reset correction to zero
      for (int ilev = l_base; ilev <= l_max; ilev++)
      {

        m_op[ilev]->incr(*(a_phi[ilev]), *(uberCorrection[ilev]), 1.0);


        m_op[ilev]->setToZero(*(uberCorrection[ilev]));
        // clean const out
      }

      //recompute residual
      rnorm = computeAMRResidual(uberResidual, a_phi, a_rhs, l_max, l_base, a_forceHomogeneous, true);
      iter++;
      if (m_verbosity >= 3) ////
      {
        pout() << "    AMREBMultiGrid:: iteration = " << iter << ", residual norm = " << rnorm;
        if (rnorm > 0.0)
        {
          pout() << ", rate = " << norm_last/rnorm;
        }
        pout() << std::endl;
      }

      goNorm = rnorm > m_normThresh;                        //keep iterating if norm is not small enough
      goRedu = rnorm > m_eps*initial_rnorm;                 //keep iterating if initial norm is not reduced enough
      goIter = iter < m_iterMax;                            //keep iterating if iter < max iteration count
      goHang = iter < m_imin || rnorm <(1-m_hang)*norm_last;//keep iterating if we didn't hang
      goMin = iter < m_iterMin ; // keep iterating if iter < min
    }
    // if ((rnorm > 10.*initial_rnorm) && (rnorm > 10.*m_eps))
    //   {
    //     pout() << "solver seems to have blown up" << endl;
    //     MayDay::Error("kaboom");
    //   }
    m_exitStatus = int(!goRedu) + int(!goIter)*2 + int(!goHang)*4 + int(!goNorm)*8;
    if (m_verbosity >= 2)
    {
      pout() << "    AMREBMultiGrid:: iteration = " << iter << ", residual norm = " << rnorm << std::endl;
    }
    if (m_verbosity > 1)
    {
      if (!goIter && goRedu && goNorm) // goRedu=T, goIter=F, goHang=?, goNorm=T
      { // m_exitStatus == 0 + 2 + 0|4 + 0 = 2|6
        pout() << "    AMREBMultiGrid:: WARNING: Exit because max iteration count exceeded" << std::endl;
      }
      if (!goHang && goRedu && goNorm) // goRedu=T, goIter=?, goHang=F, goNorm=T
      { // m_exitStatus == 0 + 0|2 + 4 + 0 = 4|6
        pout() << "    AMREBMultiGrid:: WARNING: Exit because of solver hang" << std::endl;
      }
      if (m_verbosity > 4)
      {
        pout() << "    AMREBMultiGrid:: exitStatus = " << m_exitStatus << std::endl;
      }
    }
    for (int i = lowlim; i <= l_max; i++)
    {
      m_op[i]->clear(*uberCorrection[i]);
      delete uberCorrection[i];
    }
  }


  template <class T>
  void AMREBMultiGrid<T>::clear()
  {
    for (int i = 0; i < m_op.size(); i++)
    {
      m_op[i]->clear(*m_correction[i]);
      m_op[i]->clear(*m_residual[i]);
      m_op[i]->clear(*m_resC[i]);

      delete m_correction[i];
      delete m_residual[i];
      delete m_resC[i];
      delete m_op[i];
      delete m_mg[i];

      m_correction[i] = NULL;
      m_residual[i] = NULL;
      m_resC[i] = NULL;
      m_op[i] = NULL;
      m_mg[i] = NULL;
    }
    m_hasInitBeenCalled.resize(0);
  }

  template <class T>
  void AMREBMultiGrid<T>::init(const std::vector<T*>& a_phi, const std::vector<T*>& a_rhs,
                             int l_max, int l_base)
  {
    //if (m_hasInitBeenCalled[l_max]=='t' && m_hasInitBeenCalled[l_base]=='t') return;
    BL_PROFILE("AMREBMultiGrid::init");
//   BL_ASSERT(a_phi.size()>=m_op.size());
//   BL_ASSERT(a_rhs.size()>=m_op.size());
    for (int i = l_base; i <= l_max; i++)
    {
      m_hasInitBeenCalled[i] = 't';
      AMRLevelOp<T>& op = *(m_op[i]);

      op.create(*m_correction[i], *a_phi[i]);
      op.create(*m_residual[i],   *a_rhs[i]);
      int r = getRefToCoarser(i);
      if (i!= l_base)
      {
        int intermediateDepth = 0;

        while (r>2)
        {
          r/=2;
          intermediateDepth++;
        }
        m_mg[i]->m_depth = intermediateDepth;
      }
      m_mg[i]->init(*a_phi[i],     *a_rhs[i]);

      if (i != l_base)
      {
        r = getRefToCoarser(i);
        op.createCoarsened(*m_resC[i], *a_rhs[i], r);
      }
    }
  }

  template <class T>
  void AMREBMultiGrid<T>::revert(const std::vector<T*>& a_phi, const std::vector<T*>& a_rhs,
                               int l_max, int l_base)
  {
    BL_PROFILE("AMREBMultiGrid::revert");
    for (int i = l_base; i <= l_max; i++)
    {
      if (i!= l_base)
      {
        m_mg[i]->m_depth = m_mg[i]->m_defaultDepth;
      }
    }
  }

  template <class T>
  void AMREBMultiGrid<T>::setBottomSolver(int l_max, int l_base)
  {
    //  for (int ilev = 0; ilev <= m_mg.size(); ilev++)
    for (int ilev = l_base; ilev <= l_max; ilev++)
    {
      // only do this if we haven't already set these
      if (!m_solverParamsSet)
      {
        m_mg[ilev]->m_pre = m_pre;
        m_mg[ilev]->m_post = m_post;
        m_mg[ilev]->m_bottom = m_bottom;
      }
      m_mg[ilev]->m_bottomSolver = m_nosolve;
    }

    m_mg[l_base]->setBottomSolver(m_bottomSolver);
  }

  template <class T>
  void AMREBMultiGrid<T>::setBottomSolverEpsCushion(Real a_bottomSolverEpsCushion)
  {
    BL_ASSERT((a_bottomSolverEpsCushion >  0.) &&
              (a_bottomSolverEpsCushion <= 1.));
    m_bottomSolverEpsCushion = a_bottomSolverEpsCushion;
  }

  template <class T>
  void AMREBMultiGrid<T>::define(const Box& a_coarseDomain,
                               AMRLevelOpFactory<T>& a_factory,
                               shared_ptr<LinearSolver<T> > & a_bottomSolver,
                               const vector<int>      & a_refRat,
                               bool a_enableLevelSolves,
                               int a_maxAMRLevels)
  {
    BL_PROFILE("AMREBMultiGrid::define");
    this->clear();
    m_refRat = a_refRat;
    m_enableLevelSolves = a_enableLevelSolves;
    m_op.resize( a_maxAMRLevels, NULL);
    m_mg.resize( a_maxAMRLevels, NULL);
    m_hasInitBeenCalled.resize(a_maxAMRLevels, 'f');
    m_nosolve  = shared_ptr<LinearSolver<T> >(static_cast<LinearSolver<T>*>(new NoOpSolver<T>));
    m_correction.resize( a_maxAMRLevels, NULL);
    m_residual.  resize( a_maxAMRLevels, NULL);
    m_resC.      resize( a_maxAMRLevels, NULL);
    m_bottomSolver = a_bottomSolver;

    Box current = a_coarseDomain;
    for (int i = 0; i < a_maxAMRLevels; i++)
    {
      m_correction[i] = new T();
      m_residual[i]   = new T();
      m_resC[i]       = new T();
      m_mg[i]= new EBMultiGrid<T>();
      m_op[i] = a_factory.AMRnewOp(current);
      int refToCoarser = 2;
      if(i > 0)
      {
        refToCoarser = a_refRat[i-1];
      }
      if((i==0) || (refToCoarser > 2) || a_enableLevelSolves)
      {
        m_mg[i]->define(a_factory, m_nosolve, current, m_maxDepth);
      }

      // Only do this if it will be used (avoiding a reference to invalid
      // and/or unavailable refinement ratios)
      if (i < a_maxAMRLevels-1)
      {
        current.refine(a_refRat[i]);
      }
    }
  }

  	 
  template<class T>
  void AMREBMultiGrid<T>::AMRVCycle(std::vector<T*>& a_uberCorrection,
                                  std::vector<T*>& a_uberResidual,
                                  int ilev, int l_max, int l_base)
  {
    if (ilev == l_max)
    {
      for (int level = l_base; level <= l_max; level++)
      {
        m_op[level]->assign(*m_residual[level], *a_uberResidual[level]);
        m_op[level]->setToZero(*m_correction[level]);
      }
    }

    if (l_max == l_base)
    {
      BL_ASSERT(ilev == l_base);
      m_mg[l_base]->oneCycle(*(a_uberCorrection[ilev]), *(a_uberResidual[ilev]));
    }
    else if (ilev == l_base)
    {
      m_mg[l_base]->oneCycle(*(m_correction[ilev]), *(m_residual[ilev]));
      m_op[ilev]->incr(*(a_uberCorrection[ilev]), *(m_correction[ilev]), 1.0);
    }
    else
    {
      //============= Downsweep ========================

      this->relax(*(m_correction[ilev]), *(m_residual[ilev]), ilev, m_pre);
      m_op[ilev]->incr(*(a_uberCorrection[ilev]), *(m_correction[ilev]), 1.0);

      // Set next coarser level correction to zero
      m_op[ilev-1]->setToZero(*(m_correction[ilev-1]));

      // Recompute residual on next coarser level
      //  for the valid region NOT covered by this level.
      computeAMRResidualLevel(m_residual,
                              a_uberCorrection,
                              a_uberResidual,
                              l_max, l_base, ilev-1,
                              true);

      // Compute the restriction of the residual to the coarser level resC.
      m_op[ilev]->AMRRestrict(*(m_resC[ilev]),
                              *(m_residual[ilev]),
                              *(m_correction[ilev]),
                              *(m_correction[ilev-1]));

      // Overwrite residual on the valid region of the next coarser level
      //  with coarsened residual from this level
      m_op[ilev-1]->assign(*m_residual[ilev-1], *(m_resC[ilev]));

      //============finish Compute residual for the next coarser level======

      for (int img = 0; img < m_numMG; img++)
      {
        AMRVCycle(a_uberCorrection, a_uberResidual, ilev-1, l_max, l_base);
      }

      //================= Upsweep ======================
      //increment the correction with coarser version
      //m_op[ilev]->AMRProlong(*(m_correction[ilev]), *(m_correction[ilev-1]));
      m_op[ilev]->AMRProlong(*(m_correction[ilev]), *(m_correction[ilev-1]));
      //recompute residual
      m_op[ilev]->AMRUpdateResidual(*(m_residual[ilev]), *(m_correction[ilev]), *(m_correction[ilev-1]));

      //compute correction to the correction
      T& dCorr = *(a_uberCorrection[ilev]); // user uberCorrection as holder for correction to correction
      m_op[ilev]->setToZero(dCorr);
      this->relax(dCorr, *(m_residual[ilev]), ilev, m_post);

      //correct the correction with the correction to the correction
      m_op[ilev]->incr(*(m_correction[ilev]), dCorr, 1.0);

      m_op[ilev]->assign(*(a_uberCorrection[ilev]), *(m_correction[ilev]));
    }
  }


}
#endif
