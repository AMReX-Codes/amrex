//
// $Id: MultiGrid.cpp,v 1.1 1998-03-23 23:36:17 almgren Exp $
// 

#ifdef BL_USE_NEW_HFILES
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#include <ParmParse.H>
#include <Misc.H>
#include <Utility.H>
#include <ParallelDescriptor.H>
#include <CGSolver.H>
#include <MG_F.H>
#include <MultiGrid.H>

const char NL = '\n';

bool MultiGrid::initialized = false;
int MultiGrid::def_nu_0 = 1;
int MultiGrid::def_nu_1 = 2;
int MultiGrid::def_nu_2 = 2;
int MultiGrid::def_nu_f = 8;
int MultiGrid::def_maxiter = 40;
int MultiGrid::def_numiter = -1;
int MultiGrid::def_verbose = 0;
int MultiGrid::def_usecg = 1;
Real MultiGrid::def_rtol_b = 0.01;
Real MultiGrid::def_atol_b = -1.0;
int MultiGrid::def_nu_b = 0;
int MultiGrid::def_numLevelsMAX = 1024;

void
MultiGrid::initialize()
{
    ParmParse pp("mg");
    initialized = true;
    pp.query("maxiter", def_maxiter);
    pp.query("numiter", def_numiter);
    pp.query("nu_0", def_nu_0);
    pp.query("nu_1", def_nu_1);
    pp.query("nu_2", def_nu_2);
    pp.query("nu_f", def_nu_f);
    pp.query("verbose", def_verbose);
    pp.query("v", def_verbose);
    pp.query("usecg", def_usecg);
    pp.query("rtol_b", def_rtol_b);
    pp.query("bot_atol", def_atol_b);
    pp.query("nu_b", def_nu_b);
    pp.query("numLevelsMAX", def_numLevelsMAX);
    if(def_verbose) {
      if(ParallelDescriptor::IOProcessor()) {
        cout << "def_nu_0 = " << def_nu_0 << NL;
        cout << "def_nu_1 = " << def_nu_1 << NL;
        cout << "def_nu_2 = " << def_nu_2 << NL;
        cout << "def_nu_f = " << def_nu_f << NL;
        cout << "def_maxiter = " << def_maxiter << NL;
        cout << "def_usecg = "  << def_usecg << NL;
        cout << "def_rtol_b = " << def_rtol_b << NL;
        cout << "def_atol_b = " << def_atol_b << NL;
        cout << "def_nu_b = "   << def_nu_b << NL;
        cout << "def_numLevelsMAX = "   << def_numLevelsMAX << NL;
      }
    }
}

MultiGrid::MultiGrid(LinOp &_Lp)
    : Lp(_Lp), initialsolution((MultiFab*)NULL)
{
    if(!initialized)
        initialize();

    maxiter = def_maxiter;
    numiter = def_numiter;
    nu_0 = def_nu_0;
    nu_1 = def_nu_1;
    nu_2 = def_nu_2;
    nu_f = def_nu_f;
    usecg = def_usecg;
    verbose = def_verbose;
    rtol_b = def_rtol_b;
    atol_b = def_atol_b;
    nu_b = def_nu_b;
    numLevelsMAX = def_numLevelsMAX;
    numlevels = numLevels();
}

MultiGrid::~MultiGrid()
{
    delete initialsolution;
    int i;
    for(i = 0; i < cor.length(); ++i) {
        delete res[i];
        delete rhs[i];
        delete cor[i];
    }
}

Real
MultiGrid::errorEstimate(int level, LinOp::BC_Mode bc_mode)
{
      // Get inf-norm of residual
    int p = 0;
    Lp.residual(*res[level], *rhs[level], *cor[level], level, bc_mode);
    Real restot = 0.0;
    Real resk  = 0.0;
    const BoxArray &gbox = Lp.boxArray(0);
    //int gn;
    //for(gn = 0; gn < res[level]->length(); ++gn) {
    for(MultiFabIterator resmfi(*res[level]); resmfi.isValid(); ++resmfi) {
        assert(gbox[resmfi.index()] == resmfi.validbox());
        Box tmp(resmfi.validbox());
        resk = (resmfi().norm(tmp, p));
        restot = Max(restot, resk);
    }
    ParallelDescriptor::ReduceRealMax(restot);
    // Warning:  this reduction will result in
    //           restot = Max(restot<all procs with grids>, 0.0)
    //           if there is a processor with no grids
    //           because restot is initialized to 0.0
    return restot;
}

void
MultiGrid::prepareForLevel(int level)
{
      // Build this level by allocating reqd internal MultiFabs if necessary
    if(cor.length() > level)
        return;
    res.resize(level+1, (MultiFab*)0);
    rhs.resize(level+1, (MultiFab*)0);
    cor.resize(level+1, (MultiFab*)0);
    Lp.prepareForLevel(level);
    if(cor[level] == 0) {
        res[level] = new MultiFab(Lp.boxArray(level), 1, 1, Fab_allocate);
        rhs[level] = new MultiFab(Lp.boxArray(level), 1, 1, Fab_allocate);
        cor[level] = new MultiFab(Lp.boxArray(level), 1, 1, Fab_allocate);
        if(level == 0) {
            initialsolution = new MultiFab(Lp.boxArray(0), 1, 1, Fab_allocate);
            
        }
    }
}

void
MultiGrid::residualCorrectionForm(MultiFab& resL, const MultiFab& rhsL,
                                  MultiFab& solnL, const MultiFab& inisol,
                                  LinOp::BC_Mode bc_mode, int level)
{
      // Using the linearity of the operator, Lp, we can solve this system
      // instead by solving for the correction required to the initial guess.
    initialsolution->copy(inisol);
    solnL.copy(inisol);
    Lp.residual(resL, rhsL, solnL, level, bc_mode);
    solnL.setVal(0.0);
}

void
MultiGrid::solve(MultiFab &_sol, const MultiFab &_rhs,
                 Real _eps_rel, Real _eps_abs,
                 LinOp::BC_Mode bc_mode)
{
      // Prepare memory for new level, and solve the general boundary
      // value problem to within relative error _eps_rel.  Customized
      // to solve at level=0
    int level = 0;
    prepareForLevel(level);
    residualCorrectionForm(*rhs[level], _rhs, *cor[level], _sol,
                           bc_mode, level);
    if( ! solve_(_sol, _eps_rel, _eps_abs, LinOp::Homogeneous_BC, level) ) {
        BoxLib::Error("MultiGrid:: failed to converge!");
    }
}

int
MultiGrid::solve_(MultiFab &_sol, Real eps_rel, Real eps_abs,
                  LinOp::BC_Mode bc_mode, int level)
{
      // Relax system maxiter times, stop if relative error <= _eps_rel or
      // if absolute err <= _abs_eps
    long returnVal = 0;  // should use bool for return value from this function
    Real error0 = errorEstimate(level, bc_mode);
    Real error  = error0;
    if( verbose ) {
      if(ParallelDescriptor::IOProcessor()) {
        for (int k=0; k<level; k++) cout << "   ";
        cout << "MultiGrid: Initial error (error0) = " << error0 << NL;
      }
    }

    if( eps_rel < 1.0e-16 && eps_rel > 0.0) {
      if(ParallelDescriptor::IOProcessor()) {
        cerr << "MultiGrid: Tolerance " << eps_rel
             << " < 1e-16 is probably set too low" << NL;
      }
    }
    int nit;

      // Initialize correction to zero at this level (auto-filled at levels below)
    (*cor[level]).setVal(0.0);
    
      //  Note: if eps_rel, eps_abs < 0 then that test is effectively bypassed
    for(nit = 0;
        (nit < maxiter)
            &&   (nit < numiter || numiter < 0)
            &&   (error > eps_rel*error0)
            &&   (error > eps_abs);
        ++nit) {

        relax(*cor[level], *rhs[level], level, eps_rel, eps_abs, bc_mode);
        error = errorEstimate(level, bc_mode);
        if(  verbose > 1 ||
           (((eps_rel > 0. && error < eps_rel*error0) ||
             (eps_abs > 0. && error < eps_abs)) && verbose) ) {
            for (int k=0; k<level; k++) cout << "   ";
            if(ParallelDescriptor::IOProcessor()) {
              cout << "MultiGrid: Iteration " << nit << " error/error0 "
                   << error/error0 << NL;
            }
        }
    }

    if ( nit==numiter
         || error <= eps_rel*error0
         || error <= eps_abs ) {
        
      // Omit ghost update since maybe not initialized in calling routine.
      // BoxLib_1.99 has no MultiFab::plus(MultiFab&) member, which would
      // operate only in valid regions; do explicitly.  Add to boundary
      // values stored in initialsolution.
        _sol.copy(*cor[level]);
        _sol.plus(*initialsolution,0,_sol.nComp(),0);
        returnVal = 1;
    } else {
      returnVal = 0;
    }

    // use long instead of bool--there is a bug on the t3e which
    // causes reduction problems with bool
    ParallelDescriptor::ReduceLongAnd(returnVal);

      // Otherwise, failed to solve satisfactorily
    return ((int) returnVal);
}

int
MultiGrid::getNumLevels(int _numlevels)
{
    assert(_numlevels >= 0);
    int oldnumlevels = numlevels;
    numlevels = Min(_numlevels, numLevels());
    return oldnumlevels;
}

int
MultiGrid::numLevels() const
{
    int ng = Lp.numGrids();
    int lv = numLevelsMAX;

      // the routine `falls through' since coarsening and refining
      // a unit box does not yield the initial box.
    const BoxArray &bs = Lp.boxArray(0);
    int i;
    for(i = 0; i < ng; ++i) {
        int llv = 0;
        Box tmp = bs[i];
        for(;;) {
            Box ctmp = tmp; ctmp.coarsen(2);
            Box rctmp = ctmp; rctmp.refine(2);
            if (tmp != rctmp) {
                break;
            }
            if (ctmp.numPts() == 1)
                break;
            llv++;
            tmp = ctmp;
        }
          // Set number of levels so that every box can be refined to there
        lv = (lv < llv)? lv : llv;
    }
    return lv+1;                // including coarsest
}

void
MultiGrid::relax(MultiFab& solL, MultiFab& rhsL, int level,
                 Real eps_rel, Real eps_abs, LinOp::BC_Mode bc_mode)
{
      // Recursively relax system.  Equivalent to multigrid V-cycle.
      // At coarsest grid, call coarsestSmooth
    int i;
    if(level < numlevels - 1 ) {
        for(i = preSmooth() ; i > 0 ; i--) {
            Lp.smooth(solL, rhsL, level, bc_mode);
        }
        Lp.residual(*res[level], rhsL, solL, level, bc_mode);
        prepareForLevel(level+1);
        average(*rhs[level+1], *res[level]);
        cor[level+1]->setVal(0.0);
        for(i = cntRelax(); i > 0 ; i--) {
            relax(*cor[level+1], *rhs[level+1], level+1,
                  eps_rel, eps_abs, bc_mode);
        }
        interpolate(solL, *cor[level+1]);
        for(i = postSmooth(); i > 0 ; i--) {
            Lp.smooth(solL, rhsL, level, bc_mode);
        }
    } else {
        coarsestSmooth(solL, rhsL, level, eps_rel, eps_abs, bc_mode);
    }
}

void
MultiGrid::coarsestSmooth(MultiFab& solL, MultiFab& rhsL, int level,
                          Real eps_rel, Real eps_abs, LinOp::BC_Mode bc_mode)
{
    prepareForLevel(level);
    if (usecg == 0) {
        for(int i = finalSmooth(); i > 0; i--) {
            Lp.smooth(solL, rhsL, level, bc_mode);
        }
    } else {
        bool use_mg_precond = false;
        CGSolver cg(Lp, use_mg_precond, level);
        cg.solve(solL, rhsL, rtol_b, atol_b, bc_mode);
        for(int i=0; i<nu_b; i++) {
            Lp.smooth(solL, rhsL, level, bc_mode);
        }
    }
}

void
MultiGrid::average(MultiFab &c, const MultiFab &f)
{
      // Use Fortran function to average down (restrict) f to c
    //int gn;
    //for(gn = 0; gn < f.length(); ++gn) {
    for(MultiFabIterator cmfi(c); cmfi.isValid(); ++cmfi) {
      DependentMultiFabIterator fmfi(cmfi, f);
        assert(c.boxArray().get(cmfi.index()) == cmfi.validbox());
        Box bx(cmfi.validbox());
        int nc = c.nComp();
        FORT_AVERAGE(
            cmfi().dataPtr(), ARLIM(cmfi().loVect()), ARLIM(cmfi().hiVect()),
            fmfi().dataPtr(), ARLIM(fmfi().loVect()), ARLIM(fmfi().hiVect()),
            bx.loVect(), bx.hiVect(),
            &nc);
    }
}

void
MultiGrid::interpolate(MultiFab &f, const MultiFab &c)
{
      // Use fortran function to interpolate up (prolong) c to f
      // Note: returns f=f+P(c) , i.e. ADDS interp'd c to f
    //int gn;
    //for(gn = 0; gn < f.length(); ++gn) {
    for(MultiFabIterator fmfi(f); fmfi.isValid(); ++fmfi) {
      DependentMultiFabIterator cmfi(fmfi, c);
        assert(c.boxArray().get(cmfi.index()) == cmfi.validbox());
        Box bx(cmfi.validbox());
        int nc = f.nComp();
        FORT_INTERP(
            fmfi().dataPtr(), ARLIM(fmfi().loVect()), ARLIM(fmfi().hiVect()),
            cmfi().dataPtr(), ARLIM(cmfi().loVect()), ARLIM(cmfi().hiVect()),
            bx.loVect(), bx.hiVect(),
            &nc);
    }
}
