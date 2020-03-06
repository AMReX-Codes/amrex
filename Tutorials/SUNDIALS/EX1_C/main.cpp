#include <fstream>
#include <iomanip>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts. */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver     */
#include <cvode/cvode_spils.h>         /* access to CVSpils interface */
#include <sundials/sundials_types.h>   /* definition of type realtype */
#include <sundials/sundials_math.h>    /* definition of ABS and EXP   */
#include <nvector/nvector_serial.h>
#include "myfunc_F.H"

/* Functions Called by the Solver */
static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

static void PrintFinalStats(void *cvode_mem);

/* Private function to check function return values */
static int check_flag(void *flagvalue, const char *funcname, int opt);

using namespace amrex;

/*
#pragma gpu
void device_ptr_wrapper(N_Vector u, Real* dptr)
      dptr=N_VGetDeviceArrayPointer_Cuda(u);
*/

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = ParallelDescriptor::second();

    std::cout << std::setprecision(15);

    bool test_ic;
    int n_cell, max_grid_size;
    int cvode_meth, cvode_itmeth, write_plotfile;
    bool do_tiling;

  realtype reltol, abstol, t, tout, umax;
  N_Vector u;
  SUNLinearSolver LS;
  void *cvode_mem;
  int iout, flag;
  long int nst;

  test_ic=false;
  u = NULL;
  LS = NULL;
  cvode_mem = NULL;

  reltol = 1e-6;  /* Set the tolerances */
  abstol = 1e-10;

    // inputs parameters
    {
      // ParmParse is way of reading inputs from the inputs file
      ParmParse pp;

      // We need to get n_cell from the inputs file - this is the number of
      // cells on each side of a square (or cubic) domain.
      pp.get("n_cell",n_cell);

      // Default nsteps to 0, allow us to set it to something else in the
      // inputs file
      pp.get("max_grid_size",max_grid_size);

      // Select CVODE solve method.
      //   1 = Adams (for non-stiff problems)
      //   2 = BDF (for stiff problems)
      pp.get("cvode_meth",cvode_meth);
      // Select CVODE solver iteration method.
      //   1 = Functional iteration
      //   2 = Newton iteration
      pp.get("cvode_itmeth",cvode_itmeth);

      pp.get("write_plotfile",write_plotfile);
      pp.get("do_tiling",do_tiling);
    }

    if (cvode_meth < 1)
      amrex::Abort("Unknown cvode_meth");
    if (cvode_itmeth < 1)
      amrex::Abort("Unknown cvode_itmeth");

    amrex::Print() << "This is AMReX version " << amrex::Version() << std::endl;
    amrex::Print() << "Problem domain size: nx = ny = nz = " << n_cell << std::endl;
    amrex::Print() << "Max grid size: " << max_grid_size << std::endl;
    amrex::Print() << "CVODE method: ";
    if (cvode_meth == 1) {
      amrex::Print() << "Adams (non-stiff)";
    } else if (cvode_meth == 2) {
        amrex::Print() << "BDF (stiff)";
    }
    amrex::Print() << std::endl;
    amrex::Print() << "CVODE iteration method: ";
    if (cvode_itmeth == 1) {
      amrex::Print() << "Functional";
    } else if (cvode_itmeth == 2) {
        amrex::Print() << "Newton";
    }
    amrex::Print() << std::endl;

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
      IntVect dom_lo(IntVect(D_DECL(0,0,0)));
      IntVect dom_hi(IntVect(D_DECL(n_cell-1, n_cell-1, n_cell-1)));
      Box domain(dom_lo, dom_hi);

      // Initialize the boxarray "ba" from the single box "bx"
      ba.define(domain);

      // Break up boxarray "ba" into chunks no larger than "max_grid_size"
      // along a direction
      ba.maxSize(max_grid_size);

      // This defines the physical size of the box.  Right now the box is
      // [-1,1] in each direction.
      RealBox real_box;
      for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n,-1.0);
        real_box.setHi(n, 1.0);
      }

      // This sets the boundary conditions to be doubly or triply periodic
      int is_periodic[BL_SPACEDIM];
      for (int i = 0; i < BL_SPACEDIM; i++) {
        is_periodic[i] = 1;
      }

      // This defines a Geometry object
      geom.define(domain,&real_box,CoordSys::cartesian,is_periodic);
    }

    // Ncomp = number of components for each array
    int Ncomp  = 1;

    // time = starting time in the simulation
    Real time = 0.0;

    DistributionMapping dm(ba);

    // Create MultiFab with no ghost cells.
    MultiFab mf(ba, dm, Ncomp, 0);

    
    mf.setVal(0.0);
    

#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(mf, do_tiling); mfi.isValid(); ++mfi )
    {
      t=0;
      tout=2;
      Real* dptr;
      Real* dptr_compare;

      const Box& tbx = mfi.tilebox();
      amrex::IntVect tile_size = tbx.size();
      const int* hi = tbx.hiVect();
      const int* lo = tbx.loVect();
      int long neq1=(hi[0]-lo[0]+1)*(hi[1]-lo[1]+1)*(hi[2]-lo[2]+1);
      int long neq=(tile_size[0])*(tile_size[1])*(tile_size[2]);

      if(neq>1)
	amrex::Print()<<"Integrating a box with "<<neq<<" cels"<<std::endl;

      /* Create a CUDA vector with initial values */
      u = N_VNew_Serial(neq);  /* Allocate u vector */
      N_Vector ucomp = N_VNew_Serial(neq);  /* Allocate u vector */
      if(check_flag((void*)u, "N_VNew_Serial", 0)) return(1);

      FSetIC_mfab(mf[mfi].dataPtr(),
        tbx.loVect(),
	    tbx.hiVect());  /* Initialize u vector */

      dptr=N_VGetArrayPointer_Serial(u);
      if(test_ic==true)
      dptr_compare=N_VGetArrayPointer_Serial(ucomp);
      mf[mfi].copyToMem(tbx,0,1,dptr);

      if(test_ic==true)      
	{
      FSetIC(mf[mfi].dataPtr(),
        tbx.loVect(),
	tbx.hiVect(),dptr_compare); 

      N_Vector z_result=N_VNew_Serial(neq);
      N_Vector weight=N_VNew_Serial(neq);
      N_VConst(1.0,weight);
      N_VLinearSum(1,u,-1,ucomp,z_result);
      double result=N_VWrmsNorm(z_result,weight);
      amrex::Print()<<"Difference in initial condition is: "<<result<<std::endl;
	}
      /* Call CVodeCreate to create the solver memory and specify the 
       * Backward Differentiation Formula and the use of a Newton iteration */
      #ifdef CV_NEWTON
      cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
      #else
      cvode_mem = CVodeCreate(CV_BDF);
      #endif
      if(check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

      /* Call CVodeInit to initialize the integrator memory and specify the
       * user's right hand side function in u'=f(t,u), the initial time T0, and
       * the initial dependent variable vector u. */
      flag = CVodeInit(cvode_mem, f, t, u);
      if(check_flag(&flag, "CVodeInit", 1)) return(1);

      /* Call CVodeSStolerances to specify the scalar relative tolerance
       * and scalar absolute tolerance */
      flag = CVodeSStolerances(cvode_mem, reltol, abstol);
      if (check_flag(&flag, "CVodeSStolerances", 1)) return(1);

      /* Create SPGMR solver structure without preconditioning
       * and the maximum Krylov dimension maxl */
      LS = SUNSPGMR(u, PREC_NONE, 0);
      if(check_flag(&flag, "SUNSPGMR", 1)) return(1);

      /* Set CVSpils linear solver to LS */
      flag = CVSpilsSetLinearSolver(cvode_mem, LS);
      if(check_flag(&flag, "CVSpilsSetLinearSolver", 1)) return(1);

      /* Call CVode */
      flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
      if(check_flag(&flag, "CVode", 1)) break;

      mf[mfi].copyFromMem(tbx,0,1,dptr);

      N_VDestroy(u);          /* Free the u vector */
      CVodeFree(&cvode_mem);  /* Free the integrator memory */
    
    }

    amrex::Print()<<"Maximum of repacked final solution: "<<mf.max(0,0,0)<<std::endl;
    
    if (write_plotfile)
    {
      amrex::WriteSingleLevelPlotfile("PLT_OUTPUT",
                                      mf,
                                      {"y1"},
                                      geom,
                                      time,
                                      0);
    }

    // Call the timer again and compute the maximum difference between the start time and stop time
    //   over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;
    
    amrex::Finalize();
    return 0;
}

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  //  RhsFn(t,u,udot,user_data);
  /*  int b=5;
      test<<<1,1>>>(b);*/
  N_VConst_Serial(2.0*t,udot);
  return 0;
}

static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */

  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if flag < 0 */

  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */

  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}
