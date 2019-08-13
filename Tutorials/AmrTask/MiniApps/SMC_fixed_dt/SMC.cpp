
#include <SMC.H>
#include <SMC_F.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

bool        SMC::initialized       = false;
int         SMC::ncons             = 0;
int         SMC::nprim             = 0;
int         SMC::nspec             = 0;
int         SMC::nplot             = 0;
int         SMC::ngrow             = 0;
Vector<int>  SMC::ncell             (3,128);
Vector<int>  SMC::max_grid_size     (3,64);
int         SMC::max_step          = 2;
Real        SMC::stop_time         = 3.e-3;
Vector<Real> SMC::prob_lo           (3,  0.0);
Vector<Real> SMC::prob_hi           (3,  0.1);
int         SMC::verbose           = 2;
int         SMC::cfl_int           = 10;
Real        SMC::cfl               = 0.1;
Real        SMC::init_shrink       = 0.5;
Real        SMC::fixed_dt          = 1.0e-4;
int         SMC::plot_int          = -1;

SMC::SMC ()
{
    if (!initialized) {

	init_runtime();

	init_stencil();

	init_chemistry();

	init_variables();

	initialized = true;
    }

    build_multifabs();

    init_from_scratch();

    wt_fb1 = wt_fb2 = wt_chem1 = wt_chem2 = wt_hypdiff = 0.0;
}

SMC::~SMC ()
{
    ;
}

    void
SMC::init_runtime ()
{
    ParmParse pp;
    int n;

    n = pp.countval("ncell");
    if (n == 1) {
	int nc;
	pp.query("ncell", nc);
	ncell[0] = ncell[1] = ncell[2] = nc;
    } else if (n == 3) {
	pp.queryarr("ncell", ncell, 0, 3);
    }

    n = pp.countval("max_grid_size");
    if (n == 1) {
	int mgs;
	pp.query("max_grid_size", mgs);
	max_grid_size[0] = max_grid_size[1] = max_grid_size[2] = mgs;
    } else if (n == 3) {
	pp.queryarr("max_grid_size", max_grid_size, 0, 3);
    }

    pp.query("max_step", max_step);
    pp.query("stop_time", stop_time);

    n = pp.countval("prob_lo");
    if (n == 1) {
	Real plo;
	pp.query("prob_lo", plo);
	prob_lo[0] = prob_lo[1] = prob_lo[2] = plo;
    } else {
	pp.queryarr("prob_lo", prob_lo, 0, 3);
    }

    n = pp.countval("prob_hi");
    if (n == 1) {
	Real phi;
	pp.query("prob_hi", phi);
	prob_hi[0] = prob_hi[1] = prob_hi[2] = phi;
    } else {
	pp.queryarr("prob_hi", prob_hi, 0, 3);
    }

    pp.query("verbose", verbose);

    pp.query("cfl_int", cfl_int);
    pp.query("cfl", cfl);
    pp.query("init_shrink", init_shrink);
    pp.query("fixed_dt", fixed_dt);

    pp.query("plot_int", plot_int);

}

    void
SMC::init_stencil ()
{
    derivative_stencil_init();
    ngrow = get_num_ghost_cells();
}

    void
SMC::init_chemistry ()
{
    chemistry_init();
    nspec = get_num_species();
}

    void
SMC::init_variables ()
{
    variables_init();
    ncons = get_num_cons();
    nprim = get_num_prim();
    nplot = 5 + nspec;
}


class MyAction :public Action{
    private:
	int stage, max_step;
	Real t, dt, stop_time;
    public:
	static SMC *thisSMC;
	MyAction(): stage(0), max_step(thisSMC->max_step), t(thisSMC->t), dt(thisSMC->fixed_dt), stop_time(thisSMC->stop_time){}
	void Compute(){
	    Box bx = validbox();
	    FArrayBox &U_fab= validFab(), &Uprime_fab= validFab(thisSMC->Uprime), &Utmp_fab= validFab(thisSMC->Utmp), &Q_fab= validFab(thisSMC->Q);
	    FArrayBox &mu_fab= validFab(thisSMC->mu), &xi_fab= validFab(thisSMC->xi), &lam_fab= validFab(thisSMC->lam), &Ddiag_fab= validFab(thisSMC->Ddiag);

	    switch(stage){
		case 0:
    		    Uprime_fab.setVal(0.0);
		    for(auto t: ta.tileArray) thisSMC->compute_dU0(U_fab, Uprime_fab, Q_fab, mu_fab, xi_fab, lam_fab, Ddiag_fab, bx, t, growntilebox(t, thisSMC->ngrow));
		    for(auto t: ta.tileArray) thisSMC->compute_dU1(U_fab, Uprime_fab, Q_fab, mu_fab, xi_fab, lam_fab, Ddiag_fab, t);
    		    Utmp_fab.setVal(0.0);
		    for(auto t: ta.tileArray) thisSMC->update_rk3(0.0, Utmp_fab, 1.0, U_fab, dt, Uprime_fab, stage, t);
		    for(auto t: ta.tileArray) thisSMC->reset_density(Utmp_fab, bx, t, thisSMC->Utmp.nGrow());
		    stage=1; break;

		case 1:
    		    Uprime_fab.setVal(0.0);
		    for(auto t: ta.tileArray) thisSMC->compute_dU0(Utmp_fab, Uprime_fab, Q_fab, mu_fab, xi_fab, lam_fab, Ddiag_fab, bx, t, growntilebox(t, thisSMC->ngrow));
		    for(auto t: ta.tileArray) thisSMC->compute_dU1(Utmp_fab, Uprime_fab, Q_fab, mu_fab, xi_fab, lam_fab, Ddiag_fab, t);
		    for(auto t: ta.tileArray) thisSMC->update_rk3(0.25, Utmp_fab, 0.75, U_fab, 0.25*dt, Uprime_fab, stage, t);
		    for(auto t: ta.tileArray) thisSMC->reset_density(Utmp_fab, bx, t, thisSMC->Utmp.nGrow());
		    stage=2; break;

		case 2:
    		    Uprime_fab.setVal(0.0);
		    for(auto t: ta.tileArray) thisSMC->compute_dU0(Utmp_fab, Uprime_fab, Q_fab, mu_fab, xi_fab, lam_fab, Ddiag_fab, bx, t, growntilebox(t, thisSMC->ngrow));
		    for(auto t: ta.tileArray) thisSMC->compute_dU1(Utmp_fab, Uprime_fab, Q_fab, mu_fab, xi_fab, lam_fab, Ddiag_fab, t);
		    for(auto t: ta.tileArray) thisSMC->update_rk3(1./3., U_fab, 2./3., Utmp_fab, (2./3.)*dt, Uprime_fab, stage, t);
		    for(auto t: ta.tileArray) thisSMC->reset_density(U_fab, bx, t, thisSMC->U.nGrow());
	            t+= dt;
	            if(_iter<3*max_step && t<stop_time){
		        extendIters(3);
		        stage=0;
			cout <<"completed step " << _iter << " time "<< t<< " dt "<< dt <<endl;
                    }//else this task completed
	    }
	}
};


SMC* MyAction::thisSMC;

    void
SMC::evolve ()
{
    int istep = 0;
    int last_plt_written = -1;

    if (plot_int > 0) {
	writePlotFile(istep);
	last_plt_written = istep;
    }

    int init_step = 1;

    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "\n" << "BEGIN MAIN EVOLUTION LOOP" << "\n" << std::endl;
    }

    ParallelDescriptor::Barrier();
    Real wt0 = ParallelDescriptor::second();


    /*
       if ( (max_step >= init_step) && (t < stop_time || stop_time < 0.0) ) 
       {
       for (istep = init_step; istep <= max_step; ++istep) 
       {
       if (ParallelDescriptor::IOProcessor()) {
       std::cout << "Advancing time step " << istep << " time = " << t << "\n";
       }



       if (ParallelDescriptor::IOProcessor()) {
       std::cout << "End of step " << istep << " time = " << t << "\n" << std::endl;
       }

       if (plot_int > 0 && istep%plot_int == 0) {
       writePlotFile(istep);
       last_plt_written = istep;
       }

       if (stop_time >= 0.0 && t >= stop_time) {
       break;
       }
       }

       if (istep > max_step) istep = max_step;

       if (plot_int > 0 && last_plt_written != istep) {
       writePlotFile(istep);
       }
       }
     */
    MyAction::thisSMC= this;
    AMFIter<MyAction> mfi(U, 3, geom.periodicity(), true /*do_tiling*/);
    mfi.Iterate();

#ifdef BL_USE_UPCXX
    upcxx::barrier();
#else
    ParallelDescriptor::Barrier();
#endif

    Real wt1 = ParallelDescriptor::second();

    Real wt_fb = wt_fb1 + wt_fb2;
    Real wt_chem = wt_chem1 + wt_chem2;

    ParallelDescriptor::ReduceRealMax(wt_fb     , ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealMax(wt_chem   , ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealMax(wt_hypdiff, ParallelDescriptor::IOProcessorNumber());
    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "======================================" << std::endl;
	std::cout << " Total Time            : " << wt1-wt0 << "\n";
	std::cout << "     Communication time: " << wt_fb << "\n";
	std::cout << "     Chemistry     time: " << wt_chem << "\n";
	std::cout << "     Hyp-Diff      time: " << wt_hypdiff << "\n";
	std::cout << "======================================" << std::endl;
    }
}
