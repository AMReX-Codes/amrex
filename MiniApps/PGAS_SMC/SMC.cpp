
#include <ParmParse.H>

#include <SMC.H>
#include <SMC_F.H>

bool        SMC::initialized       = false;
int         SMC::ncons             = 0;
int         SMC::nprim             = 0;
int         SMC::nspec             = 0;
int         SMC::nplot             = 0;
int         SMC::ngrow             = 0;
Array<int>  SMC::ncell             (3,128);
int         SMC::max_grid_size     = 64;
int         SMC::max_step          = 5;
Real        SMC::stop_time         = 3.e-3;
Array<Real> SMC::prob_lo           (3, -0.05);
Array<Real> SMC::prob_hi           (3,  0.05);
int         SMC::verbose           = 2;
int         SMC::cfl_int           = 10;
Real        SMC::cfl               = 0.1;
Real        SMC::init_shrink       = 0.5;
Real        SMC::fixed_dt          = -1.0e10;
int         SMC::plot_int          = -1;
int         SMC::overlap_comm_comp = 1;

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

    pp.query("max_grid_size", max_grid_size);
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

    pp.query("overlap_comm_comp", overlap_comm_comp);
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

    if ( (max_step >= init_step) && (t < stop_time || stop_time < 0.0) ) 
    {
	for (istep = init_step; istep <= max_step; ++istep) 
	{
	    if (ParallelDescriptor::IOProcessor()) {
		std::cout << "Advancing time step " << istep << " time = " << t << "\n" << std::endl;
	    }

	    advance(istep);

	    t = t + dt;

	    if (ParallelDescriptor::IOProcessor()) {
		std::cout << "End of step " << istep << " time = " << t << std::endl;
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

    // print timing results

}
