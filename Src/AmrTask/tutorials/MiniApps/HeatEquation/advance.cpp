#include "AMReX_AbstractTask.H"
#include "AMReX_TaskGraph.H"
#include "RTS.H" 
#include "AMReX_AsyncMFIter.H"
#include "myfunc.H"
#include "myfunc_F.H"

using namespace amrex;

class MyAction :public Action{
    public:
	static MultiFab *old_phi, *new_phi;
	static std::array<MultiFab, AMREX_SPACEDIM>* flux;
	static Real dt, *dx;
	void Compute(){
#if 0
	    Box& bx = validbox();
	    int idx = localFabIdx();
	    compute_flux(BL_TO_FORTRAN_BOX(bx),
		    BL_TO_FORTRAN_ANYD((*old_phi)[idx]),
		    BL_TO_FORTRAN_ANYD((*flux)[0][idx]),
		    BL_TO_FORTRAN_ANYD((*flux)[1][idx]),
		    BL_TO_FORTRAN_ANYD((*flux)[2][idx]),
		    dx);
	    update_phi(BL_TO_FORTRAN_BOX(bx),
		    BL_TO_FORTRAN_ANYD((*old_phi)[idx]),
		    BL_TO_FORTRAN_ANYD((*new_phi)[idx]),
		    BL_TO_FORTRAN_ANYD((*flux)[0][idx]),
		    BL_TO_FORTRAN_ANYD((*flux)[1][idx]),
		    BL_TO_FORTRAN_ANYD((*flux)[2][idx]),
		    dx, dt);
#endif
	}
};

MultiFab * MyAction::old_phi;
MultiFab * MyAction::new_phi;
std::array<MultiFab, AMREX_SPACEDIM>* MyAction::flux;
Real MyAction::dt;
Real* MyAction::dx;

void advance (MultiFab& old_phi, MultiFab& new_phi,
	std::array<MultiFab, AMREX_SPACEDIM>& flux,
	Real dt, const Geometry& geom, int nsteps)
{
    // Fill the ghost cells of each grid from the other grids
    // includes periodic domain boundaries
    //old_phi.FillBoundary(geom.periodicity());

    // Fill non-periodic physical boundaries
    //fill_physbc(old_phi, geom);

    Real* dx = (Real*)geom.CellSize();
    MyAction::old_phi= &old_phi;
    MyAction::new_phi= &new_phi;
    MyAction::flux= &flux;
    MyAction::dt= dt;
    MyAction::dx= dx;

    //
    // Note that this simple example is not optimized.
    // The following two MFIter loops could be merged
    // and we do not have to use flux MultiFab.
    // 

    AMFIter<MyAction> mfi(old_phi, nsteps);
    mfi.Iterate();


}
