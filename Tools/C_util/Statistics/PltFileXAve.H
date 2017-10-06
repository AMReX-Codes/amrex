#include <AMReX_DataServices.H>

void 
compute_flux_all(int  nstart,
		 int  nmax,
		 int  nfac,
		 std::string iFile,
		 Real phi);


void 
compute_flux_all(int nstart,
		 int nmax,
		 int nfac,
		 std::string iFile,
		 MultiFab& phidata);
void
compute_flux(AmrData&           amrData, 
	     int                dir, 
	     Vector<std::string> cNames,
	     Real               dt,
	     Vector<Real>&       xold,
	     Real&              flux,
	     Real               phi,
	     bool               do_init,
 	     Real*              barr=0);
