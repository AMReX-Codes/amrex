#include <DataServices.H>

void 
compute_flux_all(int  nstart,
		 int  nmax,
		 int  nfac,
		 int  dir,
		 std::string iFile,
		 Real phi);


void 
compute_flux_all(int nstart,
		 int nmax,
		 int nfac,
		 int dir,
		 std::string iFile,
		 MultiFab& phidata);
void
compute_flux(AmrData&           amrData, 
	     int                dir, 
	     Array<std::string> cNames,
	     Real               dt,
	     Array<Real>&       xold,
	     Real&              flux,
	     Real               phi,
	     bool               do_init,
 	     Array<Real>&       barr);
