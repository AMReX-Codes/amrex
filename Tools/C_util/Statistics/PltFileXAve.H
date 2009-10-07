#include "DataServices.H"

void 
compute_flux_all(int nstart,
		 int nmax,
		 int nfac,
		 std::string iFile);

void
compute_flux(AmrData&           amrData, 
	     int                dir, 
	     Array<std::string> cNames,
	     Real               dt,
	     Real               flux,
	     Real               phi,
 	     Real*              barr=0);
