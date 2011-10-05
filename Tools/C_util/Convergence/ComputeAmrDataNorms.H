#include <DataServices.H>

void
ComputeAmrDataNorms (AmrData&     amrData,
		     Array<Real>& norm0,
		     Array<Real>& norm1,
		     Array<Real>& norm2,
		     bool         verbose);

void
ComputeAmrDataInt (AmrData&     amrData,
		   Array<Real>& norm1,
		   bool         verbose);
