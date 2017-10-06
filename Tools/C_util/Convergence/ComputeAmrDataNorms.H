#include <AMReX_DataServices.H>

void
ComputeAmrDataNorms (AmrData&     amrData,
		     Vector<Real>& norm0,
		     Vector<Real>& norm1,
		     Vector<Real>& norm2,
		     bool         verbose);

void
ComputeAmrDataInt (AmrData&     amrData,
		   Vector<Real>& norm1,
		   bool         verbose);
