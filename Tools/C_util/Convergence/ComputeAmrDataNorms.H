#include <AMReX_DataServices.H>

void
ComputeAmrDataNorms (amrex::AmrData&     amrData,
		     amrex::Vector<amrex::Real>& norm0,
		     amrex::Vector<amrex::Real>& norm1,
		     amrex::Vector<amrex::Real>& norm2,
		     bool         verbose);


void
ComputeAmrDataInt (AmrData&     amrData,
		   Vector<Real>& norm1,
		   bool         verbose);
