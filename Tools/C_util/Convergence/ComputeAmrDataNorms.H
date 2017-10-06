#include <AMReX_DataServices.H>

void
ComputeAmrDataNorms (amrex::AmrData&     amrData,
		     amrex::Array<amrex::Real>& norm0,
		     amrex::Array<amrex::Real>& norm1,
		     amrex::Array<amrex::Real>& norm2,
		     bool         verbose);

void
ComputeAmrDataInt (amrex::AmrData&     amrData,
		   amrex::Array<amrex::Real>& norm1,
		   bool         verbose);
