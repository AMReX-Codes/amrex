#include <AMReX_DataServices.H>
#include <AMReX_AmrData.H>

#ifndef _computeamredatanorms__
#define _computeamredatanorms__
void
ComputeAmrDataNorms (amrex::AmrData&     amrData,
		     amrex::Vector<amrex::Real>& norm0,
		     amrex::Vector<amrex::Real>& norm1,
		     amrex::Vector<amrex::Real>& norm2,
		     bool         verbose);


void
ComputeAmrDataInt (amrex::AmrData&     amrData,
                     amrex::Vector<amrex::Real>& norm1,
		   bool         verbose);
#endif
