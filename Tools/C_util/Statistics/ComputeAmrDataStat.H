#include "DataServices.H"

void
ComputeAmrDataMeanVar (AmrData&     amrData,
		       Array<Real>& mean,
		       Array<Real>& variance,
		       int          sComp,
		       int          nComp,
		       bool         verbose);

void
ComputeAmrDataList (AmrData&         amrData,
		    Array<MultiFab*> mean,
		    Array<MultiFab*> variance,
		    int              sComp,
		    int              nComp);
