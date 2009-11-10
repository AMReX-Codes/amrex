#include "DataServices.H"

void
ComputeAmrDataMeanVar (AmrData&     amrData,
		       Array<Real>& mean,
		       Array<Real>& variance,
		       int          sComp,
		       int          nComp,
		       bool         verbose);

void
ComputeAmrDataMeanVar (AmrData&           amrData,
		       Array<std::string> cNames,
		       Array<BoxArray>    bas,
		       Array<Real>&       mean,
		       Array<Real>&       variance);

void
ComputeAmrDataList (AmrData&         amrData,
		    Array<MultiFab*> mean,
		    Array<MultiFab*> variance,
		    int              sComp,
		    int              nComp);

void
ComputeAmrDataPDF (AmrData&     amrData,
		   Real**       icount,
		   int          nBin,
		   Array<std::string> cNames);

void
ComputeAmrDataPDF (AmrData&     amrData,
		   Real**       icount,
		   int          nBin,
		   Array<std::string> cNames,
		   Array<BoxArray> bas);

void
ComputeAmrDataVAR (AmrData&           amrData,
		   int                nBin,
		   Array<std::string> cNames,
		   Array<Real>        barr,
		   std::string        oFile);

void
VariogramUniform (AmrData&             amrData,
		  Array<std::string>   cNames,
		  Array<Real>          barr,
		  Array< Array<int> >  ivoption,
		  int                  nlag,
		  int                  isill,
		  Array<Real>          sills,
		  std::string          oFile);
