#include <DataServices.H>

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
ComputeMeanVarMF (MultiFab&          mf,
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

void
VariogramCross(AmrData&             amrData,
	       Array<std::string>   cNames,
	       MultiFab&            mf,
	       Array<Real>          barr,
	       Array< Array<int> >  ivoption,
	       int                  nlag,
	       int                  isill,
	       Array<Real>          sills,
	       std::string          oFile);

void
VariogramUniformMF (const MultiFab&      mf,
		    Array<Real>          dx,
		    Array< Array<int> >  ivoption,
		    int                  nlag,
		    int                  isill,
		    Array<Real>          sills,
		    std::string          oFile);
void
VariogramUniformMFG (const MultiFab&      mf,
		     Array<Real>          dx,
		     IntVect              sm,
		     IntVect              bg,
		     Array< Array<int> >  ivoption,
		     int                  nlag,
		     int                  isill,
		     Array<Real>          sills,
		     std::string          oFile);
void
TakeDifferenceFine(AmrData&             amrDataf,
		   AmrData&             amrDatac,
		   Array<std::string>   cNames,
		   Array<Real>          barr,
		   std::string          oFile);

void
TakeDifferenceCrse(AmrData&             amrDataf,
		   AmrData&             amrDatac,
		   Array<std::string>   cNames,
		   Array<Real>          barr,
		   std::string          oFile);
void
TakeDifferenceSum(AmrData&             amrDataf,
		  AmrData&             amrDatac,
		  Array<std::string>   cNames,
		  Array<Real>          barr,
		  std::string          oFile);
void
TakeDifferenceMean(AmrData&             amrDataf,
		   Array<std::string>   cNames,
		   Array<Real>          barr,
		   Array<int>           rratio,
		   std::string          oFile);
