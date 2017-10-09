#include <AMReX_DataServices.H>

void
ComputeAmrDataMeanVar (AmrData&     amrData,
		       Vector<Real>& mean,
		       Vector<Real>& variance,
		       int          sComp,
		       int          nComp,
		       bool         verbose);

void
ComputeAmrDataMeanVar (AmrData&           amrData,
		       Vector<std::string> cNames,
		       Vector<BoxArray>    bas,
		       Vector<Real>&       mean,
		       Vector<Real>&       variance);
void
ComputeMeanVarMF (MultiFab&          mf,
		  Vector<Real>&       mean,
		  Vector<Real>&       variance);

void
ComputeAmrDataList (AmrData&         amrData,
		    Vector<MultiFab*> mean,
		    Vector<MultiFab*> variance,
		    int              sComp,
		    int              nComp);

void
ComputeAmrDataPDF (AmrData&     amrData,
		   Real**       icount,
		   int          nBin,
		   Vector<std::string> cNames);

void
ComputeAmrDataPDF (AmrData&     amrData,
		   Real**       icount,
		   int          nBin,
		   Vector<std::string> cNames,
		   Vector<BoxArray> bas);

void
ComputeAmrDataVAR (AmrData&           amrData,
		   int                nBin,
		   Vector<std::string> cNames,
		   Vector<Real>        barr,
		   std::string        oFile);

void
VariogramUniform (AmrData&             amrData,
		  Vector<std::string>   cNames,
		  Vector<Real>          barr,
		  Vector< Vector<int> >  ivoption,
		  int                  nlag,
		  int                  isill,
		  Vector<Real>          sills,
		  std::string          oFile);

void
VariogramCross(AmrData&             amrData,
	       Vector<std::string>   cNames,
	       MultiFab&            mf,
	       Vector<Real>          barr,
	       Vector< Vector<int> >  ivoption,
	       int                  nlag,
	       int                  isill,
	       Vector<Real>          sills,
	       std::string          oFile);

void
VariogramUniformMF (const MultiFab&      mf,
		    Vector<Real>          dx,
		    Vector< Vector<int> >  ivoption,
		    int                  nlag,
		    int                  isill,
		    Vector<Real>          sills,
		    std::string          oFile);
void
VariogramUniformMFG (const MultiFab&      mf,
		     Vector<Real>          dx,
		     IntVect              sm,
		     IntVect              bg,
		     Vector< Vector<int> >  ivoption,
		     int                  nlag,
		     int                  isill,
		     Vector<Real>          sills,
		     std::string          oFile);
void
TakeDifferenceFine(AmrData&             amrDataf,
		   AmrData&             amrDatac,
		   Vector<std::string>   cNames,
		   Vector<Real>          barr,
		   std::string          oFile);

void
TakeDifferenceCrse(AmrData&             amrDataf,
		   AmrData&             amrDatac,
		   Vector<std::string>   cNames,
		   Vector<Real>          barr,
		   std::string          oFile);
void
TakeDifferenceSum(AmrData&             amrDataf,
		  AmrData&             amrDatac,
		  Vector<std::string>   cNames,
		  Vector<Real>          barr,
		  std::string          oFile);
void
TakeDifferenceMean(AmrData&             amrDataf,
		   Vector<std::string>   cNames,
		   Vector<Real>          barr,
		   Vector<int>           rratio,
		   std::string          oFile);
