#include "myfunc.H"
#include "mykernel.H"

#ifdef AMREX_USE_CUPTI
#include <AMReX_ActivityTraceAsync.H>
#endif // AMREX_USE_CUPTI

void doRandomSleep (MultiFab& mf) {
  // Call device sleep function at each iteration of MFIter loop
#ifdef AMREX_USE_CUPTI
  BL_PROFILE_VAR_NS_CUPTI("randomSleep", blp_sleep);
  int iter = 0;
#else
  BL_PROFILE_VAR_NS("randomSleep", blp_sleep);
#endif // AMREX_USE_CUPTI
  for ( MFIter mfi(mf); mfi.isValid(); ++mfi )
    {
      // Code is instrumented here to capture CUPTI kernel activity
#ifdef AMREX_USE_CUPTI
      BL_PROFILE_VAR_START_CUPTI(blp_sleep);
#else
      BL_PROFILE_VAR_START(blp_sleep);
#endif // AMREX_USE_CUPTI
      randomSleep<<<1, 1, 0, amrex::Gpu::Device::gpuStream()>>>();
#ifdef AMREX_USE_CUPTI
      BL_PROFILE_VAR_STOP_CUPTI(blp_sleep, iter,
				("char ID " + std::to_string(iter)).c_str());
      iter +=1;
#else
      BL_PROFILE_VAR_STOP(blp_sleep)
#endif // AMREX_USE_CUPTI

#ifdef AMREX_USE_CUPTI
      // Kernel data captured in a vector `activityRecordUserdata`
      // Print here some information from the captured records
      unsigned long long t_start = 0;
      unsigned long long t_stop = 0;
      
      for (auto record : activityRecordUserdata) {
	CUpti_ActivityKernel4 *kernel = (CUpti_ActivityKernel4 *) record;
	t_start = (unsigned long long) (kernel->start - startTimestamp);
	t_stop = (unsigned long long) (kernel->end - startTimestamp);
	
	unsigned long long dt = 0;
	dt = (((unsigned long long)t_stop) - ((unsigned long long)t_start));
	
	printf("MF Iteration: %u; Box name: %s; t_elapsed: %llu ns\n",
	       record->getUintID(), record->getCharID(), dt);
      }
#endif // AMREX_USE_CUPTI
    }
}

void init_mf(amrex::MultiFab& mf, amrex::Geometry const& geom){
    // Initialize the MultiFab; MFIter = MultiFab Iterator
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
      {
      const Box& vbx = mfi.validbox();
      auto const& mfTmp = mf.array(mfi);
      AMREX_FOR_3D ( vbx, i, j, k,
      {
	init_mf(i,j,k,mfTmp);
      });
    }
}
