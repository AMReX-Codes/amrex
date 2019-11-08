#ifndef AMREX_CUPTI_TRACE_H_
#define AMREX_CUPTI_TRACE_H_

#ifdef AMREX_USE_CUPTI
#include <stdio.h>
#include <vector>
#include <cuda.h>
#include <cupti.h>

// CUPTI buffer size, enough for 4096 activity records in a single buffer;
// `CUpti_Activity` objects are 8 bytes long
#define BFR_SIZE (32768)

// 8-byte alignment
#define ALIGNMENT (8)

// Round down to `align` boundary; mask out last 3 bits (address therefore ends in a 0 or 8)
#define ALIGN_BFR(bfr, align) ( (uint8_t*) (((uintptr_t)bfr + (align - 1)) & ~ (uintptr_t)(align - 1)) )


namespace amrex {

void CUPTIAPI bfrRequestCallback (uint8_t* *bfr, size_t* size,
				  size_t* maxNumRecords) noexcept;
void CUPTIAPI bfrCompleteCallback (CUcontext ctx, uint32_t streamId, uint8_t* bfr,
				   size_t size, size_t validSize) noexcept;

void initCuptiTrace () noexcept;
void cuptiTraceStart () noexcept;
void cuptiTraceStop () noexcept;
void cuptiTraceStop (unsigned boxUintID) noexcept;

class CUpti_Activity_Userdata {
private:
  CUpti_Activity* record_;
  unsigned uintID_;
  const char* charID_;
  unsigned long long startTime_;
  unsigned long long endTime_;
  unsigned long long timeElapsed_;
  int streamID_;
  const char* name_;
  
public:
  void setRecord (CUpti_Activity*) noexcept;
  void setUintID (unsigned) noexcept;
  void setCharID (const char*) noexcept;
  void setStartTime (unsigned long long) noexcept;
  void setEndTime (unsigned long long) noexcept;
  void setTimeElapsed (unsigned long long) noexcept;
  void setStreamID (int) noexcept;
  void setName (const char*) noexcept;
  CUpti_Activity* getRecord () noexcept;
  unsigned getUintID () noexcept;
  const char* getCharID () noexcept;
  unsigned long long getStartTime () noexcept;
  unsigned long long getEndTime () noexcept;
  unsigned long long getTimeElapsed () noexcept;
  int getStreamID () noexcept;
  const char* getName () noexcept;
};

class CuptiTrace {
public:
  CuptiTrace () noexcept;
  ~CuptiTrace () noexcept;
  void start () noexcept;
  void stop () noexcept;
  void stop (unsigned boxUintID) noexcept;
  void cleanup () noexcept;
};    

extern std::vector<CUpti_Activity_Userdata*> activityRecordUserdata;
double computeElapsedTimeUserdata(std::vector<CUpti_Activity_Userdata*>
				  activityRecordUserdata) noexcept;
}

#endif // AMREX_USE_CUPTI

#endif
