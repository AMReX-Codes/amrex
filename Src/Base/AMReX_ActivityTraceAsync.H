#ifndef AMREX_ACTIVITY_TRACE_ASYNC_H_
#define AMREX_ACTIVITY_TRACE_ASYNC_H_

#ifdef AMREX_USE_CUPTI
#include <stdio.h>
#include <vector>
#include <cuda.h>
#include <cupti.h>

// Error checking
#define CUPTI_CALL(call)                                                   \
  do {                                                                     \
    CUptiResult _status = call;                                            \
    if (_status != CUPTI_SUCCESS) {                                        \
      const char *errstr;                                                  \
      cuptiGetResultString(_status, &errstr);                              \
      fprintf(stderr, "%s:%d: error: function %s failed with error %s.\n", \
              __FILE__, __LINE__, #call, errstr);                          \
      if(_status == CUPTI_ERROR_UM_PROFILING_NOT_SUPPORTED)                \
          exit(0);                                                         \
      else                                                                 \
          exit(-1);                                                        \
    }                                                                      \
  } while (0)

#define BUF_SIZE (32*1024)
#define ALIGN_SIZE (8)

// Keep buffer size divisible by required align size
#define ALIGN_BUFFER(buffer, align)                                        \
  (((uintptr_t) (buffer) & ((align)-1)) ? ((buffer) + (align) - ((uintptr_t) (buffer) & ((align)-1))): (buffer))


namespace amrex {

static uint64_t startTimestamp;
  
const char* getActivityOverheadKindString (CUpti_ActivityOverheadKind kind) noexcept;
const char* getActivityObjectKindString (CUpti_ActivityObjectKind kind) noexcept;
uint32_t getActivityObjectKindId (CUpti_ActivityObjectKind kind,
				  CUpti_ActivityObjectKindId *id) noexcept;

void printActivity (CUpti_Activity *record) noexcept;

void CUPTIAPI bufferRequested (uint8_t **buffer, size_t *size,
			      size_t *maxNumRecords) noexcept;
void CUPTIAPI bufferCompleted (CUcontext ctx, uint32_t streamId, uint8_t *buffer,
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
  
public:
  void setRecord (CUpti_Activity*) noexcept;
  void setUintID (unsigned) noexcept;
  void setCharID (const char*) noexcept;
  CUpti_Activity* getRecord () noexcept;
  unsigned getUintID () noexcept;
  const char* getCharID () noexcept;
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
