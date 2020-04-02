#ifndef AMREX_CUPTI_TRACE_H_
#define AMREX_CUPTI_TRACE_H_

#ifdef AMREX_USE_CUPTI
#include <stdio.h>
#include <vector>
#include <cuda.h>
#include <cupti.h>
#include <memory>


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
    unsigned uintID_;
    std::string charID_;
    unsigned long long startTime_;
    unsigned long long endTime_;
    unsigned long long timeElapsed_;
    int streamID_;
    std::string name_;
  
public:
    void setUintID (unsigned) noexcept;
    void setCharID (std::string) noexcept;
    void setStartTime (unsigned long long) noexcept;
    void setEndTime (unsigned long long) noexcept;
    void setTimeElapsed (unsigned long long) noexcept;
    void setStreamID (int) noexcept;
    void setName (std::string) noexcept;
    unsigned getUintID () noexcept;
    std::string getCharID () noexcept;
    unsigned long long getStartTime () noexcept;
    unsigned long long getEndTime () noexcept;
    unsigned long long getTimeElapsed () noexcept;
    int getStreamID () noexcept;
    std::string getName () noexcept;
};

class CuptiTrace {
public:
    CuptiTrace () noexcept;
    ~CuptiTrace () noexcept;
    void start () noexcept;
    void stop () noexcept;
    void stop (unsigned boxUintID) noexcept;
};

extern std::vector<std::unique_ptr<CUpti_Activity_Userdata>> activityRecordUserdata;
double computeElapsedTimeUserdata(const std::vector<std::unique_ptr<CUpti_Activity_Userdata>>&
                                  activityRecordUserdata) noexcept;
}

#endif
#endif
