#include <AMReX_CuptiTrace.H>
#ifdef AMREX_USE_CUPTI
#include <AMReX.H>
#include <AMReX_Print.H>
#include <stdio.h>
#include <map>
#include <memory>
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

std::vector<std::unique_ptr<CUpti_Activity_Userdata>> activityRecordUserdata;

void CuptiInitialize ()
{
    initCuptiTrace();
    amrex::ExecOnFinalize(CuptiFinalize);
}

void CuptiFinalize ()
{
    activityRecordUserdata.clear();
}

void CUPTIAPI
bfrRequestCallback (uint8_t* *bfr, size_t* size, size_t* maxNumRecords) noexcept
{
    // Allocate a buffer for use by CUPTI; activity records are stored in the buffer
    uint8_t* buffer = (uint8_t*) std::malloc(BFR_SIZE + ALIGNMENT);
    if (buffer == NULL) {
        amrex::Abort("Error: CUPTI requested a buffer but memory cannot be allocated.");
    }

    // Return buffer and buffer size
    *bfr = ALIGN_BFR(buffer, ALIGNMENT);
    *size = BFR_SIZE;

    // Controls maximum number of records to be placed in buffer; setting to zero
    // fills the buffer with as many records as possible, setting to n > 0 fills
    // with n records before buffer is returned
    *maxNumRecords = 0;
}

void CUPTIAPI
bfrCompleteCallback (CUcontext ctx, uint32_t streamId, uint8_t* bfr,
                     size_t size, size_t validSize) noexcept
{ 
    CUptiResult status;
    CUpti_Activity* record = NULL;
    
    if (validSize > 0) {
        do {
            status = cuptiActivityGetNextRecord(bfr, validSize, &record);
            if (status == CUPTI_SUCCESS) {
                std::unique_ptr<CUpti_Activity_Userdata> recordUserData;
                recordUserData.reset(new CUpti_Activity_Userdata());
                CUpti_ActivityKernel4* kernel = (CUpti_ActivityKernel4*) record;
                
                // Save record data
                recordUserData->setStartTime(kernel->start);
                recordUserData->setEndTime(kernel->end);
                recordUserData->setTimeElapsed(kernel->end - kernel->start);
                recordUserData->setStreamID(kernel->streamId);
                recordUserData->setName(kernel->name);
                activityRecordUserdata.push_back( std::move(recordUserData) );
            }
            else if (status == CUPTI_ERROR_MAX_LIMIT_REACHED) {
                // No more records in the buffer
                break;
            }
            else if (status != CUPTI_SUCCESS) {
                const char* err;
                cuptiGetResultString(status, &err);
                amrex::Abort(err);
            }
        } while (true);

        size_t dropped;
        cuptiActivityGetNumDroppedRecords(ctx, streamId, &dropped);
        if (dropped != 0 and amrex::Verbose() > 1) {
            amrex::AllPrint() << (unsigned int) dropped
                              << " activity records were dropped due to insufficient buffer space\n";
        }
    }
    std::free(bfr);
}

void
initCuptiTrace () noexcept
{
    // Enable collection of kernel activity records
    cuptiActivityEnable(CUPTI_ACTIVITY_KIND_CONCURRENT_KERNEL);

    // Register callback functions to handle request of buffers to store
    // activity records and delivery of completed buffers to CUPTI client
    cuptiActivityRegisterCallbacks(bfrRequestCallback, bfrCompleteCallback);

    if (amrex::Verbose() > 0) {
        amrex::Print() << "CUPTI initialized\n";
    }
}
    
void
cuptiTraceStart () noexcept
{
    cudaDeviceSynchronize();
    cuptiActivityFlushAll(0);  
    activityRecordUserdata.clear();
}

void
cuptiTraceStop () noexcept
{
    cudaDeviceSynchronize();
    cuptiActivityFlushAll(0);
}

void
cuptiTraceStop (unsigned boxUintID) noexcept
{
    cudaDeviceSynchronize();
    cuptiActivityFlushAll(0);
    for (auto& record : activityRecordUserdata) {
        record->setUintID(boxUintID);
        record->setCharID( std::string("CharID_") + std::to_string(boxUintID) );
    }
}

double
computeElapsedTimeUserdata (const std::vector<std::unique_ptr<CUpti_Activity_Userdata>>&
                            activityRecordUserdata) noexcept
{
    if (activityRecordUserdata.size() == 0) {
        return 0.0;
    }
  
    std::map<int, unsigned long long> streamIDToElapsedTimeMap;
    
    // Initialize tally of unique streams
    for (auto& record : activityRecordUserdata) {
        if (streamIDToElapsedTimeMap.find(record->getStreamID())
            == streamIDToElapsedTimeMap.end()) {
            // Not found
            streamIDToElapsedTimeMap[record->getStreamID()] = 0;
        } else {
            // Found
        }
    }

    // Sum kernel times in each stream
    for (auto& record : activityRecordUserdata) {
        streamIDToElapsedTimeMap[record->getStreamID()] += record->getTimeElapsed();
    }

    // Sum over streams
    unsigned long long res = 0;
    for (auto const& kv : streamIDToElapsedTimeMap) {
        res += kv.second;
    }  

    // Average time per kernel
    res /= (1.*activityRecordUserdata.size());

    // Alternative: average time per stream
    //res /= streamIDToElapsedTimeMap.size();
  
    // Default is ns, convert to sec
    return (double) res*1e-9;
}

void
CUpti_Activity_Userdata::setUintID (unsigned uintID) noexcept
{
    uintID_ = uintID;
}

void
CUpti_Activity_Userdata::setStartTime (unsigned long long startTime) noexcept
{
    startTime_ = startTime;
}

void
CUpti_Activity_Userdata::setEndTime (unsigned long long endTime) noexcept
{
    endTime_ = endTime;
}

void
CUpti_Activity_Userdata::setTimeElapsed (unsigned long long timeElapsed) noexcept
{
    timeElapsed_ = timeElapsed;
}

void
CUpti_Activity_Userdata::setStreamID (int streamID) noexcept
{
    streamID_ = streamID;
}

void
CUpti_Activity_Userdata::setName (std::string name) noexcept
{
    name_ = std::move(name);
}

unsigned
CUpti_Activity_Userdata::getUintID () const noexcept
{ 
    return uintID_;
}

void
CUpti_Activity_Userdata::setCharID (std::string charID) noexcept
{
    charID_ = std::move(charID);
}

std::string const&
CUpti_Activity_Userdata::getCharID () const noexcept
{ 
    return charID_;
}

unsigned long long
CUpti_Activity_Userdata::getStartTime () const noexcept
{
    return startTime_;
}

unsigned long long
CUpti_Activity_Userdata::getEndTime () const noexcept
{
    return endTime_;
}

unsigned long long
CUpti_Activity_Userdata::getTimeElapsed () const noexcept
{
    return timeElapsed_;
}

int
CUpti_Activity_Userdata::getStreamID () const noexcept
{
    return streamID_;
}

std::string const&
CUpti_Activity_Userdata::getName () const noexcept
{
    return name_;
}
  
CuptiTrace::CuptiTrace () noexcept
{
}

CuptiTrace::~CuptiTrace () noexcept
{
}
 
void
CuptiTrace::start () noexcept
{
    cuptiTraceStart();
}

void
CuptiTrace::stop () noexcept
{
    cuptiTraceStop();
}

void
CuptiTrace::stop (unsigned boxUintID) noexcept
{
    cuptiTraceStop(boxUintID);
}

}

#endif
