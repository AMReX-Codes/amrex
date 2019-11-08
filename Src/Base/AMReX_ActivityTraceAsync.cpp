/*
 * Downsized version of activity_trace_async.cpp in:
 *     $(CUDA_HOME)/extras/CUPTI/samples 
 *
 * Copyright 2011-2015 NVIDIA Corporation. All rights reserved
 *
 * Sample CUPTI app to print a trace of CUDA API and GPU activity
 * using asynchronous handling of activity buffers.
 *
 */

#include <AMReX_ActivityTraceAsync.H>
#ifdef AMREX_USE_CUPTI
#include <stdio.h>
#include <map>
#include <cuda.h>
#include <cupti.h>


namespace amrex {
  
std::vector<CUpti_Activity_Userdata*> activityRecordUserdata;

const char *
getActivityOverheadKindString(CUpti_ActivityOverheadKind kind) noexcept
{
  switch (kind) {
  case CUPTI_ACTIVITY_OVERHEAD_DRIVER_COMPILER:
    return "COMPILER";
  case CUPTI_ACTIVITY_OVERHEAD_CUPTI_BUFFER_FLUSH:
    return "BUFFER_FLUSH";
  case CUPTI_ACTIVITY_OVERHEAD_CUPTI_INSTRUMENTATION:
    return "INSTRUMENTATION";
  case CUPTI_ACTIVITY_OVERHEAD_CUPTI_RESOURCE:
    return "RESOURCE";
  default:
    break;
  }

  return "<unknown>";
}

const char *
getActivityObjectKindString(CUpti_ActivityObjectKind kind) noexcept
{
  switch (kind) {
  case CUPTI_ACTIVITY_OBJECT_PROCESS:
    return "PROCESS";
  case CUPTI_ACTIVITY_OBJECT_THREAD:
    return "THREAD";
  case CUPTI_ACTIVITY_OBJECT_DEVICE:
    return "DEVICE";
  case CUPTI_ACTIVITY_OBJECT_CONTEXT:
    return "CONTEXT";
  case CUPTI_ACTIVITY_OBJECT_STREAM:
    return "STREAM";
  default:
    break;
  }

  return "<unknown>";
}

uint32_t
getActivityObjectKindId(CUpti_ActivityObjectKind kind, CUpti_ActivityObjectKindId *id) noexcept
{
  switch (kind) {
  case CUPTI_ACTIVITY_OBJECT_PROCESS:
    return id->pt.processId;
  case CUPTI_ACTIVITY_OBJECT_THREAD:
    return id->pt.threadId;
  case CUPTI_ACTIVITY_OBJECT_DEVICE:
    return id->dcs.deviceId;
  case CUPTI_ACTIVITY_OBJECT_CONTEXT:
    return id->dcs.contextId;
  case CUPTI_ACTIVITY_OBJECT_STREAM:
    return id->dcs.streamId;
  default:
    break;
  }

  return 0xffffffff;
}

void printActivity(CUpti_Activity *record) noexcept
{
  switch (record->kind)
  {
    case CUPTI_ACTIVITY_KIND_KERNEL:
    case CUPTI_ACTIVITY_KIND_CONCURRENT_KERNEL:
    {
      const char* kindString = (record->kind == CUPTI_ACTIVITY_KIND_KERNEL) ? "KERNEL" : "CONC KERNEL";
      CUpti_ActivityKernel4 *kernel = (CUpti_ActivityKernel4 *) record;
      printf("%s \"%s\" [ %llu - %llu, %llu ns] device %u, context %u, stream %u, correlation %u\n\n",
             kindString,
             kernel->name,
             (unsigned long long) (kernel->start - startTimestamp),
             (unsigned long long) (kernel->end - startTimestamp),
	     (unsigned long long) (((kernel->end - startTimestamp) - (kernel->start - startTimestamp))),
             kernel->deviceId, kernel->contextId, kernel->streamId,
             kernel->correlationId);
      printf("    grid [%u,%u,%u], block [%u,%u,%u], shared memory (static %u, dynamic %u)\n",
             kernel->gridX, kernel->gridY, kernel->gridZ,
             kernel->blockX, kernel->blockY, kernel->blockZ,
             kernel->staticSharedMemory, kernel->dynamicSharedMemory);
      break;
    }
  case CUPTI_ACTIVITY_KIND_OVERHEAD:
    {
      CUpti_ActivityOverhead *overhead = (CUpti_ActivityOverhead *) record;
      printf("OVERHEAD %s [ %llu, %llu ] %llu ns, %s id %u\n",
             getActivityOverheadKindString(overhead->overheadKind),
             (unsigned long long) overhead->start - startTimestamp,
             (unsigned long long) overhead->end - startTimestamp,
	     (unsigned long long) (((overhead->end - startTimestamp) - (overhead->start - startTimestamp))),
             getActivityObjectKindString(overhead->objectKind),
             getActivityObjectKindId(overhead->objectKind, &overhead->objectId));
      break;
    }
  default:
    const char* kindString = (record->kind == CUPTI_ACTIVITY_KIND_KERNEL) ? "KERNEL" : "CONC KERNEL";
    CUpti_ActivityKernel4 *kernel = (CUpti_ActivityKernel4 *) record;
    printf("%s \"%s\" [ %llu - %llu, %llu ns] device %u, context %u, stream %u, correlation %u\n\n",
	   kindString,
	   kernel->name,
	   (unsigned long long) (kernel->start - startTimestamp),
	   (unsigned long long) (kernel->end - startTimestamp),
	   (unsigned long long) (((kernel->end - startTimestamp) - (kernel->start - startTimestamp))),
	   kernel->deviceId, kernel->contextId, kernel->streamId,
	   kernel->correlationId);
    printf("    grid [%u,%u,%u], block [%u,%u,%u], shared memory (static %u, dynamic %u)\n",
	   kernel->gridX, kernel->gridY, kernel->gridZ,
	   kernel->blockX, kernel->blockY, kernel->blockZ,
	   kernel->staticSharedMemory, kernel->dynamicSharedMemory);
    break;
    printf("  <unknown>\n");
    break;
  }
}

void CUPTIAPI bufferRequested(uint8_t **buffer, size_t *size, size_t *maxNumRecords) noexcept
{
  uint8_t *bfr = (uint8_t *) malloc(BUF_SIZE + ALIGN_SIZE);
  if (bfr == NULL) {
    printf("Error: out of memory\n");
    exit(-1);
  }

  *size = BUF_SIZE;
  *buffer = ALIGN_BUFFER(bfr, ALIGN_SIZE);
  *maxNumRecords = 0;
}

void CUPTIAPI bufferCompleted(CUcontext ctx, uint32_t streamId, uint8_t *buffer, size_t size, size_t validSize) noexcept
{ 
  CUptiResult status;
  CUpti_Activity *record = NULL;
  
  if (validSize > 0) {
    do {
      status = cuptiActivityGetNextRecord(buffer, validSize, &record);
      if (status == CUPTI_SUCCESS) {
	CUpti_Activity_Userdata* recordUserData = (CUpti_Activity_Userdata*) record;
	activityRecordUserdata.push_back(recordUserData);
	//printActivity(record);
      }
      else if (status == CUPTI_ERROR_MAX_LIMIT_REACHED)
        break;
      else {
        CUPTI_CALL(status);
      }
    } while (1);
    
    // Report any records dropped from the queue
    size_t dropped;
    CUPTI_CALL(cuptiActivityGetNumDroppedRecords(ctx, streamId, &dropped));
    if (dropped != 0) {
      printf("Dropped %u activity records\n", (unsigned int) dropped);
    }
  }

  free(buffer);
}

void
initTrace() noexcept
{
  // 1) CUPTI activity record is created when CUDA initializes, so must enable
  //    before cuInit() or any CUDA runtime call
  // 2) Register callback functions
  // 3) Get and set activity API attributes, for further control over behavior
  //    of the activity API

  size_t attrValue = 0, attrValueSize = sizeof(size_t);
  
  // Kernel activity records
  CUPTI_CALL(cuptiActivityEnable(CUPTI_ACTIVITY_KIND_KERNEL));
  
  // Overhead activity records
  //CUPTI_CALL(cuptiActivityEnable(CUPTI_ACTIVITY_KIND_OVERHEAD));

  // Register callbacks for buffer requests and for buffers completed by CUPTI
  CUPTI_CALL(cuptiActivityRegisterCallbacks(bufferRequested, bufferCompleted));

  // Get and set activity attributes; allows for control over behavior of the
  // activity API.  Some attributes require to be set before any CUDA context is
  // created to be effective, e.g. to be applied to all device buffer allocations.
  CUPTI_CALL(cuptiActivityGetAttribute(CUPTI_ACTIVITY_ATTR_DEVICE_BUFFER_SIZE, &attrValueSize, &attrValue));
  printf("%s = %llu\n", "CUPTI_ACTIVITY_ATTR_DEVICE_BUFFER_SIZE", (long long unsigned)attrValue);
  attrValue *= 2; // Just for example
  CUPTI_CALL(cuptiActivitySetAttribute(CUPTI_ACTIVITY_ATTR_DEVICE_BUFFER_SIZE, &attrValueSize, &attrValue));
  
  // Maximum number of buffers per context; buffers can be reused by the context.
  // Increasing this reduces the number of times CUPTI needs to flush the buffers.
  // Default value is 100.
  CUPTI_CALL(cuptiActivityGetAttribute(CUPTI_ACTIVITY_ATTR_DEVICE_BUFFER_POOL_LIMIT, &attrValueSize, &attrValue));
  printf("%s = %llu\n", "CUPTI_ACTIVITY_ATTR_DEVICE_BUFFER_POOL_LIMIT", (long long unsigned)attrValue);
  attrValue *= 2; // Just for example.
  CUPTI_CALL(cuptiActivitySetAttribute(CUPTI_ACTIVITY_ATTR_DEVICE_BUFFER_POOL_LIMIT, &attrValueSize, &attrValue));

  // Time at initialization, used for normalizing other printed times
  CUPTI_CALL(cuptiGetTimestamp(&startTimestamp));
}

double
computeElapsedTimeUserdata(std::vector<CUpti_Activity_Userdata*> activityRecordUserdata) noexcept
{
  unsigned long long t_elapsed = 0;
  unsigned long long t_start = 0;
  unsigned long long t_stop = 0;

  std::map<int, unsigned long long> streamIDToElapsedTimeMap;
  
  // Initialize tally of unique streams
  for (auto record : activityRecordUserdata) {
    CUpti_ActivityKernel4 *kernel = (CUpti_ActivityKernel4 *) record;
    if (streamIDToElapsedTimeMap.find(kernel->streamId) == streamIDToElapsedTimeMap.end()) {
      // Not found
      streamIDToElapsedTimeMap[kernel->streamId] = 0;
    } else {
      // Found
    }
  }

  // Sum kernel times in each stream
  for (auto record : activityRecordUserdata) {
    CUpti_ActivityKernel4 *kernel = (CUpti_ActivityKernel4 *) record;
    t_start = (unsigned long long) (kernel->start - startTimestamp);
    t_stop = (unsigned long long) (kernel->end - startTimestamp);
    t_elapsed = (((unsigned long long)t_stop) - ((unsigned long long)t_start));
    streamIDToElapsedTimeMap[kernel->streamId] += t_elapsed;
  }

  // Average over streams
  unsigned long long res = 0;
  for (auto const& kv : streamIDToElapsedTimeMap) {
    res += kv.second;
  }
  res /= streamIDToElapsedTimeMap.size();
  
  // Default is ns, convert to sec  
  return (double) res*1e-9;
    
}

void
CUpti_Activity_Userdata::setUintID (unsigned pUintID) noexcept
{
  uintID = pUintID;
}

unsigned
CUpti_Activity_Userdata::getUintID () noexcept
{
  return uintID;
}

void
CUpti_Activity_Userdata::setCharID (const char* pCharID) noexcept
{
  charID = pCharID;
}

const char*
CUpti_Activity_Userdata::getCharID () noexcept
{
  return charID;
}

}

#endif // AMREX_USE_CUPTI
