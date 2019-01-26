/*
  testMemInfo.cpp: example usage of MemInfo class; see 
                   MemInfo.h for detailed explanations. 

                     Zarija Lukic, Berkeley, June 2013
                              zarija@lbl.gov
*/


#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "MemInfo.H"

int main(int argc, char **argv) {
  MPI_Init(NULL, NULL);

  // MemInfo class:
  MemInfo* mInfo = MemInfo::GetInstance();

  // Optional: initialize it with particular file name: 
  mInfo->Init("mem_info.log");

  // Test for the amount of memory in GB (all nodes):
  float chk_mem_gb = 2.0;
  int fail = mInfo->BelowThreshold(chk_mem_gb);
  if (fail > 0)
    fprintf(stderr, "WARNING: %d cores report memory shortage!\n", fail);

  // Print one line logs:
  for (int i = 0; i < 3; ++i) {
    void* p;
    if ((p=malloc(1<<28)) == NULL) {  // allocate 256MB
      puts("malloc() failed");
      break;
    }
    memset(p, 0, (1<<28));
    char info[32];
    snprintf(info, sizeof(info), "Step %d", i);
    mInfo->LogSummary(info);    // log memory status
  }

  // Print info for all MPI ranks to stdout:
  mInfo->PrintAll(stdout);

  // Get available and total memory for "my" node [GB]:
  float avail_mem, total_mem;
  mInfo->GetMemInfo(&avail_mem, &total_mem);

  MPI_Finalize();
  return(0);
}
