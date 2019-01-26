/*
  MemInfo: Utility class for checking/logging 
       free memory on nodes in MPI runtime environment.
       See also testMemInfo.cpp for example usage.

  Systems: 
       NOTE THAT THIS CODE IS MACHINE/OS DEPENDENT!
       Should work on Linux/Unix and Mac. It certainly 
       will not work on a Windows system.

  Usage:
       0) To get class instance (singleton class!): 

        MemInfo* mInfo=MemInfo::GetInstance();

          Optional: to specify file for log summaries 
          (dafault file name: mem_info.log):

          mInfo->Init("filename");

       1) To check how many nodes have less than 
        a desired amount of GBs (threshold):

        int fail = mInfo->BelowThreshold(float threshold);

       2) To write one line summary to mem_info.log 
        file, containing min, max and average available 
        memory (info: string to be printed out, e.g. name 
        of the routine writing the log line):

        mInfo->LogSummary(char* info);

       3) To write full memory info for all MPI ranks to 
        a file fname (output is ordered by MPI ranks): 

        mInfo->PrintAll(FILE* fname);

       4) To get available/free and total memory for "my" 
        node:

        mInfo->GetMemInfo(&availMem, &totalMem);


                     Zarija Lukic, Berkeley, June 2013
                              zarija@lbl.gov
*/


#ifndef _MemInfo_H_
#define _MemInfo_H_

#include <string>

class MemInfo {
 public:
  static MemInfo* GetInstance() {
    if (!classInstance) classInstance = new MemInfo;
    return(classInstance);
  }

  void Init(const std::string& fname);
  int  BelowThreshold(float threshold);
  void LogSummary(const char* info);
  void PrintAll(FILE* fname);
  bool GetMemInfo(float* availMem, float* totalMem);

 private:
  static const int kStrLen = 512;
  static const int64_t kKB = 1024;
  static const int64_t kMB = kKB*1024;
  static const int64_t kGB = kMB*1024;

  MemInfo();
  ~MemInfo();
  explicit MemInfo(MemInfo const& copy);
  MemInfo& operator=(MemInfo const& copy);
  static MemInfo* classInstance;

  bool MPIData();
  bool GetMyHostName();
  bool GetMyPages();
  void StampSysInfo(FILE* fname);

  std::string default_fname_;
  bool class_initialized_;
  int nlen_;
  char* hostname_;
  long page_size_, total_pages_, avail_pages_;
  float global_min_, global_max_, psize_in_gb_;

  FILE* logfile_;

  int mpi_v_, mpi_subv_;                   // MPI
  int my_rank_, num_ranks_, master_rank_;  // stuff
};

#endif  // _MemInfo_H_
