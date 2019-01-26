/*
  MemInfo.cpp: Implementation of the class; see 
               MemInfo.h for detailed explanations, 
               and testMemInfo.cpp for example usage. 

                   Zarija Lukic, Berkeley, June 2013
                            zarija@lbl.gov
*/


#include <mpi.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
#include <sys/utsname.h>

#if defined __APPLE__
#include <sys/sysctl.h>
#include <mach/mach.h>
#endif

#include <cstdlib>
#include <cstdio>
#include <string>

#include <AMReX_ParallelDescriptor.H>

#include "MemInfo.H"

using namespace amrex;

MemInfo* MemInfo::classInstance = NULL;


MemInfo::MemInfo() {
  class_initialized_ = false;
  default_fname_ = "mem_info.log";
  bool success = MPIData();
  if ( ! success) {
     puts("MemInfo::MemInfo(): MPI not initialized!");
  }
  return;
}


MemInfo::~MemInfo() {
  delete hostname_;
  if (my_rank_ == master_rank_ && logfile_ != NULL)
    fclose(logfile_);
  return;
}


void MemInfo::Init(const std::string& fname) {
  if (class_initialized_)
    return;
  class_initialized_ = true;

  // Get the hostname:
  errno = 0;
  if ((nlen_ = sysconf(_SC_HOST_NAME_MAX)) == -1) {
    if (errno == 0) {
      printf("Process %d: HOST_NAME_MAX not supported\n", my_rank_);
    } else {
      fprintf(stderr, "Process %d: ", my_rank_);
      perror("MemInfo::Init(): problem with host name max");
    }
    nlen_ = 64;
  }

  hostname_ = new char[nlen_];
  bool success = GetMyHostName();
  if (!success)
     printf("Process %d: could not get the hostname\n", my_rank_);

  // Get page size and total number of pages on the node:
  errno = 0;
  if ((page_size_ = sysconf(_SC_PAGE_SIZE)) == -1) {
    if (errno == 0) {
      printf("Process %d: PAGE_SIZE not supported\n", my_rank_);
    } else {
      fprintf(stderr, "Process %d: ", my_rank_);
      perror("MemInfo::Init(): problem with page size");
    }
  }
  psize_in_gb_ = static_cast<float>(page_size_)/kGB;

#if defined __APPLE__
  int mib[2] = { CTL_HW, HW_MEMSIZE };
  uint64_t memsize;
  size_t len = sizeof(memsize);
  if (sysctl(mib, 2, &memsize, &len, NULL, 0) == -1) {
    fprintf(stderr, "Process %d: ", my_rank_);
    perror("MemInfo::Init(): problem in sysctl call");
    total_pages_ = 0;
  } else {
    total_pages_ = memsize/page_size_;
  }
#else
  errno = 0;
  if ((total_pages_ = sysconf(_SC_PHYS_PAGES)) == -1) {
    if (errno == 0) {
      printf("Process %d: PHYS_PAGES not supported\n", my_rank_);
    } else {
      fprintf(stderr, "Process %d: ", my_rank_);
      perror("MemInfo::Init(): problem with phys pages");
    }
    total_pages_ = 0;
  }
#endif  // __APPLE__

  // Find min and max across all nodes:
  long min_pag = 0, max_pag = 0;
  BL_MPI_REQUIRE( MPI_Allreduce(&total_pages_, &min_pag, 1, MPI_LONG, MPI_MIN,
                ParallelDescriptor::Communicator()) );
  BL_MPI_REQUIRE( MPI_Allreduce(&total_pages_, &max_pag, 1, MPI_LONG, MPI_MAX,
                ParallelDescriptor::Communicator()) );
  global_min_ = static_cast<float>(min_pag) * psize_in_gb_;
  global_max_ = static_cast<float>(max_pag) * psize_in_gb_;

  // Open log file:
  if (my_rank_ == master_rank_) {
    std::string filename;
    if (!fname.empty())
      filename = fname;
    else
      filename = default_fname_;

    if ((logfile_ = fopen(filename.c_str(), "w")) == NULL)
      perror("MemInfo::Init(): Cannot open memory log file");
    else
      StampSysInfo(logfile_);
  }

  return;
}


bool MemInfo::MPIData() {
  bool success = false;
  int flag = 0;
  BL_MPI_REQUIRE( MPI_Initialized(&flag) );

  if (flag) {
    BL_MPI_REQUIRE( MPI_Comm_size(ParallelDescriptor::Communicator(), &num_ranks_) );
    BL_MPI_REQUIRE( MPI_Comm_rank(ParallelDescriptor::Communicator(), &my_rank_) );
    master_rank_ = 0;
    BL_MPI_REQUIRE( MPI_Get_version(&mpi_v_, &mpi_subv_) );
    success = true;
  } else {
    my_rank_ = master_rank_ = 0;
    num_ranks_ = 1;
    mpi_v_ = -1;
    mpi_subv_ = -1;
    success = false;
  }

  return(success);
}


bool MemInfo::GetMyHostName() {
  bool success = false;

  if (gethostname(hostname_, nlen_) == 0)
    success = true;
  else
    *hostname_ = '\0';

  return(success);
}


bool MemInfo::GetMyPages() {
  bool success = true;

#if defined __APPLE__
  vm_statistics_data_t   host_info;
  mach_msg_type_number_t count = HOST_VM_INFO_COUNT;
  if (host_statistics(mach_host_self(), HOST_VM_INFO,
      reinterpret_cast<host_info_t>(&host_info), &count) != 0) {
    fprintf(stderr, "Process %d: ", my_rank_);
    perror("MemInfo::GetMyPages(): host_statistics error");
    avail_pages_ = 0;
    success = false;
  } else {
    avail_pages_ = host_info.free_count;
  }
#else
  errno = 0;
  if ((avail_pages_ = sysconf(_SC_AVPHYS_PAGES)) == -1) {
    if (errno == 0) {
      printf("Process %d: AVPHYS_PAGES not supported\n", my_rank_);
    } else {
      fprintf(stderr, "Process %d: ", my_rank_);
      perror("MemInfo::GetMyPages() error");
    }
    avail_pages_ = 0;
    success = false;
  }
#endif  // __APPLE__

  return(success);
}


void MemInfo::StampSysInfo(FILE* fname) {
  struct utsname info;
  if ((uname(&info) < 0) || (fname == NULL))
    return;

  fprintf(fname, "#################################");
  fprintf(fname, "#################################\n");

  char buf[kStrLen], t_stamp[64];
  const time_t t_now = time(NULL);
  struct tm* current_time = new tm;
  current_time = localtime_r(&t_now, current_time);
  asctime_r(current_time, t_stamp);

  int len = 0;
  len += snprintf(buf+len, kStrLen-len, "# Time: %s",      t_stamp);
  len += snprintf(buf+len, kStrLen-len, "# System: %s ",   info.sysname);
  len += snprintf(buf+len, kStrLen-len, "(%s), ",          info.release);
  len += snprintf(buf+len, kStrLen-len, "%s\n",            info.machine);
  len += snprintf(buf+len, kStrLen-len, "# Version: %s\n", info.version);
  len += snprintf(buf+len, kStrLen-len, "# Total memory on nodes: ");
  len += snprintf(buf+len, kStrLen-len, "min=%5.2fGB, max=%5.2fGB\n",
           global_min_, global_max_);

  fprintf(fname, "%s", buf);
  fprintf(fname, "#################################");
  fprintf(fname, "#################################\n\n");

  return;
}


int MemInfo::BelowThreshold(float threshold) {
  if (!class_initialized_)
    Init(default_fname_);
  bool success = GetMyPages();
  if (!success)
    return(-1);

  int i = ((avail_pages_ * psize_in_gb_) >= threshold);
  if (i == 0)
    printf("Process %d: node %s is short on memory\n",
           my_rank_, hostname_);

  int sum = 0;
  BL_MPI_REQUIRE( MPI_Allreduce(&i, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD) );

  return(num_ranks_ - sum);
}


void MemInfo::LogSummary(const char* info) {
  if ( ! class_initialized_) {
    Init(default_fname_);
  }
  bool success = GetMyPages();
  if ( ! success) {
    return;
  }

  long avg_pag = 0, min_pag = 0, max_pag = 0;
  BL_MPI_REQUIRE( MPI_Reduce(&avail_pages_, &min_pag, 1, MPI_LONG, MPI_MIN, master_rank_,
             ParallelDescriptor::Communicator()) );
  BL_MPI_REQUIRE( MPI_Reduce(&avail_pages_, &max_pag, 1, MPI_LONG, MPI_MAX, master_rank_,
             ParallelDescriptor::Communicator()) );
  BL_MPI_REQUIRE( MPI_Reduce(&avail_pages_, &avg_pag, 1, MPI_LONG, MPI_SUM, master_rank_,
             ParallelDescriptor::Communicator()) );

  if (my_rank_ == master_rank_) {
    if (logfile_ == NULL) {
      puts("Log file was not opened!");
    } else {
      float avg_gb = static_cast<float>(avg_pag)/num_ranks_ * psize_in_gb_;
      float min_gb = static_cast<float>(min_pag) * psize_in_gb_;
      float max_gb = static_cast<float>(max_pag) * psize_in_gb_;
      fprintf(logfile_, "%s: Max = %f, Min = %f, Average = %f \n",
          info, max_gb, min_gb, avg_gb);
      fflush(logfile_);
    }
  }

  return;
}


void MemInfo::PrintAll(FILE* fout) {
  if (fout == NULL) {
    return;
  }
  
  if ( ! class_initialized_) {
    Init(default_fname_);
  }

  float avail_mem, total_mem;
  bool success = GetMemInfo(&avail_mem, &total_mem);

  char buf[kStrLen];
  int len = 0;
  len += snprintf(buf, kStrLen, "Rank %6d on node %s:",
                  my_rank_, hostname_);
  if (success) {
    len += snprintf(buf+len, kStrLen-len, "%6.2fGB free, %5.2f%% of total\n",
                    avail_mem, avail_mem/total_mem*100.0);
  } else {
    len += snprintf(buf+len, kStrLen-len, "could not get memory info");
  }

  if (my_rank_ == 0) {
    StampSysInfo(fout);
    fprintf(fout, "%s", buf);
    for (int i = 1; i < num_ranks_; ++i) {
      MPI_Status stat;
      buf[0] = '\0';
      BL_MPI_REQUIRE( MPI_Recv(buf, kStrLen, MPI_CHAR, i, i, MPI_COMM_WORLD, &stat) );
      fprintf(fout, "%s", buf);
    }
    fflush(fout);
  } else {
    BL_MPI_REQUIRE( MPI_Send(buf, kStrLen, MPI_CHAR, 0, my_rank_, MPI_COMM_WORLD) );
  }

  return;
}


bool MemInfo::GetMemInfo(float* avail_mem, float* total_mem) {
  if ( ! class_initialized_) {
    Init(default_fname_);
  }
  bool success = GetMyPages();

  if (success) {
    *avail_mem = static_cast<float>(avail_pages_) * psize_in_gb_;
    *total_mem = static_cast<float>(total_pages_) * psize_in_gb_;
  } else {
    *avail_mem = *total_mem = -1.0;
  }

  return(success);
}
