#ifndef P_PERILLACONFIG_H
#define P_PERILLACONFIG_H

//class PerillaConfig
namespace perilla
{
  static const int NUM_THREAD_TEAMS=1;
  static const int NUM_THREADS_PER_TEAM=6;
  static const int NUM_COMM_THREADS=1;
  static const int NUM_PREGENERATED_PACKAGES=2;
  static const int MSG_QUEUE_DEFAULT_MAXSIZE=2; // same as num of pregenerated packages because these got swaped between pQ and recycleQ
  static const int TASK_QUEUE_DEFAULT_MAXSIZE=32;
  static const int MAX_SQRT_TAG=512;  //sqrt(512*512)
  static const bool NUMA_AWARE_MESSAGEHANDLER=true;
  static const int LOCK_GRANULARITY=1; //!0 is queue level, 1 is region level
  
  static const bool LAZY_PUSH=false;
  static const bool LAZY_PULL=false;
  static const bool PACKING_FINEGRAIN=false; //!PACKING_FINEGRAIN = TRUE is effective when LAZY_PUSH=FALSE
  static const bool UNPACKING_FINEGRAIN=false; //!is effective when LAZY_PUSH=FALSE

};
#endif
