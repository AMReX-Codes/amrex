#ifndef NODE_INFO
#define NODE_INFO

# include <cstdlib>
# include <iostream>
# include <sstream>
# include <string>
# include <unordered_map>

// Implementation provided machine primitives.

  
  struct NodeHardware {
    int thread_per_core;
    int core_per_numa;
    int numa_per_node;
    
    // cpuid indexing dimension order
    int thread_stride;
    int core_stride;
    int numa_stride;
    
    int thread_per_numa() const { return thread_per_core*core_per_numa; }
    int thread_per_node() const { return thread_per_core*core_per_numa*numa_per_node; }
    int core_per_node() const { return core_per_numa*numa_per_node; }
    
    int numa_of_core(int core) const { return core/core_per_numa; }
    int numa_of_thread(int thread) const { return thread/(thread_per_core*core_per_numa); }
    
    int core_of_thread(int thread) const { return thread/thread_per_core; }
    
    int cpuid_of(int thread, int core, int numa) const {
      return thread*thread_stride + core*core_stride + numa*numa_stride;
    }
  };


   NodeHardware query_node_hardware();
#endif
