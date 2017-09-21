//This code was developed in Rambutan to parse machine information

#include "sysInfo.H"
#include <fstream>
#include <algorithm>
#include <vector>
#include <cassert>
using namespace std;


namespace {
  bool read_line_prefix(istream &is, string &prefix) {
    int spaces = 0;
    
    prefix.clear();
    
    while(true) {
      int c = is.get();

      if(c == EOF)
        return false;
      
      if(c == '\n') {
        prefix.clear();
        spaces = 0;
        continue;
      }
      
      if(c == ':')
        return true;
      
      if(c == ' ')
        spaces += 1;
      else if(c == '\t')
        ;
      else {
        while(spaces > 0) {
          prefix.push_back(' ');
          spaces -= 1;
        }
        prefix.push_back((char)c);
      }
    }
  }
}

NodeHardware query_node_hardware() {
  NodeHardware ans;
  
  struct Entry {
    int proc;
    int phys;
    int core;
    int sib_n;
    int core_n;
  };
  
  vector<Entry> procs;
  
  { // parse /proc/cpuinfo
    ifstream f("/proc/cpuinfo", ios::in);
    
    string pre;
    
    Entry e;
    int n = 5;
    
    while(read_line_prefix(f, pre)) {
      if(pre == "processor") {
        assert(n == 5);
        f >> e.proc;
      }
      else if(pre == "physical id")
        f >> e.phys;
      else if(pre == "siblings")
        f >> e.sib_n;
      else if(pre == "core id")
        f >> e.core;
      else if(pre == "cpu cores")
        f >> e.core_n;
      else
        continue;
      
      n -= 1;
      if(n == 0) {
        procs.push_back(e);
        n = 5;
      }
    }
  }
  
  // sort procs by "proc id", weird if wasnt already sorted by OS
  std::sort(procs.begin(), procs.end(), [](Entry a, Entry b) { return a.proc < b.proc; });
  
  // assert OS numbers procs contiguously from zero
  bool ok = true;
  for(int i=0; i < (int)procs.size(); i++)
    ok = ok && procs[i].proc == i;
  assert(ok);
  
  ans.numa_per_node = procs.size() / procs[0].sib_n;
  ans.core_per_numa = procs[0].core_n;
  ans.thread_per_core = procs[0].sib_n / procs[0].core_n;
  
  if(procs.size() == 1) {
    ans.thread_stride = 0;
    ans.core_stride = 0;
    ans.numa_stride = 0;
  }
  else {
    ans.thread_stride = procs.size();
    for(int i=1; i < (int)procs.size(); i++) {
      if(procs[i].core == procs[0].core && procs[i].phys == procs[0].phys) {
        ans.thread_stride = i;
        break;
      }
    }
    
    ans.core_stride = procs.size();
    for(int i=1; i < (int)procs.size(); i++) {
      if(procs[i].core != procs[0].core) {
        ans.core_stride = i;
        break;
      }
    }
    
    ans.numa_stride = procs.size();
    for(int i=1; i < (int)procs.size(); i++) {
      if(procs[i].phys != procs[0].phys) {
        ans.numa_stride = i;
        break;
      }
    }
  }
  
  return ans;
}

