
#ifndef _BL_PROFILER_H_
#define _BL_PROFILER_H_

class Profiler 
{
public:
  class Tag
  {
  public:
    explicit Tag(const char* str) {}
  };
  Profiler(const Tag&, bool) {}
  static void Initialize() {}
  static void Finalize() {}
  void start() {}
  void stop() {}
};

#define BL_PROFILE(a)
#define BL_PROFILE_TIMER(var, a)
#define BL_PROFILE_START(var)
#define BL_PROFILE_STOP(var)
#define BL_PROFILE_THIS_NAME()

#endif
