#ifndef _BL_PROFILER_H_
#define _BL_PROFILER_H_

#include <iostream>
#include <string>
#include <typeinfo>

#ifdef BL_NAMESPACE
namespace BoxLib
#endif

class Profiler
{
public:
    class Tag;
    Profiler(const Tag& tag_, bool hold = false);
    void start();
    void stop();
    ~Profiler();
    static void on();
    static void off();
    static bool is_profiling();
    static std::string clean_name(const std::string& str);
    static void Initialize(int& argc, char**& argv);
    static void Finalize();
private:
    static void glean();
    const Tag& tag;
    bool started;
    static bool profiling;
    static bool initialized;
    void* operator new(size_t);
    Profiler(const Profiler& prf);
    Profiler& operator=(const Profiler& prf);
};

class Profiler::Tag
{
public:
    explicit Tag(const std::string& tag_);
    const std::string& name() const;
private:
    int itag;
    const std::string tag;
    static int next_itag;
};

#ifdef BL_PROFILING
#define BL_PROFILE(a)							\
static Profiler::Tag bl_prf_tag_((a));			\
Profiler bl_prf_(bl_prf_tag_, false)

#define BL_PROFILE_TIMER(var, a)					\
static Profiler::Tag bl_prf_tag_##var((a));			\
Profiler bl_prf_##var(bl_prf_tag_##var, true)

#define BL_PROFILE_START(var) bl_prf_##var.start()
#define BL_PROFILE_STOP(var) bl_prf_##var.stop()
#define BL_PROFILE_THIS_NAME() BoxLib::Profiler::clean_name(std::string(typeid(*this).name()))
#else
#define BL_PROFILE(a)
#define BL_PROFILE_TIMER(var, a)
#define BL_PROFILE_START(var)
#define BL_PROFILE_STOP(var)
#define BL_PROFILE_THIS_NAME()
#endif

#ifdef BL_NAMESPACE
}
#endif

#endif
