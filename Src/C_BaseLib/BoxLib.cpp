//
// $Id: BoxLib.cpp,v 1.30 2002-11-14 18:43:29 car Exp $
//
#include <winstd.H>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <new>

#include <BoxLib.H>
#include <BLVERSION.H>
#include <DistributionMapping.H>
#include <FArrayBox.H>
#include <ParallelDescriptor.H>
#include <ParmParse.H>
#include <Profiler.H>
#include <Utility.H>
#include <WorkQueue.H>

#define bl_str(s)  # s
#define bl_xstr(s) bl_str(s)
//
// The definition of our version string.
//    
// Takes the form:  boxlib version 2.0 built Jun 25 1996 at 14:52:36
//
const char * const version =

"boxlib version "
bl_xstr(BL_VERSION_MAJOR)
"."
bl_xstr(BL_VERSION_MINOR)
" built "
__DATE__
" at "
__TIME__;

#undef bl_str
#undef bl_xstr

//
// This is used by BoxLib::Error(), BoxLib::Abort(), and BoxLib::Assert()
// to ensure that when writing the message to stderr, that no additional
// heap-based memory is allocated.
//

static
void
write_to_stderr_without_buffering (const char* str)
{
    //
    // Flush all buffers.
    //
    fflush(NULL);

    if (str)
    {
        //
        // Add some `!'s and a newline to the string.
        //
        const char * const end = " !!!\n";
        fwrite(str, strlen(str), 1, stderr);
        fwrite(end, strlen(end), 1, stderr);
    }
}

void
BL_this_is_a_dummy_routine_to_force_version_into_executable ()
{
    write_to_stderr_without_buffering(version);    
}

static
void
write_lib_id(const char* msg)
{
    fflush(0);
    const char* const boxlib = "BoxLib::";
    fwrite(boxlib, strlen(boxlib), 1, stderr);
    if ( msg ) 
    {
	fwrite(msg, strlen(msg), 1, stderr);
	fwrite("::", 2, 1, stderr);
    }
}

void
BoxLib::Error (const char* msg)
{
    write_lib_id("Error");
    write_to_stderr_without_buffering(msg);
    ParallelDescriptor::Abort();
}

namespace
{
  const int EOS = -1;

  std::string
  Trim (const std::string& str)
  {
    int n;
    for ( n = str.size(); --n >= 0; )
      {
	if ( str[n] != ' ' ) break;
      }
    std::string result;
    for (int i = 0; i <= n; ++i )
      {
	result += str[i];
      }
    return result;
  }

  std::string
  Fint_2_string (const int* iarr, int nlen)
  {
    std::string res;
    for ( int i = 0; i < nlen && *iarr != EOS; ++i )
      {
	res += *iarr++;
      }
    return Trim(res);
  }
}

#include <BLFort.H>

BL_FORT_PROC_DECL(BL_ERROR_CPP,bl_error_cpp)
  (
   const int istr[], const int* NSTR
   )
{
  std::string res = Fint_2_string(istr, *NSTR);
  BoxLib::Error(res.c_str());
}

BL_FORT_PROC_DECL(BL_WARNING_CPP,bl_warning_cpp)
  (
   const int istr[], const int* NSTR
   )
{
  std::string res = Fint_2_string(istr, *NSTR);
  BoxLib::Warning(res.c_str());
}

BL_FORT_PROC_DECL(BL_ABORT_CPP,bl_abort_cpp)
  (
   const int istr[], const int* NSTR
   )
{
  std::string res = Fint_2_string(istr, *NSTR);
  BoxLib::Abort(res.c_str());
}

void
BoxLib::Abort (const char* msg)
{
    write_lib_id("Abort");
    write_to_stderr_without_buffering(msg);
    ParallelDescriptor::Abort();
}

void
BoxLib::Warning (const char* msg)
{
    if (msg)
    {
        std::cerr << msg << '!' << '\n';
    }
}

void
BoxLib::Assert (const char* EX,
                const char* file,
                int         line)
{
    const int DIMENSION = 1024;

    char buf[DIMENSION+1];

    sprintf(buf,
            "Assertion `%s' failed, file \"%s\", line %d",
            EX,
            file,
            line);
    //
    // Just to be a little safer :-)
    //
    buf[DIMENSION] = 0;

    write_to_stderr_without_buffering(buf);

    ParallelDescriptor::Abort();
}

void
BoxLib::OutOfMemory (const char* file,
                     int         line)
{
    BoxLib::Assert("operator new", file, line);
}

namespace
{
    Profiler* bl_prf;
}

void
BoxLib::Initialize (int& argc, char**& argv)
{
    static Profiler::Tag bl_prf_tag("BoxLib");

#ifndef WIN32
    //
    // Make sure to catch new failures.
    //
    std::set_new_handler(BoxLib::OutOfMemory);
#endif

    bl_prf = new Profiler(bl_prf_tag, true);
    bl_prf->start();

    ParallelDescriptor::StartParallel(&argc, &argv);

    if (argc == 1)
    {
        ParmParse::Initialize(0,0,0);
    }
    else
    {
        if (strchr(argv[1],'='))
        {
            ParmParse::Initialize(argc-1,argv+1,0);
        }
        else
        {
            ParmParse::Initialize(argc-2,argv+2,argv[1]);
        }
    }

    Profiler::Initialize();
    WorkQueue::Initialize();
    //
    // Initialize random seed after we're running in parallel.
    //
    BoxLib::InitRandom(ParallelDescriptor::MyProc() + 1);

    FArrayBox::Initialize();

    DistributionMapping::Initialize();

    std::cout << std::setprecision(10);
}

void
BoxLib::Finalize ()
{
    bl_prf->stop();
    delete bl_prf;
    DistributionMapping::Finalize();
    FArrayBox::Finalize();
    WorkQueue::Finalize();
    Profiler::Finalize();
    ParmParse::Finalize();
    ParallelDescriptor::EndParallel();
}
