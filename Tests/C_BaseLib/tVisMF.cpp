//BL_COPYRIGHT_NOTICE

//
// $Id: tVisMF.cpp,v 1.16 1998-04-01 17:01:59 car Exp $
//

#include <stdlib.h>

#include <VisMF.H>
#include <Utility.H>

//
// This is HCAll/preload.cpp -- too bad this couldn't be in a library.
//

#ifdef BL_USE_BSP
//
// $Id: tVisMF.cpp,v 1.16 1998-04-01 17:01:59 car Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cstring>
#else
#include <string.h>
#endif

#include "ParallelDescriptor.H"

extern int _bsp_do_stat;
extern int _bsp_do_cgprof;
extern int _bsp_do_prof;
extern int _bsp_nbuffers;
extern int _bsp_buffer_size;
extern int _bsp_slotsize_usecs;
extern int _bsp_buffer_stalls;
extern int _bsp_throttle_procs;
extern int _bsp_comm_fifo_size;
extern int _bsp_opt_contention_level;
extern int _bsp_opt_fcombine_puts;
extern int _bsp_opt_fcombine_puts_max;
extern int _bsp_opt_fcombine_puts_min;
extern int _bsp_opt_bsmp_buffer_size;
extern char *_bsp_compile_flags;
extern char *_bsp_arch;
extern char *_bsp_device;
extern char *_bsp_include_dir;
extern int _bsp_check_syncs;
extern char *_bsp_exec_file;
extern char _bsp_library_type;
extern int  _bsp_opt_flibrary_level;


extern "C" void _bsp_preload_init ();

//
// Set _bsp_include_dir from the environment else take precompiled default.
//

static
void
get_bsp_include_dir ()
{
    const char* dir = getenv("_bsp_include_dir");

    if (dir == 0 || *dir == 0)
    {
        if (_bsp_include_dir == 0)
            bsp_abort("_bsp_include_dir must be set");
    }
    else
    {
        if (!(_bsp_include_dir == 0))
            free(_bsp_include_dir);

        if ((_bsp_include_dir = (char*) malloc(strlen(dir) + 1)) == 0)
            bsp_abort("malloc() failed");

        strcpy(_bsp_include_dir, dir);
    }

    printf("Using _bsp_include_dir=%s\n", _bsp_include_dir);

    fflush(stdout);
}

//
// This function is written by BSP when configuring BSP.
// BSP expects to be able to call this function on startup.
//

#ifdef BL_T3E

void
_bsp_preload_init()
{
   _bsp_do_cgprof        = 0;
   _bsp_do_prof          = 0;
   _bsp_do_stat          = 0;
   _bsp_nbuffers         = 1;
   _bsp_buffer_size      = 1048576;
   _bsp_slotsize_usecs   = 0;
   _bsp_throttle_procs   = 0;
   _bsp_comm_fifo_size   = 100;
   _bsp_buffer_stalls    = 2;
   _bsp_opt_contention_level = 1;
   _bsp_opt_fcombine_puts= 0;
   _bsp_opt_fcombine_puts_max=10485760;
   _bsp_opt_fcombine_puts_min=10240;
   _bsp_opt_bsmp_buffer_size =-1;
   _bsp_check_syncs  =0;
   _bsp_library_type ='O';
   _bsp_opt_flibrary_level=2;

   _bsp_compile_flags  = (char*) malloc(1+strlen(" -O3 -flibrary-level 2 -fcontention-resolve 1"));
   _bsp_arch=(char*) malloc(1+strlen("CRAYT3E"));
   _bsp_device=(char*) malloc(1+strlen("DRMA_SHMEM"));
#ifdef BL_T3E_NERSC
   _bsp_include_dir=(char*) malloc(1+strlen("/u1/vince/BSP/include/"));
   _bsp_include_dir  =strcpy(_bsp_include_dir,"/u1/vince/BSP/include/");
#endif
#ifdef BL_T3E_NAVO
   _bsp_include_dir=(char*) malloc(1+strlen("/home/Cvince/BSP/include/"));
   _bsp_include_dir  =strcpy(_bsp_include_dir,"/home/Cvince/BSP/include/");
#endif
   _bsp_exec_file= (char*)malloc(1+strlen("main"));
   if (_bsp_compile_flags==NULL || _bsp_arch==NULL ||
       _bsp_include_dir==NULL || _bsp_exec_file==NULL)
     bsp_abort("{bsp_begin}: unable to malloc for compile flags");

   _bsp_compile_flags=strcpy(_bsp_compile_flags, " -O3 -flibrary-level 2 -fcontention-resolve 1");
   _bsp_arch         =strcpy(_bsp_arch,"CRAYT3E");
   _bsp_device       =strcpy(_bsp_device,"DRMA_SHMEM");
   _bsp_exec_file    =strcpy(_bsp_exec_file,"main");

    //
    // This call is not part of the original BSP code.
    // This allows us to override where BSP_INCLUDE_DIR is found.
    //
    get_bsp_include_dir();
}

#else   /* not BL_T3E */

void
_bsp_preload_init()
{
   _bsp_do_cgprof        = 0;
   _bsp_do_prof          = 0;
   _bsp_do_stat          = 0;
   _bsp_nbuffers         = 2;
   _bsp_buffer_size      = 10240;
   _bsp_slotsize_usecs   = 0;
   _bsp_throttle_procs   = 0;
   _bsp_comm_fifo_size   = 100;
   _bsp_buffer_stalls    = 2;
   _bsp_opt_contention_level = 1;
   _bsp_opt_fcombine_puts= 20480;
   _bsp_opt_fcombine_puts_max=102400;
   _bsp_opt_fcombine_puts_min=5120;
   _bsp_opt_bsmp_buffer_size =-1;
   _bsp_check_syncs  =1;
   _bsp_library_type ='O';
   _bsp_opt_flibrary_level=2;

   _bsp_compile_flags  = (char*) malloc(1+strlen(" -O3 -flibrary-level 2 -fcombine-puts-buffer 20480,102400,5120 -fcontention-resolve 1"));
   _bsp_arch=(char*) malloc(1+strlen("OSF1"));
   _bsp_device=(char*) malloc(1+strlen("SHMEM_SYSV"));
   _bsp_include_dir=(char*) malloc(1+strlen("/usr/people/vince/Parallel/BSP/BSP/include/"));
   _bsp_exec_file= (char*)malloc(1+strlen("hedgehog"));
   if (_bsp_compile_flags==NULL || _bsp_arch==NULL ||
       _bsp_include_dir==NULL || _bsp_exec_file==NULL)
     bsp_abort("{bsp_begin}: unable to malloc for compile flags");

   _bsp_compile_flags=strcpy(_bsp_compile_flags, " -O3 -flibrary-level 2 -fcombine-puts-buffer 20480,102400,5120 -fcontention-resolve 1");
   _bsp_arch         =strcpy(_bsp_arch,"OSF1");
   _bsp_device       =strcpy(_bsp_device,"SHMEM_SYSV");
   _bsp_include_dir  =strcpy(_bsp_include_dir,"/usr/people/vince/Parallel/BSP/BSP/include/");
   _bsp_exec_file    =strcpy(_bsp_exec_file,"hedgehog");

    //
    // This call is not part of the original BSP code.
    // This allows us to override where _bsp_include_dir is found.
    //
    get_bsp_include_dir();
}

#endif

#endif /*BL_USE_BSP*/



static int nProcs = 1;

static int nBoxs  = 10;

static char* the_prog_name;

static aString PerFab("PerFab");

static aString PerCPU("PerCPU");

//
// How defaults to PerCPU.
//
static aString How(PerCPU);

static
void
usage ()
{
    std::cout << "usage: "
              << the_prog_name
              << " [-how PerFab|PerCPU]"
              << " [-nprocs N]"
              << " [-nboxs N]"
              << std::endl;
    exit(1);
}

static
void
parse_args (char**& argv)
{
    while (*++argv && **argv == '-')
    {
        if (strcmp(*argv, "-nprocs") ==  0)
        {
            if (*++argv)
            {
                nProcs = atoi(*argv);

                if (nProcs <= 0)
                {
                    std::cout << "nprocs must be positive" << std::endl;
                    usage();
                }
            }
            else
            {
                std::cout << "No argument to -nprocs supplied.\n";
                usage();
            }
        }
        else if (strcmp(*argv, "-nboxs") ==  0)
        {
            if (*++argv)
            {
                nBoxs = atoi(*argv);

                if (nBoxs <= 0)
                {
                    std::cout << "nboxs must be positive" << std::endl;
                    usage();
                }
            }
            else
            {
                std::cout << "No argument to -nprocs supplied.\n";
                usage();
            }
        }
        else if (strcmp(*argv, "-how") ==  0)
        {
            if (*++argv)
            {
                How = *argv;

                if (!(How == PerCPU || How == PerFab))
                {
                    std::cout << "Invalid value for -how argument\n";
                    usage();
                }
            }
            else
            {
                std::cout << "No argument to -how supplied.\n";
                usage();
            }
        }
        else
        {
            std::cout << "Exiting, unknown option: " << *argv << std::endl;
            usage();
        }
    }
}

static
void
Write_N_Read (const MultiFab& mf,
              const aString&  mf_name,
              VisMF::How      how)
{
    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "Writing the MultiFab to disk ...\n";
    }

    double start, end;

    ParallelDescriptor::Synchronize();

    if (ParallelDescriptor::IOProcessor())
    {
        start = Utility::wsecond();
    }

    switch (how)
    {
    case VisMF::OneFilePerCPU:
        VisMF::Write(mf, mf_name, VisMF::OneFilePerCPU); break;
    case VisMF::OneFilePerFab:
        VisMF::Write(mf, mf_name, VisMF::OneFilePerFab); break;
    default:
        BoxLib::Error("Bad case in switch");
    }

    ParallelDescriptor::Synchronize();

    if (ParallelDescriptor::IOProcessor())
    {
        end = Utility::wsecond();

        std::cout << "\nWallclock time for MF write: " << (end-start) << '\n';

        std::cout << "Reading the MultiFab from disk ...\n";
    }

    VisMF vmf(mf_name);

    assert(vmf.length() == mf.boxArray().length());

    for (ConstMultiFabIterator mfi(mf); mfi.isValid(); ++mfi)
    {
        //const FArrayBox& fab = vmf[mfi.index()];
        const FArrayBox& fab = vmf.GetFab(mfi.index(), 0);

        std::cout << "\tCPU #"
                  << ParallelDescriptor::MyProc()
                  << " read FAB #"
                  << mfi.index()
                  << '\n';
    }

    ParallelDescriptor::Synchronize();

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "Building new MultiFab from disk version ....\n\n";
    }

    MultiFab new_mf;
    
    VisMF::Read(new_mf, mf_name);
}

int
main (int, char** argv)
{
    the_prog_name = argv[0];

    parse_args(argv);

    ParallelDescriptor::StartParallel(nProcs, &argc, &argv);

    BoxArray ba(nBoxs);

    ba.set(0, Box(IntVect(D_DECL(0,0,0)), IntVect(D_DECL(2,2,2))));

    for (int i = 1; i < nBoxs; i++)
    {
        ba.set(i,grow(ba[i-1],2));
    }

    MultiFab mf(ba, 2, 1);

    for (ConstMultiFabIterator mfi(mf); mfi.isValid(); ++mfi)
    {
        mf[mfi.index()].setVal(mfi.index()+1);
    }
    //
    // Set cells in ghost region to zero.
    //
    mf.setBndry(0);

    static const aString mf_name = "Spam-n-Eggs";

    Write_N_Read (mf,
                  mf_name,
                  (How==PerCPU) ? VisMF::OneFilePerCPU : VisMF::OneFilePerFab);

    ParallelDescriptor::EndParallel();
}
