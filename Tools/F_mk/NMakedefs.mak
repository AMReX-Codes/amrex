#MPI=T
#NDEBUG=T
#OMP=T
#ALL=T

FOR=f90
F77=f77

sources = $(f90sources) $(fsources) $(csources)
objects = $(f90objects) $(cobjects) $(fobjects)

!IFDEF MPI
mpi_home = C:\Program Files\Argonne National Lab\MPICH.NT.1.2.4\SDK
mpi_incl = /I"$(mpi_home)\Include"
!IFDEF NDEBUG
mpi_ldfl = /LIBPATH:"$(mpi_home)\LIB" ws2_32.lib mpich.lib
!ELSE
mpi_ldfl = /LIBPATH:"$(mpi_home)\LIB" ws2_32.lib mpichd.lib
!ENDIF
!ENDIF

tdir    = t

!IFDEF NDEBUG
mdir    = $(tdir)\obj
obj_dir = $(tdir)\obj
CFLAGS  = /nologo /MT /DWIN32
FFLAGS  = /nologo /iface:cref /fast /MT /tune:host /inline:all /optimize:5 $(mpi_incl)
FFLAGS  = /nologo /iface:cref /fast /MT /tune:host /optimize:5 $(mpi_incl)
FFLAGS  = $(FFLAGS) /module:$(mdir)
!ELSE
mdir    = $(tdir)\debug.obj
obj_dir = $(tdir)\debug.obj
CFLAGS  = /nologo /MTd /DWIN32
FFLAGS  = /nologo /iface:cref /debug /traceback /check:all /MTd $(mpi_incl)
FFLAGS  = $(FFLAGS) /module:$(mdir)
FFLAGS  = $(FFLAGS) /warn:uninitialized
FFLAGS  = $(FFLAGS) /warn:unused
FFLAGS  = $(FFLAGS) /warn:truncated_source
FFLAGS  = $(FFLAGS) /warn:uncalled
FFLAGS  = $(FFLAGS) /check:nounderflow
!ENDIF

LDFLAGS = /link $(mpi_ldfl) /stack:8000000

MODDEP=$(TOP)\scripts\moddep.pl
MKDEP=$(TOP)\scripts\mkdep.pl
