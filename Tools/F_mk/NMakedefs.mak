#MPI=T
#NDEBUG=T
#OMP=T
#ALL=T

FOR=f90
F77=f77

MULTIGRID = t
EIGEN     = t
VARDEN    =

msources =		\
	main.f90	\
	t_main.f90

msources = $(msources)	\
	BoxLib.f90	\
	parallel_stubs.f90	\
	bl_constants.f90 \
	bl_error.f90	\
	bl_space.f90	\
	bl_types.f90	\
	bl_utility.f90	\
	box.f90		\
	boxarray.f90	\
	fab.f90		\
	layout.f90	\
	multifab.f90	\
	bndry_reg.f90   \
	mboxarray.f90	\
	box_util.f90	\
	sort_box.f90	\
	omp.f90		\
	knapsack.f90	\
	vector_i.f90	\
	sort_d.f90	\
	plotfile.f90	\
	fabio.f90	\
	list_box.f90	\
	cluster.f90	\
	mt19937ar.f90	\
	omp_stubs.f90	\
	f2kcli_win32.f90	\
	ppm_util.f90

!IFDEF	EIGEN
msources = $(msources)	\
	cc_eigen.f90    \
	t_eigen.f90
!ENDIF

!IFDEF MULTIGRID
msources = $(msources)	\
	ml.f90		\
	mlmg.f90	\
	mg.f90		\
	itsol.f90       \
	mg_smoother.f90 \
	stencil.f90	\
	stencil_nodal.f90	\
	sparse_solve.f90	\
	wrapper.f90	\
	cc_multi.f90	\
	cc_single.f90	\
	nodal_multi.f90
!ENDIF

!IFDEF VARDEN
msources = $(msources)	\
	varden.f90	\
	macproject.f90	\
	hgproject.f90	\
	advance.f90	\
	estdt.f90	\
	initdata.f90	\
	mkflux.f90	\
	slope.f90	\
	bc.f90		\
	mkutrans.f90	\
	cvmg.f90	\
	setvelbc.f90
!ENDIF

!IFDEF	MULTIGRID
ssources = sk_sup.f  iters.f ilut.f
!ENDIF

!IFDEF	EIGEN
ssources = $(ssources) dsygv.f
!ENDIF

csources = f2kgetcl.c \
	timer_c.c \
	fabio_c.c \
	ppm_util_c.c

sources = $(msources) $(csources)

mobjects=$(msources:.f90=.obj)
cobjects=$(csources:.c=.obj)
sobjects=$(ssources:.f=.obj)

objects=$(mobjects:.f=.obj) $(cobjects) $(sobjects)

!IFDEF MPI
mpi_home = C:\Program Files\Argonne National Lab\MPICH.NT.1.2.4\SDK
mpi_incl = /I"$(mpi_home)\Include"
!IFDEF NDEBUG
mpi_ldfl = /LIBPATH:"$(mpi_home)\LIB" ws2_32.lib mpich.lib
!ELSE
mpi_ldfl = /LIBPATH:"$(mpi_home)\LIB" ws2_32.lib mpichd.lib
!ENDIF
!ENDIF

!IFDEF NDEBUG
mdir    = fmod_ndebug
CFLAGS  = /nologo /MT /DWIN32
FFLAGS  = /nologo /iface:cref /fast /MT /tune:host /inline:all /optimize:5 $(mpi_incl)
FFLAGS  = /nologo /iface:cref /fast /MT /tune:host /optimize:5 $(mpi_incl)
FFLAGS  = $(FFLAGS) /module:$(mdir)
!ELSE
mdir    = fmod
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

!IF	DEFINED(DEPENDS) || ! EXIST(f90.depends)
f90.depends: $(sources)
	-mkdir $(mdir)
	perl moddep.pl --objext obj $(sources) > f90.depends
	nmake main.exe
!ENDIF

main.exe: $(objects) 
	$(FOR) $(FFLAGS) $(LDFLAGS) /out:main.exe $(objects)

f2kgetcl.obj: f2kgetcl.c
	$(CC) /DUPPER $(CFLAGS) /c f2kgetcl.c

ilut.obj: SPARSKIT\ilut.f
	$(F77) -c $(FFLAGS) SPARSKIT\ilut.f

iters.obj: SPARSKIT\iters.f
	$(F77) -c $(FFLAGS) SPARSKIT\iters.f

sk_sup.obj: SPARSKIT\sk_sup.f
	$(F77) -c $(FFLAGS) SPARSKIT\sk_sup.f

dsygv.obj: dsygv.f
	$(F77) -c $(FFLAGS) /optimize:3 /check:nounderflow dsygv.f

clean:
	-del /q main.exe
	-del /q *.mod *.obj
	-del /q *.pdb
	-del /q f90.depends
	-rd /s/q fmod fnmod_ndebug

!IF	!DEFINED(DEPENDS) && EXIST(f90.depends)
!INCLUDE f90.depends
!ENDIF
