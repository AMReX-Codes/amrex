TOP = ..\..

!INCLUDE $(TOP)\mk\NMakedefs.mak

boxlib_dir=$(TOP)\boxlib
!INCLUDE $(boxlib_dir)\NPackage.mak
mg_dir=$(TOP)\mg
!INCLUDE $(mg_dir)\NPackage.mak
sparskit_dir=$(TOP)\extern\SPARSKIT
!INCLUDE $(sparskit_dir)\NPackage.mak
lapack_dir=$(TOP)\extern\LAPACK
!INCLUDE $(lapack_dir)\NPackage.mak

f90sources = $(f90sources) main.f90
f90objects = $(f90objects) main.obj

f90sources = $(f90sources) wrapper.f90
f90objects = $(f90objects) wrapper.obj

f90sources = $(f90sources) cc_multi.f90
f90objects = $(f90objects) cc_multi.obj

f90sources = $(f90sources) cc_single.f90
f90objects = $(f90objects) cc_single.obj

f90sources = $(f90sources) nodal_single.f90
f90objects = $(f90objects) nodal_single.obj

f90sources = $(f90sources) t_smoother.f90
f90objects = $(f90objects) t_smoother.obj

f90sources = $(f90sources) nodal_multi.f90
f90objects = $(f90objects) nodal_multi.obj

main.exe: $(objects)
	$(FOR) $(FFLAGS) $(LDFLAGS) /out:main.exe $(objects)

!INCLUDE $(TOP)\mk\NMakerules.mak
