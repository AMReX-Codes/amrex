TOP = ..\..

#f90sources = $(f90sources) main.f90

!INCLUDE $(TOP)\mk\NMakedefs.mak

boxlib_dir=..
!INCLUDE ..\NPackage.mak

f90objects = $(f90objects) $(obj_dir)\main.obj
f90objects = $(f90objects) $(obj_dir)\t_main.obj
f90objects = $(f90objects) $(obj_dir)\t_knapsack.obj

main.exe: $(objects)
	$(FOR) $(FFLAGS) $(LDFLAGS) /out:main.exe $(objects)

!INCLUDE $(TOP)\mk\NMakerules.mak
