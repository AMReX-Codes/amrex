!IF	DEFINED(DEPENDS) || ! EXIST(f90.depends)
f90.depends: $(f90sources) $(fsources)
	-mkdir $(mdir)
	perl $(MODDEP) --odir $(obj_dir) --objext obj $(f90sources) $(fsources) > f90.depends
	nmake main.exe
!ENDIF

$(obj_dir):
	@mkdir $(obj_dir)

clean:
	-del /q main.exe
	-del /q $(obj_dir)\*.mod $(obj_dir)\*.obj
	-del /q $(obj_dir)\*.pdb
	-del /q f90.depends
	-rd /s/q fmod fnmod_ndebug

!IF	!DEFINED(DEPENDS) && EXIST(f90.depends)
!INCLUDE f90.depends
!ENDIF


