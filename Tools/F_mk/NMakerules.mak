!IF	DEFINED(DEPENDS) || ! EXIST(f90.depends)
f90.depends: $(sources)
	-mkdir $(mdir)
	perl moddep.pl --objext obj $(sources) > f90.depends
	nmake main.exe
!ENDIF

clean:
	-del /q main.exe
	-del /q *.mod *.obj
	-del /q *.pdb
	-del /q f90.depends
	-rd /s/q fmod fnmod_ndebug

!IF	!DEFINED(DEPENDS) && EXIST(f90.depends)
!INCLUDE f90.depends
!ENDIF
