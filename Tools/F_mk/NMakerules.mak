$(tdir)\f90.depends: $(f90sources) $(fsources) "$(tdir)"
	perl $(MODDEP) --odir $(obj_dir) $(f90sources) $(fsources) > $(tdir)\f90.depends

$(tdir)\c.depends: $(csources) "$(tdir)"
	perl $(MKDEP) --odir $(obj_dir) $(csources) > $(tdir)\c.depends

"$(tdir)":
	if not exist "$(tdir)\" mkdir "$(tdir)"

"$(obj_dir)": "$(tdir)"
	if not exist "$(obj_dir)\" mkdir "$(obj_dir)"

.f{$(obj_dir)}.obj:
	@if not exist "$(obj_dir)\" mkdir "$(obj_dir)"
	$(FOR) /c $(FFLAGS) /object:$(obj_dir)\ $<

.f90{$(obj_dir)}.obj:
	@if not exist "$(obj_dir)\" mkdir "$(obj_dir)"
	$(FOR) /c $(FFLAGS) /object:$(obj_dir)\ $<

.c{$(obj_dir)}.obj:
	@if not exist "$(obj_dir)\" mkdir "$(obj_dir)"
	$(CC) /c $(CFLAGS) /Fo$(obj_dir)\ $<

clean:
	-del /q main.exe main.pdb df60.pdb
	-rd /s/q $(tdir)

depend:	$(tdir)\f90.depends $(tdir)\c.depends

!IF	EXIST($(tdir)\f90.depends)
!INCLUDE $(tdir)\f90.depends
!ENDIF
!IF	EXIST($(tdir)\c.depends)
!INCLUDE $(tdir)\c.depends
!ENDIF
