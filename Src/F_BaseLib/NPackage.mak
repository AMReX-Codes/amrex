f90sources = $(f90sources) $(boxlib_dir)\BoxLib.f90
f90objects = $(f90objects) $(obj_dir)\BoxLib.obj

f90sources = $(f90sources) $(boxlib_dir)\omp.f90
f90objects = $(f90objects) $(obj_dir)\omp.obj

f90sources = $(f90sources) $(boxlib_dir)\f2kcli_win32.f90
f90objects = $(f90objects) $(obj_dir)\f2kcli_win32.obj

f90sources = $(f90sources) $(boxlib_dir)\mt19937ar.f90
f90objects = $(f90objects) $(obj_dir)\mt19937ar.obj

f90sources = $(f90sources) $(boxlib_dir)\bl_constants.f90
f90objects = $(f90objects) $(obj_dir)\bl_constants.obj

f90sources = $(f90sources) $(boxlib_dir)\bl_error.f90
f90objects = $(f90objects) $(obj_dir)\bl_error.obj

f90sources = $(f90sources) $(boxlib_dir)\bl_IO.f90
f90objects = $(f90objects) $(obj_dir)\bl_IO.obj

f90sources = $(f90sources) $(boxlib_dir)\bl_mem_stat.f90
f90objects = $(f90objects) $(obj_dir)\bl_mem_stat.obj

f90sources = $(f90sources) $(boxlib_dir)\bl_parmparse.f90
f90objects = $(f90objects) $(obj_dir)\bl_parmparse.obj

f90sources = $(f90sources) $(boxlib_dir)\bl_space.f90
f90objects = $(f90objects) $(obj_dir)\bl_space.obj

f90sources = $(f90sources) $(boxlib_dir)\bl_stream.f90
f90objects = $(f90objects) $(obj_dir)\bl_stream.obj

f90sources = $(f90sources) $(boxlib_dir)\bl_string.f90
f90objects = $(f90objects) $(obj_dir)\bl_string.obj

f90sources = $(f90sources) $(boxlib_dir)\bl_timer.f90
f90objects = $(f90objects) $(obj_dir)\bl_timer.obj

f90sources = $(f90sources) $(boxlib_dir)\bl_types.f90
f90objects = $(f90objects) $(obj_dir)\bl_types.obj

f90sources = $(f90sources) $(boxlib_dir)\box.f90
f90objects = $(f90objects) $(obj_dir)\box.obj

f90sources = $(f90sources) $(boxlib_dir)\knapsack.f90
f90objects = $(f90objects) $(obj_dir)\knapsack.obj

f90sources = $(f90sources) $(boxlib_dir)\layout.f90
f90objects = $(f90objects) $(obj_dir)\layout.obj

f90sources = $(f90sources) $(boxlib_dir)\boxarray.f90
f90objects = $(f90objects) $(obj_dir)\boxarray.obj

f90sources = $(f90sources) $(boxlib_dir)\mboxarray.f90
f90objects = $(f90objects) $(obj_dir)\mboxarray.obj

f90sources = $(f90sources) $(boxlib_dir)\box_util.f90
f90objects = $(f90objects) $(obj_dir)\box_util.obj

f90sources = $(f90sources) $(boxlib_dir)\fab.f90
f90objects = $(f90objects) $(obj_dir)\fab.obj

f90sources = $(f90sources) $(boxlib_dir)\multifab.f90
f90objects = $(f90objects) $(obj_dir)\multifab.obj

f90sources = $(f90sources) $(boxlib_dir)\fabio.f90
f90objects = $(f90objects) $(obj_dir)\fabio.obj

f90sources = $(f90sources) $(boxlib_dir)\plotfile.f90
f90objects = $(f90objects) $(obj_dir)\plotfile.obj

f90sources = $(f90sources) $(boxlib_dir)\filler.f90
f90objects = $(f90objects) $(obj_dir)\filler.obj

f90sources = $(f90sources) $(boxlib_dir)\cluster.f90
f90objects = $(f90objects) $(obj_dir)\cluster.obj

f90sources = $(f90sources) $(boxlib_dir)\interp.f90
f90objects = $(f90objects) $(obj_dir)\interp.obj

f90sources = $(f90sources) $(boxlib_dir)\bndry_reg.f90
f90objects = $(f90objects) $(obj_dir)\bndry_reg.obj

f90sources = $(f90sources) $(boxlib_dir)\list_box.f90
f90objects = $(f90objects) $(obj_dir)\list_box.obj

f90sources = $(f90sources) $(boxlib_dir)\sort_box.f90
f90objects = $(f90objects) $(obj_dir)\sort_box.obj

f90sources = $(f90sources) $(boxlib_dir)\vector_i.f90
f90objects = $(f90objects) $(obj_dir)\vector_i.obj

f90sources = $(f90sources) $(boxlib_dir)\sort_d.f90
f90objects = $(f90objects) $(obj_dir)\sort_d.obj

f90sources = $(f90sources) $(boxlib_dir)\ppm_util.f90
f90objects = $(f90objects) $(obj_dir)\ppm_util.obj

f90sources = $(f90sources) $(boxlib_dir)\parallel_stubs.f90
f90objects = $(f90objects) $(obj_dir)\parallel_stubs.obj

f90sources = $(f90sources) $(boxlib_dir)\omp_stubs.f90
f90objects = $(f90objects) $(obj_dir)\omp_stubs.obj

csources = $(csources) $(boxlib_dir)\fabio_c.c
cobjects = $(cobjects) $(obj_dir)\fabio_c.obj

csources = $(csources) $(boxlib_dir)\timer_c.c
cobjects = $(cobjects) $(obj_dir)\timer_c.obj

csources = $(csources) $(boxlib_dir)\ppm_util_c.c
cobjects = $(cobjects) $(obj_dir)\ppm_util_c.obj

csources = $(csources) $(boxlib_dir)\f2kgetcl.c
cobjects = $(cobjects) $(obj_dir)\f2kgetcl.obj

{$(boxlib_dir)}.f{$(obj_dir)}.obj:
	@if not exist "$(obj_dir)\" mkdir "$(obj_dir)"
	$(FOR) /c $(FFLAGS) $< $(FOB)

{$(boxlib_dir)}.f90{$(obj_dir)}.obj:
	@if not exist "$(obj_dir)\" mkdir "$(obj_dir)"
	$(FOR) /c $(FFLAGS) $< $(FOB)

{$(boxlib_dir)}.c{$(obj_dir)}.obj:
	@if not exist "$(obj_dir)\" mkdir "$(obj_dir)"
	$(CC) /c $(CFLAGS) /Fo$(obj_dir)/ $<
