f90sources = $(f90sources) BoxLib.f90
f90sources = $(f90sources) omp.f90
f90sources = $(f90sources) f2kcli$_win32.f90
f90sources = $(f90sources) mt19937ar.f90

f90sources = $(f90sources) bl_constants.f90
f90sources = $(f90sources) bl_error.f90  
f90sources = $(f90sources) bl_IO.f90
f90sources = $(f90sources) bl_kiss.f90
f90sources = $(f90sources) bl_mem_stat.f90
f90sources = $(f90sources) bl_space.f90  
f90sources = $(f90sources) bl_stream.f90  
f90sources = $(f90sources) bl_string.f90  
f90sources = $(f90sources) bl_timer.f90  
f90sources = $(f90sources) bl_types.f90  

f90sources = $(f90sources) box.f90  
f90sources = $(f90sources) knapsack.f90 
f90sources = $(f90sources) layout.f90  
f90sources = $(f90sources) boxarray.f90  
f90sources = $(f90sources) mboxarray.f90
f90sources = $(f90sources) box_util.f90
f90sources = $(f90sources) fab.f90  
f90sources = $(f90sources) multifab.f90
f90sources = $(f90sources) fabio.f90
f90sources = $(f90sources) plotfile.f90
f90sources = $(f90sources) filler.f90
f90sources = $(f90sources) cluster.f90
f90sources = $(f90sources) interp.f90
f90sources = $(f90sources) bndry_reg.f90

f90sources = $(f90sources) list_box.f90
f90sources = $(f90sources) sort_box.f90
f90sources = $(f90sources) vector_i.f90 
f90sources = $(f90sources) sort_d.f90

f90sources = $(f90sources) ppm_util.f90

f90sources = $(f90sources) parallel_stubs.f90
f90sources = $(f90sources) omp_stubs.f90

csources = $(csources) fabio_c.c
csources = $(csources) timer_c.c
csources = $(csources) ppm_util_c.c
