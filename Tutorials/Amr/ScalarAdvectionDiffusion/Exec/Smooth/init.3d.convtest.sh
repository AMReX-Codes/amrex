#!/bin/csh

   ./main2d.gnu.DEBUG.ex smooth.inputs amr.n_cell="128 128 128" ; mv plt00000 _plt.init.128
   ./main2d.gnu.DEBUG.ex smooth.inputs amr.n_cell="256 256 256" ; mv plt00000 _plt.init.256
   ./main2d.gnu.DEBUG.ex smooth.inputs amr.n_cell="512 512 512" ; mv plt00000 _plt.init.512
   
   ${AMREX_HOME}/Tools/C_util/Convergence/RichardsonConvergenceTest2d.gnu.ex
   ${AMREX_HOME}/Tools/C_util/Convergence/RichardsonConvergenceTest2d.gnu.ex coarFile=_plt.init.128 mediFile=_plt.init.256 fineFile=_plt.init.512 coarError=_err.init.128 mediError=_err.init.256
