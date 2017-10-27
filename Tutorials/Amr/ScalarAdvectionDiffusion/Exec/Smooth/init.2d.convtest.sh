#!/bin/csh

   ./main2d.gnu.DEBUG.ex smooth.inputs amr.n_cell="128 128 128" ; mv plt00000 _plt.2d.init.128
   ./main2d.gnu.DEBUG.ex smooth.inputs amr.n_cell="256 256 256" ; mv plt00000 _plt.2d.init.256
   ./main2d.gnu.DEBUG.ex smooth.inputs amr.n_cell="512 512 512" ; mv plt00000 _plt.2d.init.512
   
   ${AMREX_HOME}/Tools/C_util/Convergence/RichardsonConvergenceTest2d.gnu.ex
   ${AMREX_HOME}/Tools/C_util/Convergence/RichardsonConvergenceTest2d.gnu.ex coarFile=_plt.2d.init.128 mediFile=_plt.2d.init.256 fineFile=_plt.2d.init.512 coarError=_err2d.init.128 mediError=_err2d.init.256 |& tee init_convtest.2d.out

