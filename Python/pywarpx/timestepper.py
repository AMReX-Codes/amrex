# Copyright 2017-2018 Andrew Myers, David Grote, Weiqun Zhang
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from ._libwarpx import libwarpx
from . import callbacks

class TimeStepper(object):

    def step(self, nsteps=1):
        for i in range(nsteps):
            self.onestep()

    def onestep(self):

        callbacks._beforestep()

        self.cur_time = libwarpx.warpx_gett_new(0)
        self.istep = libwarpx.warpx_getistep(0)

        #if mpi.rank == 0:
        print("\nSTEP %d starts ..."%(self.istep + 1))

        #if (ParallelDescriptor::NProcs() > 1)
        #   if (okToRegrid(step)) RegridBaseLevel();

        dt = libwarpx.warpx_getdt(0)

        # --- At the beginning, we have B^{n-1/2} and E^{n}.
        # --- Particles have p^{n-1/2} and x^{n}.
        libwarpx.warpx_FillBoundaryE()
        libwarpx.warpx_EvolveB(0.5*dt,1) # We now B^{n}

        libwarpx.warpx_FillBoundaryB()
        libwarpx.warpx_UpdateAuxilaryData()

        # --- Evolve particles to p^{n+1/2} and x^{n+1}
        # --- Depose current, j^{n+1/2}
        callbacks._particleinjection()
        callbacks._particlescraper()
        callbacks._beforedeposition()
        libwarpx.warpx_PushParticlesandDepose(self.cur_time)
        callbacks._afterdeposition()

        libwarpx.mypc_Redistribute() # Redistribute particles

        libwarpx.warpx_FillBoundaryE()
        libwarpx.warpx_EvolveB(0.5*dt,2) # We now B^{n+1/2}

        libwarpx.warpx_SyncCurrent()

        libwarpx.warpx_FillBoundaryB()
        callbacks._beforeEsolve()
        libwarpx.warpx_EvolveE(dt,0) # We now have E^{n+1}
        callbacks._afterEsolve()

        self.istep += 1

        self.cur_time += dt

        libwarpx.warpx_MoveWindow();

        #if mpi.rank == 0:
        print("STEP %d ends. TIME = %e DT = %e"%(self.istep, self.cur_time, dt))

        # --- Sync up time
        for i in range(libwarpx.warpx_finestLevel()+1):
            libwarpx.warpx_sett_new(i, self.cur_time)
            libwarpx.warpx_setistep(i, self.istep)

        max_time_reached = ((self.cur_time >= libwarpx.warpx_stopTime() - 1.e-6*dt) or (self.istep >= libwarpx.warpx_maxStep()))

        if (libwarpx.warpx_plotInt() > 0 and (self.istep+1)%libwarpx.warpx_plotInt() == 0) or max_time_reached:
            libwarpx.warpx_WritePlotFile()
        if (libwarpx.warpx_openpmdInt() > 0 and (self.istep+1)%libwarpx.warpx_openpmdInt() == 0) or max_time_reached:
            libwarpx.warpx_WriteOpenPMDFile()

        if libwarpx.warpx_checkInt() > 0 and (self.istep+1)%libwarpx.warpx_plotInt() == 0 or max_time_reached:
            libwarpx.warpx_WriteCheckPointFile()

        callbacks._afterstep()
