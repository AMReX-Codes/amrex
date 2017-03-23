import warp
from warp import top
from warp import w3d

from ._libwarpx import libwarpx

class TimeStepper(object):

    def step(self, nsteps=1):
        for i in range(nsteps):
            self.onestep()

    def onestep(self):

        self.cur_time = libwarpx.warpx_gett_new(0)
        self.istep = libwarpx.warpx_getistep(0)

        #if mpi.rank == 0:
        print "\nSTEP %d starts ..."%(self.istep + 1)

        #if (ParallelDescriptor::NProcs() > 1)
        #   if (okToRegrid(step)) RegridBaseLevel();

        libwarpx.warpx_ComputeDt()
        dt = libwarpx.warpx_getdt(0)

        # --- Advance level 0 by dt
        lev = 0

        # --- At the beginning, we have B^{n-1/2} and E^{n}.
        # --- Particles have p^{n-1/2} and x^{n}.
        libwarpx.warpx_EvolveB(lev, 0.5*dt) # We now B^{n}

        libwarpx.warpx_FillBoundaryE(lev, False)
        libwarpx.warpx_FillBoundaryB(lev, False)

        # --- Evolve particles to p^{n+1/2} and x^{n+1}
        # --- Depose current, j^{n+1/2}
        libwarpx.warpx_PushParticlesandDepose(lev, self.cur_time)

        libwarpx.mypc_Redistribute() # Redistribute particles

        libwarpx.warpx_EvolveB(lev, 0.5*dt) # We now B^{n+1/2}

        libwarpx.warpx_FillBoundaryB(lev, True)

        libwarpx.warpx_EvolveE(lev, dt) # We now have E^{n+1}

        self.istep += 1

        self.cur_time += dt

        libwarpx.warpx_MoveWindow();

        #if mpi.rank == 0:
        print "STEP %d ends. TIME = %e DT = %e"%(self.istep, self.cur_time, dt)

        # --- Sync up time
        for i in range(libwarpx.warpx_finestLevel()+1):
            libwarpx.warpx_sett_new(i, self.cur_time)
            libwarpx.warpx_setistep(i, self.istep)

        max_time_reached = ((self.cur_time >= libwarpx.warpx_stopTime() - 1.e-6*dt) or (self.istep >= libwarpx.warpx_maxStep()))

        if libwarpx.warpx_plotInt() > 0 and (self.istep+1)%libwarpx.warpx_plotInt() == 0 or max_time_reached:
            libwarpx.warpx_WritePlotFile()

        if libwarpx.warpx_checkInt() > 0 and (self.istep+1)%libwarpx.warpx_plotInt() == 0 or max_time_reached:
            libwarpx.warpx_WriteCheckPointFile()












# --- This is not used
class TimeStepperFromPICSAR(object):

    def __init__(self, package=None, solver=None, l_debug=False):
        self.package = package
        self.solver = solver
        self.l_debug = l_debug

    def setpackagestepnumber(self, it):
        if self.package is not None:
            self.package.setstepnumber(it)

    def step(self, n=1, freq_print=10, lallspecl=0):
      """
      This function performs a range of Particle-In-Cell iterations

      Inputs:
      - n: number of iterations
      - freq_print: print frequency
      """

      if (self.l_debug): print("Call step")

      for i in range(n):
          if(me == 0):
              if top.it%freq_print == 0:
                  print 'it = %g time = %g'%(top.it, top.time)

          l_first = (lallspecl or (i == 0))
          l_last = (lallspecl or (i == n-1))

          self.onestep(l_first, l_last)

      if (self.l_debug): print("End step")

    def onestep(self, l_first, l_last):
        """
        Perform a single particle-in-cell step
        """

        if (self.l_debug): print("Call onestep")

        # --- Iteration number
        self.setpackagestepnumber(top.it)

        # --- call beforestep functions
        if (self.l_debug): print("Call beforestep functions")
        warp.callbeforestepfuncs.callfuncsinlist()

        # --- gather fields from grid to particles
        if (self.l_debug): print("Call Field gathering and particle push")

        # --- push
        if l_first:
            if self.package is None:
                # --- Standard Warp advance
                for specie in warp.listofallspecies:
                    for pg in specie.flatten(specie.pgroups):
                        for js in range(pg.ns):
                            self.push_velocity_second_part(js, pg)
                            self.push_positions(js, pg)
                        warp.particleboundaries3d(pg, -1, False)

            else:
                # --- Particle pusher
                if (self.l_debug): print("Call package particle push")
                pxr.pxrpush_particles_part2()

                # --- Particle boundary consitions
                if (self.l_debug): print("Call package particle boundary conditions")
                pxr.particle_bcs()

                if (self.l_debug): print("Call aliasparticlearrays()")
                self.aliasparticlearrays()

        else:
            if self.package is None:
                # --- Standard Warp advance

                for specie in warp.listofallspecies:
                    for pg in specie.flatten(specie.pgroups):
                        for js in range(pg.ns):
                            self.push_velocity_full(js, pg)
                            self.push_positions(js, pg)

                        warp.particleboundaries3d(pg, -1, False)

            else:
                # --- Particle pusher
                if (self.l_debug): print("Call package particle pusher")
                pxr.field_gathering_plus_particle_pusher()

                # --- Particle boundary conditions
                if (self.l_debug): print("Call package particle boundary conditions")
                pxr.particle_bcs()

                self.aliasparticlearrays()

        # --- Particle sorting
        if (self.l_debug): print("Call Particle Sorting")
        if self.package is not None:
            # --- This should be a function installed before load rho
            pxr.particle_sorting_sub()

        # --- call beforeloadrho functions
        if (self.l_debug): print("Call beforeloadrho functions")
        warp.beforeloadrho.callfuncsinlist()

        pgroups = []
        for specie in warp.listofallspecies:
            pgroups += specie.flatten(specie.pgroups)

        self.pgroups = pgroups

        # --- Call user-defined injection routines
        if (self.l_debug): print("Call user-defined injection routines")
        warp.userinjection.callfuncsinlist()

        xgriddiff = w3d.xmmin - pxr.xmin
        ygriddiff = w3d.ymmin - pxr.ymin
        zgriddiff = w3d.zmmin - pxr.zmin

        if (xgriddiff != 0 or ygriddiff != 0 or zgriddiff != 0):
            pxr.pxr_move_sim_boundaries(xgriddiff, ygriddiff, zgriddiff)


        if (self.l_debug): print("Call loadrho")
        self.solver.loadrho(pgroups = pgroups)
        if (self.l_debug): print("Call loadj")
        self.solver.loadj(pgroups = pgroups)

        if (self.l_debug): print("Call dosolve")
        self.solver.dosolve()

        if self.package is None:
            for specie in warp.listofallspecies:
                for pg in specie.flatten(specie.pgroups):
                    for js in range(pg.ns):
                        self.fetcheb(js, pg)
                        if l_last:
                            self.push_velocity_first_part(js, pg)

        else:
            if l_last:
                if (self.l_debug): print("Call package push particles 1")
                pxr.pxrpush_particles_part1()

        # --- update time, time counter
        top.time += top.dt
        if top.it%top.nhist == 0:
#           zmmnt()
           minidiag(top.it, top.time, top.lspecial)
        top.it += 1

        # --- Load balance function should be installed after step

        # --- call afterstep functions
        if (self.l_debug): print("Call callafterstepfuncs.callfuncsinlist()")
        warp.callafterstepfuncs.callfuncsinlist()


    # --- The following methods use the standard Warp routines
    def fetcheb(self, js, pg=None):
        if self.l_verbose:print me, 'enter fetcheb'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np == 0: return
        il = pg.ins[js] - 1
        iu = il + pg.nps[js]
        w3d.pgroupfsapi = pg
        w3d.ipminfsapi = pg.ins[js]
        w3d.npfsapi = pg.nps[js]
        pg.ex[il:iu] = 0.
        pg.ey[il:iu] = 0.
        pg.ez[il:iu] = 0.
        pg.bx[il:iu] = 0.
        pg.by[il:iu] = 0.
        pg.bz[il:iu] = 0.
        self.fetche()
        self.fetchb()
        w3d.pgroupfsapi = None

    def push_velocity_full(self, js, pg=None):
        if self.l_verbose:print me, 'enter push_velocity_full'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np == 0: return
        il = pg.ins[js] - 1
        iu = il + pg.nps[js]
        if pg.lebcancel_pusher:
          warp.ebcancelpush3d(np, pg.uxp[il:iu], pg.uyp[il:iu], pg.uzp[il:iu], pg.gaminv[il:iu],
                              pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                              pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                              pg.sq[js], pg.sm[js], top.dt, 0)
        else:
          # --- push velocity from electric field (half step)
          warp.epush3d(np, pg.uxp[il:iu], pg.uyp[il:iu], pg.uzp[il:iu],
                       pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                       pg.sq[js], pg.sm[js], 0.5*top.dt)
          # --- update gamma
          self.set_gamma(js, pg)
          # --- push velocity from magnetic field
          warp.bpush3d(np, pg.uxp[il:iu], pg.uyp[il:iu], pg.uzp[il:iu], pg.gaminv[il:iu],
                       pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                       pg.sq[js], pg.sm[js], top.dt, top.ibpush)
          # --- push velocity from electric field (half step)
          warp.epush3d(np, pg.uxp[il:iu], pg.uyp[il:iu], pg.uzp[il:iu],
                       pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                       pg.sq[js], pg.sm[js], 0.5*top.dt)
          # --- update gamma
          self.set_gamma(js, pg)

        if self.l_verbose:print me, 'exit push_velocity_first_part'

    def push_velocity_first_part(self, js, pg=None):
        if self.l_verbose:print me, 'enter push_velocity_first_part'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np == 0: return
        il = pg.ins[js] - 1
        iu = il + pg.nps[js]
        if pg.lebcancel_pusher:
          warp.ebcancelpush3d(np, pg.uxp[il:iu], pg.uyp[il:iu], pg.uzp[il:iu], pg.gaminv[il:iu],
                              pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                              pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                              pg.sq[js], pg.sm[js], top.dt, 1)
        else:
          # --- push velocity from electric field (half step)
          warp.epush3d(np, pg.uxp[il:iu], pg.uyp[il:iu], pg.uzp[il:iu],
                      pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                      pg.sq[js], pg.sm[js], 0.5*top.dt)
          # --- update gamma
          self.set_gamma(js,pg)
          # --- push velocity from magnetic field
          warp.bpush3d(np, pg.uxp[il:iu], pg.uyp[il:iu], pg.uzp[il:iu], pg.gaminv[il:iu],
                       pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                       pg.sq[js], pg.sm[js], 0.5*top.dt, top.ibpush)

        if self.l_verbose:print me, 'exit push_velocity_first_part'

    def push_velocity_second_part(self, js, pg=None):
        if self.l_verbose:print me, 'enter push_velocity_second_part'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np == 0: return
        il = pg.ins[js] - 1
        iu = il + pg.nps[js]
        if pg.lebcancel_pusher:
          warp.ebcancelpush3d(np, pg.uxp[il:iu], pg.uyp[il:iu], pg.uzp[il:iu], pg.gaminv[il:iu],
                              pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                              pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                              pg.sq[js], pg.sm[js], top.dt, 2)
        else:
          # --- push velocity from magnetic field
          warp.bpush3d(np, pg.uxp[il:iu], pg.uyp[il:iu], pg.uzp[il:iu], pg.gaminv[il:iu],
                       pg.bx[il:iu], pg.by[il:iu], pg.bz[il:iu],
                       pg.sq[js], pg.sm[js], 0.5*top.dt, top.ibpush)
          # --- push velocity from electric field (half step)
          warp.epush3d(np, pg.uxp[il:iu], pg.uyp[il:iu], pg.uzp[il:iu],
                       pg.ex[il:iu], pg.ey[il:iu], pg.ez[il:iu],
                       pg.sq[js], pg.sm[js], 0.5*top.dt)
        # --- update gamma
        self.set_gamma(js, pg)

        if self.l_verbose:print me, 'exit push_velocity_second_part'

    def set_gamma(self, js, pg=None):
        if self.l_verbose:print me, 'enter set_gamma'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np == 0: return
        il = pg.ins[js] - 1
        iu = il + pg.nps[js]
        # --- update gamma
        warp.gammaadv(np, pg.gaminv[il:iu], pg.uxp[il:iu], pg.uyp[il:iu], pg.uzp[il:iu], top.gamadv, top.lrelativ)

        if self.l_verbose:print me, 'exit set_gamma'

    def push_positions(self, js, pg=None):
        if self.l_verbose:print me, 'enter push_positions'
        if pg is None:
            pg = top.pgroup
        np = pg.nps[js]
        if np == 0: return
        il = pg.ins[js] - 1
        iu = il + pg.nps[js]
        warp.xpush3d(np, pg.xp[il:iu], pg.yp[il:iu], pg.zp[il:iu],
                     pg.uxp[il:iu], pg.uyp[il:iu], pg.uzp[il:iu],
                     pg.gaminv[il:iu], top.dt)

        if self.l_verbose:print me, 'exit push_positions'








