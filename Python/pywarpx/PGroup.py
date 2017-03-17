import numpy as np
from . import WarpX
from . import warpxC

class PGroup(object):
    """Implements a class that has the same API as a warp ParticleGroup instance.
    """

    def __init__(self):
        self.ns = 1 # Number of species
        self.npmax = 0 # Size of data arrays
        self.npid = 0 # number of columns for pid.

        self.gallot()

    def name(self):
        return 'WarpXParticleGroup'

    def gallot(self):
        self.lebcancel_pusher = 0 # turns on/off cancellation of E+VxB within V push (logical type)
        self.lebcancel = 0 # turns on/off cancellation of E+VxB before V push (logical type)

        self.sm = np.zeros(self.ns) # Species mass [kg]
        self.sq = np.zeros(self.ns) # Species charge [C]
        self.sw = np.zeros(self.ns) # Species weight, (real particles per simulation particles)

        self.ins = np.ones(self.ns, dtype=int) # Index of first particle in species
        self.nps = np.zeros(self.ns, dtype=int) # Number of particles in species
        self.ipmax = np.zeros(self.ns+1, dtype=int) # Max extent within the arrays of each species

        self.sid = np.arange(self.ns, dtype=int) # Global species index for each species
        self.ndts = np.ones(self.ns, dtype=int) # Stride for time step advance for each species
        self.ldts = np.ones(self.ns, dtype=int) # (logical type)
        self.lvdts = np.ones(self.ns, dtype=int) # (logical type)
        self.iselfb = np.zeros(self.ns, dtype=int) # Group number for particles that are affected by
                                                   # their own magnetic field, using the 1/gamma**2
                                                   # approximation. The correction is not applied to
                                                   # group number -1.
        self.fselfb = np.zeros(self.ns) # The scaling factor, vz.
        self.l_maps = np.zeros(self.ns)
        self.dtscale = np.ones(self.ns) # Scale factor applied to time step size for each
                                        # species. Only makes sense in steaday and and
                                        # transverse slice modes.

        self.limplicit = np.zeros(self.ns, dtype=int) # Flags implicit particle species (logical type)
        self.iimplicit = np.full(self.ns, -1, dtype=int) # Group number for implicit particles
        self.ldoadvance = np.ones(self.ns, dtype=int) # Flags whether particles are time advanced (logical type)
        self.lboundaries = np.ones(self.ns, dtype=int) # Flags whether boundary conditions need to be applied (logical type)
        self.lparaxial = np.zeros(self.ns, dtype=int) # Flags to turn on/off paraxial approximation (logical type)

        self.zshift = np.zeros(self.ns)
        self.gamma_ebcancel_max = np.ones(self.ns) # maximum value allowed for ExB cancellation

        self.gaminv = np.ones(self.npmax) # inverse relativistic gamma factor
        self._xp = np.zeros(self.npmax) # X-positions of particles [m]
        self._yp = np.zeros(self.npmax) # Y-positions of particles [m]
        self._zp = np.zeros(self.npmax) # Z-positions of particles [m]
        self._uxp = np.zeros(self.npmax) # gamma * X-velocities of particles [m/s]
        self._uyp = np.zeros(self.npmax) # gamma * Y-velocities of particles [m/s]
        self._uzp = np.zeros(self.npmax) # gamma * Z-velocities of particles [m/s]
        self._ex = np.zeros(self.npmax) # Ex of particles [V/m]
        self._ey = np.zeros(self.npmax) # Ey of particles [V/m]
        self._ez = np.zeros(self.npmax) # Ez of particles [V/m]
        self._bx = np.zeros(self.npmax) # Bx of particles [T]
        self._by = np.zeros(self.npmax) # By of particles [T]
        self._bz = np.zeros(self.npmax) # Bz of particles [T]
        self._pid = np.zeros((self.npmax, self.npid)) # Particle ID - used for various purposes

    # --- Temporary fix
    gchange = gallot

    def allocated(self, name):
        return (getattr(self, name, None) is not None)

    def addspecies(self):
        pass

    def _updatelocations(self):
        warpx = WarpX.warpx.warpx
        mypc = warpx.GetPartContainer()

        xplist = []
        yplist = []
        zplist = []
        for ispecie in range(mypc.nSpecies()):
            pc = mypc.GetParticleContainer(ispecie)
            xx = pc.getLocations()
            xplist.append(xx[0,:])
            yplist.append(xx[1,:])
            zplist.append(xx[2,:])
            self.nps[ispecie] = len(xplist[-1])
            if ispecie > 0:
                self.ins[ispecie] = self.ins[ispecie-1] + self.nps[ispecie-1]
                self.ipmax[ispecie+1] = self.ins[ispecie] + self.nps[ispecie] - 1

        self._xp = np.concatenate(xplist)
        self._yp = np.concatenate(yplist)
        self._zp = np.concatenate(zplist)
        self.npmax = len(self._xp)

    def _updatevelocities(self):
        warpx = WarpX.warpx.warpx
        mypc = warpx.GetPartContainer()

        uxplist = []
        uyplist = []
        uzplist = []
        for ispecie in range(mypc.nSpecies()):
            pc = mypc.GetParticleContainer(ispecie)
            vv = pc.getData(0, 3)
            uxplist.append(vv[0,:])
            uyplist.append(vv[1,:])
            uzplist.append(vv[2,:])
            self.nps[ispecie] = len(uxplist[-1])
            if ispecie > 0:
                self.ins[ispecie] = self.ins[ispecie-1] + self.nps[ispecie-1]
                self.ipmax[ispecie+1] = self.ins[ispecie] + self.nps[ispecie] - 1

        self._uxp = np.concatenate(uxplist)
        self._uyp = np.concatenate(uyplist)
        self._uzp = np.concatenate(uzplist)
        self.npmax = len(self._xp)

    def _updatepids(self):
        warpx = WarpX.warpx.warpx
        mypc = warpx.GetPartContainer()

        pidlist = []
        for ispecie in range(mypc.nSpecies()):
            pc = mypc.GetParticleContainer(ispecie)
            self.npid = pc.nAttribs - 3
            vv = pc.getData(3, self.npid)
            pidlist.append(vv)
            self.nps[ispecie] = len(uxplist[-1])
            if ispecie > 0:
                self.ins[ispecie] = self.ins[ispecie-1] + self.nps[ispecie-1]
                self.ipmax[ispecie+1] = self.ins[ispecie] + self.nps[ispecie] - 1

        self._pid = np.concatenate(pidlist.T, axis=0)
        self.npmax = self._pid.shape[0]

    def getxp(self):
        self._updatelocations()
        return self._xp
    xp = property(getxp)

    def getyp(self):
        self._updatelocations()
        return self._yp
    yp = property(getyp)

    def getzp(self):
        self._updatelocations()
        return self._zp
    zp = property(getzp)

    def getuxp(self):
        self._updatevelocities()
        return self._uxp
    uxp = property(getuxp)

    def getuyp(self):
        self._updatevelocities()
        return self._uyp
    uyp = property(getuyp)

    def getuzp(self):
        self._updatevelocities()
        return self._uzp
    uzp = property(getuzp)

    def getpid(self):
        self._updatepids()
        return self._pid
    pid = property(getpid)

    def getgaminv(self):
        uxp = self.uxp
        uyp = self.uyp
        uzp = self.uzp
        return sqrt(1. - (uxp**2 + uyp**2 + uzp**2)/warpxC.c**2)
    gaminv = property(getgaminv)

    def getex(self):
        return np.zeros(self.npmax)
    ex = property(getex)
    def getey(self):
        return np.zeros(self.npmax)
    ey = property(getey)
    def getez(self):
        return np.zeros(self.npmax)
    ez = property(getez)
    def getbx(self):
        return np.zeros(self.npmax)
    bx = property(getbx)
    def getby(self):
        return np.zeros(self.npmax)
    by = property(getby)
    def getbz(self):
        return np.zeros(self.npmax)
    bz = property(getbz)

