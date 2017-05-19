import numpy as np
from . import WarpX
from . import _libwarpx

class PGroup(object):
    """Implements a class that has the same API as a warp ParticleGroup instance.
    """

    def __init__(self, igroup):
        self.igroup = igroup
        self.ns = 1 # Number of species

        self.gallot()

    def name(self):
        return 'WarpXParticleGroup'

    def gallot(self):
        self.lebcancel_pusher = 0 # turns on/off cancellation of E+VxB within V push (logical type)
        self.lebcancel = 0 # turns on/off cancellation of E+VxB before V push (logical type)

        self.sm = np.zeros(self.ns) # Species mass [kg]
        self.sq = np.zeros(self.ns) # Species charge [C]
        self.sw = np.zeros(self.ns) # Species weight, (real particles per simulation particles)

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

    # --- Temporary fix
    gchange = gallot

    def allocated(self, name):
        return True

    def addspecies(self):
        pass

    def getnpid(self):
        return _labwarpx.get_nattr()
    npid = property(getnpid)

    def getnps(self):
        return np.array([len(self.xp)], dtype='l')
    nps = property(getnps)

    def getins(self):
        return np.ones(self.ns, dtype='l')
    ins = property(getins)

    def getipmax(self):
        return np.array([0, len(self.xp)], dtype='l')
    ipmax = property(getipmax)

    def getnpmax(self):
        return self.nps.sum()
    npmax = property(getnpmax)

    def getxp(self, js=0):
        return _libwarpx.get_particle_x(js)[self.igroup]
    xp = property(getxp)

    def getyp(self, js=0):
        return _libwarpx.get_particle_y(js)[self.igroup]
    yp = property(getyp)

    def getzp(self, js=0):
        return _libwarpx.get_particle_z(js)[self.igroup]
    zp = property(getzp)

    def getuxp(self, js=0):
        return _libwarpx.get_particle_ux(js)[self.igroup]
    uxp = property(getuxp)

    def getuyp(self, js=0):
        return _libwarpx.get_particle_uy(js)[self.igroup]
    uyp = property(getuyp)

    def getuzp(self, js=0):
        return _libwarpx.get_particle_uz(js)[self.igroup]
    uzp = property(getuzp)

    def getpid(self, js=0):
        return _libwarpx.get_particle_id(js)[self.igroup]
    pid = property(getpid)

    def getgaminv(self, js=0):
        uxp = self.getuxp(js)
        uyp = self.getuyp(js)
        uzp = self.getuzp(js)
        return np.sqrt(1. - (uxp**2 + uyp**2 + uzp**2)/_libwarpx.clight**2)
    gaminv = property(getgaminv)

    def getex(self, js=0):
        return _libwarpx.get_particle_Ex(js)[self.igroup]
    ex = property(getex)

    def getey(self, js=0):
        return _libwarpx.get_particle_Ey(js)[self.igroup]
    ey = property(getey)

    def getez(self, js=0):
        return _libwarpx.get_particle_Ez(js)[self.igroup]
    ez = property(getez)

    def getbx(self, js=0):
        return _libwarpx.get_particle_Bx(js)[self.igroup]
    bx = property(getbx)

    def getby(self, js=0):
        return _libwarpx.get_particle_By(js)[self.igroup]
    by = property(getby)

    def getbz(self, js=0):
        return _libwarpx.get_particle_Bz(js)[self.igroup]
    bz = property(getbz)

class PGroups(object):
    def __init__(self, ispecie=0):
        self.ispecie = ispecie
        xall = _libwarpx.get_particle_x(ispecie)
        self.ngroups = len(xall)

        self._pgroups = []
        for igroup in range(self.ngroups):
            self._pgroups.append(PGroup(igroup))
        
    def __iter__(self):
        for igroup in range(self.ngroups):
            yield self._pgroups[igroup]

    def __getitem__(self, key):
        return self._pgroups[key]

