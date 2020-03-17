# Copyright 2017-2019 David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import numpy as np
from . import _libwarpx

class PGroup(object):
    """Implements a class that has the same API as a warp ParticleGroup instance.
    """

    def __init__(self, igroup, ispecie, level=0):
        self.igroup = igroup
        self.ispecie = ispecie
        self.level = level
        self.ns = 1 # Number of species

        self.gallot()

    def name(self):
        return 'WarpXParticleGroup'

    def gallot(self):
        self.lebcancel_pusher = 0 # turns on/off cancellation of E+VxB within V push (logical type)
        self.lebcancel = 0 # turns on/off cancellation of E+VxB before V push (logical type)

        self.sm = np.zeros(self.ns) # Species mass [kg]
        self.sq = np.zeros(self.ns) # Species charge [C]
        self.sw = np.ones(self.ns) # Species weight, (real particles per simulation particles)

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
        return _libwarpx.get_nattr()
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

    def getxp(self):
        return _libwarpx.get_particle_x(self.ispecie, self.level)[self.igroup]
    xp = property(getxp)

    def getyp(self):
        return _libwarpx.get_particle_y(self.ispecie, self.level)[self.igroup]
    yp = property(getyp)

    def getrp(self):
        return _libwarpx.get_particle_r(self.ispecie, self.level)[self.igroup]
    rp = property(getrp)

    def getzp(self):
        return _libwarpx.get_particle_z(self.ispecie, self.level)[self.igroup]
    zp = property(getzp)

    def getuxp(self):
        return _libwarpx.get_particle_ux(self.ispecie, self.level)[self.igroup]
    uxp = property(getuxp)

    def getuyp(self):
        return _libwarpx.get_particle_uy(self.ispecie, self.level)[self.igroup]
    uyp = property(getuyp)

    def getuzp(self):
        return _libwarpx.get_particle_uz(self.ispecie, self.level)[self.igroup]
    uzp = property(getuzp)

    def getw(self):
        return _libwarpx.get_particle_weight(self.ispecie, self.level)[self.igroup]

    def getpid(self, id):
        pid = _libwarpx.get_particle_arrays(self.ispecie, id, self.level)[self.igroup]
        return np.array([pid]).T

    def getgaminv(self):
        uxp = self.getuxp()
        uyp = self.getuyp()
        uzp = self.getuzp()
        return np.sqrt(1. - (uxp**2 + uyp**2 + uzp**2)/_libwarpx.clight**2)
    gaminv = property(getgaminv)

    def getex(self):
        return _libwarpx.get_particle_Ex(self.ispecie, self.level)[self.igroup]
    ex = property(getex)

    def getey(self):
        return _libwarpx.get_particle_Ey(self.ispecie, self.level)[self.igroup]
    ey = property(getey)

    def getez(self):
        return _libwarpx.get_particle_Ez(self.ispecie, self.level)[self.igroup]
    ez = property(getez)

    def getbx(self):
        return _libwarpx.get_particle_Bx(self.ispecie, self.level)[self.igroup]
    bx = property(getbx)

    def getby(self):
        return _libwarpx.get_particle_By(self.ispecie, self.level)[self.igroup]
    by = property(getby)

    def getbz(self):
        return _libwarpx.get_particle_Bz(self.ispecie, self.level)[self.igroup]
    bz = property(getbz)

    def gettheta(self):
        return _libwarpx.get_particle_theta(self.ispecie, self.level)[self.igroup]
    theta = property(gettheta)

class PGroups(object):
    def __init__(self, ispecie=0, level=0):
        self.ispecie = ispecie
        self.level = level

    def setuppgroups(self):
        xall = _libwarpx.get_particle_x(self.ispecie, self.level)
        self.ngroups = len(xall)

        self._pgroups = []
        for igroup in range(self.ngroups):
            self._pgroups.append(PGroup(igroup, self.ispecie, self.level))

    def __iter__(self):
        self.setuppgroups()
        for igroup in range(self.ngroups):
            yield self._pgroups[igroup]

    def __getitem__(self, key):
        self.setuppgroups()
        return self._pgroups[key]

    def __len__(self):
        self.setuppgroups()
        return len(self._pgroups)
