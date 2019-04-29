from .Bucket import Bucket
from .Constants import my_constants
from .Amr import amr
from .Geometry import geometry
from .Algo import algo
from .Langmuirwave import langmuirwave
from .Interpolation import interpolation
from .Lasers import lasers, lasers_list
from . import Particles
from .Particles import particles, particles_list


class WarpX(Bucket):
    """
    A Python wrapper for the WarpX C++ class
    """

    def create_argv_list(self):
        argv = []
        argv += warpx.attrlist()
        argv += my_constants.attrlist()
        argv += amr.attrlist()
        argv += geometry.attrlist()
        argv += algo.attrlist()
        argv += langmuirwave.attrlist()
        argv += interpolation.attrlist()
        argv += particles.attrlist()
        argv += lasers.attrlist()

        # --- Search through species_names and add any predefined particle objects in the list.
        particles_list_names = [p.instancename for p in particles_list]
        for pstring in particles.species_names.split(' '):
            if pstring in particles_list_names:
                # --- The species is already included in particles_list
                continue
            elif hasattr(Particles, pstring):
                # --- Add the predefined species to particles_list
                particles_list.append(getattr(Particles, pstring))
                particles_list_names.append(pstring)
            else:
                raise Exception('Species %s listed in species_names not defined'%pstring)

        for particle in particles_list:
            argv += particle.attrlist()

        for laser in lasers_list:
            argv += laser.attrlist()

        return argv

    def init(self):
        from . import wx
        argv = ['warpx'] + self.create_argv_list()
        wx.initialize(argv)

    def evolve(self, nsteps=-1):
        from . import wx
        wx.evolve(nsteps)

    def finalize(self, finalize_mpi=1):
        from . import wx
        wx.finalize(finalize_mpi)

    def getProbLo(self, direction):
        from . import wx
        return wx.libwarpx.warpx_getProbLo(direction)

    def getProbHi(self, direction):
        from . import wx
        return wx.libwarpx.warpx_getProbHi(direction)

    def write_inputs(self, filename='inputs', **kw):
        argv = self.create_argv_list()
        with open(filename, 'w') as ff:

            for k, v in kw.items():
                ff.write('{0} = {1}\n'.format(k, v))

            for arg in argv:
                ff.write('{0}\n'.format(arg))

warpx = WarpX('warpx')

