"""Classes following the PICMI standard
"""
from PICMI_Base import *
import numpy as np
from pywarpx import *

codename = 'WarpX'

def _args_to_string(*args):
    # --- Converts of sequence of number to a string that is appropriate for input.
    return ' '.join(map(repr, args))

class Grid(PICMI_Grid):
    def init(self, **kw):

        amr.n_cell = _args_to_string(self.nx, self.ny, self.nz)

        # Maximum allowable size of each subdomain in the problem domain;
        #    this is used to decompose the domain for parallel calculations.
        amr.max_grid_size = kw.get('max_grid_size', 32)

        # Maximum level in hierarchy (for now must be 0, i.e., one level in total)
        amr.max_level = kw.get('max_level', 0)

        # Geometry
        geometry.coord_sys = kw.get('coord_sys', 0)  # 0: Cartesian
        geometry.is_periodic = '%d %d %d'%(self.bcxmin=='periodic', self.bcymin=='periodic', self.bczmin=='periodic')  # Is periodic?
        geometry.prob_lo = _args_to_string(self.xmin, self.ymin, self.zmin)  # physical domain
        geometry.prob_hi = _args_to_string(self.xmax, self.ymax, self.zmax)

        if self.moving_window_velocity is not None and np.any(self.moving_window_velocity != 0):
            warpx.do_moving_window = 1
            if self.moving_window_velocity[0] != 0.:
                warpx.moving_window_dir = 'x'
                warpx.moving_window_v = self.moving_window_velocity[0]/clight  # in units of the speed of light
            if self.moving_window_velocity[1] != 0.:
                warpx.moving_window_dir = 'y'
                warpx.moving_window_v = self.moving_window_velocity[1]/clight  # in units of the speed of light
            if self.moving_window_velocity[2] != 0.:
                warpx.moving_window_dir = 'z'
                warpx.moving_window_v = self.moving_window_velocity[2]/clight  # in units of the speed of light

    def getmins(self, **kw):
        return np.array([warpx.getProbLo(0), warpx.getProbLo(1), warpx.getProbLo(2)])

    def getmaxs(self, **kw):
        return np.array([warpx.getProbHi(0), warpx.getProbHi(1), warpx.getProbHi(2)])

    def getxmin(self):
        return warpx.getProbLo(0)
    def getxmax(self):
        return warpx.getProbHi(0)
    def getymin(self):
        return warpx.getProbLo(1)
    def getymax(self):
        return warpx.getProbHi(1)
    def getzmin(self):
        return warpx.getProbLo(2)
    def getzmax(self):
        return warpx.getProbHi(2)


class EM_solver(PICMI_EM_solver):
    def init(self, **kw):

        if self.current_deposition_algo is not None:
            algo.current_deposition = self.current_deposition_algo
        if self.charge_deposition_algo is not None:
            algo.charge_deposition = self.charge_deposition_algo
        if self.field_gathering_algo is not None:
            algo.field_gathering = self.field_gathering_algo
        if self.particle_pusher_algo is not None:
            algo.particle_pusher = self.particle_pusher_algo


class Gaussian_laser(PICMI_Gaussian_laser):
    def init(self, **kw):

        warpx.use_laser = 1
        laser.profile = "Gaussian"
        laser.position = _args_to_string(self.antenna_x0, self.antenna_y0, self.antenna_z0)  # This point is on the laser plane
        laser.direction = _args_to_string(self.antenna_xvec, self.antenna_yvec, self.antenna_zvec)  # The plane normal direction
        laser.polarization = _args_to_string(np.cos(self.pol_angle), np.sin(self.pol_angle), 0.)  # The main polarization vector
        laser.e_max = self.E0  # Maximum amplitude of the laser field (in V/m)
        laser.profile_waist = self.waist  # The waist of the laser (in meters)
        laser.profile_duration = self.duration  # The duration of the laser (in seconds)
        laser.profile_t_peak = self.t_peak  # The time at which the laser reaches its peak (in seconds)
        laser.profile_focal_distance = self.focal_position - self.antenna_z0  # Focal distance from the antenna (in meters)
        laser.wavelength = self.wavelength  # The wavelength of the laser (in meters)


class Species(PICMI_Species):
    def init(self, **kw):

        self.species_number = particles.nspecies
        particles.nspecies = particles.nspecies + 1
        particles.species_names = particles.species_names + ' ' + self.name

        self.bucket = Bucket.Bucket(self.name, mass=self.mass, charge=self.charge, injection_style = 'python')
        Particles.particles_list.append(self.bucket)

    def add_particles(self, n=None,
                      x=None, y=None, z=None,
                      ux=None, uy=None, uz=None, w=None,
                      unique_particles=None, **kw):
        pid = np.array([w]).T
        add_particles(self.species_number, x, y, z, ux, uy, uz, pid, unique_particles)


class Simulation(PICMI_Simulation):
    def set_warpx_attr(self, warpx_obj, attr, kw):
        value = kw.get(attr, None)
        if value is not None:
            setattr(warpx_obj, attr, value)
            setattr(self, attr, value)

    def init(self, **kw):

        warpx.verbose = self.verbose
        warpx.cfl = self.cfl
        amr.plot_int = self.plot_int

        self.amrex = AMReX()
        self.amrex.init()
        warpx.init()

    def step(self, nsteps=-1):
        warpx.evolve(nsteps)

    def finalize(self):
        warpx.finalize()
        self.amrex.finalize()

