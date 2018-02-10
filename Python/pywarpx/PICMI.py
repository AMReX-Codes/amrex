"""Classes following the PICMI standard
"""
from PICMI_Base import *
import numpy as np
from pywarpx import *

codename = 'WarpX'

class Grid(PICMI_Grid):
    def init(self, **kw):

        amr.n_cell = [self.nx, self.ny, self.nz]

        # Maximum allowable size of each subdomain in the problem domain;
        #    this is used to decompose the domain for parallel calculations.
        amr.max_grid_size = kw.get('max_grid_size', 32)

        # Maximum level in hierarchy (for now must be 0, i.e., one level in total)
        amr.max_level = kw.get('max_level', 0)

        # Geometry
        geometry.coord_sys = kw.get('coord_sys', 0)  # 0: Cartesian
        geometry.is_periodic = '%d %d %d'%(self.bcxmin=='periodic', self.bcymin=='periodic', self.bczmin=='periodic')  # Is periodic?
        geometry.prob_lo = [self.xmin, self.ymin, self.zmin]  # physical domain
        geometry.prob_hi = [self.xmax, self.ymax, self.zmax]

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
        laser.wavelength = self.wavelength  # The wavelength of the laser (in meters)
        laser.e_max = self.E0  # Maximum amplitude of the laser field (in V/m)
        laser.polarization = [np.cos(self.pol_angle), np.sin(self.pol_angle), 0.]  # The main polarization vector
        laser.profile_waist = self.waist  # The waist of the laser (in meters)
        laser.profile_duration = self.duration  # The duration of the laser (in seconds)
        laser.profile_t_peak = (self.focal_position - self.z0)/clight  # The time at which the laser reaches its peak (in seconds)


class Laser_antenna(PICMI_Laser_antenna):
    def init(self, **kw):

        laser.position = [self.antenna_x0, self.antenna_y0, self.antenna_z0]  # This point is on the laser plane
        laser.direction = [self.antenna_xvec, self.antenna_yvec, self.antenna_zvec]  # The plane normal direction
        laser.profile_focal_distance = self.laser.focal_position - self.antenna_z0  # Focal distance from the antenna (in meters)


class Species(PICMI_Species):
    def init(self, **kw):

        self.species_number = particles.nspecies
        particles.nspecies = particles.nspecies + 1
        if particles.species_names is None:
            particles.species_names = self.name
        else:
            particles.species_names = particles.species_names + ' ' + self.name

        self.bucket = Bucket.Bucket(self.name, mass=self.mass, charge=self.charge, injection_style = 'python')
        Particles.particles_list.append(self.bucket)

    def add_particles(self, n=None,
                      x=None, y=None, z=None,
                      ux=None, uy=None, uz=None, w=None,
                      unique_particles=None, **kw):
        pid = np.array([w]).T
        add_particles(self.species_number, x, y, z, ux, uy, uz, pid, unique_particles)


class GaussianBeam(PICMI_GaussianBeam):
    def init(self, **kw):

        self.species.bucket.injection_style = "gaussian_beam"
        self.species.bucket.x_m = self.Xmean
        self.species.bucket.y_m = self.Ymean
        self.species.bucket.z_m = self.Zmean
        self.species.bucket.x_rms = self.Xrms
        self.species.bucket.y_rms = self.Yrms
        self.species.bucket.z_rms = self.Zrms
        self.species.bucket.npart = self.number_sim_particles
        self.species.bucket.q_tot = self.number_sim_particles*self.species.charge

        # --- These are unused but need to be set (maybe)
        self.species.bucket.profile = 'constant'
        self.species.bucket.density = 1

        # --- Momentum distribution
        if 'u_over_r' in kw:
            # --- Radial expansion
            self.species.bucket.momentum_distribution_type = "radial_expansion"
            self.species.bucket.u_over_r = kw['u_over_r']

        elif self.UXrms == 0. and self.UYrms == 0. and self.UZrms == 0.:
            # --- Constant velocity
            self.species.bucket.momentum_distribution_type = "constant"
            self.species.bucket.ux = self.UXmean
            self.species.bucket.uy = self.UYmean
            self.species.bucket.uz = self.UZmean

        else:
            # --- Gaussian velocity distribution
            self.species.bucket.momentum_distribution_type = "gaussian"
            self.species.bucket.ux_m = self.UXmean
            self.species.bucket.uy_m = self.UYmean
            self.species.bucket.uz_m = self.UZmean
            self.species.bucket.u_th = self.UZrms
            # !!! UXrms and UYrms are unused. Only an isotropic distribution is supported
            # !!! Maybe an error should be raised


class Plasma(PICMI_Plasma):
    def init(self, **kw):

        for species in self.species:
            species.bucket.injection_style = "NUniformPerCell"
            species.bucket.xmin = self.xmin
            species.bucket.xmax = self.xmax
            species.bucket.ymin = self.ymin
            species.bucket.ymax = self.ymax
            species.bucket.zmin = self.zmin
            species.bucket.zmax = self.zmax

            species.bucket.profile = 'constant'
            species.bucket.density = self.density

            if self.number_per_cell is not None:
                species.bucket.nrandompercell = self.number_per_cell
            elif self.number_per_cell_each_dim is not None:
                species.bucket.num_particles_per_cell_each_dim = self.number_per_cell_each_dim

            # --- Momentum distribution
            if 'u_over_r' in kw:
                # --- Radial expansion
                species.bucket.momentum_distribution_type = "radial_expansion"
                species.bucket.u_over_r = kw['u_over_r']

            elif self.vthx == 0. and self.vthy == 0. and self.vthz == 0.:
                # --- Constant velocity
                species.bucket.momentum_distribution_type = "constant"
                species.bucket.ux = self.vxmean
                species.bucket.uy = self.vymean
                species.bucket.uz = self.vzmean

            else:
                # --- Gaussian velocity distribution
                species.bucket.momentum_distribution_type = "gaussian"
                species.bucket.ux_m = self.vxmean
                species.bucket.uy_m = self.vymean
                species.bucket.uz_m = self.vzmean
                species.bucket.u_th = self.vthz
                # !!! vthx and vthy are unused. Only an isotropic distribution is supported
                # !!! Maybe an error should be raised


class ParticleList(PICMI_ParticleList):
    def init(self, **kw):

        assert len(self.x) == 1, "WarpX only supports initializing with a single particle"

        self.species.bucket.injection_style = "SingleParticle"
        self.species.bucket.single_particle_pos = [self.x[0], self.y[0], self.z[0]]
        self.species.bucket.single_particle_vel = [self.ux[0]/clight, self.uy[0]/clight, self.uz[0]/clight]
        self.species.bucket.single_particle_weight = self.weight


class Simulation(PICMI_Simulation):
    def set_warpx_attr(self, warpx_obj, attr, kw):
        value = kw.get(attr, None)
        if value is not None:
            setattr(warpx_obj, attr, value)
            setattr(self, attr, value)

    def init(self, **kw):

        warpx.verbose = self.verbose
        warpx.cfl = self.timestep_over_cfl
        if self.timestep == 0.:
            warpx.cfl = self.timestep_over_cfl
        else:
            warpx.const_dt = self.timestep
        amr.plot_int = self.plot_int

        self.initialized = False

    def initialize(self, inputs_name=None):
        if not self.initialized:
            self.initialized = True
            warpx.init()

    def write_inputs(self, inputs_name='inputs'):
        kw = {}
        if hasattr(self, 'max_step'):
            kw['max_step'] = self.max_step
        warpx.write_inputs(inputs_name, **kw)

    def step(self, nsteps=None):
        self.initialize()
        if nsteps is None:
            if self.max_step is not None:
                nsteps = self.max_step
            else:
                nsteps = -1
        warpx.evolve(nsteps)

    def finalize(self):
        warpx.finalize()

