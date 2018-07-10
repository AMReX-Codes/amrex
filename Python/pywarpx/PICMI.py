"""Classes following the PICMI standard
"""
import PICMI_Base
import numpy as np
import pywarpx

codename = 'WarpX'

# --- Values from WarpXConst.H
c = 299792458.
ep0 = 8.854187817e-12
mu0 = 1.2566370614359173e-06
q_e = 1.602176462e-19
m_e = 9.10938291e-31
m_p = 1.6726231e-27


class Species(PICMI_Base.PICMI_Species):
    def init(self, **kw):

        if self.type == 'electron':
            if self.charge is None: self.charge = '-q_e'
            if self.mass is None: self.mass = 'm_e'
        elif self.type == 'positron':
            if self.charge is None: self.charge = 'q_e'
            if self.mass is None: self.mass = 'm_e'
        elif self.type == 'proton':
            if self.charge is None: self.charge = 'q_e'
            if self.mass is None: self.mass = 'm_p'
        elif self.type == 'anti-proton':
            if self.charge is None: self.charge = '-q_e'
            if self.mass is None: self.mass = 'm_p'

        self.species_number = pywarpx.particles.nspecies
        pywarpx.particles.nspecies += 1

        if self.name is None:
            self.name = 'species{}'.format(self.species_number)

        if pywarpx.particles.species_names is None:
            pywarpx.particles.species_names = self.name
        else:
            pywarpx.particles.species_names += ' ' + self.name

        self.bucket = pywarpx.Bucket.Bucket(self.name, mass=self.mass, charge=self.charge, injection_style = 'python')
        pywarpx.Particles.particles_list.append(self.bucket)


class GaussianBeam(PICMI_Base.PICMI_GaussianBeam):
    def init(self, **kw):

        self.species.bucket.injection_style = "gaussian_beam"
        self.species.bucket.x_m = self.Xmean
        self.species.bucket.y_m = self.Ymean
        self.species.bucket.z_m = self.Zmean
        self.species.bucket.x_rms = self.Xrms
        self.species.bucket.y_rms = self.Yrms
        self.species.bucket.z_rms = self.Zrms
        self.species.bucket.npart = self.number_sim_particles

        # --- Calculate the total charge. Note that charge might be a string instead of a number.
        charge = self.species.bucket.charge
        if charge == 'q_e' or charge == '+q_e':
            charge = q_e
        elif charge == '-q_e':
            charge = -q_e
        self.species.bucket.q_tot = self.number_real_particles*charge

        # --- These need to be defined even though they are not used
        self.species.bucket.profile = "constant"
        self.species.bucket.density = 1

        # --- The PICMI standard doesn't yet have a way of specifying these values.
        # --- They should default to the size of the domain. They are not typically
        # --- necessary though since any particles outside the domain are rejected.
        #self.species.bucket.xmin
        #self.species.bucket.xmax
        #self.species.bucket.ymin
        #self.species.bucket.ymax
        #self.species.bucket.zmin
        #self.species.bucket.zmax

        if self.UXdiv != 0. and self.UYdiv != 0. and self.UZdiv != 0.:
            self.species.bucket.momentum_distribution_type = "radial_expansion"
            self.species.bucket.u_over_r = self.UXdiv
            #self.species.bucket.u_over_y = self.UYdiv
            #self.species.bucket.u_over_z = self.UZdiv
        elif self.UXrms != 0. or self.UYrms != 0. or self.UZrms != 0.:
            self.species.bucket.momentum_distribution_type = "gaussian"
            self.species.bucket.ux_m = self.UXmean
            self.species.bucket.uy_m = self.UYmean
            self.species.bucket.uz_m = self.UZmean
            self.species.bucket.ux_th = self.UXrms
            self.species.bucket.uy_th = self.UYrms
            self.species.bucket.uz_th = self.UZrms
        else:
            self.species.bucket.momentum_distribution_type = "constant"
            self.species.bucket.ux = self.UXmean
            self.species.bucket.uy = self.UYmean
            self.species.bucket.uz = self.UZmean


class Plasma(PICMI_Base.PICMI_Plasma):
    def init(self, **kw):

        for species in self.species:
            if self.number_per_cell_each_dim is not None:
                species.bucket.injection_style = "nuniformpercell"
                species.bucket.num_particles_per_cell_each_dim = self.number_per_cell_each_dim
            elif self.number_per_cell is not None:
                species.bucket.injection_style = "nrandompercell"
                species.bucket.num_particles_per_cell = self.number_per_cell
            else:
                raise Exception('Either nuniformpercell or nrandompercell must be specified')

            species.bucket.xmin = self.xmin
            species.bucket.xmax = self.xmax
            species.bucket.ymin = self.ymin
            species.bucket.ymax = self.ymax
            species.bucket.zmin = self.zmin
            species.bucket.zmax = self.zmax

            # --- Only constant density is supported at this time
            species.bucket.profile = "constant"
            species.bucket.density = self.density

            if self.vthx != 0. or self.vthy != 0. or self.vthz != 0.:
                species.bucket.momentum_distribution_type = "gaussian"
                species.bucket.ux_m = self.vxmean
                species.bucket.uy_m = self.vymean
                species.bucket.uz_m = self.vzmean
                species.bucket.ux_th = self.vthx
                species.bucket.uy_th = self.vthy
                species.bucket.uz_th = self.vthz
            else:
                species.bucket.momentum_distribution_type = "constant"
                species.bucket.ux = self.vxmean
                species.bucket.uy = self.vymean
                species.bucket.uz = self.vzmean

            if self.fill_in:
                pywarpx.warpx.do_plasma_injection = 1
                if not hasattr(pywarpx.warpx, 'injected_plasma_species'):
                    pywarpx.warpx.injected_plasma_species = []

                pywarpx.warpx.injected_plasma_species.append(species.species_number)
                pywarpx.warpx.num_injected_species = len(pywarpx.warpx.injected_plasma_species)


class ParticleList(PICMI_Base.PICMI_ParticleList):
    def init(self, **kw):

        if len(x) > 1:
            raise Exception('Only a single particle can be loaded')

        self.species.bucket.injection_style = "singleparticle"
        self.species.bucket.single_particle_pos = [self.x[0], self.y[0], self.z[0]]
        self.species.bucket.single_particle_vel = [self.ux[0]/c, self.uy[0]/c, self.uz[0]/c]
        self.species.bucket.single_particle_weight = self.weight

        # --- These need to be defined even though they are not used
        self.species.bucket.profile = "constant"
        self.species.bucket.density = 1
        self.species.bucket.momentum_distribution_type = 'constant'


class Grid(PICMI_Base.PICMI_Grid):
    def init(self, **kw):

        pywarpx.amr.n_cell = [self.nx, self.ny, self.nz]

        # Maximum allowable size of each subdomain in the problem domain;
        #    this is used to decompose the domain for parallel calculations.
        pywarpx.amr.max_grid_size = kw.get('max_grid_size', 32)

        # Maximum level in hierarchy (for now must be 0, i.e., one level in total)
        pywarpx.amr.max_level = kw.get('max_level', 0)

        # Geometry
        pywarpx.geometry.coord_sys = kw.get('coord_sys', 0)  # 0: Cartesian
        pywarpx.geometry.is_periodic = '%d %d %d'%(self.bcxmin=='periodic', self.bcymin=='periodic', self.bczmin=='periodic')  # Is periodic?
        pywarpx.geometry.prob_lo = [self.xmin, self.ymin, self.zmin]  # physical domain
        pywarpx.geometry.prob_hi = [self.xmax, self.ymax, self.zmax]

        if self.moving_window_velocity is not None and np.any(np.not_equal(self.moving_window_velocity, 0.)):
            pywarpx.warpx.do_moving_window = 1
            if self.moving_window_velocity[0] != 0.:
                pywarpx.warpx.moving_window_dir = 'x'
                pywarpx.warpx.moving_window_v = self.moving_window_velocity[0]/c  # in units of the speed of light
            if self.moving_window_velocity[1] != 0.:
                pywarpx.warpx.moving_window_dir = 'y'
                pywarpx.warpx.moving_window_v = self.moving_window_velocity[1]/c  # in units of the speed of light
            if self.moving_window_velocity[2] != 0.:
                pywarpx.warpx.moving_window_dir = 'z'
                pywarpx.warpx.moving_window_v = self.moving_window_velocity[2]/c  # in units of the speed of light

    def getmins(self, **kw):
        return np.array([pywarpx.warpx.getProbLo(0), pywarpx.warpx.getProbLo(1), pywarpx.warpx.getProbLo(2)])

    def getmaxs(self, **kw):
        return np.array([pywarpx.warpx.getProbHi(0), pywarpx.warpx.getProbHi(1), pywarpx.warpx.getProbHi(2)])

    def getxmin(self):
        return pywarpx.warpx.getProbLo(0)

    def getxmax(self):
        return pywarpx.warpx.getProbHi(0)

    def getymin(self):
        return pywarpx.warpx.getProbLo(1)

    def getymax(self):
        return pywarpx.warpx.getProbHi(1)

    def getzmin(self):
        return pywarpx.warpx.getProbLo(2)

    def getzmax(self):
        return pywarpx.warpx.getProbHi(2)


class EM_solver(PICMI_Base.PICMI_EM_solver):
    def init(self, **kw):

        if self.method is None:
            self.method = 'Yee'

        assert self.method in ['Yee'], Exception("Only 'Yee' FDTD is supported")

        if 'current_deposition_algo' in kw:
            pywarpx.algo.current_deposition = kw['current_deposition_algo']
        if 'charge_deposition_algo' in kw:
            pywarpx.algo.charge_deposition = kw['charge_deposition_algo']
        if 'field_gathering_algo' in kw:
            pywarpx.algo.field_gathering = kw['field_gathering_algo']
        if 'particle_pusher_algo' in kw:
            pywarpx.algo.particle_pusher = kw['particle_pusher_algo']

        pywarpx.interpolation.nox = self.norderx
        pywarpx.interpolation.noy = self.nordery
        pywarpx.interpolation.noz = self.norderz

class Simulation(PICMI_Base.PICMI_Simulation):
    def init(self, **kw):

        pywarpx.warpx.verbose = self.verbose
        pywarpx.warpx.cfl = self.timestep_over_cfl
        if self.timestep == 0.:
            pywarpx.warpx.cfl = self.timestep_over_cfl
        else:
            pywarpx.warpx.const_dt = self.timestep

        if 'plot_int' in kw:
            pywarpx.amr.plot_int = kw['plot_int']

        self.initialized = False

    def initialize(self, inputs_name=None):
        if not self.initialized:
            self.initialized = True
            pywarpx.warpx.init()

    def write_inputs(self, inputs_name='inputs'):
        kw = {}
        if self.max_step is not None:
            kw['max_step'] = self.max_step
        if self.max_time is not None:
            kw['stop_time'] = self.max_time
        pywarpx.warpx.write_inputs(inputs_name, **kw)

    def step(self, nsteps=None):
        self.initialize()
        if nsteps is None:
            if self.max_step is not None:
                nsteps = self.max_step
            else:
                nsteps = -1
        pywarpx.warpx.evolve(nsteps)

    def finalize(self):
        if self.initialized:
            self.initialized = False
            pywarpx.warpx.finalize()


class Gaussian_laser(PICMI_Base.PICMI_Gaussian_laser):
    def init(self, **kw):

        pywarpx.warpx.use_laser = 1
        pywarpx.laser.profile = "Gaussian"
        pywarpx.laser.wavelength = self.wavelength  # The wavelength of the laser (in meters)
        pywarpx.laser.e_max = self.E0  # Maximum amplitude of the laser field (in V/m)
        pywarpx.laser.polarization = [np.cos(self.pol_angle), np.sin(self.pol_angle), 0.]  # The main polarization vector
        pywarpx.laser.profile_waist = self.waist  # The waist of the laser (in meters)
        pywarpx.laser.profile_duration = self.duration  # The duration of the laser (in seconds)
        pywarpx.laser.profile_t_peak = (self.focal_position - self.z0)/c  # The time at which the laser reaches its peak (in seconds)


class Laser_antenna(PICMI_Base.PICMI_Laser_antenna):
    def init(self, **kw):

        pywarpx.laser.position = [self.antenna_x0, self.antenna_y0, self.antenna_z0]  # This point is on the laser plane
        pywarpx.laser.direction = [self.antenna_xvec, self.antenna_yvec, self.antenna_zvec]  # The plane normal direction
        pywarpx.laser.profile_focal_distance = self.laser.focal_position - self.antenna_z0  # Focal distance from the antenna (in meters)
