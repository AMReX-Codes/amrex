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

        if self.particle_type == 'electron':
            if self.charge is None: self.charge = '-q_e'
            if self.mass is None: self.mass = 'm_e'
        elif self.particle_type == 'positron':
            if self.charge is None: self.charge = 'q_e'
            if self.mass is None: self.mass = 'm_e'
        elif self.particle_type == 'proton':
            if self.charge is None: self.charge = 'q_e'
            if self.mass is None: self.mass = 'm_p'
        elif self.particle_type == 'anti-proton':
            if self.charge is None: self.charge = '-q_e'
            if self.mass is None: self.mass = 'm_p'

    def initialize_inputs(self, layout):
        self.species_number = pywarpx.particles.nspecies
        pywarpx.particles.nspecies += 1

        if self.name is None:
            self.name = 'species{}'.format(self.species_number)

        if pywarpx.particles.species_names is None:
            pywarpx.particles.species_names = self.name
        else:
            pywarpx.particles.species_names += ' ' + self.name

        self.species = pywarpx.Bucket.Bucket(self.name, mass=self.mass, charge=self.charge, injection_style = 'python')
        pywarpx.Particles.particles_list.append(self.species)

        if self.initial_distribution is not None:
            self.initial_distribution.initialize_inputs(self.species_number, layout, self.species)


PICMI_Base.PICMI_MultiSpecies.Species_class = Species
class MultiSpecies(PICMI_Base.PICMI_MultiSpecies):
    pass


class GaussianBunchDistribution(PICMI_Base.PICMI_GaussianBunchDistribution):
    def init(self, **kw):
        if self.seed is not None:
            print('Warning: WarpX does not support specifying the random number seed')

    def initialize_inputs(self, species_number, layout, species):
        species.injection_style = "gaussian_beam"
        species.x_m = self.centroid_position[0]
        species.y_m = self.centroid_position[1]
        species.z_m = self.centroid_position[2]
        species.x_rms = self.rms_bunch_size[0]
        species.y_rms = self.rms_bunch_size[1]
        species.z_rms = self.rms_bunch_size[2]

        # --- Only PseudoRandomLayout is supported
        species.npart = layout.n_macroparticles

        # --- Calculate the total charge. Note that charge might be a string instead of a number.
        charge = species.charge
        if charge == 'q_e' or charge == '+q_e':
            charge = q_e
        elif charge == '-q_e':
            charge = -q_e
        species.q_tot = self.number_real_particles*charge

        # --- These need to be defined even though they are not used
        species.profile = "constant"
        species.density = 1

        # --- The PICMI standard doesn't yet have a way of specifying these values.
        # --- They should default to the size of the domain. They are not typically
        # --- necessary though since any particles outside the domain are rejected.
        #species.xmin
        #species.xmax
        #species.ymin
        #species.ymax
        #species.zmin
        #species.zmax

        if np.any(np.not_equal(self.velocity_divergence, 0.)):
            species.momentum_distribution_type = "radial_expansion"
            species.u_over_r = self.velocity_divergence[0]
            #species.u_over_y = self.velocity_divergence[1]
            #species.u_over_z = self.velocity_divergence[2]
        elif np.any(np.not_equal(self.rms_velocity, 0.)):
            species.momentum_distribution_type = "gaussian"
            species.ux_m = self.centroid_velocity[0]
            species.uy_m = self.centroid_velocity[1]
            species.uz_m = self.centroid_velocity[2]
            species.ux_th = self.rms_velocity[0]
            species.uy_th = self.rms_velocity[1]
            species.uz_th = self.rms_velocity[2]
        else:
            species.momentum_distribution_type = "constant"
            species.ux = self.centroid_velocity[0]
            species.uy = self.centroid_velocity[1]
            species.uz = self.centroid_velocity[2]


class UniformDistribution(PICMI_Base.PICMI_UniformDistribution):

    def initialize_inputs(self, species_number, layout, species):

        if isinstance(layout, GriddedLayout):
            species.injection_style = "nuniformpercell"
            species.num_particles_per_cell_each_dim = layout.n_macroparticle_per_cell
        elif isinstance(layout, PseudoRandomLayout):
            assert (layout.n_macroparticles_per_cell is not None), Exception('WarpX only supports n_macroparticles_per_cell for the GriddedLayout with UniformDistribution')
            species.injection_style = "nrandompercell"
            species.num_particles_per_cell = layout.n_macroparticles_per_cell
        else:
            raise Exception('WarpX does not support the specified layout for UniformDistribution')

        species.xmin = self.lower_bound[0]
        species.xmax = self.upper_bound[0]
        species.ymin = self.lower_bound[1]
        species.ymax = self.upper_bound[1]
        species.zmin = self.lower_bound[2]
        species.zmax = self.upper_bound[2]

        # --- Only constant density is supported at this time
        species.profile = "constant"
        species.density = self.density

        if np.any(np.not_equal(self.rms_velocity, 0.)):
            species.momentum_distribution_type = "gaussian"
            species.ux_m = self.directed_velocity[0]
            species.uy_m = self.directed_velocity[1]
            species.uz_m = self.directed_velocity[2]
            species.ux_th = self.rms_velocity[0]
            species.uy_th = self.rms_velocity[1]
            species.uz_th = self.rms_velocity[2]
        else:
            species.momentum_distribution_type = "constant"
            species.ux = self.directed_velocity[0]
            species.uy = self.directed_velocity[1]
            species.uz = self.directed_velocity[2]

        if self.fill_in:
            pywarpx.warpx.do_plasma_injection = 1
            if not hasattr(pywarpx.warpx, 'injected_plasma_species'):
                pywarpx.warpx.injected_plasma_species = []

            pywarpx.warpx.injected_plasma_species.append(species_number)
            pywarpx.warpx.num_injected_species = len(pywarpx.warpx.injected_plasma_species)


class AnalyticDistribution(PICMI_Base.PICMI_AnalyticDistribution):

    def initialize_inputs(self, species_number, layout, species):
        raise Exception('WarpX does not support AnalyticDistribution')


class ParticleList(PICMI_Base.PICMI_ParticleList):
    def init(self, **kw):

        if len(x) > 1:
            raise Exception('Only a single particle can be loaded')

    def initialize_inputs(self, species_number, layout, species):

        species.injection_style = "singleparticle"
        species.single_particle_pos = [self.x[0], self.y[0], self.z[0]]
        species.single_particle_vel = [self.ux[0]/c, self.uy[0]/c, self.uz[0]/c]
        species.single_particle_weight = self.weight

        # --- These need to be defined even though they are not used
        species.profile = "constant"
        species.density = 1
        species.momentum_distribution_type = 'constant'


class ParticleDistributionPlanarInjector(PICMI_Base.PICMI_ParticleDistributionPlanarInjector):
    pass


class GriddedLayout(PICMI_Base.PICMI_GriddedLayout):
    pass


class PseudoRandomLayout(PICMI_Base.PICMI_PseudoRandomLayout):
    pass


class BinomialSmoother(PICMI_Base.PICMI_BinomialSmoother):
    pass


class CylindricalGrid(PICMI_Base.PICMI_CylindricalGrid):
    def init(self, **kw):
        raise Exception('WarpX does not support CylindricalGrid')


class Cartesian2DGrid(PICMI_Base.PICMI_Cartesian2DGrid):
    def init(self, **kw):
        self.max_grid_size = kw.get('max_grid_size', 32)
        self.max_level = kw.get('max_level', 0)
        self.coord_sys = kw.get('coord_sys', 0)

    def initialize_inputs(self):
        pywarpx.amr.n_cell = self.number_of_cells

        # Maximum allowable size of each subdomain in the problem domain;
        #    this is used to decompose the domain for parallel calculations.
        pywarpx.amr.max_grid_size = self.max_grid_size

        # Maximum level in hierarchy (for now must be 0, i.e., one level in total)
        pywarpx.amr.max_level = self.max_level

        # Geometry
        pywarpx.geometry.coord_sys = self.coord_sys
        pywarpx.geometry.is_periodic = '%d %d %d'%(self.bc_xmin=='periodic', self.bc_ymin=='periodic')  # Is periodic?
        pywarpx.geometry.prob_lo = self.lower_bound  # physical domain
        pywarpx.geometry.prob_hi = self.upper_bound

        if self.moving_window_velocity is not None and np.any(np.not_equal(self.moving_window_velocity, 0.)):
            pywarpx.warpx.do_moving_window = 1
            if self.moving_window_velocity[0] != 0.:
                pywarpx.warpx.moving_window_dir = 'x'
                pywarpx.warpx.moving_window_v = self.moving_window_velocity[0]/c  # in units of the speed of light
            if self.moving_window_velocity[1] != 0.:
                pywarpx.warpx.moving_window_dir = 'y'
                pywarpx.warpx.moving_window_v = self.moving_window_velocity[1]/c  # in units of the speed of light


class Cartesian3DGrid(PICMI_Base.PICMI_Cartesian3DGrid):
    def init(self, **kw):
        self.max_grid_size = kw.get('max_grid_size', 32)
        self.max_level = kw.get('max_level', 0)
        self.coord_sys = kw.get('coord_sys', 0)

    def initialize_inputs(self):
        pywarpx.amr.n_cell = self.number_of_cells

        # Maximum allowable size of each subdomain in the problem domain;
        #    this is used to decompose the domain for parallel calculations.
        pywarpx.amr.max_grid_size = self.max_grid_size

        # Maximum level in hierarchy (for now must be 0, i.e., one level in total)
        pywarpx.amr.max_level = self.max_level

        # Geometry
        pywarpx.geometry.coord_sys = self.coord_sys
        pywarpx.geometry.is_periodic = '%d %d %d'%(self.bc_xmin=='periodic', self.bc_ymin=='periodic', self.bc_zmin=='periodic')  # Is periodic?
        pywarpx.geometry.prob_lo = self.lower_bound  # physical domain
        pywarpx.geometry.prob_hi = self.upper_bound

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


class ElectromagneticSolver(PICMI_Base.PICMI_ElectromagneticSolver):
    def init(self, **kw):
        assert self.method is None or self.method in ['Yee'], Exception("Only 'Yee' FDTD is supported")

    def initialize_inputs(self):

        self.grid.initialize_inputs()

        if self.cfl is not None:
            pywarpx.warpx.cfl = self.cfl

        if self.stencil_order is not None:
            pywarpx.interpolation.nox = self.stencil_order[0]
            pywarpx.interpolation.noy = self.stencil_order[1]
            pywarpx.interpolation.noz = self.stencil_order[2]


class Electrostatic_solver(PICMI_Base.PICMI_Electrostatic_solver):
    def initialize_inputs(self):
        pass


class GaussianLaser(PICMI_Base.PICMI_GaussianLaser):

    def initialize_inputs(self):
        pywarpx.warpx.use_laser = 1
        pywarpx.laser.profile = "Gaussian"
        pywarpx.laser.wavelength = self.wavelength  # The wavelength of the laser (in meters)
        pywarpx.laser.e_max = self.E0  # Maximum amplitude of the laser field (in V/m)
        pywarpx.laser.polarization = [np.cos(self.polarization_angle), np.sin(self.polarization_angle), 0.]  # The main polarization vector
        pywarpx.laser.profile_waist = self.waist  # The waist of the laser (in meters)
        pywarpx.laser.profile_duration = self.duration  # The duration of the laser (in seconds)


class LaserAntenna(PICMI_Base.PICMI_LaserAntenna):

    def initialize_inputs(self, laser):
        pywarpx.laser.position = self.position  # This point is on the laser plane
        pywarpx.laser.direction = self.normal_vector  # The plane normal direction
        pywarpx.laser.profile_focal_distance = laser.focal_position[2] - self.position[2]  # Focal distance from the antenna (in meters)
        pywarpx.laser.profile_t_peak = (self.position[2] - laser.centroid_position[2])/c  # The time at which the laser reaches its peak (in seconds)


class Simulation(PICMI_Base.PICMI_Simulation):
    def init(self, **kw):

        self.plot_int = kw.get('plot_int', None)
        self.current_deposition_algo = kw.get('current_deposition_algo', None)
        self.charge_deposition_algo = kw.get('charge_deposition_algo', None)
        self.field_gathering_algo = kw.get('field_gathering_algo', None)
        self.particle_pusher_algo = kw.get('particle_pusher_algo', None)

        self.inputs_initialized = False
        self.warpx_initialized = False

    def initialize_inputs(self):
        if self.inputs_initialized:
            return

        self.inputs_initialized = True

        pywarpx.warpx.verbose = self.verbose
        if self.time_step_size is not None:
            pywarpx.warpx.const_dt = self.timestep

        pywarpx.amr.plot_int = self.plot_int
        pywarpx.algo.current_deposition = self.current_deposition_algo
        pywarpx.algo.charge_deposition = self.charge_deposition_algo
        pywarpx.algo.field_gathering = self.field_gathering_algo
        pywarpx.algo.particle_pusher = self.particle_pusher_algo

        self.solver.initialize_inputs()

        for i in range(len(self.species)):
            assert self.calculate_self_fields[i], Exception('WarpX does not support species without self fields')
            self.species[i].initialize_inputs(self.layouts[i])
            
        for i in range(len(self.lasers)):
            self.lasers[i].initialize_inputs()
            self.laser_injection_methods[i].initialize_inputs(self.lasers[i])

    def initialize_warpx(self, inputs_name=None):
        if self.warpx_initialized:
            return

        self.warpx_initialized = True
        pywarpx.warpx.init()

    def write_input_file(self, inputs_name='inputs'):
        self.initialize_inputs()
        kw = {}
        if self.max_steps is not None:
            kw['max_step'] = self.max_steps
        if self.max_time is not None:
            kw['stop_time'] = self.max_time
        pywarpx.warpx.write_inputs(inputs_name, **kw)

    def step(self, nsteps=None):
        self.initialize_inputs()
        self.initialize_warpx()
        if nsteps is None:
            if self.max_steps is not None:
                nsteps = self.max_steps
            else:
                nsteps = -1
        pywarpx.warpx.evolve(nsteps)

    def finalize(self):
        if self.warpx_initialized:
            self.warpx_initialized = False
            pywarpx.warpx.finalize()
