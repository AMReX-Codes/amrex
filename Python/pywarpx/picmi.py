"""Classes following the PICMI standard
"""
import picmistandard
import numpy as np
import pywarpx

codename = 'warpx'
picmistandard.register_codename(codename)

# --- Values from WarpXConst.H
c = 299792458.
ep0 = 8.854187817e-12
mu0 = 1.2566370614359173e-06
q_e = 1.602176462e-19
m_e = 9.10938291e-31
m_p = 1.6726231e-27


class Species(picmistandard.PICMI_Species):
    def init(self, kw):

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
        elif self.particle_type == 'H' and self.charge_state == 1:
            if self.charge is None: self.charge = 'q_e'
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


picmistandard.PICMI_MultiSpecies.Species_class = Species
class MultiSpecies(picmistandard.PICMI_MultiSpecies):
    pass


class GaussianBunchDistribution(picmistandard.PICMI_GaussianBunchDistribution):
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
        species.q_tot = self.n_physical_particles*charge

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

        # --- Note that WarpX takes gamma*beta as input
        if np.any(np.not_equal(self.velocity_divergence, 0.)):
            species.momentum_distribution_type = "radial_expansion"
            species.u_over_r = self.velocity_divergence[0]/c
            #species.u_over_y = self.velocity_divergence[1]/c
            #species.u_over_z = self.velocity_divergence[2]/c
        elif np.any(np.not_equal(self.rms_velocity, 0.)):
            species.momentum_distribution_type = "gaussian"
            species.ux_m = self.centroid_velocity[0]/c
            species.uy_m = self.centroid_velocity[1]/c
            species.uz_m = self.centroid_velocity[2]/c
            species.ux_th = self.rms_velocity[0]/c
            species.uy_th = self.rms_velocity[1]/c
            species.uz_th = self.rms_velocity[2]/c
        else:
            species.momentum_distribution_type = "constant"
            species.ux = self.centroid_velocity[0]/c
            species.uy = self.centroid_velocity[1]/c
            species.uz = self.centroid_velocity[2]/c


class UniformDistribution(picmistandard.PICMI_UniformDistribution):
    def initialize_inputs(self, species_number, layout, species):

        if isinstance(layout, GriddedLayout):
            # --- Note that the grid attribute of GriddedLayout is ignored
            species.injection_style = "nuniformpercell"
            species.num_particles_per_cell_each_dim = layout.n_macroparticle_per_cell
        elif isinstance(layout, PseudoRandomLayout):
            assert (layout.n_macroparticles_per_cell is not None), Exception('WarpX only supports n_macroparticles_per_cell for the PseudoRandomLayout with UniformDistribution')
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

        # --- Note that WarpX takes gamma*beta as input
        if np.any(np.not_equal(self.rms_velocity, 0.)):
            species.momentum_distribution_type = "gaussian"
            species.ux_m = self.directed_velocity[0]/c
            species.uy_m = self.directed_velocity[1]/c
            species.uz_m = self.directed_velocity[2]/c
            species.ux_th = self.rms_velocity[0]/c
            species.uy_th = self.rms_velocity[1]/c
            species.uz_th = self.rms_velocity[2]/c
        else:
            species.momentum_distribution_type = "constant"
            species.ux = self.directed_velocity[0]/c
            species.uy = self.directed_velocity[1]/c
            species.uz = self.directed_velocity[2]/c

        if self.fill_in:
            species.do_continuous_injection = 1


class AnalyticDistribution(picmistandard.PICMI_AnalyticDistribution):
    def initialize_inputs(self, species_number, layout, species):

        if isinstance(layout, GriddedLayout):
            # --- Note that the grid attribute of GriddedLayout is ignored
            species.injection_style = "nuniformpercell"
            species.num_particles_per_cell_each_dim = layout.n_macroparticle_per_cell
        elif isinstance(layout, PseudoRandomLayout):
            assert (layout.n_macroparticles_per_cell is not None), Exception('WarpX only supports n_macroparticles_per_cell for the PseudoRandomLayout with UniformDistribution')
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
        species.profile = "parse_density_function"
        species.__setattr__('density_function(x,y,z)', self.density_expression)

        for k,v in self.user_defined_kw.items():
            setattr(pywarpx.my_constants, k, v)

        # --- Note that WarpX takes gamma*beta as input
        if np.any(np.not_equal(self.rms_velocity, 0.)):
            species.momentum_distribution_type = "gaussian"
            species.ux_m = self.directed_velocity[0]/c
            species.uy_m = self.directed_velocity[1]/c
            species.uz_m = self.directed_velocity[2]/c
            species.ux_th = self.rms_velocity[0]/c
            species.uy_th = self.rms_velocity[1]/c
            species.uz_th = self.rms_velocity[2]/c
        else:
            species.momentum_distribution_type = "constant"
            species.ux = self.directed_velocity[0]/c
            species.uy = self.directed_velocity[1]/c
            species.uz = self.directed_velocity[2]/c

        if self.fill_in:
            species.do_continuous_injection = 1


class ParticleListDistribution(picmistandard.PICMI_ParticleListDistribution):
    def init(self, kw):

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


class ParticleDistributionPlanarInjector(picmistandard.PICMI_ParticleDistributionPlanarInjector):
    pass


class GriddedLayout(picmistandard.PICMI_GriddedLayout):
    pass


class PseudoRandomLayout(picmistandard.PICMI_PseudoRandomLayout):
    def init(self, kw):
        if self.seed is not None:
            print('Warning: WarpX does not support specifying the random number seed')


class BinomialSmoother(picmistandard.PICMI_BinomialSmoother):
    pass


class CylindricalGrid(picmistandard.PICMI_CylindricalGrid):
    """This assumes that WarpX was compiled with USE_RZ = TRUE
    """
    def init(self, kw):
        self.max_grid_size = kw.pop('warpx_max_grid_size', 32)
        self.blocking_factor = kw.pop('warpx_blocking_factor', None)

    def initialize_inputs(self):
        pywarpx.amr.n_cell = self.number_of_cells

        # Maximum allowable size of each subdomain in the problem domain;
        #    this is used to decompose the domain for parallel calculations.
        pywarpx.amr.max_grid_size = self.max_grid_size

        assert self.lower_bound[0] >= 0., Exception('Lower radial boundary must be >= 0.')
        assert self.bc_rmin != 'periodic' and self.bc_rmax != 'periodic', Exception('Radial boundaries can not be periodic')
        assert self.n_azimuthal_modes is None or self.n_azimuthal_modes == 1, Exception('Only one azimuthal mode supported')

        # Geometry
        pywarpx.geometry.coord_sys = 1  # RZ
        pywarpx.geometry.is_periodic = '0 %d'%(self.bc_zmin=='periodic')  # Is periodic?
        pywarpx.geometry.prob_lo = self.lower_bound  # physical domain
        pywarpx.geometry.prob_hi = self.upper_bound

        if self.moving_window_velocity is not None and np.any(np.not_equal(self.moving_window_velocity, 0.)):
            pywarpx.warpx.do_moving_window = 1
            if self.moving_window_velocity[0] != 0.:
                raise Exception('In cylindrical coordinates, a moving window in r can not be done')
            if self.moving_window_velocity[1] != 0.:
                pywarpx.warpx.moving_window_dir = 'z'
                pywarpx.warpx.moving_window_v = self.moving_window_velocity[1]/c  # in units of the speed of light

        if self.refined_regions:
            assert len(self.refined_regions) == 1, Exception('WarpX only supports one refined region.')
            assert self.refined_regions[0][0] == 1, Exception('The one refined region can only be level 1')
            pywarpx.amr.max_level = 1
            pywarpx.warpx.fine_tag_lo = self.refined_regions[0][1]
            pywarpx.warpx.fine_tag_hi = self.refined_regions[0][2]
            # The refinement_factor is ignored (assumed to be [2,2])
        else:
            pywarpx.amr.max_level = 0


class Cartesian2DGrid(picmistandard.PICMI_Cartesian2DGrid):
    def init(self, kw):
        self.max_grid_size = kw.pop('warpx_max_grid_size', 32)
        self.blocking_factor = kw.pop('warpx_blocking_factor', None)

    def initialize_inputs(self):
        pywarpx.amr.n_cell = self.number_of_cells

        # Maximum allowable size of each subdomain in the problem domain;
        #    this is used to decompose the domain for parallel calculations.
        pywarpx.amr.max_grid_size = self.max_grid_size

        # Geometry
        pywarpx.geometry.coord_sys = 0  # Cartesian
        pywarpx.geometry.is_periodic = '%d %d'%(self.bc_xmin=='periodic', self.bc_ymin=='periodic')  # Is periodic?
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

        if self.refined_regions:
            assert len(self.refined_regions) == 1, Exception('WarpX only supports one refined region.')
            assert self.refined_regions[0][0] == 1, Exception('The one refined region can only be level 1')
            pywarpx.amr.max_level = 1
            pywarpx.warpx.fine_tag_lo = self.refined_regions[0][1]
            pywarpx.warpx.fine_tag_hi = self.refined_regions[0][2]
            # The refinement_factor is ignored (assumed to be [2,2])
        else:
            pywarpx.amr.max_level = 0


class Cartesian3DGrid(picmistandard.PICMI_Cartesian3DGrid):
    def init(self, kw):
        self.max_grid_size = kw.pop('warpx_max_grid_size', 32)
        self.blocking_factor = kw.pop('warpx_blocking_factor', None)

    def initialize_inputs(self):
        pywarpx.amr.n_cell = self.number_of_cells

        # Maximum allowable size of each subdomain in the problem domain;
        #    this is used to decompose the domain for parallel calculations.
        pywarpx.amr.max_grid_size = self.max_grid_size

        pywarpx.amr.blocking_factor = self.blocking_factor

        # Geometry
        pywarpx.geometry.coord_sys = 0  # Cartesian
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

        if self.refined_regions:
            assert len(self.refined_regions) == 1, Exception('WarpX only supports one refined region.')
            assert self.refined_regions[0][0] == 1, Exception('The one refined region can only be level 1')
            pywarpx.amr.max_level = 1
            pywarpx.warpx.fine_tag_lo = self.refined_regions[0][1]
            pywarpx.warpx.fine_tag_hi = self.refined_regions[0][2]
            # The refinement_factor is ignored (assumed to be [2,2,2])
        else:
            pywarpx.amr.max_level = 0

class ElectromagneticSolver(picmistandard.PICMI_ElectromagneticSolver):
    def init(self, kw):
        assert self.method is None or self.method in ['Yee', 'CKC'], Exception("Only 'Yee' and 'CKC' FDTD are supported")

        self.do_pml = kw.pop('warpx_do_pml', None)
        self.pml_ncell = kw.pop('warpx_pml_ncell', None)

    def initialize_inputs(self):

        self.grid.initialize_inputs()

        pywarpx.warpx.do_pml = self.do_pml
        pywarpx.warpx.pml_ncell = self.pml_ncell

        # --- Same method names are used, though mapped to lower case.
        pywarpx.warpx.maxwell_fdtd_solver = self.method

        if self.cfl is not None:
            pywarpx.warpx.cfl = self.cfl

        if self.stencil_order is not None:
            pywarpx.interpolation.nox = self.stencil_order[0]
            pywarpx.interpolation.noy = self.stencil_order[1]
            pywarpx.interpolation.noz = self.stencil_order[2]


class ElectrostaticSolver(picmistandard.PICMI_ElectrostaticSolver):
    def initialize_inputs(self):
        pass


class GaussianLaser(picmistandard.PICMI_GaussianLaser):
    def initialize_inputs(self):
        self.laser_number = pywarpx.lasers.nlasers + 1
        self.name = 'laser{}'.format(self.laser_number)

        self.laser = pywarpx.Lasers.newlaser(self.name)

        self.laser.profile = "Gaussian"
        self.laser.wavelength = self.wavelength  # The wavelength of the laser (in meters)
        self.laser.e_max = self.E0  # Maximum amplitude of the laser field (in V/m)
        self.laser.polarization = [np.cos(self.polarization_angle), np.sin(self.polarization_angle), 0.]  # The main polarization vector
        self.laser.profile_waist = self.waist  # The waist of the laser (in meters)
        self.laser.profile_duration = self.duration  # The duration of the laser (in seconds)
        self.laser.zeta = self.zeta
        self.laser.beta = self.beta
        self.laser.phi2 = self.phi2


class LaserAntenna(picmistandard.PICMI_LaserAntenna):
    def initialize_inputs(self, laser):
        laser.laser.position = self.position  # This point is on the laser plane
        laser.laser.direction = self.normal_vector  # The plane normal direction
        laser.laser.profile_focal_distance = laser.focal_position[2] - self.position[2]  # Focal distance from the antenna (in meters)
        laser.laser.profile_t_peak = (self.position[2] - laser.centroid_position[2])/c  # The time at which the laser reaches its peak (in seconds)


class Simulation(picmistandard.PICMI_Simulation):
    def init(self, kw):

        self.plot_int = kw.pop('warpx_plot_int', None)
        self.plot_file = kw.pop('warpx_plot_file', None)
        self.current_deposition_algo = kw.pop('warpx_current_deposition_algo', None)
        self.charge_deposition_algo = kw.pop('warpx_charge_deposition_algo', None)
        self.field_gathering_algo = kw.pop('warpx_field_gathering_algo', None)
        self.particle_pusher_algo = kw.pop('warpx_particle_pusher_algo', None)
        self.use_filter = kw.pop('warpx_use_filter', None)
        self.serialize_ics = kw.pop('warpx_serialize_ics', None)
        self.do_dynamic_scheduling = kw.pop('warpx_do_dynamic_scheduling', None)
        self.load_balance_int = kw.pop('warpx_load_balance_int', None)
        self.load_balance_with_sfc = kw.pop('warpx_load_balance_with_sfc', None)

        self.inputs_initialized = False
        self.warpx_initialized = False

    def initialize_inputs(self):
        if self.inputs_initialized:
            return

        self.inputs_initialized = True

        pywarpx.warpx.verbose = self.verbose
        if self.time_step_size is not None:
            pywarpx.warpx.const_dt = self.timestep

        if self.gamma_boost is not None:
            pywarpx.warpx.gamma_boost = self.gamma_boost
            pywarpx.warpx.boost_direction = 'z'

        pywarpx.amr.plot_int = self.plot_int
        pywarpx.amr.plot_file = self.plot_file
        pywarpx.algo.current_deposition = self.current_deposition_algo
        pywarpx.algo.charge_deposition = self.charge_deposition_algo
        pywarpx.algo.field_gathering = self.field_gathering_algo
        pywarpx.algo.particle_pusher = self.particle_pusher_algo

        pywarpx.warpx.use_filter = self.use_filter
        pywarpx.warpx.serialize_ics = self.serialize_ics

        pywarpx.warpx.do_dynamic_scheduling = self.do_dynamic_scheduling
        pywarpx.warpx.load_balance_int = self.load_balance_int
        pywarpx.warpx.load_balance_with_sfc = self.load_balance_with_sfc

        particle_shape = self.particle_shape
        for s in self.species:
            if s.particle_shape is not None:
                assert particle_shape is None or particle_shape == s.particle_shape, Exception('WarpX only supports one particle shape for all species')
                # --- If this was set for any species, use that value.
                particle_shape = s.particle_shape

        if particle_shape is not None:
            if isinstance(particle_shape, str):
                interpolation_order = {'NGP':0, 'linear':1, 'quadratic':2, 'cubic':3}[particle_shape]
            else:
                interpolation_order = particle_shape
            pywarpx.interpolation.nox = interpolation_order
            pywarpx.interpolation.noy = interpolation_order
            pywarpx.interpolation.noz = interpolation_order

        self.solver.initialize_inputs()

        for i in range(len(self.species)):
            assert not self.initialize_self_fields[i], Exception('WarpX does not support initializing self fields')
            self.species[i].initialize_inputs(self.layouts[i])

        for i in range(len(self.lasers)):
            self.lasers[i].initialize_inputs()
            self.laser_injection_methods[i].initialize_inputs(self.lasers[i])

        for diagnostic in self.diagnostics:
            diagnostic.initialize_inputs()

    def initialize_warpx(self):
        if self.warpx_initialized:
            return

        self.warpx_initialized = True
        pywarpx.warpx.init()

    def write_input_file(self, file_name='inputs'):
        self.initialize_inputs()
        kw = {}
        if self.max_steps is not None:
            kw['max_step'] = self.max_steps
        if self.max_time is not None:
            kw['stop_time'] = self.max_time
        pywarpx.warpx.write_inputs(file_name, **kw)

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


# ----------------------------
# Simulation frame diagnostics
# ----------------------------


class FieldDiagnostic(picmistandard.PICMI_FieldDiagnostic):
    def init(self, kw):

        self.plot_raw_fields = kw.pop('warpx_plot_raw_fields', None)
        self.plot_raw_fields_guards = kw.pop('warpx_plot_raw_fields_guards', None)
        self.plot_finepatch = kw.pop('warpx_plot_finepatch', None)
        self.plot_crsepatch = kw.pop('warpx_plot_crsepatch', None)

    def initialize_inputs(self):
        # --- For now, the period must be the same as plot_int if set
        pywarpx.amr.check_consistency('plot_int', self.period, 'The period must be the same for all simulation frame diagnostics')
        pywarpx.amr.plot_int = self.period

        if 'rho' in self.data_list:
            pywarpx.warpx.plot_rho = 1
        if 'dive' in self.data_list:
            pywarpx.warpx.plot_dive = 1
        if 'divb' in self.data_list:
            pywarpx.warpx.plot_divb = 1
        if 'F' in self.data_list:
            pywarpx.warpx.plot_F = 1
        if 'proc_number' in self.data_list:
            pywarpx.warpx.plot_proc_number = 1

        pywarpx.warpx.plot_raw_fields = self.plot_raw_fields
        pywarpx.warpx.plot_raw_fields_guards = self.plot_raw_fields_guards

        pywarpx.amr.check_consistency('plot_finepatch', self.plot_finepatch, 'The fine patch flag must be the same for all simulation frame field diagnostics')
        pywarpx.amr.check_consistency('plot_crsepatch', self.plot_crsepatch, 'The coarse patch flag must be the same for all simulation frame field diagnostics')
        pywarpx.warpx.plot_finepatch = self.plot_finepatch
        pywarpx.warpx.plot_crsepatch = self.plot_crsepatch

class ElectrostaticFieldDiagnostic(picmistandard.PICMI_ElectrostaticFieldDiagnostic):
    def initialize_inputs(self):
        # --- For now, the period must be the same as plot_int if set
        pywarpx.amr.check_consistency('plot_int', self.period, 'The period must be the same for all simulation frame diagnostics')
        pywarpx.amr.plot_int = self.period


class ParticleDiagnostic(picmistandard.PICMI_ParticleDiagnostic):
    def initialize_inputs(self):
        # --- For now, the period must be the same as plot_int if set
        pywarpx.amr.check_consistency('plot_int', self.period, 'The period must be the same for all simulation frame diagnostics')
        pywarpx.amr.plot_int = self.period

        if 'part_per_cell' in self.data_list:
            pywarpx.warpx.plot_part_per_cell = 1
        if 'part_per_grid' in self.data_list:
            pywarpx.warpx.plot_part_per_grid = 1
        if 'part_per_proc' in self.data_list:
            pywarpx.warpx.plot_part_per_proc = 1


# ----------------------------
# Lab frame diagnostics
# ----------------------------


class LabFrameFieldDiagnostic(picmistandard.PICMI_LabFrameFieldDiagnostic):
    def initialize_inputs(self):

        pywarpx.warpx.check_consistency('num_snapshots_lab', self.num_snapshots, 'The number of snapshots must be the same in all lab frame diagnostics')
        pywarpx.warpx.check_consistency('dt_snapshots_lab', self.dt_snapshots, 'The time between snapshots must be the same in all lab frame diagnostics')

        pywarpx.warpx.do_boosted_frame_diagnostic = 1
        pywarpx.warpx.num_snapshots_lab = self.num_snapshots
        pywarpx.warpx.dt_snapshots_lab = self.dt_snapshots
        pywarpx.warpx.do_boosted_frame_fields = 1


class LabFrameParticleDiagnostic(picmistandard.PICMI_LabFrameParticleDiagnostic):
    def initialize_inputs(self):

        pywarpx.warpx.check_consistency('num_snapshots_lab', self.num_snapshots, 'The number of snapshots must be the same in all lab frame diagnostics')
        pywarpx.warpx.check_consistency('dt_snapshots_lab', self.dt_snapshots, 'The time between snapshots must be the same in all lab frame diagnostics')

        pywarpx.warpx.do_boosted_frame_diagnostic = 1
        pywarpx.warpx.num_snapshots_lab = self.num_snapshots
        pywarpx.warpx.dt_snapshots_lab = self.dt_snapshots
        pywarpx.warpx.do_boosted_frame_particles = 1
