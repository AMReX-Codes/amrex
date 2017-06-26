"""Classes following the PICMI standard
"""
import numpy as np
from pywarpx import *

pi = 3.14159265358979323  # ratio of a circle's circumference to its diameter
euler = 0.57721566490153386  # Euler-Masceroni constant. Base of the natural logarithm.
amu = 1.660538921e-27  # Atomic Mass Unit [kg]
clight = 2.99792458e+8  # Speed of light in vacuum (exact) [m/s]
echarge = 1.602176565e-19  # Elementary charge [C]
emass = 9.10938291e-31  # Electron mass [kg]
mu0 = 4.e-7*pi  # Permeability of free space [kg.m/(s.s.A.A)=H/m=T.m/A]
eps0 = 1./(mu0*clight*clight)  # Permittivity of free space [F/m]
boltzmann = 1.3806488e-23  # Boltzmann's constant [J/K]
avogadro = 6.02214129e23  # Avogadro's Number [atoms/mole]
planck = 6.62606957e-34  # Planck's constant [J.s]

class Grid(object):
    """
    - `Grid`
      - **type**: *object*
      - `Nx=nx` - **type**: *integer* - "Number of cells along X (Nb nodes=nx+1)."
      - `Ny=ny` - **type**: *integer* - "Number of cells along Y (Nb nodes=ny+1)."
      - `Nr=nr` - **type**: *integer* - "Number of cells along R (Nb nodes=nr+1)."
      - `Nz=nz` - **type**: *integer* - "Number of cells along Z (Nb nodes=nz+1)."
      - `Nm=nm` - **type**: *integer* - "Number of azimuthal modes."
      - `Xmin=xmin` - **type**: *double* - "Position of first node along X."
      - `Xmax=xmax` - **type**: *double* - "Position of last node along X."
      - `Ymin=ymin` - **type**: *double* - "Position of first node along Y."
      - `Ymax=ymax` - **type**: *double* - "Position of last node along Y."
      - `Rmax=rmax` - **type**: *double* - "Position of last node along R."
      - `Zmin=zmin` - **type**: *double* - "Position of first node along Z."
      - `Zmax=zmax` - **type**: *double* - "Position of last node along Z."
      - `bcxmin` - **type**: *string* - "Boundary condition at min X: periodic/open/dirichlet/neumann."
      - `bcxmax` - **type**: *string* - "Boundary condition at max X: periodic/open/dirichlet/neumann."
      - `bcymin` - **type**: *string* - "Boundary condition at min Y: periodic/open/dirichlet/neumann."
      - `bcymax` - **type**: *string* - "Boundary condition at max Y: periodic/open/dirichlet/neumann."
      - `bcrmax` - **type**: *string* - "Boundary condition at max R: open/dirichlet/neumann."
      - `bczmin` - **type**: *string* - "Boundary condition at min Z: periodic/open/dirichlet/neumann."
      - `bczmax` - **type**: *string* - "Boundary condition at max Z: periodic/open/dirichlet/neumann."

      - max_grid_size
      - max_level
      - coord_sys
    """

    def __init__(self, nx=None, ny=None, nr=None, nz=None, nm=None,
                 xmin=None, xmax=None, ymin=None, ymax=None, rmax=None, zmin=None, zmax=None,
                 bcxmin=None, bcxmax=None, bcymin=None, bcymax=None, bcrmax=None, bczmin=None, bczmax=None,
                 **kw):
        self.nx = nx
        self.ny = ny
        self.nr = nr
        self.nz = nz
        self.nm = nm
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.rmax = rmax
        self.zmin = zmin
        self.zmax = zmax
        self.bcxmin = bcxmin
        self.bcxmax = bcxmax
        self.bcymin = bcymin
        self.bcymax = bcymax
        self.bcrmax = bcrmax
        self.bczmin = bczmin
        self.bczmax = bczmax

        amr.n_cell = '%d  %d  %d'%(nx, ny, nz)

        # Maximum allowable size of each subdomain in the problem domain; 
        #    this is used to decompose the domain for parallel calculations.
        amr.max_grid_size = kw.get('max_grid_size', 32)

        # Maximum level in hierarchy (for now must be 0, i.e., one level in total)
        amr.max_level = kw.get('max_level', 0)

        # Geometry
        geometry.coord_sys = kw.get('coord_sys', 0)  # 0: Cartesian
        geometry.is_periodic = '%d %d %d'%(bcxmin=='periodic', bcymin=='periodic', bczmin=='periodic')  # Is periodic?  
        geometry.prob_lo = '%7.0e %7.0e %7.0e'%(xmin, ymin, zmin)  # physical domain
        geometry.prob_hi = '%7.0e %7.0e %7.0e'%(xmax, ymax, zmax)


class EM_solver(object):
    Methods_list = ['Yee', 'CK', 'CKC', 'Lehe', 'PSTD', 'PSATD', 'GPSTD']
    def __init__(self, Method=None,
                 norderx=None, nordery=None, norderr=None, norderz=None,
                 l_nodal=None,
                 current_deposition_algo=None, charge_deposition_algo=None,
                 field_gathering_algo=None, particle_pusher_algo=None, **kw):

        assert Method is None or Method in EM_solver.Methods_list, Exception('Method has incorrect value')

        if current_deposition_algo is not None:
            algo.current_deposition = current_deposition_algo
        if charge_deposition_algo is not None:
            algo.charge_deposition = charge_deposition_algo
        if field_gathering_algo is not None:
            algo.field_gathering = field_gathering_algo
        if particle_pusher_algo is not None:
            algo.particle_pusher = particle_pusher_algo


class Particle(object):
    def __init__(self, Charge=None, charge=None, Q=None, q=None,
                 Mass=None, mass=None, M=None, m=None,
                 Symbol=None, symbol=None, S=None, s=None,
                 Name=None, name=None, N=None, n=None, **kw):
        # --- Accept multpiple names, but use 'charge', 'mass', 'symbol', 'name' internally.
        if Charge is not None: charge = Charge
        if Q is not None: charge = Q
        if q is not None: charge = q
        if Mass is not None: mass = Mass
        if M is not None: mass = M
        if m is not None: mass = m
        if Symbol is not None: symbol = Symbol
        if S is not None: symbol = S
        if s is not None: symbol = s
        if Name is not None: name = Name
        if N is not None: name = N
        if n is not None: name = n
        self.charge = charge
        self.mass = mass
        self.symbol = symbol

Electron = Particle(q=-echarge, m=emass, symbol='e-', name='Electron')
Positron = Particle(q=echarge, m=emass, symbol='e+', name='Positron')
Proton = Particle(q=echarge, m=1.6726231e-27, symbol='p', name='Proton')
AntiProton = Particle(q=-echarge, m=1.6726231e-27, symbol='p-', name='Antiproton')
Neutron = Particle(q=0. , m=1.6749286e-27, symbol='n', name='Neutron')
Muon = Particle(q=-echarge, m=1.883531475e-28, symbol='mu-', name='Muon')
Antimuon = Particle(q=echarge, m=1.883531475e-28, symbol='mu+', name='Antimuon')
Photon = Particle(q=0., m=0., symbol='gnu', name='Photon')
    

class Species(object):
    def __init__(self,
                  Type=None, type=None,
                  Name=None, name=None,
                  Sid=None, sid=None,
                  Charge_state=None, charge_state=None,
                  Charge=None, charge=None, Q=None, q=None,
                  Mass=None, mass=None, M=None, m=None,
                  Weight=None, weight=None, W=None, w=None, **kw):
        # --- Accept multpiple names, but use 'type', 'name', 'sid', 'charge_state', 'charge', 'mass', 'weight'
        if Type is not None: type = Type
        if Name is not None: name = Name
        if Sid is not None: sid = Sid
        if Charge_state is not None: charge_state = Charge_state
        if Charge is not None: charge = Charge
        if Q is not None: charge = Q
        if q is not None: charge = q
        if Mass is not None: mass = Mass
        if M is not None: mass = M
        if m is not None: mass = m
        if Weight is not None: weight = Weight
        if W is not None: weight = W
        if w is not None: weight = w
        self.type = type
        self.name = name
        self.sid = sid
        self.charg_state = charg_state
        self.charge = charge
        self.mass = mass
        self.weight = weight

        self.species_number = particles.nspecies
        particles.nspecies = particles.nspecies + 1
        particles.species_names = particles.species_names + ' ' + name

    def add_particles(self, n=None,
                      x=None, y=None, z=None,
                      ux=None, uy=None, uz=None, w=None,
                      unique_particles=None, **kw):
        pid = np.array([w]).T
        add_particles(self.species_number, x, y, z, ux, uy, uz, pid, unique_particles)


class Simulation(object):
    def __init__(self, plot_int=None, verbose=None, cfl=None):
        amr.plot_int = plot_int
        warpx.verbose = verbose
        warpx.cfl = cfl

        self.amrex = AMReX()
        self.amrex.init()
        warpx.init()

    def step(self, nsteps=-1):
        warpx.evolve(nsteps)
        
    def finalize(self):
        warpx.finalize()
        self.amrex.finalize()

        
