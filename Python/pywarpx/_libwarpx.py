# --- This defines the wrapper functions that directly call the underlying compiled routines
import os
import sys
import ctypes
from ctypes.util import find_library as _find_library
import numpy as np
from numpy.ctypeslib import ndpointer as _ndpointer

from .Geometry import geometry

# --- Is there a better way of handling constants?
clight = 2.99792458e+8 # m/s

def _get_package_root():
    '''
    Get the path to the installation location (where libwarpx.so would be installed).
    '''
    cur = os.path.abspath(__file__)
    while True:
        name = os.path.basename(cur)
        if name == 'pywarpx':
            return cur
        elif not name:
            return ''
        cur = os.path.dirname(cur)

# --- Use geometry to determine whether to import the 2D or 3D version.
# --- This assumes that the input is setup before this module is imported,
# --- which should normally be the case.
# --- Default to 3D if geometry is not setup yet.
try:
    _prob_lo = geometry.prob_lo
except AttributeError:
    geometry_dim = 3
else:
    geometry_dim = len(_prob_lo)
    del _prob_lo

_libc = ctypes.CDLL(_find_library('c'))

try:
    libwarpx = ctypes.CDLL(os.path.join(_get_package_root(), "libwarpx%dd.so"%geometry_dim))
except OSError:
    raise Exception('libwarpx%dd.so was not installed. It can be installed by running "make" in the Python directory of WarpX'%geometry_dim)

dim = libwarpx.warpx_SpaceDim()

# our particle data type
_p_struct = [(d, 'f8') for d in 'xyz'[:dim]] + [('id', 'i4'), ('cpu', 'i4')]
_p_dtype = np.dtype(_p_struct, align=True)

_numpy_to_ctypes = {}
_numpy_to_ctypes['f8'] = ctypes.c_double
_numpy_to_ctypes['i4'] = ctypes.c_int

class Particle(ctypes.Structure):
    _fields_ = [(v[0], _numpy_to_ctypes[v[1]]) for v in _p_struct]


# some useful typenames
_LP_particle_p = ctypes.POINTER(ctypes.POINTER(Particle))
_LP_c_int = ctypes.POINTER(ctypes.c_int)
_LP_c_void_p = ctypes.POINTER(ctypes.c_void_p)
_LP_c_double = ctypes.POINTER(ctypes.c_double)
_LP_LP_c_double = ctypes.POINTER(_LP_c_double)
_LP_c_char = ctypes.POINTER(ctypes.c_char)
_LP_LP_c_char = ctypes.POINTER(_LP_c_char)

# this is a function for converting a ctypes pointer to a numpy array
def _array1d_from_pointer(pointer, dtype, size):
    if sys.version_info.major >= 3:
        # from where do I import these? this might only work for CPython...
        #PyBuf_READ  = 0x100
        PyBUF_WRITE = 0x200
        buffer_from_memory = ctypes.pythonapi.PyMemoryView_FromMemory
        buffer_from_memory.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int)
        buffer_from_memory.restype = ctypes.py_object
        buf = buffer_from_memory(pointer, dtype.itemsize*size, PyBUF_WRITE)
    else:
        buffer_from_memory = ctypes.pythonapi.PyBuffer_FromReadWriteMemory
        buffer_from_memory.restype = ctypes.py_object
        buf = buffer_from_memory(pointer, dtype.itemsize*size)
    return np.frombuffer(buf, dtype=dtype, count=size)


# set the arg and return types of the wrapped functions
libwarpx.amrex_init.argtypes = (ctypes.c_int, _LP_LP_c_char)
libwarpx.warpx_getParticleStructs.restype = _LP_particle_p
libwarpx.warpx_getParticleArrays.restype = _LP_LP_c_double
libwarpx.warpx_getEfield.restype = _LP_LP_c_double
libwarpx.warpx_getEfieldLoVects.restype = _LP_c_int
libwarpx.warpx_getEfieldCP.restype = _LP_LP_c_double
libwarpx.warpx_getEfieldCPLoVects.restype = _LP_c_int
libwarpx.warpx_getEfieldFP.restype = _LP_LP_c_double
libwarpx.warpx_getEfieldFPLoVects.restype = _LP_c_int
libwarpx.warpx_getBfield.restype = _LP_LP_c_double
libwarpx.warpx_getBfieldLoVects.restype = _LP_c_int
libwarpx.warpx_getBfieldCP.restype = _LP_LP_c_double
libwarpx.warpx_getBfieldCPLoVects.restype = _LP_c_int
libwarpx.warpx_getBfieldFP.restype = _LP_LP_c_double
libwarpx.warpx_getBfieldFPLoVects.restype = _LP_c_int
libwarpx.warpx_getCurrentDensity.restype = _LP_LP_c_double
libwarpx.warpx_getCurrentDensityLoVects.restype = _LP_c_int
libwarpx.warpx_getCurrentDensityCP.restype = _LP_LP_c_double
libwarpx.warpx_getCurrentDensityCPLoVects.restype = _LP_c_int
libwarpx.warpx_getCurrentDensityFP.restype = _LP_LP_c_double
libwarpx.warpx_getCurrentDensityFPLoVects.restype = _LP_c_int

#libwarpx.warpx_getPMLSigma.restype = _LP_c_double
#libwarpx.warpx_getPMLSigmaStar.restype = _LP_c_double
#libwarpx.warpx_ComputePMLFactors.argtypes = (ctypes.c_int, ctypes.c_double)

libwarpx.warpx_addNParticles.argtypes = (ctypes.c_int, ctypes.c_int,
                                         _ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                         _ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                         _ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                         _ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                         _ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                         _ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                         ctypes.c_int,
                                         _ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                         ctypes.c_int)

libwarpx.warpx_getProbLo.restype = ctypes.c_double
libwarpx.warpx_getProbHi.restype = ctypes.c_double
libwarpx.warpx_getistep.restype = ctypes.c_int
libwarpx.warpx_gett_new.restype = ctypes.c_double
libwarpx.warpx_getdt.restype = ctypes.c_double
libwarpx.warpx_maxStep.restype = ctypes.c_int
libwarpx.warpx_stopTime.restype = ctypes.c_double
libwarpx.warpx_checkInt.restype = ctypes.c_int
libwarpx.warpx_plotInt.restype = ctypes.c_int
libwarpx.warpx_finestLevel.restype = ctypes.c_int

libwarpx.warpx_EvolveE.argtypes = [ctypes.c_double]
libwarpx.warpx_EvolveB.argtypes = [ctypes.c_double]
libwarpx.warpx_FillBoundaryE.argtypes = []
libwarpx.warpx_FillBoundaryB.argtypes = []
libwarpx.warpx_UpdateAuxilaryData.argtypes = []
libwarpx.warpx_SyncCurrent.argtypes = []
libwarpx.warpx_PushParticlesandDepose.argtypes = [ctypes.c_double]
libwarpx.warpx_getistep.argtypes = [ctypes.c_int]
libwarpx.warpx_setistep.argtypes = [ctypes.c_int, ctypes.c_int]
libwarpx.warpx_gett_new.argtypes = [ctypes.c_int]
libwarpx.warpx_sett_new.argtypes = [ctypes.c_int, ctypes.c_double]
libwarpx.warpx_getdt.argtypes = [ctypes.c_int]

def get_nattr():
    '''

    Get the number of extra attributes.

    '''
    # --- The -3 is because the comps include the velocites
    return libwarpx.warpx_nComps() - 3

def amrex_init(argv):
    # --- Construct the ctype list of strings to pass in
    argc = len(argv)
    argvC = (_LP_c_char * (argc+1))()
    for i, arg in enumerate(argv):
        enc_arg = arg.encode('utf-8')
        argvC[i] = ctypes.create_string_buffer(enc_arg)

    libwarpx.amrex_init(argc, argvC)

def initialize(argv=None):
    '''

    Initialize WarpX and AMReX. Must be called before
    doing anything else.

    '''
    if argv is None:
        argv = sys.argv
    amrex_init(argv)
    libwarpx.warpx_ConvertLabParamsToBoost()
    libwarpx.warpx_init()


def finalize(finalize_mpi=1):
    '''

    Call finalize for WarpX and AMReX. Must be called at
    the end of your script.

    '''
    libwarpx.warpx_finalize()
    libwarpx.amrex_finalize(finalize_mpi)


def evolve(num_steps=-1):
    '''

    Evolve the simulation for num_steps steps. If num_steps=-1,
    the simulation will be run until the end as specified in the
    inputs file.

    Parameters
    ----------

    num_steps: int, the number of steps to take

    '''

    libwarpx.warpx_evolve(num_steps);


#def get_sigma(direction):
#    '''
#
#    Return the 'sigma' PML coefficients for the electric field
#    in a given direction.
#
#    '''
#
#    size = ctypes.c_int(0)
#    data = libwarpx.warpx_getPMLSigma(direction, ctypes.byref(size))
#    arr = np.ctypeslib.as_array(data, (size.value,))
#    arr.setflags(write=1)
#    return arr
#
#
#def get_sigma_star(direction):
#    '''
#
#    Return the 'sigma*' PML coefficients for the magnetic field
#    in the given direction.
#
#    '''
#
#    size = ctypes.c_int(0)
#    data = libwarpx.warpx_getPMLSigmaStar(direction, ctypes.byref(size))
#    arr = np.ctypeslib.as_array(data, (size.value,))
#    arr.setflags(write=1)
#    return arr
#
#
#def compute_pml_factors(lev, dt):
#    '''
#
#    This recomputes the PML coefficients for a given level, using the
#    time step dt. This needs to be called after modifying the coefficients
#    from Python.
#
#    '''
#
#    libwarpx.warpx_ComputePMLFactors(lev, dt)

def add_particles(species_number=0,
                  x=0., y=0., z=0., ux=0., uy=0., uz=0., attr=0.,
                  unique_particles=True):
    '''

    A function for adding particles to the WarpX simulation.

    Parameters
    ----------

    species_number   : the species to add the particle to (default = 0)
    x, y, z          : arrays or scalars of the particle positions (default = 0.)
    ux, uy, uz       : arrays or scalars of the particle momenta (default = 0.)
    attr             : a 2D numpy array or scalar with the particle attributes (default = 0.)
    unique_particles : whether the particles are unique or duplicated on
                       several processes. (default = True)

    '''

    # --- Get length of arrays, set to one for scalars
    lenx = np.size(x)
    leny = np.size(y)
    lenz = np.size(z)
    lenux = np.size(ux)
    lenuy = np.size(uy)
    lenuz = np.size(uz)
    lenattr = np.size(attr)

    if (lenx == 0 or leny == 0 or lenz == 0 or lenux == 0 or
        lenuy == 0 or lenuz == 0 or lenattr == 0):
        return

    maxlen = max(lenx, leny, lenz, lenux, lenuy, lenuz, lenattr)
    assert lenx==maxlen or lenx==1, "Length of x doesn't match len of others"
    assert leny==maxlen or leny==1, "Length of y doesn't match len of others"
    assert lenz==maxlen or lenz==1, "Length of z doesn't match len of others"
    assert lenux==maxlen or lenux==1, "Length of ux doesn't match len of others"
    assert lenuy==maxlen or lenuy==1, "Length of uy doesn't match len of others"
    assert lenuz==maxlen or lenuz==1, "Length of uz doesn't match len of others"
    assert lenattr==maxlen or lenattr==1, "Length of attr doesn't match len of others"

    if lenx == 1:
        x = np.array(x)*np.ones(maxlen)
    if leny == 1:
        y = np.array(y)*np.ones(maxlen)
    if lenz == 1:
        z = np.array(z)*np.ones(maxlen)
    if lenux == 1:
        ux = np.array(ux)*np.ones(maxlen)
    if lenuy == 1:
        uy = np.array(uy)*np.ones(maxlen)
    if lenuz == 1:
        uz = np.array(uz)*np.ones(maxlen,'d')
    if lenattr == 1:
        nattr = get_nattr()
        attr = np.array(attr)*np.ones([maxlen,nattr])

    libwarpx.warpx_addNParticles(species_number, x.size,
                                 x, y, z, ux, uy, uz,
                                 attr.shape[-1], attr, unique_particles)

def get_particle_structs(species_number):
    '''

    This returns a list of numpy arrays containing the particle struct data
    on each tile for this process. The particle data is represented as a structured
    numpy array and contains the particle 'x', 'y', 'z', 'id', and 'cpu'.

    The data for the numpy arrays are not copied, but share the underlying
    memory buffer with WarpX. The numpy arrays are fully writeable.

    Parameters
    ----------

        species_number : the species id that the data will be returned for

    Returns
    -------

        A List of numpy arrays.

    '''

    particles_per_tile = _LP_c_int()
    num_tiles = ctypes.c_int(0)
    data = libwarpx.warpx_getParticleStructs(species_number,
                                             ctypes.byref(num_tiles),
                                             ctypes.byref(particles_per_tile))

    particle_data = []
    for i in range(num_tiles.value):
        arr = _array1d_from_pointer(data[i], _p_dtype, particles_per_tile[i])
        particle_data.append(arr)

    _libc.free(particles_per_tile)
    _libc.free(data)
    return particle_data


def get_particle_arrays(species_number, comp):
    '''

    This returns a list of numpy arrays containing the particle array data
    on each tile for this process.

    The data for the numpy arrays are not copied, but share the underlying
    memory buffer with WarpX. The numpy arrays are fully writeable.

    Parameters
    ----------

        species_number : the species id that the data will be returned for
        comp           : the component of the array data that will be returned.

    Returns
    -------

        A List of numpy arrays.

    '''

    particles_per_tile = _LP_c_int()
    num_tiles = ctypes.c_int(0)
    data = libwarpx.warpx_getParticleArrays(species_number, comp,
                                            ctypes.byref(num_tiles),
                                            ctypes.byref(particles_per_tile))

    particle_data = []
    for i in range(num_tiles.value):
        arr = np.ctypeslib.as_array(data[i], (particles_per_tile[i],))
        arr.setflags(write=1)
        particle_data.append(arr)

    _libc.free(particles_per_tile)
    _libc.free(data)
    return particle_data


def get_particle_x(species_number):
    '''

    Return a list of numpy arrays containing the particle 'x'
    positions on each tile.

    '''
    structs = get_particle_structs(species_number)
    return [struct['x'] for struct in structs]


def get_particle_y(species_number):
    '''

    Return a list of numpy arrays containing the particle 'y'
    positions on each tile.

    '''
    structs = get_particle_structs(species_number)
    return [struct['y'] for struct in structs]


def get_particle_z(species_number):
    '''

    Return a list of numpy arrays containing the particle 'z'
    positions on each tile.

    '''
    structs = get_particle_structs(species_number)
    return [struct['z'] for struct in structs]


def get_particle_id(species_number):
    '''

    Return a list of numpy arrays containing the particle 'id'
    positions on each tile.

    '''
    structs = get_particle_structs(species_number)
    return [struct['id'] for struct in structs]


def get_particle_cpu(species_number):
    '''

    Return a list of numpy arrays containing the particle 'cpu'
    positions on each tile.

    '''
    structs = get_particle_structs(species_number)
    return [struct['cpu'] for struct in structs]


def get_particle_weight(species_number):
    '''

    Return a list of numpy arrays containing the particle
    weight on each tile.

    '''

    return get_particle_arrays(species_number, 0)


def get_particle_ux(species_number):
    '''

    Return a list of numpy arrays containing the particle
    x momentum on each tile.

    '''

    return get_particle_arrays(species_number, 1)


def get_particle_uy(species_number):
    '''

    Return a list of numpy arrays containing the particle
    y momentum on each tile.

    '''

    return get_particle_arrays(species_number, 2)


def get_particle_uz(species_number):
    '''

    Return a list of numpy arrays containing the particle
    z momentum on each tile.

    '''

    return get_particle_arrays(species_number, 3)


def get_particle_Ex(species_number):
    '''

    Return a list of numpy arrays containing the particle
    x electric field on each tile.

    '''

    return get_particle_arrays(species_number, 4)


def get_particle_Ey(species_number):
    '''

    Return a list of numpy arrays containing the particle
    y electric field on each tile.

    '''

    return get_particle_arrays(species_number, 5)


def get_particle_Ez(species_number):
    '''

    Return a list of numpy arrays containing the particle
    z electric field on each tile.

    '''

    return get_particle_arrays(species_number, 6)


def get_particle_Bx(species_number):
    '''

    Return a list of numpy arrays containing the particle
    x magnetic field on each tile.

    '''

    return get_particle_arrays(species_number, 7)


def get_particle_By(species_number):
    '''

    Return a list of numpy arrays containing the particle
    y magnetic field on each tile.

    '''

    return get_particle_arrays(species_number, 8)


def get_particle_Bz(species_number):
    '''

    Return a list of numpy arrays containing the particle
    z magnetic field on each tile.

    '''

    return get_particle_arrays(species_number, 9)


def get_mesh_electric_field(level, direction, include_ghosts=True):
    '''

    This returns a list of numpy arrays containing the mesh electric field
    data on each grid for this process.

    This version is for the full "auxillary" solution on the given level.

    The data for the numpy arrays are not copied, but share the underlying
    memory buffer with WarpX. The numpy arrays are fully writeable.

    Parameters
    ----------

        level          : the AMR level to get the data for
        direction      : the component of the data you want
        include_ghosts : whether to include ghost zones or not

    Returns
    -------

        A List of numpy arrays.

    '''

    assert(level == 0)

    shapes = _LP_c_int()
    size = ctypes.c_int(0)
    ngrow = ctypes.c_int(0)
    data = libwarpx.warpx_getEfield(level, direction,
                                    ctypes.byref(size), ctypes.byref(ngrow),
                                    ctypes.byref(shapes))
    ng = ngrow.value
    grid_data = []
    for i in range(size.value):
        shape = tuple([shapes[dim*i + d] for d in range(dim)])
        # --- The data is stored in Fortran order, hence shape is reversed and a transpose is taken.
        arr = np.ctypeslib.as_array(data[i], shape[::-1]).T
        arr.setflags(write=1)
        if include_ghosts:
            grid_data.append(arr)
        else:
            grid_data.append(arr[[slice(ng, -ng) for _ in range(dim)]])

    _libc.free(shapes)
    _libc.free(data)
    return grid_data


def get_mesh_electric_field_cp(level, direction, include_ghosts=True):
    '''

    This returns a list of numpy arrays containing the mesh electric field
    data on each grid for this process. This version returns the field on
    the coarse patch for the given level.

    The data for the numpy arrays are not copied, but share the underlying
    memory buffer with WarpX. The numpy arrays are fully writeable.

    Parameters
    ----------

        level          : the AMR level to get the data for
        direction      : the component of the data you want
        include_ghosts : whether to include ghost zones or not

    Returns
    -------

        A List of numpy arrays.

    '''

    assert(level == 0)

    shapes = _LP_c_int()
    size = ctypes.c_int(0)
    ngrow = ctypes.c_int(0)
    data = libwarpx.warpx_getEfieldCP(level, direction,
                                      ctypes.byref(size), ctypes.byref(ngrow),
                                      ctypes.byref(shapes))
    ng = ngrow.value
    grid_data = []
    for i in range(size.value):
        shape = tuple([shapes[dim*i + d] for d in range(dim)])
        # --- The data is stored in Fortran order, hence shape is reversed and a transpose is taken.
        arr = np.ctypeslib.as_array(data[i], shape[::-1]).T
        arr.setflags(write=1)
        if include_ghosts:
            grid_data.append(arr)
        else:
            grid_data.append(arr[[slice(ng, -ng) for _ in range(dim)]])

    _libc.free(shapes)
    _libc.free(data)
    return grid_data


def get_mesh_electric_field_fp(level, direction, include_ghosts=True):
    '''

    This returns a list of numpy arrays containing the mesh electric field
    data on each grid for this process. This version returns the field on
    the fine patch for the given level.

    The data for the numpy arrays are not copied, but share the underlying
    memory buffer with WarpX. The numpy arrays are fully writeable.

    Parameters
    ----------

        level          : the AMR level to get the data for
        direction      : the component of the data you want
        include_ghosts : whether to include ghost zones or not

    Returns
    -------

        A List of numpy arrays.

    '''

    assert(level == 0)

    shapes = _LP_c_int()
    size = ctypes.c_int(0)
    ngrow = ctypes.c_int(0)
    data = libwarpx.warpx_getEfieldFP(level, direction,
                                      ctypes.byref(size), ctypes.byref(ngrow),
                                      ctypes.byref(shapes))
    ng = ngrow.value
    grid_data = []
    for i in range(size.value):
        shape = tuple([shapes[dim*i + d] for d in range(dim)])
        # --- The data is stored in Fortran order, hence shape is reversed and a transpose is taken.
        arr = np.ctypeslib.as_array(data[i], shape[::-1]).T
        arr.setflags(write=1)
        if include_ghosts:
            grid_data.append(arr)
        else:
            grid_data.append(arr[[slice(ng, -ng) for _ in range(dim)]])

    _libc.free(shapes)
    _libc.free(data)
    return grid_data


def get_mesh_magnetic_field(level, direction, include_ghosts=True):
    '''

    This returns a list of numpy arrays containing the mesh magnetic field
    data on each grid for this process.

    This version is for the full "auxillary" solution on the given level.

    The data for the numpy arrays are not copied, but share the underlying
    memory buffer with WarpX. The numpy arrays are fully writeable.

    Parameters
    ----------

        level          : the AMR level to get the data for
        direction      : the component of the data you want
        include_ghosts : whether to include ghost zones or not

    Returns
    -------

        A List of numpy arrays.

    '''

    assert(level == 0)

    shapes = _LP_c_int()
    size = ctypes.c_int(0)
    ngrow = ctypes.c_int(0)
    data = libwarpx.warpx_getBfield(level, direction,
                                    ctypes.byref(size), ctypes.byref(ngrow),
                                    ctypes.byref(shapes))
    ng = ngrow.value
    grid_data = []
    for i in range(size.value):
        shape = tuple([shapes[dim*i + d] for d in range(dim)])
        # --- The data is stored in Fortran order, hence shape is reversed and a transpose is taken.
        arr = np.ctypeslib.as_array(data[i], shape[::-1]).T
        arr.setflags(write=1)
        if include_ghosts:
            grid_data.append(arr)
        else:
            grid_data.append(arr[[slice(ng, -ng) for _ in range(dim)]])

    _libc.free(shapes)
    _libc.free(data)
    return grid_data


def get_mesh_magnetic_field_cp(level, direction, include_ghosts=True):
    '''

    This returns a list of numpy arrays containing the mesh magnetic field
    data on each grid for this process. This version returns the field on
    the coarse patch for the given level.

    The data for the numpy arrays are not copied, but share the underlying
    memory buffer with WarpX. The numpy arrays are fully writeable.

    Parameters
    ----------

        level          : the AMR level to get the data for
        direction      : the component of the data you want
        include_ghosts : whether to include ghost zones or not

    Returns
    -------

        A List of numpy arrays.

    '''

    assert(level == 0)

    shapes = _LP_c_int()
    size = ctypes.c_int(0)
    ngrow = ctypes.c_int(0)
    data = libwarpx.warpx_getBfieldCP(level, direction,
                                      ctypes.byref(size), ctypes.byref(ngrow),
                                      ctypes.byref(shapes))
    ng = ngrow.value
    grid_data = []
    for i in range(size.value):
        shape = tuple([shapes[dim*i + d] for d in range(dim)])
        # --- The data is stored in Fortran order, hence shape is reversed and a transpose is taken.
        arr = np.ctypeslib.as_array(data[i], shape[::-1]).T
        arr.setflags(write=1)
        if include_ghosts:
            grid_data.append(arr)
        else:
            grid_data.append(arr[[slice(ng, -ng) for _ in range(dim)]])

    _libc.free(shapes)
    _libc.free(data)
    return grid_data


def get_mesh_magnetic_field_fp(level, direction, include_ghosts=True):
    '''

    This returns a list of numpy arrays containing the mesh magnetic field
    data on each grid for this process. This version returns the field on
    the fine patch for the given level.

    The data for the numpy arrays are not copied, but share the underlying
    memory buffer with WarpX. The numpy arrays are fully writeable.

    Parameters
    ----------

        level          : the AMR level to get the data for
        direction      : the component of the data you want
        include_ghosts : whether to include ghost zones or not

    Returns
    -------

        A List of numpy arrays.

    '''

    assert(level == 0)

    shapes = _LP_c_int()
    size = ctypes.c_int(0)
    ngrow = ctypes.c_int(0)
    data = libwarpx.warpx_getBfieldFP(level, direction,
                                      ctypes.byref(size), ctypes.byref(ngrow),
                                      ctypes.byref(shapes))
    ng = ngrow.value
    grid_data = []
    for i in range(size.value):
        shape = tuple([shapes[dim*i + d] for d in range(dim)])
        # --- The data is stored in Fortran order, hence shape is reversed and a transpose is taken.
        arr = np.ctypeslib.as_array(data[i], shape[::-1]).T
        arr.setflags(write=1)
        if include_ghosts:
            grid_data.append(arr)
        else:
            grid_data.append(arr[[slice(ng, -ng) for _ in range(dim)]])

    _libc.free(shapes)
    _libc.free(data)
    return grid_data


def get_mesh_current_density(level, direction, include_ghosts=True):
    '''

    This returns a list of numpy arrays containing the mesh current density
    data on each grid for this process.

    The data for the numpy arrays are not copied, but share the underlying
    memory buffer with WarpX. The numpy arrays are fully writeable.

    Parameters
    ----------

        level          : the AMR level to get the data for
        direction      : the component of the data you want
        include_ghosts : whether to include ghost zones or not

    Returns
    -------

        A List of numpy arrays.

    '''

    assert(level == 0)

    shapes = _LP_c_int()
    size = ctypes.c_int(0)
    ngrow = ctypes.c_int(0)
    data = libwarpx.warpx_getCurrentDensity(level, direction,
                                            ctypes.byref(size), ctypes.byref(ngrow),
                                            ctypes.byref(shapes))
    ng = ngrow.value
    grid_data = []
    for i in range(size.value):
        shape = tuple([shapes[dim*i + d] for d in range(dim)])
        # --- The data is stored in Fortran order, hence shape is reversed and a transpose is taken.
        arr = np.ctypeslib.as_array(data[i], shape[::-1]).T
        arr.setflags(write=1)
        if include_ghosts:
            grid_data.append(arr)
        else:
            grid_data.append(arr[[slice(ng, -ng) for _ in range(dim)]])

    _libc.free(shapes)
    _libc.free(data)
    return grid_data


def get_mesh_current_density_cp(level, direction, include_ghosts=True):
    '''

    This returns a list of numpy arrays containing the mesh current density
    data on each grid for this process. This version returns the density for
    the coarse patch on the given level.

    The data for the numpy arrays are not copied, but share the underlying
    memory buffer with WarpX. The numpy arrays are fully writeable.

    Parameters
    ----------

        level          : the AMR level to get the data for
        direction      : the component of the data you want
        include_ghosts : whether to include ghost zones or not

    Returns
    -------

        A List of numpy arrays.

    '''

    assert(level == 0)

    shapes = _LP_c_int()
    size = ctypes.c_int(0)
    ngrow = ctypes.c_int(0)
    data = libwarpx.warpx_getCurrentDensityCP(level, direction,
                                              ctypes.byref(size), ctypes.byref(ngrow),
                                              ctypes.byref(shapes))
    ng = ngrow.value
    grid_data = []
    for i in range(size.value):
        shape = tuple([shapes[dim*i + d] for d in range(dim)])
        # --- The data is stored in Fortran order, hence shape is reversed and a transpose is taken.
        arr = np.ctypeslib.as_array(data[i], shape[::-1]).T
        arr.setflags(write=1)
        if include_ghosts:
            grid_data.append(arr)
        else:
            grid_data.append(arr[[slice(ng, -ng) for _ in range(dim)]])

    _libc.free(shapes)
    _libc.free(data)
    return grid_data


def get_mesh_current_density_fp(level, direction, include_ghosts=True):
    '''

    This returns a list of numpy arrays containing the mesh current density
    data on each grid for this process. This version returns the density on
    the fine patch for the given level.

    The data for the numpy arrays are not copied, but share the underlying
    memory buffer with WarpX. The numpy arrays are fully writeable.

    Parameters
    ----------

        level          : the AMR level to get the data for
        direction      : the component of the data you want
        include_ghosts : whether to include ghost zones or not

    Returns
    -------

        A List of numpy arrays.

    '''

    assert(level == 0)

    shapes = _LP_c_int()
    size = ctypes.c_int(0)
    ngrow = ctypes.c_int(0)
    data = libwarpx.warpx_getCurrentDensityFP(level, direction,
                                              ctypes.byref(size), ctypes.byref(ngrow),
                                              ctypes.byref(shapes))
    ng = ngrow.value
    grid_data = []
    for i in range(size.value):
        shape = tuple([shapes[dim*i + d] for d in range(dim)])
        # --- The data is stored in Fortran order, hence shape is reversed and a transpose is taken.
        arr = np.ctypeslib.as_array(data[i], shape[::-1]).T
        arr.setflags(write=1)
        if include_ghosts:
            grid_data.append(arr)
        else:
            grid_data.append(arr[[slice(ng, -ng) for _ in range(dim)]])

    _libc.free(shapes)
    _libc.free(data)
    return grid_data


def _get_mesh_array_lovects(level, direction, include_ghosts=True, getarrayfunc=None):
    assert(0 <= level and level <= libwarpx.warpx_finestLevel())

    size = ctypes.c_int(0)
    ngrow = ctypes.c_int(0)
    data = getarrayfunc(level, direction, ctypes.byref(size), ctypes.byref(ngrow))

    lovects_ref = np.ctypeslib.as_array(data, (size.value, dim))

    # --- Make a copy of the data to avoid memory problems
    # --- Also, take the transpose to give shape (dims, number of grids)
    lovects = lovects_ref.copy().T

    if not include_ghosts:
        lovects += ngrow.value

    del lovects_ref
    _libc.free(data)
    return lovects


def get_mesh_electric_field_lovects(level, direction, include_ghosts=True):
    '''

    This returns a list of the lo vectors of the arrays containing the mesh electric field
    data on each grid for this process.

    This version is for the full "auxillary" solution on the given level.

    Parameters
    ----------

        level          : the AMR level to get the data for
        direction      : the component of the data you want
        include_ghosts : whether to include ghost zones or not

    Returns
    -------

        A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

    '''
    return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getEfieldLoVects)


def get_mesh_electric_field_cp_lovects(level, direction, include_ghosts=True):
    '''

    This returns a list of the lo vectors of the arrays containing the mesh electric field
    data on each grid for this process.

    Parameters
    ----------

        level          : the AMR level to get the data for
        direction      : the component of the data you want
        include_ghosts : whether to include ghost zones or not

    Returns
    -------

        A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

    '''
    return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getEfieldCPLoVects)


def get_mesh_electric_field_fp_lovects(level, direction, include_ghosts=True):
    '''

    This returns a list of the lo vectors of the arrays containing the mesh electric field
    data on each grid for this process.

    Parameters
    ----------

        level          : the AMR level to get the data for
        direction      : the component of the data you want
        include_ghosts : whether to include ghost zones or not

    Returns
    -------

        A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

    '''
    return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getEfieldFPLoVects)


def get_mesh_magnetic_field_lovects(level, direction, include_ghosts=True):
    '''

    This returns a list of the lo vectors of the arrays containing the mesh electric field
    data on each grid for this process.

    This version is for the full "auxillary" solution on the given level.

    Parameters
    ----------

        level          : the AMR level to get the data for
        direction      : the component of the data you want
        include_ghosts : whether to include ghost zones or not

    Returns
    -------

        A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

    '''
    return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getBfieldLoVects)


def get_mesh_magnetic_field_cp_lovects(level, direction, include_ghosts=True):
    '''

    This returns a list of the lo vectors of the arrays containing the mesh electric field
    data on each grid for this process.

    Parameters
    ----------

        level          : the AMR level to get the data for
        direction      : the component of the data you want
        include_ghosts : whether to include ghost zones or not

    Returns
    -------

        A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

    '''
    return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getBfieldCPLoVects)


def get_mesh_magnetic_field_fp_lovects(level, direction, include_ghosts=True):
    '''

    This returns a list of the lo vectors of the arrays containing the mesh electric field
    data on each grid for this process.

    Parameters
    ----------

        level          : the AMR level to get the data for
        direction      : the component of the data you want
        include_ghosts : whether to include ghost zones or not

    Returns
    -------

        A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

    '''
    return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getBfieldFPLoVects)


def get_mesh_current_density_lovects(level, direction, include_ghosts=True):
    '''

    This returns a list of the lo vectors of the arrays containing the mesh electric field
    data on each grid for this process.

    Parameters
    ----------

        level          : the AMR level to get the data for
        direction      : the component of the data you want
        include_ghosts : whether to include ghost zones or not

    Returns
    -------

        A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

    '''
    return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getCurrentDensityLoVects)


def get_mesh_current_density_cp_lovects(level, direction, include_ghosts=True):
    '''

    This returns a list of the lo vectors of the arrays containing the mesh electric field
    data on each grid for this process.

    Parameters
    ----------

        level          : the AMR level to get the data for
        direction      : the component of the data you want
        include_ghosts : whether to include ghost zones or not

    Returns
    -------

        A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

    '''
    return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getCurrentDensityCPLoVects)

def get_mesh_current_density_fp_lovects(level, direction, include_ghosts=True):
    '''

    This returns a list of the lo vectors of the arrays containing the mesh electric field
    data on each grid for this process.

    Parameters
    ----------

        level          : the AMR level to get the data for
        direction      : the component of the data you want
        include_ghosts : whether to include ghost zones or not

    Returns
    -------

        A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

    '''
    return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getCurrentDensityFPLoVects)
