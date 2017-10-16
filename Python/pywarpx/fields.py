"""Provides wrappers around field and current density on multiFABs

Available routines:

ExWrapper, EyWrapper, EzWrapper
BxWrapper, ByWrapper, BzWrapper
JxWrapper, JyWrapper, JzWrapper

"""
import numpy as np
from . import _libwarpx


class MultiFABWrapper(object):
    """Wrapper around field arrays at level 0
    This provides a convenient way to query and set fields that are broken up into FABs.
    The indexing is based on global indices.
     - direction: component to access, one of the values (0, 1, 2)
     - overlaps: is one along the axes where the grid boundaries overlap the neighboring grid
     - get_lovects: routine that returns the list of lo vectors
     - get_fabs: routine that returns the list of FABs
     - level=0: ignored
    """
    def __init__(self, direction, overlaps, get_lovects, get_fabs, level=0):
        self.direction = direction
        self.overlaps = np.array(overlaps)
        self.get_lovects = get_lovects
        self.get_fabs = get_fabs
        self.level = 0 #level
        self.include_ghosts = False

    def _getlovects(self):
        return self.get_lovects(self.level, self.direction, self.include_ghosts)

    def _gethivects(self):
        lovects = self._getlovects()
        fields = self._getfields()

        hivects = np.zeros_like(lovects)
        for i in range(len(fields)):
            hivects[:,i] = lovects[:,i] + np.array(fields[i].shape) - self.overlaps

        return hivects

    def _getfields(self):
        return self.get_fabs(self.level, self.direction, self.include_ghosts)

    def __len__(self):
        return lend(self._getlovects())

    def __getitem__(self, index):
        """Returns slices of a decomposed array, The shape of
      the object returned depends on the number of ix, iy and iz specified, which
      can be from none to all three. Note that the values of ix, iy and iz are
      relative to the fortran indexing, meaning that 0 is the lower boundary
      of the domain.
        """
        if index == Ellipsis:
            index = (slice(None), slice(None), slice(None))

        ix = index[0]
        iy = index[1]
        iz = index[2]

        lovects = self._getlovects()
        hivects = self._gethivects()
        fields = self._getfields()

        nx = hivects[0,:].max()
        ny = hivects[1,:].max()
        nz = hivects[2,:].max()

        if isinstance(ix, slice):
            ixstart = max(ix.start, 0)
            ixstop = min(ix.stop or nx + 1, nx + self.overlaps[0])
        else:
            ixstart = ix
            ixstop = ix + 1
        if isinstance(iy, slice):
            iystart = max(iy.start, 0)
            iystop = min(iy.stop or ny + 1, ny + self.overlaps[1])
        else:
            iystart = iy
            iystop = iy + 1
        if isinstance(iz, slice):
            izstart = max(iz.start, 0)
            izstop = min(iz.stop or nz + 1, nz + self.overlaps[2])
        else:
            izstart = iz
            izstop = iz + 1

        # --- Setup the size of the array to be returned and create it.
        sss = (max(0, ixstop - ixstart),
               max(0, iystop - iystart),
               max(0, izstop - izstart))
        resultglobal = np.zeros(sss)

        for i in range(len(fields)):

            # --- The ix1, 2 etc are relative to global indexing
            ix1 = max(ixstart, lovects[0,i])
            ix2 = min((ixstop or nx+1), lovects[0,i] + fields[i].shape[0])
            iy1 = max(iystart, lovects[1,i])
            iy2 = min((iystop or ny+1), lovects[1,i] + fields[i].shape[1])
            iz1 = max(izstart, lovects[2,i])
            iz2 = min((izstop or nz+1), lovects[2,i] + fields[i].shape[2])

            if ix1 < ix2 and iy1 < iy2 and iz1 < iz2:

                sss = (slice(ix1 - lovects[0,i], ix2 - lovects[0,i]),
                       slice(iy1 - lovects[1,i], iy2 - lovects[1,i]),
                       slice(iz1 - lovects[2,i], iz2 - lovects[2,i]))

                vslice = (slice(ix1 - ixstart, ix2 - ixstart),
                          slice(iy1 - iystart, iy2 - iystart),
                          slice(iz1 - izstart, iz2 - izstart))

                resultglobal[vslice] = fields[i][sss]

        # --- Now remove any of the reduced dimensions.
        sss = [slice(None), slice(None), slice(None)]
        if not isinstance(ix, slice):
            sss[0] = 0
        if not isinstance(iy, slice):
            sss[1] = 0
        if not isinstance(iz, slice):
            sss[2] = 0

        return resultglobal[sss]

    def __setitem__(self, index, value):
        """Sets slices of a decomposed array. The shape of
      the input object depends on the number of arguments specified, which can
      be from none to all three.
        - value: input array (must be supplied)
        """
        if index == Ellipsis:
            index = (slice(None), slice(None), slice(None))

        ix = index[0]
        iy = index[1]
        iz = index[2]

        lovects = self._getlovects()
        hivects = self._gethivects()
        fields = self._getfields()

        nx = hivects[0,:].max()
        ny = hivects[1,:].max()
        nz = hivects[2,:].max()

        # --- Add extra dimensions so that the input has the same number of
        # --- dimensions as array.
        if isinstance(value, np.ndarray):
            value3d = np.array(value, copy=False)
            sss = list(value3d.shape)
            if not isinstance(ix, slice): sss[0:0] = [1]
            if not isinstance(iy, slice): sss[1:1] = [1]
            if not isinstance(iz, slice): sss[2:2] = [1]
            value3d.shape = sss

        if isinstance(ix, slice):
            ixstart = max(ix.start, 0)
            ixstop = min(ix.stop or nx + 1, nx + self.overlaps[0])
        else:
            ixstart = ix
            ixstop = ix + 1
        if isinstance(iy, slice):
            iystart = max(iy.start, 0)
            iystop = min(iy.stop or ny + 1, ny + self.overlaps[1])
        else:
            iystart = iy
            iystop = iy + 1
        if isinstance(iz, slice):
            izstart = max(iz.start, 0)
            izstop = min(iz.stop or nz + 1, nz + self.overlaps[2])
        else:
            izstart = iz
            izstop = iz + 1

        for i in range(len(fields)):

            # --- The ix1, 2 etc are relative to global indexing
            ix1 = max(ixstart, lovects[0,i])
            ix2 = min((ixstop or nx+1), lovects[0,i] + fields[i].shape[0])
            iy1 = max(iystart, lovects[1,i])
            iy2 = min((iystop or ny+1), lovects[1,i] + fields[i].shape[1])
            iz1 = max(izstart, lovects[2,i])
            iz2 = min((izstop or nz+1), lovects[2,i] + fields[i].shape[2])

            if ix1 < ix2 and iy1 < iy2 and iz1 < iz2:

                sss = (slice(ix1 - lovects[0,i], ix2 - lovects[0,i]),
                       slice(iy1 - lovects[1,i], iy2 - lovects[1,i]),
                       slice(iz1 - lovects[2,i], iz2 - lovects[2,i]))

                if isinstance(value, np.ndarray):
                    vslice = (slice(ix1 - ixstart, ix2 - ixstart),
                              slice(iy1 - iystart, iy2 - iystart),
                              slice(iz1 - izstart, iz2 - izstart))
                    fields[i][sss] = value3d[vslice]
                else:
                    fields[i][sss] = value


def ExWrapper(level=0):
    return MultiFABWrapper(direction=0, overlaps=[0,1,1],
                           get_lovects=_libwarpx.get_mesh_electric_field_lovects,
                           get_fabs=_libwarpx.get_mesh_electric_field, level=level)

def EyWrapper(level=0):
    return MultiFABWrapper(direction=1, overlaps=[1,0,1],
                           get_lovects=_libwarpx.get_mesh_electric_field_lovects,
                           get_fabs=_libwarpx.get_mesh_electric_field, level=level)

def EzWrapper(level=0):
    return MultiFABWrapper(direction=2, overlaps=[1,1,0],
                           get_lovects=_libwarpx.get_mesh_electric_field_lovects,
                           get_fabs=_libwarpx.get_mesh_electric_field, level=level)

def BxWrapper(level=0):
    return MultiFABWrapper(direction=0, overlaps=[1,0,0],
                           get_lovects=_libwarpx.get_mesh_magnetic_field_lovects,
                           get_fabs=_libwarpx.get_mesh_magnetic_field, level=level)

def ByWrapper(level=0):
    return MultiFABWrapper(direction=1, overlaps=[0,1,0],
                           get_lovects=_libwarpx.get_mesh_magnetic_field_lovects,
                           get_fabs=_libwarpx.get_mesh_magnetic_field, level=level)

def BzWrapper(level=0):
    return MultiFABWrapper(direction=2, overlaps=[0,0,1],
                           get_lovects=_libwarpx.get_mesh_magnetic_field_lovects,
                           get_fabs=_libwarpx.get_mesh_magnetic_field, level=level)

def JxWrapper(level=0):
    return MultiFABWrapper(direction=0, overlaps=[0,1,1],
                           get_lovects=_libwarpx.get_mesh_current_density_lovects,
                           get_fabs=_libwarpx.get_mesh_current_density, level=level)

def JyWrapper(level=0):
    return MultiFABWrapper(direction=1, overlaps=[1,0,1],
                           get_lovects=_libwarpx.get_mesh_current_density_lovects,
                           get_fabs=_libwarpx.get_mesh_current_density, level=level)

def JzWrapper(level=0):
    return MultiFABWrapper(direction=2, overlaps=[1,1,0],
                           get_lovects=_libwarpx.get_mesh_current_density_lovects,
                           get_fabs=_libwarpx.get_mesh_current_density, level=level)

