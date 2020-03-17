# Copyright 2019 David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import warp
from . import fields
from pywarpx import PGroup

# The particle weight is always the first pid
warp.top.wpid = 1

def warp_species(warp_type, picmi_species, level=0):
    """Returns a Warp species that has a reference to the WarpX particles.
    """
    pgroups = PGroup.PGroups(ispecie=picmi_species.species_number, level=level)
    return warp.Species(type=warp_type, pgroups=pgroups)


class _WarpX_FIELDtype(object):
    """Mirrors part of the EM3D_FIELDtype type from Warp
    yf
    """
    def __init__(self, yf):
        self.yf = yf


class _WarpX_BLOCKtype(object):
    """Mirrors part of the EM3D_BLOCKtype type from Warp
    core
    xmin, ymin, zmin
    xmax, ymax, zmax
    dx, dy, dz
    xrbnd, yrbnd, zrbnd
    """
    def __init__(self, fields, picmi_grid):
        self.core = _WarpX_FIELDtype(fields)
        self.picmi_grid = picmi_grid

        self.xmin = self.picmi_grid.lower_bound[0]
        if self.picmi_grid.number_of_dimensions == 3:
            self.ymin = self.picmi_grid.lower_bound[1]
        else:
            self.ymin = 0.
        self.zmin = self.picmi_grid.lower_bound[-1]

        self.xmax = self.picmi_grid.upper_bound[0]
        if self.picmi_grid.number_of_dimensions == 3:
            self.ymax = self.picmi_grid.upper_bound[1]
        else:
            self.ymax = 0.
        self.zmax = self.picmi_grid.upper_bound[-1]

        self.dx = (self.xmax - self.xmin)/self.picmi_grid.number_of_cells[0]
        if self.picmi_grid.number_of_dimensions == 3:
            self.dy = (self.ymax - self.ymin)/self.picmi_grid.number_of_cells[1]
        else:
            self.dy = 1.
        self.dz = (self.zmax - self.zmin)/self.picmi_grid.number_of_cells[-1]

        self.xrbnd = 0
        self.yrbnd = 0
        self.zrbnd = 0


class _WarpX_YEEFIELDtype(object):
    """Mirrors part of the EM3D_YEEFIELDtype type from Warp
    Exp, Eyp, Ezp
    Bxp, Byp, Bzp
    Ex, Ey, Ez
    Bx, By, Bz
    Jx, Jy, Jz
    """
    def __init__(self, level=0):
        self.level = level
        self._Ex_wrap = fields.ExWrapper(level, include_ghosts=True)
        self._Ey_wrap = fields.EyWrapper(level, include_ghosts=True)
        self._Ez_wrap = fields.EzWrapper(level, include_ghosts=True)
        self._Bx_wrap = fields.BxWrapper(level, include_ghosts=True)
        self._By_wrap = fields.ByWrapper(level, include_ghosts=True)
        self._Bz_wrap = fields.BzWrapper(level, include_ghosts=True)
        self._Jx_wrap = fields.JxWrapper(level, include_ghosts=True)
        self._Jy_wrap = fields.JyWrapper(level, include_ghosts=True)
        self._Jz_wrap = fields.JzWrapper(level, include_ghosts=True)

        self._Ex_wrap._getlovects()  # --- Calculated nghosts
        self.nxguard = self._Ex_wrap.nghosts
        self.nyguard = self._Ex_wrap.nghosts
        self.nzguard = self._Ex_wrap.nghosts

    def _get_wrapped_array(self, wrapper):
        result = wrapper[...]
        if len(result.shape) == 2:
            # --- Add the middle dimension that Warp is expecting
            result.shape = [result.shape[0], 1, result.shape[1]]
        return result

    @property
    def Exp(self):
        return self._get_wrapped_array(self._Ex_wrap)
    @property
    def Eyp(self):
        return self._get_wrapped_array(self._Ey_wrap)
    @property
    def Ezp(self):
        return self._get_wrapped_array(self._Ez_wrap)
    @property
    def Bxp(self):
        return self._get_wrapped_array(self._Bx_wrap)
    @property
    def Byp(self):
        return self._get_wrapped_array(self._By_wrap)
    @property
    def Bzp(self):
        return self._get_wrapped_array(self._Bz_wrap)
    @property
    def Ex(self):
        return self._get_wrapped_array(self._Ex_wrap)
    @property
    def Ey(self):
        return self._get_wrapped_array(self._Ey_wrap)
    @property
    def Ez(self):
        return self._get_wrapped_array(self._Ez_wrap)
    @property
    def Bx(self):
        return self._get_wrapped_array(self._Bx_wrap)
    @property
    def By(self):
        return self._get_wrapped_array(self._By_wrap)
    @property
    def Bz(self):
        return self._get_wrapped_array(self._Bz_wrap)
    @property
    def Jx(self):
        return self._get_wrapped_array(self._Jx_wrap)
    @property
    def Jy(self):
        return self._get_wrapped_array(self._Jy_wrap)
    @property
    def Jz(self):
        return self._get_wrapped_array(self._Jz_wrap)


class WarpX_EM3D(warp.EM3D):
    """Mirrors part of the Warp EM3D class, mostly diagnostics.
    """
    def __init__(self, picmi_grid, level=0):
        self.picmi_grid = picmi_grid
        self.level = level

        # --- Only define what is necessary for the diagnostics
        self.fields = _WarpX_YEEFIELDtype(level)
        self.block = _WarpX_BLOCKtype(self.fields, picmi_grid)

        self.isactive = True

        self.l_1dz = (picmi_grid.number_of_dimensions == 1)
        self.l_2dxz = (picmi_grid.number_of_dimensions == 2)
        try:
            picmi_grid.nr
        except AttributeError:
            self.l_2drz = False
        else:
            self.l_2drz = True

        self.l4symtry = False
        self.l2symtry = False

        self.nx = picmi_grid.number_of_cells[0]
        if not self.l_2dxz:
            self.ny = picmi_grid.number_of_cells[1]
        else:
            self.ny = 0
        self.nz = picmi_grid.number_of_cells[-1]

        self.zgrid = 0.  # --- This should be obtained from WarpX

