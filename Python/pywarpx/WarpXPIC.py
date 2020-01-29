# Copyright 2017 David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from warp.run_modes.timestepper import PICAPI
from ._libwarpx import libwarpx

class WarpXPIC(PICAPI):

    def get_time(self):
        return libwarpx.warpx_gett_new(0)

    def set_time(self, time):
        for i in range(libwarpx.warpx_finestLevel()+1):
            libwarpx.warpx_sett_new(i, time)

    def get_step_size(self):
        libwarpx.warpx_ComputeDt()
        return libwarpx.warpx_getdt(0)

    def get_step_number(self):
        return libwarpx.warpx_getistep(0)

    def set_step_number(self, it):
        for i in range(libwarpx.warpx_finestLevel()+1):
            libwarpx.warpx_setistep(i, it)

    def push_positions(self, dt):
        libwarpx.warpx_PushX(0, dt)

    def push_velocities_withE(self, dt):
        libwarpx.warpx_EPushV(0, dt)

    def push_velocities_withB(self, dt):
        libwarpx.warpx_BPushV(0, dt)

    def get_self_fields(self):
        libwarpx.warpx_FieldGather(0)

    def calculate_source(self):
        libwarpx.warpx_CurrentDeposition(0)

    def push_Efields(self, dt):
        libwarpx.warpx_EvolveE(0, dt)
        libwarpx.warpx_FillBoundaryE(0, True)

    def push_Bfields(self, dt):
        libwarpx.warpx_EvolveB(0, dt)
        libwarpx.warpx_FillBoundaryB(0, True)

    def apply_particle_boundary_conditions(self):
        libwarpx.mypc_Redistribute() # Redistribute particles
        libwarpx.warpx_MoveWindow() # !!! not the correct place yet

