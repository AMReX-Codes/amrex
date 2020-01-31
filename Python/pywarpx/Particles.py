# Copyright 2017-2020 Andrew Myers, David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from .Bucket import Bucket

particles = Bucket('particles', nspecies=0, species_names=[])
particles_list = []

electrons = Bucket('electrons')
electrons.charge = "-q_e"
electrons.mass = "m_e"
electrons.injection_style = "python"

positrons = Bucket('positrons')
positrons.charge = "q_e"
positrons.mass = "m_e"
positrons.injection_style = "python"

protons = Bucket('protons')
protons.charge = "q_e"
protons.mass = "m_p"
protons.injection_style = "python"

particle_dict = {'electrons':electrons,
                 'positrons':positrons,
                 'protons':protons
                 }

def newspecies(name):
    result = Bucket(name)
    particles_list.append(result)
    return result
