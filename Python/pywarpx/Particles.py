from .Bucket import Bucket

particles = Bucket('particles', nspecies=0, species_names=None)
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
