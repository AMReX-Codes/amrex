# Copyright 2019-2020 Axel Huebl, Luca Fedeli, Maxence Thevenet
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

'''
This python script reads ionization tables in atomic_data.txt (generated from
the NIST website) and extracts ionization levels into C++ file
IonizationEnergiesTable.H, which contains tables + metadata.
'''

import re, os
import numpy as np

filename = os.path.join( '.', 'atomic_data.txt' )
with open(filename) as f:
    text_data = f.read()

# Read full table from file and get names, atomic numbers and offsets
# position in table of ionization energies for all species
regex_command = '\n\s+(\d+)\s+\|\s+([A-Z]+[a-z]*)\s+\w+\s+\|\s+\+*(\d+)\s+\|\s+\(*\[*(\d+\.*\d*)'
list_of_tuples = re.findall( regex_command, text_data )
ion_atom_numbers = [int(i) for i in list(dict.fromkeys( [x[0] for x in list_of_tuples] ))]
ion_names = list(dict.fromkeys( [x[1] for x in list_of_tuples] ))
ion_offsets = np.concatenate(([0], np.cumsum(np.array(ion_atom_numbers)[:-1])), axis=0)

# Head of CPP file
cpp_string = '// This script was automatically generated!\n'
cpp_string += '// Edit dev/Source/Utils/write_atomic_data_cpp.py instead!\n'
cpp_string += '#ifndef WARPX_IONIZATION_TABLE_H_\n'
cpp_string += '#define WARPX_IONIZATION_TABLE_H_\n\n'
cpp_string += '#include <AMReX_AmrCore.H>\n'
cpp_string += '#include <map>\n\n'

# Map each element to ID in table
cpp_string += 'std::map<std::string, int> ion_map_ids = {'
for count, name in enumerate(ion_names):
    cpp_string += '\n    {"' + name + '", ' + str(count) + '},'
cpp_string = cpp_string[:-1]
cpp_string += ' };\n\n'

# Atomic number of each species
cpp_string += 'const int nelements = ' + str(len(ion_names)) + ';\n\n'
cpp_string += 'const int ion_atomic_numbers[nelements] = {\n    '
for count, atom_num in enumerate(ion_atom_numbers):
    if count%10==0 and count>0: cpp_string = cpp_string[:-2] + ',\n    '
    cpp_string += str(atom_num) + ', '
cpp_string = cpp_string[:-2]
cpp_string += '};\n\n'

# Offset of each element in table of ionization energies
cpp_string += 'const int ion_energy_offsets[nelements] = {\n    '
for count, offset in enumerate(ion_offsets):
    if count%10==0 and count>0: cpp_string = cpp_string[:-2] + ',\n    '
    cpp_string += str(offset) + ', '
cpp_string = cpp_string[:-2]
cpp_string += '};\n\n'

# Table of ionization energies
cpp_string += 'const int energies_tab_length = ' + str(len(list_of_tuples)) + ';\n\n'
cpp_string += 'const amrex::Real table_ionization_energies[energies_tab_length]{'
for element in ion_names:
    cpp_string += '\n    // ' + element + '\n    '
    regex_command = \
        '\n\s+(\d+)\s+\|\s+%s\s+\w+\s+\|\s+\+*(\d+)\s+\|\s+\(*\[*(\d+\.*\d*)' \
        %element
    list_of_tuples = re.findall( regex_command, text_data )
    for count, energy in enumerate([x[2] for x in list_of_tuples]):
        if count%7==0 and count>0: cpp_string = cpp_string[:-2] + ',\n    '
        cpp_string += energy + ', '
    cpp_string = cpp_string[:-1]
cpp_string = cpp_string[:-1]
cpp_string += '\n};\n\n'

# Write the string to file
cpp_string += '#endif // #ifndef WARPX_IONIZATION_TABLE_H_\n'
f= open("IonizationEnergiesTable.H","w")
f.write(cpp_string)
f.close()
