set(_sources ElectrostaticParticleContainer.cpp ElectrostaticParticleContainer.H
   electrostatic_pic_${DIM}d.f90 electrostatic_pic_F.H main.cpp PhysConst.H)

set(_input_files inputs)

setup_tutorial(_sources _input_files)

unset(_sources)
unset(_input_files)
