import yt
import numpy as np

ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")
ad = ds.all_data()

pos = ad['particle_position']

num_particles = np.array(pos.shape[0], dtype=np.int64)
dim = np.array(pos.shape[1], dtype=np.int32)
nxtra = np.array([0], dtype=np.int32)

particle_positions = pos.reshape(pos.shape[0]*pos.shape[1])

data = particle_positions.d

with file("binary_particle_file.dat", "w") as f:
    f.write(num_particles.tobytes())
    f.write(dim.tobytes())
    f.write(nxtra.tobytes())
    f.write(data.tobytes())
