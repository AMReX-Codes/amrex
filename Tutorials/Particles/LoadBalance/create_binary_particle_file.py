import yt
import numpy as np

#ds = yt.load("enzo_tiny_cosmology/DD0020/DD0020")
ds = yt.load("Enzo_64/DD0040/data0040")
#ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
ad = ds.all_data()

pos = ad['particle_position']

num_particles = np.array(pos.shape[0], dtype=np.int64)
dim = np.array(pos.shape[1], dtype=np.int32)
nxtra = np.array([0], dtype=np.int32)

particle_positions = pos.reshape(pos.shape[0]*pos.shape[1])

data = particle_positions.d

print(num_particles)

with file("binary_particle_file.dat", "w") as f:
    f.write(num_particles.tobytes())
    f.write(dim.tobytes())
    f.write(nxtra.tobytes())
    f.write(data.tobytes())
