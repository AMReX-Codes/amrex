import numpy as np

class AMReXParticleHeader(object):

    def __init__(self, header_filename):

        self.real_component_names = []
        self.int_component_names = []
        with open(header_filename, "r") as f:
            self.version_string = f.readline().strip()

            particle_real_type = self.version_string.split('_')[-1]
            particle_real_type = self.version_string.split('_')[-1]
            if particle_real_type == 'double':
                self.real_type = np.float64
            elif particle_real_type == 'single':
                self.real_type = np.float32
            else:
                raise RuntimeError("Did not recognize particle real type.")
            self.int_type = np.int32

            self.dim = int(f.readline().strip())
            self.num_int_base = 2
            self.num_real_base = self.dim
            self.num_real_extra = int(f.readline().strip())
            for i in range(self.num_real_extra):
                self.real_component_names.append(f.readline().strip())
            self.num_int_extra = int(f.readline().strip())
            for i in range(self.num_int_extra):
                self.int_component_names.append(f.readline().strip())
            self.num_int = self.num_int_base + self.num_int_extra
            self.num_real = self.num_real_base + self.num_real_extra
            self.is_checkpoint = bool(int(f.readline().strip()))
            self.num_particles = int(f.readline().strip())
            self.max_next_id = int(f.readline().strip())
            self.finest_level = int(f.readline().strip())
            self.num_levels = self.finest_level + 1

            if not self.is_checkpoint:
                self.num_int_base = 0
                self.num_int_extra = 0
                self.num_int = 0

            self.grids_per_level = np.zeros(self.num_levels, dtype='int64')
            for level_num in range(self.num_levels):
                self.grids_per_level[level_num] = int(f.readline().strip())

            self.grids = [[]*level_num]
            for level_num in range(self.num_levels):
                for grid_num in range(self.grids_per_level[level_num]):
                    entry = [int(val) for val in f.readline().strip().split()]
                    self.grids[level_num].append(tuple(entry))

                    
def read_particle_data(fn, ptype="particle0"):
    base_fn = fn + "/" + ptype
    header = AMReXParticleHeader(base_fn + "/Header")
    
    idtype = "(%d,)i4" % header.num_int    
    if header.real_type == np.float64:
        fdtype = "(%d,)f8" % header.num_real
    elif header.real_type == np.float32:
        fdtype = "(%d,)f4" % header.num_real
    
    idata = np.empty((header.num_particles, header.num_int ))
    rdata = np.empty((header.num_particles, header.num_real))
    
    ip = 0
    for lvl, level_grids in enumerate(header.grids):
        for (which, count, where) in level_grids:
            if count == 0: continue
            fn = base_fn + "/Level_%d/DATA_%04d" % (lvl, which)

            with open(fn, 'rb') as f:
                f.seek(where)
                ints   = np.fromfile(f, dtype = idtype, count=count)
                floats = np.fromfile(f, dtype = fdtype, count=count)                    

            idata[ip] = ints
            rdata[ip] = floats            
            ip += 1
            
    return idata, rdata


if __name__ == "__main__":
    import pylab as plt
    import glob

    x = []
    y = []
    
    fn_list = glob.glob("plt?????")
    fn_list.sort()
    
    for fn in fn_list:
        idata, rdata = read_particle_data(fn)
        x.append(rdata[0][0])
        y.append(rdata[0][1])

    fig = plt.gcf()
    fig.set_size_inches(8, 8)
    plt.plot(x, y, '.')
    plt.axis((-30e-6, 30e-6, -30e-6, 30e-6))
    ax = plt.gca()
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    plt.savefig('particles.png')
