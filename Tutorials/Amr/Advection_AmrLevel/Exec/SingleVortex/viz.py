import yt
from yt.frontends.boxlib.data_structures import AMReXDataset


ds = AMReXDataset("plt00120")
print(ds.field_list)

sl = yt.SlicePlot(ds, 2, 'phi')
sl.annotate_particles(ds.domain_width[2])
sl.show()
