
from boxlib import fbl

def test_multifab():
  fbl.open()

  la = fbl.layout()
  la.create(boxes=[ [(0,0), (31,31)] ])

  mfab = fbl.multifab()
  mfab.create(la)

  fab = mfab.fab(1)
  fab[0,0] = 22.0

  assert mfab.nc == 1
  assert mfab.ng == 0

  fbl.close()


if __name__ == '__main__':
  test_multifab()
