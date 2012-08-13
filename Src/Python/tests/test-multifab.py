
from pyboxlib import *

def test_multifab():
  pybl.open()

  la = layout()
  la.create(boxes=[ [(0,0), (31,31)] ])

  mfab = multifab()
  mfab.create(la)

  assert mfab.nc == 1
  assert mfab.ng == 0


if __name__ == '__main__':
  test_multifab()
