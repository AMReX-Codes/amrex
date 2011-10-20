python module pyfboxlib

    interface
include(fboxlib.pyf)
    end interface

  usercode '''
include(boxlib_numpy.c)
  '''
  pymethoddef '''
undivert(1)
  '''

end python module pyfboxlib
