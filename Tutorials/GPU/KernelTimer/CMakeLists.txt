if (NOT (DIM EQUAL 3))
   return()
endif ()

set(_sources main.cpp )

setup_tutorial(_sources)

unset(_sources)
