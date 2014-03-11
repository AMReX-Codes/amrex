$(PYFBOXLIB): $(objects)
	$(F90) $(F90FLAGS) -shared -o $@ $^
