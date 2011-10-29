pyfboxlib.so: $(objects) $(tdir)/blobjects.f90 $(tdir)/pyfboxlib.pyf
	@echo Running f2py...
	@f2py --quiet --fcompiler=gnu95 --f90exec=$(FC) --f90flags="-I $(mdir)" \
		-c $(tdir)/pyfboxlib.pyf $(tdir)/blobjects.f90 $(PYBOXLIB)/src/boxlib_numpy.f90 $(pybl_sources) $(objects)

$(tdir)/blobjects.f90: $(PYBOXLIB)/src/blobjects.py
	cd $(tdir) && python $(PYBOXLIB)/src/blobjects.py

$(tdir)/%.pyf: %.f90
	f2py --overwrite-signature -h $@ $<

$(tdir)/pyfboxlib.pyf: $(PYBOXLIB)/src/boxlib_numpy.c $(pybl_pyfs) $(pybl_sources)
	@$(PYBOXLIB)/mkpyfboxlib $(PYBOXLIB)/src/boxlib_numpy.c $(pybl_pyfs) > $(tdir)/pyfboxlib.pyf
