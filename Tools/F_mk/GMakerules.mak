vpath %.f   . $(VPATH_LOCATIONS)
vpath %.c   . $(VPATH_LOCATIONS)
vpath %.f90 . $(VPATH_LOCATIONS)

doc:	$(html_sources)

clean:
	$(RM) ./*.o ./*.mod $(mdir)/*.mod $(odir)/*.o *.exe *~
	$(RM) $(odir)/*.il
	$(RM) $(tdir)/f90.depends
	$(RM) *.html
	$(RM) TAGS

realclean: clean
	$(RM) -fr t

TAGS:	$(sources)
	etags $^

${odir}/%.o: %.f
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	@if [ ! -d $(mdir) ]; then mkdir -p $(mdir); fi
	$(COMPILE.f) $(OUTPUT_OPTION) $<

${odir}/%.o: %.f90
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	@if [ ! -d $(mdir) ]; then mkdir -p $(mdir); fi
	$(COMPILE.f90) $(OUTPUT_OPTION) $<

${odir}/%.o: %.c
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	$(COMPILE.c) $(OUTPUT_OPTION) $<

${hdir}/%.html: %.f90
	@if [ ! -d $(hdir) ]; then mkdir -p $(hdir); fi
	$(F90DOC) $(OUTPUT_OPTION) $<

${hdir}/%.html: %.f
	@if [ ! -d $(hdir) ]; then mkdir -p $(hdir); fi
	$(F90DOC) $(OUTPUT_OPTION) $<

$(tdir)/f90.depends: $(sources) 
	@if [ ! -d $(tdir) ]; then mkdir -p $(tdir); fi
	perl $(MODDEP) --odir $(odir) $^ > $(tdir)/f90.depends 
-include $(tdir)/f90.depends
