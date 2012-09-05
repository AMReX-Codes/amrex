vpath %.f   . $(VPATH_LOCATIONS)
vpath %.c   . $(VPATH_LOCATIONS)
vpath %.h   . $(VPATH_LOCATIONS)
vpath %.f90 . $(VPATH_LOCATIONS)

doc:	$(html_sources)
	mv *.html $(hdir)

clean::
	$(RM) ./*.o ./*.mod $(mdir)/*.mod $(odir)/*.o *.$(suf).exe *~
	$(RM) $(odir)/*.il
	$(RM) $(tdir)/f90.depends $(tdir)/c.depends
	$(RM) *.html
	$(RM) TAGS deppairs tags

realclean:: clean
	$(RM) -fr t
	$(RM) *.exe

deppairs: $(f90sources) $(fsources)
	perl $(MODDEP) --tsort $^ > deppairs

TAGS:	$(sources)
	ctags -e --verbose=yes --fortran-kinds=+i $(abspath $^)

tags:	$(sources)
	ctags --verbose=yes --fortran-kinds=+i $^

# should prevent deletion of .o files
.SECONDARY: $(objects)


%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o main.$(suf).exe $< $(objects) $(libraries)
else
	@echo "Linking $@ ..."
	@$(LINK.f90) -o main.$(suf).exe $< $(objects) $(libraries)
endif

${odir}/%.o: %.f
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	@if [ ! -d $(mdir) ]; then mkdir -p $(mdir); fi
ifdef MKVERBOSE
	$(COMPILE.f) $(OUTPUT_OPTION) $<
else
	@echo "Building $< ..."
	@$(COMPILE.f) $(OUTPUT_OPTION) $<
endif

${odir}/%.o: %.f90
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	@if [ ! -d $(mdir) ]; then mkdir -p $(mdir); fi
ifdef MKVERBOSE
	$(COMPILE.f90) $(OUTPUT_OPTION) $<
else
	@echo "Building $< ..."
	@$(COMPILE.f90) $(OUTPUT_OPTION) $<
endif

${odir}/%.o: %.c
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
ifdef MKVERBOSE
	$(COMPILE.c) $(OUTPUT_OPTION) $<
else
	@echo "Building $< ..."
	@$(COMPILE.c) $(OUTPUT_OPTION) $<
endif

${hdir}/%.html: %.f90
	@if [ ! -d $(hdir) ]; then mkdir -p $(hdir); fi
	$(F90DOC) $(F90DOC_OPTION) $<

${hdir}/%.html: %.f
	@if [ ! -d $(hdir) ]; then mkdir -p $(hdir); fi
	$(F90DOC) $(F90DOC_OPTION) $<

$(tdir)/f90.depends: $(fsources) $(f90sources)
	@if [ ! -d $(tdir) ]; then mkdir -p $(tdir); fi
ifdef MKVERBOSE
	perl $(MODDEP) $(f_includes) --odir $(odir)  $^ > $(tdir)/f90.depends 
else
	@echo "Building f90/f dependency File ..."
	@perl $(MODDEP) $(f_includes) --odir $(odir) $^ > $(tdir)/f90.depends 
endif


$(tdir)/c.depends:  $(csources)
	@if [ ! -d $(tdir) ]; then mkdir -p $(tdir); fi
ifdef MKVERBOSE
	perl $(MKDEP) $(c_includes) --odir $(odir) $^ > $(tdir)/c.depends 
else
	@echo "Building c dependency File ..."
	@perl $(MKDEP) $(c_includes) --odir $(odir) $^ > $(tdir)/c.depends 
endif

ifneq ($(MAKECMDGOALS),realclean)
ifneq ($(MAKECMDGOALS),clean)
include $(tdir)/f90.depends

ifdef csources
include $(tdir)/c.depends
endif
endif
endif




