-include $(BOXLIB_HOME)/Tools/F_mk/Make.local

vpath %.f   . $(VPATH_LOCATIONS)
vpath %.c   . $(VPATH_LOCATIONS)
vpath %.cpp . $(VPATH_LOCATIONS)
vpath %.h   . $(VPATH_LOCATIONS)
vpath %.f90 . $(VPATH_LOCATIONS)
vpath %.F90 . $(VPATH_LOCATIONS)


ifeq ($(pnames),)

%.$(suf).exe: $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $(objects) $(libraries)
else
	@echo "Linking $@ ..."
	@$(LINK.f90) -o $@ $(objects) $(libraries)
endif

else

%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif

endif


doc:	$(html_sources)
	mv *.html $(hdir)

clean::
	$(RM) ./*.o ./*.mod $(mdir)/*.mod $(odir)/*.o *.$(suf).exe *~
	$(RM) ./*.optrpt $(odir)/*.optrpt
	$(RM) $(odir)/*.il
	$(RM) $(tdir)/f90.depends $(tdir)/c.depends
	$(RM) *.html
	$(RM) TAGS tags

realclean:: clean
	$(RM) -fr t
	$(RM) *.exe

file_locations:
	$(BOXLIB_HOME)/Tools/F_scripts/find_files_vpath.py --vpath "$(VPATH_LOCATIONS)" --files "$(sources)"

TAGS:	$(sources)
	ctags -e --verbose=yes --fortran-kinds=+i $(abspath $^)

tags:	$(sources)
	ctags --verbose=yes --fortran-kinds=+i $^

# should prevent deletion of .o files
.SECONDARY: $(objects)

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

# here we rely on the compiler convention that .F90 files will
# automatically be preprocessed
${odir}/%.o: %.F90
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	@if [ ! -d $(mdir) ]; then mkdir -p $(mdir); fi
ifdef MKVERBOSE
	$(COMPILE.f90) $(FPP_DEFINES) $(OUTPUT_OPTION) $<
else
	@echo "Building $< ..."
	@$(COMPILE.f90) $(FPP_DEFINES) $(OUTPUT_OPTION) $<
endif

${odir}/%.o: %.c
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
ifdef MKVERBOSE
	$(COMPILE.c) $(OUTPUT_OPTION) $<
else
	@echo "Building $< ..."
	@$(COMPILE.c) $(OUTPUT_OPTION) $<
endif

${odir}/%.o: %.cpp
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
ifdef MKVERBOSE
	$(COMPILE.cc) $(OUTPUT_OPTION) $<
else
	@echo "Building $< ..."
	@$(COMPILE.cc) $(OUTPUT_OPTION) $<
endif

${hdir}/%.html: %.f90
	@if [ ! -d $(hdir) ]; then mkdir -p $(hdir); fi
	$(F90DOC) $(F90DOC_OPTION) $<

${hdir}/%.html: %.F90
	@if [ ! -d $(hdir) ]; then mkdir -p $(hdir); fi
	$(F90DOC) $(F90DOC_OPTION) $<

${hdir}/%.html: %.f
	@if [ ! -d $(hdir) ]; then mkdir -p $(hdir); fi
	$(F90DOC) $(F90DOC_OPTION) $<

# determine the build dependencies among the source files.  At the moment,
# we do not preprocess the Fortran files before doing the dependency check.
# If this is an issue, then it is easy to change it, following what is done
# in C_mk/Make.rules
$(tdir)/f90.depends: $(fsources) $(f90sources) $(F90sources)
	@if [ ! -d $(tdir) ]; then mkdir -p $(tdir); fi
	@echo "Building f90/F90/f dependency File ..."
	$(MODDEP) --prefix $(odir) \
            --temp_dir $(tdir) \
            --cpp "cpp -E -traditional" \
            --defines "$(FPPFLAGS) $(FPP_DEFINES)" $^ > $(tdir)/f90.depends 
	@if [ $$? -ne 0 ]; then exit "make fail"; fi

$(tdir)/c.depends:  $(csources) $(cxxsources)
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
else
ifdef cxxsources
include $(tdir)/c.depends
endif
endif

endif
endif

#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make echo-Fmlocs".  This will
# print out the value.
echo-%: ; @echo $* is $($*)

