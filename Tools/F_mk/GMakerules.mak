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
	$(F90DOC) $(OUTPUT_OPTION) $<

${hdir}/%.html: %.f
	@if [ ! -d $(hdir) ]; then mkdir -p $(hdir); fi
	$(F90DOC) $(OUTPUT_OPTION) $<

$(tdir)/f90.depends: $(sources) 
	@if [ ! -d $(tdir) ]; then mkdir -p $(tdir); fi
	@echo "Building dependency File ..."
	@perl $(MODDEP) --odir $(odir) $^ > $(tdir)/f90.depends 

-include $(tdir)/f90.depends
