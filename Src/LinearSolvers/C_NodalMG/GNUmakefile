
# To set option variables, go to GNUmakefile.main

all: libproj proj

libproj:
	$(MAKE) -f GNUmakefile.main LBASE=proj

libell:
	$(MAKE) -f GNUmakefile.main LBASE=ell

driver:
	$(MAKE) -f GNUmakefile.main EBASE=driver

proj:
	$(MAKE) -f GNUmakefile.main EBASE=proj

poisson:
	$(MAKE) -f GNUmakefile.main EBASE=poisson
