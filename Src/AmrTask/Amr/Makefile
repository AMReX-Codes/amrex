include ../arch.common 

OBJECTS= AMReX_AmrTask.o 

AMRLIB= AMRTask.a

all: $(AMRLIB)

$(AMRLIB): $(OBJECTS)
	ar rv $(AMRLIB) $(OBJECTS) 


INCLUDE += -DBL_USE_MPI -DBL_USE_OMP -DBL_SPACEDIM=3 -DAMREX_SPACEDIM=3 -DBL_FORT_USE_UNDERSCORE -DBL_Linux

AMReX_AmrTask.o: AMReX_AmrTask.cpp
	$(C++) $(C++FLAGS) -I./ -I../../Base -I../../Amr -I../../AmrCore -I../../Boundary -I../graph -I$(INCLUDE) -c AMReX_AmrTask.cpp -o AMReX_AmrTask.o

.PHONY: clean

clean:
	$(RM) $(OBJECTS)
	$(RM) *.a

