include ../arch.common 

OBJECTS= AMReX_AsyncMFIter.o 

AMFIterLIB= AMReX_AsyncMFIter.a

all: $(AMFIterLIB)

$(AMFIterLIB): $(OBJECTS)
	ar rv $(AMFIterLIB) $(OBJECTS) 


INCLUDE += -DBL_USE_MPI -DBL_USE_OMP -DBL_SPACEDIM=3 -DAMREX_SPACEDIM=3 -DBL_FORT_USE_UNDERSCORE -DBL_Linux

all: $(OBJECTS)
  
AMReX_AsyncMFIter.o: AMReX_AsyncMFIter.cpp AMReX_AsyncMFIter.H
	$(C++) $(C++FLAGS) -I./ -I../../Base -I../../Amr -I../../AmrCore -I../graph -I$(INCLUDE) -c AMReX_AsyncMFIter.cpp -o AMReX_AsyncMFIter.o

.PHONY: clean

clean:
	$(RM) $(OBJECTS)
	$(RM) *.a

