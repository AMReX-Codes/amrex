include ../../arch.common 

OBJECTS= rts.o

all: $(OBJECTS)
  
$(OBJECTS): rts.C

rts.o: rts.C
	$(C++) $(C++FLAGS) -I. -I$(INCLUDE) -I../../graph -c rts.C -o rts.o

.PHONY: clean

clean:
	$(RM) $(OBJECTS)
