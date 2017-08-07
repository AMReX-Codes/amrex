include ../../arch.common 

OBJECTS= 001_TokenRing.o

all: 001_TokenRing

001_TokenRing: $(OBJECTS)
	$(C++LINK) $(C++FLAGS) 001_TokenRing.o ../../graph/graph.a ../../rts_impls/Serial/rts.o -o 001_TokenRing 
  
$(OBJECTS): 001_TokenRing.C

001_TokenRing.o: 001_TokenRing.C
	$(C++) $(C++FLAGS) -I. -I../../graph -c 001_TokenRing.C -o 001_TokenRing.o

.PHONY: clean

clean:
	$(RM) $(OBJECTS)
	$(RM) 001_TokenRing
