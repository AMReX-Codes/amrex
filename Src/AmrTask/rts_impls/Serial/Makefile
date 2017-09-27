include ../../arch.common 

OBJECTS= rts.o

RTS_LIB= rts.a

OBJECTS= rts.o

all: $(RTS_LIB)

$(RTS_LIB): $(OBJECTS)
	ar rv $(RTS_LIB) $(OBJECTS) 

$(OBJECTS): rts.C

rts.o: rts.C
	$(C++) $(C++FLAGS) -I. -I$(INCLUDE) -I../../graph -c rts.C -o rts.o

.PHONY: clean

clean:
	$(RM) $(OBJECTS)
