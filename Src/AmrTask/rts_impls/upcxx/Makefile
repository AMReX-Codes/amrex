include ../../arch.common 

RTS_LIB= rts.a

OBJECTS= rts.o

all: $(RTS_LIB)

$(RTS_LIB): $(OBJECTS)
	ar rv $(RTS_LIB) $(OBJECTS) 

#$(OBJECTS): rts.C

rts.o: PerillaRts.cpp
	$(C++) $(C++FLAGS) -I. -I../Utils/ -I$(INCLUDE) -c PerillaRts.cpp -o rts.o

.PHONY: clean

clean:
	$(RM) $(OBJECTS)
	$(RM) *.a
