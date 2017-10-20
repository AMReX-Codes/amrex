include ../../arch.common 

RTS_LIB= rts.a

OBJECTS= rts.o sysInfo.o dl_malloc.o

all: $(RTS_LIB)

$(RTS_LIB): $(OBJECTS)
	ar rv $(RTS_LIB) $(OBJECTS) 

#$(OBJECTS): rts.C

rts.o: rts.C
	$(C++) $(C++FLAGS) $(SEGSIZE) -DONLY_MSPACES=1 -I. -I../Utils/ -I$(INCLUDE) -I../../graph -c rts.C -o rts.o

sysInfo.o: ../Utils/sysInfo.C
	$(C++) $(C++FLAGS) -I../Utils/ -I$(INCLUDE) -c ../Utils/sysInfo.C -o sysInfo.o

dl_malloc.o:
	$(CC) -DONLY_MSPACES=1 -I../Utils/ -I$(INCLUDE) -O2 -c ../Utils/dl_malloc.c -o dl_malloc.o

.PHONY: clean

clean:
	$(RM) $(OBJECTS)
	$(RM) *.a
