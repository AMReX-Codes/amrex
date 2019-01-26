
PROGS = test

CXX = CC
CXXFLAGS = -Wall -O3 -I../
LDFLAGS  = 
LDLIBS   = 

########## SHOULD NOT CHANGE ANYTHING BELOW ##############

CPP_SOURCES = ../MemInfo.cpp testMemInfo.cpp
INCLUDE = ../MemInfo.H

CPP_OBJ = $(CPP_SOURCES:.cpp=.o)


default: $(CPP_SOURCES) $(PROGS)
	rm -rf *.o

.PHONY: clean tarball
clean:
	rm -rf *.o ../*.o $(PROGS)

tarball:
	tar -zcf source.tar.gz $(CPP_SOURCES) $(INCLUDE) Makefile

test: $(CPP_OBJ)
	$(CXX) $(LDFLAGS) $(CPP_OBJ) -o $@ $(LDLIBS)

.cpp.o: $(INCLUDE)
	$(CXX) -c $(CXXFLAGS) $< -o $@

