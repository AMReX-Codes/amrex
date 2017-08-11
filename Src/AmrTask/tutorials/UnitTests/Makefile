include ../../arch.common 

OBJECTS= 001_TokenRing.o 002_Jacobi_StaticGraph.o

all: 001_TokenRing 002_Jacobi_StaticGraph

001_TokenRing: 001_TokenRing.o
	$(C++LINK) $(C++FLAGS) 001_TokenRing.o ../../graph/graph.a ../../rts_impls/Serial/rts.o -o 001_TokenRing 

001_TokenRing.o: 001_TokenRing.C
	$(C++) $(C++FLAGS) -I$(INCLUDE) -I. -I../../graph -c 001_TokenRing.C -o 001_TokenRing.o

002_Jacobi_StaticGraph: 002_Jacobi_StaticGraph.o
	$(C++LINK) $(C++FLAGS) 002_Jacobi_StaticGraph.o ../../graph/graph.a ../../rts_impls/Serial/rts.o -o 002_Jacobi_StaticGraph

002_Jacobi_StaticGraph.o: 002_Jacobi_StaticGraph.C
	$(C++) $(C++FLAGS) -I$(INCLUDE) -I. -I../../graph -c 002_Jacobi_StaticGraph.C -o 002_Jacobi_StaticGraph.o

.PHONY: clean

clean:
	$(RM) $(OBJECTS)
	$(RM) 001_TokenRing
	$(RM) 002_Jacobi_StaticGraph
