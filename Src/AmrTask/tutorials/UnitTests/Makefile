include ../../arch.common 

OBJECTS= 001_TokenRing.o 002_Jacobi_StaticGraph.o 003_Jacobi_DynamicGraph.o

all: 001_TokenRing 002_Jacobi_StaticGraph 003_Jacobi_DynamicGraph

001_TokenRing: 001_TokenRing.o
	$(C++LINK) $(C++FLAGS) 001_TokenRing.o ../../graph/graph.a $(RTS_DIR)/rts.a $(LDLIBS) -o 001_TokenRing

001_TokenRing.o: 001_TokenRing.C
	$(C++) $(C++FLAGS) -I$(INCLUDE) -I. -I../../graph -c 001_TokenRing.C -o 001_TokenRing.o

002_Jacobi_StaticGraph: 002_Jacobi_StaticGraph.o
	$(C++LINK) $(C++FLAGS) 002_Jacobi_StaticGraph.o ../../graph/graph.a $(RTS_DIR)/rts.a $(LDLIBS)  -o 002_Jacobi_StaticGraph

002_Jacobi_StaticGraph.o: 002_Jacobi_StaticGraph.C
	$(C++) $(C++FLAGS) -I$(INCLUDE) -I. -I../../graph -c 002_Jacobi_StaticGraph.C -o 002_Jacobi_StaticGraph.o

003_Jacobi_DynamicGraph: 003_Jacobi_DynamicGraph.o
	$(C++LINK) $(C++FLAGS) 003_Jacobi_DynamicGraph.o ../../graph/graph.a $(RTS_DIR)/rts.a $(LDLIBS)  -o 003_Jacobi_DynamicGraph

003_Jacobi_DynamicGraph.o: 003_Jacobi_DynamicGraph.C
	$(C++) $(C++FLAGS) -I$(INCLUDE) -I. -I../../graph -c 003_Jacobi_DynamicGraph.C -o 003_Jacobi_DynamicGraph.o
.PHONY: clean

clean:
	$(RM) $(OBJECTS)
	$(RM) 001_TokenRing
	$(RM) 002_Jacobi_StaticGraph
	$(RM) 003_Jacobi_DynamicGraph
