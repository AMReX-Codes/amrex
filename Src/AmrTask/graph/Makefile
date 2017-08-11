include ../arch.common 

GRAPH_LIB= graph.a

OBJECTS= AMReX_AbstractTask.o AMReX_TaskGraph.o

all: $(GRAPH_LIB)
  
$(GRAPH_LIB): $(OBJECTS) 
	ar rv $(GRAPH_LIB) $(OBJECTS) 

$(OBJECTS): AMReX_AbstractTask.H AMReX_AbstractTask.cpp AMReX_TaskGraph.H AMReX_TaskGraph.cpp

AMReX_AbstractTask.o: AMReX_AbstractTask.cpp AMReX_AbstractTask.H
	$(C++) $(C++FLAGS) -I. -I$(INCLUDE) -c AMReX_AbstractTask.cpp -o AMReX_AbstractTask.o

AMReX_TaskGraph.o: AMReX_TaskGraph.cpp AMReX_TaskGraph.H
	$(C++) $(C++FLAGS) -I./ -I$(INCLUDE) -c AMReX_TaskGraph.cpp -o AMReX_TaskGraph.o

.PHONY: clean

clean:
	$(RM) $(OBJECTS)
	$(RM) graph.a
