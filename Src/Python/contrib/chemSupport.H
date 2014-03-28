
#include <list>
#include <ChemDriver.H>

std::list<ChemDriver::Edge>
GetPathDiagramEdges (const ChemDriver&  cd,
		     const std::string& trace_element_name,
		     bool               unmap_edge_reactions         = true,
		     bool               one_based_reaction_numbering = true,
		     int                verbose                      = false,
		     bool               hack_mech_splitting_rules    = true);
