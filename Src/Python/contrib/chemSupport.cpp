#include <iostream>

#include <chemSupport.H>

std::list<ChemDriver::Edge>
GetPathDiagramEdges (const ChemDriver&  cd,
		     const std::string& trace_element_name,
		     bool               unmap_edge_reactions,
		     bool               one_based_reaction_numbering,
		     int                verbose,
		     bool               hack_mech_splitting_rules)
{
  std::list<ChemDriver::Edge> edges;

  if (verbose) {
    std::cout << "Getting path diagram edges ("
	      << "Element: " << trace_element_name << ", " 
	      << "Unmap edge reactions: " << unmap_edge_reactions << ", " 
	      << "One-based reaction numbering: " << one_based_reaction_numbering << ", " 
	      << "Verbose edge generation: " << verbose << ", " 
	      << "Apply special reaction splitting rules: " << hack_mech_splitting_rules << ")" << std::endl;
  }

  if (unmap_edge_reactions) {

    // Build inverse reaction map
    //
    // Given:
    //     map[old_number] = new_number
    //
    // Construct
    //     map[new_number] = old_number
    //
    const Array<int>& reaction_map = cd.reactionMap();
    Array<int> new_to_old(cd.numReactions());
    for (int i=0; i<cd.numReactions(); ++i) {
      new_to_old[reaction_map[i]] = i;
    }
    
    // Translate edges from mapped reactions
    // 
    // For each reaction in the edges rate weight list, map the list that is
    // specified in terms of the new numbering into a list in terms of the new numbering
    //
    std::list<ChemDriver::Edge> edges_new = cd.getEdges(trace_element_name, verbose, hack_mech_splitting_rules);

    for (std::list<ChemDriver::Edge>::const_iterator it=edges_new.begin(), End=edges_new.end(); it!=End; ++it) {
      const Array<std::pair<int,Real> >& RWL_new = it->RateWeightList();
      int nr = RWL_new.size();
      Array<std::pair<int,Real> > RWL_old(nr);
      for (int i=0; i<nr; ++i) {
	const std::pair<int,Real>& rw_new = RWL_new[i];
	int new_rxn_id = rw_new.first;
	int old_rxn_id = new_to_old[new_rxn_id];
	if (one_based_reaction_numbering) old_rxn_id++;
	int rxn_weight = rw_new.second;
	RWL_old[i] = std::pair<int,Real>(old_rxn_id,rxn_weight);
      }
      edges.push_back(ChemDriver::Edge(it->left(),it->right(),RWL_old));
    }
  }
  else {
    edges = cd.getEdges(trace_element_name, verbose, hack_mech_splitting_rules);
  }

  return edges;
}
