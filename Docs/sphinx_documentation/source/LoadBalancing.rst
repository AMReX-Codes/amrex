.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

.. _ss:load_balancing:

Load Balancing
--------------

The process of load balancing is independent of the process of grid creation; 
the inputs to load balancing are a given set of grids with a set of weights 
assigned to each grid.  

Single-level load balancing algorithms are sequentially applied to each AMR level independently, 
and the resulting distributions are mapped onto the ranks taking into account the weights 
already assigned to them (assign heaviest set of grids to the least loaded rank)

Options supported by AMReX include:

- Knapsack: the default weight of a grid in the knapsack algorithm is the number of grid cells, 
  but AMReX supports the option to pass an array of weights – one per grid – or alternatively 
  to pass in a MultiFab of weights per cell which is used to compute the weight per grid

- SFC: enumerate grids with a space-filling Z-morton curve, then partition the 
  resulting ordering across ranks in a way that balances the load

- Round-robin: sort grids and assign them to ranks in round-robin fashion -- specifically
  FAB *i* is owned by CPU *i%N* where N is the total number of MPI ranks.
