This unit test

(1) initializes particles, 

(2) calls FillNeighbors to create particle copies in ghost cells of each grid

(3) calls BuildNeighborList to create particle-neighbor lists 

(4) calls checkNeighbors to compare 
    A.  how many particles each particle has in its neighbor list to 
        how many particles are within the cutoff radius (the same radius used to compute the neighbor list)
        -- computed using an N^2 loop over all particles "owned" by this grid
    B.  if the particle neighbor counts match, do the id's of the neighbors match?
