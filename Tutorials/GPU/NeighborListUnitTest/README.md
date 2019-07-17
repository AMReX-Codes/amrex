This unit test

(1) initializes particles, 

(2) counts how many particles with which grid id it "owns" (only for grid 0) -- here it only owns the
    particles it started with, e.g. 512 = 8^3 if num_ppc = 1

(3) calls FillNeighbors

(4) counts how many particles with which grid id it "owns" (only for grid 0) -- here it should have
    its original particles plus particles from every other grid, e.g. 8^3 , 6x8^2 , 12x8 , 8x1 from grids 
    0 to 26 (if num_ppc = 1)

(5) calls UpdateNeighbors

(6) counts how many particles with which grid id it "owns" (only for grid 0) -- answer should be identical to previous call

(7) calls reset_test_id which sets the "test_id" of every particle (including "ghosts") to the (local)grid id 

(8) counts how many particles with which grid id it "owns" (only for grid 0) -- should be 
    e.g. 1000 = 8^3 + 6x8^2 + 12x8 + 8x1  (if num_ppc = 1) from grid 0 and no others

(9) calls UpdateNeighbors

(10) counts how many particles with which grid id it "owns" (only for grid 0) -- answer should revert back to that in (4)
