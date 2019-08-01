# Matching

A match.txt file is generated by the file `test_squares_meanshift.cpp`. This file is currently set up to operate on
2D point clouds. It is built up like this:

* approxmatch_cpu takes two point clouds and generates a match matrix as well as two offset vectors, one for each point cloud
* it calls calc_offset that is defined in sort_indices.c
* in sort_indices there are two implementations, one takes fixed neighbours, the other a fixed window
* only a fixed window works, if we shift with only a fixed number of neighbours, the shift does typically not shift the entire object or it does not shift only one object, but multiple to the origin
* the disadvantage of a fixed window is that you have to guess beforehand the extend of an individual object

Build:

	mkdir -p build && cd build && cmake .. && make

Run in the `build` directory:

	./emd-suite

This will generate several `.txt` files.

## Visualization

The `display_match.m` octave script takes two point clouds and a matching file.

The matches are subsequently visualized by drawing lines of various thickness between the two point clouds.

Run in the `analyse` directory:

	./display_match.m ../build/cloud1.txt ../build/cloud2.txt ../build/match.txt

# Copyrights

Copyrights belong to A.C. van Rossum (2019)

License: LGPLv3, MIT, Apache License
