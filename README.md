# proteinPacking
Calculate packing fraction of residues in protein structures using voro++ (C. H. Rycroft, Voro++: A three-dimensional Voronoi cell library in C++, Chaos 19, 041111 (2009)).

- Compile code with voro++
 g++ -I src/ -I src/voro++/src/ [SOURCE CODE HERE] src/*.cpp src/voro++/src/voro++.cc -o pf.o
