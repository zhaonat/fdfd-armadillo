# fdfd-armadillo
implementation of fdfd in C++

## Dependencies
The only dependency is armadillo

## Structure
for a 2D solver, there is no classes that are needed, we just have a solveTE and solveTM.
There are just a few supporting functions which are needed, such as generating derivative operators, PMLs, averaging, but these are simple and can 
be unit-tested as well.

## Worfkow
Because we having to use an executable which comes out of a compiled src file. The input structure, source, and information such as the domain size and number of grid cells have to be specified via input files.

Parameters that you have to specify:
