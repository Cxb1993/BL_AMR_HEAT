BL_AMR_HEAT
===========

What's this
-----------
This is an example of Adaptive Mesh Refinement (AMR) computation using BoxLib (LLNL).
This application solves a heat diffusion problem by making use of some features of BoxLib.


How to build
------------

It is assumed that one has a bash-like environment to run Make files.

1. You have to download BoxLib first and save it as a sibling of the BL_AMR_HEAT directory;
2. `cd` into the BL_AMR_HEAT directory and `make`;


How to run
----------
1. execute the program passing it the argument `inputs_2d_amr.properties`;
2. in the shell type

        ls plt*/Headers | tee heat2d.visit
3. use `VisIt` to visualize the results

