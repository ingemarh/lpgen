This is a repo for doing lagprofiling and optionally create an input for GUISDAP initialisation.

cmake + make to make the lp libraries

Make a definition file: manda.py is an example

Generate parameter (and/or t_to_ps) file: ./par_gen.py -s site [manda]

Test the par file: ./pltest.py [-p] [manda_va.par] [1023465677.mat]
