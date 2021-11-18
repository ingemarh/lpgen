This is a repo for doing lagprofiling and optionally create an input for GUISDAP initialisation.

Generate parameter (and/or t2ps) file:
Make a definition file: manda.py is an example
	python3:	import par_gen as pg
			pg.gen('manda',['v'])
	shell:	./par_gen.py -s site /manda/

Test the lagprofiler:
	cmake . + make to make the lp libraries
	Test the par file (-p plot the result):
	Test the par file (-w save the result):
	python3:	import partest as pt
			pt.pltest('manda_va.par',['p'],['15234489.mat'])
	shell:	./partest.py [-p] [-w] /manda_va.par/ [15234489.mat]
