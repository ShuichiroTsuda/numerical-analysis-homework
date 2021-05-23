build:
	gfortran cavity.f90 -o cavity.out

run:
	./cavity.out && \
	gnuplot cavity_graph.plt