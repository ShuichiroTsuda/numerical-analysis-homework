build:
	gfortran cavity.f90 -o outputs/cavity.out

run:
	outputs/cavity.out && \
	gnuplot cavity_graph.plt