build:
	gfortran cavity.f90 -o outputs/build/cavity.out

run:
	outputs/build/cavity.out && \
	gnuplot cavity_graph.plt