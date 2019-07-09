%: %.f
	gfortran -static -O3 -o $@ $<

all: fillfoottable headonslice
