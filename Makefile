GMPINC ?= /usr/include/x86_64-linux-gnu
all: src/c/fsearch64 src/c/fsearch_gmp src/c/pepin_fft src/c/bench_sqr
src/c/fsearch64: src/c/fsearch64.c
	gcc -O3 -march=native -fopenmp -o $@ $<
src/c/fsearch_gmp: src/c/fsearch_gmp.c
	gcc -O3 -march=native -fopenmp -I$(GMPINC) -o $@ $< -lgmp
src/c/pepin_fft: src/c/pepin_fft.c
	gcc -O2 -I$(GMPINC) -o $@ $< -lgmp
test: all
	python3 -m pytest -q tests/
src/c/bench_sqr: src/c/bench_sqr.c
	gcc -O2 -I$(GMPINC) -o $@ $< -lgmp
