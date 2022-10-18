F77=gfortran
OPTIONS=-static-libgfortran -static-libgcc -ffixed-line-length-132

all: bdbin dmdbin nmabin

bdbin: bdrandom.o bd/bd2.f
	cd bd && \
	$(F77) $(OPTIONS) bd2.f bdrandom.o -o bd

bdrandom.o: bd/bdrandom.f
	cd bd && \
	$(F77) $(OPTIONS) -c bdrandom.f -o bdrandom.o

dmdbin: dmd/dmdgoopt.f
	cd dmd && \
	$(F77) $(OPTIONS) dmdgoopt.f -o dmd

nmabin: nma/lorellnma.f
	cd nma && \
	$(F77) $(OPTIONS) lorellnma.f -o nma

clean:
	rm bd/bdrandom.o bd/bd dmd/dmd nma/nma