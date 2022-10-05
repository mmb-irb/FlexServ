#F77=g77
F77=gfortran

all: diaghess.f
	(cd blas; make; cd ..)
	(cd lapack; make; cd ..)
	$(F77) -o diaghess diaghess.f lapack/double/libdouble.a lapack/util/libutil.a blas/libblas.a 

clean:
	(cd blas; make clean; cd ..)
	(cd lapack; make clean; cd ..)
	$(RM) -f diaghess *.o

