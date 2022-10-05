F77=gfortran
AR=ar

all: daxpy.f ddot.f dgemv.f dnrm2.f dswap.f dsyr2.f dtrmm.f dcopy.f dgemm.f dger.f dscal.f dsymv.f dsyr2k.f dtrmv.f
	$(F77) -c *.f
	$(AR) -r libblas.a *.o

clean:
	$(RM) -f libblas.a *.o

