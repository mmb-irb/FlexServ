F77=gfortran
AR=ar
OPTIONS=-static-libgfortran -static-libgcc
LINE_LENGTH_132=-ffixed-line-length-132
LINE_LENGTH_NONE=-ffixed-line-length-none
all: bd_ dmdgoopt_ lorellnma_ diaghess_

bd_: bdrandom.o bd/bd2.f
	cd bd && \
	$(F77) $(OPTIONS) $(LINE_LENGTH_NONE) bd2.f bdrandom.o -o bd

bdrandom.o: bd/bdrandom.f
	cd bd && \
	$(F77) $(OPTIONS) $(LINE_LENGTH_NONE) -c bdrandom.f -o bdrandom.o

dmdgoopt_: dmd/dmdgoopt.f
	cd dmd && \
	$(F77) $(OPTIONS) $(LINE_LENGTH_NONE) dmdgoopt.f -o dmdgoopt

lorellnma_: nma/lorellnma.f
	cd nma && \
	$(F77) $(OPTIONS) $(LINE_LENGTH_132) lorellnma.f -o lorellnma

diaghess_: nma/diaghess/diaghess.f libdouble.a libutil.a libblas.a
	cd nma/diaghess && \
	$(F77) $(OPTIONS) diaghess.f lapack/double/libdouble.a lapack/util/libutil.a blas/libblas.a -o diaghess

libdouble.a: nma/diaghess/lapack/double/dlanst.f nma/diaghess/lapack/double/dlarf.f nma/diaghess/lapack/double/dlascl.f nma/diaghess/lapack/double/dlassq.f nma/diaghess/lapack/double/dorgql.f nma/diaghess/lapack/double/dsterf.f nma/diaghess/lapack/double/dlae2.f nma/diaghess/lapack/double/dlansy.f nma/diaghess/lapack/double/dlarfg.f nma/diaghess/lapack/double/dlaset.f nma/diaghess/lapack/double/dlatrd.f nma/diaghess/lapack/double/dorgqr.f nma/diaghess/lapack/double/dsyev.f nma/diaghess/lapack/double/dlaev2.f nma/diaghess/lapack/double/dlapy2.f nma/diaghess/lapack/double/dlarft.f nma/diaghess/lapack/double/dlasr.f nma/diaghess/lapack/double/dorg2l.f nma/diaghess/lapack/double/dorgtr.f nma/diaghess/lapack/double/dsytd2.f nma/diaghess/lapack/double/dlamch.f nma/diaghess/lapack/double/dlarfb.f nma/diaghess/lapack/double/dlartg.f nma/diaghess/lapack/double/dlasrt.f nma/diaghess/lapack/double/dorg2r.f nma/diaghess/lapack/double/dsteqr.f nma/diaghess/lapack/double/dsytrd.f
	cd nma/diaghess/lapack/double && \
	$(F77) $(OPTIONS) -c *.f && \
	$(AR) -r libdouble.a *.o

libutil.a: nma/diaghess/lapack/util/ieeeck.f nma/diaghess/lapack/util/ilaenv.f nma/diaghess/lapack/util/lsame.f nma/diaghess/lapack/util/xerbla.f
	cd nma/diaghess/lapack/util && \
	$(F77) $(OPTIONS) -c *.f && \
	$(AR) -r libutil.a *.o

libblas.a: nma/diaghess/blas/daxpy.f nma/diaghess/blas/ddot.f nma/diaghess/blas/dgemv.f nma/diaghess/blas/dnrm2.f nma/diaghess/blas/dswap.f nma/diaghess/blas/dsyr2.f nma/diaghess/blas/dtrmm.f nma/diaghess/blas/dcopy.f nma/diaghess/blas/dgemm.f nma/diaghess/blas/dger.f nma/diaghess/blas/dscal.f nma/diaghess/blas/dsymv.f nma/diaghess/blas/dsyr2k.f nma/diaghess/blas/dtrmv.f
	cd nma/diaghess/blas && \
	$(F77) $(OPTIONS) -c *.f && \
	$(AR) -r libblas.a *.o

clean:
	rm bd/bdrandom.o bd/bd dmd/dmdgoopt nma/lorellnma nma/diaghess/diaghess nma/diaghess/lapack/double/*.o nma/diaghess/lapack/double/libdouble.a nma/diaghess/lapack/util/*.o nma/diaghess/lapack/util/libutil.a nma/diaghess/blas/*.o nma/diaghess/blas/libblas.a