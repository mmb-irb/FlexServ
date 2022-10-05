F77=gfortran
AR=ar

all: ieeeck.f ilaenv.f lsame.f xerbla.f
	$(F77) -c *.f
	$(AR) -r libutil.a *.o

clean:
	$(RM) -f libutil.a *.o

