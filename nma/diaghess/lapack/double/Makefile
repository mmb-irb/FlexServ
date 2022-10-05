F77=gfortran
AR=ar

all: dlanst.f dlarf.f dlascl.f dlassq.f dorgql.f dsterf.f dlae2.f dlansy.f dlarfg.f dlaset.f dlatrd.f dorgqr.f dsyev.f dlaev2.f dlapy2.f dlarft.f dlasr.f dorg2l.f dorgtr.f dsytd2.f dlamch.f dlarfb.f dlartg.f dlasrt.f dorg2r.f dsteqr.f dsytrd.f
	$(F77) -c *.f
	$(AR) -r libdouble.a *.o

clean:
	$(RM) -f libdouble.a *.o

