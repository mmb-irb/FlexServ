      FUNCTION dran_u()
      IMPLICIT real(a-h,o-z)
      PARAMETER(ip=1279)
      PARAMETER(iq=418)
      PARAMETER(is=ip-iq)
      PARAMETER (rmax=2147483647.0)
      INTEGER ix(ip)
      COMMON /ixx/ ix
      COMMON /icc/ ic
      ic=ic+1
      IF(ic.gt.ip) ic=1
      IF(ic.gt.iq) THEN
        ix(ic)=ieor(ix(ic),ix(ic-iq))
      ELSE
        ix(ic)=ieor(ix(ic),ix(ic+is))
      ENDIF
      dran_u=real(ix(ic))/rmax
      RETURN
      END
C
      SUBROUTINE dran_gv(u,n)
C
      IMPLICIT real*8 (a-h,o-z)
      PARAMETER(ip=1279)
      PARAMETER(iq=418)
      PARAMETER(np=14)
      PARAMETER(nbit=31)
      PARAMETER(m=2**np,np1=nbit-np,nn=2**np1-1,nn1=nn+1)
      PARAMETER(is=ip-iq)
      DIMENSION g(0:m)
      DIMENSION u(n)
      DIMENSION ix(ip)
      COMMON /ixx/ ix
      COMMON /icc/ic
      COMMON /gg/ g
C
      DO 99 k=1,n
        ic=ic+1
        IF(ic.gt.ip) ic=1
        IF(ic.gt.iq) THEN
C
          ix(ic)=ieor(ix(ic),ix(ic-iq))
        ELSE
          ix(ic)=ieor(ix(ic),ix(ic+is))
        ENDIF
        i=ishft(ix(ic),-np1)
        i2=iand(ix(ic),nn)
        u(k)=i2*g(i+1)+(nn1-i2)*g(i)
   99 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE dran_ini(iseed0)
      IMPLICIT double precision (a-h,o-z)
      PARAMETER(ip=1279)
      PARAMETER(np=14)
      PARAMETER(nbit=31)
      PARAMETER(m=2**np,np1=nbit-np,nn=2**np1-1,nn1=nn+1)
      INTEGER ix(ip)
      DIMENSION g(0:m)
C
      DATA c0,c1,c2/2.515517,0.802853,0.010328/
      DATA d1,d2,d3/1.432788,0.189269,0.001308/
C
      COMMON /ixx/ ix
      COMMON /icc/ ic
      COMMON /gg/ g
C
      dseed=iseed0
      DO 200 i=1,ip
        ix(i)=0
        DO 200 j=0,nbit-1
          IF(rand_xx(dseed).lt.0.5) ix(i)=ibset(ix(i),j)
  200 CONTINUE
      ic=0
C
      pi=4.0d0*datan(1.0d0)
      DO 1 i=m/2,m
        p=1.0-real(i+1)/(m+2)
        t=sqrt(-2.0*log(p))
        x=t-(c0+t*(c1+c2*t))/(1.0+t*(d1+t*(d2+t*d3)))
        g(i)=x
        g(m-i)=-x
    1 CONTINUE
C
      u2th=1.0-real(m+2)/m*sqrt(2.0/pi)*g(m)* exp(-g(m)*g(m)/2)
      u2th=nn1*sqrt(u2th)
      DO 856 i=0,m
  856 g(i)=g(i)/u2th
C
      RETURN
      END
C
      SUBROUTINE ranvel(nsel,npm,nrp,v,winv,ibelly,temp0,ig,iscr)
C
C Subroutine RANdom VELocities. This routine randomly selects NSEL atoms and resets their
C their velocities to random values with probabilities given by a maxwellian distribution
C about TEMP0
C
C NSEL... Number of randomly chosen atoms whose velocities will be randomly reassigned. If =
C 0, then all will be reassigned.
C NPM... number of replicates. Usually 1.
C NRP... number of atoms
C V()... velocity array
C WINV... inverse mass array
C IBEL()... if IBELLY>0, only atoms for wich IBEL(I)>0 will move
C IBELLY... >0 when belly being used.
C TEMP0... target temperature
C IG... seed for tandom number generator.
C ISCR()... scratch array at least NPM*NRP elements long.
C
      IMPLICIT double precision (a-h,o-z)
C
C Read coordinates and velocities
C
      DIMENSION v(nrp*3),winv(nrp),iscr(npm*nrp)
C
      iseed=ig
      boltz=(8.31441d-3/4.814d0)*temp0
      ntotat=npm*nrp
C
C Clear scratch array, after any velocity is reset, ISCR(I) is set to 1. ISCR(I) is set to 1
C immediately (don't assign velocity) if IBEL(I)=0 and belly is on.
C
      DO 10 i=1,ntotat
C
        iscr(i)=0
C
C IF(IBELLY.GT.0 .AND. IBEL(I).EQ.0) ISCR(I)=1
C
   10 CONTINUE
C
C Do NSEL velocity reassignments
C
      nchan=nsel
      IF(nchan.gt.ntotat.or.nchan.eq.0) nchan=ntotat
C
      DO 50 i=1,nchan
C
C Generate a random number on the interval 1->NTOTAT. Then search from that location
C forward, looking for the first still-untouched velocity. (Loop around if we get to end of
C list before we find untouched velocity).
C
        CALL amrand(y)
        isel=min(int(y*float(ntotat))+1,ntotat)
        DO 60 j=isel,isel+ntotat-1
          jpt=mod(j-1,ntotat)+1
          IF(iscr(jpt).eq.0) GOTO 65
   60   CONTINUE
C
C If we get here, no atoms left to change, RETURN
C
        RETURN
C
   65   ifill=3*(jpt-1)
        sd=sqrt(boltz*winv(jpt))
        DO 70 m=1,3
          CALL gauss(0.d0,sd,v(ifill+m))
   70   CONTINUE
        iscr(jpt)=1
   50 CONTINUE
      RETURN
      END
C
      SUBROUTINE amrand(y)
C
      DOUBLE PRECISION u(97),c,cd,cm
C
C ... real variables in Marsaglia algorithm
C
      INTEGER i97,j97
C
C ... pointers into U() in Marsaglia algorithm
C
      LOGICAL set
C
C ... true if amrset has been called
C
      COMMON /raset1/ u,c,cd,cm,i97,j97,set
C
C      OUTPUT:
C
C      Y: A random number between 0.0 and 1.0
      DOUBLE PRECISION y
      DOUBLE PRECISION uni
C
      IF(.not.set) THEN
        WRITE(6,'(a)') 'amrand not initd'
C
C call mexit(6,1)
C
      ENDIF
C
      uni = u(i97)-u(j97)
      IF (uni.lt.0.0d0) uni=uni+1.0d0
      u(i97)=uni
      i97=i97-1
      IF(i97.eq.0) i97=97
      j97=j97-1
      IF (j97.eq.0) j97=97
      c=c-cd
      IF(c.lt.0.0d0) c=c+cm
      uni=uni-c
      IF(uni.lt.0.0d0)uni=uni+1.0d0
C
      y=uni
C
      RETURN
      END
C
      SUBROUTINE gauss(am,sd,v)
C
      IMPLICIT double precision(a-h,o-z)
C
C  this is a version of amrand() that adds the constraint of a gaussian distribution. It also
C  requires amrset to have been called first, and 'uses up' the same sequence that amrand()
C  does.
C
      DOUBLE PRECISION u(97),c,cd,cm
C
C ... real variables in Marsaglia algorithm
C
      INTEGER i97,j97
C
C ... pointers into u() in Marsaglia algorithm
C
      LOGICAL set
C
C ... true if amrset has been called
C
      COMMON /raset1/u,c,cd,cm,i97,j97,set
      DATA zero,six /0.0d0,6.0d0/
C
      IF (.not.set) THEN
        WRITE(6,'(a)') 'amrand not initd'
C
C call mexit(6, 1)
C
      ENDIF
      a=zero
      DO 50 i=1,12
        uni=u(i97)-u(j97)
        IF(uni.lt.0.0d0) uni=uni+1.0d0
        u(i97)=uni
        i97=i97-1
        IF(i97.eq.0)i97=97
        j97=j97-1
        IF(j97.eq.0)j97=97
        c=c-cd
        IF(c.lt.0.0d0)c=c+cm
        uni=uni-c
        IF(uni.lt.0.0d0) uni=uni+1.0d0
        a=a+uni
   50 CONTINUE
      v=(a-six)*sd+am
      RETURN
      END
C
      SUBROUTINE amrset(iseed)
C
      INTEGER iseed
C
C ... integer seed greater than zero
C
      DOUBLE PRECISION u(97),c,cd,cm
      INTEGER i97,j97
      LOGICAL set
C
      COMMON /raset1/u,c,cd,cm,i97,j97,set
      INTEGER is1,is2
C
C ... the two internal seeds used in Marsaglia algorithm
C
      INTEGER is1max
C
C ... max value of first seed (is1), 31328
C
      INTEGER is2max
C
C ... max value of second seed (is2), 30081
C
      INTEGER i,j,k,l,m
C
C ... used in generation of u()
C
      DOUBLE PRECISION s,t
C
C ... used in generation of u()
C
      INTEGER ii,jj
C
C ... loop indices
C
      DATA is1max, is2max /31328, 30081/
C
      is1=max((iseed/is2max)+1,1)
      is1=min(is1,is1max)
C
      is2=max(1,mod(iseed,is2max)+1)
      is2=min(is2,is2max)
C
      i=mod(is1/177,177)+2
      j=mod(is1,177)+2
      k=mod(is2/169,178)+1
      l=mod(is2,169)
C
      DO 200 ii=1,97
        s=0.0d0
        t=0.5d0
        DO 100 jj=1,24
          m=mod(mod(i*j,179)*k,179)
          i=j
          j=k
          k=m
          l=mod(53*l+1,169)
          IF(mod(l*m,64).ge.32)s=s+t
          t=0.5d0*t
  100   CONTINUE
        u(ii)=s
  200 CONTINUE
C
      c=362436.0d0/16777216.0d0
      cd=7654321.0d0/16777216.0d0
      cm=16777213.0d0/16777216.0d0
C
      i97=97
      j97=33
C
      set=.true.
      RETURN
      END
C
      FUNCTION rand_xx(dseed)
      DOUBLE PRECISION a,c,xmm,rm,dseed,rand_xx
      PARAMETER (xmm=2.d0**32,rm=1.d0/xmm,a=69069.d0,c=1.d0)
      dseed=mod(dseed*a+c,xmm)
      rand_xx=dseed*rm
      RETURN
      END
