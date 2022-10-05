      PROGRAM verlet_brownian_md
C
      IMPLICIT none
C
      CHARACTER*255 gridfn1,gridfn2,gridfn3,gridfn4,gridfn5,gridfn6,gridfn7
      CHARACTER*50 atom,ca
      CHARACTER*50 fitxer
      CHARACTER*7 auxstr
C
C natomstot són el nombre d'àtoms del pdb inicial (no considero un pdb amb només els
C carbonis alfa, sinó que considero el de tota la proteïna sense treure'n res)
C
	integer natomsmax
	real pi
      PARAMETER (natomsmax=1500, pi=3.1415926)
C
      real fx(natomsmax),fy(natomsmax),fz(natomsmax)
      real*8 randx(natomsmax),randy(natomsmax),randz(natomsmax)
      real rx(natomsmax),ry(natomsmax),rz(natomsmax)
      real vx(natomsmax),vy(natomsmax),vz(natomsmax)
      real*8 v(natomsmax*3)
      real winv(natomsmax)
      real c(natomsmax,natomsmax),d(natomsmax,natomsmax)
      real d0(natomsmax,natomsmax)
      integer iscr(natomsmax)
      real vxmig(natomsmax),vymig(natomsmax),vzmig(natomsmax)
      real k0(natomsmax,natomsmax)
	real ar,br,cr,c0,c1,bvmig,pmass,t0,dt,const,r0,cv,bv
	real 	epsi, c2, tau, taup, boltz, theta, phi
c
      integer pairs(2,natomsmax*natomsmax),npairs,k
	integer i,j,n,m,natoms, itsnap, iseed, itempsmax
	real*8  dran_u
      real CutOff
c
      CutOff = 12.
C
C llegir condicions inicials d'alguna banda
C
      CALL getarg(1,gridfn1)
      CALL getarg(2,gridfn2)
      CALL getarg(3,gridfn3)
      CALL getarg(4,gridfn4)
      CALL getarg(5,gridfn5)
      CALL getarg(6,gridfn6)
      CALL getarg(7,gridfn7)
C
      READ(gridfn2,*)itempsmax
      READ(gridfn3,*)dt
      READ(gridfn4,*)itsnap
      READ(gridfn5,*)const
      READ(gridfn6,*)r0
C
	call readPDB (gridfn1, rx,ry,rz,natoms)
	call moveToOrigin (rx,ry,rz,natoms)
	call calcMatDist (natomsmax, natoms,rx,ry,rz,d0)
	call calcPairs (natomsmax,natoms,d0,CutOff,pairs,npairs)
C
      DO k=1,npairs
	 k0(pairs(1,k),pairs(2,k))=-const * (r0/d0(pairs(1,k),pairs(2,k)))**6
	 k0(pairs(2,k),pairs(1,k))=k0(pairs(1,k),pairs(2,k))
       ENDDO
C massa mitjana d'un residu és + o - 100.d0
C
      pmass=100.d0
      t0=300.d0
C
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   Velocitats inicials aleatòries a una certa temperatura T0
C
      DO i=1,natoms
        winv(i)=1.d0/pmass
      ENDDO
C
      iseed=73
      CALL amrset(iseed)
      CALL ranvel(0,1,natomsmax,v,winv,1,t0,iseed,iscr)
C
      iseed=94
      CALL dran_ini(iseed)
C
      DO i=1,natoms
        theta=pi*dran_u()
        phi=2.*pi*dran_u()
C
        vx(i)=abs(v(i))*cos(theta)*cos(phi)*1.e10
        vy(i)=abs(v(i))*cos(theta)*sin(phi)*1.e10
        vz(i)=abs(v(i))*sin(theta)*1.e10
      end do
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      t0=300.
      boltz=(8.31441e-3/4.814)
      pmass=100.*1.67e-27
      tau=0.4e-12
      taup=tau/pmass
      epsi=2.*boltz*t0/taup
C
      c0=exp(-dt/tau)
      c1=(1.-c0)*tau/dt
      c2=(1.-c1)*tau/dt
C
      ar=c1*dt
      br=c2*0.5*dt**2/pmass
      bvmig=(c1-c2)*dt/pmass
      bv=c2*dt*pmass
      cr=sqrt(dt)*(1.-c0)*tau/pmass*sqrt(epsi)
      cv=c0*sqrt(dt)/pmass*sqrt(epsi)
C
      DO i=1,natoms
        fx(i)=0.
        fy(i)=0.
        fz(i)=0.
	do j=1,natoms
	c(i,j) = 0.
	end do
	end do
C
      write (6,*) "Beginning trajectory for ", natoms, " atoms"
      OPEN(unit=34,file=gridfn7,status='REPLACE')
      WRITE(34,*) ''
      DO n=1,itempsmax
	call getRandom (natoms,randx,randy,randz)
C		
	call actR(natoms,ar,br,cr,rx,ry,rz,vx,vy,vz,fx,fy,fz,randx,randy,randz)
	call actVmig(natoms,c0,bvmig,vx,vy,vz,fx,fy,fz,vxmig,vymig,vzmig)
C
	do k=1,npairs
	   i=pairs(1,k)
           j=pairs(2,k)
	   d(i,j) = sqrt((rx(i)-rx(j))**2+(ry(i)-ry(j))**2+(rz(i)-rz(j))**2)
c   c(i,j) = -const * (r0/d(i,j))**6 * (1-d0(i,j)/d(i,j))
	   c(i,j) = k0(i,j) * (1-d0(i,j)/d(i,j))
	   c(j,i) = c(i,j)
	end do
c
	call actForces (natomsmax,natoms,rx,ry,rz,c,fx,fy,fz,pairs,npairs)
C
	DO i=1,natoms
		vx(i)=vxmig(i)+bv*fx(i)+cv*randx(i)
		vy(i)=vymig(i)+bv*fy(i)+cv*randy(i)
		vz(i)=vzmig(i)+bv*fz(i)+cv*randz(i)
	ENDDO
C
	IF(mod(n,itsnap).eq.0)THEN
		m=n/itsnap
		write (34, 114) (rx(i),ry(i),rz(i),i=1,natoms) 
	ENDIF 
	ENDDO 
	CLOSE(34)
C
C JL retoco els formats per evitar errors en PDBs raros, faltaria nomatom i nom residu a A4
  111 fORMAT(i7.7) 
  112 FORMAT(a6,i5,2x,a2,2x,a3,2x,i4,4x,3f8.3,2f6.2)
  113 FORMAT(a6,i5,2x,a2,2x,a3,2x,i4,4x,3f8.3)
114   FORMAT(10F8.3) 
C
      END
C---------------------------------------------------------------
	subroutine readPDB (fn, rx,ry,rz,natoms)
	character*(*) fn
	real rx(*),ry(*),rz(*)
	integer natoms
	integer i,n1,n2
	character*6 atom
	character*2 ca
	character*4 res
c
	OPEN(unit=1,file=fn,status='UNKNOWN')
	i=1
	DO
	 READ(1,112,end=69)atom,n1,ca,res,n2,rx(i),ry(i),rz(i)
	IF(atom.eq.'ATOM  '.and.ca.eq.'CA') i=i+1
	ENDDO
   69 CLOSE(1)
	natoms = i-1
  112 FORMAT(a6,i5,2x,a2,2x,a3,2x,i4,4x,3f8.3,2f6.2)
	end
c--------------------------------------------------------
	subroutine moveToOrigin (rx,ry,rz,natoms)
	integer natoms
	real rx(natoms),ry(natoms),rz(natoms)
	real rcmx, rcmy, rcmz
C
      rcmx=0.d0
      rcmy=0.d0
      rcmz=0.d0
      DO i=1,natoms
        rcmx=rx(i)/natoms+rcmx
        rcmy=ry(i)/natoms+rcmy
        rcmz=rz(i)/natoms+rcmz
      ENDDO
      DO i=1,natoms
        rx(i)=rx(i)-rcmx
        ry(i)=ry(i)-rcmy
        rz(i)=rz(i)-rcmz
      ENDDO
	end
c------------------------------------------------------------
	subroutine calcMatDist (natomsmax, natoms,rx,ry,rz,d)
	integer natomsmax,natoms
	real rx(natoms),ry(natoms),rz(natoms)
	real d(natomsmax,natomsmax)
	integer i,j
	DO i=1,natoms
		d(i,i) = 0.0
	DO j=i+1,natoms
		d(i,j)=sqrt((rx(i)-rx(j))**2+(ry(i)-ry(j))**2+(rz(i)-rz(j))**2)
		d(j,i)=d(i,j)
	ENDDO
	ENDDO
	end
c---------------------------------------------------------------------
	subroutine calcPairs (natomsmax, natoms,d,CutOff,pairs,npairs)
	integer natomsmax,natoms
	real d(natomsmax, natomsmax), CutOff
	integer npairs, pairs(2,natomsmax*natomsmax)
	integer i,j
C
	npairs=0
	DO i=1,natoms
	DO j=i+1,natoms
		if (d(i,j).lt.CutOff) then
			npairs=npairs+1
			pairs(1,npairs)=i
			pairs(2,npairs)=j
		end if
	ENDDO
	ENDDO
	write (6,*) "Applied",CutOff," A cutoff, ",npairs," pairs"
	end
C---------------------------------------------------------
	subroutine getRandom (natoms,randx,randy,randz)
	integer natoms
	double precision randx(natoms),randy(natoms),randz(natoms)
        CALL dran_gv(randx,natoms)
        CALL dran_gv(randy,natoms)
        CALL dran_gv(randz,natoms)
	end
c-----------------------------------------------------------------------
	subroutine  actR(natoms,ar,br,cr,rx,ry,rz,vx,vy,vz,fx,fy,fz,randx,randy,randz)
	integer natoms,i
	real ar,br,cr
	real rx(natoms),ry(natoms),rz(natoms)
	real vx(natoms),vy(natoms),vz(natoms)
	real fx(natoms),fy(natoms),fz(natoms)
	real*8 randx(natoms),randy(natoms),randz(natoms)
c
	DO i=1,natoms
		rx(i)=rx(i)+ar*vx(i)+br*fx(i)+cr*randx(i)
		ry(i)=ry(i)+ar*vy(i)+br*fy(i)+cr*randy(i)
		rz(i)=rz(i)+ar*vz(i)+br*fz(i)+cr*randz(i)
	end do
	end
C--------------------------------------------------------------------------
	subroutine actVmig(natoms,c0,bvmig,vx,vy,vz,fx,fy,fz,vxmig,vymig,vzmig)
	integer natoms,i
	real c0,bvmig
	real vxmig(natoms),vymig(natoms),vzmig(natoms)
	real vx(natoms),vy(natoms),vz(natoms)
	real fx(natoms),fy(natoms),fz(natoms)
c
        DO i=1,natoms
          vxmig(i)=c0*vx(i)+bvmig*fx(i)
          vymig(i)=c0*vy(i)+bvmig*fy(i)
          vzmig(i)=c0*vz(i)+bvmig*fz(i)
	end do
	end
c-------------------------------------------------------------------
	subroutine actForces (natomsmax,natoms,rx,ry,rz,c,fx,fy,fz,pairs,npairs)
	integer natomsmax,natoms,i,j,npairs,k
	real rx(natoms),ry(natoms),rz(natoms)
	real fx(natoms),fy(natoms),fz(natoms)
	real c(natomsmax,natomsmax)
	integer pairs(2,natomsmax*natomsmax)
	do i=1,natoms
		fx(i)=0.
		fy(i)=0.
		fz(i)=0.
	end do
	do k=1,npairs
		i=pairs(1,k)
		j=pairs(2,k)
		fx(i) = c(i,j) * (rx(i) - rx(j)) + fx(i)
		fy(i) = c(i,j) * (ry(i) - ry(j)) + fy(i)
		fz(i) = c(i,j) * (rz(i) - rz(j)) + fz(i)
		fx(j) = c(j,i) * (rx(j) - rx(i)) + fx(j)
		fy(j) = c(j,i) * (ry(j) - ry(i)) + fy(j)
		fz(j) = c(j,i) * (rz(j) - rz(i)) + fz(j)
	ENDDO
	end
