cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PROGRAMA DMDGOTL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C DISTANCIES I VELOCITATS QUE SON EN ANGSTROMS
C TEMPERATURA EN K
C ENERGIES EN kcal/mol
C NOMES CONSIDERA ELS CA
C
C AL PRINCIPI MEMORITZA ELS TEMPS DE COLISIO DE TOTES LES
C PARELLES I DESPRES NOMES ACTUALITZA ELS TEMPS DE COLISIO
C DE LES PARTICULES QUE HAN XOCAT 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	program dmdgotl
	implicit real*8(a-h,o-z)
	parameter(natmax=1500)
	common / xoc / r(natmax,3),v(natmax,3),rbmax,rbmin,dmin,tij
       	dimension rant(natmax,3),rb(natmax,natmax),drb(natmax,natmax),r0(natmax,natmax),timp(natmax,natmax)
	dimension rhc(natmax),xm(natmax),ind1(natmax),ind2(natmax),ipart(natmax),tpart(natmax)
	dimension rcm(3),vcm(3),vd(3)
	character*100 file7,file8,file9,file10,file12,file19,file20,file21
 	character*5 c1,c2,c3,atom(natmax),res(natmax)
	dimension tmpvd(natmax*3)
	namelist /input/ tsnap,rc,rca,rn,rcb,sigma,file9,temp,kkk,file7,file8,file9,file10,file19,file20,file21,rhf,ehf,ires1,ires2,nbloc,id,iread,rcut,ixerra,rcutgo,nevent,sigmago,xmres,file12,xmarge,rchf,sigmahf,isec,trnd
 1000	FORMAT (A4,3X,I4,2X,A2,2X,A3,3X,I3,5X,F7.3,1X,F7.3,1X,F7.3)
 1001	FORMAT (10F8.3)
 2000	FORMAT (A4,3X,I4,2X,A2,2X,A3,3X,I3,5X,F7.3,1X,F7.3,1X,F7.3,1X,F7.5)
 333	FORMAT('MODEL',8X,I4)
 334	FORMAT('ENDMDL')
c default values
	file7='topologia.dat'
	file8='molecula.pdb'
	file10='distancia.dat'
	file9='molecula.out'
	file20='snapshots.pdb'
	file12='masses.dat'
	file19='output.pdb'
	file21='snapshots.crd'
	xmarge=1.e-3
	a=1.d-10
	tmin=1.d-30
	kkk=2839
	temp=300.
	trnd=1.d-10
	nbloc=2
	isec=0
	rcut=20.
	rc=1.5
	rca=0.5
	rcb=2.2
	rn=1.3
	xmres=.1
	sigma=.05
	sigmago=0.1
	sigmahf=0.1
	rchf=0.1
	rcutgo=8.
	ehf=0.
	id=0
	ires1=1
	ires2=2
	iread=0
	ixerra=0
	read(5,INPUT)
c converteixo la temperatura en energia per particula
c utilitzo la constant de Boltzmann per mol
	temp=1.5*8.32*temp
c llegeix les coordenades (fitxer pdb)
	open(unit=9,file=file9)
	n=1
 50	read(9,1000,end=51)c1,j,c2,c3,k,x,y,z
	if(c2.eq.'CA')then
	atom(n)=c2
	ind2(n)=k
	ind1(n)=n
	res(n)=c3
	r(n,1)=x
	r(n,2)=y
	r(n,3)=z
	n=n+1
	endif
	goto 50
 51	natom=n-1
	close(9)

ccccccccccccccccccccccccccccccccccccccc
c llegeix la matriu de topologia
ccccccccccccccccccccccccccccccccccccccc
	do i=1,natom-1
	do j=i+1,natom
	rb(i,j)=0.
	drb(i,j)=0.
	enddo
	enddo


	ngo=0
	jmax=natom
c prepara la matriu de topologia rb
	do i=1,natom-1
c primer el seguent en la sequencia
	  rij1=r(i,1)-r(i+1,1)
	  rij2=r(i,2)-r(i+1,2)
	  rij3=r(i,3)-r(i+1,3)
	  dist=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(dist.lt.rcutgo)then
	 rb(i,i+1)=dist
	if(dist.lt.4.5)then
	 drb(i,i+1)=sigma*dist
	else
	 drb(i,i+1)=sigmago*dist
	endif
	endif
c i ara tota la resta
	if(isec.eq.1)jmax=i+4
	 do j=i+2,jmax
	  rij1=r(i,1)-r(j,1)
	  rij2=r(i,2)-r(j,2)
	  rij3=r(i,3)-r(j,3)
	  dist=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	 if(dist.lt.rcutgo)then
	 rb(i,j)=dist
	 drb(i,j)=sigmago*dist
	 ngo=ngo+1
	 endif
	 enddo
	enddo


	write(6,*)'natom,tsnap,nbloc',natom,tsnap,nbloc

c assigna velocitats aleatories a les particules
c assigna radis de hardcore en funcio de l'element
	open(unit=12,file=file12)
	xmassa=0.
	do i=1,natom
	xm(i)=xmres
	rhc(i)=rca
	rewind(12)
  130	read(12,*,end=131)c1,xmm
	if(c1.ne.res(i))goto 130
	xm(i)=1.d-3*xmm
  131	xmassa=xmassa+xm(i)
	do j=1,3
	call ran1(fi,kkk)
	v(i,j)=fi
	enddo
	enddo
	write(6,*)'massa',xmassa
	close(12)

c ajusta l'energia cinetica a la temperatura requerida
	ekin=0.
	do i=1,natom
	do j=1,3
	ekin=ekin+0.5*xm(i)*(v(i,j)*a)**2
	enddo
	enddo
	sto=sqrt(natom*temp/ekin)
	ekin0=0
	do i=1,natom
	do j=1,3
	v(i,j)=v(i,j)*sto
	ekin0=ekin0+0.5*xm(i)*(v(i,j)*a)**2
	enddo
	do j=i+1,natom
	timp(i,j)=1.d6
	enddo
	enddo
	ekin=ekin0


	do i=1,natom-2
	do j=i+2,natom
	  rij1=r(i,1)-r(j,1)
	  rij2=r(i,2)-r(j,2)
	  rij3=r(i,3)-r(j,3)
	  rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	rmin=rhc(j)+rhc(i)
	if(rij.lt.rmin)then
	write(6,*)'ATOMS MASSA APROP',i,j,rij,rmin
	stop
	endif
	enddo
	enddo

	write(6,*)'Contactes Go',ngo

	c1='ATOM'
	
c ara busca el CM
	do j=1,3
	rcm(j)=0.
c	vcm(j)=0.
	do i=1,natom
	rcm(j)=rcm(j)+xm(i)*r(i,j)
c	vcm(j)=vcm(j)+xm(i)*v(i,j)
	enddo
        rcm(j)=rcm(j)/xmassa
c     vcm(j)=vcm(j)/xmassa
	enddo

c escriu les coordenades en el SRCM
	open(unit=20,file=file20)
	write(20,333)ibloc
c JL
	open(unit=21,file=file21)
	write(21,*) ''
c JL
	open(unit=8,file=file8)
	do i=1,natom
	do j=1,3
	rant(i,j)=r(i,j)
	vd(j)=r(i,j)-rcm(j)
c JL
	tmpvd((i-1)*3+j)=vd(j)
c JL
	enddo
	write(8,2000)c1,ind1(i),atom(i),res(i),ind2(i),(vd(j),j=1,3),rhc(i)
	write(20,1000)c1,ind1(i),atom(i),res(i),ind2(i),(vd(j),j=1,3)
	enddo
c JL
	write(21,1001) (tmpvd(j),j=1,natom*3)
c JL
	close(8)
	write(20,334)

c aqui llegeix les coordenades i velocitats del calcul anterior
	if(iread.eq.1)then
	open(unit=9,file='molecula.out')
	do i=1,natom
	read(9,*)r(i,1),r(i,2),r(i,3),v(i,1),v(i,2),v(i,3)
	enddo
	close(9)
	endif

c comença a iterar per buscar l'event mes proper
c la variable ixoc indica quin tipus d'event ocorre:
c 0 -> enllaç
c 1 -> xoc entre atoms no enllaçats

	open(unit=10,file=file10)
	nback=0
	nhf=0
	iat=0
	temps=0.
	ierr=0

	tevent=0.
	temps=0

c al final de cada cicle escriu una foto

	do 200 ibloc=1,nbloc

        iev=0
	nback=0
	iback=0
        tacum=0.
	ierr=0

 20	do i=1,natom
	do j=1,3
	r(i,j)=rant(i,j)
	enddo
	enddo

	do i=1,natom-1
	do j=i+1,natom
        rij1=r(i,1)-r(j,1)
        rij2=r(i,2)-r(j,2)
        rij3=r(i,3)-r(j,3)
        r0(i,j)=sqrt(rij1**2+rij2**2+rij3**2)
	enddo
	enddo

c	write(6,*)'tevent0',tevent0

	ixoc=1

cccccccccccccccccccccccccccccccccccccccccccccccccc
c rastreja totes les parelles, buscant l'event mes proper
cccccccccccccccccccccccccccccccccccccccccccccccccc

       do i=1,natom-1
       do j=i+1,natom
       tij1=1.d6
       tij2=1.d6
c particules enllaçades
	if (rb(i,j).gt.1.d-20)then
	rbmin=rb(i,j)-drb(i,j)
	rbmax=rb(i,j)+drb(i,j)
	call colbond(i,j)
	  if(tij.gt.tmin)timp(i,j)=tij+temps


	else

c xoc frontal entre particules no enllaçades
	dmin=rhc(i)+rhc(j)
	call colnobond(i,j)
	  if(tij.gt.tmin)timp(i,j)=tij+temps
	endif


       enddo
       enddo

c	write(6,*)'llista de temps actualitzada',iev,ierr

	do i=1,natom-1
	tpart(i)=1.
	do j=i+1,natom
	if(timp(i,j).lt.temps)then
	write(6,*)'zero',i,j,timp(i,j),temps
	stop
	endif
	if(timp(i,j).lt.tpart(i))then
	tpart(i)=timp(i,j)
	ipart(i)=j
	endif
	enddo
	enddo

c aqui comença la evolucio temporal

 10	tevent0=tevent
c	if(timp(1,2).gt.1.d-12)then
c	write(6,*)ibloc,iev,temps,timp(1,2),mem1,mem2
c	stop
c	endif

c busca quina es la primera colisio
	tevent=1.
	do i=1,natom-1
	j=ipart(i)
	if(tpart(i).lt.tevent)then
	tevent=tpart(i)
	mem1=i
	mem2=j
	endif
	enddo
	
c	write(6,*)'natom',natom,mem1,mem2

	tevent1=tevent-temps

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c translacio i variacio dels temps
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       do i=1,natom
	do j=1,3
	rant(i,j)=r(i,j)
        r(i,j)=rant(i,j)+v(i,j)*tevent1
	enddo
       enddo

	if(rb(mem1,mem2).lt.1.d-20)nback=nback+1

       tacum=tacum+tevent1
	tacumrnd=tacumrnd+tevent1
	temps=temps+tevent1

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c canvi de velocitats
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	vdmod=0
	rij=0
	do j=1,3
	vdmod=vdmod+(v(mem2,j)-v(mem1,j))*(r(mem2,j)-r(mem1,j))
	rij=rij+(r(mem2,j)-r(mem1,j))**2
	enddo
	rij=sqrt(rij)
	vdmod=vdmod/rij
	xsum=0.5*(1./xm(mem1)+1./xm(mem2))
	if(ixoc.le.1)then
c modul del moment transferit en la colisio
	dp=-vdmod/xsum
	else
	sto=vdmod**2+4.*deltak*xsum
 	 if(sto.gt.0)then
c vario la velocitat
  	 dp=-vdmod+sqrt(sto)
	 dp=dp/(2.*xsum)
	 else
c la particula rebota (no supera la barrera)
	 ixoc=4
	 iat=iat+1
	 dp=-vdmod/xsum
	 endif
	endif
	do j=1,3
        v(mem1,j)=v(mem1,j)-dp/xm(mem1)*(r(mem2,j)-r(mem1,j))/rij
        v(mem2,j)=v(mem2,j)+dp/xm(mem2)*(r(mem2,j)-r(mem1,j))/rij
	enddo

       iev=iev+1

c casos en els quals ha canviat l'energia cinetica
	if(ixoc.eq.2.or.ixoc.eq.3)then

c ajusta l'energia cinetica a la temperatura requerida
	ekin=0.
	do i=1,natom
	do j=1,3
	ekin=ekin+0.5*xm(i)*(v(i,j)*a)**2
	enddo
	enddo
	sto=sqrt(natom*temp/ekin)
c	ekin0=0
	do i=1,natom
	do j=1,3
	v(i,j)=v(i,j)*sto
c	ekin0=ekin0+0.5*xm(i)*v(i,j)**2
	enddo
	enddo

	if(ixerra.gt.0)then
	write(6,*)'Event,ixoc,particules',iev,ixoc,mem1,mem2
	endif

	endif

	if(id.gt.0)then
	sto=0
	do j=1,3
	sto=sto+(r(ires1,j)-r(ires2,j))**2
	enddo
	sto=sqrt(sto)
	rbmin=(1.-xmarge)*(rb(ires1,ires2)-drb(ires1,ires2))
	rbmax=(1.+xmarge)*(rb(ires1,ires2)+drb(ires1,ires2))
	write(10,*)temps,sto,rbmin,rbmax,mem1,mem2,tevent1
	endif

	timp(mem1,mem2)=1.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c ara calcula els temps de colisio per a les dues particules que han xocat
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       do i=1,mem1-1
       j=mem1
c particules enllaçades
	if (rb(i,j).gt.1.d-20)then
	rbmin=rb(i,j)-drb(i,j)
	rbmax=rb(i,j)+drb(i,j)
	call colbond(i,j)
	if(tij.gt.tmin)timp(i,j)=tij+temps

	else

c xoc frontal entre particules no enllaçades
	dmin=rhc(i)+rhc(j)
	if(r0(i,j).gt.dmin.and.r0(i,j).lt.rcut)then
	call colnobond(i,j)
	 if(tij.gt.tmin)timp(i,j)=tij+temps
	endif

	endif

       enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       do j=mem1+1,natom
       i=mem1
	if(j.ne.mem2)then
c particules enllaçades
	if (rb(i,j).gt.1.d-20)then
	rbmin=rb(i,j)-drb(i,j)
	rbmax=rb(i,j)+drb(i,j)
	call colbond(i,j)
	if(tij.gt.tmin)timp(i,j)=tij+temps

	else

c xoc frontal entre particules no enllaçades
	dmin=rhc(i)+rhc(j)
	if(r0(i,j).gt.dmin.and.r0(i,j).lt.rcut)then
	call colnobond(i,j)
	 if(tij.gt.tmin)timp(i,j)=tij+temps
	endif

	endif
	endif

       enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       do i=1,mem2-1
	if(i.ne.mem1)then
        j=mem2
c particules enllaçades
	if (rb(i,j).gt.1.d-20)then
	rbmin=rb(i,j)-drb(i,j)
	rbmax=rb(i,j)+drb(i,j)
	call colbond(i,j)
	if(tij.gt.tmin)timp(i,j)=tij+temps

	else
c xoc frontal entre particules no enllaçades
	dmin=rhc(i)+rhc(j)
	if(r0(i,j).gt.dmin.and.r0(i,j).lt.rcut)then
	call colnobond(i,j)
	 if(tij.gt.tmin)timp(i,j)=tij+temps
	endif

	endif
	endif

       enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       do j=mem2+1,natom
       i=mem1
c particules enllaçades
	if (rb(i,j).gt.1.d-20)then
	rbmin=rb(i,j)-drb(i,j)
	rbmax=rb(i,j)+drb(i,j)
	call colbond(i,j)
	if(tij.gt.tmin)timp(i,j)=tij+temps

	else

c xoc frontal entre particules no enllaçades
	dmin=rhc(i)+rhc(j)
	if(r0(i,j).gt.dmin.and.r0(i,j).lt.rcut)then
	call colnobond(i,j)
	 if(tij.gt.tmin)timp(i,j)=tij+temps
	endif

	endif

       enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	do i=1,natom-1
	if(tpart(i).lt.tmin)then
	write(6,*)'cagada',i,mem1,mem2
	stop
	endif
	enddo

c actualitza la llista de temps
c la que ha xocat amb les de despres
	i=mem1
	tpart(i)=1.
	do j=mem1+1,natom
	if(timp(i,j).lt.tpart(i))then
	ipart(i)=j
	tpart(i)=timp(i,j)
	endif
	enddo
c amb les d'abans
	do j=1,mem1-1
	if(timp(j,i).lt.temps)then
	write(6,*)'zero',j,i,timp(j,i),mem1,mem2
	stop
	endif
	if(timp(j,i).lt.tpart(j))then
	ipart(j)=i
	tpart(j)=timp(j,i)
	endif
	enddo
c la que ha xocat amb les de despres
	i=mem2
	tpart(i)=1.
	do j=mem2+1,natom
	if(timp(i,j).lt.tpart(i))then
	ipart(i)=j
	tpart(i)=timp(i,j)
	endif
	enddo
c amb les d'abans
	do j=1,mem2-1
	if(timp(j,i).lt.temps)then
	write(6,*)'zero',j,i,timp(j,i),mem1,mem2
	stop
	endif
	if(timp(j,i).lt.tpart(j))then
	ipart(j)=i
	tpart(j)=timp(j,i)
	endif
	enddo

	do i=1,natom-1
	if(tpart(i).lt.tmin)then
	write(6,*)'zero',i,mem1,mem2,temps
	stop
	endif
	enddo
	

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	if(tacumrnd.gt.trnd)then

	tacumrnd=0
	do i=1,natom
	do j=1,3
	rant(i,j)=0.5*(r(i,j)+rant(i,j))
	call ran1(fi,kkk)
	v(i,j)=fi
	enddo
	enddo
c ajusta l'energia cinetica a la temperatura requerida
	ekin=0.
	do i=1,natom
	do j=1,3
	ekin=ekin+0.5*xm(i)*(v(i,j)*a)**2
	enddo
	enddo
	sto=sqrt(natom*temp/ekin)
	ekin0=0
	do i=1,natom
	do j=1,3
	v(i,j)=v(i,j)*sto
	ekin0=ekin0+0.5*xm(i)*(v(i,j)*a)**2
	enddo
	do j=i+1,natom
	timp(i,j)=1.d6
	enddo
	enddo
	ekin=ekin0
	goto 20

	endif


       if(tacum.lt.tsnap) goto 10

	sto0=ekin0/natom*2.389d-4
	sto=ekin/natom*2.389d-4

	do i=1,natom
	do j=1,3
	rant(i,j)=0.5*(r(i,j)+rant(i,j))
	call ran1(fi,kkk)
	v(i,j)=fi
	enddo
	enddo
c ajusta l'energia cinetica a la temperatura requerida
	ekin=0.
	do i=1,natom
	do j=1,3
	ekin=ekin+0.5*xm(i)*(v(i,j)*a)**2
	enddo
	enddo
	sto=sqrt(natom*temp/ekin)
	ekin0=0
	do i=1,natom
	do j=1,3
	v(i,j)=v(i,j)*sto
	ekin0=ekin0+0.5*xm(i)*(v(i,j)*a)**2
	enddo
	do j=i+1,natom
	timp(i,j)=1.d6
	enddo
	enddo
	ekin=ekin0

c calula la posicio del CM
	do j=1,3
	rcm(j)=0.
	do i=1,natom
	rcm(j)=rcm(j)+xm(i)*rant(i,j)
	enddo
        rcm(j)=rcm(j)/xmassa
	enddo

	write(20,333)ibloc
	do i=1,natom
	do j=1,3
	vd(j)=rant(i,j)-rcm(j)
c JL
	tmpvd((i-1)*3+j)=vd(j)
c JL
	enddo
	write(20,1000)c1,ind1(i),atom(i),res(i),ind2(i),(vd(j),j=1,3)
	enddo
c JL
	write(21,1001)(tmpvd(j),j=1,natom*3)
c JL
	inp=0
	do i=1,natom-1
	do j=i+1,natom
	dmin=rhc(i)+rhc(j)
       rij1=r(i,1)-r(j,1)
       rij2=r(i,2)-r(j,2)
       rij3=r(i,3)-r(j,3)
       rij=sqrt(rij1**2+rij2**2+rij3**2)
	if(rb(i,j).gt.1.d-10)then
	rbmin=rb(i,j)-drb(i,j)
	rbmax=rb(i,j)+drb(i,j)
	if(rij.gt.rbmax.or.rij.lt.rbmin)inp=inp+1
	else
	if(rij.lt.dmin)inp=inp+1
	endif
	enddo
	enddo
c
	write(20,334)
c	write(20,*)'************',sto,iev

	write(6,*)'Temps',temps,' Events',iev,' backbone',nback,ierr,iback,' no permesos',inp
c	write(6,*)'Bloc',ibloc,' Energia cinetica inicial i final',sto0,sto,' Potencial',epot
 200	continue

	close(10)
	close(20)
	close(21)
c calula la posicio del CM
	do j=1,3
	rcm(j)=0.
	do i=1,natom
	rcm(j)=rcm(j)+xm(i)*r(i,j)
	enddo
        rcm(j)=rcm(j)/xmassa
	enddo

c escriu posicions i velocitats
	open(unit=9,file='molecula.out')
	do i=1,natom
	do j=1,3
	vd(j)=r(i,j)-rcm(j)
	enddo
	write(9,*)vd(1),vd(2),vd(3),v(i,1),v(i,2),v(i,3)
	enddo
	close(9)

c escriu les coordenades finals en el SRCM
        open(unit=19,file=file19)
	do i=1,natom
	do j=1,3
	vd(j)=r(i,j)-rcm(j)
	enddo
	write(19,1000)c1,ind1(i),atom(i),res(i),ind2(i),(vd(j),j=1,3)
	enddo
	close(19)
c
	do i=1,natom
	do j=i+1,natom	
	rij1=r(i,1)-r(j,1)
        rij2=r(i,2)-r(j,2)
        rij3=r(i,3)-r(j,3)
        rij=sqrt(rij1**2+rij2**2+rij3**2)
	if(atom(i).eq.'CB'.and.atom(j).eq.'CB')then
	 if(rij.lt.rhf)write(6,*)atom(i),i,j,rij,'   DINTRE'
	endif
	enddo
	enddo

	niter=iev*nbloc
        write(6,*)'Temps,events,nback,nhf,iat',tacum,niter,nback,nhf,iat
	write(6,*)'Energies inicials',sto0,epot0
	write(6,*)'Energies finals',sto,epot
	write(6,*)'Distancia de tall de la interaccio hidrofobica',rhf
       stop
       end




ccccc
CC  
CC  rutina ran1
        subroutine ran1(z,idum)
        implicit real*8 (a-h,o-z)
c
        SAVE
c
        dimension r(97)
        parameter ( m1=259200, ia1=7141, ic1=54773, rm1=1./m1)
        parameter ( m2=134456, ia2=8121, ic2=28411, rm2=1./m2)
        parameter ( m3=243000, ia3=4561, ic3=51349)
c
        data iff / 0 /
c
        if(idum.lt.0 .or. iff.eq.0) then
          iff = 1
          ix1 = mod ( ic1 - idum , m1 )
          ix1 = mod ( ia1 * ix1 + ic1 , m1 )
          ix2 = mod ( ix1, m2 )
          ix1 = mod ( ia1 * ix1 + ic1 , m1 )
          ix3 = mod ( ix1, m3 )
          do j = 1, 97
             ix1 = mod ( ia1 * ix1 + ic1 , m1 )
             ix2 = mod ( ia2 * ix2 + ic2 , m2 )
             r(j) = ( float(ix1) + float(ix2) * rm2 ) * rm1
          end do
        end if
c
c       ccccccccccccccccccccccc   except when initializing, this is 
c                                 were we start.
        ix1 = mod ( ia1 * ix1 + ic1 , m1 )
        ix2 = mod ( ia2 * ix2 + ic2 , m2 )
        ix3 = mod ( ia3 * ix3 + ic3 , m3 )
        j = 1 + ( 97 * ix3 ) / m3
c
c
        if ( j.gt.97 .or. j.lt.1 ) write(6,*) ' AAAUUGGHH !!!'
c
        z = r(j)
        r(j) = ( float(ix1) + float(ix2) * rm2 ) * rm1
c------------------------------------------------------
c       do j = 1 , 97
c       write(6,*) ' j, r(j) = ', j, r(j)
c       end do
c------------------------------------------------------
c       write(6,*) 
c
        return
        end
CC  
	subroutine colbond(i,j)
	implicit real*8(a-h,o-z)
	parameter(natmax=1500)
	common / xoc / r(natmax,3),v(natmax,3),rbmax,rbmin,dmin,tij
        rij1=r(i,1)-r(j,1)
        rij2=r(i,2)-r(j,2)
        rij3=r(i,3)-r(j,3)
        rij=sqrt(rij1**2+rij2**2+rij3**2)
        vij1=v(i,1)-v(j,1)
        vij2=v(i,2)-v(j,2)
        vij3=v(i,3)-v(j,3)
        vij=sqrt(vij1**2+vij2**2+vij3**2)
        cosang=(rij1*vij1+rij2*vij2+rij3*vij3)/(rij*vij)
	if(rij.gt.rbmax.and.cosang.gt.0)then
	v(i,1)=v(i,1)-vij1
	v(i,2)=v(i,2)-vij2
	v(i,3)=v(i,3)-vij3
	v(j,1)=v(j,1)+vij1
	v(j,2)=v(j,2)+vij2
	v(j,3)=v(j,3)+vij3
c	ierr=ierr+1
	endif
	if(rij.lt.rbmin.and.cosang.lt.0)then
	v(i,1)=v(i,1)-vij1
	v(i,2)=v(i,2)-vij2
	v(i,3)=v(i,3)-vij3
	v(j,1)=v(j,1)+vij1
	v(j,2)=v(j,2)+vij2
	v(j,3)=v(j,3)+vij3
c	ierr=ierr+1
	endif
c si cosang < 0 vol dir que s'acosten
       sinang=sqrt(1.-cosang**2)
c distancia en el moment de la colisio
       rcol=sinang*rij
c xoc frontal
c	argf=(rb(i,j)*(1.-sigma))**2-rcol**2
	argf=rbmin**2-rcol**2
c rebot
c	argb=(rb(i,j)*(1.+sigma))**2-rcol**2
	argb=rbmax**2-rcol**2
	tij=0
         if(cosang.lt.0)then
c s'acosten: xoc frontal
          if(argf.gt.0)tij=(-rij*cosang-sqrt(argf))/vij
         else
c s'estan allunyant: rebot
	  tij=(-rij*cosang+sqrt(argb))/vij
         endif

	return
	end

	subroutine colnobond(i,j)
	implicit real*8(a-h,o-z)
	parameter(natmax=1500)
	common / xoc / r(natmax,3),v(natmax,3),rbmax,rbmin,dmin,tij
        rij1=r(i,1)-r(j,1)
        rij2=r(i,2)-r(j,2)
        rij3=r(i,3)-r(j,3)
        rij=sqrt(rij1**2+rij2**2+rij3**2)
        vij1=v(i,1)-v(j,1)
        vij2=v(i,2)-v(j,2)
        vij3=v(i,3)-v(j,3)
        vij=sqrt(vij1**2+vij2**2+vij3**2)
        cosang=(rij1*vij1+rij2*vij2+rij3*vij3)/(rij*vij)
	if(rij.lt.dmin.and.cosang.lt.0)then
	v(i,1)=v(i,1)-vij1
	v(i,2)=v(i,2)-vij2
	v(i,3)=v(i,3)-vij3
	v(j,1)=v(j,1)+vij1
	v(j,2)=v(j,2)+vij2
	v(j,3)=v(j,3)+vij3
c	ierr=ierr+1
	endif
c si cosang < 0 vol dir que s'acosten
       sinang=sqrt(1.-cosang**2)
c distancia en el moment de la colisio
       rcol=sinang*rij
	tij=0
         if(cosang.lt.0)then
	  arg=rbmin**2-rcol**2
c s'acosten: xoc frontal
          if(arg.gt.0)tij=(-rij*cosang-sqrt(arg))/vij
         endif

	return
	end








