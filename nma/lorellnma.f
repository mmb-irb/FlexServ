!*****************************************************************************************
! This program calculates the matrix of harmonic constants for a protein, according to sequence connectivity and size-dependent cutoff rules. The values of all parameters have been optimized by comparison with essential dynamics. 
! Input data needed: Protein Lenght (N, provided by a "data" file containing pdb code and protein size) and Ca coordinates (here provided by Nx3 matrix "coordinates" file) 
! Output data: "constants.dat" file containing a set of constants in format I,J,K(I,J)
! Author: Laura O. R., October 2008
! This code is written in Fortran95 language
!*****************************************************************************************

      PROGRAM NMA

      IMPLICIT NONE

      INTEGER :: N,I,J,P
      REAL :: CUTOFF
      REAL, PARAMETER :: SL=3, IN=2.8  
      REAL*4 Cont,Seq,S,X,Y,Z,D,K
      REAL*4 , allocatable  :: C(:,:)
      REAL, PARAMETER :: Back= 0       
      REAL, PARAMETER :: Rc=6, Rs=60
      INTEGER, PARAMETER :: Ps=3, E1=2, E2=6   
      CHARACTER (len=255) :: arg, coords, fcte
    !JL
      character (len=1), allocatable :: chain(:)
      integer, allocatable :: numres(:)
      character (len=4) :: dumres

    ! Set order2 square matrices dimensions NxN, which are given by protein sequence lenght...

    ! Retrieve the number of atoms and coordinate and output file name
    ! from the command line parameters
      CALL GETARG (1, arg)
      CALL GETARG (2, coords)
      CALL GETARG (3, fcte)
      READ (arg, '(I8)') N 

      ALLOCATE (C(N,N))
      ALLOCATE (chain(N))
      ALLOCATE (numres(N))
      C=0.0

    ! Set automatically cutoff for calculations according to protein lenght, N 
     
      CUTOFF=INT(SL*LOG(REAL(N))-IN)

   ! The program needs as input the Ca coordinates parsed from PDB files...
   ! JL Inclou lectura de cadena i num de res

      OPEN (UNIT=1,FILE=coords,STATUS='OLD')
      DO I=1,N 
        READ (1,'(a4,a1,i4,2x,3F8.3)') dumres, chain(i),numres(i),C(I,1),C(I,2),C(I,3)           
      ENDDO
      CLOSE (UNIT=1) 


    ! Calculate Ca-Ca distances and fill the Ca-Ca distances upper matrix, then calculate the corresponding constants matrix...
      OPEN(UNIT=35,FILE=fcte,STATUS='UNKNOWN') 
      DO I=1,N
       DO J=I+1,N
         X=(C(I,1)-C(J,1))**2
         Y=(C(I,2)-C(J,2))**2
         Z=(C(I,3)-C(J,3))**2
         D=SQRT(X+Y+Z)         

         S=J-I
         Seq=Rs/(S**E1)
         Cont=(Rc/D)**E2

   	IF(S.LE.Ps) THEN
           K=Seq                 
         ELSE
           IF (D.LE.CUTOFF) THEN
             K=Cont
           ELSE
             K=Back
           ENDIF
         ENDIF   

         WRITE (35,*) I,J,K
       ENDDO
      ENDDO 
      CLOSE (UNIT=35) 

      DEALLOCATE (C)

      WRITE (6,*) ' Work done !! '

      END

