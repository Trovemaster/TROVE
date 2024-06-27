!
! This unit is for a user defined potential: N2O PES 
!
module pot_user
  use accuracy
  use moltype
  use pot_xy2, only : MLloc2pqr_xyz

  implicit none

  public MLdipole,MLpoten,ML_MEP,MLpoten_name

  private
 
  integer(ik), parameter :: verbose     = 4                          ! Verbosity level
  !
  integer(ik) :: npol(5),order(5,5)
  real(ark)   :: gam1(5,5),x0(5,5),r0(5,5),beta(5,5)
  integer(ik) :: np(5)
  !
  contains
  !
  function ML_MEP(dim,rho)  result(f)

   integer(ik),intent(in) ::  dim
   real(ark),intent(in)   ::  rho
   real(ark)              ::  f(dim)
   !
   if (dim/=3) stop 'Illegal size of the function - must be 3'
   !
   f(:) = molec%local_eq(:)
   f(molec%Ncoords) = rho

  end function ML_MEP

 ! Check the potential name 
 subroutine MLpoten_name(name)
   !
   character(len=cl),intent(in) ::  name
   character(len=cl),parameter ::  poten_name = 'GENERAL'
   ! 
   if (poten_name/=trim(name)) then
     write(out,"(a,a,a,a)") 'Wrong Potential ',trim(name),'; should be ',trim(poten_name)
   endif
   !
   write(out,"(a,a)") '  Using USER-type PES ',trim(poten_name)
   !
 end subroutine MLpoten_name
 !


 !
 ! Defining potential energy function (built for SO2)

 function MLpoten(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force( :)
   real(ark)              ::  f,f1,f2
   integer(ik)            ::  nparam,ntot
   !
   !
   nparam = int(force(1))
   ntot = size(force)
   !
   f = MLpoten_triatomic(ncoords,natoms,local,xyz,force)
   !
 end function MLpoten
 !


 recursive subroutine MLdipole(rank,ncoords,natoms,local,xyz,f)
   !
   integer(ik),intent(in) ::  rank,ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
   real(ark),intent(out)  ::  f(rank)
    !
    real(ark) :: xyz_(natoms,3), r1, r2,  n1(3), n2(3), n3(3), tmat(3,3), CN_CM(3),&
                 x(natoms,3),xyz0(natoms,3),cos_theta,alpha,mu(3),u1(3),u2(3),u3(3),bigr,smallr

    !
    call MLdmspq_triat(rank,ncoords,natoms,local,xyz,f)
    !
  end subroutine MLdipole





  function MLpoten_triatomic(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   real(ark) :: rr2,rr3,alpha,xcos,v

   if (verbose>=6) write(out,"('MLpoten_triatomic/start')")

     rr2 = local(1); rr3 = local(2) ; alpha = local(3)
     xcos = cos(alpha)

     CALL potvv(v,rr2,rr3,xcos,force)

         f = v

    if (verbose>=6) write(out,"('MLpoten_triatomic/end')")

 end function MLpoten_triatomic


!######################################################################################################
      SUBROUTINE potvv(V_ark,r2_ark,r3_ark,xcos_ark,force)

      IMPLICIT none 

      real(ark),intent(in) :: r2_ark,r3_ark,xcos_ark
      real(ark),intent(out) :: V_ark
      real(ark) :: r2,r3,xcos,v
      real(ark),intent(in)   ::  force(:)
      integer(ik) :: NC,N

      r3=r3_ark
      r2=r2_ark
      xcos=xcos_ark

      N=size(force)
      NC=N-11

         CALL TRIAT(r2,r3,xcos,v,force,NC)

         v_ark = v

      end SUBROUTINE potvv


!######################################################################################################
      SUBROUTINE TRIAT(rr1,rr2,rr3,POT,force,NC)
      IMPLICIT none 
      real(ark),intent(in)   :: rr1,rr2,rr3
      real(ark),intent(in)   ::  force(:)
      real(ark),intent(out) :: POT
      CHARACTER(LEN=3) :: MOLTYP
      integer(ik), PARAMETER :: NX=3
      integer(ik) :: I,J,K,L,M,ID
      integer(ik) :: POLORDER,MTYPE,NC
      real(ark), DIMENSION(NC)   :: C
      real(ark), ALLOCATABLE, DIMENSION(:) :: POL
      real(ark), DIMENSION(NX) :: Y
      integer(ik) :: TOTNUM,MNUM
      integer(ik), DIMENSION(3,6) :: P
      real(ark) :: A1,A2,R1EQ,R2EQ,ALPHAEQ
      real(ark) :: PI,DEGS,B1,B2,G1,G2,R1

      ALLOCATE(POL(NC))

      MTYPE = force(1)
      POLORDER = force(2)
      A1 = force(3)
      A2 = force(4)
      R1EQ = force(5)
      R2EQ = force(6)
      ALPHAEQ = force(7)
      B1 = force(8)
      B2 = force(9)
      G1 = force(10)
      G2 = force(11) 

      IF (MTYPE .EQ. 1) THEN
         MOLTYP = 'ABC'
      ELSE IF (MTYPE .EQ. 2) THEN
         MOLTYP = 'AB2'
      ELSE IF (MTYPE .EQ. 3) THEN
         MOLTYP = 'A3'
      ELSE

         STOP 'Illegal moltype of PES'

      END IF


      DO I=1,NC
        C(I)=force(I+11)
      END DO     

      R1=SQRT(rr1**2+rr2**2-2*rr1*rr2*rr3)

      PI=3.11415926535D+00
      DEGS=180.0D+00

      DO I=1,NC
        POL(I)=0.00D+00
      END DO

      DO I=1,NX
        Y(I)=0.0D+00
      END DO

      Y(2)=1-EXP(-A1*(rr1-R1EQ))
      Y(3)=1-EXP(-A2*(rr2-R2EQ))
      Y(1)=rr3-COS(ALPHAEQ*PI/DEGS)

      TOTNUM=0
      DO M=0,POLORDER
        MNUM=0
        DO I=0,M
          DO J=0,(M-I)
            K=M-I-J
!

              IF (MOLTYP.EQ.'ABC') THEN  

                TOTNUM=TOTNUM+1
                MNUM=MNUM+1

                CALL PERMUTABC(I,J,K,P,ID)

                DO L=1,ID

                POL(TOTNUM)=POL(TOTNUM)+ (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))
  
                END DO

                POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)

              ELSE IF (MOLTYP.EQ.'AB2') THEN 

                IF (J.LE.K) THEN 

                  TOTNUM=TOTNUM+1
                  MNUM=MNUM+1

                  CALL PERMUTAB2(I,J,K,P,ID)
                   
                  DO L=1,ID

                    POL(TOTNUM)=POL(TOTNUM)+ (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))

                  END DO

                  POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)

                END IF    
         
              ELSE IF (MOLTYP.EQ.'A3') THEN

                IF (I.LE.J .AND. J.LE.K) THEN

                  TOTNUM=TOTNUM+1
                  MNUM=MNUM+1

                  CALL PERMUTA3(I,J,K,P,ID)
                   
                  DO L=1,ID

                   POL(TOTNUM)=POL(TOTNUM)+ (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))

                  END DO

                  POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)

                END IF

              END IF

          END DO
        END DO

      END DO

         POT=SUM(POL)+B1*EXP(-G1*R1)+B2*EXP(-G2*R1**2)


      DEALLOCATE(POL)

      RETURN
      END SUBROUTINE TRIAT
      
!################################################################################ 
 recursive subroutine MLdmspq_triat(rank,ncoords,natoms,local,xyz,f)

    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
    real(ark),intent(out)  ::  f(rank)
    !
    integer(ik) :: I,J,K,L,M,ID
    integer(ik) :: POLORDER1,POLORDER2,MTYPE1,MTYPE2,NC1,NC2
    integer(ik), PARAMETER :: NX=3
    real(ark), DIMENSION(NX) :: Y
    integer(ik) :: TOTNUM,MNUM
    integer(ik), DIMENSION(3,6) :: P
    real(ark), ALLOCATABLE, DIMENSION(:) :: POLP,POLQ 
    real(ark)             :: mu(3),ux(3),uy(3),uz(3),tmat(3,3),n1(3),n2(3),x(2,3),r1,r2,alpha,re1,re2,ae
    real(ark)             :: cp(1:extF%nterms(1)-5),cq(1:extF%nterms(2)-2),v0,v1,v2,v3,v4,v5,v6,v7,v8,xyz0(natoms,3)
    !
    ! xyz are undefined for the local case
    !
    !write(out,"('MLdms2pqr_xyz_coeff is temporally deactivated as a bisector frame, use DIPOLE_PQR_XYZ_Z-FRAME instead')")
    !stop 'MLdms2pqr_xyz_coeff is temporally deactivated as a bisector frame, use DIPOLE_PQR_XYZ_Z-FRAME instead'
    !
    if (all(abs(xyz)<small_)) then 
      !
      xyz0 = MLloc2pqr_xyz(local)
      !
      x(1,:) = xyz0(2,:) - xyz0(1,:)
      x(2,:) = xyz0(3,:) - xyz0(1,:)
      !
    else
      !
      x(1,:) = xyz(2,:) - xyz(1,:)
      x(2,:) = xyz(3,:) - xyz(1,:)
      !
    endif
    !
    r1 = sqrt(sum(x(1,:)**2))
    r2 = sqrt(sum(x(2,:)**2))
    !
    n1 = x(1,:) / r1
    n2 = x(2,:) / r2
    !
    alpha = acos(sum(n1*n2))
    !
 
    select case(trim(molec%coords_transform))
 
    case('R1-Z-R2-RHO') 
! 
    uz = n1
    uz = uz / sqrt(sum(uz(:)**2))
    !
    uy = MLvector_product(uz,n2)
    !
    ! assume that ab initio embedding is with x = -x and y = - y 
    !
    if (sum(uy(:)**2)<sqrt(small_)) uy = (/0.0_ark,-1.0_ark,0.0_ark/)
    !
    uy = uy / sqrt(sum(uy(:)**2))
    !
    ux = MLvector_product(uy,uz)
    !
    ux = ux / sqrt(sum(ux(:)**2))
!
    case('R-RHO-Z')
!
    ux = n1 + n2
    if (sum(ux(:)**2)<sqrt(small_)) ux = (/-1.0_ark,0.0_ark,0.0_ark/)
    !
    ux = ux / sqrt(sum(ux(:)**2))
    !
    uy = MLvector_product(n1,ux)
    !
    uy = uy / sqrt(sum(uy(:)**2))
    !
    uz = MLvector_product(ux,uy)
!
    end select
 
    tmat(1, :) = ux
    tmat(2, :) = uy
    tmat(3, :) = uz
    !
    MTYPE1 = extF%coef(1,1)
    POLORDER1 = extF%coef(2,1)
    re1 = extF%coef(3,1)
    re2 = extF%coef(4,1)
    ae = extF%coef(5,1)*pi/180.0_ark
    !
      DO I=1,NX
        Y(I)=0.0D+00
      END DO
   
    Y(2) = (r1 - re1)
    Y(3) = (r2 - re2)
    Y(1) = cos(alpha) - cos(ae)
    !
    NC1 = extF%nterms(1)-5
    ALLOCATE(POLP(NC1))
    !
     DO I=1,NC1
        cp(I)=extF%coef(I+5,1)
     END DO 
    !
      DO I=1,NC1
        POLP(I)=0.00D+00
      END DO  
   
      TOTNUM=0
      DO M=0,POLORDER1
        MNUM=0
        DO I=0,M
          DO J=0,(M-I)
            K=M-I-J
!

              IF (MTYPE1 .EQ. 1) THEN  

                TOTNUM=TOTNUM+1
                MNUM=MNUM+1

                CALL PERMUTABC(I,J,K,P,ID)

                DO L=1,ID

                POLP(TOTNUM)=POLP(TOTNUM)+ (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))
  
                END DO

                POLP(TOTNUM)=cp(TOTNUM)*POLP(TOTNUM)/DBLE(ID)

              ELSE IF (MTYPE1 .EQ. 2) THEN 

                IF (J.LE.K) THEN 

                  TOTNUM=TOTNUM+1
                  MNUM=MNUM+1

                  CALL PERMUTAB2(I,J,K,P,ID)
                   
                  DO L=1,ID

                    POLP(TOTNUM)=POLP(TOTNUM)+ (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))

                  END DO

                  POLP(TOTNUM)=cp(TOTNUM)*POLP(TOTNUM)/DBLE(ID)

                END IF    
         
              ELSE IF (MTYPE1 .EQ. 3) THEN

                IF (I.LE.J .AND. J.LE.K) THEN

                  TOTNUM=TOTNUM+1
                  MNUM=MNUM+1

                  CALL PERMUTA3(I,J,K,P,ID)
                   
                  DO L=1,ID

                   POLP(TOTNUM)=POLP(TOTNUM)+ (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))

                  END DO

                  POLP(TOTNUM)=cp(TOTNUM)*POLP(TOTNUM)/DBLE(ID)

                END IF

              END IF

          END DO
        END DO

      END DO
    !
    mu(3) = SUM(POLP)
    !
    MTYPE2 = extF%coef(1,2)
    POLORDER2 = extF%coef(2,2) 
    !
    NC2 = extF%nterms(2)-2 
    ALLOCATE(POLQ(NC2))
    !
     DO I=1,NC2
        cq(I)=extF%coef(I+2,1)
     END DO 
    !
      DO I=1,NC2
        POLQ(I)=0.00D+00
      END DO  
   
      TOTNUM=0
      DO M=0,POLORDER2
        MNUM=0
        DO I=0,M
          DO J=0,(M-I)
            K=M-I-J
!

              IF (MTYPE2 .EQ. 1) THEN  

                TOTNUM=TOTNUM+1
                MNUM=MNUM+1

                CALL PERMUTABC(I,J,K,P,ID)

                DO L=1,ID

                POLQ(TOTNUM)=POLQ(TOTNUM)+ (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))
  
                END DO

                POLQ(TOTNUM)=cq(TOTNUM)*POLQ(TOTNUM)/DBLE(ID)

              ELSE IF (MTYPE2 .EQ. 2) THEN 

                IF (J.LE.K) THEN 

                  TOTNUM=TOTNUM+1
                  MNUM=MNUM+1

                  CALL PERMUTAB2(I,J,K,P,ID)
                   
                  DO L=1,ID

                    POLQ(TOTNUM)=POLQ(TOTNUM)+ (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))

                  END DO

                  POLQ(TOTNUM)=cq(TOTNUM)*POLQ(TOTNUM)/DBLE(ID)

                END IF    
         
              ELSE IF (MTYPE2 .EQ. 3) THEN

                IF (I.LE.J .AND. J.LE.K) THEN

                  TOTNUM=TOTNUM+1
                  MNUM=MNUM+1

                  CALL PERMUTA3(I,J,K,P,ID)
                   
                  DO L=1,ID

                   POLQ(TOTNUM)=POLQ(TOTNUM)+ (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))

                  END DO

                  POLQ(TOTNUM)=cq(TOTNUM)*POLQ(TOTNUM)/DBLE(ID)

                END IF

              END IF

          END DO
        END DO

      END DO
    !
    mu(1) = SUM(POLQ)*sin(pi-alpha)
    !
    mu(2) = 0
    !
    f(1:3) = matmul(tmat,mu)
 
 DEALLOCATE(POLP,POLQ)
    !
 end subroutine MLdmspq_triat

!######################################################################################################
! NH2 GROUND STATE THREE-BODY TERM
!######################################################################################################
      SUBROUTINE CHIPR_TRIAT(rr1,rr2,rr3,POT,C)
      IMPLICIT none
      real(ark),intent(in)   :: rr1,rr2,rr3
      real(ark),intent(out) :: POT
      CHARACTER(LEN=3), PARAMETER :: MOLTYP="AB2"
      integer(ik), PARAMETER :: DEG=2
      integer(ik), PARAMETER :: NC=160
      integer(ik), PARAMETER :: NX=3
      integer(ik) :: I,J,K,L,M,S,O,ID
      integer(ik) :: POLORDER
      integer(ik), DIMENSION(DEG) :: BSORDER,NCBAS
      real(ark), DIMENSION(NC) :: C
      real(ark), ALLOCATABLE, DIMENSION(:) :: POL,BS
      real(ark), DIMENSION(3) :: R
      real(ark) :: REPD,test
      integer(ik) :: NCPOL,NCTOTAL,SUMC
      real(ark), DIMENSION(NX) :: Y
      integer(ik) :: TOTNUM,MNUM
      integer(ik), DIMENSION(3,6) :: P

      R(1)=rr1
      R(2)=rr2
      R(3)=rr3

      BSORDER(1)=  4
      BSORDER(2)=  4

      POLORDER=  10

      DO I=1,DEG
        NCBAS(I)=0
        NCBAS(I)=2*BSORDER(I)+2
      END DO

      CALL POLNC(POLORDER,MOLTYP,NCPOL)

      IF (NC.NE.(NCPOL+SUM(NCBAS))) THEN 
        WRITE(*,*) "PROBLEMS IN DEFINING PROPER NUMBER OF COEFFICIENTS"
        WRITE(*,*) "TOTAL NUMBER OF COEFFS ARE NOT SUMMING UP CORRECTLY"
        STOP
      END IF

      ALLOCATE(POL(NCPOL))

      DO I=1,NCPOL
        POL(I)=0.00D+00

      END DO

      DO I=1,NX
        Y(I)=0.0D+00
  
      END DO

      SUMC=NCPOL
      DO I=1,DEG
        ALLOCATE(BS(NCBAS(I)))
        DO J=1,NCBAS(I)
          BS(J)=0.0D+00
        END DO
        K=SUMC
        DO J=1,NCBAS(I)
          K=K+1
          BS(J)=C(K)
        END DO

!######
! AB2-TYPE
!######
        IF (DEG.EQ.2) THEN
          IF (I.EQ.1) THEN
!
!           B-B BASIS
!
            CALL BASIS_CONTRACT(1,BSORDER(1),BS,R(1),Y(1))

          ELSE
!
!           A-B BASIS
!
            DO O=2,3 
              CALL BASIS_CONTRACT(2,BSORDER(2),BS,R(O),Y(O))

!######
            END DO
          END IF            
        END IF

        SUMC=SUMC+NCBAS(I)
        DEALLOCATE(BS)
      END DO

      TOTNUM=0
      S=0
      DO M=0,POLORDER
        MNUM=0
        DO I=0,M
          DO J=0,(M-I)
            K=M-I-J
            S=I+J+K
            IF (S.NE.I .AND. S.NE.J .AND. S.NE.K) THEN

              IF (MOLTYP.EQ."ABC") THEN  

                TOTNUM=TOTNUM+1
                MNUM=MNUM+1

                CALL PERMUTABC(I,J,K,P,ID)

                DO L=1,ID

                POL(TOTNUM)=POL(TOTNUM)+ (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))
  
                END DO

                POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)

              ELSE IF (MOLTYP.EQ."AB2") THEN 

                IF (J.LE.K) THEN 

                  TOTNUM=TOTNUM+1
                  MNUM=MNUM+1

                  CALL PERMUTAB2(I,J,K,P,ID)
                   
                  DO L=1,ID

                    POL(TOTNUM)=POL(TOTNUM)+ (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))

                  END DO

                  POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)

                END IF    
         
              ELSE IF (MOLTYP.EQ."A3") THEN

                IF (I.LE.J .AND. J.LE.K) THEN

                  TOTNUM=TOTNUM+1
                  MNUM=MNUM+1

                  CALL PERMUTA3(I,J,K,P,ID)
                   
                  DO L=1,ID

                   POL(TOTNUM)=POL(TOTNUM)+ (Y(1)**(P(1,L))*Y(2)**(P(2,L))*Y(3)**(P(3,L)))

                  END DO

                  POL(TOTNUM)=C(TOTNUM)*POL(TOTNUM)/DBLE(ID)

                END IF

              END IF

            END IF
          END DO
        END DO

      END DO

     CALL REPDAMP(NX,R,test)
     repd=test
      POT=SUM(POL)*repd


      DEALLOCATE(POL)

      RETURN
      END SUBROUTINE CHIPR_TRIAT

!################################################################################
      SUBROUTINE REPDAMP(NX,R,PEPDA)
      IMPLICIT none
      integer(ik) :: I,NX
      real(ark), DIMENSION(NX) :: R,H
      real(ark) :: KAPPA, XI, R0,PEPDA
      R0=0.5D+00
      KAPPA=100.0D+00
      XI=10.0D+00
      DO I=1,NX
        H(I)=0.5D+00*(1.00D+00+TANH(KAPPA*(R(I)-R0)))
      END DO
      PEPDA=(PRODUCT(H))**(XI)    
      RETURN  
      END SUBROUTINE REPDAMP
!################################################################################           
      SUBROUTINE PERMUTABC(I,J,K,P,ID)
      IMPLICIT real(ark) (A-H,O-Y),logical(z)
      integer(ik) :: I,J,K
      integer(ik) :: L,M,ID
      integer(ik), DIMENSION(3) :: INTER
      integer(ik), DIMENSION(3,6) :: P

      DO L=1,3 
        DO M=1,6
          P(L,M)=0
        END DO
      END DO

      INTER(1)=I
      INTER(2)=J
      INTER(3)=K

      ID=1

      P(1,1)=INTER(1)
      P(2,1)=INTER(2)
      P(3,1)=INTER(3) 

      RETURN
      END SUBROUTINE PERMUTABC
!################################################################################           
      SUBROUTINE PERMUTAB2(I,J,K,P,ID)
      IMPLICIT none
      integer(ik) :: I,J,K
      integer(ik) :: L,M,ID
      integer(ik), DIMENSION(3) :: INTER
      integer(ik), DIMENSION(3,6) :: P

      DO L=1,3 
        DO M=1,6
          P(L,M)=0
        END DO
      END DO
      
      INTER(1)=I
      INTER(2)=J
      INTER(3)=K
      
      IF (J.NE.K) THEN  
 
        ID=2

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3)      

        P(1,2)=INTER(1)
        P(2,2)=INTER(3)
        P(3,2)=INTER(2) 

       ELSE 

        ID=1

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3)

       END IF

      RETURN 
      END SUBROUTINE PERMUTAB2
!################################################################################          
      SUBROUTINE PERMUTA3(I,J,K,P,ID)
      IMPLICIT none
      integer(ik) :: I,J,K
      integer(ik) :: L,M,ID
      integer(ik), DIMENSION(3) :: INTER
      integer(ik), DIMENSION(3,6) :: P

      DO L=1,3 
        DO M=1,6
          P(L,M)=0
        END DO
      END DO
      
      INTER(1)=I
      INTER(2)=J
      INTER(3)=K
      
      IF (I.EQ.J.AND.J.EQ.K.AND.I.EQ.K) THEN
      
        ID=1

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3)  
        
      ELSE IF (I.NE.J.AND.J.NE.K.AND.I.NE.K) THEN   

        ID=6

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3)      

        P(1,2)=INTER(2)
        P(2,2)=INTER(1)
        P(3,2)=INTER(3)  

        P(1,3)=INTER(3)
        P(2,3)=INTER(2)
        P(3,3)=INTER(1) 

        P(1,4)=INTER(1)
        P(2,4)=INTER(3)
        P(3,4)=INTER(2) 

        P(1,5)=P(1,2)
        P(2,5)=P(3,2)
        P(3,5)=P(2,2)  

        P(1,6)=P(1,3)
        P(2,6)=P(3,3)
        P(3,6)=P(2,3)  
        
      ELSE IF (I.EQ.J) THEN       
     
        ID=3     

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3) 

        P(1,2)=INTER(3)
        P(2,2)=INTER(2)
        P(3,2)=INTER(1) 

        P(1,3)=INTER(1)
        P(2,3)=INTER(3)
        P(3,3)=INTER(2)
        
      ELSE IF (J.EQ.K) THEN       
                    
        ID=3     

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3)      

        P(1,2)=INTER(2)
        P(2,2)=INTER(1)
        P(3,2)=INTER(3)  

        P(1,3)=INTER(3)
        P(2,3)=INTER(2)
        P(3,3)=INTER(1) 
        
      ELSE IF (I.EQ.K) THEN    
  
        ID=3  

        P(1,1)=INTER(1)
        P(2,1)=INTER(2)
        P(3,1)=INTER(3)      

        P(1,2)=INTER(2)
        P(2,2)=INTER(1)
        P(3,2)=INTER(3)  

        P(1,3)=INTER(1)
        P(2,3)=INTER(3)
        P(3,3)=INTER(2)        
    
      END IF

      RETURN 
      END  SUBROUTINE PERMUTA3
!################################################################################ 
      SUBROUTINE POLNC(ORDER,MOLTYP,NC)
      IMPLICIT none
      integer(ik) :: I,J,K,L,M,S
      integer(ik) :: NC,MNUM,ORDER
      CHARACTER(LEN=3) :: MOLTYP

      NC=0
      S=0

      DO M=0,ORDER
        MNUM=0
        DO I=0,M
          DO J=0,(M-I)
            K=M-I-J
            S=I+J+K
            IF (S.NE.I .AND. S.NE.J .AND. S.NE.K) THEN
 
              IF (MOLTYP.EQ."ABC") THEN  

                NC=NC+1
                MNUM=MNUM+1

              ELSE IF (MOLTYP.EQ."AB2") THEN 

                IF (J.LE.K) THEN 

                  NC=NC+1
                  MNUM=MNUM+1

                END IF    
        
              ELSE IF (MOLTYP.EQ."A3") THEN

                IF (I.LE.J .AND. J.LE.K) THEN

                  NC=NC+1
                  MNUM=MNUM+1

                END IF

              END IF

            END IF
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE POLNC
!######################################################################################################
      SUBROUTINE BASIS_CONTRACT(DEG,M,C,R,YVAL)
      IMPLICIT none
      integer(ik) :: I,J,DEG
      integer(ik) :: M
      real(ark) :: RREF0,ZETA
      real(ark), DIMENSION(2*M+2) :: C
      real(ark), DIMENSION(M) :: GAMA,VAL
      real(ark) :: R,YVAL

      DO I=1,M
        GAMA(I)=C(M+I)
      END DO

      RREF0=C(2*M+1)
      ZETA=C(2*M+2)

      DO I=1,M-1
        CALL PHISECBASIS(DEG,I,RREF0,ZETA,GAMA(I),1,R,VAL(I))
      END DO

      DO I=M,M
        CALL PHICSECBASIS(DEG,I,RREF0,ZETA,GAMA(I),1,6,R,VAL(I))
      END DO

      YVAL=0.00D+00

      DO J=1,M
        YVAL=YVAL+C(J)*VAL(J)
      END DO

      RETURN
      END SUBROUTINE BASIS_CONTRACT
!######################################################################################################
      SUBROUTINE PHISECBASIS(DEG,IND,RREF0,ZETA,GAMA,ETA,R,VAL)
      IMPLICIT none
      integer(ik) :: DEG, IND, ETA
      real(ark) :: RREF0, ZETA, RREFIND 
      real(ark) :: GAMA, R, VAL, RHO
      RREFIND=ZETA*(RREF0)**(DBLE(IND)-1.0D+00)
      RHO=(R-RREFIND)
      VAL=(1.00D+00/(COSH(GAMA*RHO)))**(DBLE(ETA))
      RETURN
      END SUBROUTINE PHISECBASIS
!######################################################################################################
       SUBROUTINE PHICSECBASIS(DEG,IND,RREF0,ZETA,GAMA,ETA,LR,R,VAL)
      IMPLICIT none
      integer(ik) :: DEG, IND, LR, ETA
      real(ark) :: RREF0, ZETA, RREFIND 
      real(ark) :: GAMA, R, VAL, RHO
      real(ark) :: BETA, FAC
      BETA=1.00D+00/5.0D+00
      RREFIND=ZETA*(RREF0)**(DBLE(IND)-1.0D+00)
      RHO=(R-RREFIND)
      FAC=(TANH(BETA*R)/R)**(DBLE(LR))
      VAL=FAC*(1.00D+00/(COSH(GAMA*RHO)))**(DBLE(ETA))
      RETURN
      END SUBROUTINE PHICSECBASIS



end module pot_user
