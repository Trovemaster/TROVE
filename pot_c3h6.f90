! This unit defines all specific routines for a nine-atomic molecule of propene-type

module pot_c3h6
use accuracy
use moltype

implicit none

public MLpoten_c3h6_harmtest,MLpoten_c3h6_sym_II

private

integer(ik), parameter :: verbose = 4 ! Verbosity level


contains


function MLpoten_c3h6_harmtest(ncoords, natoms, local, xyz, force) result(f)

  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f

  real(ark) :: xi(21)
  real(ark) :: chi(21,6),term,rad
  integer(ik) :: ioper,ipower(21),i

  real(ark) :: r_12,r_23,r_14,r_15,r_16,r_38
  real(ark) :: r_39,r_27,alp_321,alp_214,alp_215
  real(ark) :: alp_216,alp_238,alp_239,alp_127
  real(ark) :: di_hcch,di_hcch2,di_hcch3,di_ccch
  real(ark) :: di_hcch4,di_hcchg
  real(ark) :: r_12e,r_23e,r_14e,r_15e,r_16e,r_38e
  real(ark) :: r_39e,r_27e,alp_321e,alp_214e,alp_215e
  real(ark) :: alp_216e,alp_238e,alp_239e,alp_127e
  real(ark) :: di_hcche,di_hcch2e,di_hcch3e,di_ccche
  real(ark) :: di_hcch4e,di_hcchge,tol
  real(ark) :: a1,a2,a3,a4,a5,a6,a7,a8,Ax1,Ax2,tbar,T1,T2,T3,T1_

  tol = 1.0d-11
  rad = pi/180.0_ark

  ! expansion functions

!   write(6,*) !BPM
!   write(6,*) local !BPM


  r_12      = local(1)
  r_23      = local(2)
  r_14      = local(3)
  r_15      = local(4)
  r_16      = local(5)
  r_38      = local(6)
  r_39      = local(7)
  r_27      = local(8)

  alp_321   = local(9)
  alp_214   = local(10)
  alp_215   = local(11)
  alp_216   = local(12)
  alp_238   = local(13)
  alp_239   = local(14)
  alp_127   = local(15)

  di_hcch   = local(16)
  di_hcch2  = local(17)
  di_hcch3  = local(18)
  di_ccch   = local(19)
  di_hcch4  = local(20)
  di_hcchg  = local(21)

  r_12e = force(1)
  r_23e = force(2)
  r_14e = force(3)
  r_15e = force(4)
  r_16e = force(5)
  r_38e = force(6)
  r_39e = force(7)
  r_27e = force(8)


  alp_321e   = force(9)*rad
  alp_214e   = force(10)*rad
  alp_215e   = force(11)*rad
  alp_216e   = force(12)*rad
  alp_238e   = force(13)*rad
  alp_239e   = force(14)*rad
  alp_127e   = force(15)*rad


  di_hcche   = force(16)*rad
  di_hcch2e  = force(17)*rad
  di_hcch3e  = force(18)*rad
  di_ccche   = force(19)*rad
  di_hcch4e  = force(20)*rad
  di_hcchge  = force(21)*rad

  a1 = force(22)
  a2 = force(23)
  a3 = force(24)
  a4 = force(25)
  a5 = force(26)
  a6 = force(27)
  a7 = force(28)
  a8 = force(29)


  select case(trim(molec%coords_transform))
    !
  case default
    !
    write(out, '(/a,1x,a,1x,a)') &
    'MLpoten_c3h6_harmtest error', trim(molec%coords_transform), 'is unknown'
    stop 'MLpoten_c3h6_harmtest error error: bad coordinate type'

  case('ZMAT_7ALF_5TAU','C3H6_7ALF_5TAU')

    xi(1)=1.0_ark-exp(-a1*(r_12-r_12e))
    xi(2)=1.0_ark-exp(-a2*(r_23-r_23e))
    xi(3)=1.0_ark-exp(-a3*(r_14-r_14e))
    xi(4)=1.0_ark-exp(-a4*(r_15-r_15e))
    XI(5)=1.0_ark-EXP(-a5*(r_16-r_16e))
    XI(6)=1.0_ark-EXP(-a7*(r_38-r_38e))
    XI(7)=1.0_ark-EXP(-a8*(r_39-r_39e))
    XI(8)=1.0_ark-EXP(-a6*(r_27-r_27e))
  !
    xi(9) =  alp_321   -  alp_321e
    xi(10) = alp_214   -  alp_214e
    xi(11) = alp_215   -  alp_215e
    xi(12) = alp_216   -  alp_216e
    xi(13) = alp_238   -  alp_238e
    xi(14) = alp_239   -  alp_239e
    xi(15) = alp_127   -  alp_127e

    xi(16) =  di_ccch   -  di_ccche 
    !To ensure always near zero:
    !if( abs( xi(16)-2.0_ark*pi).lt.1d-4 ) then
    if(xi(16).gt.6.d0) then
    xi(16) = xi(16) - 2.0_ark*pi
    end if
    xi(17) =  di_hcch4  -  di_hcch4e
    xi(18) =  di_hcchg  -  di_hcchge


    T1 = di_hcch


!       ALWAYS WANT TO MEASURE ANGLES IN SAME WAY
!       TO GET CONSISTENT TRANSFORMS, DEFINE RANGES FOR ANGLES ELSE ADD/SUB 2*PI

    IF(T1.LT.0.0.AND.ABS(T1).GT.TOL) THEN
    T1 = T1 + 2.0_ark*PI
    END IF
    IF(T1.GT.(2.0_ark*PI).AND.ABS(T1-2.0_ARK*PI).GT.TOL) THEN
    T1 = T1 - 2.0_ark*PI
    END IF


!    T2 = di_hcch2 + di_hcch
!    T3 = di_hcch3 + di_hcch

    T2 = di_hcch2 + T1
    T3 = di_hcch3 + T1

!    IF(T2.LT.2.0_ark*Pi/3.0_ark.AND.abs(T2-2.0_ark*Pi/3.0_ark).GT.tol) THEN
!    T2 = T2 + 2.0_ark*PI
!    END IF
!    IF(T2.GT.8.0_ark*Pi/3.0_ark.AND.(T2-8.0_ark*Pi/3.0_ark).GT.tol) THEN
!    T2 = T2 - 2.0_ark*PI
!    END IF

!    IF(T3.LT.(4.0_ark*PI/3.0_ark).AND.ABS(T3-4.0_ark*PI/3.0_ark).GT.TOL) THEN
!    T3 = T3 + 2.0_ark*PI
!    END IF
!    IF(T3.GT.(10.0_ark*PI/3.0_ark).AND.ABS(T3-10.0_ark*PI/3.0_ark).GT.TOL) THEN
!    T3 = T3 - 2.0_ark*PI
!    END IF

    T1 = T1 - force(16)*rad
    T2 = T2 - force(17)*rad
    T3 = T3 - force(18)*rad

!       SUBTRACT EQUILBRIUM THETA VALUES TO MAKE A1/A2 ZERO AT EQUILIBRIUM
!       AND ENSURE CONSISTENT TRANSFROMS


    Ax1  = (2.0_ark*T1 - T2 - T3)/sqrt(6.0_ark)
    Ax2  = (             T2 - T3)/sqrt(2.0_ark)
    tbar = (T1 + T2 + T3)/3.0_ark

    xi(19) = Ax1 
    xi(20) = Ax2 
    xi(21) = (1.0_ark - cos(3.d0*tbar) )
    !
  case('7ALF_TAU_1RHO')
    !
    xi(1)=1.0_ark-exp(-a1*(r_12-r_12e))
    xi(2)=1.0_ark-exp(-a2*(r_23-r_23e))
    xi(3)=1.0_ark-exp(-a3*(r_14-r_14e))
    xi(4)=1.0_ark-exp(-a4*(r_15-r_15e))
    XI(5)=1.0_ark-EXP(-a5*(r_16-r_16e))
    XI(6)=1.0_ark-EXP(-a7*(r_38-r_38e))
    XI(7)=1.0_ark-EXP(-a8*(r_39-r_39e))
    XI(8)=1.0_ark-EXP(-a6*(r_27-r_27e))
  !
    xi(9) =  alp_321   -  alp_321e
    xi(10) = alp_214   -  alp_214e
    xi(11) = alp_215   -  alp_215e
    xi(12) = alp_216   -  alp_216e
    xi(13) = alp_238   -  alp_238e
    xi(14) = alp_239   -  alp_239e
    xi(15) = alp_127   -  alp_127e

    xi(16) =  di_ccch   -  di_ccche 
    !
    !To ensure always near zero:
    !if( abs( xi(16)-2.0_ark*pi).lt.1d-4 ) then
    if(xi(16).gt.6.d0) then
       xi(16) = xi(16) - 2.0_ark*pi
    end if
    xi(17) =  di_hcch4  -  di_hcch4e
    xi(18) =  di_hcchg  -  di_hcchge
    !
    T1 = di_hcch
    !
    !       ALWAYS WANT TO MEASURE ANGLES IN SAME WAY
    !       TO GET CONSISTENT TRANSFORMS, DEFINE RANGES FOR ANGLES ELSE ADD/SUB 2*PI
    !
    IF(T1.LT.0.0.AND.ABS(T1).GT.TOL) THEN
      T1 = T1 + 2.0_ark*PI
    END IF
    IF(T1.GT.(2.0_ark*PI).AND.ABS(T1-2.0_ARK*PI).GT.TOL) THEN
      T1 = T1 - 2.0_ark*PI
    END IF
    !
    T2 = di_hcch2
    T3 = di_hcch3
    !
    !IF(T2.LT.2.0_ark*Pi/3.0_ark.AND.abs(T2-2.0_ark*Pi/3.0_ark).GT.tol) THEN
    !T2 = T2 + 2.0_ark*PI
    !END IF
    !IF(T2.GT.8.0_ark*Pi/3.0_ark.AND.(T2-8.0_ark*Pi/3.0_ark).GT.tol) THEN
    !T2 = T2 - 2.0_ark*PI
    !END IF

    !IF(T3.LT.(4.0_ark*PI/3.0_ark).AND.ABS(T3-4.0_ark*PI/3.0_ark).GT.TOL) THEN
    !T3 = T3 + 2.0_ark*PI
    !END IF
    !IF(T3.GT.(10.0_ark*PI/3.0_ark).AND.ABS(T3-10.0_ark*PI/3.0_ark).GT.TOL) THEN
    !T3 = T3 - 2.0_ark*PI
    !END IF
    !
    T1 = T1 - force(16)*rad
    T2 = T2 - force(17)*rad
    T3 = T3 - force(18)*rad
    !
    Ax1  = (T2 + T3)*3.0_ark/sqrt(6.0_ark)   
    Ax2  = (T2 - T3)/sqrt(2.0_ark)   
    !
    !tbar = (src(16) + src(17) + src(18))/3.0_ark 
    !
    T2 = T2 + T1
    T3 =-T3 + T1
    !
    tbar = (T1 + T2 + T3)/3.0_ark 


!       SUBTRACT EQUILBRIUM THETA VALUES TO MAKE A1/A2 ZERO AT EQUILIBRIUM
!       AND ENSURE CONSISTENT TRANSFROMS
    !

    xi(19) = Ax1 
    xi(20) = Ax2 
    xi(21) = (1.0_ark - cos(3.d0*tbar) )
    !
   end select

      !  write(6,*) !BPM
      !  write(6,*) xi !BPM


  f = 0
  !
!  do ioper = 1,12
    !
!    call ML_symmetry_transformation_XY3_II(ioper,xi,chi(:,ioper),18)
    !
!  enddo
  ! 
  do i = 30, molec%parmax

    ipower(1:21) = molec%pot_ind(1:21,i)

!!    term = 0 

!    do ioper = 1,12
   
      term = product(xi(1:21)**ipower(1:21))   !no symmetry invoked here,
                                                      !just test of harmonic pot

!    end do

!    term = term/12.0_ark

    f = f + term*force(i)

  end do

      !  write(6,*) !BPM
      !  write(6,*) f!BPM
      !  stop ! BPM
        if(f.gt.1e6) then 
          write(out,"('MLpoten_c3h6_harmtest error, V is too large:',e12.5)") f
          write(out,"(21f14.8)") xi
          !stop 'problem here'
        endif


  !
  !
  !f = force(6)+force(7)*xi(1)**2+&
  !    force(8)*sum(xi(2:7)**2)+force(8)*xi(18)**2+&
  !    force(9)*sum(xi(8:13)**2)+force(10)*sum(xi(14:17)**2)+&
  !    force(11)*( xi(14)*xi(16)+xi(15)*xi(17) )
  !!
  !f = f/12.0_ark
  !  

end function  MLpoten_c3h6_harmtest



function MLpoten_c3h6_sym_II(ncoords, natoms, local, xyz, force) result(f)

  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f

  real(ark) :: xi(21)
  real(ark) :: chi(21,6),term,deg
  integer(ik) :: ioper,ipower(21),i,nparam_eq
  !
  real(ark) :: rC1,rC2,rH1 ,rH2 ,rH3 ,rH4 ,rH5 ,rH6 ,rC1e,rC2e,rH1e,rH2e,rH4e,rH5e,rH6e,&
                      alpha1e,alpha2e,alpha3e ,alpha4e ,alpha5e,alpha6e,taue,delta1e ,delta2e ,delta3e ,&
                      delta4e ,delta5e ,delta6e 
  real(ark) :: b1,b2,b3,b4,b5,theta1,theta2,theta3,a1,a2,tau
  !
  nparam_eq = 22
  deg = pi/180.0_ark
  !
  ! expansion functions
  !
  rC1e      = force(1)
  rC2e      = force(2)
  rH1e      = force(3)
  rH4e      = force(4)
  rH5e      = force(5)
  rH6e      = force(6)
  !
  alpha1e    = force( 7)*deg
  alpha2e    = force( 8)*deg
  alpha4e    = force( 9)*deg
  alpha5e    = force(10)*deg
  alpha6e    = force(11)*deg
  !
  delta1e    = force(12)*deg
  delta2e    = force(13)*deg
  delta4e    = force(14)*deg
  delta5e    = force(15)*deg
  delta6e    = force(16)*deg
  !
  a1       = force(17)
  a2       = force(18)
  b1       = force(19)
  b3       = force(20)
  b4       = force(21)
  b5       = force(22)


  select case(trim(molec%coords_transform))
    !
  case default
    !
    write(out, '(/a,1x,a,1x,a)') &
    'MLpoten_c3h6_harmtest error', trim(molec%coords_transform), 'is unknown'
    stop 'MLpoten_c3h6_harmtest error error: bad coordinate type'

      !
  case('7ALF_TAU_1RHO')
    !
    rC1      = local(1)
    rC2      = local(2)
    rH1      = local(3)
    rH2      = local(4)
    rH3      = local(5)
    rH4      = local(6)
    rH5      = local(7)
    rH6      = local(8)
    !
    ! C-C
    xi(1)=1.0d+00-exp(-a1*(rC1-rC1e))
    xi(2)=1.0d+00-exp(-a2*(rC2-rC2e))
    !CH3
    xi(3)=1.0d+00-exp(-b1*(rH1-rH1e))
    xi(4)=1.0d+00-exp(-b1*(rH2-rH1e))
    xi(5)=1.0d+00-exp(-b1*(rH3-rH1e))
    !CHs
    xi(6)=1.0d+00-exp(-b3*(rH4-rH4e))
    xi(7)=1.0d+00-exp(-b4*(rH5-rH5e))
    xi(8)=1.0d+00-exp(-b5*(rH6-rH6e))
    !
    !CCH1
    xi(9) = local(9)- alpha1e
    !xi(9) = cos(local(9))- cos(alpha1e)
    !
    ! alphas
    !H2-C-H3,H1-C-H2,H1-C-H3
    xi(10) = local(10)- alpha2e
    xi(11) = local(11)- alpha2e
    xi(12) = local(12)- alpha2e
    !
    !other CHs
    xi(13) = local(13)- alpha4e
    xi(14) = local(14)- alpha5e
    xi(15) = local(15)- alpha6e
    !
    xi(16) = local(19)- delta4e
    xi(17) = local(20)- delta5e
    xi(18) = local(21)- delta6e
    !
    theta2 = local(17)- delta2e
    theta3 = local(18)- delta2e
    !
    theta1 = 2.0_ark*pi-theta2-theta3
    !
    xi(19)  = ( theta2 + theta3 )*3.d0/sqrt(6.0_ark)
    xi(20)  = ( theta2 - theta3 )/sqrt(2.0_ark)
    !
    theta1 = local(16)-delta1e
    !
    theta2 = theta2 + theta1
    theta3 =-theta3 + theta1
    !
    tau = (theta1 + theta2 + theta3)/3.0_ark
    !
    xi(21) = 1.0_ark - cos(3.0_ark*tau)
    !
  end select
  !
  f = 0
  !
  call ML_symmetry_transformation_XY3_II(6,xi,chi,21)
  ! 
  do i = nparam_eq+1, molec%parmax
    !
    ipower(1:21) = molec%pot_ind(1:21,i)
    !
    term = 0 
    !
    do ioper = 1,6
      !
      term = term + product(chi(1:21,ioper)**ipower(1:21))
      !
    end do
    !
    term = term/6.0_ark
    !
    f = f + term*force(i)
    !
  end do
  !  
end function  MLpoten_c3h6_sym_II







subroutine ML_symmetry_transformation_XY3_II(nsym,src,dst,ndeg)
    implicit none
    !
    integer(rk),intent(in)    :: nsym,ndeg  ! group operation  
    real(ark),intent(in)      :: src(1:ndeg)
    real(ark),intent(out)     :: dst(1:ndeg,6)
    !
    real(ark)         :: repres(6,ndeg,ndeg),a,b,e,o
    integer(rk)       :: ioper
    !
    if (verbose>=5) write(out,"('ML_symmetry_transformation_XY3/start')")
    !
    a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
    !
    repres = 0
    !
    repres(:,1,1) = 1.0_ark
    repres(:,2,2) = 1.0_ark
    !
    repres(:,6,6) = 1.0_ark
    repres(:,7,7) = 1.0_ark
    repres(:,8,8) = 1.0_ark
    repres(:,9,9) = 1.0_ark
    !
    repres(:,13,13) = 1.0_ark
    repres(:,14,14) = 1.0_ark
    repres(:,15,15) = 1.0_ark
    !
    repres(:,16,16) = 1.0_ark
    repres(:,17,17) = 1.0_ark
    repres(:,18,18) = 1.0_ark
    repres(:,21,21) = 1.0_ark

    !
    ! E
    ! r123
    repres(1,3,3) = 1.0_ark
    repres(1,4,4) = 1.0_ark
    repres(1,5,5) = 1.0_ark
    ! a123
    repres(1,10,10) = 1.0_ark
    repres(1,11,11) = 1.0_ark
    repres(1,12,12) = 1.0_ark
    !d12
    repres(1,16,16) = 1.0_ark
    repres(1,17,17) = 1.0_ark
    !
    !C3+/(123)
    repres(2,3,4) = 1.0_ark
    repres(2,4,5) = 1.0_ark
    repres(2,5,3) = 1.0_ark
    !
    repres(2,10,11) = 1.0_ark
    repres(2,11,12) = 1.0_ark
    repres(2,12,10) = 1.0_ark
    !
    repres(2,19,19) = -a
    repres(2,20,19) = -b
    repres(2,19,20) =  b
    repres(2,20,20) = -a
    !
    !C3-/(132)
    !
    repres(3,3,5) = 1.0_ark
    repres(3,4,3) = 1.0_ark
    repres(3,5,4) = 1.0_ark
    !
    repres(3,10,12) =  1.0_ark
    repres(3,11,10) = 1.0_ark
    repres(3,12,11) = 1.0_ark
    !
    repres(3,19,19) = -a
    repres(3,20,19) =  b
    repres(3,19,20) = -b
    repres(3,20,20) = -a
    !
    !C2/(23)->(45)
    !
    repres(4,3,3) = 1.0_ark
    repres(4,4,5) = 1.0_ark
    repres(4,5,4) = 1.0_ark
    !
    repres(4,10,10) = 1.0_ark
    repres(4,11,12) =  1.0_ark
    repres(4,12,11) = 1.0_ark
    !
    repres(4,19,19) =  1.0_ark
    repres(4,20,20) = -1._ark
    !
    !C2'/(12)->(34)
    repres(5,3,4) = 1.0_ark
    repres(5,4,3) = 1.0_ark
    repres(5,5,5) = 1.0_ark
    !
    repres(5,9,10) = 1.0_ark
    repres(5,10,9) = 1.0_ark
    repres(5,11,11) =  1.0_ark
    !
    repres(5,19,19) = -a
    repres(5,19,20) = b
    repres(5,20,19) = b
    repres(5,20,20) = a
    !
    !(13)->(35)
    repres(6,3,5) = 1.0_ark
    repres(6,4,4) = 1.0_ark
    repres(6,5,3) = 1.0_ark
    !
    repres(6,10,12) = 1.0_ark
    repres(6,11,11) = 1.0_ark
    repres(6,12,10) =  1.0_ark
    repres(6,19,19) = -a
    repres(6,19,20) = -b
    repres(6,20,19) = -b
    repres(6,20,20) = a
    !
    do ioper = 1,nsym
      dst(:,ioper) = matmul(repres(ioper,:,:),src) 
    enddo

  end subroutine ML_symmetry_transformation_XY3_II





recursive subroutine ML_dipole_c3h6_4m_dummy(rank,ncoords,natoms,local,xyz0,f)
    !
    implicit none
    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords),xyz0(natoms,3)
    real(ark)              ::  f(rank)
    !
    !
    stop 'dipole_c3h6_4m_dummy is not implemented'
    !
end subroutine ML_dipole_c3h6_4m_dummy



end module pot_c3h6
