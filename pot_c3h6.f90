! This unit defines all specific routines for a nine-atomic molecule of propene-type

module pot_c3h6
use accuracy
use moltype

implicit none

public MLpoten_c3h6_harmtest

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
  real(ark) :: a1,a2,a3,a4,a5,a6,a7,a8,Ax1,Ax2,tbar,T1,T2,T3

  tol = 1.0d-11
  rad = pi/180.0_ark

  ! expansion functions


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
  di_hcch4   = local(20)
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
  di_hcch4e   = force(20)*rad
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
    !
  case('ZMAT_7ALF_5TAU','C3H6_7ALF_5TAU')
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
   end select


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












subroutine ML_symmetry_transformation_XY3_II(ioper,src,dst,NDEG)
!       THIS COULD ALSO BE MADE AS GENERAL AS POSSIBLE BUT WILL ALWAYS NEED TO
!       DEFINE MATRICES FOR EACH MOLECULE BPM
    implicit none
    !
    integer,intent(in)    :: ioper,NDEG  ! group operation  
    real(ark),intent(in)      :: src(1:NDEG)
    real(ark),intent(out)     :: dst(1:NDEG)
    !
    double precision         :: repres(12,NDEG,NDEG),a,b,e,o
    INTEGER NSYM
    !
    if (verbose>=5) write(out,"('ML_symmetry_transformation_XY3/start')")
    !
    a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
    NSYM = 12
    !
    repres = 0
    !
    ! E
    repres(1,1,1) = 1.0_ark
    repres(1,2,2) = 1.0_ark
    repres(1,3,3) = 1.0_ark
    repres(1,4,4) = 1.0_ark
    repres(1,5,5) = 1.0_ark
    repres(1,6,6) = 1.0_ark
    repres(1,7,7) = 1.0_ark
    repres(1,8,8) = 1.0_ark
    repres(1,9,9) = 1.0_ark
    REPRES(1,10,10) = 1.0_ark
    REPRES(1,11,11) = 1.0_ark
    REPRES(1,12,12) = 1.0_ark
    REPRES(1,13,13) = 1.0_ark
    REPRES(1,14,14) = 1.0_ark
    REPRES(1,15,15) = 1.0_ark
    REPRES(1,16,16) = 1.0_ark
    REPRES(1,17,17) = 1.0_ark
    REPRES(1,18,18)  = 1.0_ark
!
    !C3+/(123)(465)
    repres(2,1,1) = 1.0_ark
    repres(2,4,2) = 1.0_ark
    repres(2,2,3) = 1.0_ark
    repres(2,3,4) = 1.0_ark
    repres(2,6,5) = 1.0_ark
    repres(2,7,6) = 1.0_ark
    repres(2,5,7) = 1.0_ark
    repres(2,10,8) = 1.0_ark
    repres(2,8,9) = 1.0_ark
    repres(2,9,10) =  1.0_ark
    repres(2,12,11) = 1.0_ark
    repres(2,13,12) = 1.0_ark
    repres(2,11,13) = 1.0_ark
    repres(2,14,14) = -a
    repres(2,15,14) = b
    repres(2,14,15) = -b
    repres(2,15,15) = -a
    repres(2,16,16) = -a
    repres(2,17,16) = b
    repres(2,16,17) = -b
    repres(2,17,17) = -a
    REPRES(2,18,18)  = 1.0_ark
    !
!C3-/(123)(465)
    repres(3,1,1) = 1.0_ark
    repres(3,3,2) = 1.0_ark
    repres(3,4,3) = 1.0_ark
    repres(3,2,4) = 1.0_ark
    repres(3,7,5) = 1.0_ark
    repres(3,5,6) = 1.0_ark
    repres(3,6,7) = 1.0_ark
    repres(3,9,8) = 1.0_ark
    repres(3,10,9) = 1.0_ark
    repres(3,8,10) =  1.0_ark
    repres(3,13,11) = 1.0_ark
    repres(3,11,12) = 1.0_ark
    repres(3,12,13) = 1.0_ark
    repres(3,14,14) = -a
    repres(3,15,14) = -b
    repres(3,14,15) = b
    repres(3,15,15) = -a
    repres(3,16,16) = -a
    repres(3,17,16) = -b
    repres(3,16,17) = b
    repres(3,17,17) = -a
    REPRES(3,18,18)  = 1.0_ark
    !
    !C2/(16)(24)(35)(78)
    repres(4,1,1) = 1.0_ark
    repres(4,7,2) = 1.0_ark
    repres(4,5,3) = 1.0_ark
    repres(4,6,4) = 1.0_ark
    repres(4,3,5) = 1.0_ark
    repres(4,4,6) = 1.0_ark
    repres(4,2,7) = 1.0_ark
    repres(4,13,8) = 1.0_ark
    repres(4,11,9) = 1.0_ark
    repres(4,12,10) =  1.0_ark
    repres(4,9,11) = 1.0_ark
    repres(4,10,12) = 1.0_ark
    repres(4,8,13) = 1.0_ark
    repres(4,16,14) = 1._ark
    repres(4,17,15) = -1._ark
    repres(4,14,16) = 1._ark
    repres(4,15,17) = -1._ark
    REPRES(4,18,18)  = 1.0_ark
    !
    !
    !C2'/(16)(24)(35)(78)
    repres(5,1,1) = 1.0_ark
    repres(5,6,2) = 1.0_ark
    repres(5,7,3) = 1.0_ark
    repres(5,5,4) = 1.0_ark
    repres(5,4,5) = 1.0_ark
    repres(5,2,6) = 1.0_ark
    repres(5,3,7) = 1.0_ark
    repres(5,12,8) = 1.0_ark
    repres(5,13,9) = 1.0_ark
    repres(5,11,10) =  1.0_ark
    repres(5,10,11) = 1.0_ark
    repres(5,8,12) = 1.0_ark
    repres(5,9,13) = 1.0_ark
    repres(5,16,14) = -a
    repres(5,17,14) = -b
    repres(5,16,15) = -b
    repres(5,17,15) = a
    repres(5,14,16) = -a
    repres(5,15,16) = -b
    repres(5,14,17) = -b
    repres(5,15,17) =  a
    REPRES(5,18,18)  = 1.0_ark
    !
    !C2''/(16)(24)(35)(78)
    repres(6,1,1) = 1.0_ark
    repres(6,5,2) = 1.0_ark
    repres(6,6,3) = 1.0_ark
    repres(6,7,4) = 1.0_ark
    repres(6,2,5) = 1.0_ark
    repres(6,3,6) = 1.0_ark
    repres(6,4,7) = 1.0_ark
    repres(6,11,8) = 1.0_ark
    repres(6,12,9) = 1.0_ark
    repres(6,13,10) =  1.0_ark
    repres(6,8,11) = 1.0_ark
    repres(6,9,12) = 1.0_ark
    repres(6,10,13) = 1.0_ark
    repres(6,16,14) = -a
    repres(6,17,14) = b
    repres(6,16,15) = b
    repres(6,17,15) = a
    repres(6,14,16) = -a
    repres(6,14,16) = -a
    repres(6,15,16) = b
    repres(6,14,17) = b
    repres(6,15,17) = a
    REPRES(6,18,18)  = 1.0_ark
    !
    !i/(14)(26)(35)(78)*
    repres(7,1,1) = 1.0_ark
    repres(7,5,2) = 1.0_ark
    repres(7,7,3) = 1.0_ark
    repres(7,6,4) = 1.0_ark
    repres(7,2,5) = 1.0_ark
    repres(7,4,6) = 1.0_ark
    repres(7,3,7) = 1.0_ark
    repres(7,11,8) = 1.0_ark
    repres(7,13,9) = 1.0_ark
    repres(7,12,10) =  1.0_ark
    repres(7,8,11) = 1.0_ark
    repres(7,10,12) = 1.0_ark
    repres(7,9,13) = 1.0_ark
    repres(7,16,14) = 1._ark
    repres(7,17,15) = 1._ark
    repres(7,14,16) = 1._ark
    repres(7,15,17) = 1._ark
    REPRES(7,18,18)  = -1.0_ark
    !
    !S6/(163425)(78)*
    repres(8,1,1) = 1.0_ark
    repres(8,6,2) = 1.0_ark
    repres(8,5,3) = 1.0_ark
    repres(8,7,4) = 1.0_ark
    repres(8,4,5) = 1.0_ark
    repres(8,3,6) = 1.0_ark
    repres(8,2,7) = 1.0_ark
    repres(8,12,8) = 1.0_ark
    repres(8,11,9) = 1.0_ark
    repres(8,13,10) =  1.0_ark
    repres(8,10,11) = 1.0_ark
    repres(8,9,12) = 1.0_ark
    repres(8,8,13) = 1.0_ark
    repres(8,16,14) = -a
    repres(8,17,14) = b
    repres(8,16,15) = -b
    repres(8,17,15) = -a
    repres(8,14,16) = -a
    repres(8,15,16) = b
    repres(8,14,17) = -b
    repres(8,15,17) = -a
    REPRES(8,18,18)  = -1.0_ark
    !
    !S6'/(14)(26)(35)(78)*
    repres(9,1,1) = 1.0_ark
    repres(9,7,2) = 1.0_ark
    repres(9,6,3) = 1.0_ark
    repres(9,5,4) = 1.0_ark
    repres(9,3,5) = 1.0_ark
    repres(9,2,6) = 1.0_ark
    repres(9,4,7) = 1.0_ark
    repres(9,13,8) = 1.0_ark
    repres(9,12,9) = 1.0_ark
    repres(9,11,10) =  1.0_ark
    repres(9,9,11) = 1.0_ark
    repres(9,8,12) = 1.0_ark
    repres(9,10,13) = 1.0_ark
    repres(9,16,14) = -a
    repres(9,17,14) = -b
    repres(9,16,15) = b
    repres(9,17,15) = -a
    repres(9,14,16) = -a
    repres(9,15,16) = -b
    repres(9,14,17) = b
    repres(9,15,17) = -a
    REPRES(9,18,18)  = -1.0_ark
    !
 !sigmad/(12)(46)*
    repres(10,1,1) = 1.0_ark
    repres(10,4,2) = 1.0_ark
    repres(10,3,3) = 1.0_ark
    repres(10,2,4) = 1.0_ark
    repres(10,6,5) = 1.0_ark
    repres(10,5,6) = 1.0_ark
    repres(10,7,7) = 1.0_ark
    repres(10,10,8) = 1.0_ark
    repres(10,9,9) = 1.0_ark
    repres(10,8,10) =  1.0_ark
    repres(10,12,11) = 1.0_ark
    repres(10,11,12) = 1.0_ark
    repres(10,13,13) = 1.0_ark
    repres(10,14,14) = -a
    repres(10,15,14) = -b
    repres(10,14,15) = -b
    repres(10,15,15) = a
    repres(10,16,16) = -a
    repres(10,17,16) = -b
    repres(10,16,17) = -b
    repres(10,17,17) = a
    REPRES(10,18,18)  = -1.0_ark
!sigmad'/(12)(46)*
    repres(11,1,1) = 1.0_ark
    repres(11,2,2) = 1.0_ark
    repres(11,4,3) = 1.0_ark
    repres(11,3,4) = 1.0_ark
    repres(11,5,5) = 1.0_ark
    repres(11,7,6) = 1.0_ark
    repres(11,6,7) = 1.0_ark
    repres(11,8,8) = 1.0_ark
    repres(11,10,9) = 1.0_ark
    repres(11,9,10) =  1.0_ark
    repres(11,11,11) = 1.0_ark
    repres(11,13,12) = 1.0_ark
    repres(11,12,13) = 1.0_ark
    repres(11,14,14) =  -a
    repres(11,15,14) =  b
    repres(11,14,15) =  b
    repres(11,15,15) =  a
    repres(11,16,16) = -a
    repres(11,17,16) =  b
    repres(11,16,17) =  b
    repres(11,17,17) =  a
    REPRES(11,18,18)  = -1.0_ark
    !
!sigmad''/(12)(46)*
    repres(12,1,1) = 1.0_ark
    repres(12,3,2) = 1.0_ark
    repres(12,2,3) = 1.0_ark
    repres(12,4,4) = 1.0_ark
    repres(12,7,5) = 1.0_ark
    repres(12,6,6) = 1.0_ark
    repres(12,5,7) = 1.0_ark
    repres(12,9,8) = 1.0_ark
    repres(12,8,9) = 1.0_ark
    repres(12,10,10) =  1.0_ark
    repres(12,13,11) = 1.0_ark
    repres(12,12,12) = 1.0_ark
    repres(12,11,13) = 1.0_ark
    repres(12,14,14) = 1.0_ark
    repres(12,15,15) =  -1.0_ark
    repres(12,16,16) = 1.0_ark
    repres(12,17,17) =  -1.0_ark
    REPRES(12,18,18)  = -1.0_ark
    !
    if (ioper<0.or.ioper>NSYM) then
      write (6,"('symmetry_transformation_local: operation ',i8,' unknown')")
      stop 'symmetry_transformation_local - bad operation. type'
    endif
    !
    ! BPM : REVERSED MULTIPLICATION ORDER TO AGREE WITH MY NOTES
    dst = matmul(src,repres(ioper,:,:))

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
