! This unit defines all specific routines for a six-atomic molecule of ethane-type

module pot_c2h6
use accuracy
use moltype

implicit none

public MLpoten_c2h6_88,ML_dipole_c2h6_4m_dummy

private

integer(ik), parameter :: verbose = 4 ! Verbosity level


contains

function MLpoten_c2h6_88(ncoords, natoms, local, xyz, force) result(f)

  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f

  real(ark) :: xi(18),r1,r2,r3,r4,r5,r6,r7,r1e,r2e,betae,a,b
  real(ark) :: chi(18,12),term,rad
  integer(ik) :: ioper,ipower(18),i

  real(ark) :: tau14,tau24,tau25,tau35,tau36,theta12,theta23,theta13,theta56,theta45,theta46


  rad = pi/180.0_ark

  ! expansion functions

! TRANSFORM HERE FROM Z-MATRIX COORDINATES TO SYMMETRIZED COORDS
! AS FITTING WAS CARRIED OUT USING THESE COORDINATES

  r1      = local(1)
  r2      = local(2)
  r3      = local(3)
  r4      = local(4)
  r5      = local(5)
  r6      = local(6)
  r7      = local(7)

  r1e = force(1)
  r2e = force(2)
  betae = force(3)*rad
  a = force(4)
  b = force(5)


  select case(trim(molec%coords_transform))
    !
  case default
    !
    write(out, '(/a,1x,a,1x,a)') &
    'MLpoten_c2h6_88 error', trim(molec%coords_transform), 'is unknown'
    stop 'MLpoten_c2h6_88 error error: bad coordinate type'
    !
  case('ZMAT_4BETA_1TAU','C2H6_4BETA_1TAU')
    ! ORDER CHANGED HERE AS WANT TBAR TO BE 18TH COORDINATE BUT NEEDS
    !  TO BE 16TH COORDINATE IN Z-MAT
    !

    xi(14) = (3.0_ark*local(14)  - 2.0_ark*pi)/sqrt(6.0_ark)
    xi(15) = ( 2.0_ark*local(15) - 2.0_ark*pi + local(14))/sqrt(2.0_ark)
    xi(16) = (3.0_ark*local(18)  - 2.0_ark*pi)/sqrt(6.0_ark)
    xi(17) = ( 2.0_ark*local(17) - 2.0_ark*pi + local(18))/sqrt(2.0_ark)
    xi(18) = ( ( 3.0_ark*local(16) + local(14) - local(15) + local(17) - local(18))/3.0_ark ) - pi

    !
  case('R-R16-BETA16-THETA-TAU')
      !
      tau14 = local(14)
      tau24 = local(15)
      tau25 = local(16)
      tau35 = local(17)
      tau36 = local(18)
      !
      theta12 = tau14-tau24
      theta23 = tau25-tau35
      theta13 = 2.0_ark*pi-theta12-theta23
      !
      theta56 = tau36-tau35
      theta45 = tau25-tau24
      theta46 = 2.0_ark*pi-theta56-theta45
      !
      xi(14)  = ( 2.0_ark*theta12 - theta13 - theta23 )/sqrt(6.0_ark)
      xi(15)  = (                   theta13 - theta23 )/sqrt(2.0_ark)
      !
      xi(16)  = ( 2.0_ark*theta46 - theta45 - theta56 )/sqrt(6.0_ark)
      xi(17)  = (                   theta45 - theta56 )/sqrt(2.0_ark)
      !
      xi(18)  = ( tau14+tau25+tau36 )/sqrt(3.0_ark)-3.0_ark*pi/sqrt(3.0_ark)
    !
  end select


  xi(1)=1.0_ark-exp(-a*(r1-r1e))
  xi(2)=1.0_ark-exp(-b*(r2-r2e))
  xi(3)=1.0_ark-exp(-b*(r3-r2e))
  xi(4)=1.0_ark-exp(-b*(r4-r2e))
  XI(5)=1.0_ark-EXP(-B*(R5-R2E))
  XI(6)=1.0_ark-EXP(-B*(R6-R2E))
  XI(7)=1.0_ark-EXP(-B*(R7-R2E))
  !
  xi(8) = local(8)   - betae
  xi(9) = local(9)   - betae
  xi(10) = local(10) - betae
  XI(11) = LOCAL(11) - BETAE
  XI(12) = LOCAL(12) - BETAE
  XI(13) = LOCAL(13) - BETAE


  !write(out,*) xi, 'pot coords'  !BPM

  f = 0
  !
  do ioper = 1,12
    !
    call ML_symmetry_transformation_XY3_II(ioper,xi,chi(:,ioper),18)
    !
  enddo
  ! 
  do i = 6, molec%parmax

    ipower(1:18) = molec%pot_ind(1:18,i)

    term = 0 

    do ioper = 1,12
   
      term = term + product(chi(1:18,ioper)**ipower(1:18))
  
    end do

    term = term/12.0_ark

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

end function MLpoten_c2h6_88

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





recursive subroutine ML_dipole_c2h6_4m_dummy(rank,ncoords,natoms,local,xyz0,f)
    !
    implicit none
    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords),xyz0(natoms,3)
    real(ark)              ::  f(rank)
    !
    !
    stop 'dipole_c2h6_4m_dummy is not implemented'
    !
end subroutine ML_dipole_c2h6_4m_dummy



end module pot_c2h6
