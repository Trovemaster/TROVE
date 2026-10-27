! This unit defines all specific routines for a six-atomic molecule of ethane-type

module pot_c2h6
use accuracy
use moltype
use mol_c2h6

implicit none

public MLpoten_c2h6_88,ML_dipole_c2h6_4m_dummy,MLpoten_c2h6_88_cos3tau,&
       MLpoten_c2h6_88_cos3tau_142536,MLpoten_c2h6_88_cos3tau_sym,MLpoten_c2h6_Duncan
public MLpoten_c2h6_88_cos3tau_G36,ML_alpha_C2H6_zero_order,MLpoten_c2h6_88_cos3tau_sin3tau_G36,&
       MLpoten_c2h6_explct_M2_P6,MLpoten_c2h6_Morse_displace_fourier_G36,MLdms2xyz_c2h6_g36

private

integer(ik), parameter :: verbose = 4 ! Verbosity level
!
! G_36(EM) symmetry-adapted dipole moment surface of C2H6 (DMS_C2H6_G36, MLdms2xyz_c2h6_g36): numbers of
! the expansion parameters of the parallel mu_z (A_4s) and perpendicular (mu_x,mu_y) (E_2d) components;
! PARAM_DEG_Z, PARAM_DEG_XY: polynomial degree of each basis function in xi(1:17);
! PARAM_M_Z: torsional harmonic m of each mu_z function, cos/sin(m*tau);
! PARAM_HFREQ_XY: torsional harmonic K of each (mu_x,mu_y) function in units of tau/2, cos/sin(K*tau/2).
! The basis functions are defined in eval_c2h6_dms_basis at the end of this module.
!
integer, parameter :: N_PARAMS_Z  = 70
integer, parameter :: N_PARAMS_XY = 109
integer, parameter :: N_PARAMS_DMS = 179
integer, parameter :: PARAM_DEG_Z(N_PARAMS_Z) = [ &
    1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
    2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 4, 4, 4, 4, 4, 4, 4, 4, 4 ]
integer, parameter :: PARAM_M_Z(N_PARAMS_Z) = [ &
    0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, &
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
integer, parameter :: PARAM_DEG_XY(N_PARAMS_XY) = [ &
    1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 ]
integer, parameter :: PARAM_HFREQ_XY(N_PARAMS_XY) = [ &
    1, 1, 1, 5, 5, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 1, 1, &
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, &
    3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3 ]


contains

function MLpoten_c2h6_88(ncoords, natoms, local, xyz, force) result(f)

  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f

  real(ark) :: xi(18),r1,r2,r3,r4,r5,r6,r7,r1e,r2e,betae,a,b
  real(ark) :: beta1,beta2,beta3,beta4,beta5,beta6
  real(ark) :: chi(18,12),term,rad
  integer(ik) :: ioper,ipower(18),i

  real(ark) :: tau14,tau24,tau25,tau35,tau36,theta12,theta23,theta13,theta56,theta45,theta46


  rad = pi/180.0_ark

  ! expansion functions

! TRANSFORM HERE FROM Z-MATRIX COORDINATES TO SYMMETRIZED COORDS
! AS FITTING WAS CARRIED OUT USING THESE COORDINATES

  r1      = local(1)
  r2      = local(2)
  r3      = local(4)
  r4      = local(6)
  r5      = local(3)
  r6      = local(5)
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
  case('R-R16-BETA16-THETA-TAU-2')
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


      !r1=    1.53580000   
      !r2=    1.08770000   
      !r3=    1.08770000   
      !r4=    1.08770000   
      !r5=    1.08770000   
      !r6=    1.08770000   
      !r7=    1.08770000   
      
      !beta1 =  111.22000000*rad   
      !beta2 = 111.22000000*rad   
      !beta3 = 111.22000000*rad   
      !beta4 = 111.22000000*rad   
      !beta5 = 111.22000000*rad   
      !beta6 = 111.22000000*rad   
      !
      !theta45 =  120.00000000*rad   
      !theta46 =  120.00000000*rad    ! - 1
      !tau14   =  145.00000000*rad   
      !theta12 =  120.00000000*rad   
      !theta13 =  120.00000000*rad    ! -1
      !theta56 = 2.0_ark*pi-theta46-theta45
      !theta13 = 2.0_ark*pi-theta12-theta23
      !
      !tau24 = tau14-theta12 
      !tau25 = theta45+tau24
      !tau35 = tau25-theta23
      !tau36 = theta56+tau35

      !
      xi(14)  = ( 2.0_ark*theta12 - theta13 - theta23 )/sqrt(6.0_ark)
      xi(15)  = (                   theta13 - theta23 )/sqrt(2.0_ark)
      !
      xi(16)  = ( 2.0_ark*theta46 - theta45 - theta56 )/sqrt(6.0_ark)
      xi(17)  = (                   theta45 - theta56 )/sqrt(2.0_ark)
      !
      !!!xi(18)  = ( tau14+tau25+tau36 )/sqrt(3.0_ark)-3.0_ark*pi/sqrt(3.0_ark)

      xi(18)  = ( tau14+tau25+tau36 )/3.0_ark-pi


    !
  end select


  xi(1)=1.0_ark-exp(-a*(r1-r1e))
  xi(2)=1.0_ark-exp(-b*(r2-r2e))
  xi(3)=1.0_ark-exp(-b*(r3-r2e))
  xi(4)=1.0_ark-exp(-b*(r4-r2e))
  xi(5)=1.0_ark-exp(-b*(r5-r2e))
  xi(6)=1.0_ark-exp(-b*(r6-r2e))
  xi(7)=1.0_ark-exp(-b*(r7-r2e))
  !
  xi(8) = local(8)   - betae
  xi(9) = local(10)   - betae
  xi(10) = local(12) - betae
  xi(11) = local(9) - betae
  xi(12) = local(11) - betae
  xi(13) = local(13) - betae

  !xi(8)  = beta1 - betae
  !xi(9)  = beta2 - betae
  !xi(10) = beta3 - betae
  !xi(11) = beta4 - betae
  !xi(12) = beta5 - betae
  !xi(13) = beta6 - betae
  !xi(18) = xi(18)/sqrt(3.0_ark)

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
    real(ark) :: repres(12,NDEG,NDEG),a,b,e,o
    integer(ik) :: NSYM
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



function MLpoten_c2h6_88_cos3tau(ncoords, natoms, local, xyz, force) result(f)

  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f

  real(ark) :: xi(18),r1,r2,r3,r4,r5,r6,r7,r1e,r2e,betae,a,b
  real(ark) :: beta1,beta2,beta3,beta4,beta5,beta6
  real(ark) :: chi(18,12),term,rad,rhobar
  integer(ik) :: ioper,ipower(18),i

  real(ark) :: tau14,tau24,tau25,tau35,tau36,theta12,theta23,theta13,theta56,theta45,theta46


  rad = pi/180.0_ark

  ! expansion functions

! TRANSFORM HERE FROM Z-MATRIX COORDINATES TO SYMMETRIZED COORDS
! AS FITTING WAS CARRIED OUT USING THESE COORDINATES

  r1      = local(1)
  r2      = local(2)
  r3      = local(4)
  r4      = local(6)
  r5      = local(3)
  r6      = local(7)
  r7      = local(5)

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
  case('R-R16-BETA16-THETA-TAU-33')
    !
    tau14 = local(14)
    tau24 = local(15)
    tau25 = local(16)
    tau35 = local(17)
    tau36 = local(18)
    !
    theta12 = mod(tau14-tau24+2.0_ark*pi,2.0_ark*pi)
    theta23 = mod(tau25-tau35+2.0_ark*pi,2.0_ark*pi)
    theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
    !

    tau14 = local(14)
    tau24 = local(15)
    tau25 = local(16)
    tau35 = local(17)
    tau36 = local(18)
    !
    theta12 = mod(tau14-tau24+2.0_ark*pi,2.0_ark*pi)
    theta23 = mod(tau25-tau35+2.0_ark*pi,2.0_ark*pi)
    theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
    !
    theta56 = mod(tau36-tau35+2.0_ark*pi,2.0_ark*pi)
    theta45 = mod(tau25-tau24+2.0_ark*pi,2.0_ark*pi)
    theta46 = mod(2.0_ark*pi-theta56-theta45+2.0_ark*pi,2.0_ark*pi)
    !
    xi(14)  = ( 2.0_ark*theta23 - theta13 - theta12 )/sqrt(6.0_ark)
    xi(15)  = (                   theta13 - theta12 )/sqrt(2.0_ark)
    !
    xi(16)  = ( 2.0_ark*theta46 - theta45 - theta56 )/sqrt(6.0_ark)
    xi(17)  = (                   theta45 - theta56 )/sqrt(2.0_ark)
    !
    !!!xi(18)  = ( tau14+tau25+tau36 )/sqrt(3.0_ark)-3.0_ark*pi/sqrt(3.0_ark)

    rhobar = ( tau14+tau25+tau36 )/3.0_ark

    xi(18) = 1.0_ark + cos(3.0_ark*rhobar)
    !
  case('R-R16-BETA16-THETA-TAU-2')
    !
    tau14 = local(14)
    tau24 = local(15)
    tau25 = local(16)
    tau35 = local(17)
    tau36 = local(18)
    !
    theta12 = mod(tau14-tau24+2.0_ark*pi,2.0_ark*pi)
    theta23 = mod(tau25-tau35+2.0_ark*pi,2.0_ark*pi)
    theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
    !

    tau14 = local(14)
    tau24 = local(15)
    tau25 = local(16)
    tau35 = local(17)
    tau36 = local(18)
    !
    theta12 = mod(tau14-tau24+2.0_ark*pi,2.0_ark*pi)
    theta23 = mod(tau25-tau35+2.0_ark*pi,2.0_ark*pi)
    theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
    !
    theta56 = mod(tau36-tau35+2.0_ark*pi,2.0_ark*pi)
    theta45 = mod(tau25-tau24+2.0_ark*pi,2.0_ark*pi)
    theta46 = mod(2.0_ark*pi-theta56-theta45+2.0_ark*pi,2.0_ark*pi)


    !r1=    1.53580000   
    !r2=    1.08770000   
    !r3=    1.08770000   
    !r4=    1.08770000   
    !r5=    1.08770000   
    !r6=    1.08770000   
    !r7=    1.08770000   
    
    !beta1 =  111.22000000*rad   
    !beta2 = 111.22000000*rad   
    !beta3 = 111.22000000*rad   
    !beta4 = 111.22000000*rad   
    !beta5 = 111.22000000*rad   
    !beta6 = 111.22000000*rad   
    !
    !theta45 =  120.00000000*rad   
    !theta46 =  120.00000000*rad    ! - 1
    !tau14   =  145.00000000*rad   
    !theta12 =  120.00000000*rad   
    !theta13 =  120.00000000*rad    ! -1
    !theta56 = 2.0_ark*pi-theta46-theta45
    !theta13 = 2.0_ark*pi-theta12-theta23
    !
    !tau24 = tau14-theta12 
    !tau25 = theta45+tau24
    !tau35 = tau25-theta23
    !tau36 = theta56+tau35
    !
    xi(14)  = ( 2.0_ark*theta12 - theta13 - theta23 )/sqrt(6.0_ark)
    xi(15)  = (                   theta13 - theta23 )/sqrt(2.0_ark)
    !
    xi(16)  = ( 2.0_ark*theta46 - theta45 - theta56 )/sqrt(6.0_ark)
    xi(17)  = (                   theta45 - theta56 )/sqrt(2.0_ark)
    !
    !!!xi(18)  = ( tau14+tau25+tau36 )/sqrt(3.0_ark)-3.0_ark*pi/sqrt(3.0_ark)

    rhobar = ( tau14+tau25+tau36 )/3.0_ark

    xi(18) = 1.0_ark + cos(3.0_ark*rhobar)
    !
  case('R-R16-BETA16-THETA-TAU','R-R16-BETA16-THETA-TAU-4','R-R16-BETA16-THETA-TAU-5','R-R16-BETA16-THETA-TAU-6',&
       'R-R16-BETA16-THETA-TAU-7')
    !
    tau14 = local(14)
    tau24 = local(15)
    tau25 = local(16)
    tau35 = local(17)
    tau36 = local(18)
    !
    theta12 = mod(tau14-tau24+2.0_ark*pi,2.0_ark*pi)
    theta23 = mod(tau25-tau35+2.0_ark*pi,2.0_ark*pi)
    theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
    !

    tau14 = local(14)
    tau24 = local(15)
    tau25 = local(16)
    tau35 = local(17)
    tau36 = local(18)
    !
    theta12 = mod(tau14-tau24+2.0_ark*pi,2.0_ark*pi)
    theta23 = mod(tau25-tau35+2.0_ark*pi,2.0_ark*pi)
    theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
    !
    theta56 = mod(tau36-tau35+2.0_ark*pi,2.0_ark*pi)
    theta45 = mod(tau25-tau24+2.0_ark*pi,2.0_ark*pi)
    theta46 = mod(2.0_ark*pi-theta56-theta45+2.0_ark*pi,2.0_ark*pi)


    !r1=    1.53580000   
    !r2=    1.08770000   
    !r3=    1.08770000   
    !r4=    1.08770000   
    !r5=    1.08770000   
    !r6=    1.08770000   
    !r7=    1.08770000   
    
    !beta1 =  111.22000000*rad   
    !beta2 = 111.22000000*rad   
    !beta3 = 111.22000000*rad   
    !beta4 = 111.22000000*rad   
    !beta5 = 111.22000000*rad   
    !beta6 = 111.22000000*rad   
    !
    !theta45 =  120.00000000*rad   
    !theta46 =  120.00000000*rad    ! - 1
    !tau14   =  145.00000000*rad   
    !theta12 =  120.00000000*rad   
    !theta13 =  120.00000000*rad    ! -1
    !theta56 = 2.0_ark*pi-theta46-theta45
    !theta13 = 2.0_ark*pi-theta12-theta23
    !
    !tau24 = tau14-theta12 
    !tau25 = theta45+tau24
    !tau35 = tau25-theta23
    !tau36 = theta56+tau35

    !
    xi(14)  = ( 2.0_ark*theta12 - theta13 - theta23 )/sqrt(6.0_ark)
    xi(15)  = (                   theta13 - theta23 )/sqrt(2.0_ark)
    !
    xi(16)  = ( 2.0_ark*theta46 - theta45 - theta56 )/sqrt(6.0_ark)
    xi(17)  = (                   theta45 - theta56 )/sqrt(2.0_ark)
    !
    !!!xi(18)  = ( tau14+tau25+tau36 )/sqrt(3.0_ark)-3.0_ark*pi/sqrt(3.0_ark)

    rhobar = ( tau14+tau25+tau36 )/3.0_ark

    xi(18) = 1.0_ark + cos(3.0_ark*rhobar)
    !
  end select
  !
  xi(1)=1.0_ark-exp(-a*(r1-r1e))
  xi(2)=1.0_ark-exp(-b*(r2-r2e))
  xi(3)=1.0_ark-exp(-b*(r3-r2e))
  xi(4)=1.0_ark-exp(-b*(r4-r2e))
  xi(5)=1.0_ark-exp(-b*(r5-r2e))
  xi(6)=1.0_ark-exp(-b*(r6-r2e))
  xi(7)=1.0_ark-exp(-b*(r7-r2e))
  !
  xi(8) = local(8)   - betae
  xi(9) = local(10)   - betae
  xi(10) = local(12) - betae
  xi(11) = local(9) - betae
  xi(12) = local(13) - betae
  xi(13) = local(11) - betae

  !xi(8)  = beta1 - betae
  !xi(9)  = beta2 - betae
  !xi(10) = beta3 - betae
  !xi(11) = beta4 - betae
  !xi(12) = beta5 - betae
  !xi(13) = beta6 - betae
  !xi(18) = xi(18)/sqrt(3.0_ark)

  !write(out,*) xi, 'pot coords'  !BPM

  f = 0
  !
  do ioper = 1,12
    !
    call ML_symmetry_transformation_XY3_III(ioper,xi,chi(:,ioper),18)
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
  !    force(8)*sum(xi(2:7)**2)+&
  !    force(9)*sum(xi(8:13)**2)+force(10)*sum(xi(14:17)**2)
  !    !force(11)*( xi(14)*xi(16)+xi(15)*xi(17) )
  !
  !f = f +sum(force(11:16)*xi(18)**molec%pot_ind(18,11:16))
  !
  f = f/12.0_ark
  !  

end function MLpoten_c2h6_88_cos3tau


subroutine ML_symmetry_transformation_XY3_III(ioper,src,dst,NDEG)
!       THIS COULD ALSO BE MADE AS GENERAL AS POSSIBLE BUT WILL ALWAYS NEED TO
!       DEFINE MATRICES FOR EACH MOLECULE BPM
    implicit none
    !
    integer,intent(in)    :: ioper,NDEG  ! group operation  
    real(ark),intent(in)      :: src(1:NDEG)
    real(ark),intent(out)     :: dst(1:NDEG)
    !
    real(ark)                 :: repres(12,NDEG,NDEG),a,b,e,o
    INTEGER NSYM
    !
    if (verbose>=5) write(out,"('ML_symmetry_transformation_XY3/start')")
    !
    a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
    Nsym = 12
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
    repres(1,10,10) = 1.0_ark
    repres(1,11,11) = 1.0_ark
    repres(1,12,12) = 1.0_ark
    repres(1,13,13) = 1.0_ark
    repres(1,14,14) = 1.0_ark
    repres(1,15,15) = 1.0_ark
    repres(1,16,16) = 1.0_ark
    repres(1,17,17) = 1.0_ark
    repres(1,18,18)  = 1.0_ark
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
    repres(2,18,18)  = 1.0_ark
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
    repres(3,18,18)  = 1.0_ark
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
    repres(4,18,18)  = 1.0_ark
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
    repres(5,18,18)  = 1.0_ark
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
    repres(6,15,16) = b
    repres(6,14,17) = b
    repres(6,15,17) = a
    repres(6,18,18)  = 1.0_ark
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
    repres(7,18,18)  = 1.0_ark
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
    repres(8,18,18)  = 1.0_ark
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
    repres(9,18,18)  = 1.0_ark
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
    repres(10,18,18)  = 1.0_ark
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
    repres(11,18,18)  = 1.0_ark
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
    repres(12,18,18)  = 1.0_ark
    !
    if (ioper<0.or.ioper>NSYM) then
      write (6,"('symmetry_transformation_local: operation ',i8,' unknown')")
      stop 'symmetry_transformation_local - bad operation. type'
    endif
    !
    ! BPM : REVERSED MULTIPLICATION ORDER TO AGREE WITH MY NOTES
    dst = matmul(src,repres(ioper,:,:))

  end subroutine ML_symmetry_transformation_XY3_III



function MLpoten_c2h6_88_cos3tau_sym(ncoords, natoms, local, xyz, force) result(f)

  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f

  real(ark) :: xi(18),r1,r2,r3,r4,r5,r6,r7,r1e,r2e,betae,a,b
  real(ark) :: beta1,beta2,beta3,beta4,beta5,beta6
  real(ark) :: chi(18,12),term,rad,rhobar
  integer(ik) :: ioper,ipower(18),i

  real(ark) :: tau14,tau24,tau25,tau35,tau36,theta12,theta23,theta13,theta56,theta45,theta46,xi_A,xi_B,xi_C,xi_D


  rad = pi/180.0_ark

  ! expansion functions

! TRANSFORM HERE FROM Z-MATRIX COORDINATES TO SYMMETRIZED COORDS
! AS FITTING WAS CARRIED OUT USING THESE COORDINATES

  r1      = local(1)
  r2      = local(2)
  r3      = local(4)
  r4      = local(6)
  r5      = local(3)
  r6      = local(5)
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
  case('R-R16-BETA16-THETA-TAU-4','R-R16-BETA16-THETA-TAU-5','R-R16-BETA16-THETA-TAU-6',&
       'R-R16-BETA16-THETA-TAU-7','R-R16-BETA16-THETA-TAU-8','R-R16-BETA16-THETA-TAU-9')
      !
      tau14 = local(14)
      tau24 = local(15)
      tau25 = local(16)
      tau35 = local(17)
      tau36 = local(18)
      !
      theta12 = mod(tau14-tau24+2.0_ark*pi,2.0_ark*pi)
      theta23 = mod(tau25-tau35+2.0_ark*pi,2.0_ark*pi)
      theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
      !

      tau14 = local(14)
      tau24 = local(15)
      tau25 = local(16)
      tau35 = local(17)
      tau36 = local(18)
      !
      theta12 = mod(tau14-tau24+2.0_ark*pi,2.0_ark*pi)
      theta23 = mod(tau25-tau35+2.0_ark*pi,2.0_ark*pi)
      theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
      !
      theta56 = mod(tau36-tau35+2.0_ark*pi,2.0_ark*pi)
      theta45 = mod(tau25-tau24+2.0_ark*pi,2.0_ark*pi)
      theta46 = mod(2.0_ark*pi-theta56-theta45+2.0_ark*pi,2.0_ark*pi)

      xi_A  = ( 2.0_ark*theta23 - theta13 - theta12 )/sqrt(12.0_ark)
      xi_B  = (                   theta13 - theta12 )/(2.0_ark)
      !
      xi_C  = ( 2.0_ark*theta56 - theta46 - theta45 )/sqrt(12.0_ark)
      xi_D  = (                   theta46 - theta45 )/(2.0_ark)
      !
      xi(14) = xi_A + xi_C 
      xi(15) = xi_B + xi_D
      xi(16) = xi_A - xi_C 
      xi(17) = xi_B - xi_D 
      !
      rhobar = ( tau14+tau25+tau36 )/3.0_ark
      !
      xi(8)  = local(8)  - betae
      xi(9)  = local(10) - betae
      xi(10) = local(12) - betae
      xi(11) = local(9)  - betae
      xi(12) = local(11) - betae
      xi(13) = local(13) - betae
      !
      xi(18) = 1.0_ark + cos(3.0_ark*rhobar)
      !
  case('R-R16-BETA16-THETA-TAU-1','R-R16-BETA16-THETA-TAU-2','R-R16-BETA16-THETA-TAU-3',&
       'R-R16-BETA16-THETA-TAU-11','R-R16-BETA16-THETA-TAU-12','R-R16-BETA16-THETA-TAU-17')
      !
      r1 = local(1)
      r2 = local(2)
      r3 = local(4)
      r4 = local(6)
      r5 = local(3)
      r6 = local(7)
      r7 = local(5)
      !
      xi(8)  = local(8)  - betae
      xi(9)  = local(10) - betae
      xi(10) = local(12) - betae
      xi(11) = local(9)  - betae
      xi(12) = local(13) - betae
      xi(13) = local(11) - betae
      !
      tau14 = mod(local(14)+4.0_ark*pi,4.0_ark*pi)
      tau24 = mod(local(15)+2.0_ark*pi,2.0_ark*pi)
      tau25 = mod(local(16)+2.0_ark*pi,2.0_ark*pi)
      tau35 = mod(local(17)+2.0_ark*pi,2.0_ark*pi)
      tau36 = mod(local(18)+2.0_ark*pi,2.0_ark*pi)
      !
      if (tau14>2.0_ark*pi) then 
         tau25 = tau25 + 2.0_ark*pi
         tau36 = tau36 + 2.0_ark*pi
      endif
      !
      rhobar  = ( tau14+tau25+tau36 )/(3.0_ark)
      !
      tau14 = mod(local(14)+2.0_ark*pi,2.0_ark*pi)
      tau24 = mod(local(15)+2.0_ark*pi,2.0_ark*pi)
      tau25 = mod(local(16)+2.0_ark*pi,2.0_ark*pi)
      tau35 = mod(local(17)+2.0_ark*pi,2.0_ark*pi)
      tau36 = mod(local(18)+2.0_ark*pi,2.0_ark*pi)
      !
      theta12 = mod(tau14-tau24+2.0_ark*pi,2.0_ark*pi)
      theta23 = mod(tau25-tau35+2.0_ark*pi,2.0_ark*pi)
      theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
      !
      theta56 = mod(tau36-tau35+2.0_ark*pi,2.0_ark*pi)
      theta45 = mod(tau25-tau24+2.0_ark*pi,2.0_ark*pi)
      theta46 = mod(2.0_ark*pi-theta56-theta45+2.0_ark*pi,2.0_ark*pi)
      !
      xi_A  = ( 2.0_ark*theta23 - theta13 - theta12 )/sqrt(6.0_ark)
      xi_B  = (                   theta13 - theta12 )/sqrt(2.0_ark)
      !
      xi_C  = ( 2.0_ark*theta56 - theta45 - theta46 )/sqrt(6.0_ark)
      xi_D  = (                   theta45 - theta46 )/sqrt(2.0_ark)
      !
      xi(14) = xi_A + xi_C 
      xi(15) = xi_B + xi_D
      xi(16) =-xi_A + xi_C 
      xi(17) =-xi_B + xi_D       
      !
      rhobar = ( tau14+tau25+tau36 )/3.0_ark
      !
      xi(18) = 1.0_ark + cos(3.0_ark*rhobar)
      !
  case('R-R16-BETA16-THETA-TAU')
      !
      tau14 = local(14)
      tau24 = local(15)
      tau25 = local(16)
      tau35 = local(17)
      tau36 = local(18)
      !
      tau14 = mod(tau14+4.0_ark*pi,4.0_ark*pi)
      tau24 = mod(tau24+4.0_ark*pi,4.0_ark*pi)
      tau25 = mod(tau25+4.0_ark*pi,4.0_ark*pi)
      tau35 = mod(tau35+4.0_ark*pi,4.0_ark*pi)
      tau36 = mod(tau36+4.0_ark*pi,4.0_ark*pi)
      !
      rhobar  = ( tau14+tau25+tau36 )/(3.0_ark)-pi
      !
      tau14 = mod(tau14+2.0_ark*pi,2.0_ark*pi)
      tau24 = mod(tau24+2.0_ark*pi,2.0_ark*pi)
      tau25 = mod(tau25+2.0_ark*pi,2.0_ark*pi)
      tau35 = mod(tau35+2.0_ark*pi,2.0_ark*pi)
      tau36 = mod(tau36+2.0_ark*pi,2.0_ark*pi)
      !
      theta12 = mod(tau14-tau24+2.0_ark*pi,2.0_ark*pi)
      theta23 = mod(tau25-tau35+2.0_ark*pi,2.0_ark*pi)
      theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
      !
      theta56 = mod(tau36-tau35+2.0_ark*pi,2.0_ark*pi)
      theta45 = mod(tau25-tau24+2.0_ark*pi,2.0_ark*pi)
      theta46 = mod(2.0_ark*pi-theta56-theta45+2.0_ark*pi,2.0_ark*pi)
      !
      xi_A  = ( 2.0_ark*theta23 - theta13 - theta12 )/sqrt(12.0_ark)
      xi_B  = (                   theta13 - theta12 )/(2.0_ark)
      !
      xi_C  = ( 2.0_ark*theta56 - theta46 - theta45 )/sqrt(12.0_ark)
      xi_D  = (                   theta46 - theta45 )/(2.0_ark)
      !
      xi(14) = xi_A + xi_C 
      xi(15) = xi_B + xi_D
      xi(16) = xi_A - xi_C 
      xi(17) = xi_B - xi_D 
      !
      !rhobar = ( tau14+tau25+tau36 )/3.0_ark
      !
      xi(18) = 1.0_ark + cos(3.0_ark*rhobar)

      xi(8)  = local(8)  - betae
      xi(9)  = local(10) - betae
      xi(10) = local(12) - betae
      xi(11) = local(9)  - betae
      xi(12) = local(11) - betae
      xi(13) = local(13) - betae


      !
  end select
  !
  xi(1)=1.0_ark-exp(-a*(r1-r1e))
  xi(2)=1.0_ark-exp(-b*(r2-r2e))
  xi(3)=1.0_ark-exp(-b*(r3-r2e))
  xi(4)=1.0_ark-exp(-b*(r4-r2e))
  xi(5)=1.0_ark-exp(-b*(r5-r2e))
  xi(6)=1.0_ark-exp(-b*(r6-r2e))
  xi(7)=1.0_ark-exp(-b*(r7-r2e))
  !
  !
  f = force(6)*xi(1)**2+force(7)*sum(xi(2:7)**2)/6.0_ark+force(8)*sum(xi(8:13)**2)/6.0_ark+&
      force(9)*(xi(14)**2+xi(15)**2)/2.0_ark+force(10)*(xi(16)**2+xi(17)**2)/2.0_ark+&
      force(11)*xi(18)+force(12)*xi(18)**2+force(13)*xi(18)**3+force(14)*xi(18)**4+force(15)*xi(18)**5+force(16)*xi(18)**6
   !  

end function MLpoten_c2h6_88_cos3tau_sym



function MLpoten_c2h6_Duncan(ncoords, natoms, local, xyz, force) result(f)

  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f

  real(ark) :: xi(18),r1,r2,r3,r4,r5,r6,r7,r1e,r2e,betae,a,b
  real(ark) :: beta1,beta2,beta3,beta4,beta5,beta6
  real(ark) :: chi(18,12),term,rad,rhobar
  integer(ik) :: ioper,ipower(18),i

  real(ark) :: tau14,tau24,tau25,tau35,tau36,theta12,theta23,theta13,theta56,theta45,theta46,xi_A,xi_B,xi_C,xi_D,alphae
  real(ark) :: F11,F12,F13,F22,F23,F33,F44,F55,F56,F66,F77,F78,F79,F88,F89,F99,F1010,F1011,F1012,F1111,F1112,F1212
  real(ark) :: kappa,S1,S2,S3,S4,S5,S6,S7a,S7b,S8a,S8b,S9a,S9b,S10a,S10b,S11a,S11b,S12a,S12b,S2a,S2b,S6a,S6b
  real(ark) :: c12,c13,c23,c45,c46,c56,thetae,cosae,alpha12,alpha13,alpha23,alpha45,alpha46,alpha56

  rad = pi/180.0_ark

  ! expansion functions

! TRANSFORM HERE FROM Z-MATRIX COORDINATES TO SYMMETRIZED COORDS
! AS FITTING WAS CARRIED OUT USING THESE COORDINATES

  r1      = local(1)
  r2      = local(2)
  r3      = local(4)
  r4      = local(6)
  r5      = local(3)
  r6      = local(5)
  r7      = local(7)

  r1e = force(1)
  r2e = force(2)
  thetae = 2.0_ark*pi/3.0_ark
  betae = force(3)*rad
  !
  cosae = cos(betae)*cos(betae)+cos(thetae)*sin(betae)*sin(betae)
  !
  alphae = acos(cosae)
  !
  a = force(4)
  b = force(5)
  !
  F11 = force( 6)
  F12 = force( 7)
  F13 = force( 8)
  F22 = force( 9)
  F23 = force(10)
  F33 = force(11)
  F44 = force(12)
  F55 = force(13)
  F56 = force(14)
  F66 = force(15)
  F77 = force(16)
  F78 = force(17)
  F79 = force(18)
  F88 = force(19)
  F89 = force(20)
  F99 = force(21)
  F1010 = force(22)
  F1011 = force(23)
  F1012 = force(24)
  F1111 = force(25)
  F1112 = force(26)
  F1212 = force(27)
  !
  select case(trim(molec%coords_transform))
    !
  case default
    !
    write(out, '(/a,1x,a,1x,a)') &
    'MLpoten_c2h6_88 error', trim(molec%coords_transform), 'is unknown'
    stop 'MLpoten_c2h6_88 error error: bad coordinate type'
    !
  case('R-R16-BETA16-THETA-TAU')
      !
      tau14 = local(14)
      tau24 = local(15)
      tau25 = local(16)
      tau35 = local(17)
      tau36 = local(18)
      !
      theta12 = mod(tau14-tau24+2.0_ark*pi,2.0_ark*pi)
      theta23 = mod(tau25-tau35+2.0_ark*pi,2.0_ark*pi)
      theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
      !
      theta56 = mod(tau36-tau35+2.0_ark*pi,2.0_ark*pi)
      theta45 = mod(tau25-tau24+2.0_ark*pi,2.0_ark*pi)
      theta46 = mod(2.0_ark*pi-theta56-theta45+2.0_ark*pi,2.0_ark*pi)
      !
      beta1 = local(8)
      beta2 = local(10)
      beta3 = local(12)
      beta4 = local(9)
      beta5 = local(11)
      beta6 = local(13)
      !
      c12 = cos(beta1)*cos(beta2)+cos(theta12)*sin(beta1)*sin(beta2)
      c13 = cos(beta1)*cos(beta3)+cos(theta13)*sin(beta1)*sin(beta3)
      c23 = cos(beta2)*cos(beta3)+cos(theta23)*sin(beta2)*sin(beta3)
      !
      c45 = cos(beta4)*cos(beta5)+cos(theta45)*sin(beta4)*sin(beta5)
      c46 = cos(beta4)*cos(beta6)+cos(theta46)*sin(beta4)*sin(beta6)
      c56 = cos(beta5)*cos(beta6)+cos(theta56)*sin(beta5)*sin(beta6)
      !  
      alpha12 = acos(c12)
      alpha13 = acos(c13)
      alpha23 = acos(c23)
      !              
      alpha45 = acos(c45)
      alpha46 = acos(c46)
      alpha56 = acos(c56)
      !
      xi_A  = ( 2.0_ark*alpha23 - alpha13 - alpha12 )
      xi_B  = (                   alpha13 - alpha12 )
      !
      xi_C  = ( 2.0_ark*alpha56 - alpha46 - alpha45 )
      xi_D  = (                   alpha46 - alpha45 )
      !
      xi(14) = xi_A + xi_C 
      xi(15) = xi_B + xi_D
      xi(16) = xi_A - xi_C 
      xi(17) = xi_B - xi_D 

      rhobar = ( tau14+tau25+tau36 )/3.0_ark

      xi(18) = 1.0_ark + cos(3.0_ark*rhobar)
      !
  end select
  !
  xi(1)=r1-r1e
  xi(2)=r2-r2e
  xi(3)=r3-r2e
  xi(4)=r4-r2e
  xi(5)=r5-r2e
  xi(6)=r6-r2e
  xi(7)=r7-r2e
  !
  xi(8)  = local(8)  - betae
  xi(9)  = local(10) - betae
  xi(10) = local(12) - betae
  xi(11) = local(9)  - betae
  xi(12) = local(11) - betae
  xi(13) = local(13) - betae

  S1 = 1.0_ark/sqrt(6.0_ark)*(sum(xi(2:7)))
  !
  kappa = 3.0_ark*sin(betae)*cos(betae)/sin(alphae)
  !
  S2a = (alpha23 + alpha13 + alpha12 + alpha56 + alpha46 + alpha45)-6.0_ark*alphae
  S2b = sum(xi(8:13))
  !
  S2 = -(kappa*S2a+S2b)/sqrt(6.0_ark*(1.0_ark+kappa**2))
  !
  S3 = xi(1)
  !
  S4 = rhobar-pi
  !
  S5 = ( r2+r3+r4-(r5+r6+r7) )/sqrt(6.0_ark)
  !
  S6a = (alpha23 + alpha13 + alpha12 - ( alpha56 + alpha46 + alpha45) )
  S6b = ( beta1+beta2+beta3-(beta4+beta5+beta6) )
  !
  S6 = -(kappa*S6a+S6b)/sqrt( 6.0_ark*( 1.0_ark+kappa**2) )
  !
  S7a = ( 2.0_ark*r2-r3-r4-(2.0_ark*r5-r6-r7) )/sqrt(12.0_ark)
  S7b = (            r3-r4-(           r6-r7) )/sqrt(4.0_ark)
  !
  S8a = ( xi_A-xi_C )/sqrt(12.0_ark)
  S8b = ( xi_B-xi_D )/sqrt(4.0_ark)
  !
  S9a = ( 2.0_ark*xi(8)-xi(9)-xi(10)-(2.0_ark*xi(11)-xi(12)-xi(13)) )/sqrt(12.0_ark)
  S9b = (               xi(9)-xi(10)-(               xi(12)-xi(13)) )/sqrt(4.0_ark)
  !
  S10a = ( 2.0_ark*r2-r3-r4+(2.0_ark*r5-r6-r7) )/sqrt(12.0_ark)
  S10b = (            r3-r4+(           r6-r7) )/sqrt(4.0_ark)
  !
  S11a = ( xi_A+xi_C )/sqrt(12.0_ark)
  S11b = ( xi_B+xi_D )/sqrt(4.0_ark)
  !
  S12a = ( 2.0_ark*xi(8)-xi(9)-xi(10)+(2.0_ark*xi(11)-xi(12)-xi(13)) )/sqrt(12.0_ark)
  S12b = (               xi(9)-xi(10)+(               xi(12)-xi(13)) )/sqrt(4.0_ark)
  !
  f = F11*S1**2+2.0_ark*(F12*S1*S2+F13*S1*S3)+&
      F22*S2**2+2.0_ark*F23*S2*S3+&
      F33*S3**2+&
      F44*S4**2+&
      F55*S5**2+2.0_ark*F56*S5*S6+&
      F66*S6**2+&
      F77*(S7a**2+S7b**2)+F88*(S8a**2+S8b**2)+F99*(S9a**2+S9b**2)+&
      2.0_ark*( F78*(S7a*S8a+S7b*S8b )+F79*(S7a*S9a+S7b*S9b)+F89*( S8a*S9a+S8b*S9b ) )+&
      F1010*( S10a**2+S10b**2 )+F1111*( S11a**2+S11b**2 )+F1212*( S12a**2+S12b**2 )+&
      2.0_ark*( F1112*( S11a*S12a+S11b*S12b )+F1011*( S10a*S11a+S10b*S11b )+F1012*( S10a*S12a+S10b*S12b ) )
      !  
  f = f/2.0_ark*1.0e-11/planck/vellgt    
  !

end function MLpoten_c2h6_Duncan


function MLpoten_c2h6_88_cos3tau_142536(ncoords, natoms, local, xyz, force) result(f)

  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f

  real(ark) :: xi(18),r1,r2,r3,r4,r5,r6,r7,r1e,r2e,betae,a,b
  real(ark) :: beta1,beta2,beta3,beta4,beta5,beta6
  real(ark) :: chi(18,12),term,rad,rhobar
  integer(ik) :: ioper,ipower(18),i

  real(ark) :: tau14,tau24,tau25,tau35,tau36,theta12,theta23,theta13,theta56,theta45,theta46


  rad = pi/180.0_ark

  ! expansion functions

! TRANSFORM HERE FROM Z-MATRIX COORDINATES TO SYMMETRIZED COORDS
! AS FITTING WAS CARRIED OUT USING THESE COORDINATES

  r1      = local(1)
  r2      = local(2)
  r3      = local(4)
  r4      = local(6)
  r5      = local(3)
  r6      = local(5)
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
      theta12 = mod(tau14-tau24+2.0_ark*pi,2.0_ark*pi)
      theta23 = mod(tau25-tau35+2.0_ark*pi,2.0_ark*pi)
      theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
      !

      tau14 = local(14)
      tau24 = local(15)
      tau25 = local(16)
      tau35 = local(17)
      tau36 = local(18)
      !
      theta12 = mod(tau14-tau24+2.0_ark*pi,2.0_ark*pi)
      theta23 = mod(tau25-tau35+2.0_ark*pi,2.0_ark*pi)
      theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
      !
      theta56 = mod(tau36-tau35+2.0_ark*pi,2.0_ark*pi)
      theta45 = mod(tau25-tau24+2.0_ark*pi,2.0_ark*pi)
      theta46 = mod(2.0_ark*pi-theta56-theta45+2.0_ark*pi,2.0_ark*pi)


      !r1=    1.53580000   
      !r2=    1.08770000   
      !r3=    1.08770000   
      !r4=    1.08770000   
      !r5=    1.08770000   
      !r6=    1.08770000   
      !r7=    1.08770000   
      
      !beta1 =  111.22000000*rad   
      !beta2 = 111.22000000*rad   
      !beta3 = 111.22000000*rad   
      !beta4 = 111.22000000*rad   
      !beta5 = 111.22000000*rad   
      !beta6 = 111.22000000*rad   
      !
      !theta45 =  120.00000000*rad   
      !theta46 =  120.00000000*rad    ! - 1
      !tau14   =  145.00000000*rad   
      !theta12 =  120.00000000*rad   
      !theta13 =  120.00000000*rad    ! -1
      !theta56 = 2.0_ark*pi-theta46-theta45
      !theta13 = 2.0_ark*pi-theta12-theta23
      !
      !tau24 = tau14-theta12 
      !tau25 = theta45+tau24
      !tau35 = tau25-theta23
      !tau36 = theta56+tau35

      !
      xi(14)  = ( 2.0_ark*theta23 - theta13 - theta12 )/sqrt(6.0_ark)
      xi(15)  = (                   theta13 - theta12 )/sqrt(2.0_ark)
      !
      xi(16)  = ( 2.0_ark*theta56 - theta46 - theta45 )/sqrt(6.0_ark)
      xi(17)  = (                   theta46 - theta45 )/sqrt(2.0_ark)
      !
      !!!xi(18)  = ( tau14+tau25+tau36 )/sqrt(3.0_ark)-3.0_ark*pi/sqrt(3.0_ark)

      rhobar = ( tau14+tau25+tau36 )/3.0_ark

      xi(18) = 1.0_ark + cos(3.0_ark*rhobar)
      !
  end select
  !
  xi(1)=1.0_ark-exp(-a*(r1-r1e))
  xi(2)=1.0_ark-exp(-b*(r2-r2e))
  xi(3)=1.0_ark-exp(-b*(r3-r2e))
  xi(4)=1.0_ark-exp(-b*(r4-r2e))
  xi(5)=1.0_ark-exp(-b*(r5-r2e))
  xi(6)=1.0_ark-exp(-b*(r6-r2e))
  xi(7)=1.0_ark-exp(-b*(r7-r2e))
  !
  xi(8)  = local(8)   - betae
  xi(9)  = local(10)   - betae
  xi(10) = local(12) - betae
  xi(11) = local(9) - betae
  xi(12) = local(11) - betae
  xi(13) = local(13) - betae

  !xi(8)  = beta1 - betae
  !xi(9)  = beta2 - betae
  !xi(10) = beta3 - betae
  !xi(11) = beta4 - betae
  !xi(12) = beta5 - betae
  !xi(13) = beta6 - betae
  !xi(18) = xi(18)/sqrt(3.0_ark)

  !write(out,*) xi, 'pot coords'  !BPM

  f = 0
  !
  do ioper = 1,12
    !
    !call ML_symmetry_transformation_XY3_IV(ioper,xi,chi(:,ioper),18)
    !
    call ML_symmetry_transformation_C2H6(ioper, 18, xi, chi(:,ioper))
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

end function MLpoten_c2h6_88_cos3tau_142536


subroutine ML_symmetry_transformation_XY3_IV(ioper,src,dst,NDEG)
!       THIS COULD ALSO BE MADE AS GENERAL AS POSSIBLE BUT WILL ALWAYS NEED TO
!       DEFINE MATRICES FOR EACH MOLECULE BPM
    implicit none
    !
    integer,intent(in)    :: ioper,NDEG  ! group operation  
    real(ark),intent(in)      :: src(1:NDEG)
    real(ark),intent(out)     :: dst(1:NDEG)
    !
    real(ark)                 :: repres(12,NDEG,NDEG),a,b,e,o
    INTEGER NSYM
    !
    if (verbose>=5) write(out,"('ML_symmetry_transformation_XY3/start')")
    !
    a = 0.5_ark ; b = 0.5_ark*sqrt(3.0_ark) ; e = 1.0_ark ; o = 0.0_ark
    Nsym = 12
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
    repres(1,10,10) = 1.0_ark
    repres(1,11,11) = 1.0_ark
    repres(1,12,12) = 1.0_ark
    repres(1,13,13) = 1.0_ark
    repres(1,14,14) = 1.0_ark
    repres(1,15,15) = 1.0_ark
    repres(1,16,16) = 1.0_ark
    repres(1,17,17) = 1.0_ark
    repres(1,18,18)  = 1.0_ark
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
    repres(2,18,18)  = 1.0_ark
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
    repres(3,18,18)  = 1.0_ark
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
    repres(4,18,18)  = 1.0_ark
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
    repres(5,18,18)  = 1.0_ark
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
    repres(6,15,16) = b
    repres(6,14,17) = b
    repres(6,15,17) = a
    repres(6,18,18)  = 1.0_ark
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
    repres(7,18,18)  = 1.0_ark
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
    repres(8,18,18)  = 1.0_ark
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
    repres(9,18,18)  = 1.0_ark
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
    repres(10,18,18)  = 1.0_ark
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
    repres(11,18,18)  = 1.0_ark
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
    repres(12,18,18)  = 1.0_ark
    !
    if (ioper<0.or.ioper>NSYM) then
      write (6,"('symmetry_transformation_local: operation ',i8,' unknown')")
      stop 'symmetry_transformation_local - bad operation. type'
    endif
    !
    ! BPM : REVERSED MULTIPLICATION ORDER TO AGREE WITH MY NOTES
    dst = matmul(src,repres(ioper,:,:))

  end subroutine ML_symmetry_transformation_XY3_IV



   recursive subroutine ML_symmetry_transformation_C2H6_G36(ioper, nmodes, src, dst, s18, d18)
    !
    ! Symmetry transformation rules of coordinates
    !
    integer(ik), intent(in)  ::  ioper
    integer(ik), intent(in)  ::  nmodes
    real(ark), intent(in)    ::  src(1:nmodes)
    real(ark), intent(out)   ::  dst(1:nmodes)
    real(ark), optional,intent(in)    ::  s18
    real(ark), optional,intent(out)   ::  d18
    !
    real(ark) :: a,b,e,o,g(1:4,1:4),c123(2,2),c132(2,2),a123(3,3),a132(3,3),sxy(2,2),i(3,3),i2(3,3)
    real(ark) :: s18_,d18_
    !
    real(ark),dimension(size(src)) :: tmp
    !
    integer(ik)  :: tn(72,2), temp(144)
    integer(ik) :: nsrc
    !
    temp(1:36)   = (/0, 0, 2, 0, 6, 4, 0, 7, 2, 3, 2, 3, 4, 5, 6, 4, 5, 6, 0,21,19, 2, 3, 2, 3, 2, 3, 4, 5, 6, 4, 5, 6, 4, 5, 6/)
    temp(73:108) = (/0, 0, 2, 0, 2, 2, 0, 7, 7, 7, 8, 8, 7, 7, 7, 8, 8, 8, 0, 7, 7,19,19,20,20,21,21,19,19,19,20,20,20,21,21,21/)
    !
    tn = reshape( temp, (/ 72, 2/))
    !
    s18_ = 0 ; d18_ = 0
    !
    if (present(s18)) then 
       s18_ = s18
    endif 
    !
    a = 0.5_ark
    b = 0.5_ark*sqrt(3.0_ark)
    !
    e = 1.0_ark
    o = 0.0_ark
    !!
    c123 = transpose(reshape( (/ -a, -b, &
                                  b, -a/), (/ 2, 2/)))
    c132  = matmul(c123,c123)
    !
    a123 = transpose(reshape( (/ o, o, e,&
                                 e, o, o,&
                                 o, e, o/), (/ 3, 3/)))
    a132  = matmul(a123,a123)
    !
    i = transpose(reshape( (/ e, o, o,&
                              o, e, o,&
                              o, o, e/), (/ 3, 3/)))
    !
    i2= transpose(reshape( (/ e, o, o,&
                              o, o, e,&
                              o, e, o/), (/ 3, 3/)))
    !
    sxy = transpose(reshape( (/ e,  o, &
                                o, -e /), (/ 2, 2/)))
    !
    select case(ioper)
      !
    case (1) ! E
      !
      dst(1:18) = src(1:18)
      !
    case (2) !C(+)/(123)(456)
      !
      dst(1) = src(1)
      !
      dst(2:4) = matmul(a123,src(2:4))
      dst(5:7) = matmul(a123,src(5:7))
      !
      dst(8:10) = matmul(a123,src(8:10))
      dst(11:13) = matmul(a123,src(11:13))
      !!
      dst(14:15) = matmul(c123,src(14:15))
      dst(16:17) = matmul(c123,src(16:17))
      !
      !!!! remember that 
      !dst(18) = src(18)  - 4.0_ark/3.0_ark*pi
      !!!! src(18) = 1 + cos(3 tau)
      !
      dst(18) = src(18)
      !    
    case (4) !sxy(+)/(14)(26)(35)(ab)* 
      !
      dst(1) = src(1)
      !
      dst(2:4) = matmul(i2,src(5:7))
      dst(5:7) = matmul(i2,src(2:4))
      !
      dst(8:10) = matmul(i2,src(11:13))
      dst(11:13) = matmul(i2,src(8:10))
      !
      !!
      dst(14:15) = matmul(sxy,src(16:17))
      dst(16:17) = matmul(sxy,src(14:15))
      !
      !!!!!
      !dst(18) =  4.0_ark*pi - src(18)
      !!!!!
      dst(18) = src(18)
      !
      d18_ = -s18_
      !
    case (7) ! C(-)/(132)(456)
      !
      dst(1) = src(1)
      !
      dst(2:4) = matmul(a132,src(2:4))
      dst(5:7) = matmul(a123,src(5:7))
      !
      dst(8:10) = matmul(a132,src(8:10))
      dst(11:13) = matmul(a123,src(11:13))
      !!
      dst(14:15) = matmul(c132,src(14:15))
      dst(16:17) = matmul(c123,src(16:17))
      !
      !!
      dst(18) = src(18)
      !
    case (19) !sxy(-)/(14)(25)(36)(ab)
      !
      dst(1) = src(1)
      !
      dst(2:4) = matmul(i,src(5:7))
      dst(5:7) = matmul(i,src(2:4))
      !
      dst(8:10) = matmul(i,src(11:13))
      dst(11:13) = matmul(i,src(8:10))
      !!
      dst(14:15) = src(16:17)
      dst(16:17) = src(14:15)
      !!
      dst(18) = src(18)
      !
    case(37) !E'
      !
      dst(1:18) = src(1:18)
      !dst(18) = src(18) + 2.0_ark*pi
      !
    end select
    !
    if (all(tn(ioper,:)/=0)) then
        call ML_symmetry_transformation_C2H6_G36(tn(ioper,1),nmodes,src,tmp,s18_,d18_)
        call ML_symmetry_transformation_C2H6_G36(tn(ioper,2),nmodes,tmp,dst,s18_,d18_)
    endif 
    !
    if (present(d18)) then 
       d18 = d18_
    endif 
    !
  end subroutine ML_symmetry_transformation_C2H6_G36


   recursive subroutine ML_symmetry_transformation_C2H6_G36_full(ioper, nmodes, src, dst)
    !
    ! Symmetry transformation rules of coordinates
    !
    integer(ik), intent(in)  ::  ioper
    integer(ik), intent(in)  ::  nmodes
    real(ark), intent(in)    ::  src(1:nmodes)
    real(ark), intent(out)   ::  dst(1:nmodes)
    !
    real(ark) :: a,b,e,o,g(1:4,1:4),c123(2,2),c132(2,2),a123(3,3),a132(3,3),sxy(2,2),i(3,3),i2(3,3)
    !
    real(ark),dimension(size(src)) :: tmp
    !
    integer(ik)  :: tn(72,2), temp(144)
    integer(ik) :: nsrc
    !
    temp(1:36)   = (/0, 0, 2, 0, 6, 4, 0, 7, 2, 3, 2, 3, 4, 5, 6, 4, 5, 6, 0,21,19, 2, 3, 2, 3, 2, 3, 4, 5, 6, 4, 5, 6, 4, 5, 6/)
    temp(73:108) = (/0, 0, 2, 0, 2, 2, 0, 7, 7, 7, 8, 8, 7, 7, 7, 8, 8, 8, 0, 7, 7,19,19,20,20,21,21,19,19,19,20,20,20,21,21,21/)
    !
    tn = reshape( temp, (/ 72, 2/))
    !
    a = 0.5_ark
    b = 0.5_ark*sqrt(3.0_ark)
    !
    e = 1.0_ark
    o = 0.0_ark
    !!
    c123 = transpose(reshape( (/ -a, -b, &
                                  b, -a/), (/ 2, 2/)))
    c132  = matmul(c123,c123)
    !
    a123 = transpose(reshape( (/ o, o, e,&
                                 e, o, o,&
                                 o, e, o/), (/ 3, 3/)))
    a132  = matmul(a123,a123)
    !
    i = transpose(reshape( (/ e, o, o,&
                              o, e, o,&
                              o, o, e/), (/ 3, 3/)))
    !
    i2= transpose(reshape( (/ e, o, o,&
                              o, o, e,&
                              o, e, o/), (/ 3, 3/)))
    !
    sxy = transpose(reshape( (/ e,  o, &
                                o, -e /), (/ 2, 2/)))
    !
    select case(ioper)
      !
    case (1) ! E
      !
      dst(1:18) = src(1:18)
      !
    case (2) !C(+)/(123)(456)
      !
      dst(1) = src(1)
      !
      dst(2:4) = matmul(a123,src(2:4))
      dst(5:7) = matmul(a123,src(5:7))
      !
      dst(8:10) = matmul(a123,src(8:10))
      dst(11:13) = matmul(a123,src(11:13))
      !!
      dst(14:15) = matmul(c123,src(14:15))
      dst(16:17) = matmul(c123,src(16:17))
      !
      dst(18) = src(18)  - 4.0_ark/3.0_ark*pi
      !
      do while(dst(18) < 0.0_ark) 
            dst(18) = dst(18) + 4.0_ark*pi
      enddo
      !
      do while(dst(18) > 4.0_ark*pi) 
            dst(18) = dst(18) - 4.0_ark*pi
      enddo
      !    
    case (4) !sxy(+)/(14)(26)(35)(ab)* 
      !
      dst(1) = src(1)
      !
      dst(2:4) = matmul(i2,src(5:7))
      dst(5:7) = matmul(i2,src(2:4))
      !
      dst(8:10) = matmul(i2,src(11:13))
      dst(11:13) = matmul(i2,src(8:10))
      !
      !!
      dst(14:15) = matmul(sxy,src(16:17))
      dst(16:17) = matmul(sxy,src(14:15))
      !
      dst(18) =  4.0_ark*pi - src(18)
      !
      do while(dst(18) < 0.0_ark) 
            dst(18) = dst(18) + 4.0_ark*pi
      enddo
      do while(dst(18) > 4.0_ark*pi) 
            dst(18) = dst(18) - 4.0_ark*pi
      enddo
      !
    case (7) ! C(-)/(132)(456)
      !
      dst(1) = src(1)
      !
      dst(2:4) = matmul(a132,src(2:4))
      dst(5:7) = matmul(a123,src(5:7))
      !
      dst(8:10) = matmul(a132,src(8:10))
      dst(11:13) = matmul(a123,src(11:13))
      !!
      dst(14:15) = matmul(c132,src(14:15))
      dst(16:17) = matmul(c123,src(16:17))
      !
      dst(18) = src(18)
      !
    case (19) !sxy(-)/(14)(25)(36)(ab)
      !
      dst(1) = src(1)
      !
      dst(2:4) = matmul(i,src(5:7))
      dst(5:7) = matmul(i,src(2:4))
      !
      dst(8:10) = matmul(i,src(11:13))
      dst(11:13) = matmul(i,src(8:10))
      !!
      dst(14:15) = src(16:17)
      dst(16:17) = src(14:15)
      !!
      dst(18) = src(18)
      !
   case(37) !E'
      !
      dst(1:17) = src(1:17)
      dst(18) = src(18) + 2.0_ark*pi
      do while(dst(18) < 0.0_ark) 
           dst(18) = dst(18) + 4.0_ark*pi
      enddo
      do while(dst(18) > 4.0_ark*pi) 
           dst(18) = dst(18) - 4.0_ark*pi
      enddo
      !
    end select
    !
    if (all(tn(ioper,:)/=0)) then
        call ML_symmetry_transformation_C2H6_G36_full(tn(ioper,1),nmodes,src,tmp)
        call ML_symmetry_transformation_C2H6_G36_full(tn(ioper,2),nmodes,tmp,dst)
    endif 
    !
  end subroutine ML_symmetry_transformation_C2H6_G36_full


function MLpoten_c2h6_88_cos3tau_G36(ncoords, natoms, local, xyz, force) result(f)
  !
  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f
  !
  real(ark) :: xi(18),r1,r2,r3,r4,r5,r6,r7,r1e,r2e,betae,a,b
  real(ark) :: beta1,beta2,beta3,beta4,beta5,beta6
  real(ark) :: chi(18,36),term,rad,rhobar
  integer(ik) :: ioper,ipower(18),i

  real(ark) :: tau14,tau24,tau25,tau35,tau36,theta12,theta23,theta13,theta56,theta45,theta46,xi_A,xi_B,xi_C,xi_D
  real(ark) :: tau41,tau51,tau52,tau62,tau63,theta31,theta64,tau16,tau53,tau34,tau26
  !
  rad = pi/180.0_ark
  !
  !r1      = local(1)
  !r2      = local(2)
  !r3      = local(4)
  !r4      = local(6)
  !r5      = local(3)
  !r6      = local(5)
  !r7      = local(7)
  !
  r1e = force(1)
  r2e = force(2)
  betae = force(3)*rad
  a = force(4)
  b = force(5)
  !
  select case(trim(molec%coords_transform))
    !
  case default
    !
    write(out, '(/a,1x,a,1x,a)') &
    'MLpoten_c2h6_88 error', trim(molec%coords_transform), 'is unknown'
    stop 'MLpoten_c2h6_88 error error: bad coordinate type'
    !
  case('R-R16-BETA16-THETA-TAU-1')
    !
    r1 = local(1)
    r2 = local(2)
    r3 = local(6)
    r4 = local(4)
    r5 = local(3)
    r6 = local(5)
    r7 = local(7)
    !
    xi(1)=1.0_ark-exp(-a*(r1-r1e))
    xi(2)=1.0_ark-exp(-b*(r2-r2e))
    xi(3)=1.0_ark-exp(-b*(r3-r2e))
    xi(4)=1.0_ark-exp(-b*(r4-r2e))
    xi(5)=1.0_ark-exp(-b*(r5-r2e))
    xi(6)=1.0_ark-exp(-b*(r6-r2e))
    xi(7)=1.0_ark-exp(-b*(r7-r2e))
    !
    xi(8)  = local(8)  - betae
    xi(9)  = local(12) - betae
    xi(10) = local(10) - betae
    xi(11) = local(9)  - betae
    xi(12) = local(11) - betae
    xi(13) = local(13) - betae


    tau14 = mod(local(14)+4.0_ark*pi,4.0_ark*pi)
    tau34 = mod(local(15)+2.0_ark*pi,2.0_ark*pi)
    tau36 = mod(local(16)+2.0_ark*pi,2.0_ark*pi)
    tau26 = mod(local(17)+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod(local(18)+2.0_ark*pi,2.0_ark*pi)
    !
    ! assuming this is the 404-type (0..720) for tau14, tau25 and tau36 are extended to 0-720 as well
    if (tau14>2.0_ark*pi) then 
       tau25 = tau25 + 2.0_ark*pi
       tau36 = tau36 + 2.0_ark*pi
    endif
    !
    ! for oher dihedral modes the extension is not needed and removed by mod(2 pi)
    !
    tau14 = mod(tau14+2.0_ark*pi,2.0_ark*pi)
    tau34 = mod(tau34+2.0_ark*pi,2.0_ark*pi)
    tau36 = mod(tau36+2.0_ark*pi,2.0_ark*pi)
    tau26 = mod(tau26+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod(tau25+2.0_ark*pi,2.0_ark*pi)
    !
    theta13 = mod(tau14-tau34+2.0_ark*pi,2.0_ark*pi)
    theta23 = mod(tau36-tau26+2.0_ark*pi,2.0_ark*pi)
    theta12 = mod(2.0_ark*pi-theta13-theta23+2.0_ark*pi,2.0_ark*pi)
    !
    theta46 = mod(tau36-tau34+2.0_ark*pi,2.0_ark*pi)
    theta56 = mod(tau25-tau26+2.0_ark*pi,2.0_ark*pi)
    theta45 = mod(2.0_ark*pi-theta46-theta56+2.0_ark*pi,2.0_ark*pi)
    !
    xi(14)  = ( 2.0_ark*theta23 - theta13 - theta12 )/sqrt(6.0_ark)
    xi(15)  = (                   theta13 - theta12 )/sqrt(2.0_ark)
    !
    xi(16)  = ( 2.0_ark*theta56 - theta45 - theta46 )/sqrt(6.0_ark)
    xi(17)  = (                   theta45 - theta46 )/sqrt(2.0_ark)
    !
    rhobar = ( tau14+tau25+tau36 )/3.0_ark
    !
    xi(18) = 1.0_ark + cos(3.0_ark*rhobar)
    !
  case('R-R16-BETA16-THETA-TAU-2')
    !
    r1 = local(1)
    r2 = local(2)
    r3 = local(6)
    r4 = local(4)
    r5 = local(3)
    r6 = local(5)
    r7 = local(7)
    !
    xi(1)=1.0_ark-exp(-a*(r1-r1e))
    xi(2)=1.0_ark-exp(-b*(r2-r2e))
    xi(3)=1.0_ark-exp(-b*(r3-r2e))
    xi(4)=1.0_ark-exp(-b*(r4-r2e))
    xi(5)=1.0_ark-exp(-b*(r5-r2e))
    xi(6)=1.0_ark-exp(-b*(r6-r2e))
    xi(7)=1.0_ark-exp(-b*(r7-r2e))
    !
    xi(8)  = local(8)  - betae
    xi(9)  = local(12) - betae
    xi(10) = local(10) - betae
    xi(11) = local(9)  - betae
    xi(12) = local(11) - betae
    xi(13) = local(13) - betae

    tau14 = mod(local(14)+4.0_ark*pi,4.0_ark*pi)
    tau34 = mod(local(15)+2.0_ark*pi,2.0_ark*pi)
    tau36 = mod(local(16)+2.0_ark*pi,2.0_ark*pi)
    tau26 = mod(local(17)+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod(local(18)+2.0_ark*pi,2.0_ark*pi)
    !
    ! assuming this is the 404-type (0..720) for tau14, tau25 and tau36 are extended to 0-720 as well
    if (tau14>2.0_ark*pi) then 
       tau25 = tau25 + 2.0_ark*pi
       tau36 = tau36 + 2.0_ark*pi
    endif
    !
    ! for oher dihedral modes the extension is not needed and removed by mod(2 pi)
    !
    tau14 = mod(tau14+2.0_ark*pi,2.0_ark*pi)
    tau34 = mod(tau34+2.0_ark*pi,2.0_ark*pi)
    tau36 = mod(tau36+2.0_ark*pi,2.0_ark*pi)
    tau26 = mod(tau26+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod(tau25+2.0_ark*pi,2.0_ark*pi)
    !
    theta13 = mod(tau14-tau34+2.0_ark*pi,2.0_ark*pi)
    theta23 = mod(tau36-tau26+2.0_ark*pi,2.0_ark*pi)
    theta12 = mod(2.0_ark*pi-theta13-theta23+2.0_ark*pi,2.0_ark*pi)
    !
    theta46 = mod(tau36-tau34+2.0_ark*pi,2.0_ark*pi)
    theta56 = mod(tau25-tau26+2.0_ark*pi,2.0_ark*pi)
    theta45 = mod(2.0_ark*pi-theta46-theta56+2.0_ark*pi,2.0_ark*pi)
    !
    xi(14)  = ( 2.0_ark*theta23 - theta12 - theta13 )/sqrt(6.0_ark)
    xi(15)  = (                   theta12 - theta13 )/sqrt(2.0_ark)
    !
    xi(16)  = ( 2.0_ark*theta56 - theta46 - theta45 )/sqrt(6.0_ark)
    xi(17)  = (                   theta46 - theta45 )/sqrt(2.0_ark)
    !
    rhobar = ( tau14+tau25+tau36 )/3.0_ark
    !
    xi(18) = 1.0_ark + cos(3.0_ark*rhobar)
    !
  case('R-R16-BETA16-THETA-TAU-3','R-R16-BETA16-THETA-TAU-11','R-R16-BETA16-THETA-TAU-17','R-R16-BETA16-THETA-TAU-18')
    !
    r1 = local(1)
    r2 = local(2)
    r3 = local(4)
    r4 = local(6)
    r5 = local(3)
    r6 = local(5)
    r7 = local(7)
    !
    xi(1)=1.0_ark-exp(-a*(r1-r1e))
    xi(2)=1.0_ark-exp(-b*(r2-r2e))
    xi(3)=1.0_ark-exp(-b*(r3-r2e))
    xi(4)=1.0_ark-exp(-b*(r4-r2e))
    xi(5)=1.0_ark-exp(-b*(r5-r2e))
    xi(6)=1.0_ark-exp(-b*(r6-r2e))
    xi(7)=1.0_ark-exp(-b*(r7-r2e))
    !
    xi(8)  = local(8)  - betae
    xi(9)  = local(10) - betae
    xi(10) = local(12) - betae
    xi(11) = local(9)  - betae
    xi(12) = local(11) - betae
    xi(13) = local(13) - betae
    !
    tau14 = mod(local(14)+4.0_ark*pi,4.0_ark*pi)
    tau24 = mod(local(15)+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod(local(16)+2.0_ark*pi,2.0_ark*pi)
    tau35 = mod(local(17)+2.0_ark*pi,2.0_ark*pi)
    tau36 = mod(local(18)+2.0_ark*pi,2.0_ark*pi)
    !
    if (tau14>2.0_ark*pi) then 
       tau25 = tau25 + 2.0_ark*pi
       tau36 = tau36 + 2.0_ark*pi
    endif
    !
    !rhobar  = ( tau14+tau25+tau36 )/(3.0_ark)
    !
    tau14 = mod(local(14)+2.0_ark*pi,2.0_ark*pi)
    tau24 = mod(local(15)+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod(local(16)+2.0_ark*pi,2.0_ark*pi)
    tau35 = mod(local(17)+2.0_ark*pi,2.0_ark*pi)
    tau36 = mod(local(18)+2.0_ark*pi,2.0_ark*pi)
    !
    theta12 = mod(tau14-tau24+2.0_ark*pi,2.0_ark*pi)
    theta23 = mod(tau25-tau35+2.0_ark*pi,2.0_ark*pi)
    theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
    !
    theta56 = mod(tau36-tau35+2.0_ark*pi,2.0_ark*pi)
    theta45 = mod(tau25-tau24+2.0_ark*pi,2.0_ark*pi)
    theta46 = mod(2.0_ark*pi-theta56-theta45+2.0_ark*pi,2.0_ark*pi)
    !
    xi(14)  = ( 2.0_ark*theta23 - theta13 - theta12 )/sqrt(6.0_ark)
    xi(15)  = (                   theta13 - theta12 )/sqrt(2.0_ark)
    xi(16)  = ( 2.0_ark*theta56 - theta45 - theta46 )/sqrt(6.0_ark)
    xi(17)  = (                   theta46 - theta45 )/sqrt(2.0_ark)
    !
    rhobar = ( tau14+tau25+tau36 )/3.0_ark
    !
    xi(18) = 1.0_ark + cos(3.0_ark*rhobar)    !
  case('R-R16-BETA16-THETA-TAU-12','R-R16-BETA16-THETA-TAU-13','R-R16-BETA16-THETA-TAU-14',&
       'R-R16-BETA16-THETA-TAU-16','R-R16-BETA16-THETA-TAU-19')
    !
    r1 = local(1)
    r2 = local(3)
    r3 = local(5)
    r4 = local(7)
    r5 = local(2)
    r6 = local(6)
    r7 = local(4)
    !
    xi(1)=1.0_ark-exp(-a*(r1-r1e))
    xi(2)=1.0_ark-exp(-b*(r2-r2e))
    xi(3)=1.0_ark-exp(-b*(r3-r2e))
    xi(4)=1.0_ark-exp(-b*(r4-r2e))
    xi(5)=1.0_ark-exp(-b*(r5-r2e))
    xi(6)=1.0_ark-exp(-b*(r6-r2e))
    xi(7)=1.0_ark-exp(-b*(r7-r2e))
    !
    xi(8)  = local(9)  - betae
    xi(9)  = local(11) - betae
    xi(10) = local(13) - betae
    xi(11) = local(8)  - betae
    xi(12) = local(12) - betae
    xi(13) = local(10) - betae
    !
    tau41 = mod(local(14)+4.0_ark*pi,4.0_ark*pi)
    tau16 = mod(local(15)+2.0_ark*pi,2.0_ark*pi)
    tau62 = mod(local(16)+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod(local(17)+2.0_ark*pi,2.0_ark*pi)
    tau53 = mod(local(18)+2.0_ark*pi,2.0_ark*pi)
    !
    ! assuming this is the 404-type (0..720) for tau14, tau25 and tau36 are extended to 0-720 as well
    if (tau41>2.0_ark*pi) then 
       tau53 = tau53 + 2.0_ark*pi
       tau62 = tau62 + 2.0_ark*pi
    endif
    !
    rhobar  = ( tau41+tau53+tau62)/(3.0_ark)
    !
    !
    tau41 = mod(tau41+2.0_ark*pi,2.0_ark*pi)
    tau16 = mod(tau16+2.0_ark*pi,2.0_ark*pi)
    tau62 = mod(tau62+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod(tau25+2.0_ark*pi,2.0_ark*pi)
    tau53 = mod(tau53+2.0_ark*pi,2.0_ark*pi)
    !
    theta12 = mod(tau62-tau16+2.0_ark*pi,2.0_ark*pi)
    theta23 = mod(tau53-tau25+2.0_ark*pi,2.0_ark*pi)
    theta31 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
    !
    theta56 = mod(tau25-tau62+2.0_ark*pi,2.0_ark*pi)
    theta64 = mod(tau16-tau41+2.0_ark*pi,2.0_ark*pi)
    theta45 = mod(2.0_ark*pi-theta56-theta64+2.0_ark*pi,2.0_ark*pi)
    !
    xi(14)  = ( 2.0_ark*theta23 - theta31 - theta12 )/sqrt(6.0_ark)
    xi(15)  = (                   theta31 - theta12 )/sqrt(2.0_ark)
    !
    xi(16)  = ( 2.0_ark*theta56 - theta64 - theta45 )/sqrt(6.0_ark)
    xi(17)  = (                   theta64 - theta45 )/sqrt(2.0_ark)
    !
    xi(18) = 1.0_ark + cos(3.0_ark*rhobar)
    !  
  case('R-R16-BETA16-THETA-TAU')
    !
    r1 = local(1)
    r2 = local(2)
    r3 = local(4)
    r4 = local(6)
    r5 = local(3)
    r6 = local(5)
    r7 = local(7)
    !
    xi(1)=1.0_ark-exp(-a*(r1-r1e))
    xi(2)=1.0_ark-exp(-b*(r2-r2e))
    xi(3)=1.0_ark-exp(-b*(r3-r2e))
    xi(4)=1.0_ark-exp(-b*(r4-r2e))
    xi(5)=1.0_ark-exp(-b*(r5-r2e))
    xi(6)=1.0_ark-exp(-b*(r6-r2e))
    xi(7)=1.0_ark-exp(-b*(r7-r2e))
    !
    tau14 = local(14)
    tau24 = local(15)
    tau25 = local(16)
    tau35 = local(17)
    tau36 = local(18)
    !
    tau14 = mod(tau14+4.0_ark*pi,4.0_ark*pi)
    tau24 = mod(tau24+4.0_ark*pi,4.0_ark*pi)
    tau25 = mod(tau25+4.0_ark*pi,4.0_ark*pi)
    tau35 = mod(tau35+4.0_ark*pi,4.0_ark*pi)
    tau36 = mod(tau36+4.0_ark*pi,4.0_ark*pi)
    !
    rhobar  = ( tau14+tau25+tau36 )/(3.0_ark)-pi
    !
    tau14 = mod(tau14+2.0_ark*pi,2.0_ark*pi)
    tau24 = mod(tau24+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod(tau25+2.0_ark*pi,2.0_ark*pi)
    tau35 = mod(tau35+2.0_ark*pi,2.0_ark*pi)
    tau36 = mod(tau36+2.0_ark*pi,2.0_ark*pi)
    !
    theta12 = mod(tau14-tau24+2.0_ark*pi,2.0_ark*pi)
    theta23 = mod(tau25-tau35+2.0_ark*pi,2.0_ark*pi)
    theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
    !
    theta56 = mod(tau36-tau35+2.0_ark*pi,2.0_ark*pi)
    theta45 = mod(tau25-tau24+2.0_ark*pi,2.0_ark*pi)
    theta46 = mod(2.0_ark*pi-theta56-theta45+2.0_ark*pi,2.0_ark*pi)
    !
    xi(14)  = ( 2.0_ark*theta23 - theta13 - theta12 )/sqrt(12.0_ark)
    xi(15)  = (                   theta13 - theta12 )/(2.0_ark)
    xi(16)  = ( 2.0_ark*theta56 - theta46 - theta45 )/sqrt(12.0_ark)
    xi(17)  = (                   theta46 - theta45 )/(2.0_ark)
    !
    xi(18) = 1.0_ark + cos(3.0_ark*rhobar)

    xi(8)  = local(8)  - betae
    xi(9)  = local(10) - betae
    xi(10) = local(12) - betae
    xi(11) = local(9)  - betae
    xi(12) = local(11) - betae
    xi(13) = local(13) - betae
    !
  end select
  !
  f = 0
  !
  do ioper = 1,36
    !
    ! for xi(18) = 1.0_ark + cos(3.0_ark*rhobar) 
    ! all operations on xi18->xi18
    !
    call ML_symmetry_transformation_C2H6_G36(ioper,18,xi,chi(:,ioper))
    !
  enddo
  ! 
  do i = 6, molec%parmax
    !
    ipower(1:18) = molec%pot_ind(1:18,i)
    !
    term = 0 
    !
    do ioper = 1,36
      !
      term = term + product(chi(1:18,ioper)**ipower(1:18))
      !
    end do
    !
    term = term/36.0_ark
    !
    f = f + term*force(i)
    !
  enddo
  !  
end function MLpoten_c2h6_88_cos3tau_G36


function MLpoten_c2h6_Morse_displace_fourier_G36(ncoords, natoms, local, xyz, force) result(f)
  !
  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f
  !
  real(ark) :: xi(18),x(18),r1,r2,r3,r4,r5,r6,r7,r1e,r2e,betae,a,b
  real(ark) :: beta1,beta2,beta3,beta4,beta5,beta6
  real(ark) :: chi(18,36),term,rad,rhobar
  integer(ik) :: ioper,ipower(18),i,nmodes,Nsin,Ncos
  !
  real(ark) :: tau14,tau24,tau25,tau35,tau36,theta12,theta23,theta13,theta56,theta45,theta46
  real(ark) :: tau41,tau51,tau52,tau62,tau63,theta31,theta64,tau16,tau53
  !
  rad = pi/180.0_ark
  !
  r1e = force(1)
  r2e = force(2)
  betae = force(3)*rad
  a = force(4)
  b = force(5)
  Ncos = force(6)
  Nsin = force(6+Ncos+1)
  !
  nmodes = 18
  !
  call coordinate_transformation(ncoords,nmodes,local,xi)
  !
  do ioper = 1,36
    !
    call ML_symmetry_transformation_C2H6_G36_full(ioper,18,xi,chi(:,ioper))
    !
  enddo
  !
  f = 0
  ! 
  do i = 7, 7+Ncos
    !
    ipower(1:18) = molec%pot_ind(1:18,i)
    !
    term = 0 
    do ioper = 1,36
      x( 1   )= 1.0_ark   - exp(-a*(chi(1,ioper)-r1e))
      x( 2: 7)= 1.0_ark   - exp(-b*(chi(2:7,ioper)-r2e))
      x( 8:13)= chi( 8:13,ioper) - betae
      x(14:18)= chi(14:18,ioper)
      term = term + product(x(1:17)**ipower(1:17))*cos(ipower(18)*x(18))
      !
    end do
    !
    term = term/36.0_ark
    !
    f = f + term*force(i)
    !
  enddo
  ! 
  do i =7+Ncos+1,7+Ncos+Nsin
    !
    ipower(1:18) = molec%pot_ind(1:18,i)
    !
    term = 0 
    do ioper = 1,36
      x( 1   )= 1.0_ark   - exp(-a*(chi(1,ioper)-r1e))
      x( 2: 7)= 1.0_ark   - exp(-b*(chi(2:7,ioper)-r2e))
      x( 8:13)= chi( 8:13,ioper) - betae
      x(14:18)= chi(14:18,ioper)
      term = term + product(x(1:17)**ipower(1:17))*sin(ipower(18)*x(18))
      !
    end do
    !
    term = term/36.0_ark
    !
    f = f + term*force(i)
    !
  enddo
  !  
end function MLpoten_c2h6_Morse_displace_fourier_G36



function MLpoten_c2h6_explct_M2_P6(ncoords, natoms, local, xyz, force) result(f)
  !
  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f
  !
  real(ark) :: xi(18),r1e,r2e,betae,a,b,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,s1
  real(ark) :: rad,g(molec%parmax)
  integer(ik) :: i,nmodes
  !
  rad = pi/180.0_ark
  !
  !r1      = local(1)
  !r2      = local(2)
  !r3      = local(4)
  !r4      = local(6)
  !r5      = local(3)
  !r6      = local(5)
  !r7      = local(7)
  !
  r1e = force(1)
  r2e = force(2)
  betae = force(3)*rad
  a = force(4)
  b = force(5)
  !
  nmodes = 18
  !
  call coordinate_transformation(ncoords,nmodes,local,xi)
  !
  y1  = 1.0_ark-exp(-a*(xi( 1)-r1e))
  y2  = 1.0_ark-exp(-b*(xi( 2)-r2e))
  y3  = 1.0_ark-exp(-b*(xi( 3)-r2e))
  y4  = 1.0_ark-exp(-b*(xi( 4)-r2e))
  y5  = 1.0_ark-exp(-b*(xi( 5)-r2e))
  y6  = 1.0_ark-exp(-b*(xi( 6)-r2e))
  y7  = 1.0_ark-exp(-b*(xi( 7)-r2e))
  
  
  y8  = xi( 8) - betae
  y9  = xi( 9) - betae
  y10 = xi(10) - betae
  y11 = xi(11) - betae
  y12 = xi(12) - betae
  y13 = xi(13) - betae
  y14 = xi(14)
  y15 = xi(15)
  y16 = xi(16)
  y17 = xi(17)
  y18 = xi(18)

  g = 0
  !
  g(9) = 1._ark
  g(10) = cos(3._ark*y18)
  g(11) = cos(6.0_ark*y18)
  g(12) = cos(9.0_ark*y18)
  g(13) = cos(12.0_ark*y18)


  g(14) = cos(9.0_ark*y18)*(y12+y13+y11+y10+y9+y8)
  g(15) = y12+y13+y11+y10+y9+y8
  g(16) = cos(3.0_ark*y18)*(y12+y13+y11+y10+y9+y8)
  g(17) = cos(6.0_ark*y18)*(y12+y13+y11+y10+y9+y8)
  g(18) = cos(12.0_ark*y18)*(y12+y13+y11+y10+y9+y8)
  g(19) = cos(3.0_ark*y18)*(y7+y5+y6+y4+y3+y2)
  g(20) = cos(9.0_ark*y18)*(y7+y5+y6+y4+y3+y2)
  g(21) = y7+y5+y6+y4+y3+y2
  g(22) = cos(6.0_ark*y18)*(y7+y5+y6+y4+y3+y2)
  g(23) = cos(12.0_ark*y18)*(y7+y5+y6+y4+y3+y2)
  g(24) = y1
  g(25) = y1*cos(3.0_ark*y18)
  g(26) = y1*cos(6.0_ark*y18)
  g(27) = y1*cos(9.0_ark*y18)
  g(28) = y1*cos(12.0_ark*y18)
  g(29) = cos(3.0_ark*y18)*(y17**2+y15**2+y16**2+y14**2)
  g(30) = cos(6.0_ark*y18)*(y17**2+y15**2+y16**2+y14**2)
  g(31) = y17**2+y15**2+y16**2+y14**2
  g(32) = y15*y17*cos(y18)+y15*y16*sin(y18)-y14*y16*cos(y18)+y14*y17*sin(y18)
  g(33) = y15*y17*cos(2.0_ark*y18)-y15*y16*sin(2.0_ark*y18)-y14*y16*cos(2.0_ark*y18)-y14*y17*sin(2.0_ark*y18)
  g(34) = y15*y16*sin(4.0_ark*y18)+y14*y17*sin(4.0_ark*y18)+y15*y17*cos(4.0_ark*y18)-y14*y16*cos(4.0_ark*y18)
  g(35) = y15*y16*sin(5.0_ark*y18)+y14*y16*cos(5.0_ark*y18)+y14*y17*sin(5.0_ark*y18)-y15*y17*cos(5.0_ark*y18)
  g(36) = cos(6.0_ark*y18)*(y11**2+y12**2+y10**2+y9**2+y8**2+y13**2)
  g(37) = y11**2+y12**2+y10**2+y9**2+y8**2+y13**2
  g(38) = cos(3.0_ark*y18)*(y11**2+y12**2+y10**2+y9**2+y8**2+y13**2)
  g(39) = y7**2+y6**2+y2**2+y5**2+y4**2+y3**2
  g(40) = cos(3.0_ark*y18)*(y7**2+y6**2+y2**2+y5**2+y4**2+y3**2)
  g(41) = cos(6.0_ark*y18)*(y7**2+y6**2+y2**2+y5**2+y4**2+y3**2)
  g(42) = y1**2
  g(43) = y1**2*cos(3.0_ark*y18)
  g(44) = y1**2*cos(6.0_ark*y18)
  g(45) = -y15*y17**2*sin(y18)-y14**2*y16*cos(y18)-y14*y16**2*cos(y18)-y15**2*y17*sin(y18)+y15*y16**2*sin(y18)-&
          2.0_ark*y15*y16*y17*cos(y18)-2.0_ark*y14*y15*y17*cos(y18)-2.0_ark*y14*y15*y16*sin(y18)+y14*y17**2*cos(y18)-&
          2.0_ark*y14*y16*y17*sin(y18)+y15**2*y16*cos(y18)+y14**2*y17*sin(y18)
  g(46) = -y16**3/3.0_ark-y14**3/3.0_ark+y16*y17**2+y14*y15**2
  g(47) = y9**3+y8**3+y12**3+y13**3+y11**3+y10**3
  g(48) = y4**3+y6**3+y7**3+y2**3+y5**3+y3**3
  g(49) = y1**3
  g(50) = y14*y15**2*y17*sin(y18)+y15**3*y16*sin(y18)+y15**3*y17*cos(y18)+y14**3*y17*sin(y18)+y14**2*y15*y16*sin(y18)-&
  y14*y16*y17**2*cos(y18)-y14*y16**3*cos(y18)-y14*y15**2*y16*cos(y18)-y14**3*y16*cos(y18)+y15*y16**2*y17*cos(y18)+&
  y15*y16**3*sin(y18)+y15*y16*y17**2*sin(y18)+y15*y17**3*cos(y18)+y14*y17**3*sin(y18)+y14*y16**2*y17*sin(y18)+&
     y14**2*y15*y17*cos(y18)
  g(51) = -y15**2*y17**2*cos(y18)+y14**2*y17**2*cos(y18)+2.0_ark*y15**2*y16*y17*sin(y18)+4.0_ark*y14*y15*y16*y17*cos(y18)+&
        2.0_ark*y14*y15*y17**2*sin(y18)-y14**2*y16**2*cos(y18)+y15**2*y16**2*cos(y18)-2.0_ark*y14*y15*y16**2*sin(y18)-&
        2.0_ark*y14**2*y16*y17*sin(y18)
  g(52) = y14**2*y16**2+y15**2*y17**2+y15**2*y16**2+y14**2*y17**2
  g(53) = 2.0_ark*y16**2*y17**2+y15**4+y17**4+2.0_ark*y14**2*y15**2+y16**4+y14**4
  g(54) = y8**4+y9**4+y11**4+y10**4+y13**4+y12**4
  g(55) = y6**4+y7**4+y2**4+y3**4+y5**4+y4**4
  g(56) = y1**4
  g(57) = -y14**2*y16**2*y17*sin(y18)-y14**2*y17**3*sin(y18)-y14**3*y17**2*cos(y18)+2.0_ark*y14**2*y15*y16*y17*cos(y18)+&
        y14**3*y16**2*cos(y18)+2.0_ark*y14**3*y16*y17*sin(y18)+y14**2*y15*y17**2*sin(y18)+y14**2*y16**3*cos(y18)+&
        y15**2*y17**3*sin(y18)+y14*y15**2*y16**2*cos(y18)-y15**2*y16**3*cos(y18)-y15**2*y16*y17**2*cos(y18)-&
        y15**3*y16**2*sin(y18)+2.0_ark*y15**3*y16*y17*cos(y18)+2.0_ark*y14*y15*y17**3*cos(y18)+2.0_ark*y14*y15*y16**2*y17*cos(y18)+&
        2.0_ark*y14*y15*y16*y17**2*sin(y18)+2.0_ark*y14*y15*y16**3*sin(y18)+y14**2*y16*y17**2*cos(y18)-y14*y15**2*y17**2*cos(y18)+&
        2.0_ark*y14*y15**2*y16*y17*sin(y18)-y14**2*y15*y16**2*sin(y18)+y15**2*y16**2*y17*sin(y18)+y15**3*y17**2*sin(y18)
  g(58) = -3.0_ark*y14**3*y15*y16*sin(y18)+y15*y16*y17**3*cos(y18)+y14*y16*y17**3*sin(y18)+y14*y17**4*cos(y18)+y15**4*y16*cos(y18)-&
        y15**4*y17*sin(y18)+y14*y15**3*y16*sin(y18)+3.0_ark*y14**2*y15**2*y17*sin(y18)-3.0_ark*y14**2*y15**2*y16*cos(y18)-&
        3.0_ark*y14**3*y15*y17*cos(y18)-3.0_ark*y15*y16**3*y17*cos(y18)+3.0_ark*y15*y16**2*y17**2*sin(y18)-y15*y17**4*sin(y18)-&
        3.0_ark*y14*y16**2*y17**2*cos(y18)-3.0_ark*y14*y16**3*y17*sin(y18)+y14*y15**3*y17*cos(y18)
  g(59) = -y14**2*y16**3/3.0_ark-y14**3*y17**2/3.0_ark+y14**2*y16*y17**2+y14*y15**2*y17**2+y15**2*y16*y17**2-y14**3*y16**2/3.0_ark-&
          y15**2*y16**3/3.0_ark+y14*y15**2*y16**2
  g(60) = y14**3*y15**2-y14**5/2.0_ark-y16**5/2.0_ark+y16**3*y17**2+3.0_ark/2.0_ark*y14*y15**4+3.0_ark/2.0_ark*y16*y17**4
  g(61) = -y14**4*y17*sin(y18)+8.0_ark*y14**3*y15*y16*sin(y18)+y14**4*y16*cos(y18)+y14*y16**4*cos(y18)-3.0_ark*y14*y17**4*cos(y18)-&
           3.0_ark*y15**4*y16*cos(y18)+3.0_ark*y15**4*y17*sin(y18)-6.0_ark*y14**2*y15**2*y17*sin(y18)+&
           6.0_ark*y14**2*y15**2*y16*cos(y18)+8.0_ark*y14**3*y15*y17*cos(y18)-y15*y16**4*sin(y18)+8.0_ark*y15*y16**3*y17*cos(y18)-&
           6.0_ark*y15*y16**2*y17**2*sin(y18)+3.0_ark*y15*y17**4*sin(y18)+6.0_ark*y14*y16**2*y17**2*cos(y18)+&
           8.0_ark*y14*y16**3*y17*sin(y18)
  g(62) = y8**5+y9**5+y11**5+y13**5+y10**5+y12**5
  g(63) = y7**5+y2**5+y3**5+y4**5+y6**5+y5**5
  g(64) = y1**5
  g(65) = -6.0_ark*y14**2*y15**4+9.0_ark*y14**4*y15**2-6.0_ark*y16**2*y17**4+y15**6+9.0_ark*y16**4*y17**2+y17**6
  g(66) = 9.0_ark*y14**2*y15**4-6.0_ark*y14**4*y15**2+9.0_ark*y16**2*y17**4+y16**6-6.0_ark*y16**4*y17**2+y14**6
  g(67) = -4.0_ark*y14**2*y15**3*y16*sin(y18)+3.0_ark*y14**4*y15*y16*sin(y18)-6.0_ark*y14**3*y15**2*y16*cos(y18)+&
          2.0_ark*y14*y16*y17**4*cos(y18)+y14*y17**5*sin(y18)+y15**5*y17*cos(y18)+3.0_ark*y15*y16**4*y17*cos(y18)+&
          6.0_ark*y15*y16**3*y17**2*sin(y18)+2.0_ark*y14*y15**4*y16*cos(y18)-4.0_ark*y14**2*y15**3*y17*cos(y18)+&
          6.0_ark*y14**3*y15**2*y17*sin(y18)+3.0_ark*y14**4*y15*y17*cos(y18)+y15*y17**5*cos(y18)- &
          2.0_ark*y15*y16*y17**4*sin(y18)-4.0_ark*y15*y16**2*y17**3*cos(y18)+y15**5*y16*sin(y18)-&
          6.0_ark*y14*y16**3*y17**2*cos(y18)-4.0_ark*y14*y16**2*y17**3*sin(y18)+3.0_ark*y14*y16**4*y17*sin(y18)-&
          2.0_ark*y14*y15**4*y17*sin(y18)
  g(68) = 6.0_ark*y14**2*y15**3*y16*sin(y18)-y14**5*y16*cos(y18)-2.0_ark*y14**4*y15*y16*sin(y18)+&
          4.0_ark*y14**3*y15**2*y16*cos(y18)-3.0_ark*y14*y16*y17**4*cos(y18)-2.0_ark*y15*y16**4*y17*cos(y18)-&
          4.0_ark*y15*y16**3*y17**2*sin(y18)-3.0_ark*y14*y15**4*y16*cos(y18)+ &
          6.0_ark*y14**2*y15**3*y17*cos(y18)-y14*y16**5*cos(y18)-4.0_ark*y14**3*y15**2*y17*sin(y18)-&
          2.0_ark*y14**4*y15*y17*cos(y18)+3.0_ark*y15*y16*y17**4*sin(y18)+6.0_ark*y15*y16**2*y17**3*cos(y18)+&
          y15*y16**5*sin(y18)+y14**5*y17*sin(y18)+ &
          4.0_ark*y14*y16**3*y17**2*cos(y18)+6.0_ark*y14*y16**2*y17**3*sin(y18)-2.0_ark*y14*y16**4*y17*sin(y18)+&
          3.0_ark*y14*y15**4*y17*sin(y18)
  g(69) = -y14**3*y15*y16**2*sin(y18)+2.0_ark*y14**3*y15*y16*y17*cos(y18)+y15**4*y16**2*cos(y18)/3.0_ark+&
          y15**2*y16**2*y17**2*cos(y18)- y15**2*y16*y17**3*sin(y18)/3.0_ark-y14*y15**3*y17**2*sin(y18)/3.0_ark-&
          2.0_ark/3.0_ark*y14*y15*y16*y17**3*cos(y18)+ 2.0_ark/3.0_ark*y14*y15*y17**4*sin(y18)-y14**2*y16**2*y17**2*cos(y18)+&
          y14**2*y16*y17**3*sin(y18)/3.0_ark+y14**2*y17**4*cos(y18)/3.0_ark+2.0_ark*y14*y15*y16**3*y17*cos(y18)- &
          y14**2*y16**3*y17*sin(y18)+y14**3*y15*y17**2*sin(y18)-y15**2*y17**4*cos(y18)/3.0_ark+y15**2*y16**3*y17*sin(y18)- &
          y15**4*y17**2*cos(y18)/3.0_ark+2.0_ark/3.0_ark*y15**4*y16*y17*sin(y18)-2.0_ark*y14*y15*y16**2*y17**2*sin(y18)+&
          y14*y15**3*y16**2*sin(y18)/3.0_ark-&
          2.0_ark/3.0_ark*y14*y15**3*y16*y17*cos(y18)+y14**2*y15**2*y17**2*cos(y18)-2.0_ark*y14**2*y15**2*y16*y17*sin(y18)-&
          y14**2*y15**2*y16**2*cos(y18)
  g(70) = 2.0_ark*y15**2*y16**2*y17**2+y14**4*y17**2+y15**4*y16**2+y15**2*y17**4+y14**2*y17**4+y14**2*y16**4+&
          2.0_ark*y14**2*y16**2*y17**2+2.0_ark*y14**2*y15**2*y17**2+2.0_ark*y14**2*y15**2*y16**2+y14**4*y16**2+&
          y15**2*y16**4+y15**4*y17**2
  g(71) = -y14**4*y16**2*cos(y18)-2.0_ark*y14**4*y16*y17*sin(y18)+y15**4*y16**2*cos(y18)/3.0_ark-&
          2.0_ark*y15**2*y16**2*y17**2*cos(y18)+8.0_ark/3.0_ark*y15**2*y16*y17**3*sin(y18)+&
          8.0_ark/3.0_ark*y14*y15**3*y17**2*sin(y18)+16.0_ark/3.0_ark*y14*y15*y16*y17**3*cos(y18)+&
          2.0_ark/3.0_ark*y14*y15*y17**4*sin(y18)+2.0_ark*y14**2*y16**2*y17**2*cos(y18)- &
          8.0_ark/3.0_ark*y14**2*y16*y17**3*sin(y18)+y14**2*y17**4*cos(y18)/3.0_ark-y14**2*y16**4*cos(y18)+&
          y14**4*y17**2*cos(y18)-y15**2*y17**4*cos(y18)/3.0_ark+y15**2*y16**4*cos(y18)-y15**4*y17**2*cos(y18)/3.0_ark+&
          2.0_ark/3.0_ark*y15**4*y16*y17*sin(y18)+ 4.0_ark*y14*y15*y16**2*y17**2*sin(y18)-2.0_ark*y14*y15*y16**4*sin(y18)-&
          8.0_ark/3.0_ark*y14*y15**3*y16**2*sin(y18)+16.0_ark/3.0_ark*y14*y15**3*y16*y17*cos(y18)-&
          2.0_ark*y14**2*y15**2*y17**2*cos(y18)+4.0_ark*y14**2*y15**2*y16*y17*sin(y18)+2.0_ark*y14**2*y15**2*y16**2*cos(y18)
  g(72) = -3.0_ark*y15**3*y16**2*y17-3.0_ark*y14**2*y15*y17**3+9.0_ark*y14**2*y15*y16**2*y17+y15**3*y17**3
  g(73) = -y14**3*y16*y17**2*cos(y18)+y15**3*y16**3*sin(y18)+y15**3*y16*y17**2*sin(y18)+y14*y15**2*y16**2*y17*sin(y18)+ &
           y14**2*y15*y16*y17**2*sin(y18)-y14*y15**2*y16**3*cos(y18)-y14*y15**2*y16*y17**2*cos(y18)+y14**3*y17**3*sin(y18)&
          -y14**3*y16**3*cos(y18)+y14**3*y16**2*y17*sin(y18)+y15**3*y17**3*cos(y18)+y15**3*y16**2*y17*cos(y18)+&
          y14*y15**2*y17**3*sin(y18)+y14**2*y15*y17**3*cos(y18)+y14**2*y15*y16**2*y17*cos(y18)+y14**2*y15*y16**3*sin(y18)
  g(74) = -3.0_ark*y14**3*y16*y17**2-3.0_ark*y14*y15**2*y16**3+9.0_ark*y14*y15**2*y16*y17**2+y14**3*y16**3
  g(75) = y9**6+y10**6+y12**6+y8**6+y13**6+y11**6
  g(76) = y5**6+y7**6+y6**6+y3**6+y4**6+y2**6
  g(77) = y1**6
  !
  f = 0 
  ! 
  !do i = 9,molec%parmax
  !  !
  !  f = f + g(i)*force(i)
  !  !
  !enddo
  !
  f = sum(g(9:)*force(9:))
  !  
end function MLpoten_c2h6_explct_M2_P6


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



!define cartesian components of the polarizabilityin space-fixed system
 !
 !
 recursive subroutine ML_alpha_C2H6_zero_order(rank,ncoords,natoms,local,xyz,f)
    !
    implicit none
    !
    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
    real(ark),intent(out)  ::  f(rank)
    !
    real(ark)              :: e1(3),e2(3),e3(3),rad,r1e,r2e,betae,r1,r2,r3,r4,r5,r6,r7,xi(18)
    real(ark)              :: tau14,tau24,tau25,tau35,tau36,theta12,theta23,theta13,theta56,theta45,theta46
    real(ark)              :: rhobar

    !
    e1(:) = xyz(1,:)-xyz(2,:)
    e2(:) = xyz(3,:)-xyz(1,:)
    e3(:) = xyz(4,:)-xyz(2,:)
    !
    rad = pi/180.0_ark
    !
    r1e   = extF%coef(1, 1)
    r2e   = extF%coef(2, 1)
    betae = extF%coef(3, 1)*rad
    !
    select case(trim(molec%coords_transform))
      !
    case default
      !
      write(out, '(/a,1x,a,1x,a)') &
      'ML_alpha_C2H6_zero_order error', trim(molec%coords_transform), 'is unknown'
      stop 'ML_alpha_C2H6_zero_order error: bad coordinate type'
      !
    case('R-R16-BETA16-THETA-TAU-1','R-R16-BETA16-THETA-TAU-2','R-R16-BETA16-THETA-TAU-3',&
         'R-R16-BETA16-THETA-TAU-11','R-R16-BETA16-THETA-TAU-17','R-R16-BETA16-THETA-TAU-18')
      !
      r1 = local(1)
      r2 = local(2)
      r3 = local(4)
      r4 = local(6)
      r5 = local(3)
      r6 = local(7)
      r7 = local(5)
      !
      xi(1)=(r1-r1e)
      xi(2)=(r2-r2e)
      xi(3)=(r3-r2e)
      xi(4)=(r4-r2e)
      xi(5)=(r5-r2e)
      xi(6)=(r6-r2e)
      xi(7)=(r7-r2e)
      !
      xi(8)  = local(8)  - betae
      xi(9)  = local(10) - betae
      xi(10) = local(12) - betae
      xi(11) = local(9)  - betae
      xi(12) = local(13) - betae
      xi(13) = local(11) - betae
      !
      tau14 = mod(local(14)+4.0_ark*pi,4.0_ark*pi)
      tau24 = mod(local(15)+2.0_ark*pi,2.0_ark*pi)
      tau25 = mod(local(16)+2.0_ark*pi,2.0_ark*pi)
      tau35 = mod(local(17)+2.0_ark*pi,2.0_ark*pi)
      tau36 = mod(local(18)+2.0_ark*pi,2.0_ark*pi)
      !
      if (tau14>2.0_ark*pi) then 
         tau25 = tau25 + 2.0_ark*pi
         tau36 = tau36 + 2.0_ark*pi
      endif
      !
      rhobar  = ( tau14+tau25+tau36 )/(3.0_ark)
      !
      tau14 = mod(local(14)+2.0_ark*pi,2.0_ark*pi)
      tau24 = mod(local(15)+2.0_ark*pi,2.0_ark*pi)
      tau25 = mod(local(16)+2.0_ark*pi,2.0_ark*pi)
      tau35 = mod(local(17)+2.0_ark*pi,2.0_ark*pi)
      tau36 = mod(local(18)+2.0_ark*pi,2.0_ark*pi)
      !
      theta12 = mod(tau14-tau24+2.0_ark*pi,2.0_ark*pi)
      theta23 = mod(tau25-tau35+2.0_ark*pi,2.0_ark*pi)
      theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
      !
      theta56 = mod(tau36-tau35+2.0_ark*pi,2.0_ark*pi)
      theta45 = mod(tau25-tau24+2.0_ark*pi,2.0_ark*pi)
      theta46 = mod(2.0_ark*pi-theta56-theta45+2.0_ark*pi,2.0_ark*pi)
      !
      xi(14)  = ( 2.0_ark*theta23 - theta13 - theta12 )/sqrt(6.0_ark)
      xi(15)  = (                   theta13 - theta12 )/sqrt(2.0_ark)
      xi(16)  = ( 2.0_ark*theta56 - theta45 - theta46 )/sqrt(6.0_ark)
      xi(17)  = (                   theta45 - theta46 )/sqrt(2.0_ark)
      !
      rhobar = ( tau14+tau25+tau36 )/3.0_ark
      !
      xi(18) = cos(3.0_ark*rhobar)
      !
    end select
    !
    f(1) = extF%coef(4,1)+extF%coef(5,1)*xi(18) !(1,1)
    f(2) = extF%coef(1,2)+extF%coef(2,2)*xi(18) !(1,2)
    f(3) = extF%coef(1,3)+extF%coef(2,3)*xi(18) !(1,3)
    f(4) = extF%coef(1,4)+extF%coef(2,4)*xi(18) !(2,2)
    f(5) = extF%coef(1,5)+extF%coef(2,5)*xi(18) !(2,3)
    f(6) = extF%coef(1,6)+extF%coef(2,6)*xi(18) !(3,3)
    !
  end subroutine ML_alpha_C2H6_zero_order


function MLpoten_c2h6_88_cos3tau_sin3tau_G36(ncoords, natoms, local, xyz, force) result(f)
  !
  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f
  !
  real(ark) :: xi(18),r1,r2,r3,r4,r5,r6,r7,r1e,r2e,betae,a,b
  real(ark) :: beta1,beta2,beta3,beta4,beta5,beta6
  real(ark) :: chi(18,36),term,rad,rhobar,d18,s18
  integer(ik) :: ioper,ipower(18),i

  real(ark) :: tau14,tau24,tau25,tau35,tau36,theta12,theta23,theta13,theta56,theta45,theta46,xi_A,xi_B,xi_C,xi_D
  real(ark) :: tau41,tau51,tau52,tau62,tau63,theta31,theta64,tau16,tau53
  !
  rad = pi/180.0_ark
  !
  !r1      = local(1)
  !r2      = local(2)
  !r3      = local(4)
  !r4      = local(6)
  !r5      = local(3)
  !r6      = local(5)
  !r7      = local(7)
  !
  r1e = force(1)
  r2e = force(2)
  betae = force(3)*rad
  a = force(4)
  b = force(5)
  !
  select case(trim(molec%coords_transform))
    !
  case default
    !
    write(out, '(/a,1x,a,1x,a)') &
    'ML_alpha_C2H6_zero_order error', trim(molec%coords_transform), 'is unknown'
    stop 'MLpoten_c2h6_88 error error: bad coordinate type'
    !
  case('R-R16-BETA16-THETA-TAU-1','R-R16-BETA16-THETA-TAU-2','R-R16-BETA16-THETA-TAU-11',&
       'R-R16-BETA16-THETA-TAU-17','R-R16-BETA16-THETA-TAU-18','R-R16-BETA16-THETA-TAU-3')
    !
    r1 = local(1)
    r2 = local(2)
    r3 = local(4)
    r4 = local(6)
    r5 = local(3)
    r6 = local(5)
    r7 = local(7)
    !
    xi(1)=1.0_ark-exp(-a*(r1-r1e))
    xi(2)=1.0_ark-exp(-b*(r2-r2e))
    xi(3)=1.0_ark-exp(-b*(r3-r2e))
    xi(4)=1.0_ark-exp(-b*(r4-r2e))
    xi(5)=1.0_ark-exp(-b*(r5-r2e))
    xi(6)=1.0_ark-exp(-b*(r6-r2e))
    xi(7)=1.0_ark-exp(-b*(r7-r2e))
    !
    xi(8)  = local(8)  - betae
    xi(9)  = local(10) - betae
    xi(10) = local(12) - betae
    xi(11) = local(9)  - betae
    xi(12) = local(11) - betae
    xi(13) = local(13) - betae
    !
    tau14 = mod(local(14)+4.0_ark*pi,4.0_ark*pi)
    tau24 = mod(local(15)+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod(local(16)+2.0_ark*pi,2.0_ark*pi)
    tau35 = mod(local(17)+2.0_ark*pi,2.0_ark*pi)
    tau36 = mod(local(18)+2.0_ark*pi,2.0_ark*pi)
    !
    if (tau14>2.0_ark*pi) then 
       tau25 = tau25 + 2.0_ark*pi
       tau36 = tau36 + 2.0_ark*pi
    endif
    !
    !rhobar  = ( tau14+tau25+tau36 )/(3.0_ark)
    !
    tau14 = mod(local(14)+2.0_ark*pi,2.0_ark*pi)
    tau24 = mod(local(15)+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod(local(16)+2.0_ark*pi,2.0_ark*pi)
    tau35 = mod(local(17)+2.0_ark*pi,2.0_ark*pi)
    tau36 = mod(local(18)+2.0_ark*pi,2.0_ark*pi)
    !
    theta12 = mod(tau14-tau24+2.0_ark*pi,2.0_ark*pi)
    theta23 = mod(tau25-tau35+2.0_ark*pi,2.0_ark*pi)
    theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
    !
    theta56 = mod(tau36-tau35+2.0_ark*pi,2.0_ark*pi)
    theta45 = mod(tau25-tau24+2.0_ark*pi,2.0_ark*pi)
    theta46 = mod(2.0_ark*pi-theta56-theta45+2.0_ark*pi,2.0_ark*pi)
    !
    xi(14)  = ( 2.0_ark*theta23 - theta13 - theta12 )/sqrt(6.0_ark)
    xi(15)  = (                   theta13 - theta12 )/sqrt(2.0_ark)
    xi(16)  = ( 2.0_ark*theta56 - theta45 - theta46 )/sqrt(6.0_ark)
    xi(17)  = (                   theta46 - theta45 )/sqrt(2.0_ark)
    !
    rhobar = ( tau14+tau25+tau36 )/3.0_ark
    !
    xi(18) = 1.0_ark + cos(3.0_ark*rhobar)
    !
    s18 = sin(3.0_ark*rhobar)
    !
  end select
  !
  f = 0
  !
  do ioper = 1,36
    !
    ! for xi(18) = 1.0_ark + cos(3.0_ark*rhobar) 
    ! all operations on xi18->xi18
    !
    call ML_symmetry_transformation_C2H6_G36(ioper,18,xi,chi(:,ioper),s18,d18)
    !
  enddo
  ! 
  do i = 7, molec%parmax
    !
    ipower(1:18) = molec%pot_ind(1:18,i)
    !
    term = 0 
    !
    do ioper = 1,36
      !
      term = term + product(chi(1:18,ioper)**ipower(1:18))
      !
    end do
    !
    ! parameter 6 contains the number of even expansion terms. All parameters after are odd. 
    !
    if (i>int(force(6))+6) term = term*d18
    !
    term = term/36.0_ark
    !
    f = f + term*force(i)
    !
  enddo
  !  
end function MLpoten_c2h6_88_cos3tau_sin3tau_G36



subroutine coordinate_transformation(ncoords,nmodes,local,xi)

  integer(ik),intent(in) :: ncoords,nmodes 
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(out)   :: xi(nmodes)
  !
  real(ark) :: rhobar
  real(ark) :: tau14,tau24,tau25,tau35,tau36,theta12,theta23,theta13,theta56,theta45,theta46
  real(ark) :: tau41,tau51,tau52,tau62,tau63,theta31,theta64,tau16,tau53,tau34,tau26
  !
  select case(trim(molec%coords_transform))
    !
  case default
    !
    write(out, '(/a,1x,a,1x,a)') &
    'MLpoten_c2h6_88 error', trim(molec%coords_transform), 'is unknown'
    stop 'MLpoten_c2h6_88 error error: bad coordinate type'
    !
  case('R-R16-BETA16-THETA-TAU-3','R-R16-BETA16-THETA-TAU-11','R-R16-BETA16-THETA-TAU-17','R-R16-BETA16-THETA-TAU-18')
    !
    xi(1)=local(1)
    xi(2)=local(2)
    xi(3)=local(4)
    xi(4)=local(6)
    xi(5)=local(3)
    xi(6)=local(5)
    xi(7)=local(7)
    !
    xi(8)  = local(8) 
    xi(9)  = local(10)
    xi(10) = local(12)
    xi(11) = local(9) 
    xi(12) = local(11)
    xi(13) = local(13)
    !
    tau14 = mod(local(14)+4.0_ark*pi,4.0_ark*pi)
    tau24 = mod(local(15)+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod(local(16)+2.0_ark*pi,2.0_ark*pi)
    tau35 = mod(local(17)+2.0_ark*pi,2.0_ark*pi)
    tau36 = mod(local(18)+2.0_ark*pi,2.0_ark*pi)
    !
    if (tau14>2.0_ark*pi) then 
       tau25 = tau25 + 2.0_ark*pi
       tau36 = tau36 + 2.0_ark*pi
    endif
    !
    rhobar  = ( tau14+tau25+tau36 )/(3.0_ark)
    !
    tau14 = mod(local(14)+2.0_ark*pi,2.0_ark*pi)
    tau24 = mod(local(15)+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod(local(16)+2.0_ark*pi,2.0_ark*pi)
    tau35 = mod(local(17)+2.0_ark*pi,2.0_ark*pi)
    tau36 = mod(local(18)+2.0_ark*pi,2.0_ark*pi)
    !
    theta12 = mod(tau14-tau24+2.0_ark*pi,2.0_ark*pi)
    theta23 = mod(tau25-tau35+2.0_ark*pi,2.0_ark*pi)
    theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
    !
    theta56 = mod(tau36-tau35+2.0_ark*pi,2.0_ark*pi)
    theta45 = mod(tau25-tau24+2.0_ark*pi,2.0_ark*pi)
    theta46 = mod(2.0_ark*pi-theta56-theta45+2.0_ark*pi,2.0_ark*pi)
    !
    xi(14)  = ( 2.0_ark*theta23 - theta13 - theta12 )/sqrt(6.0_ark)
    xi(15)  = (                   theta13 - theta12 )/sqrt(2.0_ark)
    xi(16)  = ( 2.0_ark*theta56 - theta45 - theta46 )/sqrt(6.0_ark)
    xi(17)  = (                   theta46 - theta45 )/sqrt(2.0_ark)
    !
    !rhobar = ( tau14+tau25+tau36 )/3.0_ark
    !
    xi(18) = rhobar
    !
  case('R-R16-BETA16-THETA-TAU-1')
    !
    xi(1)=local(1)
    xi(2)=local(2)
    xi(3)=local(6)
    xi(4)=local(4)
    xi(5)=local(3)
    xi(6)=local(5)
    xi(7)=local(7)
    !
    xi(8)  = local(8) 
    xi(9)  = local(12)
    xi(10) = local(10)
    xi(11) = local(9) 
    xi(12) = local(11)
    xi(13) = local(13)
    !
    tau14 = mod(local(14)+4.0_ark*pi,4.0_ark*pi)
    tau24 = mod(local(15)+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod(local(16)+2.0_ark*pi,2.0_ark*pi)
    tau35 = mod(local(17)+2.0_ark*pi,2.0_ark*pi)
    tau36 = mod(local(18)+2.0_ark*pi,2.0_ark*pi)
    !
    if (tau14>2.0_ark*pi) then 
       tau25 = tau25 + 2.0_ark*pi
       tau36 = tau36 + 2.0_ark*pi
    endif
    !
    !rhobar  = ( tau14+tau25+tau36 )/(3.0_ark)
    !
    tau14 = mod(local(14)+2.0_ark*pi,2.0_ark*pi)
    tau34 = mod(local(15)+2.0_ark*pi,2.0_ark*pi)
    tau36 = mod(local(16)+2.0_ark*pi,2.0_ark*pi)
    tau26 = mod(local(17)+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod(local(18)+2.0_ark*pi,2.0_ark*pi)
    !
    theta13 = mod(tau14-tau34+2.0_ark*pi,2.0_ark*pi)
    theta23 = mod(tau36-tau26+2.0_ark*pi,2.0_ark*pi)
    theta12 = mod(2.0_ark*pi-theta13-theta23+2.0_ark*pi,2.0_ark*pi)
    !
    theta46 = mod(tau36-tau34+2.0_ark*pi,2.0_ark*pi)
    theta56 = mod(tau25-tau26+2.0_ark*pi,2.0_ark*pi)
    theta45 = mod(2.0_ark*pi-theta46-theta56+2.0_ark*pi,2.0_ark*pi)
    !
    xi(14)  = ( 2.0_ark*theta23 - theta13 - theta12 )/sqrt(6.0_ark)
    xi(15)  = (                   theta13 - theta12 )/sqrt(2.0_ark)
    xi(16)  = ( 2.0_ark*theta56 - theta45 - theta46 )/sqrt(6.0_ark)
    xi(17)  = (                   theta45 - theta46 )/sqrt(2.0_ark)
    !
    rhobar = ( tau14+tau25+tau36 )/3.0_ark
    !
    xi(18) = rhobar
    !
  case('R-R16-BETA16-THETA-TAU-2')
    !
    xi(1)=local(1)
    xi(2)=local(2)
    xi(3)=local(6)
    xi(4)=local(4)
    xi(5)=local(3)
    xi(6)=local(5)
    xi(7)=local(7)
    !
    xi(8)  = local(8) 
    xi(9)  = local(12)
    xi(10) = local(10)
    xi(11) = local(9) 
    xi(12) = local(11)
    xi(13) = local(13)
    !
    tau14 = mod( local(14)+4.0_ark*pi,4.0_ark*pi)
    tau24 = mod( local(15)+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod( local(16)+2.0_ark*pi,2.0_ark*pi)
    tau35 = mod( local(17)+2.0_ark*pi,2.0_ark*pi)
    tau36 = mod( local(18)+2.0_ark*pi,2.0_ark*pi)
    !
    if (tau14>2.0_ark*pi) then 
       tau25 = tau25 + 2.0_ark*pi
       tau36 = tau36 + 2.0_ark*pi
    endif
    !
    !rhobar  = ( tau14+tau25+tau36 )/(3.0_ark)
    !
    tau14 = mod(local(14)+2.0_ark*pi,2.0_ark*pi)
    tau34 = mod(local(15)+2.0_ark*pi,2.0_ark*pi)
    tau36 = mod(local(16)+2.0_ark*pi,2.0_ark*pi)
    tau26 = mod(local(17)+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod(local(18)+2.0_ark*pi,2.0_ark*pi)
    !
    theta13 = mod(tau14-tau34+2.0_ark*pi,2.0_ark*pi)
    theta23 = mod(tau36-tau26+2.0_ark*pi,2.0_ark*pi)
    theta12 = mod(2.0_ark*pi-theta13-theta23+2.0_ark*pi,2.0_ark*pi)
    !
    theta46 = mod(tau36-tau34+2.0_ark*pi,2.0_ark*pi)
    theta56 = mod(tau25-tau26+2.0_ark*pi,2.0_ark*pi)
    theta45 = mod(2.0_ark*pi-theta46-theta56+2.0_ark*pi,2.0_ark*pi)
    !
    xi(14)  = ( 2.0_ark*theta23 - theta12 - theta13 )/sqrt(6.0_ark)
    xi(15)  = (                   theta12 - theta13 )/sqrt(2.0_ark)
    xi(16)  = ( 2.0_ark*theta56 - theta46 - theta45 )/sqrt(6.0_ark)
    xi(17)  = (                   theta46 - theta45 )/sqrt(2.0_ark)
    !
    rhobar = ( tau14+tau25+tau36 )/3.0_ark
    !
    xi(18) = rhobar
    !
  case('R-R16-BETA16-THETA-TAU-12','R-R16-BETA16-THETA-TAU-13','R-R16-BETA16-THETA-TAU-14',&
       'R-R16-BETA16-THETA-TAU-16','R-R16-BETA16-THETA-TAU-19')
    !
    xi(1)=local(1)
    xi(2)=local(3)
    xi(3)=local(5)
    xi(4)=local(7)
    xi(5)=local(2)
    xi(6)=local(6)
    xi(7)=local(4)
    !
    xi(8)  = local(9)  
    xi(9)  = local(11) 
    xi(10) = local(13) 
    xi(11) = local(8)  
    xi(12) = local(12) 
    xi(13) = local(10) 
    !
    tau41 = mod(local(14)+4.0_ark*pi,4.0_ark*pi)
    tau16 = mod(local(15)+2.0_ark*pi,2.0_ark*pi)
    tau62 = mod(local(16)+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod(local(17)+2.0_ark*pi,2.0_ark*pi)
    tau53 = mod(local(18)+2.0_ark*pi,2.0_ark*pi)
    !
    ! assuming this is the 404-type (0..720) for tau14, tau25 and tau36 are extended to 0-720 as well
    if (tau41>2.0_ark*pi) then 
       tau53 = tau53 + 2.0_ark*pi
       tau62 = tau62 + 2.0_ark*pi
    endif
    !
    rhobar  = ( tau41+tau53+tau62)/(3.0_ark)
    !
    !
    tau41 = mod(tau41+2.0_ark*pi,2.0_ark*pi)
    tau16 = mod(tau16+2.0_ark*pi,2.0_ark*pi)
    tau62 = mod(tau62+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod(tau25+2.0_ark*pi,2.0_ark*pi)
    tau53 = mod(tau53+2.0_ark*pi,2.0_ark*pi)
    !
    theta12 = mod(tau62-tau16+2.0_ark*pi,2.0_ark*pi)
    theta23 = mod(tau53-tau25+2.0_ark*pi,2.0_ark*pi)
    theta31 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
    !
    theta56 = mod(tau25-tau62+2.0_ark*pi,2.0_ark*pi)
    theta64 = mod(tau16-tau41+2.0_ark*pi,2.0_ark*pi)
    theta45 = mod(2.0_ark*pi-theta56-theta64+2.0_ark*pi,2.0_ark*pi)
    !
    xi(14)  = ( 2.0_ark*theta23 - theta31 - theta12 )/sqrt(6.0_ark)
    xi(15)  = (                   theta31 - theta12 )/sqrt(2.0_ark)
    !
    xi(16)  = ( 2.0_ark*theta56 - theta64 - theta45 )/sqrt(6.0_ark)
    xi(17)  = (                   theta64 - theta45 )/sqrt(2.0_ark)
    !
    xi(18) = rhobar
    !  
  case('R-R16-BETA16-THETA-TAU')
    !
    xi(1)=local(1)
    xi(2)=local(2)
    xi(3)=local(4)
    xi(4)=local(6)
    xi(5)=local(3)
    xi(6)=local(5)
    xi(7)=local(7)
    !
    tau14 = local(14)
    tau24 = local(15)
    tau25 = local(16)
    tau35 = local(17)
    tau36 = local(18)
    !
    tau14 = mod(tau14+4.0_ark*pi,4.0_ark*pi)
    tau24 = mod(tau24+4.0_ark*pi,4.0_ark*pi)
    tau25 = mod(tau25+4.0_ark*pi,4.0_ark*pi)
    tau35 = mod(tau35+4.0_ark*pi,4.0_ark*pi)
    tau36 = mod(tau36+4.0_ark*pi,4.0_ark*pi)
    !
    rhobar  = ( tau14+tau25+tau36 )/(3.0_ark)-pi
    !
    tau14 = mod(tau14+2.0_ark*pi,2.0_ark*pi)
    tau24 = mod(tau24+2.0_ark*pi,2.0_ark*pi)
    tau25 = mod(tau25+2.0_ark*pi,2.0_ark*pi)
    tau35 = mod(tau35+2.0_ark*pi,2.0_ark*pi)
    tau36 = mod(tau36+2.0_ark*pi,2.0_ark*pi)
    !
    theta12 = mod(tau14-tau24+2.0_ark*pi,2.0_ark*pi)
    theta23 = mod(tau25-tau35+2.0_ark*pi,2.0_ark*pi)
    theta13 = mod(2.0_ark*pi-theta12-theta23+2.0_ark*pi,2.0_ark*pi)
    !
    theta56 = mod(tau36-tau35+2.0_ark*pi,2.0_ark*pi)
    theta45 = mod(tau25-tau24+2.0_ark*pi,2.0_ark*pi)
    theta46 = mod(2.0_ark*pi-theta56-theta45+2.0_ark*pi,2.0_ark*pi)
    !
    xi(14)  = ( 2.0_ark*theta23 - theta13 - theta12 )/sqrt(12.0_ark)
    xi(15)  = (                   theta13 - theta12 )/(2.0_ark)
    xi(16)  = ( 2.0_ark*theta56 - theta46 - theta45 )/sqrt(12.0_ark)
    xi(17)  = (                   theta46 - theta45 )/(2.0_ark)
    !
    xi(18) = rhobar

    xi(8)  = local(8) 
    xi(9)  = local(10)
    xi(10) = local(12)
    xi(11) = local(9) 
    xi(12) = local(11)
    xi(13) = local(13)
    !
  end select
  !
  end subroutine coordinate_transformation


  recursive subroutine MLdms2xyz_c2h6_g36(rank, ncoords, natoms, local, xyz, f)
    !
    ! G_36(EM) symmetry-adapted electric dipole moment surface for C2H6
    ! (70 A_4s parallel z parameters + 109 E_2d perpendicular x,y parameters)
    !
    integer(ik), intent(in) :: rank, ncoords, natoms
    real(ark), intent(in)   :: local(ncoords), xyz(natoms, 3)
    real(ark), intent(out)  :: f(rank)
    !
    real(ark) :: xi(18), r1e, r2e, betae, rad_loc
    real(rk)  :: xi_rk(18), params_z(N_PARAMS_Z), params_xy(N_PARAMS_XY), mu_rk(3)
    integer(ik) :: n1, n2
    !
    rad_loc = pi / 180.0_ark
    r1e   = extF%coef(1, 1)
    r2e   = extF%coef(2, 1)
    betae = extF%coef(3, 1) * rad_loc
    !
    call coordinate_transformation(ncoords, 18_ik, local, xi)
    !
    xi(1)    = xi(1)    - r1e
    xi(2:7)  = xi(2:7)  - r2e
    xi(8:13) = xi(8:13) - betae
    !
    xi_rk(1:18) = real(xi(1:18), kind=rk)
    params_z  = 0.0_rk
    params_xy = 0.0_rk
    n1 = size(extF%coef, 1)
    n2 = size(extF%coef, 2)
    !
    if (n1 >= 3 + N_PARAMS_DMS) then
      params_z(1:N_PARAMS_Z)   = real(extF%coef(4:3+N_PARAMS_Z, 1), kind=rk)
      params_xy(1:N_PARAMS_XY) = real(extF%coef(4+N_PARAMS_Z:3+N_PARAMS_DMS, 1), kind=rk)
    else
      if (n1 >= 3 + N_PARAMS_XY) then
        params_xy(1:N_PARAMS_XY) = real(extF%coef(4:3+N_PARAMS_XY, 1), kind=rk)
      endif
      if (n2 >= 3 .and. n1 >= 3 + N_PARAMS_Z) then
        params_z(1:N_PARAMS_Z) = real(extF%coef(4:3+N_PARAMS_Z, 3), kind=rk)
      endif
    endif
    !
    call eval_c2h6_dms(xi_rk, params_z, params_xy, mu_rk)
    !
    f = 0.0_ark
    f(1) = real(mu_rk(1), kind=ark)
    f(2) = real(mu_rk(2), kind=ark)
    f(3) = real(mu_rk(3), kind=ark)
    !
  end subroutine MLdms2xyz_c2h6_g36


  ! =============================================================================
  ! Automatically generated G_36(EM) symmetry-adapted DMS for C2H6 (DMS_C2H6_G36),
  ! formerly the standalone module dms_c2h6_g36 (dms_c2h6_g36.f90)
  ! Frame: KEO Half-Torsion Body-Fixed Frame (ML_b0_C2H6 / A2B6.nb / disp_to_cartesian)
  ! Coordinates: xi(1..17) = vibrational displacements, xi(18) = tau (radians)
  ! real(rk) is the double-precision kind rk of the module accuracy (the same kind as
  ! selected_real_kind(15, 307) of the generated module)
  ! =============================================================================

  subroutine eval_c2h6_dms_basis(xi, phi_z, phi_x, phi_y)
    real(rk), intent(in)  :: xi(18)
    real(rk), intent(out) :: phi_z(N_PARAMS_Z)
    real(rk), intent(out) :: phi_x(N_PARAMS_XY), phi_y(N_PARAMS_XY)
    real(rk) :: tau, half_tau
    real(rk) :: c_1, s_1
    real(rk) :: c_2, s_2
    real(rk) :: c_3, s_3
    real(rk) :: c_4, s_4
    real(rk) :: c_5, s_5
    real(rk) :: c_6, s_6

    tau = xi(18)
    half_tau = 0.5_rk * tau
    c_1 = cos(1.0_rk * half_tau)
    s_1 = sin(1.0_rk * half_tau)
    c_2 = cos(2.0_rk * half_tau)
    s_2 = sin(2.0_rk * half_tau)
    c_3 = cos(3.0_rk * half_tau)
    s_3 = sin(3.0_rk * half_tau)
    c_4 = cos(4.0_rk * half_tau)
    s_4 = sin(4.0_rk * half_tau)
    c_5 = cos(5.0_rk * half_tau)
    s_5 = sin(5.0_rk * half_tau)
    c_6 = cos(6.0_rk * half_tau)
    s_6 = sin(6.0_rk * half_tau)

    ! --- Parallel dipole basis functions phi_z (A_4s) ---
    phi_z(1) = (+1.000000000000000e+00_rk) * xi(7) &
         + (+1.000000000000000e+00_rk) * xi(6) &
         + (+1.000000000000000e+00_rk) * xi(5) &
         + (-1.000000000000000e+00_rk) * xi(4) &
         + (-1.000000000000000e+00_rk) * xi(3) &
         + (-1.000000000000000e+00_rk) * xi(2)
    phi_z(2) = (+1.000000000000000e+00_rk) * xi(13) &
         + (+1.000000000000000e+00_rk) * xi(12) &
         + (+1.000000000000000e+00_rk) * xi(11) &
         + (-1.000000000000000e+00_rk) * xi(10) &
         + (-1.000000000000000e+00_rk) * xi(9) &
         + (-1.000000000000000e+00_rk) * xi(8)
    phi_z(3) = (+1.000000000000000e+00_rk) * xi(7) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(6) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(5) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(4) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(3) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(2) * c_6
    phi_z(4) = (+1.000000000000000e+00_rk) * xi(13) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(12) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(11) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(10) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(9) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(8) * c_6
    phi_z(5) = (+1.000000000000000e+00_rk) * xi(7)**2 &
         + (+1.000000000000000e+00_rk) * xi(6)**2 &
         + (+1.000000000000000e+00_rk) * xi(5)**2 &
         + (-1.000000000000000e+00_rk) * xi(4)**2 &
         + (-1.000000000000000e+00_rk) * xi(3)**2 &
         + (-1.000000000000000e+00_rk) * xi(2)**2
    phi_z(6) = (+1.000000000000000e+00_rk) * xi(13)**2 &
         + (+1.000000000000000e+00_rk) * xi(12)**2 &
         + (+1.000000000000000e+00_rk) * xi(11)**2 &
         + (-1.000000000000000e+00_rk) * xi(10)**2 &
         + (-1.000000000000000e+00_rk) * xi(9)**2 &
         + (-1.000000000000000e+00_rk) * xi(8)**2
    phi_z(7) = (+1.000000000000000e+00_rk) * xi(17)**2 &
         + (+1.000000000000000e+00_rk) * xi(16)**2 &
         + (-1.000000000000000e+00_rk) * xi(15)**2 &
         + (-1.000000000000000e+00_rk) * xi(14)**2
    phi_z(8) = (+1.000000000000000e+00_rk) * xi(1) * xi(7) &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(6) &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(5) &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(4) &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(3) &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(2)
    phi_z(9) = (+1.000000000000000e+00_rk) * xi(1) * xi(13) &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(12) &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(11) &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(10) &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(9) &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(8)
    phi_z(10) = (+1.000000000000000e+00_rk) * xi(6) * xi(7) &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(7) &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(6) &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(4) &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(4) &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(3)
    phi_z(11) = (+1.000000000000000e+00_rk) * xi(12) * xi(13) &
         + (+1.000000000000000e+00_rk) * xi(11) * xi(13) &
         + (+1.000000000000000e+00_rk) * xi(11) * xi(12) &
         + (-1.000000000000000e+00_rk) * xi(9) * xi(10) &
         + (-1.000000000000000e+00_rk) * xi(8) * xi(10) &
         + (-1.000000000000000e+00_rk) * xi(8) * xi(9)
    phi_z(12) = (+1.000000000000000e+00_rk) * xi(7) * xi(13) &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(12) &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(11) &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(10) &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(9) &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(8)
    phi_z(13) = (+1.000000000000000e+00_rk) * xi(7) * xi(12) &
         + (+1.000000000000000e+00_rk) * xi(7) * xi(11) &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(13) &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(11) &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(13) &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(12) &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(9) &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(8) &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(10) &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(8) &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(10) &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(9)
    phi_z(14) = (+1.000000000000000e+00_rk) * xi(7) * xi(10) &
         + (+1.000000000000000e+00_rk) * xi(7) * xi(9) &
         + (+1.000000000000000e+00_rk) * xi(7) * xi(8) &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(10) &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(9) &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(8) &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(10) &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(9) &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(8) &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(13) &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(12) &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(11) &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(13) &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(12) &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(11) &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(13) &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(12) &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(11)
    phi_z(15) = (-8.660254037844385e-01_rk) * xi(7) * xi(17) &
         + (-4.999999999999999e-01_rk) * xi(7) * xi(16) &
         + (+8.660254037844385e-01_rk) * xi(6) * xi(17) &
         + (-5.000000000000001e-01_rk) * xi(6) * xi(16) &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(16) &
         + (+8.660254037844385e-01_rk) * xi(4) * xi(15) &
         + (+5.000000000000001e-01_rk) * xi(4) * xi(14) &
         + (-8.660254037844386e-01_rk) * xi(3) * xi(15) &
         + (+4.999999999999999e-01_rk) * xi(3) * xi(14) &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(14)
    phi_z(16) = (-8.660254037844385e-01_rk) * xi(13) * xi(17) &
         + (-4.999999999999999e-01_rk) * xi(13) * xi(16) &
         + (+8.660254037844385e-01_rk) * xi(12) * xi(17) &
         + (-5.000000000000001e-01_rk) * xi(12) * xi(16) &
         + (+1.000000000000000e+00_rk) * xi(11) * xi(16) &
         + (+8.660254037844385e-01_rk) * xi(10) * xi(15) &
         + (+5.000000000000001e-01_rk) * xi(10) * xi(14) &
         + (-8.660254037844386e-01_rk) * xi(9) * xi(15) &
         + (+4.999999999999999e-01_rk) * xi(9) * xi(14) &
         + (-1.000000000000000e+00_rk) * xi(8) * xi(14)
    phi_z(17) = (+1.000000000000000e+00_rk) * xi(7) * xi(10) * c_2 &
         + (-4.999999999999997e-01_rk) * xi(7) * xi(9) * c_2 &
         + (-8.660254037844387e-01_rk) * xi(7) * xi(9) * s_2 &
         + (-5.000000000000003e-01_rk) * xi(7) * xi(8) * c_2 &
         + (+8.660254037844385e-01_rk) * xi(7) * xi(8) * s_2 &
         + (-5.000000000000003e-01_rk) * xi(6) * xi(10) * c_2 &
         + (+8.660254037844385e-01_rk) * xi(6) * xi(10) * s_2 &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(9) * c_2 &
         + (-4.999999999999997e-01_rk) * xi(6) * xi(8) * c_2 &
         + (-8.660254037844387e-01_rk) * xi(6) * xi(8) * s_2 &
         + (-4.999999999999997e-01_rk) * xi(5) * xi(10) * c_2 &
         + (-8.660254037844387e-01_rk) * xi(5) * xi(10) * s_2 &
         + (-5.000000000000003e-01_rk) * xi(5) * xi(9) * c_2 &
         + (+8.660254037844385e-01_rk) * xi(5) * xi(9) * s_2 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(8) * c_2 &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(13) * c_2 &
         + (+5.000000000000003e-01_rk) * xi(4) * xi(12) * c_2 &
         + (-8.660254037844385e-01_rk) * xi(4) * xi(12) * s_2 &
         + (+4.999999999999997e-01_rk) * xi(4) * xi(11) * c_2 &
         + (+8.660254037844387e-01_rk) * xi(4) * xi(11) * s_2 &
         + (+4.999999999999994e-01_rk) * xi(3) * xi(13) * c_2 &
         + (+8.660254037844388e-01_rk) * xi(3) * xi(13) * s_2 &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(12) * c_2 &
         + (+5.000000000000003e-01_rk) * xi(3) * xi(11) * c_2 &
         + (-8.660254037844385e-01_rk) * xi(3) * xi(11) * s_2 &
         + (+5.000000000000003e-01_rk) * xi(2) * xi(13) * c_2 &
         + (-8.660254037844385e-01_rk) * xi(2) * xi(13) * s_2 &
         + (+4.999999999999997e-01_rk) * xi(2) * xi(12) * c_2 &
         + (+8.660254037844387e-01_rk) * xi(2) * xi(12) * s_2 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(11) * c_2
    phi_z(18) = (-8.660254037844385e-01_rk) * xi(7) * xi(15) * c_2 &
         + (-5.000000000000001e-01_rk) * xi(7) * xi(15) * s_2 &
         + (-5.000000000000002e-01_rk) * xi(7) * xi(14) * c_2 &
         + (+8.660254037844383e-01_rk) * xi(7) * xi(14) * s_2 &
         + (+8.660254037844386e-01_rk) * xi(6) * xi(15) * c_2 &
         + (-4.999999999999998e-01_rk) * xi(6) * xi(15) * s_2 &
         + (-4.999999999999996e-01_rk) * xi(6) * xi(14) * c_2 &
         + (-8.660254037844386e-01_rk) * xi(6) * xi(14) * s_2 &
         + (+9.999999999999998e-01_rk) * xi(5) * xi(15) * s_2 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(14) * c_2 &
         + (+8.660254037844386e-01_rk) * xi(4) * xi(17) * c_2 &
         + (-4.999999999999998e-01_rk) * xi(4) * xi(17) * s_2 &
         + (+4.999999999999996e-01_rk) * xi(4) * xi(16) * c_2 &
         + (+8.660254037844386e-01_rk) * xi(4) * xi(16) * s_2 &
         + (-8.660254037844383e-01_rk) * xi(3) * xi(17) * c_2 &
         + (-5.000000000000002e-01_rk) * xi(3) * xi(17) * s_2 &
         + (+5.000000000000003e-01_rk) * xi(3) * xi(16) * c_2 &
         + (-8.660254037844385e-01_rk) * xi(3) * xi(16) * s_2 &
         + (+9.999999999999998e-01_rk) * xi(2) * xi(17) * s_2 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(16) * c_2
    phi_z(19) = (-8.660254037844385e-01_rk) * xi(13) * xi(15) * c_2 &
         + (-5.000000000000001e-01_rk) * xi(13) * xi(15) * s_2 &
         + (-5.000000000000002e-01_rk) * xi(13) * xi(14) * c_2 &
         + (+8.660254037844383e-01_rk) * xi(13) * xi(14) * s_2 &
         + (+8.660254037844386e-01_rk) * xi(12) * xi(15) * c_2 &
         + (-4.999999999999998e-01_rk) * xi(12) * xi(15) * s_2 &
         + (-4.999999999999996e-01_rk) * xi(12) * xi(14) * c_2 &
         + (-8.660254037844386e-01_rk) * xi(12) * xi(14) * s_2 &
         + (+9.999999999999998e-01_rk) * xi(11) * xi(15) * s_2 &
         + (+1.000000000000000e+00_rk) * xi(11) * xi(14) * c_2 &
         + (+8.660254037844386e-01_rk) * xi(10) * xi(17) * c_2 &
         + (-4.999999999999998e-01_rk) * xi(10) * xi(17) * s_2 &
         + (+4.999999999999996e-01_rk) * xi(10) * xi(16) * c_2 &
         + (+8.660254037844386e-01_rk) * xi(10) * xi(16) * s_2 &
         + (-8.660254037844383e-01_rk) * xi(9) * xi(17) * c_2 &
         + (-5.000000000000002e-01_rk) * xi(9) * xi(17) * s_2 &
         + (+5.000000000000003e-01_rk) * xi(9) * xi(16) * c_2 &
         + (-8.660254037844385e-01_rk) * xi(9) * xi(16) * s_2 &
         + (+9.999999999999998e-01_rk) * xi(8) * xi(17) * s_2 &
         + (-1.000000000000000e+00_rk) * xi(8) * xi(16) * c_2
    phi_z(20) = (+1.000000000000000e+00_rk) * xi(7) * xi(10) * c_4 &
         + (-5.000000000000006e-01_rk) * xi(7) * xi(9) * c_4 &
         + (+8.660254037844385e-01_rk) * xi(7) * xi(9) * s_4 &
         + (-4.999999999999994e-01_rk) * xi(7) * xi(8) * c_4 &
         + (-8.660254037844389e-01_rk) * xi(7) * xi(8) * s_4 &
         + (-4.999999999999994e-01_rk) * xi(6) * xi(10) * c_4 &
         + (-8.660254037844389e-01_rk) * xi(6) * xi(10) * s_4 &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(9) * c_4 &
         + (-5.000000000000006e-01_rk) * xi(6) * xi(8) * c_4 &
         + (+8.660254037844385e-01_rk) * xi(6) * xi(8) * s_4 &
         + (-5.000000000000006e-01_rk) * xi(5) * xi(10) * c_4 &
         + (+8.660254037844385e-01_rk) * xi(5) * xi(10) * s_4 &
         + (-4.999999999999994e-01_rk) * xi(5) * xi(9) * c_4 &
         + (-8.660254037844389e-01_rk) * xi(5) * xi(9) * s_4 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(8) * c_4 &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(13) * c_4 &
         + (+4.999999999999994e-01_rk) * xi(4) * xi(12) * c_4 &
         + (+8.660254037844389e-01_rk) * xi(4) * xi(12) * s_4 &
         + (+5.000000000000006e-01_rk) * xi(4) * xi(11) * c_4 &
         + (-8.660254037844385e-01_rk) * xi(4) * xi(11) * s_4 &
         + (+5.000000000000012e-01_rk) * xi(3) * xi(13) * c_4 &
         + (-8.660254037844379e-01_rk) * xi(3) * xi(13) * s_4 &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(12) * c_4 &
         + (+4.999999999999994e-01_rk) * xi(3) * xi(11) * c_4 &
         + (+8.660254037844389e-01_rk) * xi(3) * xi(11) * s_4 &
         + (+4.999999999999994e-01_rk) * xi(2) * xi(13) * c_4 &
         + (+8.660254037844389e-01_rk) * xi(2) * xi(13) * s_4 &
         + (+5.000000000000006e-01_rk) * xi(2) * xi(12) * c_4 &
         + (-8.660254037844385e-01_rk) * xi(2) * xi(12) * s_4 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(11) * c_4
    phi_z(21) = (+8.660254037844389e-01_rk) * xi(7) * xi(15) * c_4 &
         + (-4.999999999999997e-01_rk) * xi(7) * xi(15) * s_4 &
         + (+4.999999999999994e-01_rk) * xi(7) * xi(14) * c_4 &
         + (+8.660254037844388e-01_rk) * xi(7) * xi(14) * s_4 &
         + (-8.660254037844382e-01_rk) * xi(6) * xi(15) * c_4 &
         + (-5.000000000000002e-01_rk) * xi(6) * xi(15) * s_4 &
         + (+5.000000000000006e-01_rk) * xi(6) * xi(14) * c_4 &
         + (-8.660254037844385e-01_rk) * xi(6) * xi(14) * s_4 &
         + (+9.999999999999998e-01_rk) * xi(5) * xi(15) * s_4 &
         + (-1.000000000000000e+00_rk) * xi(5) * xi(14) * c_4 &
         + (-8.660254037844383e-01_rk) * xi(4) * xi(17) * c_4 &
         + (-5.000000000000002e-01_rk) * xi(4) * xi(17) * s_4 &
         + (-5.000000000000006e-01_rk) * xi(4) * xi(16) * c_4 &
         + (+8.660254037844385e-01_rk) * xi(4) * xi(16) * s_4 &
         + (+8.660254037844393e-01_rk) * xi(3) * xi(17) * c_4 &
         + (-4.999999999999996e-01_rk) * xi(3) * xi(17) * s_4 &
         + (-4.999999999999992e-01_rk) * xi(3) * xi(16) * c_4 &
         + (-8.660254037844386e-01_rk) * xi(3) * xi(16) * s_4 &
         + (+9.999999999999998e-01_rk) * xi(2) * xi(17) * s_4 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(16) * c_4
    phi_z(22) = (+8.660254037844389e-01_rk) * xi(13) * xi(15) * c_4 &
         + (-4.999999999999997e-01_rk) * xi(13) * xi(15) * s_4 &
         + (+4.999999999999994e-01_rk) * xi(13) * xi(14) * c_4 &
         + (+8.660254037844388e-01_rk) * xi(13) * xi(14) * s_4 &
         + (-8.660254037844382e-01_rk) * xi(12) * xi(15) * c_4 &
         + (-5.000000000000002e-01_rk) * xi(12) * xi(15) * s_4 &
         + (+5.000000000000006e-01_rk) * xi(12) * xi(14) * c_4 &
         + (-8.660254037844385e-01_rk) * xi(12) * xi(14) * s_4 &
         + (+9.999999999999998e-01_rk) * xi(11) * xi(15) * s_4 &
         + (-1.000000000000000e+00_rk) * xi(11) * xi(14) * c_4 &
         + (-8.660254037844383e-01_rk) * xi(10) * xi(17) * c_4 &
         + (-5.000000000000002e-01_rk) * xi(10) * xi(17) * s_4 &
         + (-5.000000000000006e-01_rk) * xi(10) * xi(16) * c_4 &
         + (+8.660254037844385e-01_rk) * xi(10) * xi(16) * s_4 &
         + (+8.660254037844393e-01_rk) * xi(9) * xi(17) * c_4 &
         + (-4.999999999999996e-01_rk) * xi(9) * xi(17) * s_4 &
         + (-4.999999999999992e-01_rk) * xi(9) * xi(16) * c_4 &
         + (-8.660254037844386e-01_rk) * xi(9) * xi(16) * s_4 &
         + (+9.999999999999998e-01_rk) * xi(8) * xi(17) * s_4 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(16) * c_4
    phi_z(23) = (+1.000000000000000e+00_rk) * xi(7)**2 * c_6 &
         + (+1.000000000000000e+00_rk) * xi(6)**2 * c_6 &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * c_6 &
         + (-1.000000000000000e+00_rk) * xi(4)**2 * c_6 &
         + (-1.000000000000000e+00_rk) * xi(3)**2 * c_6 &
         + (-1.000000000000000e+00_rk) * xi(2)**2 * c_6
    phi_z(24) = (+1.000000000000000e+00_rk) * xi(13)**2 * c_6 &
         + (+1.000000000000000e+00_rk) * xi(12)**2 * c_6 &
         + (+1.000000000000000e+00_rk) * xi(11)**2 * c_6 &
         + (-1.000000000000000e+00_rk) * xi(10)**2 * c_6 &
         + (-1.000000000000000e+00_rk) * xi(9)**2 * c_6 &
         + (-1.000000000000000e+00_rk) * xi(8)**2 * c_6
    phi_z(25) = (+1.000000000000000e+00_rk) * xi(17)**2 * c_6 &
         + (+1.000000000000000e+00_rk) * xi(16)**2 * c_6 &
         + (-1.000000000000000e+00_rk) * xi(15)**2 * c_6 &
         + (-1.000000000000000e+00_rk) * xi(14)**2 * c_6
    phi_z(26) = (+1.000000000000000e+00_rk) * xi(1) * xi(7) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(6) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(5) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(4) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(3) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(2) * c_6
    phi_z(27) = (+1.000000000000000e+00_rk) * xi(1) * xi(13) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(12) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(11) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(10) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(9) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(8) * c_6
    phi_z(28) = (+1.000000000000000e+00_rk) * xi(6) * xi(7) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(7) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(6) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(4) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(4) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(3) * c_6
    phi_z(29) = (+1.000000000000000e+00_rk) * xi(12) * xi(13) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(11) * xi(13) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(11) * xi(12) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(9) * xi(10) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(8) * xi(10) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(8) * xi(9) * c_6
    phi_z(30) = (+1.000000000000000e+00_rk) * xi(7) * xi(13) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(12) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(11) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(10) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(9) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(8) * c_6
    phi_z(31) = (+1.000000000000000e+00_rk) * xi(7) * xi(12) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(7) * xi(11) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(13) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(11) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(13) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(12) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(9) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(8) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(10) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(8) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(10) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(9) * c_6
    phi_z(32) = (+1.000000000000000e+00_rk) * xi(7) * xi(12) * s_6 &
         + (-1.000000000000000e+00_rk) * xi(7) * xi(11) * s_6 &
         + (-1.000000000000000e+00_rk) * xi(6) * xi(13) * s_6 &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(11) * s_6 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(13) * s_6 &
         + (-1.000000000000000e+00_rk) * xi(5) * xi(12) * s_6 &
         + (+1.000000000000000e+00_rk) * xi(4) * xi(9) * s_6 &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(8) * s_6 &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(10) * s_6 &
         + (+1.000000000000000e+00_rk) * xi(3) * xi(8) * s_6 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(10) * s_6 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(9) * s_6
    phi_z(33) = (+1.000000000000000e+00_rk) * xi(7) * xi(10) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(7) * xi(9) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(7) * xi(8) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(10) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(9) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(8) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(10) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(9) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(8) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(13) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(12) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(11) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(13) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(12) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(11) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(13) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(12) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(11) * c_6
    phi_z(34) = (-8.660254037844385e-01_rk) * xi(7) * xi(17) * c_6 &
         + (-4.999999999999999e-01_rk) * xi(7) * xi(16) * c_6 &
         + (+8.660254037844385e-01_rk) * xi(6) * xi(17) * c_6 &
         + (-5.000000000000001e-01_rk) * xi(6) * xi(16) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(16) * c_6 &
         + (+8.660254037844385e-01_rk) * xi(4) * xi(15) * c_6 &
         + (+5.000000000000001e-01_rk) * xi(4) * xi(14) * c_6 &
         + (-8.660254037844386e-01_rk) * xi(3) * xi(15) * c_6 &
         + (+4.999999999999999e-01_rk) * xi(3) * xi(14) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(14) * c_6
    phi_z(35) = (-5.000000000000000e-01_rk) * xi(7) * xi(17) * s_6 &
         + (+8.660254037844386e-01_rk) * xi(7) * xi(16) * s_6 &
         + (-5.000000000000001e-01_rk) * xi(6) * xi(17) * s_6 &
         + (-8.660254037844386e-01_rk) * xi(6) * xi(16) * s_6 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(17) * s_6 &
         + (-5.000000000000001e-01_rk) * xi(4) * xi(15) * s_6 &
         + (+8.660254037844386e-01_rk) * xi(4) * xi(14) * s_6 &
         + (-5.000000000000000e-01_rk) * xi(3) * xi(15) * s_6 &
         + (-8.660254037844388e-01_rk) * xi(3) * xi(14) * s_6 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(15) * s_6
    phi_z(36) = (-8.660254037844385e-01_rk) * xi(13) * xi(17) * c_6 &
         + (-4.999999999999999e-01_rk) * xi(13) * xi(16) * c_6 &
         + (+8.660254037844385e-01_rk) * xi(12) * xi(17) * c_6 &
         + (-5.000000000000001e-01_rk) * xi(12) * xi(16) * c_6 &
         + (+1.000000000000000e+00_rk) * xi(11) * xi(16) * c_6 &
         + (+8.660254037844385e-01_rk) * xi(10) * xi(15) * c_6 &
         + (+5.000000000000001e-01_rk) * xi(10) * xi(14) * c_6 &
         + (-8.660254037844386e-01_rk) * xi(9) * xi(15) * c_6 &
         + (+4.999999999999999e-01_rk) * xi(9) * xi(14) * c_6 &
         + (-1.000000000000000e+00_rk) * xi(8) * xi(14) * c_6
    phi_z(37) = (-5.000000000000000e-01_rk) * xi(13) * xi(17) * s_6 &
         + (+8.660254037844386e-01_rk) * xi(13) * xi(16) * s_6 &
         + (-5.000000000000001e-01_rk) * xi(12) * xi(17) * s_6 &
         + (-8.660254037844386e-01_rk) * xi(12) * xi(16) * s_6 &
         + (+1.000000000000000e+00_rk) * xi(11) * xi(17) * s_6 &
         + (-5.000000000000001e-01_rk) * xi(10) * xi(15) * s_6 &
         + (+8.660254037844386e-01_rk) * xi(10) * xi(14) * s_6 &
         + (-5.000000000000000e-01_rk) * xi(9) * xi(15) * s_6 &
         + (-8.660254037844388e-01_rk) * xi(9) * xi(14) * s_6 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(15) * s_6
    phi_z(38) = (+1.000000000000000e+00_rk) * xi(7)**3 &
         + (+1.000000000000000e+00_rk) * xi(6)**3 &
         + (+1.000000000000000e+00_rk) * xi(5)**3 &
         + (-1.000000000000000e+00_rk) * xi(4)**3 &
         + (-1.000000000000000e+00_rk) * xi(3)**3 &
         + (-1.000000000000000e+00_rk) * xi(2)**3
    phi_z(39) = (+1.000000000000000e+00_rk) * xi(13)**3 &
         + (+1.000000000000000e+00_rk) * xi(12)**3 &
         + (+1.000000000000000e+00_rk) * xi(11)**3 &
         + (-1.000000000000000e+00_rk) * xi(10)**3 &
         + (-1.000000000000000e+00_rk) * xi(9)**3 &
         + (-1.000000000000000e+00_rk) * xi(8)**3
    phi_z(40) = (+1.000000000000000e+00_rk) * xi(16) * xi(17)**2 &
         + (-3.333333333333334e-01_rk) * xi(16)**3 &
         + (-1.000000000000000e+00_rk) * xi(14) * xi(15)**2 &
         + (+3.333333333333334e-01_rk) * xi(14)**3
    phi_z(41) = (+1.000000000000000e+00_rk) * xi(1)**2 * xi(7) &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(6) &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(5) &
         + (-1.000000000000000e+00_rk) * xi(1)**2 * xi(4) &
         + (-1.000000000000000e+00_rk) * xi(1)**2 * xi(3) &
         + (-1.000000000000000e+00_rk) * xi(1)**2 * xi(2)
    phi_z(42) = (+1.000000000000000e+00_rk) * xi(1) * xi(7)**2 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(6)**2 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(5)**2 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(4)**2 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(3)**2 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(2)**2
    phi_z(43) = (+1.000000000000000e+00_rk) * xi(1)**2 * xi(13) &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(12) &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(11) &
         + (-1.000000000000000e+00_rk) * xi(1)**2 * xi(10) &
         + (-1.000000000000000e+00_rk) * xi(1)**2 * xi(9) &
         + (-1.000000000000000e+00_rk) * xi(1)**2 * xi(8)
    phi_z(44) = (+1.000000000000000e+00_rk) * xi(1) * xi(13)**2 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(12)**2 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(11)**2 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(10)**2 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(9)**2 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(8)**2
    phi_z(45) = (+1.000000000000000e+00_rk) * xi(1) * xi(17)**2 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(16)**2 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(15)**2 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(14)**2
    phi_z(46) = (+1.000000000000000e+00_rk) * xi(6) * xi(7)**2 &
         + (+1.000000000000000e+00_rk) * xi(6)**2 * xi(7) &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(7)**2 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(6)**2 &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * xi(7) &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * xi(6) &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(4)**2 &
         + (-1.000000000000000e+00_rk) * xi(3)**2 * xi(4) &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(4)**2 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(3)**2 &
         + (-1.000000000000000e+00_rk) * xi(2)**2 * xi(4) &
         + (-1.000000000000000e+00_rk) * xi(2)**2 * xi(3)
    phi_z(47) = (+1.000000000000000e+00_rk) * xi(4) * xi(7)**2 &
         + (+1.000000000000000e+00_rk) * xi(4) * xi(6)**2 &
         + (+1.000000000000000e+00_rk) * xi(4) * xi(5)**2 &
         + (-1.000000000000000e+00_rk) * xi(4)**2 * xi(7) &
         + (-1.000000000000000e+00_rk) * xi(4)**2 * xi(6) &
         + (-1.000000000000000e+00_rk) * xi(4)**2 * xi(5) &
         + (+1.000000000000000e+00_rk) * xi(3) * xi(7)**2 &
         + (+1.000000000000000e+00_rk) * xi(3) * xi(6)**2 &
         + (+1.000000000000000e+00_rk) * xi(3) * xi(5)**2 &
         + (-1.000000000000000e+00_rk) * xi(3)**2 * xi(7) &
         + (-1.000000000000000e+00_rk) * xi(3)**2 * xi(6) &
         + (-1.000000000000000e+00_rk) * xi(3)**2 * xi(5) &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(7)**2 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(6)**2 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(5)**2 &
         + (-1.000000000000000e+00_rk) * xi(2)**2 * xi(7) &
         + (-1.000000000000000e+00_rk) * xi(2)**2 * xi(6) &
         + (-1.000000000000000e+00_rk) * xi(2)**2 * xi(5)
    phi_z(48) = (+1.000000000000000e+00_rk) * xi(12) * xi(13)**2 &
         + (+1.000000000000000e+00_rk) * xi(12)**2 * xi(13) &
         + (+1.000000000000000e+00_rk) * xi(11) * xi(13)**2 &
         + (+1.000000000000000e+00_rk) * xi(11) * xi(12)**2 &
         + (+1.000000000000000e+00_rk) * xi(11)**2 * xi(13) &
         + (+1.000000000000000e+00_rk) * xi(11)**2 * xi(12) &
         + (-1.000000000000000e+00_rk) * xi(9) * xi(10)**2 &
         + (-1.000000000000000e+00_rk) * xi(9)**2 * xi(10) &
         + (-1.000000000000000e+00_rk) * xi(8) * xi(10)**2 &
         + (-1.000000000000000e+00_rk) * xi(8) * xi(9)**2 &
         + (-1.000000000000000e+00_rk) * xi(8)**2 * xi(10) &
         + (-1.000000000000000e+00_rk) * xi(8)**2 * xi(9)
    phi_z(49) = (+1.000000000000000e+00_rk) * xi(10) * xi(13)**2 &
         + (+1.000000000000000e+00_rk) * xi(10) * xi(12)**2 &
         + (+1.000000000000000e+00_rk) * xi(10) * xi(11)**2 &
         + (-1.000000000000000e+00_rk) * xi(10)**2 * xi(13) &
         + (-1.000000000000000e+00_rk) * xi(10)**2 * xi(12) &
         + (-1.000000000000000e+00_rk) * xi(10)**2 * xi(11) &
         + (+1.000000000000000e+00_rk) * xi(9) * xi(13)**2 &
         + (+1.000000000000000e+00_rk) * xi(9) * xi(12)**2 &
         + (+1.000000000000000e+00_rk) * xi(9) * xi(11)**2 &
         + (-1.000000000000000e+00_rk) * xi(9)**2 * xi(13) &
         + (-1.000000000000000e+00_rk) * xi(9)**2 * xi(12) &
         + (-1.000000000000000e+00_rk) * xi(9)**2 * xi(11) &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(13)**2 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(12)**2 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(11)**2 &
         + (-1.000000000000000e+00_rk) * xi(8)**2 * xi(13) &
         + (-1.000000000000000e+00_rk) * xi(8)**2 * xi(12) &
         + (-1.000000000000000e+00_rk) * xi(8)**2 * xi(11)
    phi_z(50) = (+1.000000000000000e+00_rk) * xi(7)**2 * xi(13) &
         + (+1.000000000000000e+00_rk) * xi(6)**2 * xi(12) &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * xi(11) &
         + (-1.000000000000000e+00_rk) * xi(4)**2 * xi(10) &
         + (-1.000000000000000e+00_rk) * xi(3)**2 * xi(9) &
         + (-1.000000000000000e+00_rk) * xi(2)**2 * xi(8)
    phi_z(51) = (+1.000000000000000e+00_rk) * xi(7) * xi(13)**2 &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(12)**2 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(11)**2 &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(10)**2 &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(9)**2 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(8)**2
    phi_z(52) = (+1.000000000000000e+00_rk) * xi(7)**2 * xi(12) &
         + (+1.000000000000000e+00_rk) * xi(7)**2 * xi(11) &
         + (+1.000000000000000e+00_rk) * xi(6)**2 * xi(13) &
         + (+1.000000000000000e+00_rk) * xi(6)**2 * xi(11) &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * xi(13) &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * xi(12) &
         + (-1.000000000000000e+00_rk) * xi(4)**2 * xi(9) &
         + (-1.000000000000000e+00_rk) * xi(4)**2 * xi(8) &
         + (-1.000000000000000e+00_rk) * xi(3)**2 * xi(10) &
         + (-1.000000000000000e+00_rk) * xi(3)**2 * xi(8) &
         + (-1.000000000000000e+00_rk) * xi(2)**2 * xi(10) &
         + (-1.000000000000000e+00_rk) * xi(2)**2 * xi(9)
    phi_z(53) = (+1.000000000000000e+00_rk) * xi(7) * xi(12)**2 &
         + (+1.000000000000000e+00_rk) * xi(7) * xi(11)**2 &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(13)**2 &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(11)**2 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(13)**2 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(12)**2 &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(9)**2 &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(8)**2 &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(10)**2 &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(8)**2 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(10)**2 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(9)**2
    phi_z(54) = (-8.660254037844385e-01_rk) * xi(7)**2 * xi(17) &
         + (-4.999999999999999e-01_rk) * xi(7)**2 * xi(16) &
         + (+8.660254037844385e-01_rk) * xi(6)**2 * xi(17) &
         + (-5.000000000000001e-01_rk) * xi(6)**2 * xi(16) &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * xi(16) &
         + (+8.660254037844385e-01_rk) * xi(4)**2 * xi(15) &
         + (+5.000000000000001e-01_rk) * xi(4)**2 * xi(14) &
         + (-8.660254037844386e-01_rk) * xi(3)**2 * xi(15) &
         + (+4.999999999999999e-01_rk) * xi(3)**2 * xi(14) &
         + (-1.000000000000000e+00_rk) * xi(2)**2 * xi(14)
    phi_z(55) = (+7.499999999999998e-01_rk) * xi(7) * xi(17)**2 &
         + (+8.660254037844385e-01_rk) * xi(7) * xi(16) * xi(17) &
         + (+2.500000000000001e-01_rk) * xi(7) * xi(16)**2 &
         + (+7.499999999999998e-01_rk) * xi(6) * xi(17)**2 &
         + (-8.660254037844386e-01_rk) * xi(6) * xi(16) * xi(17) &
         + (+2.500000000000000e-01_rk) * xi(6) * xi(16)**2 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(16)**2 &
         + (-7.499999999999998e-01_rk) * xi(4) * xi(15)**2 &
         + (-8.660254037844386e-01_rk) * xi(4) * xi(14) * xi(15) &
         + (-2.499999999999999e-01_rk) * xi(4) * xi(14)**2 &
         + (-7.499999999999998e-01_rk) * xi(3) * xi(15)**2 &
         + (+8.660254037844386e-01_rk) * xi(3) * xi(14) * xi(15) &
         + (-2.500000000000001e-01_rk) * xi(3) * xi(14)**2 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(14)**2
    phi_z(56) = (+2.500000000000001e-01_rk) * xi(7) * xi(17)**2 &
         + (-8.660254037844386e-01_rk) * xi(7) * xi(16) * xi(17) &
         + (+7.499999999999999e-01_rk) * xi(7) * xi(16)**2 &
         + (+2.500000000000001e-01_rk) * xi(6) * xi(17)**2 &
         + (+8.660254037844387e-01_rk) * xi(6) * xi(16) * xi(17) &
         + (+7.499999999999999e-01_rk) * xi(6) * xi(16)**2 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(17)**2 &
         + (-2.500000000000000e-01_rk) * xi(4) * xi(15)**2 &
         + (+8.660254037844387e-01_rk) * xi(4) * xi(14) * xi(15) &
         + (-7.499999999999999e-01_rk) * xi(4) * xi(14)**2 &
         + (-2.500000000000001e-01_rk) * xi(3) * xi(15)**2 &
         + (-8.660254037844387e-01_rk) * xi(3) * xi(14) * xi(15) &
         + (-7.499999999999999e-01_rk) * xi(3) * xi(14)**2 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(15)**2
    phi_z(57) = (-8.660254037844385e-01_rk) * xi(13)**2 * xi(17) &
         + (-4.999999999999999e-01_rk) * xi(13)**2 * xi(16) &
         + (+8.660254037844385e-01_rk) * xi(12)**2 * xi(17) &
         + (-5.000000000000001e-01_rk) * xi(12)**2 * xi(16) &
         + (+1.000000000000000e+00_rk) * xi(11)**2 * xi(16) &
         + (+8.660254037844385e-01_rk) * xi(10)**2 * xi(15) &
         + (+5.000000000000001e-01_rk) * xi(10)**2 * xi(14) &
         + (-8.660254037844386e-01_rk) * xi(9)**2 * xi(15) &
         + (+4.999999999999999e-01_rk) * xi(9)**2 * xi(14) &
         + (-1.000000000000000e+00_rk) * xi(8)**2 * xi(14)
    phi_z(58) = (+7.499999999999998e-01_rk) * xi(13) * xi(17)**2 &
         + (+8.660254037844385e-01_rk) * xi(13) * xi(16) * xi(17) &
         + (+2.500000000000001e-01_rk) * xi(13) * xi(16)**2 &
         + (+7.499999999999998e-01_rk) * xi(12) * xi(17)**2 &
         + (-8.660254037844386e-01_rk) * xi(12) * xi(16) * xi(17) &
         + (+2.500000000000000e-01_rk) * xi(12) * xi(16)**2 &
         + (+1.000000000000000e+00_rk) * xi(11) * xi(16)**2 &
         + (-7.499999999999998e-01_rk) * xi(10) * xi(15)**2 &
         + (-8.660254037844386e-01_rk) * xi(10) * xi(14) * xi(15) &
         + (-2.499999999999999e-01_rk) * xi(10) * xi(14)**2 &
         + (-7.499999999999998e-01_rk) * xi(9) * xi(15)**2 &
         + (+8.660254037844386e-01_rk) * xi(9) * xi(14) * xi(15) &
         + (-2.500000000000001e-01_rk) * xi(9) * xi(14)**2 &
         + (-1.000000000000000e+00_rk) * xi(8) * xi(14)**2
    phi_z(59) = (+2.500000000000001e-01_rk) * xi(13) * xi(17)**2 &
         + (-8.660254037844386e-01_rk) * xi(13) * xi(16) * xi(17) &
         + (+7.499999999999999e-01_rk) * xi(13) * xi(16)**2 &
         + (+2.500000000000001e-01_rk) * xi(12) * xi(17)**2 &
         + (+8.660254037844387e-01_rk) * xi(12) * xi(16) * xi(17) &
         + (+7.499999999999999e-01_rk) * xi(12) * xi(16)**2 &
         + (+1.000000000000000e+00_rk) * xi(11) * xi(17)**2 &
         + (-2.500000000000000e-01_rk) * xi(10) * xi(15)**2 &
         + (+8.660254037844387e-01_rk) * xi(10) * xi(14) * xi(15) &
         + (-7.499999999999999e-01_rk) * xi(10) * xi(14)**2 &
         + (-2.500000000000001e-01_rk) * xi(9) * xi(15)**2 &
         + (-8.660254037844387e-01_rk) * xi(9) * xi(14) * xi(15) &
         + (-7.499999999999999e-01_rk) * xi(9) * xi(14)**2 &
         + (-1.000000000000000e+00_rk) * xi(8) * xi(15)**2
    phi_z(60) = (+1.000000000000000e+00_rk) * xi(4) * xi(7)**2 * c_2 &
         + (-5.000000000000003e-01_rk) * xi(4) * xi(6)**2 * c_2 &
         + (+8.660254037844385e-01_rk) * xi(4) * xi(6)**2 * s_2 &
         + (-4.999999999999997e-01_rk) * xi(4) * xi(5)**2 * c_2 &
         + (-8.660254037844387e-01_rk) * xi(4) * xi(5)**2 * s_2 &
         + (-1.000000000000000e+00_rk) * xi(4)**2 * xi(7) * c_2 &
         + (+5.000000000000003e-01_rk) * xi(4)**2 * xi(6) * c_2 &
         + (-8.660254037844385e-01_rk) * xi(4)**2 * xi(6) * s_2 &
         + (+4.999999999999997e-01_rk) * xi(4)**2 * xi(5) * c_2 &
         + (+8.660254037844387e-01_rk) * xi(4)**2 * xi(5) * s_2 &
         + (-4.999999999999997e-01_rk) * xi(3) * xi(7)**2 * c_2 &
         + (-8.660254037844387e-01_rk) * xi(3) * xi(7)**2 * s_2 &
         + (+1.000000000000000e+00_rk) * xi(3) * xi(6)**2 * c_2 &
         + (-5.000000000000003e-01_rk) * xi(3) * xi(5)**2 * c_2 &
         + (+8.660254037844385e-01_rk) * xi(3) * xi(5)**2 * s_2 &
         + (+4.999999999999994e-01_rk) * xi(3)**2 * xi(7) * c_2 &
         + (+8.660254037844388e-01_rk) * xi(3)**2 * xi(7) * s_2 &
         + (-1.000000000000000e+00_rk) * xi(3)**2 * xi(6) * c_2 &
         + (+5.000000000000003e-01_rk) * xi(3)**2 * xi(5) * c_2 &
         + (-8.660254037844385e-01_rk) * xi(3)**2 * xi(5) * s_2 &
         + (-5.000000000000003e-01_rk) * xi(2) * xi(7)**2 * c_2 &
         + (+8.660254037844385e-01_rk) * xi(2) * xi(7)**2 * s_2 &
         + (-4.999999999999997e-01_rk) * xi(2) * xi(6)**2 * c_2 &
         + (-8.660254037844387e-01_rk) * xi(2) * xi(6)**2 * s_2 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(5)**2 * c_2 &
         + (+5.000000000000003e-01_rk) * xi(2)**2 * xi(7) * c_2 &
         + (-8.660254037844385e-01_rk) * xi(2)**2 * xi(7) * s_2 &
         + (+4.999999999999997e-01_rk) * xi(2)**2 * xi(6) * c_2 &
         + (+8.660254037844387e-01_rk) * xi(2)**2 * xi(6) * s_2 &
         + (-1.000000000000000e+00_rk) * xi(2)**2 * xi(5) * c_2
    phi_z(61) = (+1.000000000000000e+00_rk) * xi(10) * xi(13)**2 * c_2 &
         + (-5.000000000000003e-01_rk) * xi(10) * xi(12)**2 * c_2 &
         + (+8.660254037844385e-01_rk) * xi(10) * xi(12)**2 * s_2 &
         + (-4.999999999999997e-01_rk) * xi(10) * xi(11)**2 * c_2 &
         + (-8.660254037844387e-01_rk) * xi(10) * xi(11)**2 * s_2 &
         + (-1.000000000000000e+00_rk) * xi(10)**2 * xi(13) * c_2 &
         + (+5.000000000000003e-01_rk) * xi(10)**2 * xi(12) * c_2 &
         + (-8.660254037844385e-01_rk) * xi(10)**2 * xi(12) * s_2 &
         + (+4.999999999999997e-01_rk) * xi(10)**2 * xi(11) * c_2 &
         + (+8.660254037844387e-01_rk) * xi(10)**2 * xi(11) * s_2 &
         + (-4.999999999999997e-01_rk) * xi(9) * xi(13)**2 * c_2 &
         + (-8.660254037844387e-01_rk) * xi(9) * xi(13)**2 * s_2 &
         + (+1.000000000000000e+00_rk) * xi(9) * xi(12)**2 * c_2 &
         + (-5.000000000000003e-01_rk) * xi(9) * xi(11)**2 * c_2 &
         + (+8.660254037844385e-01_rk) * xi(9) * xi(11)**2 * s_2 &
         + (+4.999999999999994e-01_rk) * xi(9)**2 * xi(13) * c_2 &
         + (+8.660254037844388e-01_rk) * xi(9)**2 * xi(13) * s_2 &
         + (-1.000000000000000e+00_rk) * xi(9)**2 * xi(12) * c_2 &
         + (+5.000000000000003e-01_rk) * xi(9)**2 * xi(11) * c_2 &
         + (-8.660254037844385e-01_rk) * xi(9)**2 * xi(11) * s_2 &
         + (-5.000000000000003e-01_rk) * xi(8) * xi(13)**2 * c_2 &
         + (+8.660254037844385e-01_rk) * xi(8) * xi(13)**2 * s_2 &
         + (-4.999999999999997e-01_rk) * xi(8) * xi(12)**2 * c_2 &
         + (-8.660254037844387e-01_rk) * xi(8) * xi(12)**2 * s_2 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(11)**2 * c_2 &
         + (+5.000000000000003e-01_rk) * xi(8)**2 * xi(13) * c_2 &
         + (-8.660254037844385e-01_rk) * xi(8)**2 * xi(13) * s_2 &
         + (+4.999999999999997e-01_rk) * xi(8)**2 * xi(12) * c_2 &
         + (+8.660254037844387e-01_rk) * xi(8)**2 * xi(12) * s_2 &
         + (-1.000000000000000e+00_rk) * xi(8)**2 * xi(11) * c_2
    phi_z(62) = (+1.000000000000000e+00_rk) * xi(7)**4 &
         + (+1.000000000000000e+00_rk) * xi(6)**4 &
         + (+1.000000000000000e+00_rk) * xi(5)**4 &
         + (-1.000000000000000e+00_rk) * xi(4)**4 &
         + (-1.000000000000000e+00_rk) * xi(3)**4 &
         + (-1.000000000000000e+00_rk) * xi(2)**4
    phi_z(63) = (+1.000000000000000e+00_rk) * xi(13)**4 &
         + (+1.000000000000000e+00_rk) * xi(12)**4 &
         + (+1.000000000000000e+00_rk) * xi(11)**4 &
         + (-1.000000000000000e+00_rk) * xi(10)**4 &
         + (-1.000000000000000e+00_rk) * xi(9)**4 &
         + (-1.000000000000000e+00_rk) * xi(8)**4
    phi_z(64) = (+4.999999999999999e-01_rk) * xi(17)**4 &
         + (+9.999999999999999e-01_rk) * xi(16)**2 * xi(17)**2 &
         + (+4.999999999999998e-01_rk) * xi(16)**4 &
         + (-4.999999999999999e-01_rk) * xi(15)**4 &
         + (-1.000000000000000e+00_rk) * xi(14)**2 * xi(15)**2 &
         + (-4.999999999999997e-01_rk) * xi(14)**4
    phi_z(65) = (+1.000000000000000e+00_rk) * xi(1)**2 * xi(7)**2 &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(6)**2 &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(5)**2 &
         + (-1.000000000000000e+00_rk) * xi(1)**2 * xi(4)**2 &
         + (-1.000000000000000e+00_rk) * xi(1)**2 * xi(3)**2 &
         + (-1.000000000000000e+00_rk) * xi(1)**2 * xi(2)**2
    phi_z(66) = (+1.000000000000000e+00_rk) * xi(1)**2 * xi(13)**2 &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(12)**2 &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(11)**2 &
         + (-1.000000000000000e+00_rk) * xi(1)**2 * xi(10)**2 &
         + (-1.000000000000000e+00_rk) * xi(1)**2 * xi(9)**2 &
         + (-1.000000000000000e+00_rk) * xi(1)**2 * xi(8)**2
    phi_z(67) = (+1.000000000000000e+00_rk) * xi(6)**2 * xi(7)**2 &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * xi(7)**2 &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * xi(6)**2 &
         + (-1.000000000000000e+00_rk) * xi(3)**2 * xi(4)**2 &
         + (-1.000000000000000e+00_rk) * xi(2)**2 * xi(4)**2 &
         + (-1.000000000000000e+00_rk) * xi(2)**2 * xi(3)**2
    phi_z(68) = (+1.000000000000000e+00_rk) * xi(12)**2 * xi(13)**2 &
         + (+1.000000000000000e+00_rk) * xi(11)**2 * xi(13)**2 &
         + (+1.000000000000000e+00_rk) * xi(11)**2 * xi(12)**2 &
         + (-1.000000000000000e+00_rk) * xi(9)**2 * xi(10)**2 &
         + (-1.000000000000000e+00_rk) * xi(8)**2 * xi(10)**2 &
         + (-1.000000000000000e+00_rk) * xi(8)**2 * xi(9)**2
    phi_z(69) = (+1.000000000000000e+00_rk) * xi(7)**2 * xi(13)**2 &
         + (+1.000000000000000e+00_rk) * xi(6)**2 * xi(12)**2 &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * xi(11)**2 &
         + (-1.000000000000000e+00_rk) * xi(4)**2 * xi(10)**2 &
         + (-1.000000000000000e+00_rk) * xi(3)**2 * xi(9)**2 &
         + (-1.000000000000000e+00_rk) * xi(2)**2 * xi(8)**2
    phi_z(70) = (+1.000000000000000e+00_rk) * xi(7)**2 * xi(12)**2 &
         + (+1.000000000000000e+00_rk) * xi(7)**2 * xi(11)**2 &
         + (+1.000000000000000e+00_rk) * xi(6)**2 * xi(13)**2 &
         + (+1.000000000000000e+00_rk) * xi(6)**2 * xi(11)**2 &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * xi(13)**2 &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * xi(12)**2 &
         + (-1.000000000000000e+00_rk) * xi(4)**2 * xi(9)**2 &
         + (-1.000000000000000e+00_rk) * xi(4)**2 * xi(8)**2 &
         + (-1.000000000000000e+00_rk) * xi(3)**2 * xi(10)**2 &
         + (-1.000000000000000e+00_rk) * xi(3)**2 * xi(8)**2 &
         + (-1.000000000000000e+00_rk) * xi(2)**2 * xi(10)**2 &
         + (-1.000000000000000e+00_rk) * xi(2)**2 * xi(9)**2

    ! --- Perpendicular dipole basis functions phi_x, phi_y (E_2d) ---
    phi_x(1) = (-5.000000000000000e-01_rk) * xi(7) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(7) * s_1 &
         + (-4.999999999999997e-01_rk) * xi(6) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(6) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(5) * c_1 &
         + (-4.999999999999997e-01_rk) * xi(4) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(4) * s_1 &
         + (-5.000000000000002e-01_rk) * xi(3) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(3) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2) * c_1
    phi_y(1) = (+8.660254037844385e-01_rk) * xi(7) * c_1 &
         + (+5.000000000000000e-01_rk) * xi(7) * s_1 &
         + (-8.660254037844386e-01_rk) * xi(6) * c_1 &
         + (+4.999999999999999e-01_rk) * xi(6) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(5) * s_1 &
         + (+8.660254037844386e-01_rk) * xi(4) * c_1 &
         + (-4.999999999999999e-01_rk) * xi(4) * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3) * c_1 &
         + (-5.000000000000000e-01_rk) * xi(3) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2) * s_1
    phi_x(2) = (-5.000000000000000e-01_rk) * xi(13) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(13) * s_1 &
         + (-4.999999999999997e-01_rk) * xi(12) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(12) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(11) * c_1 &
         + (-4.999999999999997e-01_rk) * xi(10) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(10) * s_1 &
         + (-5.000000000000002e-01_rk) * xi(9) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(9) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8) * c_1
    phi_y(2) = (+8.660254037844385e-01_rk) * xi(13) * c_1 &
         + (+5.000000000000000e-01_rk) * xi(13) * s_1 &
         + (-8.660254037844386e-01_rk) * xi(12) * c_1 &
         + (+4.999999999999999e-01_rk) * xi(12) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(11) * s_1 &
         + (+8.660254037844386e-01_rk) * xi(10) * c_1 &
         + (-4.999999999999999e-01_rk) * xi(10) * s_1 &
         + (-8.660254037844385e-01_rk) * xi(9) * c_1 &
         + (-5.000000000000000e-01_rk) * xi(9) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8) * s_1
    phi_x(3) = (+9.999999999999996e-01_rk) * xi(17) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(16) * c_1 &
         + (-9.999999999999996e-01_rk) * xi(15) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(14) * c_1
    phi_y(3) = (+9.999999999999997e-01_rk) * xi(17) * c_1 &
         + (+9.999999999999996e-01_rk) * xi(16) * s_1 &
         + (+9.999999999999997e-01_rk) * xi(15) * c_1 &
         + (-9.999999999999996e-01_rk) * xi(14) * s_1
    phi_x(4) = (-4.999999999999994e-01_rk) * xi(7) * c_5 &
         + (-8.660254037844388e-01_rk) * xi(7) * s_5 &
         + (-5.000000000000006e-01_rk) * xi(6) * c_5 &
         + (+8.660254037844385e-01_rk) * xi(6) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(5) * c_5 &
         + (-5.000000000000006e-01_rk) * xi(4) * c_5 &
         + (+8.660254037844385e-01_rk) * xi(4) * s_5 &
         + (-4.999999999999989e-01_rk) * xi(3) * c_5 &
         + (-8.660254037844389e-01_rk) * xi(3) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(2) * c_5
    phi_y(4) = (+8.660254037844389e-01_rk) * xi(7) * c_5 &
         + (-4.999999999999997e-01_rk) * xi(7) * s_5 &
         + (-8.660254037844385e-01_rk) * xi(6) * c_5 &
         + (-5.000000000000001e-01_rk) * xi(6) * s_5 &
         + (+9.999999999999998e-01_rk) * xi(5) * s_5 &
         + (+8.660254037844385e-01_rk) * xi(4) * c_5 &
         + (+5.000000000000001e-01_rk) * xi(4) * s_5 &
         + (-8.660254037844389e-01_rk) * xi(3) * c_5 &
         + (+4.999999999999997e-01_rk) * xi(3) * s_5 &
         + (-9.999999999999998e-01_rk) * xi(2) * s_5
    phi_x(5) = (-4.999999999999994e-01_rk) * xi(13) * c_5 &
         + (-8.660254037844388e-01_rk) * xi(13) * s_5 &
         + (-5.000000000000006e-01_rk) * xi(12) * c_5 &
         + (+8.660254037844385e-01_rk) * xi(12) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(11) * c_5 &
         + (-5.000000000000006e-01_rk) * xi(10) * c_5 &
         + (+8.660254037844385e-01_rk) * xi(10) * s_5 &
         + (-4.999999999999989e-01_rk) * xi(9) * c_5 &
         + (-8.660254037844389e-01_rk) * xi(9) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(8) * c_5
    phi_y(5) = (+8.660254037844389e-01_rk) * xi(13) * c_5 &
         + (-4.999999999999997e-01_rk) * xi(13) * s_5 &
         + (-8.660254037844385e-01_rk) * xi(12) * c_5 &
         + (-5.000000000000001e-01_rk) * xi(12) * s_5 &
         + (+9.999999999999998e-01_rk) * xi(11) * s_5 &
         + (+8.660254037844385e-01_rk) * xi(10) * c_5 &
         + (+5.000000000000001e-01_rk) * xi(10) * s_5 &
         + (-8.660254037844389e-01_rk) * xi(9) * c_5 &
         + (+4.999999999999997e-01_rk) * xi(9) * s_5 &
         + (-9.999999999999998e-01_rk) * xi(8) * s_5
    phi_x(6) = (-9.999999999999998e-01_rk) * xi(17) * s_5 &
         + (-1.000000000000000e+00_rk) * xi(16) * c_5 &
         + (+1.000000000000000e+00_rk) * xi(15) * s_5 &
         + (-1.000000000000000e+00_rk) * xi(14) * c_5
    phi_y(6) = (+9.999999999999997e-01_rk) * xi(17) * c_5 &
         + (-9.999999999999996e-01_rk) * xi(16) * s_5 &
         + (+9.999999999999997e-01_rk) * xi(15) * c_5 &
         + (+9.999999999999996e-01_rk) * xi(14) * s_5
    phi_x(7) = (-5.000000000000000e-01_rk) * xi(7)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(7)**2 * s_1 &
         + (-4.999999999999997e-01_rk) * xi(6)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(6)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * c_1 &
         + (-4.999999999999997e-01_rk) * xi(4)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(4)**2 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(3)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(3)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**2 * c_1
    phi_y(7) = (+8.660254037844385e-01_rk) * xi(7)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(7)**2 * s_1 &
         + (-8.660254037844386e-01_rk) * xi(6)**2 * c_1 &
         + (+4.999999999999999e-01_rk) * xi(6)**2 * s_1 &
         + (-1.000000000000000e+00_rk) * xi(5)**2 * s_1 &
         + (+8.660254037844386e-01_rk) * xi(4)**2 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(4)**2 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(3)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**2 * s_1
    phi_x(8) = (-5.000000000000000e-01_rk) * xi(13)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(13)**2 * s_1 &
         + (-4.999999999999997e-01_rk) * xi(12)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(12)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(11)**2 * c_1 &
         + (-4.999999999999997e-01_rk) * xi(10)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(10)**2 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(9)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(9)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8)**2 * c_1
    phi_y(8) = (+8.660254037844385e-01_rk) * xi(13)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(13)**2 * s_1 &
         + (-8.660254037844386e-01_rk) * xi(12)**2 * c_1 &
         + (+4.999999999999999e-01_rk) * xi(12)**2 * s_1 &
         + (-1.000000000000000e+00_rk) * xi(11)**2 * s_1 &
         + (+8.660254037844386e-01_rk) * xi(10)**2 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(10)**2 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(9)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(9)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8)**2 * s_1
    phi_x(9) = (-5.000000000000000e-01_rk) * xi(17)**2 * c_1 &
         + (+9.999999999999998e-01_rk) * xi(16) * xi(17) * s_1 &
         + (+4.999999999999999e-01_rk) * xi(16)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(15)**2 * c_1 &
         + (-9.999999999999998e-01_rk) * xi(14) * xi(15) * s_1 &
         + (+4.999999999999999e-01_rk) * xi(14)**2 * c_1
    phi_y(9) = (+4.999999999999998e-01_rk) * xi(17)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(16) * xi(17) * c_1 &
         + (-4.999999999999998e-01_rk) * xi(16)**2 * s_1 &
         + (-4.999999999999998e-01_rk) * xi(15)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(14) * xi(15) * c_1 &
         + (+4.999999999999998e-01_rk) * xi(14)**2 * s_1
    phi_x(10) = (-5.000000000000000e-01_rk) * xi(1) * xi(7) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(1) * xi(7) * s_1 &
         + (-4.999999999999997e-01_rk) * xi(1) * xi(6) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(1) * xi(6) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(5) * c_1 &
         + (-4.999999999999997e-01_rk) * xi(1) * xi(4) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(1) * xi(4) * s_1 &
         + (-5.000000000000002e-01_rk) * xi(1) * xi(3) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(1) * xi(3) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(2) * c_1
    phi_y(10) = (+8.660254037844385e-01_rk) * xi(1) * xi(7) * c_1 &
         + (+5.000000000000000e-01_rk) * xi(1) * xi(7) * s_1 &
         + (-8.660254037844386e-01_rk) * xi(1) * xi(6) * c_1 &
         + (+4.999999999999999e-01_rk) * xi(1) * xi(6) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(5) * s_1 &
         + (+8.660254037844386e-01_rk) * xi(1) * xi(4) * c_1 &
         + (-4.999999999999999e-01_rk) * xi(1) * xi(4) * s_1 &
         + (-8.660254037844385e-01_rk) * xi(1) * xi(3) * c_1 &
         + (-5.000000000000000e-01_rk) * xi(1) * xi(3) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(2) * s_1
    phi_x(11) = (-5.000000000000000e-01_rk) * xi(1) * xi(13) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(1) * xi(13) * s_1 &
         + (-4.999999999999997e-01_rk) * xi(1) * xi(12) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(1) * xi(12) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(11) * c_1 &
         + (-4.999999999999997e-01_rk) * xi(1) * xi(10) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(1) * xi(10) * s_1 &
         + (-5.000000000000002e-01_rk) * xi(1) * xi(9) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(1) * xi(9) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(8) * c_1
    phi_y(11) = (+8.660254037844385e-01_rk) * xi(1) * xi(13) * c_1 &
         + (+5.000000000000000e-01_rk) * xi(1) * xi(13) * s_1 &
         + (-8.660254037844386e-01_rk) * xi(1) * xi(12) * c_1 &
         + (+4.999999999999999e-01_rk) * xi(1) * xi(12) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(11) * s_1 &
         + (+8.660254037844386e-01_rk) * xi(1) * xi(10) * c_1 &
         + (-4.999999999999999e-01_rk) * xi(1) * xi(10) * s_1 &
         + (-8.660254037844385e-01_rk) * xi(1) * xi(9) * c_1 &
         + (-5.000000000000000e-01_rk) * xi(1) * xi(9) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(8) * s_1
    phi_x(12) = (+9.999999999999996e-01_rk) * xi(1) * xi(17) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(16) * c_1 &
         + (-9.999999999999996e-01_rk) * xi(1) * xi(15) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(14) * c_1
    phi_y(12) = (+9.999999999999997e-01_rk) * xi(1) * xi(17) * c_1 &
         + (+9.999999999999996e-01_rk) * xi(1) * xi(16) * s_1 &
         + (+9.999999999999997e-01_rk) * xi(1) * xi(15) * c_1 &
         + (-9.999999999999996e-01_rk) * xi(1) * xi(14) * s_1
    phi_x(13) = (+9.999999999999991e-01_rk) * xi(6) * xi(7) * c_1 &
         + (-4.999999999999995e-01_rk) * xi(5) * xi(7) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(5) * xi(7) * s_1 &
         + (-4.999999999999999e-01_rk) * xi(5) * xi(6) * c_1 &
         + (+8.660254037844386e-01_rk) * xi(5) * xi(6) * s_1 &
         + (+9.999999999999994e-01_rk) * xi(3) * xi(4) * c_1 &
         + (-4.999999999999999e-01_rk) * xi(2) * xi(4) * c_1 &
         + (+8.660254037844386e-01_rk) * xi(2) * xi(4) * s_1 &
         + (-4.999999999999994e-01_rk) * xi(2) * xi(3) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(2) * xi(3) * s_1
    phi_y(13) = (-1.000000000000000e+00_rk) * xi(6) * xi(7) * s_1 &
         + (-8.660254037844383e-01_rk) * xi(5) * xi(7) * c_1 &
         + (+5.000000000000001e-01_rk) * xi(5) * xi(7) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(5) * xi(6) * c_1 &
         + (+4.999999999999999e-01_rk) * xi(5) * xi(6) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(3) * xi(4) * s_1 &
         + (-8.660254037844382e-01_rk) * xi(2) * xi(4) * c_1 &
         + (-4.999999999999999e-01_rk) * xi(2) * xi(4) * s_1 &
         + (+8.660254037844383e-01_rk) * xi(2) * xi(3) * c_1 &
         + (-5.000000000000001e-01_rk) * xi(2) * xi(3) * s_1
    phi_x(14) = (-4.999999999999998e-01_rk) * xi(4) * xi(7) * c_1 &
         + (-4.999999999999998e-01_rk) * xi(4) * xi(6) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(4) * xi(6) * s_1 &
         + (+2.500000000000000e-01_rk) * xi(4) * xi(5) * c_1 &
         + (-4.330127018922191e-01_rk) * xi(4) * xi(5) * s_1 &
         + (-5.000000000000000e-01_rk) * xi(3) * xi(7) * c_1 &
         + (+8.660254037844383e-01_rk) * xi(3) * xi(7) * s_1 &
         + (-4.999999999999998e-01_rk) * xi(3) * xi(6) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(3) * xi(5) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(3) * xi(5) * s_1 &
         + (+2.499999999999999e-01_rk) * xi(2) * xi(7) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(2) * xi(7) * s_1 &
         + (+2.500000000000000e-01_rk) * xi(2) * xi(6) * c_1 &
         + (-4.330127018922191e-01_rk) * xi(2) * xi(6) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(5) * c_1
    phi_y(14) = (+8.660254037844385e-01_rk) * xi(4) * xi(7) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(4) * xi(5) * c_1 &
         + (-7.499999999999997e-01_rk) * xi(4) * xi(5) * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3) * xi(6) * c_1 &
         + (-4.330127018922191e-01_rk) * xi(3) * xi(5) * c_1 &
         + (-7.499999999999998e-01_rk) * xi(3) * xi(5) * s_1 &
         + (+4.330127018922191e-01_rk) * xi(2) * xi(7) * c_1 &
         + (+7.499999999999998e-01_rk) * xi(2) * xi(7) * s_1 &
         + (-4.330127018922193e-01_rk) * xi(2) * xi(6) * c_1 &
         + (+7.499999999999997e-01_rk) * xi(2) * xi(6) * s_1
    phi_x(15) = (+9.999999999999991e-01_rk) * xi(12) * xi(13) * c_1 &
         + (-4.999999999999995e-01_rk) * xi(11) * xi(13) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(11) * xi(13) * s_1 &
         + (-4.999999999999999e-01_rk) * xi(11) * xi(12) * c_1 &
         + (+8.660254037844386e-01_rk) * xi(11) * xi(12) * s_1 &
         + (+9.999999999999994e-01_rk) * xi(9) * xi(10) * c_1 &
         + (-4.999999999999999e-01_rk) * xi(8) * xi(10) * c_1 &
         + (+8.660254037844386e-01_rk) * xi(8) * xi(10) * s_1 &
         + (-4.999999999999994e-01_rk) * xi(8) * xi(9) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(8) * xi(9) * s_1
    phi_y(15) = (-1.000000000000000e+00_rk) * xi(12) * xi(13) * s_1 &
         + (-8.660254037844383e-01_rk) * xi(11) * xi(13) * c_1 &
         + (+5.000000000000001e-01_rk) * xi(11) * xi(13) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(11) * xi(12) * c_1 &
         + (+4.999999999999999e-01_rk) * xi(11) * xi(12) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(9) * xi(10) * s_1 &
         + (-8.660254037844382e-01_rk) * xi(8) * xi(10) * c_1 &
         + (-4.999999999999999e-01_rk) * xi(8) * xi(10) * s_1 &
         + (+8.660254037844383e-01_rk) * xi(8) * xi(9) * c_1 &
         + (-5.000000000000001e-01_rk) * xi(8) * xi(9) * s_1
    phi_x(16) = (-4.999999999999998e-01_rk) * xi(10) * xi(13) * c_1 &
         + (-4.999999999999998e-01_rk) * xi(10) * xi(12) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(10) * xi(12) * s_1 &
         + (+2.500000000000000e-01_rk) * xi(10) * xi(11) * c_1 &
         + (-4.330127018922191e-01_rk) * xi(10) * xi(11) * s_1 &
         + (-5.000000000000000e-01_rk) * xi(9) * xi(13) * c_1 &
         + (+8.660254037844383e-01_rk) * xi(9) * xi(13) * s_1 &
         + (-4.999999999999998e-01_rk) * xi(9) * xi(12) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(9) * xi(11) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(9) * xi(11) * s_1 &
         + (+2.499999999999999e-01_rk) * xi(8) * xi(13) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(8) * xi(13) * s_1 &
         + (+2.500000000000000e-01_rk) * xi(8) * xi(12) * c_1 &
         + (-4.330127018922191e-01_rk) * xi(8) * xi(12) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(11) * c_1
    phi_y(16) = (+8.660254037844385e-01_rk) * xi(10) * xi(13) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(10) * xi(11) * c_1 &
         + (-7.499999999999997e-01_rk) * xi(10) * xi(11) * s_1 &
         + (-8.660254037844385e-01_rk) * xi(9) * xi(12) * c_1 &
         + (-4.330127018922191e-01_rk) * xi(9) * xi(11) * c_1 &
         + (-7.499999999999998e-01_rk) * xi(9) * xi(11) * s_1 &
         + (+4.330127018922191e-01_rk) * xi(8) * xi(13) * c_1 &
         + (+7.499999999999998e-01_rk) * xi(8) * xi(13) * s_1 &
         + (-4.330127018922193e-01_rk) * xi(8) * xi(12) * c_1 &
         + (+7.499999999999997e-01_rk) * xi(8) * xi(12) * s_1
    phi_x(17) = (-5.000000000000000e-01_rk) * xi(7) * xi(13) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(7) * xi(13) * s_1 &
         + (-4.999999999999997e-01_rk) * xi(6) * xi(12) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(6) * xi(12) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(11) * c_1 &
         + (-4.999999999999997e-01_rk) * xi(4) * xi(10) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(4) * xi(10) * s_1 &
         + (-5.000000000000002e-01_rk) * xi(3) * xi(9) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(3) * xi(9) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(8) * c_1
    phi_y(17) = (+8.660254037844385e-01_rk) * xi(7) * xi(13) * c_1 &
         + (+5.000000000000000e-01_rk) * xi(7) * xi(13) * s_1 &
         + (-8.660254037844386e-01_rk) * xi(6) * xi(12) * c_1 &
         + (+4.999999999999999e-01_rk) * xi(6) * xi(12) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(5) * xi(11) * s_1 &
         + (+8.660254037844386e-01_rk) * xi(4) * xi(10) * c_1 &
         + (-4.999999999999999e-01_rk) * xi(4) * xi(10) * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3) * xi(9) * c_1 &
         + (-5.000000000000000e-01_rk) * xi(3) * xi(9) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(8) * s_1
    phi_x(18) = (-4.999999999999995e-01_rk) * xi(7) * xi(12) * c_1 &
         + (+8.660254037844387e-01_rk) * xi(7) * xi(12) * s_1 &
         + (-5.000000000000006e-01_rk) * xi(7) * xi(11) * c_1 &
         + (+8.660254037844382e-01_rk) * xi(7) * xi(11) * s_1 &
         + (-4.999999999999993e-01_rk) * xi(6) * xi(13) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(6) * xi(13) * s_1 &
         + (-5.000000000000002e-01_rk) * xi(6) * xi(11) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(6) * xi(11) * s_1 &
         + (+9.999999999999999e-01_rk) * xi(5) * xi(13) * c_1 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(12) * c_1 &
         + (-4.999999999999993e-01_rk) * xi(4) * xi(9) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(4) * xi(9) * s_1 &
         + (-5.000000000000003e-01_rk) * xi(4) * xi(8) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(4) * xi(8) * s_1 &
         + (-4.999999999999997e-01_rk) * xi(3) * xi(10) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(3) * xi(10) * s_1 &
         + (-5.000000000000004e-01_rk) * xi(3) * xi(8) * c_1 &
         + (+8.660254037844382e-01_rk) * xi(3) * xi(8) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(10) * c_1 &
         + (+9.999999999999999e-01_rk) * xi(2) * xi(9) * c_1
    phi_y(18) = (+8.660254037844385e-01_rk) * xi(7) * xi(12) * c_1 &
         + (+4.999999999999998e-01_rk) * xi(7) * xi(12) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(7) * xi(11) * c_1 &
         + (+5.000000000000000e-01_rk) * xi(7) * xi(11) * s_1 &
         + (-8.660254037844388e-01_rk) * xi(6) * xi(13) * c_1 &
         + (+4.999999999999997e-01_rk) * xi(6) * xi(13) * s_1 &
         + (-8.660254037844382e-01_rk) * xi(6) * xi(11) * c_1 &
         + (+5.000000000000000e-01_rk) * xi(6) * xi(11) * s_1 &
         + (-9.999999999999999e-01_rk) * xi(5) * xi(13) * s_1 &
         + (-9.999999999999997e-01_rk) * xi(5) * xi(12) * s_1 &
         + (+8.660254037844388e-01_rk) * xi(4) * xi(9) * c_1 &
         + (-4.999999999999997e-01_rk) * xi(4) * xi(9) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(4) * xi(8) * c_1 &
         + (-5.000000000000000e-01_rk) * xi(4) * xi(8) * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3) * xi(10) * c_1 &
         + (-4.999999999999998e-01_rk) * xi(3) * xi(10) * s_1 &
         + (-8.660254037844382e-01_rk) * xi(3) * xi(8) * c_1 &
         + (-5.000000000000000e-01_rk) * xi(3) * xi(8) * s_1 &
         + (+9.999999999999997e-01_rk) * xi(2) * xi(10) * s_1 &
         + (+9.999999999999999e-01_rk) * xi(2) * xi(9) * s_1
    phi_x(19) = (-8.660254037844385e-01_rk) * xi(7) * xi(12) * c_1 &
         + (-4.999999999999998e-01_rk) * xi(7) * xi(12) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(7) * xi(11) * c_1 &
         + (+5.000000000000000e-01_rk) * xi(7) * xi(11) * s_1 &
         + (-8.660254037844388e-01_rk) * xi(6) * xi(13) * c_1 &
         + (+4.999999999999997e-01_rk) * xi(6) * xi(13) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(6) * xi(11) * c_1 &
         + (-5.000000000000000e-01_rk) * xi(6) * xi(11) * s_1 &
         + (+9.999999999999999e-01_rk) * xi(5) * xi(13) * s_1 &
         + (-9.999999999999997e-01_rk) * xi(5) * xi(12) * s_1 &
         + (-8.660254037844388e-01_rk) * xi(4) * xi(9) * c_1 &
         + (+4.999999999999997e-01_rk) * xi(4) * xi(9) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(4) * xi(8) * c_1 &
         + (-5.000000000000000e-01_rk) * xi(4) * xi(8) * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3) * xi(10) * c_1 &
         + (-4.999999999999998e-01_rk) * xi(3) * xi(10) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(3) * xi(8) * c_1 &
         + (+5.000000000000000e-01_rk) * xi(3) * xi(8) * s_1 &
         + (-9.999999999999997e-01_rk) * xi(2) * xi(10) * s_1 &
         + (+9.999999999999999e-01_rk) * xi(2) * xi(9) * s_1
    phi_y(19) = (-4.999999999999995e-01_rk) * xi(7) * xi(12) * c_1 &
         + (+8.660254037844387e-01_rk) * xi(7) * xi(12) * s_1 &
         + (+5.000000000000006e-01_rk) * xi(7) * xi(11) * c_1 &
         + (-8.660254037844382e-01_rk) * xi(7) * xi(11) * s_1 &
         + (+4.999999999999993e-01_rk) * xi(6) * xi(13) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(6) * xi(13) * s_1 &
         + (-5.000000000000002e-01_rk) * xi(6) * xi(11) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(6) * xi(11) * s_1 &
         + (+9.999999999999999e-01_rk) * xi(5) * xi(13) * c_1 &
         + (-1.000000000000000e+00_rk) * xi(5) * xi(12) * c_1 &
         + (-4.999999999999993e-01_rk) * xi(4) * xi(9) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(4) * xi(9) * s_1 &
         + (+5.000000000000003e-01_rk) * xi(4) * xi(8) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(4) * xi(8) * s_1 &
         + (+4.999999999999997e-01_rk) * xi(3) * xi(10) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(3) * xi(10) * s_1 &
         + (-5.000000000000004e-01_rk) * xi(3) * xi(8) * c_1 &
         + (+8.660254037844382e-01_rk) * xi(3) * xi(8) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(10) * c_1 &
         + (-9.999999999999999e-01_rk) * xi(2) * xi(9) * c_1
    phi_x(20) = (-5.000000000000000e-01_rk) * xi(7) * xi(10) * c_1 &
         + (-5.000000000000001e-01_rk) * xi(7) * xi(9) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(7) * xi(9) * s_1 &
         + (+2.499999999999999e-01_rk) * xi(7) * xi(8) * c_1 &
         + (+4.330127018922194e-01_rk) * xi(7) * xi(8) * s_1 &
         + (-4.999999999999999e-01_rk) * xi(6) * xi(10) * c_1 &
         + (-8.660254037844387e-01_rk) * xi(6) * xi(10) * s_1 &
         + (-5.000000000000000e-01_rk) * xi(6) * xi(9) * c_1 &
         + (+2.500000000000001e-01_rk) * xi(6) * xi(8) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(6) * xi(8) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(5) * xi(10) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(5) * xi(10) * s_1 &
         + (+2.499999999999999e-01_rk) * xi(5) * xi(9) * c_1 &
         + (+4.330127018922194e-01_rk) * xi(5) * xi(9) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(8) * c_1 &
         + (-5.000000000000000e-01_rk) * xi(4) * xi(13) * c_1 &
         + (-4.999999999999999e-01_rk) * xi(4) * xi(12) * c_1 &
         + (-8.660254037844387e-01_rk) * xi(4) * xi(12) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(4) * xi(11) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(4) * xi(11) * s_1 &
         + (-5.000000000000002e-01_rk) * xi(3) * xi(13) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(3) * xi(13) * s_1 &
         + (-5.000000000000000e-01_rk) * xi(3) * xi(12) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(3) * xi(11) * c_1 &
         + (+4.330127018922194e-01_rk) * xi(3) * xi(11) * s_1 &
         + (+2.499999999999999e-01_rk) * xi(2) * xi(13) * c_1 &
         + (+4.330127018922194e-01_rk) * xi(2) * xi(13) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(2) * xi(12) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(2) * xi(12) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(11) * c_1
    phi_y(20) = (+8.660254037844386e-01_rk) * xi(7) * xi(10) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(7) * xi(8) * c_1 &
         + (+7.500000000000000e-01_rk) * xi(7) * xi(8) * s_1 &
         + (-8.660254037844386e-01_rk) * xi(6) * xi(9) * c_1 &
         + (-4.330127018922194e-01_rk) * xi(6) * xi(8) * c_1 &
         + (+7.499999999999999e-01_rk) * xi(6) * xi(8) * s_1 &
         + (+4.330127018922194e-01_rk) * xi(5) * xi(10) * c_1 &
         + (-7.499999999999999e-01_rk) * xi(5) * xi(10) * s_1 &
         + (-4.330127018922192e-01_rk) * xi(5) * xi(9) * c_1 &
         + (-7.500000000000000e-01_rk) * xi(5) * xi(9) * s_1 &
         + (+8.660254037844386e-01_rk) * xi(4) * xi(13) * c_1 &
         + (+4.330127018922194e-01_rk) * xi(4) * xi(11) * c_1 &
         + (-7.499999999999999e-01_rk) * xi(4) * xi(11) * s_1 &
         + (-8.660254037844386e-01_rk) * xi(3) * xi(12) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(3) * xi(11) * c_1 &
         + (-7.500000000000000e-01_rk) * xi(3) * xi(11) * s_1 &
         + (+4.330127018922192e-01_rk) * xi(2) * xi(13) * c_1 &
         + (+7.500000000000000e-01_rk) * xi(2) * xi(13) * s_1 &
         + (-4.330127018922194e-01_rk) * xi(2) * xi(12) * c_1 &
         + (+7.499999999999999e-01_rk) * xi(2) * xi(12) * s_1
    phi_x(21) = (-8.660254037844386e-01_rk) * xi(7) * xi(10) * s_1 &
         + (+7.500000000000000e-01_rk) * xi(7) * xi(8) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(7) * xi(8) * s_1 &
         + (+8.660254037844386e-01_rk) * xi(6) * xi(9) * s_1 &
         + (+7.499999999999999e-01_rk) * xi(6) * xi(8) * c_1 &
         + (+4.330127018922194e-01_rk) * xi(6) * xi(8) * s_1 &
         + (-7.499999999999999e-01_rk) * xi(5) * xi(10) * c_1 &
         + (-4.330127018922194e-01_rk) * xi(5) * xi(10) * s_1 &
         + (-7.500000000000000e-01_rk) * xi(5) * xi(9) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(5) * xi(9) * s_1 &
         + (+8.660254037844386e-01_rk) * xi(4) * xi(13) * s_1 &
         + (+7.499999999999999e-01_rk) * xi(4) * xi(11) * c_1 &
         + (+4.330127018922194e-01_rk) * xi(4) * xi(11) * s_1 &
         + (-8.660254037844386e-01_rk) * xi(3) * xi(12) * s_1 &
         + (+7.500000000000000e-01_rk) * xi(3) * xi(11) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(3) * xi(11) * s_1 &
         + (-7.500000000000000e-01_rk) * xi(2) * xi(13) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(2) * xi(13) * s_1 &
         + (-7.499999999999999e-01_rk) * xi(2) * xi(12) * c_1 &
         + (-4.330127018922194e-01_rk) * xi(2) * xi(12) * s_1
    phi_y(21) = (-5.000000000000000e-01_rk) * xi(7) * xi(10) * s_1 &
         + (-8.660254037844385e-01_rk) * xi(7) * xi(9) * c_1 &
         + (-5.000000000000001e-01_rk) * xi(7) * xi(9) * s_1 &
         + (-4.330127018922194e-01_rk) * xi(7) * xi(8) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(7) * xi(8) * s_1 &
         + (+8.660254037844387e-01_rk) * xi(6) * xi(10) * c_1 &
         + (-4.999999999999999e-01_rk) * xi(6) * xi(10) * s_1 &
         + (-5.000000000000000e-01_rk) * xi(6) * xi(9) * s_1 &
         + (+4.330127018922192e-01_rk) * xi(6) * xi(8) * c_1 &
         + (+2.500000000000001e-01_rk) * xi(6) * xi(8) * s_1 &
         + (+4.330127018922192e-01_rk) * xi(5) * xi(10) * c_1 &
         + (+2.500000000000001e-01_rk) * xi(5) * xi(10) * s_1 &
         + (-4.330127018922194e-01_rk) * xi(5) * xi(9) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(5) * xi(9) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(8) * s_1 &
         + (+5.000000000000000e-01_rk) * xi(4) * xi(13) * s_1 &
         + (-8.660254037844387e-01_rk) * xi(4) * xi(12) * c_1 &
         + (+4.999999999999999e-01_rk) * xi(4) * xi(12) * s_1 &
         + (-4.330127018922192e-01_rk) * xi(4) * xi(11) * c_1 &
         + (-2.500000000000001e-01_rk) * xi(4) * xi(11) * s_1 &
         + (+8.660254037844385e-01_rk) * xi(3) * xi(13) * c_1 &
         + (+5.000000000000002e-01_rk) * xi(3) * xi(13) * s_1 &
         + (+5.000000000000000e-01_rk) * xi(3) * xi(12) * s_1 &
         + (+4.330127018922194e-01_rk) * xi(3) * xi(11) * c_1 &
         + (-2.499999999999999e-01_rk) * xi(3) * xi(11) * s_1 &
         + (+4.330127018922194e-01_rk) * xi(2) * xi(13) * c_1 &
         + (-2.499999999999999e-01_rk) * xi(2) * xi(13) * s_1 &
         + (-4.330127018922192e-01_rk) * xi(2) * xi(12) * c_1 &
         + (-2.500000000000001e-01_rk) * xi(2) * xi(12) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(11) * s_1
    phi_x(22) = (+4.330127018922195e-01_rk) * xi(7) * xi(17) * c_1 &
         + (-7.499999999999998e-01_rk) * xi(7) * xi(17) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(7) * xi(16) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(7) * xi(16) * s_1 &
         + (-4.330127018922191e-01_rk) * xi(6) * xi(17) * c_1 &
         + (-7.500000000000000e-01_rk) * xi(6) * xi(17) * s_1 &
         + (+2.499999999999999e-01_rk) * xi(6) * xi(16) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(6) * xi(16) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(16) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(4) * xi(15) * c_1 &
         + (+7.500000000000000e-01_rk) * xi(4) * xi(15) * s_1 &
         + (+2.499999999999999e-01_rk) * xi(4) * xi(14) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(4) * xi(14) * s_1 &
         + (-4.330127018922196e-01_rk) * xi(3) * xi(15) * c_1 &
         + (+7.499999999999998e-01_rk) * xi(3) * xi(15) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(3) * xi(14) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(3) * xi(14) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(14) * c_1
    phi_y(22) = (-7.500000000000000e-01_rk) * xi(7) * xi(17) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(7) * xi(17) * s_1 &
         + (-4.330127018922192e-01_rk) * xi(7) * xi(16) * c_1 &
         + (-2.500000000000000e-01_rk) * xi(7) * xi(16) * s_1 &
         + (-7.500000000000001e-01_rk) * xi(6) * xi(17) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(6) * xi(17) * s_1 &
         + (+4.330127018922193e-01_rk) * xi(6) * xi(16) * c_1 &
         + (-2.499999999999999e-01_rk) * xi(6) * xi(16) * s_1 &
         + (-9.999999999999998e-01_rk) * xi(5) * xi(16) * s_1 &
         + (-7.500000000000001e-01_rk) * xi(4) * xi(15) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(4) * xi(15) * s_1 &
         + (-4.330127018922193e-01_rk) * xi(4) * xi(14) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(4) * xi(14) * s_1 &
         + (-7.500000000000000e-01_rk) * xi(3) * xi(15) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(3) * xi(15) * s_1 &
         + (+4.330127018922192e-01_rk) * xi(3) * xi(14) * c_1 &
         + (+2.500000000000000e-01_rk) * xi(3) * xi(14) * s_1 &
         + (+9.999999999999998e-01_rk) * xi(2) * xi(14) * s_1
    phi_x(23) = (+4.330127018922192e-01_rk) * xi(7) * xi(17) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(7) * xi(17) * s_1 &
         + (-7.500000000000000e-01_rk) * xi(7) * xi(16) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(7) * xi(16) * s_1 &
         + (-4.330127018922193e-01_rk) * xi(6) * xi(17) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(6) * xi(17) * s_1 &
         + (-7.500000000000001e-01_rk) * xi(6) * xi(16) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(6) * xi(16) * s_1 &
         + (+9.999999999999998e-01_rk) * xi(5) * xi(17) * s_1 &
         + (+4.330127018922193e-01_rk) * xi(4) * xi(15) * c_1 &
         + (-2.499999999999999e-01_rk) * xi(4) * xi(15) * s_1 &
         + (-7.500000000000001e-01_rk) * xi(4) * xi(14) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(4) * xi(14) * s_1 &
         + (-4.330127018922192e-01_rk) * xi(3) * xi(15) * c_1 &
         + (-2.499999999999999e-01_rk) * xi(3) * xi(15) * s_1 &
         + (-7.500000000000000e-01_rk) * xi(3) * xi(14) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(3) * xi(14) * s_1 &
         + (-9.999999999999998e-01_rk) * xi(2) * xi(15) * s_1
    phi_y(23) = (+2.500000000000001e-01_rk) * xi(7) * xi(17) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(7) * xi(17) * s_1 &
         + (-4.330127018922195e-01_rk) * xi(7) * xi(16) * c_1 &
         + (+7.499999999999998e-01_rk) * xi(7) * xi(16) * s_1 &
         + (+2.499999999999999e-01_rk) * xi(6) * xi(17) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(6) * xi(17) * s_1 &
         + (+4.330127018922191e-01_rk) * xi(6) * xi(16) * c_1 &
         + (+7.500000000000000e-01_rk) * xi(6) * xi(16) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(17) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(4) * xi(15) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(4) * xi(15) * s_1 &
         + (-4.330127018922193e-01_rk) * xi(4) * xi(14) * c_1 &
         + (-7.500000000000000e-01_rk) * xi(4) * xi(14) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(3) * xi(15) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(3) * xi(15) * s_1 &
         + (+4.330127018922196e-01_rk) * xi(3) * xi(14) * c_1 &
         + (-7.499999999999998e-01_rk) * xi(3) * xi(14) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(15) * c_1
    phi_x(24) = (-9.999999999999998e-01_rk) * xi(7) * xi(15) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(7) * xi(14) * c_1 &
         + (-1.000000000000000e+00_rk) * xi(6) * xi(15) * s_1 &
         + (-9.999999999999998e-01_rk) * xi(6) * xi(14) * c_1 &
         + (-9.999999999999998e-01_rk) * xi(5) * xi(15) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(5) * xi(14) * c_1 &
         + (+1.000000000000000e+00_rk) * xi(4) * xi(17) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(16) * c_1 &
         + (+9.999999999999997e-01_rk) * xi(3) * xi(17) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(16) * c_1 &
         + (+9.999999999999998e-01_rk) * xi(2) * xi(17) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(16) * c_1
    phi_y(24) = (+9.999999999999998e-01_rk) * xi(7) * xi(15) * c_1 &
         + (-9.999999999999998e-01_rk) * xi(7) * xi(14) * s_1 &
         + (+9.999999999999998e-01_rk) * xi(6) * xi(15) * c_1 &
         + (-9.999999999999997e-01_rk) * xi(6) * xi(14) * s_1 &
         + (+9.999999999999998e-01_rk) * xi(5) * xi(15) * c_1 &
         + (-9.999999999999998e-01_rk) * xi(5) * xi(14) * s_1 &
         + (+9.999999999999998e-01_rk) * xi(4) * xi(17) * c_1 &
         + (+9.999999999999997e-01_rk) * xi(4) * xi(16) * s_1 &
         + (+9.999999999999998e-01_rk) * xi(3) * xi(17) * c_1 &
         + (+9.999999999999998e-01_rk) * xi(3) * xi(16) * s_1 &
         + (+9.999999999999998e-01_rk) * xi(2) * xi(17) * c_1 &
         + (+9.999999999999998e-01_rk) * xi(2) * xi(16) * s_1
    phi_x(25) = (+4.330127018922195e-01_rk) * xi(13) * xi(17) * c_1 &
         + (-7.499999999999998e-01_rk) * xi(13) * xi(17) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(13) * xi(16) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(13) * xi(16) * s_1 &
         + (-4.330127018922191e-01_rk) * xi(12) * xi(17) * c_1 &
         + (-7.500000000000000e-01_rk) * xi(12) * xi(17) * s_1 &
         + (+2.499999999999999e-01_rk) * xi(12) * xi(16) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(12) * xi(16) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(11) * xi(16) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(10) * xi(15) * c_1 &
         + (+7.500000000000000e-01_rk) * xi(10) * xi(15) * s_1 &
         + (+2.499999999999999e-01_rk) * xi(10) * xi(14) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(10) * xi(14) * s_1 &
         + (-4.330127018922196e-01_rk) * xi(9) * xi(15) * c_1 &
         + (+7.499999999999998e-01_rk) * xi(9) * xi(15) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(9) * xi(14) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(9) * xi(14) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(14) * c_1
    phi_y(25) = (-7.500000000000000e-01_rk) * xi(13) * xi(17) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(13) * xi(17) * s_1 &
         + (-4.330127018922192e-01_rk) * xi(13) * xi(16) * c_1 &
         + (-2.500000000000000e-01_rk) * xi(13) * xi(16) * s_1 &
         + (-7.500000000000001e-01_rk) * xi(12) * xi(17) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(12) * xi(17) * s_1 &
         + (+4.330127018922193e-01_rk) * xi(12) * xi(16) * c_1 &
         + (-2.499999999999999e-01_rk) * xi(12) * xi(16) * s_1 &
         + (-9.999999999999998e-01_rk) * xi(11) * xi(16) * s_1 &
         + (-7.500000000000001e-01_rk) * xi(10) * xi(15) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(10) * xi(15) * s_1 &
         + (-4.330127018922193e-01_rk) * xi(10) * xi(14) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(10) * xi(14) * s_1 &
         + (-7.500000000000000e-01_rk) * xi(9) * xi(15) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(9) * xi(15) * s_1 &
         + (+4.330127018922192e-01_rk) * xi(9) * xi(14) * c_1 &
         + (+2.500000000000000e-01_rk) * xi(9) * xi(14) * s_1 &
         + (+9.999999999999998e-01_rk) * xi(8) * xi(14) * s_1
    phi_x(26) = (+4.330127018922192e-01_rk) * xi(13) * xi(17) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(13) * xi(17) * s_1 &
         + (-7.500000000000000e-01_rk) * xi(13) * xi(16) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(13) * xi(16) * s_1 &
         + (-4.330127018922193e-01_rk) * xi(12) * xi(17) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(12) * xi(17) * s_1 &
         + (-7.500000000000001e-01_rk) * xi(12) * xi(16) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(12) * xi(16) * s_1 &
         + (+9.999999999999998e-01_rk) * xi(11) * xi(17) * s_1 &
         + (+4.330127018922193e-01_rk) * xi(10) * xi(15) * c_1 &
         + (-2.499999999999999e-01_rk) * xi(10) * xi(15) * s_1 &
         + (-7.500000000000001e-01_rk) * xi(10) * xi(14) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(10) * xi(14) * s_1 &
         + (-4.330127018922192e-01_rk) * xi(9) * xi(15) * c_1 &
         + (-2.499999999999999e-01_rk) * xi(9) * xi(15) * s_1 &
         + (-7.500000000000000e-01_rk) * xi(9) * xi(14) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(9) * xi(14) * s_1 &
         + (-9.999999999999998e-01_rk) * xi(8) * xi(15) * s_1
    phi_y(26) = (+2.500000000000001e-01_rk) * xi(13) * xi(17) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(13) * xi(17) * s_1 &
         + (-4.330127018922195e-01_rk) * xi(13) * xi(16) * c_1 &
         + (+7.499999999999998e-01_rk) * xi(13) * xi(16) * s_1 &
         + (+2.499999999999999e-01_rk) * xi(12) * xi(17) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(12) * xi(17) * s_1 &
         + (+4.330127018922191e-01_rk) * xi(12) * xi(16) * c_1 &
         + (+7.500000000000000e-01_rk) * xi(12) * xi(16) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(11) * xi(17) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(10) * xi(15) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(10) * xi(15) * s_1 &
         + (-4.330127018922193e-01_rk) * xi(10) * xi(14) * c_1 &
         + (-7.500000000000000e-01_rk) * xi(10) * xi(14) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(9) * xi(15) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(9) * xi(15) * s_1 &
         + (+4.330127018922196e-01_rk) * xi(9) * xi(14) * c_1 &
         + (-7.499999999999998e-01_rk) * xi(9) * xi(14) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(15) * c_1
    phi_x(27) = (-9.999999999999998e-01_rk) * xi(13) * xi(15) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(13) * xi(14) * c_1 &
         + (-1.000000000000000e+00_rk) * xi(12) * xi(15) * s_1 &
         + (-9.999999999999998e-01_rk) * xi(12) * xi(14) * c_1 &
         + (-9.999999999999998e-01_rk) * xi(11) * xi(15) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(11) * xi(14) * c_1 &
         + (+1.000000000000000e+00_rk) * xi(10) * xi(17) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(10) * xi(16) * c_1 &
         + (+9.999999999999997e-01_rk) * xi(9) * xi(17) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(9) * xi(16) * c_1 &
         + (+9.999999999999998e-01_rk) * xi(8) * xi(17) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(8) * xi(16) * c_1
    phi_y(27) = (+9.999999999999998e-01_rk) * xi(13) * xi(15) * c_1 &
         + (-9.999999999999998e-01_rk) * xi(13) * xi(14) * s_1 &
         + (+9.999999999999998e-01_rk) * xi(12) * xi(15) * c_1 &
         + (-9.999999999999997e-01_rk) * xi(12) * xi(14) * s_1 &
         + (+9.999999999999998e-01_rk) * xi(11) * xi(15) * c_1 &
         + (-9.999999999999998e-01_rk) * xi(11) * xi(14) * s_1 &
         + (+9.999999999999998e-01_rk) * xi(10) * xi(17) * c_1 &
         + (+9.999999999999997e-01_rk) * xi(10) * xi(16) * s_1 &
         + (+9.999999999999998e-01_rk) * xi(9) * xi(17) * c_1 &
         + (+9.999999999999998e-01_rk) * xi(9) * xi(16) * s_1 &
         + (+9.999999999999998e-01_rk) * xi(8) * xi(17) * c_1 &
         + (+9.999999999999998e-01_rk) * xi(8) * xi(16) * s_1
    phi_x(28) = (-4.999999999999998e-01_rk) * xi(4) * xi(7) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(4) * xi(6) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(4) * xi(5) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(3) * xi(7) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(3) * xi(6) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(3) * xi(5) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(2) * xi(7) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(2) * xi(6) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(5) * c_3
    phi_y(28) = (+8.660254037844385e-01_rk) * xi(4) * xi(7) * c_3 &
         + (-8.660254037844385e-01_rk) * xi(4) * xi(5) * c_3 &
         + (-8.660254037844385e-01_rk) * xi(3) * xi(6) * c_3 &
         + (+8.660254037844385e-01_rk) * xi(3) * xi(5) * c_3 &
         + (-8.660254037844385e-01_rk) * xi(2) * xi(7) * c_3 &
         + (+8.660254037844385e-01_rk) * xi(2) * xi(6) * c_3
    phi_x(29) = (-4.999999999999998e-01_rk) * xi(10) * xi(13) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(10) * xi(12) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(10) * xi(11) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(9) * xi(13) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(9) * xi(12) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(9) * xi(11) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(8) * xi(13) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(8) * xi(12) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(11) * c_3
    phi_y(29) = (+8.660254037844385e-01_rk) * xi(10) * xi(13) * c_3 &
         + (-8.660254037844385e-01_rk) * xi(10) * xi(11) * c_3 &
         + (-8.660254037844385e-01_rk) * xi(9) * xi(12) * c_3 &
         + (+8.660254037844385e-01_rk) * xi(9) * xi(11) * c_3 &
         + (-8.660254037844385e-01_rk) * xi(8) * xi(13) * c_3 &
         + (+8.660254037844385e-01_rk) * xi(8) * xi(12) * c_3
    phi_x(30) = (-5.000000000000000e-01_rk) * xi(7) * xi(10) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(7) * xi(9) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(7) * xi(8) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(10) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(6) * xi(9) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(6) * xi(8) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(5) * xi(10) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(5) * xi(9) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(8) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(4) * xi(13) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(4) * xi(12) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(4) * xi(11) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(3) * xi(13) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(3) * xi(12) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(3) * xi(11) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(2) * xi(13) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(2) * xi(12) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(11) * c_3
    phi_y(30) = (+8.660254037844386e-01_rk) * xi(7) * xi(10) * c_3 &
         + (-8.660254037844386e-01_rk) * xi(7) * xi(8) * c_3 &
         + (-8.660254037844386e-01_rk) * xi(6) * xi(9) * c_3 &
         + (+8.660254037844386e-01_rk) * xi(6) * xi(8) * c_3 &
         + (-8.660254037844386e-01_rk) * xi(5) * xi(10) * c_3 &
         + (+8.660254037844386e-01_rk) * xi(5) * xi(9) * c_3 &
         + (+8.660254037844386e-01_rk) * xi(4) * xi(13) * c_3 &
         + (-8.660254037844386e-01_rk) * xi(4) * xi(11) * c_3 &
         + (-8.660254037844386e-01_rk) * xi(3) * xi(12) * c_3 &
         + (+8.660254037844386e-01_rk) * xi(3) * xi(11) * c_3 &
         + (-8.660254037844386e-01_rk) * xi(2) * xi(13) * c_3 &
         + (+8.660254037844386e-01_rk) * xi(2) * xi(12) * c_3
    phi_x(31) = (-8.660254037844386e-01_rk) * xi(7) * xi(10) * s_3 &
         + (+8.660254037844386e-01_rk) * xi(7) * xi(8) * s_3 &
         + (+8.660254037844386e-01_rk) * xi(6) * xi(9) * s_3 &
         + (-8.660254037844386e-01_rk) * xi(6) * xi(8) * s_3 &
         + (+8.660254037844386e-01_rk) * xi(5) * xi(10) * s_3 &
         + (-8.660254037844386e-01_rk) * xi(5) * xi(9) * s_3 &
         + (+8.660254037844386e-01_rk) * xi(4) * xi(13) * s_3 &
         + (-8.660254037844386e-01_rk) * xi(4) * xi(11) * s_3 &
         + (-8.660254037844386e-01_rk) * xi(3) * xi(12) * s_3 &
         + (+8.660254037844386e-01_rk) * xi(3) * xi(11) * s_3 &
         + (-8.660254037844386e-01_rk) * xi(2) * xi(13) * s_3 &
         + (+8.660254037844386e-01_rk) * xi(2) * xi(12) * s_3
    phi_y(31) = (-5.000000000000000e-01_rk) * xi(7) * xi(10) * s_3 &
         + (+1.000000000000000e+00_rk) * xi(7) * xi(9) * s_3 &
         + (-5.000000000000000e-01_rk) * xi(7) * xi(8) * s_3 &
         + (+1.000000000000000e+00_rk) * xi(6) * xi(10) * s_3 &
         + (-5.000000000000000e-01_rk) * xi(6) * xi(9) * s_3 &
         + (-5.000000000000000e-01_rk) * xi(6) * xi(8) * s_3 &
         + (-5.000000000000000e-01_rk) * xi(5) * xi(10) * s_3 &
         + (-5.000000000000000e-01_rk) * xi(5) * xi(9) * s_3 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(8) * s_3 &
         + (+5.000000000000000e-01_rk) * xi(4) * xi(13) * s_3 &
         + (-1.000000000000000e+00_rk) * xi(4) * xi(12) * s_3 &
         + (+5.000000000000000e-01_rk) * xi(4) * xi(11) * s_3 &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(13) * s_3 &
         + (+5.000000000000000e-01_rk) * xi(3) * xi(12) * s_3 &
         + (+5.000000000000000e-01_rk) * xi(3) * xi(11) * s_3 &
         + (+5.000000000000000e-01_rk) * xi(2) * xi(13) * s_3 &
         + (+5.000000000000000e-01_rk) * xi(2) * xi(12) * s_3 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(11) * s_3
    phi_x(32) = (+8.660254037844384e-01_rk) * xi(7) * xi(15) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(7) * xi(14) * c_3 &
         + (-8.660254037844385e-01_rk) * xi(6) * xi(15) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(6) * xi(14) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(14) * c_3 &
         + (+8.660254037844385e-01_rk) * xi(4) * xi(17) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(4) * xi(16) * c_3 &
         + (-8.660254037844385e-01_rk) * xi(3) * xi(17) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(3) * xi(16) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(16) * c_3
    phi_y(32) = (-4.999999999999998e-01_rk) * xi(7) * xi(15) * c_3 &
         + (-8.660254037844384e-01_rk) * xi(7) * xi(14) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(6) * xi(15) * c_3 &
         + (+8.660254037844384e-01_rk) * xi(6) * xi(14) * c_3 &
         + (+9.999999999999997e-01_rk) * xi(5) * xi(15) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(4) * xi(17) * c_3 &
         + (-8.660254037844384e-01_rk) * xi(4) * xi(16) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(3) * xi(17) * c_3 &
         + (+8.660254037844384e-01_rk) * xi(3) * xi(16) * c_3 &
         + (+9.999999999999997e-01_rk) * xi(2) * xi(17) * c_3
    phi_x(33) = (-4.999999999999998e-01_rk) * xi(7) * xi(15) * s_3 &
         + (-8.660254037844384e-01_rk) * xi(7) * xi(14) * s_3 &
         + (-4.999999999999998e-01_rk) * xi(6) * xi(15) * s_3 &
         + (+8.660254037844384e-01_rk) * xi(6) * xi(14) * s_3 &
         + (+9.999999999999997e-01_rk) * xi(5) * xi(15) * s_3 &
         + (+4.999999999999998e-01_rk) * xi(4) * xi(17) * s_3 &
         + (+8.660254037844384e-01_rk) * xi(4) * xi(16) * s_3 &
         + (+4.999999999999998e-01_rk) * xi(3) * xi(17) * s_3 &
         + (-8.660254037844384e-01_rk) * xi(3) * xi(16) * s_3 &
         + (-9.999999999999997e-01_rk) * xi(2) * xi(17) * s_3
    phi_y(33) = (-8.660254037844384e-01_rk) * xi(7) * xi(15) * s_3 &
         + (+4.999999999999998e-01_rk) * xi(7) * xi(14) * s_3 &
         + (+8.660254037844385e-01_rk) * xi(6) * xi(15) * s_3 &
         + (+4.999999999999998e-01_rk) * xi(6) * xi(14) * s_3 &
         + (-1.000000000000000e+00_rk) * xi(5) * xi(14) * s_3 &
         + (+8.660254037844385e-01_rk) * xi(4) * xi(17) * s_3 &
         + (-4.999999999999998e-01_rk) * xi(4) * xi(16) * s_3 &
         + (-8.660254037844385e-01_rk) * xi(3) * xi(17) * s_3 &
         + (-4.999999999999998e-01_rk) * xi(3) * xi(16) * s_3 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(16) * s_3
    phi_x(34) = (+8.660254037844384e-01_rk) * xi(13) * xi(15) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(13) * xi(14) * c_3 &
         + (-8.660254037844385e-01_rk) * xi(12) * xi(15) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(12) * xi(14) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(11) * xi(14) * c_3 &
         + (+8.660254037844385e-01_rk) * xi(10) * xi(17) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(10) * xi(16) * c_3 &
         + (-8.660254037844385e-01_rk) * xi(9) * xi(17) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(9) * xi(16) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(16) * c_3
    phi_y(34) = (-4.999999999999998e-01_rk) * xi(13) * xi(15) * c_3 &
         + (-8.660254037844384e-01_rk) * xi(13) * xi(14) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(12) * xi(15) * c_3 &
         + (+8.660254037844384e-01_rk) * xi(12) * xi(14) * c_3 &
         + (+9.999999999999997e-01_rk) * xi(11) * xi(15) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(10) * xi(17) * c_3 &
         + (-8.660254037844384e-01_rk) * xi(10) * xi(16) * c_3 &
         + (-4.999999999999998e-01_rk) * xi(9) * xi(17) * c_3 &
         + (+8.660254037844384e-01_rk) * xi(9) * xi(16) * c_3 &
         + (+9.999999999999997e-01_rk) * xi(8) * xi(17) * c_3
    phi_x(35) = (-4.999999999999998e-01_rk) * xi(13) * xi(15) * s_3 &
         + (-8.660254037844384e-01_rk) * xi(13) * xi(14) * s_3 &
         + (-4.999999999999998e-01_rk) * xi(12) * xi(15) * s_3 &
         + (+8.660254037844384e-01_rk) * xi(12) * xi(14) * s_3 &
         + (+9.999999999999997e-01_rk) * xi(11) * xi(15) * s_3 &
         + (+4.999999999999998e-01_rk) * xi(10) * xi(17) * s_3 &
         + (+8.660254037844384e-01_rk) * xi(10) * xi(16) * s_3 &
         + (+4.999999999999998e-01_rk) * xi(9) * xi(17) * s_3 &
         + (-8.660254037844384e-01_rk) * xi(9) * xi(16) * s_3 &
         + (-9.999999999999997e-01_rk) * xi(8) * xi(17) * s_3
    phi_y(35) = (-8.660254037844384e-01_rk) * xi(13) * xi(15) * s_3 &
         + (+4.999999999999998e-01_rk) * xi(13) * xi(14) * s_3 &
         + (+8.660254037844385e-01_rk) * xi(12) * xi(15) * s_3 &
         + (+4.999999999999998e-01_rk) * xi(12) * xi(14) * s_3 &
         + (-1.000000000000000e+00_rk) * xi(11) * xi(14) * s_3 &
         + (+8.660254037844385e-01_rk) * xi(10) * xi(17) * s_3 &
         + (-4.999999999999998e-01_rk) * xi(10) * xi(16) * s_3 &
         + (-8.660254037844385e-01_rk) * xi(9) * xi(17) * s_3 &
         + (-4.999999999999998e-01_rk) * xi(9) * xi(16) * s_3 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(16) * s_3
    phi_x(36) = (+1.000000000000000e+00_rk) * xi(15) * xi(17) * c_3 &
         + (-9.999999999999994e-01_rk) * xi(14) * xi(16) * c_3
    phi_y(36) = (-9.999999999999999e-01_rk) * xi(15) * xi(16) * c_3 &
         + (-1.000000000000000e+00_rk) * xi(14) * xi(17) * c_3
    phi_x(37) = (-4.999999999999994e-01_rk) * xi(7)**2 * c_5 &
         + (-8.660254037844388e-01_rk) * xi(7)**2 * s_5 &
         + (-5.000000000000006e-01_rk) * xi(6)**2 * c_5 &
         + (+8.660254037844385e-01_rk) * xi(6)**2 * s_5 &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * c_5 &
         + (-5.000000000000006e-01_rk) * xi(4)**2 * c_5 &
         + (+8.660254037844385e-01_rk) * xi(4)**2 * s_5 &
         + (-4.999999999999989e-01_rk) * xi(3)**2 * c_5 &
         + (-8.660254037844389e-01_rk) * xi(3)**2 * s_5 &
         + (+1.000000000000000e+00_rk) * xi(2)**2 * c_5
    phi_y(37) = (+8.660254037844389e-01_rk) * xi(7)**2 * c_5 &
         + (-4.999999999999997e-01_rk) * xi(7)**2 * s_5 &
         + (-8.660254037844385e-01_rk) * xi(6)**2 * c_5 &
         + (-5.000000000000001e-01_rk) * xi(6)**2 * s_5 &
         + (+9.999999999999998e-01_rk) * xi(5)**2 * s_5 &
         + (+8.660254037844385e-01_rk) * xi(4)**2 * c_5 &
         + (+5.000000000000001e-01_rk) * xi(4)**2 * s_5 &
         + (-8.660254037844389e-01_rk) * xi(3)**2 * c_5 &
         + (+4.999999999999997e-01_rk) * xi(3)**2 * s_5 &
         + (-9.999999999999998e-01_rk) * xi(2)**2 * s_5
    phi_x(38) = (-4.999999999999994e-01_rk) * xi(13)**2 * c_5 &
         + (-8.660254037844388e-01_rk) * xi(13)**2 * s_5 &
         + (-5.000000000000006e-01_rk) * xi(12)**2 * c_5 &
         + (+8.660254037844385e-01_rk) * xi(12)**2 * s_5 &
         + (+1.000000000000000e+00_rk) * xi(11)**2 * c_5 &
         + (-5.000000000000006e-01_rk) * xi(10)**2 * c_5 &
         + (+8.660254037844385e-01_rk) * xi(10)**2 * s_5 &
         + (-4.999999999999989e-01_rk) * xi(9)**2 * c_5 &
         + (-8.660254037844389e-01_rk) * xi(9)**2 * s_5 &
         + (+1.000000000000000e+00_rk) * xi(8)**2 * c_5
    phi_y(38) = (+8.660254037844389e-01_rk) * xi(13)**2 * c_5 &
         + (-4.999999999999997e-01_rk) * xi(13)**2 * s_5 &
         + (-8.660254037844385e-01_rk) * xi(12)**2 * c_5 &
         + (-5.000000000000001e-01_rk) * xi(12)**2 * s_5 &
         + (+9.999999999999998e-01_rk) * xi(11)**2 * s_5 &
         + (+8.660254037844385e-01_rk) * xi(10)**2 * c_5 &
         + (+5.000000000000001e-01_rk) * xi(10)**2 * s_5 &
         + (-8.660254037844389e-01_rk) * xi(9)**2 * c_5 &
         + (+4.999999999999997e-01_rk) * xi(9)**2 * s_5 &
         + (-9.999999999999998e-01_rk) * xi(8)**2 * s_5
    phi_x(39) = (-4.999999999999999e-01_rk) * xi(17)**2 * c_5 &
         + (-9.999999999999998e-01_rk) * xi(16) * xi(17) * s_5 &
         + (+5.000000000000000e-01_rk) * xi(16)**2 * c_5 &
         + (-4.999999999999997e-01_rk) * xi(15)**2 * c_5 &
         + (+1.000000000000000e+00_rk) * xi(14) * xi(15) * s_5 &
         + (+5.000000000000001e-01_rk) * xi(14)**2 * c_5
    phi_y(39) = (-4.999999999999997e-01_rk) * xi(17)**2 * s_5 &
         + (+9.999999999999999e-01_rk) * xi(16) * xi(17) * c_5 &
         + (+4.999999999999998e-01_rk) * xi(16)**2 * s_5 &
         + (+4.999999999999997e-01_rk) * xi(15)**2 * s_5 &
         + (+9.999999999999999e-01_rk) * xi(14) * xi(15) * c_5 &
         + (-4.999999999999998e-01_rk) * xi(14)**2 * s_5
    phi_x(40) = (-4.999999999999994e-01_rk) * xi(1) * xi(7) * c_5 &
         + (-8.660254037844388e-01_rk) * xi(1) * xi(7) * s_5 &
         + (-5.000000000000006e-01_rk) * xi(1) * xi(6) * c_5 &
         + (+8.660254037844385e-01_rk) * xi(1) * xi(6) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(5) * c_5 &
         + (-5.000000000000006e-01_rk) * xi(1) * xi(4) * c_5 &
         + (+8.660254037844385e-01_rk) * xi(1) * xi(4) * s_5 &
         + (-4.999999999999989e-01_rk) * xi(1) * xi(3) * c_5 &
         + (-8.660254037844389e-01_rk) * xi(1) * xi(3) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(2) * c_5
    phi_y(40) = (+8.660254037844389e-01_rk) * xi(1) * xi(7) * c_5 &
         + (-4.999999999999997e-01_rk) * xi(1) * xi(7) * s_5 &
         + (-8.660254037844385e-01_rk) * xi(1) * xi(6) * c_5 &
         + (-5.000000000000001e-01_rk) * xi(1) * xi(6) * s_5 &
         + (+9.999999999999998e-01_rk) * xi(1) * xi(5) * s_5 &
         + (+8.660254037844385e-01_rk) * xi(1) * xi(4) * c_5 &
         + (+5.000000000000001e-01_rk) * xi(1) * xi(4) * s_5 &
         + (-8.660254037844389e-01_rk) * xi(1) * xi(3) * c_5 &
         + (+4.999999999999997e-01_rk) * xi(1) * xi(3) * s_5 &
         + (-9.999999999999998e-01_rk) * xi(1) * xi(2) * s_5
    phi_x(41) = (-4.999999999999994e-01_rk) * xi(1) * xi(13) * c_5 &
         + (-8.660254037844388e-01_rk) * xi(1) * xi(13) * s_5 &
         + (-5.000000000000006e-01_rk) * xi(1) * xi(12) * c_5 &
         + (+8.660254037844385e-01_rk) * xi(1) * xi(12) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(11) * c_5 &
         + (-5.000000000000006e-01_rk) * xi(1) * xi(10) * c_5 &
         + (+8.660254037844385e-01_rk) * xi(1) * xi(10) * s_5 &
         + (-4.999999999999989e-01_rk) * xi(1) * xi(9) * c_5 &
         + (-8.660254037844389e-01_rk) * xi(1) * xi(9) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(8) * c_5
    phi_y(41) = (+8.660254037844389e-01_rk) * xi(1) * xi(13) * c_5 &
         + (-4.999999999999997e-01_rk) * xi(1) * xi(13) * s_5 &
         + (-8.660254037844385e-01_rk) * xi(1) * xi(12) * c_5 &
         + (-5.000000000000001e-01_rk) * xi(1) * xi(12) * s_5 &
         + (+9.999999999999998e-01_rk) * xi(1) * xi(11) * s_5 &
         + (+8.660254037844385e-01_rk) * xi(1) * xi(10) * c_5 &
         + (+5.000000000000001e-01_rk) * xi(1) * xi(10) * s_5 &
         + (-8.660254037844389e-01_rk) * xi(1) * xi(9) * c_5 &
         + (+4.999999999999997e-01_rk) * xi(1) * xi(9) * s_5 &
         + (-9.999999999999998e-01_rk) * xi(1) * xi(8) * s_5
    phi_x(42) = (-9.999999999999998e-01_rk) * xi(1) * xi(17) * s_5 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(16) * c_5 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(15) * s_5 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(14) * c_5
    phi_y(42) = (+9.999999999999997e-01_rk) * xi(1) * xi(17) * c_5 &
         + (-9.999999999999996e-01_rk) * xi(1) * xi(16) * s_5 &
         + (+9.999999999999997e-01_rk) * xi(1) * xi(15) * c_5 &
         + (+9.999999999999996e-01_rk) * xi(1) * xi(14) * s_5
    phi_x(43) = (+1.000000000000000e+00_rk) * xi(6) * xi(7) * c_5 &
         + (-5.000000000000011e-01_rk) * xi(5) * xi(7) * c_5 &
         + (+8.660254037844356e-01_rk) * xi(5) * xi(7) * s_5 &
         + (-4.999999999999988e-01_rk) * xi(5) * xi(6) * c_5 &
         + (-8.660254037844353e-01_rk) * xi(5) * xi(6) * s_5 &
         + (+9.999999999999987e-01_rk) * xi(3) * xi(4) * c_5 &
         + (-4.999999999999988e-01_rk) * xi(2) * xi(4) * c_5 &
         + (-8.660254037844353e-01_rk) * xi(2) * xi(4) * s_5 &
         + (-5.000000000000011e-01_rk) * xi(2) * xi(3) * c_5 &
         + (+8.660254037844356e-01_rk) * xi(2) * xi(3) * s_5
    phi_y(43) = (+9.999999999999963e-01_rk) * xi(6) * xi(7) * s_5 &
         + (-8.660254037844383e-01_rk) * xi(5) * xi(7) * c_5 &
         + (-4.999999999999977e-01_rk) * xi(5) * xi(7) * s_5 &
         + (+8.660254037844396e-01_rk) * xi(5) * xi(6) * c_5 &
         + (-4.999999999999985e-01_rk) * xi(5) * xi(6) * s_5 &
         + (-9.999999999999963e-01_rk) * xi(3) * xi(4) * s_5 &
         + (-8.660254037844396e-01_rk) * xi(2) * xi(4) * c_5 &
         + (+4.999999999999985e-01_rk) * xi(2) * xi(4) * s_5 &
         + (+8.660254037844383e-01_rk) * xi(2) * xi(3) * c_5 &
         + (+4.999999999999977e-01_rk) * xi(2) * xi(3) * s_5
    phi_x(44) = (-4.999999999999998e-01_rk) * xi(4) * xi(7) * c_5 &
         + (-5.000000000000006e-01_rk) * xi(4) * xi(6) * c_5 &
         + (+8.660254037844380e-01_rk) * xi(4) * xi(6) * s_5 &
         + (+2.499999999999997e-01_rk) * xi(4) * xi(5) * c_5 &
         + (+4.330127018922193e-01_rk) * xi(4) * xi(5) * s_5 &
         + (-4.999999999999990e-01_rk) * xi(3) * xi(7) * c_5 &
         + (-8.660254037844389e-01_rk) * xi(3) * xi(7) * s_5 &
         + (-4.999999999999998e-01_rk) * xi(3) * xi(6) * c_5 &
         + (+2.500000000000003e-01_rk) * xi(3) * xi(5) * c_5 &
         + (-4.330127018922190e-01_rk) * xi(3) * xi(5) * s_5 &
         + (+2.500000000000003e-01_rk) * xi(2) * xi(7) * c_5 &
         + (-4.330127018922189e-01_rk) * xi(2) * xi(7) * s_5 &
         + (+2.499999999999997e-01_rk) * xi(2) * xi(6) * c_5 &
         + (+4.330127018922193e-01_rk) * xi(2) * xi(6) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(5) * c_5
    phi_y(44) = (+8.660254037844385e-01_rk) * xi(4) * xi(7) * c_5 &
         + (+4.330127018922188e-01_rk) * xi(4) * xi(5) * c_5 &
         + (+7.500000000000001e-01_rk) * xi(4) * xi(5) * s_5 &
         + (-8.660254037844385e-01_rk) * xi(3) * xi(6) * c_5 &
         + (-4.330127018922199e-01_rk) * xi(3) * xi(5) * c_5 &
         + (+7.499999999999993e-01_rk) * xi(3) * xi(5) * s_5 &
         + (+4.330127018922199e-01_rk) * xi(2) * xi(7) * c_5 &
         + (-7.499999999999993e-01_rk) * xi(2) * xi(7) * s_5 &
         + (-4.330127018922188e-01_rk) * xi(2) * xi(6) * c_5 &
         + (-7.500000000000001e-01_rk) * xi(2) * xi(6) * s_5
    phi_x(45) = (+1.000000000000000e+00_rk) * xi(12) * xi(13) * c_5 &
         + (-5.000000000000011e-01_rk) * xi(11) * xi(13) * c_5 &
         + (+8.660254037844356e-01_rk) * xi(11) * xi(13) * s_5 &
         + (-4.999999999999988e-01_rk) * xi(11) * xi(12) * c_5 &
         + (-8.660254037844353e-01_rk) * xi(11) * xi(12) * s_5 &
         + (+9.999999999999987e-01_rk) * xi(9) * xi(10) * c_5 &
         + (-4.999999999999988e-01_rk) * xi(8) * xi(10) * c_5 &
         + (-8.660254037844353e-01_rk) * xi(8) * xi(10) * s_5 &
         + (-5.000000000000011e-01_rk) * xi(8) * xi(9) * c_5 &
         + (+8.660254037844356e-01_rk) * xi(8) * xi(9) * s_5
    phi_y(45) = (+9.999999999999963e-01_rk) * xi(12) * xi(13) * s_5 &
         + (-8.660254037844383e-01_rk) * xi(11) * xi(13) * c_5 &
         + (-4.999999999999977e-01_rk) * xi(11) * xi(13) * s_5 &
         + (+8.660254037844396e-01_rk) * xi(11) * xi(12) * c_5 &
         + (-4.999999999999985e-01_rk) * xi(11) * xi(12) * s_5 &
         + (-9.999999999999963e-01_rk) * xi(9) * xi(10) * s_5 &
         + (-8.660254037844396e-01_rk) * xi(8) * xi(10) * c_5 &
         + (+4.999999999999985e-01_rk) * xi(8) * xi(10) * s_5 &
         + (+8.660254037844383e-01_rk) * xi(8) * xi(9) * c_5 &
         + (+4.999999999999977e-01_rk) * xi(8) * xi(9) * s_5
    phi_x(46) = (-4.999999999999998e-01_rk) * xi(10) * xi(13) * c_5 &
         + (-5.000000000000006e-01_rk) * xi(10) * xi(12) * c_5 &
         + (+8.660254037844380e-01_rk) * xi(10) * xi(12) * s_5 &
         + (+2.499999999999997e-01_rk) * xi(10) * xi(11) * c_5 &
         + (+4.330127018922193e-01_rk) * xi(10) * xi(11) * s_5 &
         + (-4.999999999999990e-01_rk) * xi(9) * xi(13) * c_5 &
         + (-8.660254037844389e-01_rk) * xi(9) * xi(13) * s_5 &
         + (-4.999999999999998e-01_rk) * xi(9) * xi(12) * c_5 &
         + (+2.500000000000003e-01_rk) * xi(9) * xi(11) * c_5 &
         + (-4.330127018922190e-01_rk) * xi(9) * xi(11) * s_5 &
         + (+2.500000000000003e-01_rk) * xi(8) * xi(13) * c_5 &
         + (-4.330127018922189e-01_rk) * xi(8) * xi(13) * s_5 &
         + (+2.499999999999997e-01_rk) * xi(8) * xi(12) * c_5 &
         + (+4.330127018922193e-01_rk) * xi(8) * xi(12) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(11) * c_5
    phi_y(46) = (+8.660254037844385e-01_rk) * xi(10) * xi(13) * c_5 &
         + (+4.330127018922188e-01_rk) * xi(10) * xi(11) * c_5 &
         + (+7.500000000000001e-01_rk) * xi(10) * xi(11) * s_5 &
         + (-8.660254037844385e-01_rk) * xi(9) * xi(12) * c_5 &
         + (-4.330127018922199e-01_rk) * xi(9) * xi(11) * c_5 &
         + (+7.499999999999993e-01_rk) * xi(9) * xi(11) * s_5 &
         + (+4.330127018922199e-01_rk) * xi(8) * xi(13) * c_5 &
         + (-7.499999999999993e-01_rk) * xi(8) * xi(13) * s_5 &
         + (-4.330127018922188e-01_rk) * xi(8) * xi(12) * c_5 &
         + (-7.500000000000001e-01_rk) * xi(8) * xi(12) * s_5
    phi_x(47) = (-4.999999999999994e-01_rk) * xi(7) * xi(13) * c_5 &
         + (-8.660254037844388e-01_rk) * xi(7) * xi(13) * s_5 &
         + (-5.000000000000006e-01_rk) * xi(6) * xi(12) * c_5 &
         + (+8.660254037844385e-01_rk) * xi(6) * xi(12) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(11) * c_5 &
         + (-5.000000000000006e-01_rk) * xi(4) * xi(10) * c_5 &
         + (+8.660254037844385e-01_rk) * xi(4) * xi(10) * s_5 &
         + (-4.999999999999989e-01_rk) * xi(3) * xi(9) * c_5 &
         + (-8.660254037844389e-01_rk) * xi(3) * xi(9) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(8) * c_5
    phi_y(47) = (+8.660254037844389e-01_rk) * xi(7) * xi(13) * c_5 &
         + (-4.999999999999997e-01_rk) * xi(7) * xi(13) * s_5 &
         + (-8.660254037844385e-01_rk) * xi(6) * xi(12) * c_5 &
         + (-5.000000000000001e-01_rk) * xi(6) * xi(12) * s_5 &
         + (+9.999999999999998e-01_rk) * xi(5) * xi(11) * s_5 &
         + (+8.660254037844385e-01_rk) * xi(4) * xi(10) * c_5 &
         + (+5.000000000000001e-01_rk) * xi(4) * xi(10) * s_5 &
         + (-8.660254037844389e-01_rk) * xi(3) * xi(9) * c_5 &
         + (+4.999999999999997e-01_rk) * xi(3) * xi(9) * s_5 &
         + (-9.999999999999998e-01_rk) * xi(2) * xi(8) * s_5
    phi_x(48) = (-5.000000000000021e-01_rk) * xi(7) * xi(12) * c_5 &
         + (-8.660254037844374e-01_rk) * xi(7) * xi(12) * s_5 &
         + (-4.999999999999963e-01_rk) * xi(7) * xi(11) * c_5 &
         + (-8.660254037844390e-01_rk) * xi(7) * xi(11) * s_5 &
         + (-5.000000000000029e-01_rk) * xi(6) * xi(13) * c_5 &
         + (+8.660254037844379e-01_rk) * xi(6) * xi(13) * s_5 &
         + (-4.999999999999978e-01_rk) * xi(6) * xi(11) * c_5 &
         + (+8.660254037844383e-01_rk) * xi(6) * xi(11) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(13) * c_5 &
         + (+9.999999999999992e-01_rk) * xi(5) * xi(12) * c_5 &
         + (-5.000000000000029e-01_rk) * xi(4) * xi(9) * c_5 &
         + (+8.660254037844379e-01_rk) * xi(4) * xi(9) * s_5 &
         + (-4.999999999999978e-01_rk) * xi(4) * xi(8) * c_5 &
         + (+8.660254037844383e-01_rk) * xi(4) * xi(8) * s_5 &
         + (-5.000000000000009e-01_rk) * xi(3) * xi(10) * c_5 &
         + (-8.660254037844380e-01_rk) * xi(3) * xi(10) * s_5 &
         + (-4.999999999999963e-01_rk) * xi(3) * xi(8) * c_5 &
         + (-8.660254037844390e-01_rk) * xi(3) * xi(8) * s_5 &
         + (+9.999999999999992e-01_rk) * xi(2) * xi(10) * c_5 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(9) * c_5
    phi_y(48) = (+8.660254037844374e-01_rk) * xi(7) * xi(12) * c_5 &
         + (-5.000000000000001e-01_rk) * xi(7) * xi(12) * s_5 &
         + (+8.660254037844399e-01_rk) * xi(7) * xi(11) * c_5 &
         + (-4.999999999999988e-01_rk) * xi(7) * xi(11) * s_5 &
         + (-8.660254037844359e-01_rk) * xi(6) * xi(13) * c_5 &
         + (-5.000000000000010e-01_rk) * xi(6) * xi(13) * s_5 &
         + (-8.660254037844399e-01_rk) * xi(6) * xi(11) * c_5 &
         + (-4.999999999999988e-01_rk) * xi(6) * xi(11) * s_5 &
         + (+9.999999999999990e-01_rk) * xi(5) * xi(13) * s_5 &
         + (+9.999999999999999e-01_rk) * xi(5) * xi(12) * s_5 &
         + (+8.660254037844359e-01_rk) * xi(4) * xi(9) * c_5 &
         + (+5.000000000000010e-01_rk) * xi(4) * xi(9) * s_5 &
         + (+8.660254037844399e-01_rk) * xi(4) * xi(8) * c_5 &
         + (+4.999999999999988e-01_rk) * xi(4) * xi(8) * s_5 &
         + (-8.660254037844374e-01_rk) * xi(3) * xi(10) * c_5 &
         + (+5.000000000000001e-01_rk) * xi(3) * xi(10) * s_5 &
         + (-8.660254037844399e-01_rk) * xi(3) * xi(8) * c_5 &
         + (+4.999999999999988e-01_rk) * xi(3) * xi(8) * s_5 &
         + (-9.999999999999999e-01_rk) * xi(2) * xi(10) * s_5 &
         + (-9.999999999999990e-01_rk) * xi(2) * xi(9) * s_5
    phi_x(49) = (-8.660254037844374e-01_rk) * xi(7) * xi(12) * c_5 &
         + (+5.000000000000001e-01_rk) * xi(7) * xi(12) * s_5 &
         + (+8.660254037844399e-01_rk) * xi(7) * xi(11) * c_5 &
         + (-4.999999999999988e-01_rk) * xi(7) * xi(11) * s_5 &
         + (-8.660254037844359e-01_rk) * xi(6) * xi(13) * c_5 &
         + (-5.000000000000010e-01_rk) * xi(6) * xi(13) * s_5 &
         + (+8.660254037844399e-01_rk) * xi(6) * xi(11) * c_5 &
         + (+4.999999999999988e-01_rk) * xi(6) * xi(11) * s_5 &
         + (-9.999999999999990e-01_rk) * xi(5) * xi(13) * s_5 &
         + (+9.999999999999999e-01_rk) * xi(5) * xi(12) * s_5 &
         + (-8.660254037844359e-01_rk) * xi(4) * xi(9) * c_5 &
         + (-5.000000000000010e-01_rk) * xi(4) * xi(9) * s_5 &
         + (+8.660254037844399e-01_rk) * xi(4) * xi(8) * c_5 &
         + (+4.999999999999988e-01_rk) * xi(4) * xi(8) * s_5 &
         + (-8.660254037844374e-01_rk) * xi(3) * xi(10) * c_5 &
         + (+5.000000000000001e-01_rk) * xi(3) * xi(10) * s_5 &
         + (+8.660254037844399e-01_rk) * xi(3) * xi(8) * c_5 &
         + (-4.999999999999988e-01_rk) * xi(3) * xi(8) * s_5 &
         + (+9.999999999999999e-01_rk) * xi(2) * xi(10) * s_5 &
         + (-9.999999999999990e-01_rk) * xi(2) * xi(9) * s_5
    phi_y(49) = (-5.000000000000021e-01_rk) * xi(7) * xi(12) * c_5 &
         + (-8.660254037844374e-01_rk) * xi(7) * xi(12) * s_5 &
         + (+4.999999999999963e-01_rk) * xi(7) * xi(11) * c_5 &
         + (+8.660254037844390e-01_rk) * xi(7) * xi(11) * s_5 &
         + (+5.000000000000029e-01_rk) * xi(6) * xi(13) * c_5 &
         + (-8.660254037844379e-01_rk) * xi(6) * xi(13) * s_5 &
         + (-4.999999999999978e-01_rk) * xi(6) * xi(11) * c_5 &
         + (+8.660254037844383e-01_rk) * xi(6) * xi(11) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(13) * c_5 &
         + (-9.999999999999992e-01_rk) * xi(5) * xi(12) * c_5 &
         + (-5.000000000000029e-01_rk) * xi(4) * xi(9) * c_5 &
         + (+8.660254037844379e-01_rk) * xi(4) * xi(9) * s_5 &
         + (+4.999999999999978e-01_rk) * xi(4) * xi(8) * c_5 &
         + (-8.660254037844383e-01_rk) * xi(4) * xi(8) * s_5 &
         + (+5.000000000000009e-01_rk) * xi(3) * xi(10) * c_5 &
         + (+8.660254037844380e-01_rk) * xi(3) * xi(10) * s_5 &
         + (-4.999999999999963e-01_rk) * xi(3) * xi(8) * c_5 &
         + (-8.660254037844390e-01_rk) * xi(3) * xi(8) * s_5 &
         + (+9.999999999999992e-01_rk) * xi(2) * xi(10) * c_5 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(9) * c_5
    phi_x(50) = (-5.000000000000000e-01_rk) * xi(7) * xi(10) * c_5 &
         + (-4.999999999999996e-01_rk) * xi(7) * xi(9) * c_5 &
         + (-8.660254037844388e-01_rk) * xi(7) * xi(9) * s_5 &
         + (+2.500000000000004e-01_rk) * xi(7) * xi(8) * c_5 &
         + (-4.330127018922191e-01_rk) * xi(7) * xi(8) * s_5 &
         + (-5.000000000000008e-01_rk) * xi(6) * xi(10) * c_5 &
         + (+8.660254037844383e-01_rk) * xi(6) * xi(10) * s_5 &
         + (-5.000000000000000e-01_rk) * xi(6) * xi(9) * c_5 &
         + (+2.499999999999998e-01_rk) * xi(6) * xi(8) * c_5 &
         + (+4.330127018922194e-01_rk) * xi(6) * xi(8) * s_5 &
         + (+2.499999999999998e-01_rk) * xi(5) * xi(10) * c_5 &
         + (+4.330127018922194e-01_rk) * xi(5) * xi(10) * s_5 &
         + (+2.500000000000003e-01_rk) * xi(5) * xi(9) * c_5 &
         + (-4.330127018922191e-01_rk) * xi(5) * xi(9) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(8) * c_5 &
         + (-5.000000000000000e-01_rk) * xi(4) * xi(13) * c_5 &
         + (-5.000000000000008e-01_rk) * xi(4) * xi(12) * c_5 &
         + (+8.660254037844383e-01_rk) * xi(4) * xi(12) * s_5 &
         + (+2.499999999999998e-01_rk) * xi(4) * xi(11) * c_5 &
         + (+4.330127018922194e-01_rk) * xi(4) * xi(11) * s_5 &
         + (-4.999999999999986e-01_rk) * xi(3) * xi(13) * c_5 &
         + (-8.660254037844394e-01_rk) * xi(3) * xi(13) * s_5 &
         + (-5.000000000000000e-01_rk) * xi(3) * xi(12) * c_5 &
         + (+2.500000000000004e-01_rk) * xi(3) * xi(11) * c_5 &
         + (-4.330127018922191e-01_rk) * xi(3) * xi(11) * s_5 &
         + (+2.500000000000003e-01_rk) * xi(2) * xi(13) * c_5 &
         + (-4.330127018922191e-01_rk) * xi(2) * xi(13) * s_5 &
         + (+2.499999999999998e-01_rk) * xi(2) * xi(12) * c_5 &
         + (+4.330127018922194e-01_rk) * xi(2) * xi(12) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(11) * c_5
    phi_y(50) = (+8.660254037844386e-01_rk) * xi(7) * xi(10) * c_5 &
         + (+4.330127018922200e-01_rk) * xi(7) * xi(8) * c_5 &
         + (-7.499999999999997e-01_rk) * xi(7) * xi(8) * s_5 &
         + (-8.660254037844386e-01_rk) * xi(6) * xi(9) * c_5 &
         + (-4.330127018922190e-01_rk) * xi(6) * xi(8) * c_5 &
         + (-7.500000000000001e-01_rk) * xi(6) * xi(8) * s_5 &
         + (+4.330127018922190e-01_rk) * xi(5) * xi(10) * c_5 &
         + (+7.500000000000001e-01_rk) * xi(5) * xi(10) * s_5 &
         + (-4.330127018922200e-01_rk) * xi(5) * xi(9) * c_5 &
         + (+7.499999999999997e-01_rk) * xi(5) * xi(9) * s_5 &
         + (+8.660254037844386e-01_rk) * xi(4) * xi(13) * c_5 &
         + (+4.330127018922190e-01_rk) * xi(4) * xi(11) * c_5 &
         + (+7.500000000000001e-01_rk) * xi(4) * xi(11) * s_5 &
         + (-8.660254037844386e-01_rk) * xi(3) * xi(12) * c_5 &
         + (-4.330127018922200e-01_rk) * xi(3) * xi(11) * c_5 &
         + (+7.499999999999997e-01_rk) * xi(3) * xi(11) * s_5 &
         + (+4.330127018922200e-01_rk) * xi(2) * xi(13) * c_5 &
         + (-7.499999999999997e-01_rk) * xi(2) * xi(13) * s_5 &
         + (-4.330127018922190e-01_rk) * xi(2) * xi(12) * c_5 &
         + (-7.500000000000001e-01_rk) * xi(2) * xi(12) * s_5
    phi_x(51) = (-8.660254037844386e-01_rk) * xi(7) * xi(10) * s_5 &
         + (-7.499999999999997e-01_rk) * xi(7) * xi(8) * c_5 &
         + (-4.330127018922200e-01_rk) * xi(7) * xi(8) * s_5 &
         + (+8.660254037844386e-01_rk) * xi(6) * xi(9) * s_5 &
         + (-7.500000000000001e-01_rk) * xi(6) * xi(8) * c_5 &
         + (+4.330127018922190e-01_rk) * xi(6) * xi(8) * s_5 &
         + (+7.500000000000001e-01_rk) * xi(5) * xi(10) * c_5 &
         + (-4.330127018922190e-01_rk) * xi(5) * xi(10) * s_5 &
         + (+7.499999999999997e-01_rk) * xi(5) * xi(9) * c_5 &
         + (+4.330127018922200e-01_rk) * xi(5) * xi(9) * s_5 &
         + (+8.660254037844386e-01_rk) * xi(4) * xi(13) * s_5 &
         + (-7.500000000000001e-01_rk) * xi(4) * xi(11) * c_5 &
         + (+4.330127018922190e-01_rk) * xi(4) * xi(11) * s_5 &
         + (-8.660254037844386e-01_rk) * xi(3) * xi(12) * s_5 &
         + (-7.499999999999997e-01_rk) * xi(3) * xi(11) * c_5 &
         + (-4.330127018922200e-01_rk) * xi(3) * xi(11) * s_5 &
         + (+7.499999999999997e-01_rk) * xi(2) * xi(13) * c_5 &
         + (+4.330127018922200e-01_rk) * xi(2) * xi(13) * s_5 &
         + (+7.500000000000001e-01_rk) * xi(2) * xi(12) * c_5 &
         + (-4.330127018922190e-01_rk) * xi(2) * xi(12) * s_5
    phi_y(51) = (-5.000000000000000e-01_rk) * xi(7) * xi(10) * s_5 &
         + (+8.660254037844388e-01_rk) * xi(7) * xi(9) * c_5 &
         + (-4.999999999999996e-01_rk) * xi(7) * xi(9) * s_5 &
         + (+4.330127018922191e-01_rk) * xi(7) * xi(8) * c_5 &
         + (+2.500000000000004e-01_rk) * xi(7) * xi(8) * s_5 &
         + (-8.660254037844383e-01_rk) * xi(6) * xi(10) * c_5 &
         + (-5.000000000000008e-01_rk) * xi(6) * xi(10) * s_5 &
         + (-5.000000000000000e-01_rk) * xi(6) * xi(9) * s_5 &
         + (-4.330127018922194e-01_rk) * xi(6) * xi(8) * c_5 &
         + (+2.499999999999998e-01_rk) * xi(6) * xi(8) * s_5 &
         + (-4.330127018922194e-01_rk) * xi(5) * xi(10) * c_5 &
         + (+2.499999999999998e-01_rk) * xi(5) * xi(10) * s_5 &
         + (+4.330127018922191e-01_rk) * xi(5) * xi(9) * c_5 &
         + (+2.500000000000003e-01_rk) * xi(5) * xi(9) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(8) * s_5 &
         + (+5.000000000000000e-01_rk) * xi(4) * xi(13) * s_5 &
         + (+8.660254037844383e-01_rk) * xi(4) * xi(12) * c_5 &
         + (+5.000000000000008e-01_rk) * xi(4) * xi(12) * s_5 &
         + (+4.330127018922194e-01_rk) * xi(4) * xi(11) * c_5 &
         + (-2.499999999999998e-01_rk) * xi(4) * xi(11) * s_5 &
         + (-8.660254037844394e-01_rk) * xi(3) * xi(13) * c_5 &
         + (+4.999999999999986e-01_rk) * xi(3) * xi(13) * s_5 &
         + (+5.000000000000000e-01_rk) * xi(3) * xi(12) * s_5 &
         + (-4.330127018922191e-01_rk) * xi(3) * xi(11) * c_5 &
         + (-2.500000000000004e-01_rk) * xi(3) * xi(11) * s_5 &
         + (-4.330127018922191e-01_rk) * xi(2) * xi(13) * c_5 &
         + (-2.500000000000003e-01_rk) * xi(2) * xi(13) * s_5 &
         + (+4.330127018922194e-01_rk) * xi(2) * xi(12) * c_5 &
         + (-2.499999999999998e-01_rk) * xi(2) * xi(12) * s_5 &
         + (-1.000000000000000e+00_rk) * xi(2) * xi(11) * s_5
    phi_x(52) = (+4.330127018922189e-01_rk) * xi(7) * xi(17) * c_5 &
         + (+7.499999999999998e-01_rk) * xi(7) * xi(17) * s_5 &
         + (+2.499999999999997e-01_rk) * xi(7) * xi(16) * c_5 &
         + (+4.330127018922193e-01_rk) * xi(7) * xi(16) * s_5 &
         + (-4.330127018922197e-01_rk) * xi(6) * xi(17) * c_5 &
         + (+7.499999999999997e-01_rk) * xi(6) * xi(17) * s_5 &
         + (+2.500000000000003e-01_rk) * xi(6) * xi(16) * c_5 &
         + (-4.330127018922191e-01_rk) * xi(6) * xi(16) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(16) * c_5 &
         + (+4.330127018922198e-01_rk) * xi(4) * xi(15) * c_5 &
         + (-7.499999999999997e-01_rk) * xi(4) * xi(15) * s_5 &
         + (+2.500000000000003e-01_rk) * xi(4) * xi(14) * c_5 &
         + (-4.330127018922191e-01_rk) * xi(4) * xi(14) * s_5 &
         + (-4.330127018922184e-01_rk) * xi(3) * xi(15) * c_5 &
         + (-7.500000000000002e-01_rk) * xi(3) * xi(15) * s_5 &
         + (+2.499999999999994e-01_rk) * xi(3) * xi(14) * c_5 &
         + (+4.330127018922195e-01_rk) * xi(3) * xi(14) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(14) * c_5
    phi_y(52) = (-7.500000000000003e-01_rk) * xi(7) * xi(17) * c_5 &
         + (+4.330127018922190e-01_rk) * xi(7) * xi(17) * s_5 &
         + (-4.330127018922194e-01_rk) * xi(7) * xi(16) * c_5 &
         + (+2.499999999999998e-01_rk) * xi(7) * xi(16) * s_5 &
         + (-7.499999999999997e-01_rk) * xi(6) * xi(17) * c_5 &
         + (-4.330127018922193e-01_rk) * xi(6) * xi(17) * s_5 &
         + (+4.330127018922191e-01_rk) * xi(6) * xi(16) * c_5 &
         + (+2.500000000000000e-01_rk) * xi(6) * xi(16) * s_5 &
         + (+9.999999999999997e-01_rk) * xi(5) * xi(16) * s_5 &
         + (-7.499999999999997e-01_rk) * xi(4) * xi(15) * c_5 &
         + (-4.330127018922193e-01_rk) * xi(4) * xi(15) * s_5 &
         + (-4.330127018922191e-01_rk) * xi(4) * xi(14) * c_5 &
         + (-2.500000000000000e-01_rk) * xi(4) * xi(14) * s_5 &
         + (-7.500000000000003e-01_rk) * xi(3) * xi(15) * c_5 &
         + (+4.330127018922190e-01_rk) * xi(3) * xi(15) * s_5 &
         + (+4.330127018922194e-01_rk) * xi(3) * xi(14) * c_5 &
         + (-2.499999999999998e-01_rk) * xi(3) * xi(14) * s_5 &
         + (-9.999999999999997e-01_rk) * xi(2) * xi(14) * s_5
    phi_x(53) = (+4.330127018922194e-01_rk) * xi(7) * xi(17) * c_5 &
         + (-2.499999999999998e-01_rk) * xi(7) * xi(17) * s_5 &
         + (-7.500000000000003e-01_rk) * xi(7) * xi(16) * c_5 &
         + (+4.330127018922190e-01_rk) * xi(7) * xi(16) * s_5 &
         + (-4.330127018922191e-01_rk) * xi(6) * xi(17) * c_5 &
         + (-2.499999999999999e-01_rk) * xi(6) * xi(17) * s_5 &
         + (-7.499999999999997e-01_rk) * xi(6) * xi(16) * c_5 &
         + (-4.330127018922193e-01_rk) * xi(6) * xi(16) * s_5 &
         + (-9.999999999999997e-01_rk) * xi(5) * xi(17) * s_5 &
         + (+4.330127018922191e-01_rk) * xi(4) * xi(15) * c_5 &
         + (+2.499999999999999e-01_rk) * xi(4) * xi(15) * s_5 &
         + (-7.499999999999997e-01_rk) * xi(4) * xi(14) * c_5 &
         + (-4.330127018922193e-01_rk) * xi(4) * xi(14) * s_5 &
         + (-4.330127018922194e-01_rk) * xi(3) * xi(15) * c_5 &
         + (+2.499999999999998e-01_rk) * xi(3) * xi(15) * s_5 &
         + (-7.500000000000003e-01_rk) * xi(3) * xi(14) * c_5 &
         + (+4.330127018922190e-01_rk) * xi(3) * xi(14) * s_5 &
         + (+9.999999999999997e-01_rk) * xi(2) * xi(15) * s_5
    phi_y(53) = (+2.499999999999997e-01_rk) * xi(7) * xi(17) * c_5 &
         + (+4.330127018922193e-01_rk) * xi(7) * xi(17) * s_5 &
         + (-4.330127018922189e-01_rk) * xi(7) * xi(16) * c_5 &
         + (-7.499999999999998e-01_rk) * xi(7) * xi(16) * s_5 &
         + (+2.500000000000003e-01_rk) * xi(6) * xi(17) * c_5 &
         + (-4.330127018922191e-01_rk) * xi(6) * xi(17) * s_5 &
         + (+4.330127018922197e-01_rk) * xi(6) * xi(16) * c_5 &
         + (-7.499999999999997e-01_rk) * xi(6) * xi(16) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(17) * c_5 &
         + (+2.500000000000003e-01_rk) * xi(4) * xi(15) * c_5 &
         + (-4.330127018922191e-01_rk) * xi(4) * xi(15) * s_5 &
         + (-4.330127018922198e-01_rk) * xi(4) * xi(14) * c_5 &
         + (+7.499999999999997e-01_rk) * xi(4) * xi(14) * s_5 &
         + (+2.499999999999994e-01_rk) * xi(3) * xi(15) * c_5 &
         + (+4.330127018922195e-01_rk) * xi(3) * xi(15) * s_5 &
         + (+4.330127018922184e-01_rk) * xi(3) * xi(14) * c_5 &
         + (+7.500000000000002e-01_rk) * xi(3) * xi(14) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(15) * c_5
    phi_x(54) = (+9.999999999999993e-01_rk) * xi(7) * xi(15) * s_5 &
         + (-9.999999999999992e-01_rk) * xi(7) * xi(14) * c_5 &
         + (+9.999999999999987e-01_rk) * xi(6) * xi(15) * s_5 &
         + (-9.999999999999993e-01_rk) * xi(6) * xi(14) * c_5 &
         + (+9.999999999999989e-01_rk) * xi(5) * xi(15) * s_5 &
         + (-9.999999999999992e-01_rk) * xi(5) * xi(14) * c_5 &
         + (-9.999999999999987e-01_rk) * xi(4) * xi(17) * s_5 &
         + (-9.999999999999993e-01_rk) * xi(4) * xi(16) * c_5 &
         + (-1.000000000000000e+00_rk) * xi(3) * xi(17) * s_5 &
         + (-9.999999999999986e-01_rk) * xi(3) * xi(16) * c_5 &
         + (-9.999999999999989e-01_rk) * xi(2) * xi(17) * s_5 &
         + (-9.999999999999992e-01_rk) * xi(2) * xi(16) * c_5
    phi_y(54) = (+9.999999999999990e-01_rk) * xi(7) * xi(15) * c_5 &
         + (+9.999999999999983e-01_rk) * xi(7) * xi(14) * s_5 &
         + (+9.999999999999990e-01_rk) * xi(6) * xi(15) * c_5 &
         + (+9.999999999999993e-01_rk) * xi(6) * xi(14) * s_5 &
         + (+9.999999999999993e-01_rk) * xi(5) * xi(15) * c_5 &
         + (+9.999999999999989e-01_rk) * xi(5) * xi(14) * s_5 &
         + (+9.999999999999990e-01_rk) * xi(4) * xi(17) * c_5 &
         + (-9.999999999999993e-01_rk) * xi(4) * xi(16) * s_5 &
         + (+9.999999999999990e-01_rk) * xi(3) * xi(17) * c_5 &
         + (-9.999999999999983e-01_rk) * xi(3) * xi(16) * s_5 &
         + (+9.999999999999993e-01_rk) * xi(2) * xi(17) * c_5 &
         + (-9.999999999999989e-01_rk) * xi(2) * xi(16) * s_5
    phi_x(55) = (+4.330127018922189e-01_rk) * xi(13) * xi(17) * c_5 &
         + (+7.499999999999998e-01_rk) * xi(13) * xi(17) * s_5 &
         + (+2.499999999999997e-01_rk) * xi(13) * xi(16) * c_5 &
         + (+4.330127018922193e-01_rk) * xi(13) * xi(16) * s_5 &
         + (-4.330127018922197e-01_rk) * xi(12) * xi(17) * c_5 &
         + (+7.499999999999997e-01_rk) * xi(12) * xi(17) * s_5 &
         + (+2.500000000000003e-01_rk) * xi(12) * xi(16) * c_5 &
         + (-4.330127018922191e-01_rk) * xi(12) * xi(16) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(11) * xi(16) * c_5 &
         + (+4.330127018922198e-01_rk) * xi(10) * xi(15) * c_5 &
         + (-7.499999999999997e-01_rk) * xi(10) * xi(15) * s_5 &
         + (+2.500000000000003e-01_rk) * xi(10) * xi(14) * c_5 &
         + (-4.330127018922191e-01_rk) * xi(10) * xi(14) * s_5 &
         + (-4.330127018922184e-01_rk) * xi(9) * xi(15) * c_5 &
         + (-7.500000000000002e-01_rk) * xi(9) * xi(15) * s_5 &
         + (+2.499999999999994e-01_rk) * xi(9) * xi(14) * c_5 &
         + (+4.330127018922195e-01_rk) * xi(9) * xi(14) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(14) * c_5
    phi_y(55) = (-7.500000000000003e-01_rk) * xi(13) * xi(17) * c_5 &
         + (+4.330127018922190e-01_rk) * xi(13) * xi(17) * s_5 &
         + (-4.330127018922194e-01_rk) * xi(13) * xi(16) * c_5 &
         + (+2.499999999999998e-01_rk) * xi(13) * xi(16) * s_5 &
         + (-7.499999999999997e-01_rk) * xi(12) * xi(17) * c_5 &
         + (-4.330127018922193e-01_rk) * xi(12) * xi(17) * s_5 &
         + (+4.330127018922191e-01_rk) * xi(12) * xi(16) * c_5 &
         + (+2.500000000000000e-01_rk) * xi(12) * xi(16) * s_5 &
         + (+9.999999999999997e-01_rk) * xi(11) * xi(16) * s_5 &
         + (-7.499999999999997e-01_rk) * xi(10) * xi(15) * c_5 &
         + (-4.330127018922193e-01_rk) * xi(10) * xi(15) * s_5 &
         + (-4.330127018922191e-01_rk) * xi(10) * xi(14) * c_5 &
         + (-2.500000000000000e-01_rk) * xi(10) * xi(14) * s_5 &
         + (-7.500000000000003e-01_rk) * xi(9) * xi(15) * c_5 &
         + (+4.330127018922190e-01_rk) * xi(9) * xi(15) * s_5 &
         + (+4.330127018922194e-01_rk) * xi(9) * xi(14) * c_5 &
         + (-2.499999999999998e-01_rk) * xi(9) * xi(14) * s_5 &
         + (-9.999999999999997e-01_rk) * xi(8) * xi(14) * s_5
    phi_x(56) = (+4.330127018922194e-01_rk) * xi(13) * xi(17) * c_5 &
         + (-2.499999999999998e-01_rk) * xi(13) * xi(17) * s_5 &
         + (-7.500000000000003e-01_rk) * xi(13) * xi(16) * c_5 &
         + (+4.330127018922190e-01_rk) * xi(13) * xi(16) * s_5 &
         + (-4.330127018922191e-01_rk) * xi(12) * xi(17) * c_5 &
         + (-2.499999999999999e-01_rk) * xi(12) * xi(17) * s_5 &
         + (-7.499999999999997e-01_rk) * xi(12) * xi(16) * c_5 &
         + (-4.330127018922193e-01_rk) * xi(12) * xi(16) * s_5 &
         + (-9.999999999999997e-01_rk) * xi(11) * xi(17) * s_5 &
         + (+4.330127018922191e-01_rk) * xi(10) * xi(15) * c_5 &
         + (+2.499999999999999e-01_rk) * xi(10) * xi(15) * s_5 &
         + (-7.499999999999997e-01_rk) * xi(10) * xi(14) * c_5 &
         + (-4.330127018922193e-01_rk) * xi(10) * xi(14) * s_5 &
         + (-4.330127018922194e-01_rk) * xi(9) * xi(15) * c_5 &
         + (+2.499999999999998e-01_rk) * xi(9) * xi(15) * s_5 &
         + (-7.500000000000003e-01_rk) * xi(9) * xi(14) * c_5 &
         + (+4.330127018922190e-01_rk) * xi(9) * xi(14) * s_5 &
         + (+9.999999999999997e-01_rk) * xi(8) * xi(15) * s_5
    phi_y(56) = (+2.499999999999997e-01_rk) * xi(13) * xi(17) * c_5 &
         + (+4.330127018922193e-01_rk) * xi(13) * xi(17) * s_5 &
         + (-4.330127018922189e-01_rk) * xi(13) * xi(16) * c_5 &
         + (-7.499999999999998e-01_rk) * xi(13) * xi(16) * s_5 &
         + (+2.500000000000003e-01_rk) * xi(12) * xi(17) * c_5 &
         + (-4.330127018922191e-01_rk) * xi(12) * xi(17) * s_5 &
         + (+4.330127018922197e-01_rk) * xi(12) * xi(16) * c_5 &
         + (-7.499999999999997e-01_rk) * xi(12) * xi(16) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(11) * xi(17) * c_5 &
         + (+2.500000000000003e-01_rk) * xi(10) * xi(15) * c_5 &
         + (-4.330127018922191e-01_rk) * xi(10) * xi(15) * s_5 &
         + (-4.330127018922198e-01_rk) * xi(10) * xi(14) * c_5 &
         + (+7.499999999999997e-01_rk) * xi(10) * xi(14) * s_5 &
         + (+2.499999999999994e-01_rk) * xi(9) * xi(15) * c_5 &
         + (+4.330127018922195e-01_rk) * xi(9) * xi(15) * s_5 &
         + (+4.330127018922184e-01_rk) * xi(9) * xi(14) * c_5 &
         + (+7.500000000000002e-01_rk) * xi(9) * xi(14) * s_5 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(15) * c_5
    phi_x(57) = (+9.999999999999993e-01_rk) * xi(13) * xi(15) * s_5 &
         + (-9.999999999999992e-01_rk) * xi(13) * xi(14) * c_5 &
         + (+9.999999999999987e-01_rk) * xi(12) * xi(15) * s_5 &
         + (-9.999999999999993e-01_rk) * xi(12) * xi(14) * c_5 &
         + (+9.999999999999989e-01_rk) * xi(11) * xi(15) * s_5 &
         + (-9.999999999999992e-01_rk) * xi(11) * xi(14) * c_5 &
         + (-9.999999999999987e-01_rk) * xi(10) * xi(17) * s_5 &
         + (-9.999999999999993e-01_rk) * xi(10) * xi(16) * c_5 &
         + (-1.000000000000000e+00_rk) * xi(9) * xi(17) * s_5 &
         + (-9.999999999999986e-01_rk) * xi(9) * xi(16) * c_5 &
         + (-9.999999999999989e-01_rk) * xi(8) * xi(17) * s_5 &
         + (-9.999999999999992e-01_rk) * xi(8) * xi(16) * c_5
    phi_y(57) = (+9.999999999999990e-01_rk) * xi(13) * xi(15) * c_5 &
         + (+9.999999999999983e-01_rk) * xi(13) * xi(14) * s_5 &
         + (+9.999999999999990e-01_rk) * xi(12) * xi(15) * c_5 &
         + (+9.999999999999993e-01_rk) * xi(12) * xi(14) * s_5 &
         + (+9.999999999999993e-01_rk) * xi(11) * xi(15) * c_5 &
         + (+9.999999999999989e-01_rk) * xi(11) * xi(14) * s_5 &
         + (+9.999999999999990e-01_rk) * xi(10) * xi(17) * c_5 &
         + (-9.999999999999993e-01_rk) * xi(10) * xi(16) * s_5 &
         + (+9.999999999999990e-01_rk) * xi(9) * xi(17) * c_5 &
         + (-9.999999999999983e-01_rk) * xi(9) * xi(16) * s_5 &
         + (+9.999999999999993e-01_rk) * xi(8) * xi(17) * c_5 &
         + (-9.999999999999989e-01_rk) * xi(8) * xi(16) * s_5
    phi_x(58) = (-5.000000000000000e-01_rk) * xi(7)**3 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(7)**3 * s_1 &
         + (-4.999999999999997e-01_rk) * xi(6)**3 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(6)**3 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(5)**3 * c_1 &
         + (-4.999999999999997e-01_rk) * xi(4)**3 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(4)**3 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(3)**3 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(3)**3 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**3 * c_1
    phi_y(58) = (+8.660254037844385e-01_rk) * xi(7)**3 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(7)**3 * s_1 &
         + (-8.660254037844386e-01_rk) * xi(6)**3 * c_1 &
         + (+4.999999999999999e-01_rk) * xi(6)**3 * s_1 &
         + (-1.000000000000000e+00_rk) * xi(5)**3 * s_1 &
         + (+8.660254037844386e-01_rk) * xi(4)**3 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(4)**3 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3)**3 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(3)**3 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**3 * s_1
    phi_x(59) = (-5.000000000000000e-01_rk) * xi(13)**3 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(13)**3 * s_1 &
         + (-4.999999999999997e-01_rk) * xi(12)**3 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(12)**3 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(11)**3 * c_1 &
         + (-4.999999999999997e-01_rk) * xi(10)**3 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(10)**3 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(9)**3 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(9)**3 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8)**3 * c_1
    phi_y(59) = (+8.660254037844385e-01_rk) * xi(13)**3 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(13)**3 * s_1 &
         + (-8.660254037844386e-01_rk) * xi(12)**3 * c_1 &
         + (+4.999999999999999e-01_rk) * xi(12)**3 * s_1 &
         + (-1.000000000000000e+00_rk) * xi(11)**3 * s_1 &
         + (+8.660254037844386e-01_rk) * xi(10)**3 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(10)**3 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(9)**3 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(9)**3 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8)**3 * s_1
    phi_x(60) = (+9.999999999999994e-01_rk) * xi(17)**3 * s_1 &
         + (-9.999999999999997e-01_rk) * xi(16) * xi(17)**2 * c_1 &
         + (+9.999999999999996e-01_rk) * xi(16)**2 * xi(17) * s_1 &
         + (-9.999999999999997e-01_rk) * xi(16)**3 * c_1 &
         + (-9.999999999999994e-01_rk) * xi(15)**3 * s_1 &
         + (-1.000000000000000e+00_rk) * xi(14) * xi(15)**2 * c_1 &
         + (-9.999999999999996e-01_rk) * xi(14)**2 * xi(15) * s_1 &
         + (-9.999999999999997e-01_rk) * xi(14)**3 * c_1
    phi_y(60) = (+9.999999999999996e-01_rk) * xi(17)**3 * c_1 &
         + (+9.999999999999991e-01_rk) * xi(16) * xi(17)**2 * s_1 &
         + (+9.999999999999996e-01_rk) * xi(16)**2 * xi(17) * c_1 &
         + (+9.999999999999990e-01_rk) * xi(16)**3 * s_1 &
         + (+9.999999999999996e-01_rk) * xi(15)**3 * c_1 &
         + (-9.999999999999991e-01_rk) * xi(14) * xi(15)**2 * s_1 &
         + (+9.999999999999996e-01_rk) * xi(14)**2 * xi(15) * c_1 &
         + (-9.999999999999991e-01_rk) * xi(14)**3 * s_1
    phi_x(61) = (-5.000000000000000e-01_rk) * xi(1)**2 * xi(7) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(1)**2 * xi(7) * s_1 &
         + (-4.999999999999997e-01_rk) * xi(1)**2 * xi(6) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(1)**2 * xi(6) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(5) * c_1 &
         + (-4.999999999999997e-01_rk) * xi(1)**2 * xi(4) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(1)**2 * xi(4) * s_1 &
         + (-5.000000000000002e-01_rk) * xi(1)**2 * xi(3) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(1)**2 * xi(3) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(2) * c_1
    phi_y(61) = (+8.660254037844385e-01_rk) * xi(1)**2 * xi(7) * c_1 &
         + (+5.000000000000000e-01_rk) * xi(1)**2 * xi(7) * s_1 &
         + (-8.660254037844386e-01_rk) * xi(1)**2 * xi(6) * c_1 &
         + (+4.999999999999999e-01_rk) * xi(1)**2 * xi(6) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(1)**2 * xi(5) * s_1 &
         + (+8.660254037844386e-01_rk) * xi(1)**2 * xi(4) * c_1 &
         + (-4.999999999999999e-01_rk) * xi(1)**2 * xi(4) * s_1 &
         + (-8.660254037844385e-01_rk) * xi(1)**2 * xi(3) * c_1 &
         + (-5.000000000000000e-01_rk) * xi(1)**2 * xi(3) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(2) * s_1
    phi_x(62) = (-5.000000000000000e-01_rk) * xi(1) * xi(7)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(1) * xi(7)**2 * s_1 &
         + (-4.999999999999997e-01_rk) * xi(1) * xi(6)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(1) * xi(6)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(5)**2 * c_1 &
         + (-4.999999999999997e-01_rk) * xi(1) * xi(4)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(1) * xi(4)**2 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(1) * xi(3)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(1) * xi(3)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(2)**2 * c_1
    phi_y(62) = (+8.660254037844385e-01_rk) * xi(1) * xi(7)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(1) * xi(7)**2 * s_1 &
         + (-8.660254037844386e-01_rk) * xi(1) * xi(6)**2 * c_1 &
         + (+4.999999999999999e-01_rk) * xi(1) * xi(6)**2 * s_1 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(5)**2 * s_1 &
         + (+8.660254037844386e-01_rk) * xi(1) * xi(4)**2 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(1) * xi(4)**2 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(1) * xi(3)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(1) * xi(3)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(2)**2 * s_1
    phi_x(63) = (-5.000000000000000e-01_rk) * xi(1)**2 * xi(13) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(1)**2 * xi(13) * s_1 &
         + (-4.999999999999997e-01_rk) * xi(1)**2 * xi(12) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(1)**2 * xi(12) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(11) * c_1 &
         + (-4.999999999999997e-01_rk) * xi(1)**2 * xi(10) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(1)**2 * xi(10) * s_1 &
         + (-5.000000000000002e-01_rk) * xi(1)**2 * xi(9) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(1)**2 * xi(9) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(8) * c_1
    phi_y(63) = (+8.660254037844385e-01_rk) * xi(1)**2 * xi(13) * c_1 &
         + (+5.000000000000000e-01_rk) * xi(1)**2 * xi(13) * s_1 &
         + (-8.660254037844386e-01_rk) * xi(1)**2 * xi(12) * c_1 &
         + (+4.999999999999999e-01_rk) * xi(1)**2 * xi(12) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(1)**2 * xi(11) * s_1 &
         + (+8.660254037844386e-01_rk) * xi(1)**2 * xi(10) * c_1 &
         + (-4.999999999999999e-01_rk) * xi(1)**2 * xi(10) * s_1 &
         + (-8.660254037844385e-01_rk) * xi(1)**2 * xi(9) * c_1 &
         + (-5.000000000000000e-01_rk) * xi(1)**2 * xi(9) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(8) * s_1
    phi_x(64) = (-5.000000000000000e-01_rk) * xi(1) * xi(13)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(1) * xi(13)**2 * s_1 &
         + (-4.999999999999997e-01_rk) * xi(1) * xi(12)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(1) * xi(12)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(11)**2 * c_1 &
         + (-4.999999999999997e-01_rk) * xi(1) * xi(10)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(1) * xi(10)**2 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(1) * xi(9)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(1) * xi(9)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(8)**2 * c_1
    phi_y(64) = (+8.660254037844385e-01_rk) * xi(1) * xi(13)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(1) * xi(13)**2 * s_1 &
         + (-8.660254037844386e-01_rk) * xi(1) * xi(12)**2 * c_1 &
         + (+4.999999999999999e-01_rk) * xi(1) * xi(12)**2 * s_1 &
         + (-1.000000000000000e+00_rk) * xi(1) * xi(11)**2 * s_1 &
         + (+8.660254037844386e-01_rk) * xi(1) * xi(10)**2 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(1) * xi(10)**2 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(1) * xi(9)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(1) * xi(9)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(8)**2 * s_1
    phi_x(65) = (+9.999999999999996e-01_rk) * xi(1)**2 * xi(17) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(1)**2 * xi(16) * c_1 &
         + (-9.999999999999996e-01_rk) * xi(1)**2 * xi(15) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(1)**2 * xi(14) * c_1
    phi_y(65) = (+9.999999999999997e-01_rk) * xi(1)**2 * xi(17) * c_1 &
         + (+9.999999999999996e-01_rk) * xi(1)**2 * xi(16) * s_1 &
         + (+9.999999999999997e-01_rk) * xi(1)**2 * xi(15) * c_1 &
         + (-9.999999999999996e-01_rk) * xi(1)**2 * xi(14) * s_1
    phi_x(66) = (-5.000000000000000e-01_rk) * xi(1) * xi(17)**2 * c_1 &
         + (+9.999999999999998e-01_rk) * xi(1) * xi(16) * xi(17) * s_1 &
         + (+4.999999999999999e-01_rk) * xi(1) * xi(16)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(1) * xi(15)**2 * c_1 &
         + (-9.999999999999998e-01_rk) * xi(1) * xi(14) * xi(15) * s_1 &
         + (+4.999999999999999e-01_rk) * xi(1) * xi(14)**2 * c_1
    phi_y(66) = (+4.999999999999998e-01_rk) * xi(1) * xi(17)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(16) * xi(17) * c_1 &
         + (-4.999999999999998e-01_rk) * xi(1) * xi(16)**2 * s_1 &
         + (-4.999999999999998e-01_rk) * xi(1) * xi(15)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1) * xi(14) * xi(15) * c_1 &
         + (+4.999999999999998e-01_rk) * xi(1) * xi(14)**2 * s_1
    phi_x(67) = (-4.999999999999995e-01_rk) * xi(6) * xi(7)**2 * c_1 &
         + (+8.660254037844387e-01_rk) * xi(6) * xi(7)**2 * s_1 &
         + (-4.999999999999993e-01_rk) * xi(6)**2 * xi(7) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(6)**2 * xi(7) * s_1 &
         + (-5.000000000000006e-01_rk) * xi(5) * xi(7)**2 * c_1 &
         + (+8.660254037844382e-01_rk) * xi(5) * xi(7)**2 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(5) * xi(6)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(5) * xi(6)**2 * s_1 &
         + (+9.999999999999999e-01_rk) * xi(5)**2 * xi(7) * c_1 &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * xi(6) * c_1 &
         + (-4.999999999999993e-01_rk) * xi(3) * xi(4)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(3) * xi(4)**2 * s_1 &
         + (-4.999999999999997e-01_rk) * xi(3)**2 * xi(4) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(3)**2 * xi(4) * s_1 &
         + (-5.000000000000003e-01_rk) * xi(2) * xi(4)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(2) * xi(4)**2 * s_1 &
         + (-5.000000000000004e-01_rk) * xi(2) * xi(3)**2 * c_1 &
         + (+8.660254037844382e-01_rk) * xi(2) * xi(3)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**2 * xi(4) * c_1 &
         + (+9.999999999999999e-01_rk) * xi(2)**2 * xi(3) * c_1
    phi_y(67) = (+8.660254037844385e-01_rk) * xi(6) * xi(7)**2 * c_1 &
         + (+4.999999999999998e-01_rk) * xi(6) * xi(7)**2 * s_1 &
         + (-8.660254037844388e-01_rk) * xi(6)**2 * xi(7) * c_1 &
         + (+4.999999999999997e-01_rk) * xi(6)**2 * xi(7) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(5) * xi(7)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(5) * xi(7)**2 * s_1 &
         + (-8.660254037844382e-01_rk) * xi(5) * xi(6)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(5) * xi(6)**2 * s_1 &
         + (-9.999999999999999e-01_rk) * xi(5)**2 * xi(7) * s_1 &
         + (-9.999999999999997e-01_rk) * xi(5)**2 * xi(6) * s_1 &
         + (+8.660254037844388e-01_rk) * xi(3) * xi(4)**2 * c_1 &
         + (-4.999999999999997e-01_rk) * xi(3) * xi(4)**2 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3)**2 * xi(4) * c_1 &
         + (-4.999999999999998e-01_rk) * xi(3)**2 * xi(4) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(2) * xi(4)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(2) * xi(4)**2 * s_1 &
         + (-8.660254037844382e-01_rk) * xi(2) * xi(3)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(2) * xi(3)**2 * s_1 &
         + (+9.999999999999997e-01_rk) * xi(2)**2 * xi(4) * s_1 &
         + (+9.999999999999999e-01_rk) * xi(2)**2 * xi(3) * s_1
    phi_x(68) = (-8.660254037844385e-01_rk) * xi(6) * xi(7)**2 * c_1 &
         + (-4.999999999999998e-01_rk) * xi(6) * xi(7)**2 * s_1 &
         + (-8.660254037844388e-01_rk) * xi(6)**2 * xi(7) * c_1 &
         + (+4.999999999999997e-01_rk) * xi(6)**2 * xi(7) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(5) * xi(7)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(5) * xi(7)**2 * s_1 &
         + (+8.660254037844382e-01_rk) * xi(5) * xi(6)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(5) * xi(6)**2 * s_1 &
         + (+9.999999999999999e-01_rk) * xi(5)**2 * xi(7) * s_1 &
         + (-9.999999999999997e-01_rk) * xi(5)**2 * xi(6) * s_1 &
         + (-8.660254037844388e-01_rk) * xi(3) * xi(4)**2 * c_1 &
         + (+4.999999999999997e-01_rk) * xi(3) * xi(4)**2 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3)**2 * xi(4) * c_1 &
         + (-4.999999999999998e-01_rk) * xi(3)**2 * xi(4) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(2) * xi(4)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(2) * xi(4)**2 * s_1 &
         + (+8.660254037844382e-01_rk) * xi(2) * xi(3)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(2) * xi(3)**2 * s_1 &
         + (-9.999999999999997e-01_rk) * xi(2)**2 * xi(4) * s_1 &
         + (+9.999999999999999e-01_rk) * xi(2)**2 * xi(3) * s_1
    phi_y(68) = (-4.999999999999995e-01_rk) * xi(6) * xi(7)**2 * c_1 &
         + (+8.660254037844387e-01_rk) * xi(6) * xi(7)**2 * s_1 &
         + (+4.999999999999993e-01_rk) * xi(6)**2 * xi(7) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(6)**2 * xi(7) * s_1 &
         + (+5.000000000000006e-01_rk) * xi(5) * xi(7)**2 * c_1 &
         + (-8.660254037844382e-01_rk) * xi(5) * xi(7)**2 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(5) * xi(6)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(5) * xi(6)**2 * s_1 &
         + (+9.999999999999999e-01_rk) * xi(5)**2 * xi(7) * c_1 &
         + (-1.000000000000000e+00_rk) * xi(5)**2 * xi(6) * c_1 &
         + (-4.999999999999993e-01_rk) * xi(3) * xi(4)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(3) * xi(4)**2 * s_1 &
         + (+4.999999999999997e-01_rk) * xi(3)**2 * xi(4) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(3)**2 * xi(4) * s_1 &
         + (+5.000000000000003e-01_rk) * xi(2) * xi(4)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(2) * xi(4)**2 * s_1 &
         + (-5.000000000000004e-01_rk) * xi(2) * xi(3)**2 * c_1 &
         + (+8.660254037844382e-01_rk) * xi(2) * xi(3)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**2 * xi(4) * c_1 &
         + (-9.999999999999999e-01_rk) * xi(2)**2 * xi(3) * c_1
    phi_x(69) = (-5.000000000000000e-01_rk) * xi(4) * xi(7)**2 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(4) * xi(6)**2 * c_1 &
         + (-8.660254037844387e-01_rk) * xi(4) * xi(6)**2 * s_1 &
         + (+2.500000000000001e-01_rk) * xi(4) * xi(5)**2 * c_1 &
         + (-4.330127018922192e-01_rk) * xi(4) * xi(5)**2 * s_1 &
         + (-5.000000000000000e-01_rk) * xi(4)**2 * xi(7) * c_1 &
         + (-4.999999999999999e-01_rk) * xi(4)**2 * xi(6) * c_1 &
         + (-8.660254037844387e-01_rk) * xi(4)**2 * xi(6) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(4)**2 * xi(5) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(4)**2 * xi(5) * s_1 &
         + (-5.000000000000001e-01_rk) * xi(3) * xi(7)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(3) * xi(7)**2 * s_1 &
         + (-5.000000000000000e-01_rk) * xi(3) * xi(6)**2 * c_1 &
         + (+2.499999999999999e-01_rk) * xi(3) * xi(5)**2 * c_1 &
         + (+4.330127018922194e-01_rk) * xi(3) * xi(5)**2 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(3)**2 * xi(7) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(3)**2 * xi(7) * s_1 &
         + (-5.000000000000000e-01_rk) * xi(3)**2 * xi(6) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(3)**2 * xi(5) * c_1 &
         + (+4.330127018922194e-01_rk) * xi(3)**2 * xi(5) * s_1 &
         + (+2.499999999999999e-01_rk) * xi(2) * xi(7)**2 * c_1 &
         + (+4.330127018922194e-01_rk) * xi(2) * xi(7)**2 * s_1 &
         + (+2.500000000000001e-01_rk) * xi(2) * xi(6)**2 * c_1 &
         + (-4.330127018922192e-01_rk) * xi(2) * xi(6)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(5)**2 * c_1 &
         + (+2.499999999999999e-01_rk) * xi(2)**2 * xi(7) * c_1 &
         + (+4.330127018922194e-01_rk) * xi(2)**2 * xi(7) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(2)**2 * xi(6) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(2)**2 * xi(6) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**2 * xi(5) * c_1
    phi_y(69) = (+8.660254037844386e-01_rk) * xi(4) * xi(7)**2 * c_1 &
         + (+4.330127018922194e-01_rk) * xi(4) * xi(5)**2 * c_1 &
         + (-7.499999999999999e-01_rk) * xi(4) * xi(5)**2 * s_1 &
         + (+8.660254037844386e-01_rk) * xi(4)**2 * xi(7) * c_1 &
         + (+4.330127018922194e-01_rk) * xi(4)**2 * xi(5) * c_1 &
         + (-7.499999999999999e-01_rk) * xi(4)**2 * xi(5) * s_1 &
         + (-8.660254037844386e-01_rk) * xi(3) * xi(6)**2 * c_1 &
         + (-4.330127018922192e-01_rk) * xi(3) * xi(5)**2 * c_1 &
         + (-7.500000000000000e-01_rk) * xi(3) * xi(5)**2 * s_1 &
         + (-8.660254037844386e-01_rk) * xi(3)**2 * xi(6) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(3)**2 * xi(5) * c_1 &
         + (-7.500000000000000e-01_rk) * xi(3)**2 * xi(5) * s_1 &
         + (+4.330127018922192e-01_rk) * xi(2) * xi(7)**2 * c_1 &
         + (+7.500000000000000e-01_rk) * xi(2) * xi(7)**2 * s_1 &
         + (-4.330127018922194e-01_rk) * xi(2) * xi(6)**2 * c_1 &
         + (+7.499999999999999e-01_rk) * xi(2) * xi(6)**2 * s_1 &
         + (+4.330127018922192e-01_rk) * xi(2)**2 * xi(7) * c_1 &
         + (+7.500000000000000e-01_rk) * xi(2)**2 * xi(7) * s_1 &
         + (-4.330127018922194e-01_rk) * xi(2)**2 * xi(6) * c_1 &
         + (+7.499999999999999e-01_rk) * xi(2)**2 * xi(6) * s_1
    phi_x(70) = (-8.660254037844386e-01_rk) * xi(4) * xi(7)**2 * s_1 &
         + (-7.499999999999999e-01_rk) * xi(4) * xi(5)**2 * c_1 &
         + (-4.330127018922194e-01_rk) * xi(4) * xi(5)**2 * s_1 &
         + (+8.660254037844386e-01_rk) * xi(4)**2 * xi(7) * s_1 &
         + (+7.499999999999999e-01_rk) * xi(4)**2 * xi(5) * c_1 &
         + (+4.330127018922194e-01_rk) * xi(4)**2 * xi(5) * s_1 &
         + (+8.660254037844386e-01_rk) * xi(3) * xi(6)**2 * s_1 &
         + (-7.500000000000000e-01_rk) * xi(3) * xi(5)**2 * c_1 &
         + (+4.330127018922192e-01_rk) * xi(3) * xi(5)**2 * s_1 &
         + (-8.660254037844386e-01_rk) * xi(3)**2 * xi(6) * s_1 &
         + (+7.500000000000000e-01_rk) * xi(3)**2 * xi(5) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(3)**2 * xi(5) * s_1 &
         + (+7.500000000000000e-01_rk) * xi(2) * xi(7)**2 * c_1 &
         + (-4.330127018922192e-01_rk) * xi(2) * xi(7)**2 * s_1 &
         + (+7.499999999999999e-01_rk) * xi(2) * xi(6)**2 * c_1 &
         + (+4.330127018922194e-01_rk) * xi(2) * xi(6)**2 * s_1 &
         + (-7.500000000000000e-01_rk) * xi(2)**2 * xi(7) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(2)**2 * xi(7) * s_1 &
         + (-7.499999999999999e-01_rk) * xi(2)**2 * xi(6) * c_1 &
         + (-4.330127018922194e-01_rk) * xi(2)**2 * xi(6) * s_1
    phi_y(70) = (-5.000000000000000e-01_rk) * xi(4) * xi(7)**2 * s_1 &
         + (+8.660254037844387e-01_rk) * xi(4) * xi(6)**2 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(4) * xi(6)**2 * s_1 &
         + (+4.330127018922192e-01_rk) * xi(4) * xi(5)**2 * c_1 &
         + (+2.500000000000001e-01_rk) * xi(4) * xi(5)**2 * s_1 &
         + (+5.000000000000000e-01_rk) * xi(4)**2 * xi(7) * s_1 &
         + (-8.660254037844387e-01_rk) * xi(4)**2 * xi(6) * c_1 &
         + (+4.999999999999999e-01_rk) * xi(4)**2 * xi(6) * s_1 &
         + (-4.330127018922192e-01_rk) * xi(4)**2 * xi(5) * c_1 &
         + (-2.500000000000001e-01_rk) * xi(4)**2 * xi(5) * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3) * xi(7)**2 * c_1 &
         + (-5.000000000000001e-01_rk) * xi(3) * xi(7)**2 * s_1 &
         + (-5.000000000000000e-01_rk) * xi(3) * xi(6)**2 * s_1 &
         + (-4.330127018922194e-01_rk) * xi(3) * xi(5)**2 * c_1 &
         + (+2.499999999999999e-01_rk) * xi(3) * xi(5)**2 * s_1 &
         + (+8.660254037844385e-01_rk) * xi(3)**2 * xi(7) * c_1 &
         + (+5.000000000000002e-01_rk) * xi(3)**2 * xi(7) * s_1 &
         + (+5.000000000000000e-01_rk) * xi(3)**2 * xi(6) * s_1 &
         + (+4.330127018922194e-01_rk) * xi(3)**2 * xi(5) * c_1 &
         + (-2.499999999999999e-01_rk) * xi(3)**2 * xi(5) * s_1 &
         + (-4.330127018922194e-01_rk) * xi(2) * xi(7)**2 * c_1 &
         + (+2.499999999999999e-01_rk) * xi(2) * xi(7)**2 * s_1 &
         + (+4.330127018922192e-01_rk) * xi(2) * xi(6)**2 * c_1 &
         + (+2.500000000000001e-01_rk) * xi(2) * xi(6)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(5)**2 * s_1 &
         + (+4.330127018922194e-01_rk) * xi(2)**2 * xi(7) * c_1 &
         + (-2.499999999999999e-01_rk) * xi(2)**2 * xi(7) * s_1 &
         + (-4.330127018922192e-01_rk) * xi(2)**2 * xi(6) * c_1 &
         + (-2.500000000000001e-01_rk) * xi(2)**2 * xi(6) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(2)**2 * xi(5) * s_1
    phi_x(71) = (-4.999999999999995e-01_rk) * xi(12) * xi(13)**2 * c_1 &
         + (+8.660254037844387e-01_rk) * xi(12) * xi(13)**2 * s_1 &
         + (-4.999999999999993e-01_rk) * xi(12)**2 * xi(13) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(12)**2 * xi(13) * s_1 &
         + (-5.000000000000006e-01_rk) * xi(11) * xi(13)**2 * c_1 &
         + (+8.660254037844382e-01_rk) * xi(11) * xi(13)**2 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(11) * xi(12)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(11) * xi(12)**2 * s_1 &
         + (+9.999999999999999e-01_rk) * xi(11)**2 * xi(13) * c_1 &
         + (+1.000000000000000e+00_rk) * xi(11)**2 * xi(12) * c_1 &
         + (-4.999999999999993e-01_rk) * xi(9) * xi(10)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(9) * xi(10)**2 * s_1 &
         + (-4.999999999999997e-01_rk) * xi(9)**2 * xi(10) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(9)**2 * xi(10) * s_1 &
         + (-5.000000000000003e-01_rk) * xi(8) * xi(10)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(8) * xi(10)**2 * s_1 &
         + (-5.000000000000004e-01_rk) * xi(8) * xi(9)**2 * c_1 &
         + (+8.660254037844382e-01_rk) * xi(8) * xi(9)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8)**2 * xi(10) * c_1 &
         + (+9.999999999999999e-01_rk) * xi(8)**2 * xi(9) * c_1
    phi_y(71) = (+8.660254037844385e-01_rk) * xi(12) * xi(13)**2 * c_1 &
         + (+4.999999999999998e-01_rk) * xi(12) * xi(13)**2 * s_1 &
         + (-8.660254037844388e-01_rk) * xi(12)**2 * xi(13) * c_1 &
         + (+4.999999999999997e-01_rk) * xi(12)**2 * xi(13) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(11) * xi(13)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(11) * xi(13)**2 * s_1 &
         + (-8.660254037844382e-01_rk) * xi(11) * xi(12)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(11) * xi(12)**2 * s_1 &
         + (-9.999999999999999e-01_rk) * xi(11)**2 * xi(13) * s_1 &
         + (-9.999999999999997e-01_rk) * xi(11)**2 * xi(12) * s_1 &
         + (+8.660254037844388e-01_rk) * xi(9) * xi(10)**2 * c_1 &
         + (-4.999999999999997e-01_rk) * xi(9) * xi(10)**2 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(9)**2 * xi(10) * c_1 &
         + (-4.999999999999998e-01_rk) * xi(9)**2 * xi(10) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(8) * xi(10)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(8) * xi(10)**2 * s_1 &
         + (-8.660254037844382e-01_rk) * xi(8) * xi(9)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(8) * xi(9)**2 * s_1 &
         + (+9.999999999999997e-01_rk) * xi(8)**2 * xi(10) * s_1 &
         + (+9.999999999999999e-01_rk) * xi(8)**2 * xi(9) * s_1
    phi_x(72) = (-8.660254037844385e-01_rk) * xi(12) * xi(13)**2 * c_1 &
         + (-4.999999999999998e-01_rk) * xi(12) * xi(13)**2 * s_1 &
         + (-8.660254037844388e-01_rk) * xi(12)**2 * xi(13) * c_1 &
         + (+4.999999999999997e-01_rk) * xi(12)**2 * xi(13) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(11) * xi(13)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(11) * xi(13)**2 * s_1 &
         + (+8.660254037844382e-01_rk) * xi(11) * xi(12)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(11) * xi(12)**2 * s_1 &
         + (+9.999999999999999e-01_rk) * xi(11)**2 * xi(13) * s_1 &
         + (-9.999999999999997e-01_rk) * xi(11)**2 * xi(12) * s_1 &
         + (-8.660254037844388e-01_rk) * xi(9) * xi(10)**2 * c_1 &
         + (+4.999999999999997e-01_rk) * xi(9) * xi(10)**2 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(9)**2 * xi(10) * c_1 &
         + (-4.999999999999998e-01_rk) * xi(9)**2 * xi(10) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(8) * xi(10)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(8) * xi(10)**2 * s_1 &
         + (+8.660254037844382e-01_rk) * xi(8) * xi(9)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(8) * xi(9)**2 * s_1 &
         + (-9.999999999999997e-01_rk) * xi(8)**2 * xi(10) * s_1 &
         + (+9.999999999999999e-01_rk) * xi(8)**2 * xi(9) * s_1
    phi_y(72) = (-4.999999999999995e-01_rk) * xi(12) * xi(13)**2 * c_1 &
         + (+8.660254037844387e-01_rk) * xi(12) * xi(13)**2 * s_1 &
         + (+4.999999999999993e-01_rk) * xi(12)**2 * xi(13) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(12)**2 * xi(13) * s_1 &
         + (+5.000000000000006e-01_rk) * xi(11) * xi(13)**2 * c_1 &
         + (-8.660254037844382e-01_rk) * xi(11) * xi(13)**2 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(11) * xi(12)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(11) * xi(12)**2 * s_1 &
         + (+9.999999999999999e-01_rk) * xi(11)**2 * xi(13) * c_1 &
         + (-1.000000000000000e+00_rk) * xi(11)**2 * xi(12) * c_1 &
         + (-4.999999999999993e-01_rk) * xi(9) * xi(10)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(9) * xi(10)**2 * s_1 &
         + (+4.999999999999997e-01_rk) * xi(9)**2 * xi(10) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(9)**2 * xi(10) * s_1 &
         + (+5.000000000000003e-01_rk) * xi(8) * xi(10)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(8) * xi(10)**2 * s_1 &
         + (-5.000000000000004e-01_rk) * xi(8) * xi(9)**2 * c_1 &
         + (+8.660254037844382e-01_rk) * xi(8) * xi(9)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8)**2 * xi(10) * c_1 &
         + (-9.999999999999999e-01_rk) * xi(8)**2 * xi(9) * c_1
    phi_x(73) = (-5.000000000000000e-01_rk) * xi(10) * xi(13)**2 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(10) * xi(12)**2 * c_1 &
         + (-8.660254037844387e-01_rk) * xi(10) * xi(12)**2 * s_1 &
         + (+2.500000000000001e-01_rk) * xi(10) * xi(11)**2 * c_1 &
         + (-4.330127018922192e-01_rk) * xi(10) * xi(11)**2 * s_1 &
         + (-5.000000000000000e-01_rk) * xi(10)**2 * xi(13) * c_1 &
         + (-4.999999999999999e-01_rk) * xi(10)**2 * xi(12) * c_1 &
         + (-8.660254037844387e-01_rk) * xi(10)**2 * xi(12) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(10)**2 * xi(11) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(10)**2 * xi(11) * s_1 &
         + (-5.000000000000001e-01_rk) * xi(9) * xi(13)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(9) * xi(13)**2 * s_1 &
         + (-5.000000000000000e-01_rk) * xi(9) * xi(12)**2 * c_1 &
         + (+2.499999999999999e-01_rk) * xi(9) * xi(11)**2 * c_1 &
         + (+4.330127018922194e-01_rk) * xi(9) * xi(11)**2 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(9)**2 * xi(13) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(9)**2 * xi(13) * s_1 &
         + (-5.000000000000000e-01_rk) * xi(9)**2 * xi(12) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(9)**2 * xi(11) * c_1 &
         + (+4.330127018922194e-01_rk) * xi(9)**2 * xi(11) * s_1 &
         + (+2.499999999999999e-01_rk) * xi(8) * xi(13)**2 * c_1 &
         + (+4.330127018922194e-01_rk) * xi(8) * xi(13)**2 * s_1 &
         + (+2.500000000000001e-01_rk) * xi(8) * xi(12)**2 * c_1 &
         + (-4.330127018922192e-01_rk) * xi(8) * xi(12)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(11)**2 * c_1 &
         + (+2.499999999999999e-01_rk) * xi(8)**2 * xi(13) * c_1 &
         + (+4.330127018922194e-01_rk) * xi(8)**2 * xi(13) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(8)**2 * xi(12) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(8)**2 * xi(12) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8)**2 * xi(11) * c_1
    phi_y(73) = (+8.660254037844386e-01_rk) * xi(10) * xi(13)**2 * c_1 &
         + (+4.330127018922194e-01_rk) * xi(10) * xi(11)**2 * c_1 &
         + (-7.499999999999999e-01_rk) * xi(10) * xi(11)**2 * s_1 &
         + (+8.660254037844386e-01_rk) * xi(10)**2 * xi(13) * c_1 &
         + (+4.330127018922194e-01_rk) * xi(10)**2 * xi(11) * c_1 &
         + (-7.499999999999999e-01_rk) * xi(10)**2 * xi(11) * s_1 &
         + (-8.660254037844386e-01_rk) * xi(9) * xi(12)**2 * c_1 &
         + (-4.330127018922192e-01_rk) * xi(9) * xi(11)**2 * c_1 &
         + (-7.500000000000000e-01_rk) * xi(9) * xi(11)**2 * s_1 &
         + (-8.660254037844386e-01_rk) * xi(9)**2 * xi(12) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(9)**2 * xi(11) * c_1 &
         + (-7.500000000000000e-01_rk) * xi(9)**2 * xi(11) * s_1 &
         + (+4.330127018922192e-01_rk) * xi(8) * xi(13)**2 * c_1 &
         + (+7.500000000000000e-01_rk) * xi(8) * xi(13)**2 * s_1 &
         + (-4.330127018922194e-01_rk) * xi(8) * xi(12)**2 * c_1 &
         + (+7.499999999999999e-01_rk) * xi(8) * xi(12)**2 * s_1 &
         + (+4.330127018922192e-01_rk) * xi(8)**2 * xi(13) * c_1 &
         + (+7.500000000000000e-01_rk) * xi(8)**2 * xi(13) * s_1 &
         + (-4.330127018922194e-01_rk) * xi(8)**2 * xi(12) * c_1 &
         + (+7.499999999999999e-01_rk) * xi(8)**2 * xi(12) * s_1
    phi_x(74) = (-8.660254037844386e-01_rk) * xi(10) * xi(13)**2 * s_1 &
         + (-7.499999999999999e-01_rk) * xi(10) * xi(11)**2 * c_1 &
         + (-4.330127018922194e-01_rk) * xi(10) * xi(11)**2 * s_1 &
         + (+8.660254037844386e-01_rk) * xi(10)**2 * xi(13) * s_1 &
         + (+7.499999999999999e-01_rk) * xi(10)**2 * xi(11) * c_1 &
         + (+4.330127018922194e-01_rk) * xi(10)**2 * xi(11) * s_1 &
         + (+8.660254037844386e-01_rk) * xi(9) * xi(12)**2 * s_1 &
         + (-7.500000000000000e-01_rk) * xi(9) * xi(11)**2 * c_1 &
         + (+4.330127018922192e-01_rk) * xi(9) * xi(11)**2 * s_1 &
         + (-8.660254037844386e-01_rk) * xi(9)**2 * xi(12) * s_1 &
         + (+7.500000000000000e-01_rk) * xi(9)**2 * xi(11) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(9)**2 * xi(11) * s_1 &
         + (+7.500000000000000e-01_rk) * xi(8) * xi(13)**2 * c_1 &
         + (-4.330127018922192e-01_rk) * xi(8) * xi(13)**2 * s_1 &
         + (+7.499999999999999e-01_rk) * xi(8) * xi(12)**2 * c_1 &
         + (+4.330127018922194e-01_rk) * xi(8) * xi(12)**2 * s_1 &
         + (-7.500000000000000e-01_rk) * xi(8)**2 * xi(13) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(8)**2 * xi(13) * s_1 &
         + (-7.499999999999999e-01_rk) * xi(8)**2 * xi(12) * c_1 &
         + (-4.330127018922194e-01_rk) * xi(8)**2 * xi(12) * s_1
    phi_y(74) = (-5.000000000000000e-01_rk) * xi(10) * xi(13)**2 * s_1 &
         + (+8.660254037844387e-01_rk) * xi(10) * xi(12)**2 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(10) * xi(12)**2 * s_1 &
         + (+4.330127018922192e-01_rk) * xi(10) * xi(11)**2 * c_1 &
         + (+2.500000000000001e-01_rk) * xi(10) * xi(11)**2 * s_1 &
         + (+5.000000000000000e-01_rk) * xi(10)**2 * xi(13) * s_1 &
         + (-8.660254037844387e-01_rk) * xi(10)**2 * xi(12) * c_1 &
         + (+4.999999999999999e-01_rk) * xi(10)**2 * xi(12) * s_1 &
         + (-4.330127018922192e-01_rk) * xi(10)**2 * xi(11) * c_1 &
         + (-2.500000000000001e-01_rk) * xi(10)**2 * xi(11) * s_1 &
         + (-8.660254037844385e-01_rk) * xi(9) * xi(13)**2 * c_1 &
         + (-5.000000000000001e-01_rk) * xi(9) * xi(13)**2 * s_1 &
         + (-5.000000000000000e-01_rk) * xi(9) * xi(12)**2 * s_1 &
         + (-4.330127018922194e-01_rk) * xi(9) * xi(11)**2 * c_1 &
         + (+2.499999999999999e-01_rk) * xi(9) * xi(11)**2 * s_1 &
         + (+8.660254037844385e-01_rk) * xi(9)**2 * xi(13) * c_1 &
         + (+5.000000000000002e-01_rk) * xi(9)**2 * xi(13) * s_1 &
         + (+5.000000000000000e-01_rk) * xi(9)**2 * xi(12) * s_1 &
         + (+4.330127018922194e-01_rk) * xi(9)**2 * xi(11) * c_1 &
         + (-2.499999999999999e-01_rk) * xi(9)**2 * xi(11) * s_1 &
         + (-4.330127018922194e-01_rk) * xi(8) * xi(13)**2 * c_1 &
         + (+2.499999999999999e-01_rk) * xi(8) * xi(13)**2 * s_1 &
         + (+4.330127018922192e-01_rk) * xi(8) * xi(12)**2 * c_1 &
         + (+2.500000000000001e-01_rk) * xi(8) * xi(12)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(11)**2 * s_1 &
         + (+4.330127018922194e-01_rk) * xi(8)**2 * xi(13) * c_1 &
         + (-2.499999999999999e-01_rk) * xi(8)**2 * xi(13) * s_1 &
         + (-4.330127018922192e-01_rk) * xi(8)**2 * xi(12) * c_1 &
         + (-2.500000000000001e-01_rk) * xi(8)**2 * xi(12) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(8)**2 * xi(11) * s_1
    phi_x(75) = (-5.000000000000000e-01_rk) * xi(7)**2 * xi(13) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(7)**2 * xi(13) * s_1 &
         + (-4.999999999999997e-01_rk) * xi(6)**2 * xi(12) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(6)**2 * xi(12) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * xi(11) * c_1 &
         + (-4.999999999999997e-01_rk) * xi(4)**2 * xi(10) * c_1 &
         + (-8.660254037844386e-01_rk) * xi(4)**2 * xi(10) * s_1 &
         + (-5.000000000000002e-01_rk) * xi(3)**2 * xi(9) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(3)**2 * xi(9) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**2 * xi(8) * c_1
    phi_y(75) = (+8.660254037844385e-01_rk) * xi(7)**2 * xi(13) * c_1 &
         + (+5.000000000000000e-01_rk) * xi(7)**2 * xi(13) * s_1 &
         + (-8.660254037844386e-01_rk) * xi(6)**2 * xi(12) * c_1 &
         + (+4.999999999999999e-01_rk) * xi(6)**2 * xi(12) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(5)**2 * xi(11) * s_1 &
         + (+8.660254037844386e-01_rk) * xi(4)**2 * xi(10) * c_1 &
         + (-4.999999999999999e-01_rk) * xi(4)**2 * xi(10) * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3)**2 * xi(9) * c_1 &
         + (-5.000000000000000e-01_rk) * xi(3)**2 * xi(9) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**2 * xi(8) * s_1
    phi_x(76) = (-5.000000000000000e-01_rk) * xi(7) * xi(13)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(7) * xi(13)**2 * s_1 &
         + (-4.999999999999997e-01_rk) * xi(6) * xi(12)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(6) * xi(12)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(11)**2 * c_1 &
         + (-4.999999999999997e-01_rk) * xi(4) * xi(10)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(4) * xi(10)**2 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(3) * xi(9)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(3) * xi(9)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(8)**2 * c_1
    phi_y(76) = (+8.660254037844385e-01_rk) * xi(7) * xi(13)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(7) * xi(13)**2 * s_1 &
         + (-8.660254037844386e-01_rk) * xi(6) * xi(12)**2 * c_1 &
         + (+4.999999999999999e-01_rk) * xi(6) * xi(12)**2 * s_1 &
         + (-1.000000000000000e+00_rk) * xi(5) * xi(11)**2 * s_1 &
         + (+8.660254037844386e-01_rk) * xi(4) * xi(10)**2 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(4) * xi(10)**2 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3) * xi(9)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(3) * xi(9)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(8)**2 * s_1
    phi_x(77) = (-4.999999999999995e-01_rk) * xi(7)**2 * xi(12) * c_1 &
         + (+8.660254037844387e-01_rk) * xi(7)**2 * xi(12) * s_1 &
         + (-5.000000000000006e-01_rk) * xi(7)**2 * xi(11) * c_1 &
         + (+8.660254037844382e-01_rk) * xi(7)**2 * xi(11) * s_1 &
         + (-4.999999999999993e-01_rk) * xi(6)**2 * xi(13) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(6)**2 * xi(13) * s_1 &
         + (-5.000000000000002e-01_rk) * xi(6)**2 * xi(11) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(6)**2 * xi(11) * s_1 &
         + (+9.999999999999999e-01_rk) * xi(5)**2 * xi(13) * c_1 &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * xi(12) * c_1 &
         + (-4.999999999999993e-01_rk) * xi(4)**2 * xi(9) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(4)**2 * xi(9) * s_1 &
         + (-5.000000000000003e-01_rk) * xi(4)**2 * xi(8) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(4)**2 * xi(8) * s_1 &
         + (-4.999999999999997e-01_rk) * xi(3)**2 * xi(10) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(3)**2 * xi(10) * s_1 &
         + (-5.000000000000004e-01_rk) * xi(3)**2 * xi(8) * c_1 &
         + (+8.660254037844382e-01_rk) * xi(3)**2 * xi(8) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**2 * xi(10) * c_1 &
         + (+9.999999999999999e-01_rk) * xi(2)**2 * xi(9) * c_1
    phi_y(77) = (+8.660254037844385e-01_rk) * xi(7)**2 * xi(12) * c_1 &
         + (+4.999999999999998e-01_rk) * xi(7)**2 * xi(12) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(7)**2 * xi(11) * c_1 &
         + (+5.000000000000000e-01_rk) * xi(7)**2 * xi(11) * s_1 &
         + (-8.660254037844388e-01_rk) * xi(6)**2 * xi(13) * c_1 &
         + (+4.999999999999997e-01_rk) * xi(6)**2 * xi(13) * s_1 &
         + (-8.660254037844382e-01_rk) * xi(6)**2 * xi(11) * c_1 &
         + (+5.000000000000000e-01_rk) * xi(6)**2 * xi(11) * s_1 &
         + (-9.999999999999999e-01_rk) * xi(5)**2 * xi(13) * s_1 &
         + (-9.999999999999997e-01_rk) * xi(5)**2 * xi(12) * s_1 &
         + (+8.660254037844388e-01_rk) * xi(4)**2 * xi(9) * c_1 &
         + (-4.999999999999997e-01_rk) * xi(4)**2 * xi(9) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(4)**2 * xi(8) * c_1 &
         + (-5.000000000000000e-01_rk) * xi(4)**2 * xi(8) * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3)**2 * xi(10) * c_1 &
         + (-4.999999999999998e-01_rk) * xi(3)**2 * xi(10) * s_1 &
         + (-8.660254037844382e-01_rk) * xi(3)**2 * xi(8) * c_1 &
         + (-5.000000000000000e-01_rk) * xi(3)**2 * xi(8) * s_1 &
         + (+9.999999999999997e-01_rk) * xi(2)**2 * xi(10) * s_1 &
         + (+9.999999999999999e-01_rk) * xi(2)**2 * xi(9) * s_1
    phi_x(78) = (-8.660254037844385e-01_rk) * xi(7)**2 * xi(12) * c_1 &
         + (-4.999999999999998e-01_rk) * xi(7)**2 * xi(12) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(7)**2 * xi(11) * c_1 &
         + (+5.000000000000000e-01_rk) * xi(7)**2 * xi(11) * s_1 &
         + (-8.660254037844388e-01_rk) * xi(6)**2 * xi(13) * c_1 &
         + (+4.999999999999997e-01_rk) * xi(6)**2 * xi(13) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(6)**2 * xi(11) * c_1 &
         + (-5.000000000000000e-01_rk) * xi(6)**2 * xi(11) * s_1 &
         + (+9.999999999999999e-01_rk) * xi(5)**2 * xi(13) * s_1 &
         + (-9.999999999999997e-01_rk) * xi(5)**2 * xi(12) * s_1 &
         + (-8.660254037844388e-01_rk) * xi(4)**2 * xi(9) * c_1 &
         + (+4.999999999999997e-01_rk) * xi(4)**2 * xi(9) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(4)**2 * xi(8) * c_1 &
         + (-5.000000000000000e-01_rk) * xi(4)**2 * xi(8) * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3)**2 * xi(10) * c_1 &
         + (-4.999999999999998e-01_rk) * xi(3)**2 * xi(10) * s_1 &
         + (+8.660254037844382e-01_rk) * xi(3)**2 * xi(8) * c_1 &
         + (+5.000000000000000e-01_rk) * xi(3)**2 * xi(8) * s_1 &
         + (-9.999999999999997e-01_rk) * xi(2)**2 * xi(10) * s_1 &
         + (+9.999999999999999e-01_rk) * xi(2)**2 * xi(9) * s_1
    phi_y(78) = (-4.999999999999995e-01_rk) * xi(7)**2 * xi(12) * c_1 &
         + (+8.660254037844387e-01_rk) * xi(7)**2 * xi(12) * s_1 &
         + (+5.000000000000006e-01_rk) * xi(7)**2 * xi(11) * c_1 &
         + (-8.660254037844382e-01_rk) * xi(7)**2 * xi(11) * s_1 &
         + (+4.999999999999993e-01_rk) * xi(6)**2 * xi(13) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(6)**2 * xi(13) * s_1 &
         + (-5.000000000000002e-01_rk) * xi(6)**2 * xi(11) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(6)**2 * xi(11) * s_1 &
         + (+9.999999999999999e-01_rk) * xi(5)**2 * xi(13) * c_1 &
         + (-1.000000000000000e+00_rk) * xi(5)**2 * xi(12) * c_1 &
         + (-4.999999999999993e-01_rk) * xi(4)**2 * xi(9) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(4)**2 * xi(9) * s_1 &
         + (+5.000000000000003e-01_rk) * xi(4)**2 * xi(8) * c_1 &
         + (+8.660254037844385e-01_rk) * xi(4)**2 * xi(8) * s_1 &
         + (+4.999999999999997e-01_rk) * xi(3)**2 * xi(10) * c_1 &
         + (-8.660254037844385e-01_rk) * xi(3)**2 * xi(10) * s_1 &
         + (-5.000000000000004e-01_rk) * xi(3)**2 * xi(8) * c_1 &
         + (+8.660254037844382e-01_rk) * xi(3)**2 * xi(8) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**2 * xi(10) * c_1 &
         + (-9.999999999999999e-01_rk) * xi(2)**2 * xi(9) * c_1
    phi_x(79) = (-4.999999999999995e-01_rk) * xi(7) * xi(12)**2 * c_1 &
         + (+8.660254037844387e-01_rk) * xi(7) * xi(12)**2 * s_1 &
         + (-5.000000000000006e-01_rk) * xi(7) * xi(11)**2 * c_1 &
         + (+8.660254037844382e-01_rk) * xi(7) * xi(11)**2 * s_1 &
         + (-4.999999999999993e-01_rk) * xi(6) * xi(13)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(6) * xi(13)**2 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(6) * xi(11)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(6) * xi(11)**2 * s_1 &
         + (+9.999999999999999e-01_rk) * xi(5) * xi(13)**2 * c_1 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(12)**2 * c_1 &
         + (-4.999999999999993e-01_rk) * xi(4) * xi(9)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(4) * xi(9)**2 * s_1 &
         + (-5.000000000000003e-01_rk) * xi(4) * xi(8)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(4) * xi(8)**2 * s_1 &
         + (-4.999999999999997e-01_rk) * xi(3) * xi(10)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(3) * xi(10)**2 * s_1 &
         + (-5.000000000000004e-01_rk) * xi(3) * xi(8)**2 * c_1 &
         + (+8.660254037844382e-01_rk) * xi(3) * xi(8)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(10)**2 * c_1 &
         + (+9.999999999999999e-01_rk) * xi(2) * xi(9)**2 * c_1
    phi_y(79) = (+8.660254037844385e-01_rk) * xi(7) * xi(12)**2 * c_1 &
         + (+4.999999999999998e-01_rk) * xi(7) * xi(12)**2 * s_1 &
         + (+8.660254037844382e-01_rk) * xi(7) * xi(11)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(7) * xi(11)**2 * s_1 &
         + (-8.660254037844388e-01_rk) * xi(6) * xi(13)**2 * c_1 &
         + (+4.999999999999997e-01_rk) * xi(6) * xi(13)**2 * s_1 &
         + (-8.660254037844382e-01_rk) * xi(6) * xi(11)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(6) * xi(11)**2 * s_1 &
         + (-9.999999999999999e-01_rk) * xi(5) * xi(13)**2 * s_1 &
         + (-9.999999999999997e-01_rk) * xi(5) * xi(12)**2 * s_1 &
         + (+8.660254037844388e-01_rk) * xi(4) * xi(9)**2 * c_1 &
         + (-4.999999999999997e-01_rk) * xi(4) * xi(9)**2 * s_1 &
         + (+8.660254037844382e-01_rk) * xi(4) * xi(8)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(4) * xi(8)**2 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3) * xi(10)**2 * c_1 &
         + (-4.999999999999998e-01_rk) * xi(3) * xi(10)**2 * s_1 &
         + (-8.660254037844382e-01_rk) * xi(3) * xi(8)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(3) * xi(8)**2 * s_1 &
         + (+9.999999999999997e-01_rk) * xi(2) * xi(10)**2 * s_1 &
         + (+9.999999999999999e-01_rk) * xi(2) * xi(9)**2 * s_1
    phi_x(80) = (-8.660254037844385e-01_rk) * xi(7) * xi(12)**2 * c_1 &
         + (-4.999999999999998e-01_rk) * xi(7) * xi(12)**2 * s_1 &
         + (+8.660254037844382e-01_rk) * xi(7) * xi(11)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(7) * xi(11)**2 * s_1 &
         + (-8.660254037844388e-01_rk) * xi(6) * xi(13)**2 * c_1 &
         + (+4.999999999999997e-01_rk) * xi(6) * xi(13)**2 * s_1 &
         + (+8.660254037844382e-01_rk) * xi(6) * xi(11)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(6) * xi(11)**2 * s_1 &
         + (+9.999999999999999e-01_rk) * xi(5) * xi(13)**2 * s_1 &
         + (-9.999999999999997e-01_rk) * xi(5) * xi(12)**2 * s_1 &
         + (-8.660254037844388e-01_rk) * xi(4) * xi(9)**2 * c_1 &
         + (+4.999999999999997e-01_rk) * xi(4) * xi(9)**2 * s_1 &
         + (+8.660254037844382e-01_rk) * xi(4) * xi(8)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(4) * xi(8)**2 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3) * xi(10)**2 * c_1 &
         + (-4.999999999999998e-01_rk) * xi(3) * xi(10)**2 * s_1 &
         + (+8.660254037844382e-01_rk) * xi(3) * xi(8)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(3) * xi(8)**2 * s_1 &
         + (-9.999999999999997e-01_rk) * xi(2) * xi(10)**2 * s_1 &
         + (+9.999999999999999e-01_rk) * xi(2) * xi(9)**2 * s_1
    phi_y(80) = (-4.999999999999995e-01_rk) * xi(7) * xi(12)**2 * c_1 &
         + (+8.660254037844387e-01_rk) * xi(7) * xi(12)**2 * s_1 &
         + (+5.000000000000006e-01_rk) * xi(7) * xi(11)**2 * c_1 &
         + (-8.660254037844382e-01_rk) * xi(7) * xi(11)**2 * s_1 &
         + (+4.999999999999993e-01_rk) * xi(6) * xi(13)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(6) * xi(13)**2 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(6) * xi(11)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(6) * xi(11)**2 * s_1 &
         + (+9.999999999999999e-01_rk) * xi(5) * xi(13)**2 * c_1 &
         + (-1.000000000000000e+00_rk) * xi(5) * xi(12)**2 * c_1 &
         + (-4.999999999999993e-01_rk) * xi(4) * xi(9)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(4) * xi(9)**2 * s_1 &
         + (+5.000000000000003e-01_rk) * xi(4) * xi(8)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(4) * xi(8)**2 * s_1 &
         + (+4.999999999999997e-01_rk) * xi(3) * xi(10)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(3) * xi(10)**2 * s_1 &
         + (-5.000000000000004e-01_rk) * xi(3) * xi(8)**2 * c_1 &
         + (+8.660254037844382e-01_rk) * xi(3) * xi(8)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(10)**2 * c_1 &
         + (-9.999999999999999e-01_rk) * xi(2) * xi(9)**2 * c_1
    phi_x(81) = (+4.330127018922195e-01_rk) * xi(7)**2 * xi(17) * c_1 &
         + (-7.499999999999998e-01_rk) * xi(7)**2 * xi(17) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(7)**2 * xi(16) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(7)**2 * xi(16) * s_1 &
         + (-4.330127018922191e-01_rk) * xi(6)**2 * xi(17) * c_1 &
         + (-7.500000000000000e-01_rk) * xi(6)**2 * xi(17) * s_1 &
         + (+2.499999999999999e-01_rk) * xi(6)**2 * xi(16) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(6)**2 * xi(16) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * xi(16) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(4)**2 * xi(15) * c_1 &
         + (+7.500000000000000e-01_rk) * xi(4)**2 * xi(15) * s_1 &
         + (+2.499999999999999e-01_rk) * xi(4)**2 * xi(14) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(4)**2 * xi(14) * s_1 &
         + (-4.330127018922196e-01_rk) * xi(3)**2 * xi(15) * c_1 &
         + (+7.499999999999998e-01_rk) * xi(3)**2 * xi(15) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(3)**2 * xi(14) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(3)**2 * xi(14) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**2 * xi(14) * c_1
    phi_y(81) = (-7.500000000000000e-01_rk) * xi(7)**2 * xi(17) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(7)**2 * xi(17) * s_1 &
         + (-4.330127018922192e-01_rk) * xi(7)**2 * xi(16) * c_1 &
         + (-2.500000000000000e-01_rk) * xi(7)**2 * xi(16) * s_1 &
         + (-7.500000000000001e-01_rk) * xi(6)**2 * xi(17) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(6)**2 * xi(17) * s_1 &
         + (+4.330127018922193e-01_rk) * xi(6)**2 * xi(16) * c_1 &
         + (-2.499999999999999e-01_rk) * xi(6)**2 * xi(16) * s_1 &
         + (-9.999999999999998e-01_rk) * xi(5)**2 * xi(16) * s_1 &
         + (-7.500000000000001e-01_rk) * xi(4)**2 * xi(15) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(4)**2 * xi(15) * s_1 &
         + (-4.330127018922193e-01_rk) * xi(4)**2 * xi(14) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(4)**2 * xi(14) * s_1 &
         + (-7.500000000000000e-01_rk) * xi(3)**2 * xi(15) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(3)**2 * xi(15) * s_1 &
         + (+4.330127018922192e-01_rk) * xi(3)**2 * xi(14) * c_1 &
         + (+2.500000000000000e-01_rk) * xi(3)**2 * xi(14) * s_1 &
         + (+9.999999999999998e-01_rk) * xi(2)**2 * xi(14) * s_1
    phi_x(82) = (+4.330127018922192e-01_rk) * xi(7)**2 * xi(17) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(7)**2 * xi(17) * s_1 &
         + (-7.500000000000000e-01_rk) * xi(7)**2 * xi(16) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(7)**2 * xi(16) * s_1 &
         + (-4.330127018922193e-01_rk) * xi(6)**2 * xi(17) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(6)**2 * xi(17) * s_1 &
         + (-7.500000000000001e-01_rk) * xi(6)**2 * xi(16) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(6)**2 * xi(16) * s_1 &
         + (+9.999999999999998e-01_rk) * xi(5)**2 * xi(17) * s_1 &
         + (+4.330127018922193e-01_rk) * xi(4)**2 * xi(15) * c_1 &
         + (-2.499999999999999e-01_rk) * xi(4)**2 * xi(15) * s_1 &
         + (-7.500000000000001e-01_rk) * xi(4)**2 * xi(14) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(4)**2 * xi(14) * s_1 &
         + (-4.330127018922192e-01_rk) * xi(3)**2 * xi(15) * c_1 &
         + (-2.499999999999999e-01_rk) * xi(3)**2 * xi(15) * s_1 &
         + (-7.500000000000000e-01_rk) * xi(3)**2 * xi(14) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(3)**2 * xi(14) * s_1 &
         + (-9.999999999999998e-01_rk) * xi(2)**2 * xi(15) * s_1
    phi_y(82) = (+2.500000000000001e-01_rk) * xi(7)**2 * xi(17) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(7)**2 * xi(17) * s_1 &
         + (-4.330127018922195e-01_rk) * xi(7)**2 * xi(16) * c_1 &
         + (+7.499999999999998e-01_rk) * xi(7)**2 * xi(16) * s_1 &
         + (+2.499999999999999e-01_rk) * xi(6)**2 * xi(17) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(6)**2 * xi(17) * s_1 &
         + (+4.330127018922191e-01_rk) * xi(6)**2 * xi(16) * c_1 &
         + (+7.500000000000000e-01_rk) * xi(6)**2 * xi(16) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * xi(17) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(4)**2 * xi(15) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(4)**2 * xi(15) * s_1 &
         + (-4.330127018922193e-01_rk) * xi(4)**2 * xi(14) * c_1 &
         + (-7.500000000000000e-01_rk) * xi(4)**2 * xi(14) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(3)**2 * xi(15) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(3)**2 * xi(15) * s_1 &
         + (+4.330127018922196e-01_rk) * xi(3)**2 * xi(14) * c_1 &
         + (-7.499999999999998e-01_rk) * xi(3)**2 * xi(14) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**2 * xi(15) * c_1
    phi_x(83) = (-3.750000000000001e-01_rk) * xi(7) * xi(17)**2 * c_1 &
         + (+6.495190528383287e-01_rk) * xi(7) * xi(17)**2 * s_1 &
         + (-4.330127018922195e-01_rk) * xi(7) * xi(16) * xi(17) * c_1 &
         + (+7.499999999999998e-01_rk) * xi(7) * xi(16) * xi(17) * s_1 &
         + (-1.250000000000000e-01_rk) * xi(7) * xi(16)**2 * c_1 &
         + (+2.165063509461096e-01_rk) * xi(7) * xi(16)**2 * s_1 &
         + (-3.749999999999998e-01_rk) * xi(6) * xi(17)**2 * c_1 &
         + (-6.495190528383288e-01_rk) * xi(6) * xi(17)**2 * s_1 &
         + (+4.330127018922191e-01_rk) * xi(6) * xi(16) * xi(17) * c_1 &
         + (+7.499999999999998e-01_rk) * xi(6) * xi(16) * xi(17) * s_1 &
         + (-1.250000000000000e-01_rk) * xi(6) * xi(16)**2 * c_1 &
         + (-2.165063509461096e-01_rk) * xi(6) * xi(16)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(16)**2 * c_1 &
         + (-3.749999999999999e-01_rk) * xi(4) * xi(15)**2 * c_1 &
         + (-6.495190528383288e-01_rk) * xi(4) * xi(15)**2 * s_1 &
         + (-4.330127018922192e-01_rk) * xi(4) * xi(14) * xi(15) * c_1 &
         + (-7.499999999999998e-01_rk) * xi(4) * xi(14) * xi(15) * s_1 &
         + (-1.250000000000000e-01_rk) * xi(4) * xi(14)**2 * c_1 &
         + (-2.165063509461096e-01_rk) * xi(4) * xi(14)**2 * s_1 &
         + (-3.750000000000002e-01_rk) * xi(3) * xi(15)**2 * c_1 &
         + (+6.495190528383287e-01_rk) * xi(3) * xi(15)**2 * s_1 &
         + (+4.330127018922196e-01_rk) * xi(3) * xi(14) * xi(15) * c_1 &
         + (-7.499999999999997e-01_rk) * xi(3) * xi(14) * xi(15) * s_1 &
         + (-1.250000000000001e-01_rk) * xi(3) * xi(14)**2 * c_1 &
         + (+2.165063509461096e-01_rk) * xi(3) * xi(14)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(14)**2 * c_1
    phi_y(83) = (+6.495190528383288e-01_rk) * xi(7) * xi(17)**2 * c_1 &
         + (+3.749999999999999e-01_rk) * xi(7) * xi(17)**2 * s_1 &
         + (+7.499999999999998e-01_rk) * xi(7) * xi(16) * xi(17) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(7) * xi(16) * xi(17) * s_1 &
         + (+2.165063509461096e-01_rk) * xi(7) * xi(16)**2 * c_1 &
         + (+1.250000000000000e-01_rk) * xi(7) * xi(16)**2 * s_1 &
         + (-6.495190528383289e-01_rk) * xi(6) * xi(17)**2 * c_1 &
         + (+3.749999999999998e-01_rk) * xi(6) * xi(17)**2 * s_1 &
         + (+7.500000000000000e-01_rk) * xi(6) * xi(16) * xi(17) * c_1 &
         + (-4.330127018922191e-01_rk) * xi(6) * xi(16) * xi(17) * s_1 &
         + (-2.165063509461096e-01_rk) * xi(6) * xi(16)**2 * c_1 &
         + (+1.249999999999999e-01_rk) * xi(6) * xi(16)**2 * s_1 &
         + (-9.999999999999997e-01_rk) * xi(5) * xi(16)**2 * s_1 &
         + (+6.495190528383289e-01_rk) * xi(4) * xi(15)**2 * c_1 &
         + (-3.749999999999998e-01_rk) * xi(4) * xi(15)**2 * s_1 &
         + (+7.500000000000000e-01_rk) * xi(4) * xi(14) * xi(15) * c_1 &
         + (-4.330127018922191e-01_rk) * xi(4) * xi(14) * xi(15) * s_1 &
         + (+2.165063509461096e-01_rk) * xi(4) * xi(14)**2 * c_1 &
         + (-1.249999999999999e-01_rk) * xi(4) * xi(14)**2 * s_1 &
         + (-6.495190528383288e-01_rk) * xi(3) * xi(15)**2 * c_1 &
         + (-3.749999999999999e-01_rk) * xi(3) * xi(15)**2 * s_1 &
         + (+7.499999999999998e-01_rk) * xi(3) * xi(14) * xi(15) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(3) * xi(14) * xi(15) * s_1 &
         + (-2.165063509461096e-01_rk) * xi(3) * xi(14)**2 * c_1 &
         + (-1.250000000000000e-01_rk) * xi(3) * xi(14)**2 * s_1 &
         + (+9.999999999999997e-01_rk) * xi(2) * xi(14)**2 * s_1
    phi_x(84) = (-1.250000000000001e-01_rk) * xi(7) * xi(17)**2 * c_1 &
         + (+2.165063509461096e-01_rk) * xi(7) * xi(17)**2 * s_1 &
         + (+4.330127018922195e-01_rk) * xi(7) * xi(16) * xi(17) * c_1 &
         + (-7.499999999999998e-01_rk) * xi(7) * xi(16) * xi(17) * s_1 &
         + (-3.750000000000001e-01_rk) * xi(7) * xi(16)**2 * c_1 &
         + (+6.495190528383287e-01_rk) * xi(7) * xi(16)**2 * s_1 &
         + (-1.250000000000000e-01_rk) * xi(6) * xi(17)**2 * c_1 &
         + (-2.165063509461096e-01_rk) * xi(6) * xi(17)**2 * s_1 &
         + (-4.330127018922191e-01_rk) * xi(6) * xi(16) * xi(17) * c_1 &
         + (-7.499999999999998e-01_rk) * xi(6) * xi(16) * xi(17) * s_1 &
         + (-3.749999999999998e-01_rk) * xi(6) * xi(16)**2 * c_1 &
         + (-6.495190528383288e-01_rk) * xi(6) * xi(16)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(5) * xi(17)**2 * c_1 &
         + (-1.250000000000000e-01_rk) * xi(4) * xi(15)**2 * c_1 &
         + (-2.165063509461096e-01_rk) * xi(4) * xi(15)**2 * s_1 &
         + (+4.330127018922192e-01_rk) * xi(4) * xi(14) * xi(15) * c_1 &
         + (+7.499999999999998e-01_rk) * xi(4) * xi(14) * xi(15) * s_1 &
         + (-3.749999999999999e-01_rk) * xi(4) * xi(14)**2 * c_1 &
         + (-6.495190528383288e-01_rk) * xi(4) * xi(14)**2 * s_1 &
         + (-1.250000000000001e-01_rk) * xi(3) * xi(15)**2 * c_1 &
         + (+2.165063509461096e-01_rk) * xi(3) * xi(15)**2 * s_1 &
         + (-4.330127018922196e-01_rk) * xi(3) * xi(14) * xi(15) * c_1 &
         + (+7.499999999999997e-01_rk) * xi(3) * xi(14) * xi(15) * s_1 &
         + (-3.750000000000002e-01_rk) * xi(3) * xi(14)**2 * c_1 &
         + (+6.495190528383287e-01_rk) * xi(3) * xi(14)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(15)**2 * c_1
    phi_y(84) = (+2.165063509461096e-01_rk) * xi(7) * xi(17)**2 * c_1 &
         + (+1.250000000000000e-01_rk) * xi(7) * xi(17)**2 * s_1 &
         + (-7.499999999999998e-01_rk) * xi(7) * xi(16) * xi(17) * c_1 &
         + (-4.330127018922191e-01_rk) * xi(7) * xi(16) * xi(17) * s_1 &
         + (+6.495190528383288e-01_rk) * xi(7) * xi(16)**2 * c_1 &
         + (+3.749999999999999e-01_rk) * xi(7) * xi(16)**2 * s_1 &
         + (-2.165063509461096e-01_rk) * xi(6) * xi(17)**2 * c_1 &
         + (+1.249999999999999e-01_rk) * xi(6) * xi(17)**2 * s_1 &
         + (-7.500000000000000e-01_rk) * xi(6) * xi(16) * xi(17) * c_1 &
         + (+4.330127018922191e-01_rk) * xi(6) * xi(16) * xi(17) * s_1 &
         + (-6.495190528383289e-01_rk) * xi(6) * xi(16)**2 * c_1 &
         + (+3.749999999999998e-01_rk) * xi(6) * xi(16)**2 * s_1 &
         + (-9.999999999999993e-01_rk) * xi(5) * xi(17)**2 * s_1 &
         + (+2.165063509461096e-01_rk) * xi(4) * xi(15)**2 * c_1 &
         + (-1.249999999999999e-01_rk) * xi(4) * xi(15)**2 * s_1 &
         + (-7.500000000000000e-01_rk) * xi(4) * xi(14) * xi(15) * c_1 &
         + (+4.330127018922191e-01_rk) * xi(4) * xi(14) * xi(15) * s_1 &
         + (+6.495190528383289e-01_rk) * xi(4) * xi(14)**2 * c_1 &
         + (-3.749999999999998e-01_rk) * xi(4) * xi(14)**2 * s_1 &
         + (-2.165063509461096e-01_rk) * xi(3) * xi(15)**2 * c_1 &
         + (-1.250000000000000e-01_rk) * xi(3) * xi(15)**2 * s_1 &
         + (-7.499999999999998e-01_rk) * xi(3) * xi(14) * xi(15) * c_1 &
         + (-4.330127018922191e-01_rk) * xi(3) * xi(14) * xi(15) * s_1 &
         + (-6.495190528383288e-01_rk) * xi(3) * xi(14)**2 * c_1 &
         + (-3.749999999999999e-01_rk) * xi(3) * xi(14)**2 * s_1 &
         + (+9.999999999999993e-01_rk) * xi(2) * xi(15)**2 * s_1
    phi_x(85) = (+4.330127018922195e-01_rk) * xi(13)**2 * xi(17) * c_1 &
         + (-7.499999999999998e-01_rk) * xi(13)**2 * xi(17) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(13)**2 * xi(16) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(13)**2 * xi(16) * s_1 &
         + (-4.330127018922191e-01_rk) * xi(12)**2 * xi(17) * c_1 &
         + (-7.500000000000000e-01_rk) * xi(12)**2 * xi(17) * s_1 &
         + (+2.499999999999999e-01_rk) * xi(12)**2 * xi(16) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(12)**2 * xi(16) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(11)**2 * xi(16) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(10)**2 * xi(15) * c_1 &
         + (+7.500000000000000e-01_rk) * xi(10)**2 * xi(15) * s_1 &
         + (+2.499999999999999e-01_rk) * xi(10)**2 * xi(14) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(10)**2 * xi(14) * s_1 &
         + (-4.330127018922196e-01_rk) * xi(9)**2 * xi(15) * c_1 &
         + (+7.499999999999998e-01_rk) * xi(9)**2 * xi(15) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(9)**2 * xi(14) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(9)**2 * xi(14) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8)**2 * xi(14) * c_1
    phi_y(85) = (-7.500000000000000e-01_rk) * xi(13)**2 * xi(17) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(13)**2 * xi(17) * s_1 &
         + (-4.330127018922192e-01_rk) * xi(13)**2 * xi(16) * c_1 &
         + (-2.500000000000000e-01_rk) * xi(13)**2 * xi(16) * s_1 &
         + (-7.500000000000001e-01_rk) * xi(12)**2 * xi(17) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(12)**2 * xi(17) * s_1 &
         + (+4.330127018922193e-01_rk) * xi(12)**2 * xi(16) * c_1 &
         + (-2.499999999999999e-01_rk) * xi(12)**2 * xi(16) * s_1 &
         + (-9.999999999999998e-01_rk) * xi(11)**2 * xi(16) * s_1 &
         + (-7.500000000000001e-01_rk) * xi(10)**2 * xi(15) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(10)**2 * xi(15) * s_1 &
         + (-4.330127018922193e-01_rk) * xi(10)**2 * xi(14) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(10)**2 * xi(14) * s_1 &
         + (-7.500000000000000e-01_rk) * xi(9)**2 * xi(15) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(9)**2 * xi(15) * s_1 &
         + (+4.330127018922192e-01_rk) * xi(9)**2 * xi(14) * c_1 &
         + (+2.500000000000000e-01_rk) * xi(9)**2 * xi(14) * s_1 &
         + (+9.999999999999998e-01_rk) * xi(8)**2 * xi(14) * s_1
    phi_x(86) = (+4.330127018922192e-01_rk) * xi(13)**2 * xi(17) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(13)**2 * xi(17) * s_1 &
         + (-7.500000000000000e-01_rk) * xi(13)**2 * xi(16) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(13)**2 * xi(16) * s_1 &
         + (-4.330127018922193e-01_rk) * xi(12)**2 * xi(17) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(12)**2 * xi(17) * s_1 &
         + (-7.500000000000001e-01_rk) * xi(12)**2 * xi(16) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(12)**2 * xi(16) * s_1 &
         + (+9.999999999999998e-01_rk) * xi(11)**2 * xi(17) * s_1 &
         + (+4.330127018922193e-01_rk) * xi(10)**2 * xi(15) * c_1 &
         + (-2.499999999999999e-01_rk) * xi(10)**2 * xi(15) * s_1 &
         + (-7.500000000000001e-01_rk) * xi(10)**2 * xi(14) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(10)**2 * xi(14) * s_1 &
         + (-4.330127018922192e-01_rk) * xi(9)**2 * xi(15) * c_1 &
         + (-2.499999999999999e-01_rk) * xi(9)**2 * xi(15) * s_1 &
         + (-7.500000000000000e-01_rk) * xi(9)**2 * xi(14) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(9)**2 * xi(14) * s_1 &
         + (-9.999999999999998e-01_rk) * xi(8)**2 * xi(15) * s_1
    phi_y(86) = (+2.500000000000001e-01_rk) * xi(13)**2 * xi(17) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(13)**2 * xi(17) * s_1 &
         + (-4.330127018922195e-01_rk) * xi(13)**2 * xi(16) * c_1 &
         + (+7.499999999999998e-01_rk) * xi(13)**2 * xi(16) * s_1 &
         + (+2.499999999999999e-01_rk) * xi(12)**2 * xi(17) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(12)**2 * xi(17) * s_1 &
         + (+4.330127018922191e-01_rk) * xi(12)**2 * xi(16) * c_1 &
         + (+7.500000000000000e-01_rk) * xi(12)**2 * xi(16) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(11)**2 * xi(17) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(10)**2 * xi(15) * c_1 &
         + (+4.330127018922193e-01_rk) * xi(10)**2 * xi(15) * s_1 &
         + (-4.330127018922193e-01_rk) * xi(10)**2 * xi(14) * c_1 &
         + (-7.500000000000000e-01_rk) * xi(10)**2 * xi(14) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(9)**2 * xi(15) * c_1 &
         + (-4.330127018922192e-01_rk) * xi(9)**2 * xi(15) * s_1 &
         + (+4.330127018922196e-01_rk) * xi(9)**2 * xi(14) * c_1 &
         + (-7.499999999999998e-01_rk) * xi(9)**2 * xi(14) * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8)**2 * xi(15) * c_1
    phi_x(87) = (-3.750000000000001e-01_rk) * xi(13) * xi(17)**2 * c_1 &
         + (+6.495190528383287e-01_rk) * xi(13) * xi(17)**2 * s_1 &
         + (-4.330127018922195e-01_rk) * xi(13) * xi(16) * xi(17) * c_1 &
         + (+7.499999999999998e-01_rk) * xi(13) * xi(16) * xi(17) * s_1 &
         + (-1.250000000000000e-01_rk) * xi(13) * xi(16)**2 * c_1 &
         + (+2.165063509461096e-01_rk) * xi(13) * xi(16)**2 * s_1 &
         + (-3.749999999999998e-01_rk) * xi(12) * xi(17)**2 * c_1 &
         + (-6.495190528383288e-01_rk) * xi(12) * xi(17)**2 * s_1 &
         + (+4.330127018922191e-01_rk) * xi(12) * xi(16) * xi(17) * c_1 &
         + (+7.499999999999998e-01_rk) * xi(12) * xi(16) * xi(17) * s_1 &
         + (-1.250000000000000e-01_rk) * xi(12) * xi(16)**2 * c_1 &
         + (-2.165063509461096e-01_rk) * xi(12) * xi(16)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(11) * xi(16)**2 * c_1 &
         + (-3.749999999999999e-01_rk) * xi(10) * xi(15)**2 * c_1 &
         + (-6.495190528383288e-01_rk) * xi(10) * xi(15)**2 * s_1 &
         + (-4.330127018922192e-01_rk) * xi(10) * xi(14) * xi(15) * c_1 &
         + (-7.499999999999998e-01_rk) * xi(10) * xi(14) * xi(15) * s_1 &
         + (-1.250000000000000e-01_rk) * xi(10) * xi(14)**2 * c_1 &
         + (-2.165063509461096e-01_rk) * xi(10) * xi(14)**2 * s_1 &
         + (-3.750000000000002e-01_rk) * xi(9) * xi(15)**2 * c_1 &
         + (+6.495190528383287e-01_rk) * xi(9) * xi(15)**2 * s_1 &
         + (+4.330127018922196e-01_rk) * xi(9) * xi(14) * xi(15) * c_1 &
         + (-7.499999999999997e-01_rk) * xi(9) * xi(14) * xi(15) * s_1 &
         + (-1.250000000000001e-01_rk) * xi(9) * xi(14)**2 * c_1 &
         + (+2.165063509461096e-01_rk) * xi(9) * xi(14)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(14)**2 * c_1
    phi_y(87) = (+6.495190528383288e-01_rk) * xi(13) * xi(17)**2 * c_1 &
         + (+3.749999999999999e-01_rk) * xi(13) * xi(17)**2 * s_1 &
         + (+7.499999999999998e-01_rk) * xi(13) * xi(16) * xi(17) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(13) * xi(16) * xi(17) * s_1 &
         + (+2.165063509461096e-01_rk) * xi(13) * xi(16)**2 * c_1 &
         + (+1.250000000000000e-01_rk) * xi(13) * xi(16)**2 * s_1 &
         + (-6.495190528383289e-01_rk) * xi(12) * xi(17)**2 * c_1 &
         + (+3.749999999999998e-01_rk) * xi(12) * xi(17)**2 * s_1 &
         + (+7.500000000000000e-01_rk) * xi(12) * xi(16) * xi(17) * c_1 &
         + (-4.330127018922191e-01_rk) * xi(12) * xi(16) * xi(17) * s_1 &
         + (-2.165063509461096e-01_rk) * xi(12) * xi(16)**2 * c_1 &
         + (+1.249999999999999e-01_rk) * xi(12) * xi(16)**2 * s_1 &
         + (-9.999999999999997e-01_rk) * xi(11) * xi(16)**2 * s_1 &
         + (+6.495190528383289e-01_rk) * xi(10) * xi(15)**2 * c_1 &
         + (-3.749999999999998e-01_rk) * xi(10) * xi(15)**2 * s_1 &
         + (+7.500000000000000e-01_rk) * xi(10) * xi(14) * xi(15) * c_1 &
         + (-4.330127018922191e-01_rk) * xi(10) * xi(14) * xi(15) * s_1 &
         + (+2.165063509461096e-01_rk) * xi(10) * xi(14)**2 * c_1 &
         + (-1.249999999999999e-01_rk) * xi(10) * xi(14)**2 * s_1 &
         + (-6.495190528383288e-01_rk) * xi(9) * xi(15)**2 * c_1 &
         + (-3.749999999999999e-01_rk) * xi(9) * xi(15)**2 * s_1 &
         + (+7.499999999999998e-01_rk) * xi(9) * xi(14) * xi(15) * c_1 &
         + (+4.330127018922192e-01_rk) * xi(9) * xi(14) * xi(15) * s_1 &
         + (-2.165063509461096e-01_rk) * xi(9) * xi(14)**2 * c_1 &
         + (-1.250000000000000e-01_rk) * xi(9) * xi(14)**2 * s_1 &
         + (+9.999999999999997e-01_rk) * xi(8) * xi(14)**2 * s_1
    phi_x(88) = (-1.250000000000001e-01_rk) * xi(13) * xi(17)**2 * c_1 &
         + (+2.165063509461096e-01_rk) * xi(13) * xi(17)**2 * s_1 &
         + (+4.330127018922195e-01_rk) * xi(13) * xi(16) * xi(17) * c_1 &
         + (-7.499999999999998e-01_rk) * xi(13) * xi(16) * xi(17) * s_1 &
         + (-3.750000000000001e-01_rk) * xi(13) * xi(16)**2 * c_1 &
         + (+6.495190528383287e-01_rk) * xi(13) * xi(16)**2 * s_1 &
         + (-1.250000000000000e-01_rk) * xi(12) * xi(17)**2 * c_1 &
         + (-2.165063509461096e-01_rk) * xi(12) * xi(17)**2 * s_1 &
         + (-4.330127018922191e-01_rk) * xi(12) * xi(16) * xi(17) * c_1 &
         + (-7.499999999999998e-01_rk) * xi(12) * xi(16) * xi(17) * s_1 &
         + (-3.749999999999998e-01_rk) * xi(12) * xi(16)**2 * c_1 &
         + (-6.495190528383288e-01_rk) * xi(12) * xi(16)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(11) * xi(17)**2 * c_1 &
         + (-1.250000000000000e-01_rk) * xi(10) * xi(15)**2 * c_1 &
         + (-2.165063509461096e-01_rk) * xi(10) * xi(15)**2 * s_1 &
         + (+4.330127018922192e-01_rk) * xi(10) * xi(14) * xi(15) * c_1 &
         + (+7.499999999999998e-01_rk) * xi(10) * xi(14) * xi(15) * s_1 &
         + (-3.749999999999999e-01_rk) * xi(10) * xi(14)**2 * c_1 &
         + (-6.495190528383288e-01_rk) * xi(10) * xi(14)**2 * s_1 &
         + (-1.250000000000001e-01_rk) * xi(9) * xi(15)**2 * c_1 &
         + (+2.165063509461096e-01_rk) * xi(9) * xi(15)**2 * s_1 &
         + (-4.330127018922196e-01_rk) * xi(9) * xi(14) * xi(15) * c_1 &
         + (+7.499999999999997e-01_rk) * xi(9) * xi(14) * xi(15) * s_1 &
         + (-3.750000000000002e-01_rk) * xi(9) * xi(14)**2 * c_1 &
         + (+6.495190528383287e-01_rk) * xi(9) * xi(14)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(15)**2 * c_1
    phi_y(88) = (+2.165063509461096e-01_rk) * xi(13) * xi(17)**2 * c_1 &
         + (+1.250000000000000e-01_rk) * xi(13) * xi(17)**2 * s_1 &
         + (-7.499999999999998e-01_rk) * xi(13) * xi(16) * xi(17) * c_1 &
         + (-4.330127018922191e-01_rk) * xi(13) * xi(16) * xi(17) * s_1 &
         + (+6.495190528383288e-01_rk) * xi(13) * xi(16)**2 * c_1 &
         + (+3.749999999999999e-01_rk) * xi(13) * xi(16)**2 * s_1 &
         + (-2.165063509461096e-01_rk) * xi(12) * xi(17)**2 * c_1 &
         + (+1.249999999999999e-01_rk) * xi(12) * xi(17)**2 * s_1 &
         + (-7.500000000000000e-01_rk) * xi(12) * xi(16) * xi(17) * c_1 &
         + (+4.330127018922191e-01_rk) * xi(12) * xi(16) * xi(17) * s_1 &
         + (-6.495190528383289e-01_rk) * xi(12) * xi(16)**2 * c_1 &
         + (+3.749999999999998e-01_rk) * xi(12) * xi(16)**2 * s_1 &
         + (-9.999999999999993e-01_rk) * xi(11) * xi(17)**2 * s_1 &
         + (+2.165063509461096e-01_rk) * xi(10) * xi(15)**2 * c_1 &
         + (-1.249999999999999e-01_rk) * xi(10) * xi(15)**2 * s_1 &
         + (-7.500000000000000e-01_rk) * xi(10) * xi(14) * xi(15) * c_1 &
         + (+4.330127018922191e-01_rk) * xi(10) * xi(14) * xi(15) * s_1 &
         + (+6.495190528383289e-01_rk) * xi(10) * xi(14)**2 * c_1 &
         + (-3.749999999999998e-01_rk) * xi(10) * xi(14)**2 * s_1 &
         + (-2.165063509461096e-01_rk) * xi(9) * xi(15)**2 * c_1 &
         + (-1.250000000000000e-01_rk) * xi(9) * xi(15)**2 * s_1 &
         + (-7.499999999999998e-01_rk) * xi(9) * xi(14) * xi(15) * c_1 &
         + (-4.330127018922191e-01_rk) * xi(9) * xi(14) * xi(15) * s_1 &
         + (-6.495190528383288e-01_rk) * xi(9) * xi(14)**2 * c_1 &
         + (-3.749999999999999e-01_rk) * xi(9) * xi(14)**2 * s_1 &
         + (+9.999999999999993e-01_rk) * xi(8) * xi(15)**2 * s_1
    phi_x(89) = (-5.000000000000000e-01_rk) * xi(4) * xi(7)**2 * c_3 &
         + (+1.000000000000000e+00_rk) * xi(4) * xi(6)**2 * c_3 &
         + (-5.000000000000000e-01_rk) * xi(4) * xi(5)**2 * c_3 &
         + (-5.000000000000000e-01_rk) * xi(4)**2 * xi(7) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(4)**2 * xi(6) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(4)**2 * xi(5) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(3) * xi(7)**2 * c_3 &
         + (-5.000000000000000e-01_rk) * xi(3) * xi(6)**2 * c_3 &
         + (-5.000000000000000e-01_rk) * xi(3) * xi(5)**2 * c_3 &
         + (+1.000000000000000e+00_rk) * xi(3)**2 * xi(7) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(3)**2 * xi(6) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(3)**2 * xi(5) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(2) * xi(7)**2 * c_3 &
         + (-5.000000000000000e-01_rk) * xi(2) * xi(6)**2 * c_3 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(5)**2 * c_3 &
         + (-5.000000000000000e-01_rk) * xi(2)**2 * xi(7) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(2)**2 * xi(6) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(2)**2 * xi(5) * c_3
    phi_y(89) = (+8.660254037844386e-01_rk) * xi(4) * xi(7)**2 * c_3 &
         + (-8.660254037844386e-01_rk) * xi(4) * xi(5)**2 * c_3 &
         + (+8.660254037844386e-01_rk) * xi(4)**2 * xi(7) * c_3 &
         + (-8.660254037844386e-01_rk) * xi(4)**2 * xi(5) * c_3 &
         + (-8.660254037844386e-01_rk) * xi(3) * xi(6)**2 * c_3 &
         + (+8.660254037844386e-01_rk) * xi(3) * xi(5)**2 * c_3 &
         + (-8.660254037844386e-01_rk) * xi(3)**2 * xi(6) * c_3 &
         + (+8.660254037844386e-01_rk) * xi(3)**2 * xi(5) * c_3 &
         + (-8.660254037844386e-01_rk) * xi(2) * xi(7)**2 * c_3 &
         + (+8.660254037844386e-01_rk) * xi(2) * xi(6)**2 * c_3 &
         + (-8.660254037844386e-01_rk) * xi(2)**2 * xi(7) * c_3 &
         + (+8.660254037844386e-01_rk) * xi(2)**2 * xi(6) * c_3
    phi_x(90) = (-8.660254037844386e-01_rk) * xi(4) * xi(7)**2 * s_3 &
         + (+8.660254037844386e-01_rk) * xi(4) * xi(5)**2 * s_3 &
         + (+8.660254037844386e-01_rk) * xi(4)**2 * xi(7) * s_3 &
         + (-8.660254037844386e-01_rk) * xi(4)**2 * xi(5) * s_3 &
         + (+8.660254037844386e-01_rk) * xi(3) * xi(6)**2 * s_3 &
         + (-8.660254037844386e-01_rk) * xi(3) * xi(5)**2 * s_3 &
         + (-8.660254037844386e-01_rk) * xi(3)**2 * xi(6) * s_3 &
         + (+8.660254037844386e-01_rk) * xi(3)**2 * xi(5) * s_3 &
         + (+8.660254037844386e-01_rk) * xi(2) * xi(7)**2 * s_3 &
         + (-8.660254037844386e-01_rk) * xi(2) * xi(6)**2 * s_3 &
         + (-8.660254037844386e-01_rk) * xi(2)**2 * xi(7) * s_3 &
         + (+8.660254037844386e-01_rk) * xi(2)**2 * xi(6) * s_3
    phi_y(90) = (-5.000000000000000e-01_rk) * xi(4) * xi(7)**2 * s_3 &
         + (+1.000000000000000e+00_rk) * xi(4) * xi(6)**2 * s_3 &
         + (-5.000000000000000e-01_rk) * xi(4) * xi(5)**2 * s_3 &
         + (+5.000000000000000e-01_rk) * xi(4)**2 * xi(7) * s_3 &
         + (-1.000000000000000e+00_rk) * xi(4)**2 * xi(6) * s_3 &
         + (+5.000000000000000e-01_rk) * xi(4)**2 * xi(5) * s_3 &
         + (+1.000000000000000e+00_rk) * xi(3) * xi(7)**2 * s_3 &
         + (-5.000000000000000e-01_rk) * xi(3) * xi(6)**2 * s_3 &
         + (-5.000000000000000e-01_rk) * xi(3) * xi(5)**2 * s_3 &
         + (-1.000000000000000e+00_rk) * xi(3)**2 * xi(7) * s_3 &
         + (+5.000000000000000e-01_rk) * xi(3)**2 * xi(6) * s_3 &
         + (+5.000000000000000e-01_rk) * xi(3)**2 * xi(5) * s_3 &
         + (-5.000000000000000e-01_rk) * xi(2) * xi(7)**2 * s_3 &
         + (-5.000000000000000e-01_rk) * xi(2) * xi(6)**2 * s_3 &
         + (+1.000000000000000e+00_rk) * xi(2) * xi(5)**2 * s_3 &
         + (+5.000000000000000e-01_rk) * xi(2)**2 * xi(7) * s_3 &
         + (+5.000000000000000e-01_rk) * xi(2)**2 * xi(6) * s_3 &
         + (-1.000000000000000e+00_rk) * xi(2)**2 * xi(5) * s_3
    phi_x(91) = (-5.000000000000000e-01_rk) * xi(10) * xi(13)**2 * c_3 &
         + (+1.000000000000000e+00_rk) * xi(10) * xi(12)**2 * c_3 &
         + (-5.000000000000000e-01_rk) * xi(10) * xi(11)**2 * c_3 &
         + (-5.000000000000000e-01_rk) * xi(10)**2 * xi(13) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(10)**2 * xi(12) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(10)**2 * xi(11) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(9) * xi(13)**2 * c_3 &
         + (-5.000000000000000e-01_rk) * xi(9) * xi(12)**2 * c_3 &
         + (-5.000000000000000e-01_rk) * xi(9) * xi(11)**2 * c_3 &
         + (+1.000000000000000e+00_rk) * xi(9)**2 * xi(13) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(9)**2 * xi(12) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(9)**2 * xi(11) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(8) * xi(13)**2 * c_3 &
         + (-5.000000000000000e-01_rk) * xi(8) * xi(12)**2 * c_3 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(11)**2 * c_3 &
         + (-5.000000000000000e-01_rk) * xi(8)**2 * xi(13) * c_3 &
         + (-5.000000000000000e-01_rk) * xi(8)**2 * xi(12) * c_3 &
         + (+1.000000000000000e+00_rk) * xi(8)**2 * xi(11) * c_3
    phi_y(91) = (+8.660254037844386e-01_rk) * xi(10) * xi(13)**2 * c_3 &
         + (-8.660254037844386e-01_rk) * xi(10) * xi(11)**2 * c_3 &
         + (+8.660254037844386e-01_rk) * xi(10)**2 * xi(13) * c_3 &
         + (-8.660254037844386e-01_rk) * xi(10)**2 * xi(11) * c_3 &
         + (-8.660254037844386e-01_rk) * xi(9) * xi(12)**2 * c_3 &
         + (+8.660254037844386e-01_rk) * xi(9) * xi(11)**2 * c_3 &
         + (-8.660254037844386e-01_rk) * xi(9)**2 * xi(12) * c_3 &
         + (+8.660254037844386e-01_rk) * xi(9)**2 * xi(11) * c_3 &
         + (-8.660254037844386e-01_rk) * xi(8) * xi(13)**2 * c_3 &
         + (+8.660254037844386e-01_rk) * xi(8) * xi(12)**2 * c_3 &
         + (-8.660254037844386e-01_rk) * xi(8)**2 * xi(13) * c_3 &
         + (+8.660254037844386e-01_rk) * xi(8)**2 * xi(12) * c_3
    phi_x(92) = (-8.660254037844386e-01_rk) * xi(10) * xi(13)**2 * s_3 &
         + (+8.660254037844386e-01_rk) * xi(10) * xi(11)**2 * s_3 &
         + (+8.660254037844386e-01_rk) * xi(10)**2 * xi(13) * s_3 &
         + (-8.660254037844386e-01_rk) * xi(10)**2 * xi(11) * s_3 &
         + (+8.660254037844386e-01_rk) * xi(9) * xi(12)**2 * s_3 &
         + (-8.660254037844386e-01_rk) * xi(9) * xi(11)**2 * s_3 &
         + (-8.660254037844386e-01_rk) * xi(9)**2 * xi(12) * s_3 &
         + (+8.660254037844386e-01_rk) * xi(9)**2 * xi(11) * s_3 &
         + (+8.660254037844386e-01_rk) * xi(8) * xi(13)**2 * s_3 &
         + (-8.660254037844386e-01_rk) * xi(8) * xi(12)**2 * s_3 &
         + (-8.660254037844386e-01_rk) * xi(8)**2 * xi(13) * s_3 &
         + (+8.660254037844386e-01_rk) * xi(8)**2 * xi(12) * s_3
    phi_y(92) = (-5.000000000000000e-01_rk) * xi(10) * xi(13)**2 * s_3 &
         + (+1.000000000000000e+00_rk) * xi(10) * xi(12)**2 * s_3 &
         + (-5.000000000000000e-01_rk) * xi(10) * xi(11)**2 * s_3 &
         + (+5.000000000000000e-01_rk) * xi(10)**2 * xi(13) * s_3 &
         + (-1.000000000000000e+00_rk) * xi(10)**2 * xi(12) * s_3 &
         + (+5.000000000000000e-01_rk) * xi(10)**2 * xi(11) * s_3 &
         + (+1.000000000000000e+00_rk) * xi(9) * xi(13)**2 * s_3 &
         + (-5.000000000000000e-01_rk) * xi(9) * xi(12)**2 * s_3 &
         + (-5.000000000000000e-01_rk) * xi(9) * xi(11)**2 * s_3 &
         + (-1.000000000000000e+00_rk) * xi(9)**2 * xi(13) * s_3 &
         + (+5.000000000000000e-01_rk) * xi(9)**2 * xi(12) * s_3 &
         + (+5.000000000000000e-01_rk) * xi(9)**2 * xi(11) * s_3 &
         + (-5.000000000000000e-01_rk) * xi(8) * xi(13)**2 * s_3 &
         + (-5.000000000000000e-01_rk) * xi(8) * xi(12)**2 * s_3 &
         + (+1.000000000000000e+00_rk) * xi(8) * xi(11)**2 * s_3 &
         + (+5.000000000000000e-01_rk) * xi(8)**2 * xi(13) * s_3 &
         + (+5.000000000000000e-01_rk) * xi(8)**2 * xi(12) * s_3 &
         + (-1.000000000000000e+00_rk) * xi(8)**2 * xi(11) * s_3
    phi_x(93) = (-5.000000000000000e-01_rk) * xi(7)**4 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(7)**4 * s_1 &
         + (-4.999999999999997e-01_rk) * xi(6)**4 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(6)**4 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(5)**4 * c_1 &
         + (-4.999999999999997e-01_rk) * xi(4)**4 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(4)**4 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(3)**4 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(3)**4 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**4 * c_1
    phi_y(93) = (+8.660254037844385e-01_rk) * xi(7)**4 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(7)**4 * s_1 &
         + (-8.660254037844386e-01_rk) * xi(6)**4 * c_1 &
         + (+4.999999999999999e-01_rk) * xi(6)**4 * s_1 &
         + (-1.000000000000000e+00_rk) * xi(5)**4 * s_1 &
         + (+8.660254037844386e-01_rk) * xi(4)**4 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(4)**4 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3)**4 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(3)**4 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**4 * s_1
    phi_x(94) = (-5.000000000000000e-01_rk) * xi(13)**4 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(13)**4 * s_1 &
         + (-4.999999999999997e-01_rk) * xi(12)**4 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(12)**4 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(11)**4 * c_1 &
         + (-4.999999999999997e-01_rk) * xi(10)**4 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(10)**4 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(9)**4 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(9)**4 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8)**4 * c_1
    phi_y(94) = (+8.660254037844385e-01_rk) * xi(13)**4 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(13)**4 * s_1 &
         + (-8.660254037844386e-01_rk) * xi(12)**4 * c_1 &
         + (+4.999999999999999e-01_rk) * xi(12)**4 * s_1 &
         + (-1.000000000000000e+00_rk) * xi(11)**4 * s_1 &
         + (+8.660254037844386e-01_rk) * xi(10)**4 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(10)**4 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(9)**4 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(9)**4 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8)**4 * s_1
    phi_x(95) = (-2.500000000000001e-01_rk) * xi(17)**4 * c_1 &
         + (+9.999999999999997e-01_rk) * xi(16) * xi(17)**3 * s_1 &
         + (-5.000000000000001e-01_rk) * xi(16)**2 * xi(17)**2 * c_1 &
         + (+3.333333333333333e-01_rk) * xi(16)**3 * xi(17) * s_1 &
         + (+4.166666666666669e-01_rk) * xi(16)**4 * c_1 &
         + (-2.500000000000001e-01_rk) * xi(15)**4 * c_1 &
         + (-1.000000000000000e+00_rk) * xi(14) * xi(15)**3 * s_1 &
         + (-5.000000000000001e-01_rk) * xi(14)**2 * xi(15)**2 * c_1 &
         + (-3.333333333333332e-01_rk) * xi(14)**3 * xi(15) * s_1 &
         + (+4.166666666666669e-01_rk) * xi(14)**4 * c_1
    phi_y(95) = (+2.499999999999999e-01_rk) * xi(17)**4 * s_1 &
         + (+9.999999999999999e-01_rk) * xi(16) * xi(17)**3 * c_1 &
         + (+4.999999999999997e-01_rk) * xi(16)**2 * xi(17)**2 * s_1 &
         + (+3.333333333333333e-01_rk) * xi(16)**3 * xi(17) * c_1 &
         + (-4.166666666666665e-01_rk) * xi(16)**4 * s_1 &
         + (-2.499999999999999e-01_rk) * xi(15)**4 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(14) * xi(15)**3 * c_1 &
         + (-4.999999999999997e-01_rk) * xi(14)**2 * xi(15)**2 * s_1 &
         + (+3.333333333333333e-01_rk) * xi(14)**3 * xi(15) * c_1 &
         + (+4.166666666666665e-01_rk) * xi(14)**4 * s_1
    phi_x(96) = (-4.166666666666669e-01_rk) * xi(17)**4 * c_1 &
         + (+3.333333333333332e-01_rk) * xi(16) * xi(17)**3 * s_1 &
         + (+5.000000000000001e-01_rk) * xi(16)**2 * xi(17)**2 * c_1 &
         + (+9.999999999999997e-01_rk) * xi(16)**3 * xi(17) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(16)**4 * c_1 &
         + (-4.166666666666669e-01_rk) * xi(15)**4 * c_1 &
         + (-3.333333333333331e-01_rk) * xi(14) * xi(15)**3 * s_1 &
         + (+5.000000000000001e-01_rk) * xi(14)**2 * xi(15)**2 * c_1 &
         + (-9.999999999999997e-01_rk) * xi(14)**3 * xi(15) * s_1 &
         + (+2.500000000000001e-01_rk) * xi(14)**4 * c_1
    phi_y(96) = (+4.166666666666663e-01_rk) * xi(17)**4 * s_1 &
         + (+3.333333333333333e-01_rk) * xi(16) * xi(17)**3 * c_1 &
         + (-4.999999999999997e-01_rk) * xi(16)**2 * xi(17)**2 * s_1 &
         + (+9.999999999999999e-01_rk) * xi(16)**3 * xi(17) * c_1 &
         + (-2.499999999999999e-01_rk) * xi(16)**4 * s_1 &
         + (-4.166666666666662e-01_rk) * xi(15)**4 * s_1 &
         + (+3.333333333333333e-01_rk) * xi(14) * xi(15)**3 * c_1 &
         + (+4.999999999999997e-01_rk) * xi(14)**2 * xi(15)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(14)**3 * xi(15) * c_1 &
         + (+2.499999999999999e-01_rk) * xi(14)**4 * s_1
    phi_x(97) = (-5.000000000000000e-01_rk) * xi(1)**2 * xi(7)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(1)**2 * xi(7)**2 * s_1 &
         + (-4.999999999999997e-01_rk) * xi(1)**2 * xi(6)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(1)**2 * xi(6)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(5)**2 * c_1 &
         + (-4.999999999999997e-01_rk) * xi(1)**2 * xi(4)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(1)**2 * xi(4)**2 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(1)**2 * xi(3)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(1)**2 * xi(3)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(2)**2 * c_1
    phi_y(97) = (+8.660254037844385e-01_rk) * xi(1)**2 * xi(7)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(1)**2 * xi(7)**2 * s_1 &
         + (-8.660254037844386e-01_rk) * xi(1)**2 * xi(6)**2 * c_1 &
         + (+4.999999999999999e-01_rk) * xi(1)**2 * xi(6)**2 * s_1 &
         + (-1.000000000000000e+00_rk) * xi(1)**2 * xi(5)**2 * s_1 &
         + (+8.660254037844386e-01_rk) * xi(1)**2 * xi(4)**2 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(1)**2 * xi(4)**2 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(1)**2 * xi(3)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(1)**2 * xi(3)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(2)**2 * s_1
    phi_x(98) = (-5.000000000000000e-01_rk) * xi(1)**2 * xi(13)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(1)**2 * xi(13)**2 * s_1 &
         + (-4.999999999999997e-01_rk) * xi(1)**2 * xi(12)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(1)**2 * xi(12)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(11)**2 * c_1 &
         + (-4.999999999999997e-01_rk) * xi(1)**2 * xi(10)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(1)**2 * xi(10)**2 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(1)**2 * xi(9)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(1)**2 * xi(9)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(8)**2 * c_1
    phi_y(98) = (+8.660254037844385e-01_rk) * xi(1)**2 * xi(13)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(1)**2 * xi(13)**2 * s_1 &
         + (-8.660254037844386e-01_rk) * xi(1)**2 * xi(12)**2 * c_1 &
         + (+4.999999999999999e-01_rk) * xi(1)**2 * xi(12)**2 * s_1 &
         + (-1.000000000000000e+00_rk) * xi(1)**2 * xi(11)**2 * s_1 &
         + (+8.660254037844386e-01_rk) * xi(1)**2 * xi(10)**2 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(1)**2 * xi(10)**2 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(1)**2 * xi(9)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(1)**2 * xi(9)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(1)**2 * xi(8)**2 * s_1
    phi_x(99) = (+9.999999999999991e-01_rk) * xi(6)**2 * xi(7)**2 * c_1 &
         + (-4.999999999999995e-01_rk) * xi(5)**2 * xi(7)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(5)**2 * xi(7)**2 * s_1 &
         + (-4.999999999999999e-01_rk) * xi(5)**2 * xi(6)**2 * c_1 &
         + (+8.660254037844386e-01_rk) * xi(5)**2 * xi(6)**2 * s_1 &
         + (+9.999999999999994e-01_rk) * xi(3)**2 * xi(4)**2 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(2)**2 * xi(4)**2 * c_1 &
         + (+8.660254037844386e-01_rk) * xi(2)**2 * xi(4)**2 * s_1 &
         + (-4.999999999999994e-01_rk) * xi(2)**2 * xi(3)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(2)**2 * xi(3)**2 * s_1
    phi_y(99) = (-1.000000000000000e+00_rk) * xi(6)**2 * xi(7)**2 * s_1 &
         + (-8.660254037844383e-01_rk) * xi(5)**2 * xi(7)**2 * c_1 &
         + (+5.000000000000001e-01_rk) * xi(5)**2 * xi(7)**2 * s_1 &
         + (+8.660254037844382e-01_rk) * xi(5)**2 * xi(6)**2 * c_1 &
         + (+4.999999999999999e-01_rk) * xi(5)**2 * xi(6)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(3)**2 * xi(4)**2 * s_1 &
         + (-8.660254037844382e-01_rk) * xi(2)**2 * xi(4)**2 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(2)**2 * xi(4)**2 * s_1 &
         + (+8.660254037844383e-01_rk) * xi(2)**2 * xi(3)**2 * c_1 &
         + (-5.000000000000001e-01_rk) * xi(2)**2 * xi(3)**2 * s_1
    phi_x(100) = (-4.999999999999998e-01_rk) * xi(4)**2 * xi(7)**2 * c_1 &
         + (-4.999999999999998e-01_rk) * xi(4)**2 * xi(6)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(4)**2 * xi(6)**2 * s_1 &
         + (+2.500000000000000e-01_rk) * xi(4)**2 * xi(5)**2 * c_1 &
         + (-4.330127018922191e-01_rk) * xi(4)**2 * xi(5)**2 * s_1 &
         + (-5.000000000000000e-01_rk) * xi(3)**2 * xi(7)**2 * c_1 &
         + (+8.660254037844383e-01_rk) * xi(3)**2 * xi(7)**2 * s_1 &
         + (-4.999999999999998e-01_rk) * xi(3)**2 * xi(6)**2 * c_1 &
         + (+2.499999999999999e-01_rk) * xi(3)**2 * xi(5)**2 * c_1 &
         + (+4.330127018922192e-01_rk) * xi(3)**2 * xi(5)**2 * s_1 &
         + (+2.499999999999999e-01_rk) * xi(2)**2 * xi(7)**2 * c_1 &
         + (+4.330127018922192e-01_rk) * xi(2)**2 * xi(7)**2 * s_1 &
         + (+2.500000000000000e-01_rk) * xi(2)**2 * xi(6)**2 * c_1 &
         + (-4.330127018922191e-01_rk) * xi(2)**2 * xi(6)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**2 * xi(5)**2 * c_1
    phi_y(100) = (+8.660254037844385e-01_rk) * xi(4)**2 * xi(7)**2 * c_1 &
         + (+4.330127018922193e-01_rk) * xi(4)**2 * xi(5)**2 * c_1 &
         + (-7.499999999999997e-01_rk) * xi(4)**2 * xi(5)**2 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3)**2 * xi(6)**2 * c_1 &
         + (-4.330127018922191e-01_rk) * xi(3)**2 * xi(5)**2 * c_1 &
         + (-7.499999999999998e-01_rk) * xi(3)**2 * xi(5)**2 * s_1 &
         + (+4.330127018922191e-01_rk) * xi(2)**2 * xi(7)**2 * c_1 &
         + (+7.499999999999998e-01_rk) * xi(2)**2 * xi(7)**2 * s_1 &
         + (-4.330127018922193e-01_rk) * xi(2)**2 * xi(6)**2 * c_1 &
         + (+7.499999999999997e-01_rk) * xi(2)**2 * xi(6)**2 * s_1
    phi_x(101) = (+9.999999999999991e-01_rk) * xi(12)**2 * xi(13)**2 * c_1 &
         + (-4.999999999999995e-01_rk) * xi(11)**2 * xi(13)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(11)**2 * xi(13)**2 * s_1 &
         + (-4.999999999999999e-01_rk) * xi(11)**2 * xi(12)**2 * c_1 &
         + (+8.660254037844386e-01_rk) * xi(11)**2 * xi(12)**2 * s_1 &
         + (+9.999999999999994e-01_rk) * xi(9)**2 * xi(10)**2 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(8)**2 * xi(10)**2 * c_1 &
         + (+8.660254037844386e-01_rk) * xi(8)**2 * xi(10)**2 * s_1 &
         + (-4.999999999999994e-01_rk) * xi(8)**2 * xi(9)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(8)**2 * xi(9)**2 * s_1
    phi_y(101) = (-1.000000000000000e+00_rk) * xi(12)**2 * xi(13)**2 * s_1 &
         + (-8.660254037844383e-01_rk) * xi(11)**2 * xi(13)**2 * c_1 &
         + (+5.000000000000001e-01_rk) * xi(11)**2 * xi(13)**2 * s_1 &
         + (+8.660254037844382e-01_rk) * xi(11)**2 * xi(12)**2 * c_1 &
         + (+4.999999999999999e-01_rk) * xi(11)**2 * xi(12)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(9)**2 * xi(10)**2 * s_1 &
         + (-8.660254037844382e-01_rk) * xi(8)**2 * xi(10)**2 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(8)**2 * xi(10)**2 * s_1 &
         + (+8.660254037844383e-01_rk) * xi(8)**2 * xi(9)**2 * c_1 &
         + (-5.000000000000001e-01_rk) * xi(8)**2 * xi(9)**2 * s_1
    phi_x(102) = (-4.999999999999998e-01_rk) * xi(10)**2 * xi(13)**2 * c_1 &
         + (-4.999999999999998e-01_rk) * xi(10)**2 * xi(12)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(10)**2 * xi(12)**2 * s_1 &
         + (+2.500000000000000e-01_rk) * xi(10)**2 * xi(11)**2 * c_1 &
         + (-4.330127018922191e-01_rk) * xi(10)**2 * xi(11)**2 * s_1 &
         + (-5.000000000000000e-01_rk) * xi(9)**2 * xi(13)**2 * c_1 &
         + (+8.660254037844383e-01_rk) * xi(9)**2 * xi(13)**2 * s_1 &
         + (-4.999999999999998e-01_rk) * xi(9)**2 * xi(12)**2 * c_1 &
         + (+2.499999999999999e-01_rk) * xi(9)**2 * xi(11)**2 * c_1 &
         + (+4.330127018922192e-01_rk) * xi(9)**2 * xi(11)**2 * s_1 &
         + (+2.499999999999999e-01_rk) * xi(8)**2 * xi(13)**2 * c_1 &
         + (+4.330127018922192e-01_rk) * xi(8)**2 * xi(13)**2 * s_1 &
         + (+2.500000000000000e-01_rk) * xi(8)**2 * xi(12)**2 * c_1 &
         + (-4.330127018922191e-01_rk) * xi(8)**2 * xi(12)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(8)**2 * xi(11)**2 * c_1
    phi_y(102) = (+8.660254037844385e-01_rk) * xi(10)**2 * xi(13)**2 * c_1 &
         + (+4.330127018922193e-01_rk) * xi(10)**2 * xi(11)**2 * c_1 &
         + (-7.499999999999997e-01_rk) * xi(10)**2 * xi(11)**2 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(9)**2 * xi(12)**2 * c_1 &
         + (-4.330127018922191e-01_rk) * xi(9)**2 * xi(11)**2 * c_1 &
         + (-7.499999999999998e-01_rk) * xi(9)**2 * xi(11)**2 * s_1 &
         + (+4.330127018922191e-01_rk) * xi(8)**2 * xi(13)**2 * c_1 &
         + (+7.499999999999998e-01_rk) * xi(8)**2 * xi(13)**2 * s_1 &
         + (-4.330127018922193e-01_rk) * xi(8)**2 * xi(12)**2 * c_1 &
         + (+7.499999999999997e-01_rk) * xi(8)**2 * xi(12)**2 * s_1
    phi_x(103) = (-5.000000000000000e-01_rk) * xi(7)**2 * xi(13)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(7)**2 * xi(13)**2 * s_1 &
         + (-4.999999999999997e-01_rk) * xi(6)**2 * xi(12)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(6)**2 * xi(12)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * xi(11)**2 * c_1 &
         + (-4.999999999999997e-01_rk) * xi(4)**2 * xi(10)**2 * c_1 &
         + (-8.660254037844386e-01_rk) * xi(4)**2 * xi(10)**2 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(3)**2 * xi(9)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(3)**2 * xi(9)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**2 * xi(8)**2 * c_1
    phi_y(103) = (+8.660254037844385e-01_rk) * xi(7)**2 * xi(13)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(7)**2 * xi(13)**2 * s_1 &
         + (-8.660254037844386e-01_rk) * xi(6)**2 * xi(12)**2 * c_1 &
         + (+4.999999999999999e-01_rk) * xi(6)**2 * xi(12)**2 * s_1 &
         + (-1.000000000000000e+00_rk) * xi(5)**2 * xi(11)**2 * s_1 &
         + (+8.660254037844386e-01_rk) * xi(4)**2 * xi(10)**2 * c_1 &
         + (-4.999999999999999e-01_rk) * xi(4)**2 * xi(10)**2 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3)**2 * xi(9)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(3)**2 * xi(9)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**2 * xi(8)**2 * s_1
    phi_x(104) = (-4.999999999999995e-01_rk) * xi(7)**2 * xi(12)**2 * c_1 &
         + (+8.660254037844387e-01_rk) * xi(7)**2 * xi(12)**2 * s_1 &
         + (-5.000000000000006e-01_rk) * xi(7)**2 * xi(11)**2 * c_1 &
         + (+8.660254037844382e-01_rk) * xi(7)**2 * xi(11)**2 * s_1 &
         + (-4.999999999999993e-01_rk) * xi(6)**2 * xi(13)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(6)**2 * xi(13)**2 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(6)**2 * xi(11)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(6)**2 * xi(11)**2 * s_1 &
         + (+9.999999999999999e-01_rk) * xi(5)**2 * xi(13)**2 * c_1 &
         + (+1.000000000000000e+00_rk) * xi(5)**2 * xi(12)**2 * c_1 &
         + (-4.999999999999993e-01_rk) * xi(4)**2 * xi(9)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(4)**2 * xi(9)**2 * s_1 &
         + (-5.000000000000003e-01_rk) * xi(4)**2 * xi(8)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(4)**2 * xi(8)**2 * s_1 &
         + (-4.999999999999997e-01_rk) * xi(3)**2 * xi(10)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(3)**2 * xi(10)**2 * s_1 &
         + (-5.000000000000004e-01_rk) * xi(3)**2 * xi(8)**2 * c_1 &
         + (+8.660254037844382e-01_rk) * xi(3)**2 * xi(8)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**2 * xi(10)**2 * c_1 &
         + (+9.999999999999999e-01_rk) * xi(2)**2 * xi(9)**2 * c_1
    phi_y(104) = (+8.660254037844385e-01_rk) * xi(7)**2 * xi(12)**2 * c_1 &
         + (+4.999999999999998e-01_rk) * xi(7)**2 * xi(12)**2 * s_1 &
         + (+8.660254037844382e-01_rk) * xi(7)**2 * xi(11)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(7)**2 * xi(11)**2 * s_1 &
         + (-8.660254037844388e-01_rk) * xi(6)**2 * xi(13)**2 * c_1 &
         + (+4.999999999999997e-01_rk) * xi(6)**2 * xi(13)**2 * s_1 &
         + (-8.660254037844382e-01_rk) * xi(6)**2 * xi(11)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(6)**2 * xi(11)**2 * s_1 &
         + (-9.999999999999999e-01_rk) * xi(5)**2 * xi(13)**2 * s_1 &
         + (-9.999999999999997e-01_rk) * xi(5)**2 * xi(12)**2 * s_1 &
         + (+8.660254037844388e-01_rk) * xi(4)**2 * xi(9)**2 * c_1 &
         + (-4.999999999999997e-01_rk) * xi(4)**2 * xi(9)**2 * s_1 &
         + (+8.660254037844382e-01_rk) * xi(4)**2 * xi(8)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(4)**2 * xi(8)**2 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3)**2 * xi(10)**2 * c_1 &
         + (-4.999999999999998e-01_rk) * xi(3)**2 * xi(10)**2 * s_1 &
         + (-8.660254037844382e-01_rk) * xi(3)**2 * xi(8)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(3)**2 * xi(8)**2 * s_1 &
         + (+9.999999999999997e-01_rk) * xi(2)**2 * xi(10)**2 * s_1 &
         + (+9.999999999999999e-01_rk) * xi(2)**2 * xi(9)**2 * s_1
    phi_x(105) = (-8.660254037844385e-01_rk) * xi(7)**2 * xi(12)**2 * c_1 &
         + (-4.999999999999998e-01_rk) * xi(7)**2 * xi(12)**2 * s_1 &
         + (+8.660254037844382e-01_rk) * xi(7)**2 * xi(11)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(7)**2 * xi(11)**2 * s_1 &
         + (-8.660254037844388e-01_rk) * xi(6)**2 * xi(13)**2 * c_1 &
         + (+4.999999999999997e-01_rk) * xi(6)**2 * xi(13)**2 * s_1 &
         + (+8.660254037844382e-01_rk) * xi(6)**2 * xi(11)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(6)**2 * xi(11)**2 * s_1 &
         + (+9.999999999999999e-01_rk) * xi(5)**2 * xi(13)**2 * s_1 &
         + (-9.999999999999997e-01_rk) * xi(5)**2 * xi(12)**2 * s_1 &
         + (-8.660254037844388e-01_rk) * xi(4)**2 * xi(9)**2 * c_1 &
         + (+4.999999999999997e-01_rk) * xi(4)**2 * xi(9)**2 * s_1 &
         + (+8.660254037844382e-01_rk) * xi(4)**2 * xi(8)**2 * c_1 &
         + (-5.000000000000000e-01_rk) * xi(4)**2 * xi(8)**2 * s_1 &
         + (-8.660254037844385e-01_rk) * xi(3)**2 * xi(10)**2 * c_1 &
         + (-4.999999999999998e-01_rk) * xi(3)**2 * xi(10)**2 * s_1 &
         + (+8.660254037844382e-01_rk) * xi(3)**2 * xi(8)**2 * c_1 &
         + (+5.000000000000000e-01_rk) * xi(3)**2 * xi(8)**2 * s_1 &
         + (-9.999999999999997e-01_rk) * xi(2)**2 * xi(10)**2 * s_1 &
         + (+9.999999999999999e-01_rk) * xi(2)**2 * xi(9)**2 * s_1
    phi_y(105) = (-4.999999999999995e-01_rk) * xi(7)**2 * xi(12)**2 * c_1 &
         + (+8.660254037844387e-01_rk) * xi(7)**2 * xi(12)**2 * s_1 &
         + (+5.000000000000006e-01_rk) * xi(7)**2 * xi(11)**2 * c_1 &
         + (-8.660254037844382e-01_rk) * xi(7)**2 * xi(11)**2 * s_1 &
         + (+4.999999999999993e-01_rk) * xi(6)**2 * xi(13)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(6)**2 * xi(13)**2 * s_1 &
         + (-5.000000000000002e-01_rk) * xi(6)**2 * xi(11)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(6)**2 * xi(11)**2 * s_1 &
         + (+9.999999999999999e-01_rk) * xi(5)**2 * xi(13)**2 * c_1 &
         + (-1.000000000000000e+00_rk) * xi(5)**2 * xi(12)**2 * c_1 &
         + (-4.999999999999993e-01_rk) * xi(4)**2 * xi(9)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(4)**2 * xi(9)**2 * s_1 &
         + (+5.000000000000003e-01_rk) * xi(4)**2 * xi(8)**2 * c_1 &
         + (+8.660254037844385e-01_rk) * xi(4)**2 * xi(8)**2 * s_1 &
         + (+4.999999999999997e-01_rk) * xi(3)**2 * xi(10)**2 * c_1 &
         + (-8.660254037844385e-01_rk) * xi(3)**2 * xi(10)**2 * s_1 &
         + (-5.000000000000004e-01_rk) * xi(3)**2 * xi(8)**2 * c_1 &
         + (+8.660254037844382e-01_rk) * xi(3)**2 * xi(8)**2 * s_1 &
         + (+1.000000000000000e+00_rk) * xi(2)**2 * xi(10)**2 * c_1 &
         + (-9.999999999999999e-01_rk) * xi(2)**2 * xi(9)**2 * c_1
    phi_x(106) = (+1.000000000000000e+00_rk) * xi(15)**2 * xi(17)**2 * c_1 &
         + (-9.999999999999998e-01_rk) * xi(15)**2 * xi(16) * xi(17) * s_1 &
         + (+9.999999999999998e-01_rk) * xi(14) * xi(15) * xi(17)**2 * s_1 &
         + (+9.999999999999996e-01_rk) * xi(14) * xi(15) * xi(16)**2 * s_1 &
         + (-9.999999999999994e-01_rk) * xi(14)**2 * xi(16) * xi(17) * s_1 &
         + (-1.000000000000000e+00_rk) * xi(14)**2 * xi(16)**2 * c_1
    phi_y(106) = (-1.000000000000000e+00_rk) * xi(15)**2 * xi(16) * xi(17) * c_1 &
         + (+9.999999999999992e-01_rk) * xi(15)**2 * xi(16)**2 * s_1 &
         + (-1.000000000000000e+00_rk) * xi(14) * xi(15) * xi(17)**2 * c_1 &
         + (-9.999999999999998e-01_rk) * xi(14) * xi(15) * xi(16)**2 * c_1 &
         + (-9.999999999999993e-01_rk) * xi(14)**2 * xi(17)**2 * s_1 &
         + (-9.999999999999998e-01_rk) * xi(14)**2 * xi(16) * xi(17) * c_1
    phi_x(107) = (-4.999999999999998e-01_rk) * xi(4)**2 * xi(7)**2 * c_3 &
         + (+1.000000000000000e+00_rk) * xi(4)**2 * xi(6)**2 * c_3 &
         + (-4.999999999999998e-01_rk) * xi(4)**2 * xi(5)**2 * c_3 &
         + (+1.000000000000000e+00_rk) * xi(3)**2 * xi(7)**2 * c_3 &
         + (-4.999999999999998e-01_rk) * xi(3)**2 * xi(6)**2 * c_3 &
         + (-4.999999999999998e-01_rk) * xi(3)**2 * xi(5)**2 * c_3 &
         + (-4.999999999999998e-01_rk) * xi(2)**2 * xi(7)**2 * c_3 &
         + (-4.999999999999998e-01_rk) * xi(2)**2 * xi(6)**2 * c_3 &
         + (+1.000000000000000e+00_rk) * xi(2)**2 * xi(5)**2 * c_3
    phi_y(107) = (+8.660254037844385e-01_rk) * xi(4)**2 * xi(7)**2 * c_3 &
         + (-8.660254037844385e-01_rk) * xi(4)**2 * xi(5)**2 * c_3 &
         + (-8.660254037844385e-01_rk) * xi(3)**2 * xi(6)**2 * c_3 &
         + (+8.660254037844385e-01_rk) * xi(3)**2 * xi(5)**2 * c_3 &
         + (-8.660254037844385e-01_rk) * xi(2)**2 * xi(7)**2 * c_3 &
         + (+8.660254037844385e-01_rk) * xi(2)**2 * xi(6)**2 * c_3
    phi_x(108) = (-4.999999999999998e-01_rk) * xi(10)**2 * xi(13)**2 * c_3 &
         + (+1.000000000000000e+00_rk) * xi(10)**2 * xi(12)**2 * c_3 &
         + (-4.999999999999998e-01_rk) * xi(10)**2 * xi(11)**2 * c_3 &
         + (+1.000000000000000e+00_rk) * xi(9)**2 * xi(13)**2 * c_3 &
         + (-4.999999999999998e-01_rk) * xi(9)**2 * xi(12)**2 * c_3 &
         + (-4.999999999999998e-01_rk) * xi(9)**2 * xi(11)**2 * c_3 &
         + (-4.999999999999998e-01_rk) * xi(8)**2 * xi(13)**2 * c_3 &
         + (-4.999999999999998e-01_rk) * xi(8)**2 * xi(12)**2 * c_3 &
         + (+1.000000000000000e+00_rk) * xi(8)**2 * xi(11)**2 * c_3
    phi_y(108) = (+8.660254037844385e-01_rk) * xi(10)**2 * xi(13)**2 * c_3 &
         + (-8.660254037844385e-01_rk) * xi(10)**2 * xi(11)**2 * c_3 &
         + (-8.660254037844385e-01_rk) * xi(9)**2 * xi(12)**2 * c_3 &
         + (+8.660254037844385e-01_rk) * xi(9)**2 * xi(11)**2 * c_3 &
         + (-8.660254037844385e-01_rk) * xi(8)**2 * xi(13)**2 * c_3 &
         + (+8.660254037844385e-01_rk) * xi(8)**2 * xi(12)**2 * c_3
    phi_x(109) = (-2.499999999999999e-01_rk) * xi(15)**2 * xi(17)**2 * c_3 &
         + (+2.499999999999999e-01_rk) * xi(15)**2 * xi(16)**2 * c_3 &
         + (+1.000000000000000e+00_rk) * xi(14) * xi(15) * xi(16) * xi(17) * c_3 &
         + (+2.499999999999999e-01_rk) * xi(14)**2 * xi(17)**2 * c_3 &
         + (-2.500000000000002e-01_rk) * xi(14)**2 * xi(16)**2 * c_3
    phi_y(109) = (-5.000000000000001e-01_rk) * xi(15)**2 * xi(16) * xi(17) * c_3 &
         + (-5.000000000000001e-01_rk) * xi(14) * xi(15) * xi(17)**2 * c_3 &
         + (+4.999999999999997e-01_rk) * xi(14) * xi(15) * xi(16)**2 * c_3 &
         + (+4.999999999999997e-01_rk) * xi(14)**2 * xi(16) * xi(17) * c_3
  end subroutine eval_c2h6_dms_basis

  subroutine eval_c2h6_dms(xi, params_z, params_xy, mu)
    real(rk), intent(in)  :: xi(18)
    real(rk), intent(in)  :: params_z(N_PARAMS_Z)
    real(rk), intent(in)  :: params_xy(N_PARAMS_XY)
    real(rk), intent(out) :: mu(3)
    real(rk) :: phi_z(N_PARAMS_Z), phi_x(N_PARAMS_XY), phi_y(N_PARAMS_XY)
    call eval_c2h6_dms_basis(xi, phi_z, phi_x, phi_y)
    mu(1) = dot_product(params_xy, phi_x)
    mu(2) = dot_product(params_xy, phi_y)
    mu(3) = dot_product(params_z,  phi_z)
  end subroutine eval_c2h6_dms


end module pot_c2h6
