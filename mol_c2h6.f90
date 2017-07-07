! This unit defines all specific routines for a eight-atomic molecule of ethane-type

module mol_c2h6
  use accuracy
  use moltype

  implicit none

  public ML_coordinate_transform_C2H6, ML_b0_C2H6, ML_symmetry_transformation_C2H6, ML_rotsymmetry_C2H6
  private

  integer(ik), parameter :: verbose = 3 ! Verbosity level

  !--------------------------
  !      ZMAT_4BETA_1TAU
  !--------------------------
  !
  !ZMAT
  !C   0  0  0  0  12.00000000
  !C   1  0  0  0  12.00000000
  !H   1  2  0  0   1.00782505
  !H   1  2  3  2   1.00782505
  !H   1  2  3 -2   1.00782505
  !H   2  1  3  2   1.00782505
  !H   2  1  5  2   1.00782505
  !H   2  1  5 -2   1.00782505
  !end
  ! 1  r_12        0         1.52579576
  ! 2  r_31        0         1.09074923
  ! 3  r_41        0         1.09074923
  ! 4  r_51        0         1.09074923
  ! 5  r_62        0         1.09074923
  ! 6  r_72        0         1.09074923
  ! 7  r_82        0         1.09074923
  ! 8   alpha_312   0         111.0 DEG
  ! 9   alpha_412   0         111.0 DEG
  ! 10  alpha_512   0         111.0 DEG
  ! 11  alpha_621   0         111.0 DEG
  ! 12  alpha_721   0         111.0 DEG
  ! 13  alpha_821   0         111.0 DEG
  ! 14  thet12      0         120.00 DEG
  ! 15  thet13      0         120.00 DEG
  ! 16  thet45      0         120.00 DEG
  ! 17  thet64      0         120.00 DEG
  ! 18  t14         0         180.00 DEG
  !end




  contains


  function ML_coordinate_transform_C2H6(src,ndst,direct) result (dst)
    !
    ! Transformtation from Z-matrix to TROVE coords
    !
    real(ark),intent(in)  :: src(:)
    integer(ik),intent(in) :: ndst
    logical,intent(in):: direct
    !
    real(ark),dimension(ndst) :: dst
    integer(ik) :: nsrc
    real(ark) :: tau4213,tau5124,tau6213,tau5126
    real(ark) :: tau1, tau2, b1_zmat, b2_zmat, tau1_zmat, dtau, db1, db2
    !
    real(ark) :: tau14,tau24,tau25,tau35,tau36,theta12,theta23,theta13,theta56,theta45,theta46,&
                 S1,S2,S3,S4,S5,S6,taubar,tau

    if (verbose>=5) write(out, '(/a)') 'ML_coordinate_transform_C2H6/start'
    !
    nsrc = size(src)
    !
    select case(trim(molec%coords_transform))
      !
    case default
      !
      write(out, '(/a,1x,a,1x,a)') &
      'ML_coordinate_transform_C2H6 error: coordinate type =', trim(molec%coords_transform), 'is unknown'
      stop 'ML_coordinate_transform_C2H6 error: bad coordinate type'
      !
    case('ZMAT_4BETA_1TAU','C2H6_4BETA_1TAU')
      ! ORDER CHANGED HERE AS WANT TBAR TO BE 18TH COORDINATE BUT NEEDS
      !  TO BE 16TH COORDINATE IN Z-MAT
      !
      if (direct) then ! transform from Z-matrix coords to TROVE coords

        dst(1:13) = src(1:13)-molec%local_eq(1:13)

        dst(14)  = (3.0_ark*src(14) - 2.0_ark*pi)/sqrt(6.0_ark)
        dst(15)  = ( 2.0_ark*src(15) - 2.0_ark*pi + src(14))/sqrt(2.0_ark)
        dst(16) = (3.0_ark*src(18) - 2.0_ark*pi)/sqrt(6.0_ark)
        dst(17) = ( 2.0_ark*src(17) - 2.0_ark*pi + src(18))/sqrt(2.0_ark)
        dst(18) = ( ( 3.0_ark*src(16) + src(14) - src(15) + src(17) - src(18))/3.0_ark ) - molec%local_eq(16)
        !
      else !  transform from TROVE coords to Z-matrix coords
        !
        dst(1:13) = src(1:13)+molec%local_eq(1:13)

        dst(14) = ((sqrt(6.0_ark))*src(14) + 2.0_ark*pi)/3.0_ark
        dst(15) = ( (sqrt(2.0_ark))*src(15) - dst(14) + 2.0_ark*pi)/2.0_ark
        dst(18) = ((sqrt(6.0_ark))*src(16) + 2.0_ark*pi)/3.0_ark
        dst(17) = ( (sqrt(2.0_ark))*src(17) - dst(18) + 2.0_ark*pi)/2.0_ark
        dst(16) = ( (3.0_ark*src(18) - dst(14) + dst(15) - dst(17) + dst(18) )/3.0_ark  ) + molec%local_eq(16)
        !
      endif
      !
    case('R-R16-BETA16-THETA-TAU')
      ! ORDER CHANGED HERE AS WANT TBAR TO BE 18TH COORDINATE BUT NEEDS
      !  TO BE 16TH COORDINATE IN Z-MAT
      !
      if (direct) then ! transform from Z-matrix coords to TROVE coords

        dst(1:13) = src(1:13)-molec%local_eq(1:13)
        !
        tau14 = src(14)
        tau24 = src(15)
        tau25 = src(16)
        tau35 = src(17)
        tau36 = src(18)
        !
        theta12 = tau14-tau24
        theta23 = tau25-tau35
        theta13 = 2.0_ark*pi-theta12-theta23
        !
        theta56 = tau36-tau35
        theta45 = tau25-tau24
        theta46 = 2.0_ark*pi-theta56-theta45
        !
        dst(14)  = ( 2.0_ark*theta23 - theta13 - theta12 )/sqrt(6.0_ark)
        dst(15)  = (                   theta13 - theta12 )/sqrt(2.0_ark)
        !
        dst(16)  = ( 2.0_ark*theta56 - theta46 - theta45 )/sqrt(6.0_ark)
        dst(17)  = (                   theta46 - theta45 )/sqrt(2.0_ark)
        !
        dst(18)  = ( tau14+tau25+tau36 )/sqrt(3.0_ark)-3.0_ark*pi/sqrt(3.0_ark)
        !
      else !  transform from TROVE coords to Z-matrix coords
        !
        dst(1:13) = src(1:13)+molec%local_eq(1:13)
        !
        S1 = src(14)
        S2 = src(15)
        S3 = 2.0_ark*pi
        S4 = src(16)
        S5 = src(17)
        S6 = 2.0_ark*pi
        taubar = src(18)+3.0_ark*pi/sqrt(3.0_ark)
        Tau = taubar*sqrt(3.0_ark)
        !
        theta23 = 1.0_ark/3.0_ark*(sqrt(6.0_ark)*S1+S3)
        theta13 = sqrt(6.0_ark)/18.0_ark*( 3.0_ark*sqrt(3.0_ark)*S2+sqrt(6.0_ark)*S3-3.0_ark*S1)
        theta12 = sqrt(6.0_ark)/18.0_ark*(-3.0_ark*sqrt(3.0_ark)*S2+sqrt(6.0_ark)*S3-3.0_ark*S1)
        !
        theta56 = 1.0_ark/3.0_ark*(sqrt(6.0_ark)*S4+S6)
        theta46 = sqrt(6.0_ark)/18.0_ark*( 3.0_ark*sqrt(3.0_ark)*S5+sqrt(6.0_ark)*S6-3.0_ark*S4)
        theta45 = sqrt(6.0_ark)/18.0_ark*(-3.0_ark*sqrt(3.0_ark)*S5+sqrt(6.0_ark)*S6-3.0_ark*S4)

        tau36 = 1.0_ark/3.0_ark*(-2.0_ark*theta23+Tau+2.0_ark*theta56+        theta45-        theta12)
        tau35 = 1.0_ark/3.0_ark*(-2.0_ark*theta23+Tau-        theta56+        theta45-        theta12)
        tau14 = 1.0_ark/3.0_ark*(         theta23+Tau-        theta56-2.0_ark*theta45+2.0_ark*theta12)
        tau24 = 1.0_ark/3.0_ark*(         theta23+Tau-        theta56-2.0_ark*theta45-        theta12)
        tau25 = 1.0_ark/3.0_ark*(         theta23+Tau-        theta56+        theta45-        theta12)
        !
        dst(14) = tau14
        dst(15) = tau24
        dst(16) = tau25
        dst(17) = tau35
        dst(18) = tau36
        !
      endif
      !
    end select
    !
    if (verbose>=5) write(out, '(/a)') 'ML_coordinate_transform_C2H6/end'
    !
  end function ML_coordinate_transform_C2H6



  subroutine ML_b0_C2H6(Npoints,Natoms,b0,rho_i,rho_ref,rho_borders)
    !
    ! Reference structure
    !
    integer(ik),intent(in) :: Npoints, Natoms
    real(ark),intent(out) :: b0(Natoms,3,0:Npoints)
    real(ark),intent(inout),optional :: rho_i(0:Npoints)
    real(ark),intent(out),optional :: rho_ref
    real(ark),intent(in),optional :: rho_borders(2)  ! rhomim, rhomax - borders
    !
    real(ark) :: a0(molec%Natoms,3),CM_shift,tau,alpha0,alpha,theta,r,r12
    integer(ik) :: i, n, iatom, ix
    !
    if (verbose>=5) write(out, '(/a)') 'ML_b0_C2H6/start'
    !
    r12 = molec%req(1)
    r   = molec%req(2)
    alpha = molec%alphaeq(1)
    theta = 2.0_ark*pi/3.0_ark
    !
    if (any(molec%req(2:7)/=r)) then
      write(out,"('ML_b0_C2H6 error: eq-m r2-r7 are not all the same:',6f12.5)") molec%req(2:7)
      stop 'ML_b0_C2H6 error: eq-m r2-r7 are not all the same'
    endif
    !
    if (any(molec%alphaeq(1:6)/=alpha)) then
      write(out,"('ML_b0_C2H6 error: eq-m alphas are not all the same:',6f12.5)") molec%alphaeq(1:6)
      stop 'ML_b0_C2H6 error: eq-m alphas are not all the same'
    endif
    !
    !
    a0 = 0
    !
    a0(3,1) = r*sin(alpha)
    a0(3,2) = 0
    a0(3,3) = r*cos(alpha)
    !
    a0(4,1) =-r*sin(alpha)
    a0(4,2) = 0
    a0(4,3) =-r*cos(alpha)
    !
    a0(5,1) = r*sin(alpha)*cos(theta)
    a0(5,2) =-r*sin(alpha)*sin(theta)
    a0(5,3) = r*cos(alpha)
    !
    a0(6,1) =-r*sin(alpha)*cos(theta)
    a0(6,2) = r*sin(alpha)*sin(theta)
    a0(6,3) =-r*cos(alpha)
    !
    a0(7,1) = r*sin(alpha)*cos(theta)
    a0(7,2) = r*sin(alpha)*sin(theta)
    a0(7,3) = r*cos(alpha)
    !
    a0(8,1) =-r*sin(alpha)*cos(theta)
    a0(8,2) =-r*sin(alpha)*sin(theta)
    a0(8,3) =-r*cos(alpha)
    !
    a0(1:7:2,3) =a0(1:7:2,3)-r12*0.5_ark
    a0(2:8:2,3) =a0(2:8:2,3)+r12*0.5_ark
    !
    do ix=1, 3
      CM_shift = sum(a0(:,ix)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
      a0(:,ix) = a0(:,ix) - CM_shift
    enddo
    !
    if (verbose>=3) then
      write(out, '(1x,a,1x,3(1x,es16.8))') 'C', a0(1,1:3)
      write(out, '(1x,a,1x,3(1x,es16.8))') 'C', a0(2,1:3)
      do iatom=3, Natoms
        write(out, '(1x,a,1x,3(1x,es16.8))') 'H', a0(iatom,1:3)
      enddo
    endif
    !
    b0(:,:,0) = a0(:,:)
    !
    if (Npoints/=0) then
      rho_ref = 0
      write(out,"('ML_b0_C2H6 error: the non-rigid part has not been implemebted yet')") 
      stop 'ML_b0_C2H6 error: the non-rigid part has not been implemebted yet'
    end if

    !
    if (verbose>=5) write(out, '(/a)') 'ML_b0_C2H6/end'
    !
  end subroutine ML_b0_C2H6



  subroutine ML_symmetry_transformation_C2H6(ioper, nmodes, src, dst)
    !
    ! Symmetry transformation rules of coordinates
    !
    integer(ik), intent(in)  ::  ioper
    integer(ik), intent(in)  ::  nmodes
    real(ark), intent(in)    ::  src(1:nmodes)
    real(ark), intent(out)   ::  dst(1:nmodes)
    real(ark) :: a,b
    !
    a = 0.5_ark
    b = 0.5_ark*sqrt(3.0_ark)
    if (verbose>=5) write(out, '(/a)') 'ML_symmetry_transformation_C2H6/start'
    !
    select case(trim(molec%coords_transform))
      !
    case default
      !
      write(out, '(/a,1x,a,1x,a)') &
      'ML_symmetry_transformation_C2H6 error: coordinate type =', trim(molec%coords_transform), 'is unknown'
      stop 'ML_symmetry_transformation_C2H6 error: bad coordinate type'
      !
    case('ZMAT_4BETA_1TAU','C2H6_4BETA_1TAU')
      !
      select case(trim(molec%symmetry))
        !
      case default
        !
        write(out, '(/a,1x,a,1x,a)') &
        'ML_symmetry_transformation_C2H6 error: symmetry =', trim(molec%symmetry), 'is unknown'
        stop
        !
      case('D3D(M)')
        !
        select case(ioper)
          !
        case default
          !
          write(out, '(/a,1x,i3,1x,a)') &
          'ML_symmetry_transformation_C2H6 error: symmetry operation ', ioper, 'is unknown'
          stop
          !
        case (1) ! E
          !
          dst(1:18) = src(1:18)
          !
        case (2) ! !C3+/(123)(465)
          !
          dst(1) = src(1)
          dst(2) = src(4)
          dst(3) = src(2)
          dst(4) = src(3)
          dst(5) = src(6)
          dst(6) = src(7)
          dst(7) = src(5)
          dst(8) = src(10)
          dst(9) = src(8)
          dst(10) = src(9)
          dst(11) = src(12)
          dst(12) = src(13)
          dst(13) = src(11)
          dst(14) = -a*src(14) + b*src(15)
          dst(15) = -b*src(14) - a*src(15) 
          dst(16) = -a*src(16) + b*src(17)
          dst(17) = -b*src(16) - a*src(17)
          dst(18) = src(18)
          !
        case (3) !C3-/(123)(465)
          !
          dst(1) = src(1)
          dst(2) = src(3)
          dst(3) = src(4)
          dst(4) = src(2)
          dst(5) = src(7)
          dst(6) = src(5)
          dst(7) = src(6)
          dst(8) = src(9)
          dst(9) = src(10)
          dst(10) = src(8)
          dst(11) = src(13)
          dst(12) = src(11)
          dst(13) = src(12)
          dst(14) = -a*src(14) - b*src(15)
          dst(15) =  b*src(14) - a*src(15) 
          dst(16) = -a*src(16) - b*src(17)
          dst(17) =  b*src(16) - a*src(17)
          dst(18) = src(18)
          !
        case (4) ! C2/(16)(24)(35)(78)
          !
          dst(1) = src(1)
          dst(2) = src(7)
          dst(3) = src(5)
          dst(4) = src(6)
          dst(5) = src(3)
          dst(6) = src(4)
          dst(7) = src(2)
          dst(8) = src(13)
          dst(9) = src(11)
          dst(10) = src(12)
          dst(11) = src(9)
          dst(12) = src(10)
          dst(13) = src(8)
          dst(14) = src(16)
          dst(15) = -src(17) 
          dst(16) = src(14)
          dst(17) = -src(15)
          dst(18) = src(18)
          !
        case (5) !C2'/(16)(24)(35)(78)
          !
          dst(1) = src(1)
          dst(2) = src(6)
          dst(3) = src(7)
          dst(4) = src(5)
          dst(5) = src(4)
          dst(6) = src(2)
          dst(7) = src(3)
          dst(8) = src(12)
          dst(9) = src(13)
          dst(10) = src(11)
          dst(11) = src(10)
          dst(12) = src(8)
          dst(13) = src(9)
          dst(14) = -a*src(16) - b*src(17)
          dst(15) =  b*src(16) + a*src(17) 
          dst(16) = -a*src(14) - b*src(16)
          dst(17) = -b*src(14) + a*src(15)
          dst(18) = src(18)
          !
        case (6) ! C2''(16)(24)(35)(78)
          !
          dst(1) = src(1)
          dst(2) = src(5)
          dst(3) = src(6)
          dst(4) = src(7)
          dst(5) = src(2)
          dst(6) = src(3)
          dst(7) = src(4)
          dst(8) = src(11)
          dst(9) = src(12)
          dst(10) = src(13)
          dst(11) = src(8)
          dst(12) = src(9)
          dst(13) = src(10)
          dst(14) = -a*src(16) + b*src(17)
          dst(15) =  b*src(16) + a*src(17) 
          dst(16) = -a*src(14) + b*src(16)
          dst(17) = b*src(14) + a*src(15)
          dst(18) = src(18)
          !
        case (7) ! i/(14)(26)(35)(78)*
          !
          dst(1) = src(1)
          dst(2) = src(5)
          dst(3) = src(7)
          dst(4) = src(6)
          dst(5) = src(2)
          dst(6) = src(4)
          dst(7) = src(3)
          dst(8) = src(11)
          dst(9) = src(13)
          dst(10) = src(12)
          dst(11) = src(8)
          dst(12) = src(10)
          dst(13) = src(9)
          dst(14) = src(16)
          dst(15) = src(17) 
          dst(16) = src(14)
          dst(17) = src(15)
          dst(18) = -src(18)
          !
        case (8) ! S6/(163425)(78)*
          !
          dst(1) = src(1)
          dst(2) = src(6)
          dst(3) = src(5)
          dst(4) = src(7)
          dst(5) = src(4)
          dst(6) = src(3)
          dst(7) = src(2)
          dst(8) = src(12)
          dst(9) = src(11)
          dst(10) = src(13)
          dst(11) = src(10)
          dst(12) = src(9)
          dst(13) = src(8)
          dst(14) = -a*src(16) + b*src(17)
          dst(15) =  -b*src(16) - a*src(17) 
          dst(16) = -a*src(14) + b*src(16)
          dst(17) = -b*src(14) - a*src(15)
          dst(18) = -src(18)
          !
        case (9) !S6'/(14)(26)(35)(78)*
          !
          dst(1) = src(1)
          dst(2) = src(7)
          dst(3) = src(6)
          dst(4) = src(5)
          dst(5) = src(3)
          dst(6) = src(2)
          dst(7) = src(4)
          dst(8) = src(13)
          dst(9) = src(12)
          dst(10) = src(11)
          dst(11) = src(9)
          dst(12) = src(8)
          dst(13) = src(10)
          dst(14) = -a*src(16) - b*src(17)
          dst(15) =  b*src(16) - a*src(17) 
          dst(16) = -a*src(14) - b*src(16)
          dst(17) =  b*src(14) - a*src(15)
          dst(18) = -src(18)
          !
        case (10) !sigmad/(12)(46)*
          !
          dst(1) = src(1)
          dst(2) = src(4)
          dst(3) = src(3)
          dst(4) = src(2)
          dst(5) = src(6)
          dst(6) = src(5)
          dst(7) = src(7)
          dst(8) = src(10)
          dst(9) = src(9)
          dst(10) = src(8)
          dst(11) = src(12)
          dst(12) = src(11)
          dst(13) = src(13)
          dst(14) = -a*src(14) - b*src(15)
          dst(15) =  -b*src(14) + a*src(15) 
          dst(16) = -a*src(16) - b*src(17)
          dst(17) =  -b*src(16) + a*src(17)
          dst(18) = -src(18)
          !
        case (11) !sigmad'/(12)(46)*
          !
          dst(1) = src(1)
          dst(2) = src(2)
          dst(3) = src(4)
          dst(4) = src(3)
          dst(5) = src(5)
          dst(6) = src(7)
          dst(7) = src(6)
          dst(8) = src(8)
          dst(9) = src(10)
          dst(10) = src(9)
          dst(11) = src(11)
          dst(12) = src(13)
          dst(13) = src(12)
          dst(14) = -a*src(14) + b*src(15)
          dst(15) =  b*src(14) + a*src(15) 
          dst(16) = -a*src(16) + b*src(17)
          dst(17) =  b*src(16) + a*src(17)
          dst(18) = -src(18)
          !
        case (12) !sigmad''/(12)(46)*
          !
          dst(1) = src(1)
          dst(2) = src(3)
          dst(3) = src(2)
          dst(4) = src(4)
          dst(5) = src(7)
          dst(6) = src(6)
          dst(7) = src(5)
          dst(8) = src(9)
          dst(9) = src(8)
          dst(10) = src(10)
          dst(11) = src(13)
          dst(12) = src(12)
          dst(13) = src(11)
          dst(14) = src(14)
          dst(15) = -src(15) 
          dst(16) =  src(16)
          dst(17) =  -src(17)
          dst(18) = -src(18)
        end select
        !
        !
      end select
      !
    end select
    !
    if (verbose>=5) write(out, '(/a)') 'ML_symmetry_transformation_C2H6/end'
    !
  end subroutine ML_symmetry_transformation_C2H6



  subroutine ML_rotsymmetry_C2H6(J,K,tau,gamma,ideg)
    !
    ! Symmetry transformation rules of TROVE rotational functions
    !
    integer(ik),intent(in) :: J,K,tau
    integer(ik),intent(out) :: gamma,ideg
    !
    if (verbose>=5) write(out, '(/a)') 'ML_rotsymmetry_C2H6/start'
    !
    select case(trim(molec%coords_transform))
      !
      !
    case default
      !
      write(out, '(/a,1x,a,1x,a)') &
      'ML_rotsymmetry_C2H6 error: coordinate type =', trim(molec%coords_transform), 'is unknown'
      stop 'ML_rotsymmetry_C2H6 error: bad coordinate type'
      !
      !
    case('ZMAT_4BETA_1TAU','C2H6_4BETA_1TAU')
      !
      select case(trim(molec%symmetry))
        !
      case default
        !
        write(out, '(/a,1x,a,1x,a)') &
        'ML_rotsymmetry_C2H6 error: symmetry =', trim(molec%symmetry), 'is unknown'
        stop 'ML_rotsymmetry_C2H6 error: bad symmetry type'
        !
      case('D2H(M)')
        !
        gamma = 0
        ideg = 1
        !
       ! if (mod(K+2,2)==0.and.tau==0) gamma = 1 ! Ag
       ! if (mod(K+2,2)==0.and.tau==1) gamma = 3 ! B1g
       ! if (mod(K+2,2)/=0.and.tau==1) gamma = 7 ! B3g
       ! if (mod(K+2,2)/=0.and.tau==0) gamma = 5 ! B2g
        if (mod(K+2,2)==0.and.tau==0) gamma = 1 ! Ag
        if (mod(K+2,2)==0.and.tau==1) gamma = 7 ! B1g
        if (mod(K+2,2)/=0.and.tau==1) gamma = 3 ! B3g
        if (mod(K+2,2)/=0.and.tau==0) gamma = 5 ! B2g
!       CHANGING TO TRY TO GET TO WORK: BARRY
        !
      case('D2H')
        !
        gamma = 0
        ideg = 1
        !
        if (mod(K+2,2)==0.and.tau==0) gamma = 1 ! Ag
        if (mod(K+2,2)==0.and.tau==1) gamma = 7 ! B3g
        if (mod(K+2,2)/=0.and.tau==1) gamma = 3 ! B1g
        if (mod(K+2,2)/=0.and.tau==0) gamma = 5 ! B2g
        !
      end select
      !
      !
    case('R_ALPHA_4TAU')
      !
      select case(trim(molec%symmetry))
        !
     case default
        !
        write(out, '(/a,1x,a,1x,a)') &
        'ML_rotsymmetry_C2H6 error: symmetry =', trim(molec%symmetry), 'is unknown'
        stop 'ML_rotsymmetry_C2H6 error: bad symmetry type'
        !
      case('D2H','D2H(M)')
        !
        gamma = 0
        ideg = 1
        !
        if (mod(K+2,2)==0.and.tau==0) gamma = 1 ! Ag
        if (mod(K+2,2)==0.and.tau==1) gamma = 3 ! B1g
        if (mod(K+2,2)/=0.and.tau==1) gamma = 5 ! B2g
        if (mod(K+2,2)/=0.and.tau==0) gamma = 7 ! B3g
        !
      case('D2H(S)')
        !
        gamma = 0
        ideg = 1
        !
        !if (mod(K+2,2)==0.and.tau==0) gamma = 1 ! Ag
        !if (mod(K+2,2)==0.and.tau==1) gamma = 3 ! B1g
        !if (mod(K+2,2)/=0.and.tau==1) gamma = 6 ! B2u
        !if (mod(K+2,2)/=0.and.tau==0) gamma = 8 ! B3u
        !
        if (mod(K+2,2)==0.and.tau==0) gamma = 1 ! Ag
        if (mod(K+2,2)==0.and.tau==1) gamma = 3 ! B1g
        if (mod(K+2,2)/=0.and.tau==1) gamma = 7 ! B3g
        if (mod(K+2,2)/=0.and.tau==0) gamma = 5 ! B2g
        !
      end select
      !
    case('R-R16-BETA16-THETA-TAU')
      !
      select case(trim(molec%symmetry))
      !
      case('C','C(M)')
        !
        gamma = 1
        ideg = 1
        !
     case default
        !
        write(out, '(/a,1x,a,1x,a)') &
        'ML_rotsymmetry_C2H6 error: symmetry =', trim(molec%symmetry), 'is unknown'
        stop 'ML_rotsymmetry_C2H6 error: bad symmetry type'
        !
      end select
      !
    end select
    !
    if (verbose>=5) write(out, '(/a)') 'ML_rotsymmetry_C2H6/end'
    !
  end subroutine ML_rotsymmetry_C2H6


end module mol_c2h6
