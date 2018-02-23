! This unit defines all specific routines for a eight-atomic molecule of ethane-type

module mol_c3h6
  use accuracy
  use moltype

  implicit none

  public ML_coordinate_transform_C3H6, ML_b0_C3H6, ML_symmetry_transformation_C2H6, ML_rotsymmetry_C2H6
  private

  integer(ik), parameter :: verbose = 3 ! Verbosity level

  !--------------------------
  !
  !ZMAT
  !1C;
  !2C, 1, r2,
  !3C, 2, r3, 1, a3,
  !4H, 1, r4, 2, a4, 3, d4;
  !5H, 1, r5, 2, a5, 4, d5;
  !6H, 1, r6, 2, a6, 4, d6;
  !7H, 2, r7, 1, a7, 3, d7;
  !8H, 3, r8, 2, a8, 1, d8;
  !9H, 3, r9, 2, a9, 8, d9;
  !
  !Variables:
  !r2= 1.4456
  !r3= 1.3380
  !a3= 124.64
  !r4= 1.0927
  !a4= 113.32
  !d4=   0.03
  !r5= 1.0913
  !a5= 109.48
  !d5= 120.50
  !r6= 1.0913
  !a6= 109.48
  !d6= 239.50
  !r7= 0.9952
  !a7= 117.31
  !d7= 179.97
  !r8= 0.9921
  !a8= 120.66
  !d8= 179.97
  !r9= 0.9923
  !a9= 121.51
  !d9= 180.0
  !end




  contains


  function ML_coordinate_transform_C3H6(src,ndst,direct) result (dst)
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
                 S1,S2,S3,S4,S5,S6,taubar,tau,S14,S15,S16,S17,S18

    if (verbose>=5) write(out, '(/a)') 'ML_coordinate_transform_C3H6/start'
    !
    nsrc = size(src)
    !
    select case(trim(molec%coords_transform))
      !
    case default
      !
      write(out, '(/a,1x,a,1x,a)') &
      'ML_coordinate_transform_C3H6 error: coordinate type =', trim(molec%coords_transform), 'is unknown'
      stop 'ML_coordinate_transform_C3H6 error: bad coordinate type'
      !
    case('ZMAT_R-A-D_TAU')
      ! ORDER CHANGED HERE AS WANT TBAR TO BE 18TH COORDINATE BUT NEEDS
      !  TO BE 16TH COORDINATE IN Z-MAT
      !
      if (direct) then ! transform from Z-matrix coords to TROVE coords

        dst(1:15) = src(1:13)-molec%local_eq(1:15)
        dst(16:20) = src(17:21)-molec%local_eq(17:21)
        dst(21) = src(16)-molec%local_eq(16)
        !
      else !  transform from TROVE coords to Z-matrix coords
        !
        dst(1:15) = src(1:15)+molec%local_eq(1:15)
        dst(16) = src(21)+molec%local_eq(16)
        dst(17:21) = src(16:20)+molec%local_eq(17:21)
        !
      endif
      !
    end select
    !
    if (verbose>=5) write(out, '(/a)') 'ML_coordinate_transform_C3H6/end'
    !
  end function ML_coordinate_transform_C3H6



  subroutine ML_b0_C3H6(Npoints,Natoms,b0,rho_i,rho_ref,rho_borders)
    !
    ! Reference structure
    !
    integer(ik),intent(in) :: Npoints, Natoms
    real(ark),intent(out) :: b0(Natoms,3,0:Npoints)
    real(ark),intent(inout),optional :: rho_i(0:Npoints)
    real(ark),intent(out),optional :: rho_ref
    real(ark),intent(in),optional :: rho_borders(2)  ! rhomim, rhomax - borders
    !
    real(ark) :: a0(molec%Natoms,3),CM_shift,tau,alpha0,alpha,theta,r,r12,tau14,tau15,tau16
    real(ark) :: transform(3,3),phi
    integer(ik) :: i, n, iatom, ix
    !
    if (verbose>=5) write(out, '(/a)') 'ML_b0_C3H6/start'
    !
    rC1e      = molec%req(1)
    rC2e      = molec%req(2)
    rH1e      = molec%req(3)
    rH2e      = molec%req(4)
    rH3e      = molec%req(5)
    rH4e      = molec%req(6)
    rH5e      = molec%req(7)
    rH6e      = molec%req(8)
    !
    alpha1e    = molec%alphaeq(1)*deg
    alpha2e    = molec%alphaeq(2)*deg
    alpha3e    = molec%alphaeq(3)*deg
    alpha4e    = molec%alphaeq(4)*deg
    alpha5e    = molec%alphaeq(5)*deg
    alpha6e    = molec%alphaeq(6)*deg
    !
    delta1e    = molec%taueq(1)*deg
    delta2e    = molec%taueq(2)*deg
    delta3e    = molec%taueq(3)*deg
    delta4e    = molec%taueq(4)*deg
    delta5e    = molec%taueq(5)*deg
    taue       = molec%taueq(6)*deg

    !
    !if (any(molec%req(2:7)/=r)) then
    !  write(out,"('ML_b0_C3H6 error: eq-m r2-r7 are not all the same:',6f12.5)") molec%req(2:7)
    !  stop 'ML_b0_C3H6 error: eq-m r2-r7 are not all the same'
    !endif
    !
    !if (any(molec%alphaeq(1:6)/=alpha)) then
    !  write(out,"('ML_b0_C3H6 error: eq-m alphas are not all the same:',6f12.5)") molec%alphaeq(1:6)
    !  stop 'ML_b0_C3H6 error: eq-m alphas are not all the same'
    !endif
    !
    !
    a0 = 0
    !
    a0(1,:) = 0
    a0(2,:) = 0
    !
    a0(3,1) = r*sin(alpha)
    a0(3,2) = 0
    a0(3,3) = r*cos(alpha)
    !
    a0(4,1) = r*sin(alpha)*cos(tau14)
    a0(4,2) = r*sin(alpha)*sin(tau14)
    a0(4,3) =-r*cos(alpha)
    !
    a0(5,1) = r*sin(alpha)*cos(theta)
    a0(5,2) =-r*sin(alpha)*sin(theta)
    a0(5,3) = r*cos(alpha)
    !
    a0(6,1) = r*sin(alpha)*cos(tau15)
    a0(6,2) = r*sin(alpha)*sin(tau15)
    a0(6,3) =-r*cos(alpha)
    !
    a0(7,1) = r*sin(alpha)*cos(theta)
    a0(7,2) = r*sin(alpha)*sin(theta)
    a0(7,3) = r*cos(alpha)
    !
    a0(8,1) = r*sin(alpha)*cos(tau16)
    a0(8,2) = r*sin(alpha)*sin(tau16)
    a0(8,3) =-r*cos(alpha)
    !
    a0(1:7:2,3) =a0(1:7:2,3)-r12*0.5_ark
    a0(2:8:2,3) =a0(2:8:2,3)+r12*0.5_ark
    !
    !
    do ix=1, 3
      CM_shift = sum(a0(:,ix)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
      a0(:,ix) = a0(:,ix) - CM_shift
    enddo
    !
    !
    phi = tau*0.5_ark
    !
    transform = 0 
    transform(3,3) = 1.0_ark
    transform(1,1) = cos(phi)
    transform(1,2) = sin(phi)
    transform(2,1) = -sin(phi)
    transform(2,2) = cos(phi)
    !
    do n = 1,Natoms
      a0(n,:) = matmul(transform,a0(n,:))
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
       ! 
       if (.not.present(rho_borders).or..not.present(rho_ref)) then  
          write(out,"('ML_b0_ABCD: rho_borders and rho_ref must be presented if Npoints ne 0 ')") 
          stop 'ML_b0_ABCD: rho_borders or rho_ref not specified '
       endif
       !
       rho_ref = molec%taueq(1)
       !
       do i = 0,npoints
          !
          tau = rho_i(i)
          !
          tau14 = tau
          tau15 = tau+theta
          tau16 = tau-theta

          b0(1,:,i) = 0
          b0(2,:,i) = 0
          !
          b0(3,1,i) = r*sin(alpha)
          b0(3,2,i) = 0
          b0(3,3,i) = r*cos(alpha)
          !
          b0(4,1,i) = r*sin(alpha)*cos(tau14)
          b0(4,2,i) = r*sin(alpha)*sin(tau14)
          b0(4,3,i) =-r*cos(alpha)
          !
          b0(5,1,i) = r*sin(alpha)*cos(theta)
          b0(5,2,i) =-r*sin(alpha)*sin(theta)
          b0(5,3,i) = r*cos(alpha)
          !
          b0(6,1,i) = r*sin(alpha)*cos(tau15)
          b0(6,2,i) = r*sin(alpha)*sin(tau15)
          b0(6,3,i) =-r*cos(alpha)
          !
          b0(7,1,i) = r*sin(alpha)*cos(theta)
          b0(7,2,i) = r*sin(alpha)*sin(theta)
          b0(7,3,i) = r*cos(alpha)
          !
          b0(8,1,i) = r*sin(alpha)*cos(tau16)
          b0(8,2,i) = r*sin(alpha)*sin(tau16)
          b0(8,3,i) =-r*cos(alpha)
          !
          b0(1:7:2,3,i) =b0(1:7:2,3,i)-r12*0.5_ark
          b0(2:8:2,3,i) =b0(2:8:2,3,i)+r12*0.5_ark
          !
          ! Find center of mass
          !
          do n = 1,3 
            CM_shift = sum(b0(:,n,i)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
            b0(:,n,i) = b0(:,n,i) - CM_shift
          enddo
          !
          !
          phi = tau*0.5_ark
          !
          transform = 0 
          transform(3,3) = 1.0_ark
          transform(1,1) = cos(phi)
          transform(1,2) = sin(phi)
          transform(2,1) = -sin(phi)
          transform(2,2) = cos(phi)
          !
          do n = 1,Natoms
            b0(n,:,i) = matmul(transform,b0(n,:,i))
          enddo
          !
          !call MLorienting_a0(molec%Natoms,molec%AtomMasses,b0(:,:,0),transform)
          !
          if (verbose>=3) then
            write(out, '(i5)') 8
            write(out,'(a)') ""
            write(out, '(1x,a,1x,3(1x,es16.8))') 'C', b0(1,1:3,i)
            write(out, '(1x,a,1x,3(1x,es16.8))') 'C', b0(2,1:3,i)
            do iatom=3, Natoms
              write(out, '(1x,a,1x,3(1x,es16.8))') 'H', b0(iatom,1:3,i)
            enddo
          endif
          !
       enddo
       !
    endif
    !
    if (verbose>=5) write(out, '(/a)') 'ML_b0_C3H6/end'
    !
  end subroutine ML_b0_C3H6



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
        case (2) ! !C3+/(132)(465)
          !
          dst(1) = src(1)
          dst(2) = src(4)
          dst(3) = src(2)
          dst(4) = src(3)
          dst(5) = src(7)
          dst(6) = src(5)
          dst(7) = src(6)
          dst(8) = src(10)
          dst(9) = src(8)
          dst(10) = src(9)
          dst(11) = src(13)
          dst(12) = src(11)
          dst(13) = src(12)
          dst(14) = -a*src(14) - b*src(15)
          dst(15) =  b*src(14) - a*src(15) 
          dst(16) = -a*src(16) - b*src(17)
          dst(17) =  b*src(16) - a*src(17)
          dst(18) = src(18)
          !
        case (3) !C3-/(123)(465)
          !
          dst(1) = src(1)
          dst(2) = src(3)
          dst(3) = src(4)
          dst(4) = src(2)
          dst(5) = src(6)
          dst(6) = src(7)
          dst(7) = src(5)
          dst(8) = src(9)
          dst(9) = src(10)
          dst(10) = src(8)
          dst(11) = src(12)
          dst(12) = src(13)
          dst(13) = src(11)
          dst(14) = -a*src(14) + b*src(15)
          dst(15) = -b*src(14) - a*src(15) 
          dst(16) = -a*src(16) + b*src(17)
          dst(17) = -b*src(16) - a*src(17)
          dst(18) =  src(18)
          !
        case (4) ! C2/(16)(24)(35)(78)
          !
          dst(1) = src(1)
          dst(2) = src(6)
          dst(3) = src(5)
          dst(4) = src(7)
          dst(5) = src(3)
          dst(6) = src(2)
          dst(7) = src(4)
          dst(8) = src(12)
          dst(9) = src(11)
          dst(10) = src(13)
          dst(11) = src(9)
          dst(12) = src(8)
          dst(13) = src(10)
          dst(14) = -a*src(16) + b*src(17)
          dst(15) = -b*src(16) - a*src(17) 
          dst(16) = -a*src(14) + b*src(15)
          dst(17) =  b*src(14) + a*src(15)
          dst(18) = src(18)
          !
        case (5) !C2'/(16)(24)(35)(78)
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
          dst(14) = src(17)
          dst(15) = src(18) 
          dst(16) = src(14)
          dst(17) = src(15)
          dst(18) = src(18)
          !
        case (6) ! C2''(16)(24)(35)(78)
          !
          dst(1) = src(1)
          dst(2) = src(7)
          dst(3) = src(6)
          dst(4) = src(5)
          dst(5) = src(4)
          dst(6) = src(3)
          dst(7) = src(2)
          dst(8) = src(13)
          dst(9) = src(12)
          dst(10) = src(11)
          dst(11) = src(10)
          dst(12) = src(9)
          dst(13) = src(8)
          dst(14) = -a*src(16) - b*src(17)
          dst(15) =  b*src(16) - a*src(17) 
          dst(16) = -a*src(14) - b*src(16)
          dst(17) = -b*src(14) + a*src(15)
          dst(18) = src(18)
          !
        case (7) ! i/(14)(26)(35)(78)*
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
          dst(14) = src(16)
          dst(15) = -src(17) 
          dst(16) = src(14)
          dst(17) = -src(15)
          dst(18) = -src(18)
          !
        case (8) ! S6/(163425)(78)*
          !
          dst(1) = src(1)
          dst(2) = src(6)
          dst(3) = src(7)
          dst(4) = src(5)
          dst(5) = src(3)
          dst(6) = src(4)
          dst(7) = src(2)
          dst(8) = src(12)
          dst(9) = src(13)
          dst(10) = src(11)
          dst(11) = src(9)
          dst(12) = src(10)
          dst(13) = src(8)
          dst(14) = -a*src(16) + b*src(17)
          dst(15) =  b*src(16) + a*src(17) 
          dst(16) = -a*src(14) + b*src(16)
          dst(17) =  -b*src(14) - a*src(15)
          dst(18) = -src(18)
          !
        case (9) !S6'/(14)(26)(35)(78)*
          !
          dst(1) = src(1)
          dst(2) = src(7)
          dst(3) = src(5)
          dst(4) = src(6)
          dst(5) = src(4)
          dst(6) = src(2)
          dst(7) = src(3)
          dst(8) = src(13)
          dst(9) = src(11)
          dst(10) = src(12)
          dst(11) = src(10)
          dst(12) = src(8)
          dst(13) = src(9)
          dst(14) = -a*src(16) - b*src(17)
          dst(15) = -b*src(16) + a*src(17) 
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
          dst(5) = src(7)
          dst(6) = src(6)
          dst(7) = src(5)
          dst(8) = src(10)
          dst(9) = src(9)
          dst(10) = src(8)
          dst(11) = src(13)
          dst(12) = src(12)
          dst(13) = src(11)
          dst(14) = -a*src(14) - b*src(15) 
          dst(15) = -b*src(14) - a*src(15)
          dst(16) = -a*src(16) - b*src(17)
          dst(17) = -b*src(16) + a*src(17)
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
          dst(14) = src(14) 
          dst(15) = -src(15) 
          dst(16) = src(16)
          dst(17) = -src(17)
          dst(18) = -src(18)
          !
        case (12) !sigmad''/(12)(46)*
          !
          dst(1) = src(1)
          dst(2) = src(3)
          dst(3) = src(2)
          dst(4) = src(4)
          dst(5) = src(6)
          dst(6) = src(5)
          dst(7) = src(7)
          dst(8) = src(9)
          dst(9) = src(8)
          dst(10) = src(10)
          dst(11) = src(12)
          dst(12) = src(11)
          dst(13) = src(13)
          dst(14) = -a*src(14) + b*src(15)
          dst(15) = -b*src(14) - a*src(15) 
          dst(16) = -a*src(16) + b*src(17)
          dst(17) =  b*src(16) + a*src(17)
          dst(18) = -src(18)
        end select
        !
      end select
      !
    case('R-R16-BETA16-THETA-TAU')
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
        case (2) ! !C3+/(132)(465)
          !
          dst(1) = src(1)
          dst(2) = src(4)
          dst(3) = src(2)
          dst(4) = src(3)
          dst(5) = src(7)
          dst(6) = src(5)
          dst(7) = src(6)
          dst(8) = src(10)
          dst(9) = src(8)
          dst(10) = src(9)
          dst(11) = src(13)
          dst(12) = src(11)
          dst(13) = src(12)
          dst(14) = -a*src(14) - b*src(15)
          dst(15) =  b*src(14) - a*src(15) 
          dst(16) = -a*src(16) - b*src(17)
          dst(17) =  b*src(16) - a*src(17)
          dst(18) = src(18)
          !
        case (3) !C3-/(123)(465)
          !
          dst(1) = src(1)
          dst(2) = src(3)
          dst(3) = src(4)
          dst(4) = src(2)
          dst(5) = src(6)
          dst(6) = src(7)
          dst(7) = src(5)
          dst(8) = src(9)
          dst(9) = src(10)
          dst(10) = src(8)
          dst(11) = src(12)
          dst(12) = src(13)
          dst(13) = src(11)
          dst(14) = -a*src(14) + b*src(15)
          dst(15) = -b*src(14) - a*src(15) 
          dst(16) = -a*src(16) + b*src(17)
          dst(17) = -b*src(16) - a*src(17)
          dst(18) =  src(18)
          !
        case (4) ! C2/(16)(24)(35)(78)
          !
          dst(1) = src(1)
          dst(2) = src(6)
          dst(3) = src(5)
          dst(4) = src(7)
          dst(5) = src(3)
          dst(6) = src(2)
          dst(7) = src(4)
          dst(8) = src(12)
          dst(9) = src(11)
          dst(10) = src(13)
          dst(11) = src(9)
          dst(12) = src(8)
          dst(13) = src(10)
          dst(14) = -a*src(16) + b*src(17)
          dst(15) =  b*src(16) + a*src(17) 
          dst(16) = -a*src(14) + b*src(15)
          dst(17) =  b*src(14) + a*src(15)
          dst(18) = src(18)
          !
        case (5) !C2'/(16)(24)(35)(78)
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
          dst(15) =-src(17) 
          dst(16) = src(14)
          dst(17) =-src(15)
          dst(18) = src(18)
          !
        case (6) ! C2''(16)(24)(35)(78)
          !
          dst(1) = src(1)
          dst(2) = src(7)
          dst(3) = src(6)
          dst(4) = src(5)
          dst(5) = src(4)
          dst(6) = src(3)
          dst(7) = src(2)
          dst(8) = src(13)
          dst(9) = src(12)
          dst(10) = src(11)
          dst(11) = src(10)
          dst(12) = src(9)
          dst(13) = src(8)
          dst(14) = -a*src(16) - b*src(17)
          dst(15) = -b*src(16) + a*src(17) 
          dst(16) = -a*src(14) - b*src(15)
          dst(17) = -b*src(14) + a*src(15)
          dst(18) = src(18)
          !
        case (7) ! i/(14)(26)(35)(78)*
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
          dst(3) = src(7)
          dst(4) = src(5)
          dst(5) = src(3)
          dst(6) = src(4)
          dst(7) = src(2)
          dst(8) = src(12)
          dst(9) = src(13)
          dst(10) = src(11)
          dst(11) = src(9)
          dst(12) = src(10)
          dst(13) = src(8)
          dst(14) = -a*src(16) + b*src(17)
          dst(15) = -b*src(16) - a*src(17) 
          dst(16) = -a*src(14) + b*src(15)
          dst(17) = -b*src(14) - a*src(15)
          dst(18) = -src(18)
          !
        case (9) !S6'/(14)(26)(35)(78)*
          !
          dst(1) = src(1)
          dst(2) = src(7)
          dst(3) = src(5)
          dst(4) = src(6)
          dst(5) = src(4)
          dst(6) = src(2)
          dst(7) = src(3)
          dst(8) = src(13)
          dst(9) = src(11)
          dst(10) = src(12)
          dst(11) = src(10)
          dst(12) = src(8)
          dst(13) = src(9)
          dst(14) = -a*src(16) - b*src(17)
          dst(15) =  b*src(16) - a*src(17) 
          dst(16) = -a*src(14) - b*src(15)
          dst(17) =  b*src(14) - a*src(15)
          dst(18) = -src(18)
          !
        case (10) !sigmad/(12)(46)*
          !
          dst(1) = src(1)
          dst(2) = src(4)
          dst(3) = src(3)
          dst(4) = src(2)
          dst(5) = src(7)
          dst(6) = src(6)
          dst(7) = src(5)
          dst(8) = src(10)
          dst(9) = src(9)
          dst(10) = src(8)
          dst(11) = src(13)
          dst(12) = src(12)
          dst(13) = src(11)
          dst(14) = -a*src(14) - b*src(15) 
          dst(15) = -b*src(14) + a*src(15)
          dst(16) = -a*src(16) - b*src(17)
          dst(17) = -b*src(16) + a*src(17)
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
          dst(14) = src(14) 
          dst(15) =-src(15) 
          dst(16) = src(16)
          dst(17) =-src(17)
          dst(18) =-src(18)
          !
        case (12) !sigmad''/(12)(46)*
          !
          dst(1) = src(1)
          dst(2) = src(3)
          dst(3) = src(2)
          dst(4) = src(4)
          dst(5) = src(6)
          dst(6) = src(5)
          dst(7) = src(7)
          dst(8) = src(9)
          dst(9) = src(8)
          dst(10) = src(10)
          dst(11) = src(12)
          dst(12) = src(11)
          dst(13) = src(13)
          dst(14) = -a*src(14) + b*src(15)
          dst(15) =  b*src(14) + a*src(15) 
          dst(16) = -a*src(16) + b*src(17)
          dst(17) =  b*src(16) + a*src(17)
          dst(18) = -src(18)
          !
        end select
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
      case('D3D(M)')
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
      case('D3D(M)')
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


end module mol_c3h6
