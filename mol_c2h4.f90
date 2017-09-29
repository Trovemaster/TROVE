! This unit defines all specific routines for a six-atomic molecule of ethylene-type

module mol_c2h4
  use accuracy
  use moltype

  implicit none

  public ML_coordinate_transform_C2H4, ML_b0_C2H4, ML_symmetry_transformation_C2H4, ML_rotsymmetry_C2H4
  private

  integer(ik), parameter :: verbose = 3 ! Verbosity level

  !--------------------------
  !      ZMAT_2BETA_1TAU
  !--------------------------
  !
  ! ZMAT
  !  1 C   0  0  0  0  12.00000000
  !  2 C   1  0  0  0  12.00000000
  !  3 H   1  2  0  0   1.00782505
  !  4 H   1  2  3 -2   1.00782505
  !  5 H   2  1  3  2   1.00782505
  !  6 H   2  1  5 -2   1.00782505
  ! end
  !
  ! ( 3       5 )
  ! (  \     /  )
  ! (   1 - 2   )
  ! (  /     \  )
  ! ( 4       6 )
  !
  ! order of coordinates:
  ! 1    r_12
  ! 2    r_31
  ! 3    r_41
  ! 4    r_52
  ! 5    r_62
  ! 6    alpha_312
  ! 7    alpha_412
  ! 8    alpha_521
  ! 9    alpha_621
  ! 10   beta_4123
  ! 11   tau_5213
  ! 12   beta_6215
  !--------------------------

  !--------------------------
  !      R_ALPHA_4TAU
  !--------------------------
  !
  ! ZMAT
  !  C   0  0  0  0  12.00000000
  !  C   1  0  0  0  12.00000000
  !  H   1  2  0  0   1.00782505
  !  H   2  1  3  2   1.00782505
  !  H   1  2  4  2   1.00782505
  !  H   2  1  5  2   1.00782505
  ! end
  !
  ! ( 3       4 )
  ! (  \     /  )
  ! (   1 - 2   )
  ! (  /     \  )
  ! ( 5       6 )
  !
  ! order of coordinates:
  ! zmat1  = r_12
  ! zmat2  = r_31
  ! zmat3  = r_42
  ! zmat4  = r_51
  ! zmat5  = r_62
  ! zmat6  = alpha_312
  ! zmat7  = alpha_421
  ! zmat8  = alpha_512
  ! zmat9  = alpha_621
  ! zmat10 = -tau_4213
  ! zmat11 = tau_5124
  ! zmat12 = -tau_6215
  ! redundant tau_6213 = 2 pi-tau5124-tau4213-tau5126
  !--------------------------



  contains


  function ML_coordinate_transform_C2H4(src,ndst,direct) result (dst)
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
    if (verbose>=5) write(out, '(/a)') 'ML_coordinate_transform_C2H4/start'
    !
    nsrc = size(src)
    !
    select case(trim(molec%coords_transform))
      !
    case default
      !
      write(out, '(/a,1x,a,1x,a)') &
      'ML_coordinate_transform_C2H4 error: coordinate type =', trim(molec%coords_transform), 'is unknown'
      stop 'ML_coordinate_transform_C2H4 error: bad coordinate type'
      !
    case('ZMAT_2BETA_1TAU','C2H4_2BETA_1TAU')
      !
      if (direct) then ! transform from Z-matrix coords to TROVE coords
        !
        b1_zmat   = src(10)
        tau1_zmat = src(11)
        b2_zmat   = src(12)
        !
        db1 = b1_zmat - pi
        db2 = b2_zmat - pi
        !
        if (tau1_zmat>pi) then
          tau1 = 2.0_ark*pi - tau1_zmat
        elseif (tau1_zmat<pi) then
          tau1 = -tau1_zmat
        endif
        !
        tau2 = db1 + db2 + tau1
        dtau = tau1 + tau2
        !
        dst(1:9) = src(1:9)-molec%local_eq(1:9)
        dst(10:12) = (/db1, db2, dtau/)
        !
      else !  transform from TROVE coords to Z-matrix coords
        !
        db1  = src(10)
        db2  = src(11)
        dtau = src(12)
        !
        b1_zmat = db1 + pi
        b2_zmat = db2 + pi
        !
        tau1 = (dtau - db1 - db2)*0.5_ark
        if (sign(1.0_ark,tau1)>0.0) then
          tau1_zmat = 2.0_ark*pi - tau1
        elseif (sign(1.0_ark,tau1)<0.0) then
          tau1_zmat = -tau1
        else
          tau1_zmat = tau1
        endif
        !
        dst(1:9) = src(1:9)+molec%local_eq(1:9)
        dst(10:12) = (/b1_zmat, tau1_zmat, b2_zmat/)
        !
      endif
      !
    case('R_ALPHA_4TAU')
      !
      if (direct) then
        !
        dst(1:9) = src(1:9) - molec%local_eq(1:9)
        !
        tau4213 =-src(10)
        if (tau4213>pi) tau4213 = tau4213 - 2.0_ark*pi

        tau5124 = src(11)
        tau5126 =-src(12)
        if (tau5126>pi) tau5126 = tau5126 - 2.0_ark*pi
        !
        tau6213 = 2.0_ark*pi-tau5124-tau4213-tau5126
        tau6213 = mod(tau6213+2.0_ark*pi,2.0_ark*pi)
        !
        dst(10) = tau4213-molec%local_eq(10)
        dst(11) = tau5126-molec%local_eq(12)
        dst(12) = tau6213-tau5124
        !
      else ! not direct
        !
        dst(1:9) = src(1:9) + molec%local_eq(1:9)
        !
        tau4213 = src(10) + molec%local_eq(10)
        tau4213 = mod(tau4213+2.0_ark*pi,2.0_ark*pi)
        !
        tau5126 = src(11) + molec%local_eq(12)
        tau5126 = mod(tau5126+2.0_ark*pi,2.0_ark*pi)
        !
        tau5124 = 0.5_ark*(2.0_ark*pi-tau5126-tau4213-src(12))
        tau5124 = mod(tau5124+2.0_ark*pi,2.0_ark*pi)
        !
        dst(10) =-tau4213
        dst(11) = tau5124
        dst(12) =-tau5126
        !
      endif
      !
    end select
    !
    if (verbose>=5) write(out, '(/a)') 'ML_coordinate_transform_C2H4/end'
    !
  end function ML_coordinate_transform_C2H4



  subroutine ML_b0_C2H4(Npoints,Natoms,b0,rho_i,rho_ref,rho_borders)
    !
    ! Reference structure
    !
    integer(ik),intent(in) :: Npoints, Natoms
    real(ark),intent(out) :: b0(Natoms,3,0:Npoints)
    real(ark),intent(in),optional :: rho_i(0:Npoints)
    real(ark),intent(out),optional :: rho_ref
    real(ark),intent(in),optional :: rho_borders(2)  ! rhomim, rhomax - borders
    !
    real(ark) :: a0(molec%Natoms,3),CM_shift,tau,alpha0
    integer(ik) :: i, n, iatom
    !
    if (verbose>=5) write(out, '(/a)') 'ML_b0_C2H4/start'
    !
    a0 = 0.0_ark
    !
    if (all(molec%dihedtype(1:3)==(/-2_ik,2_ik,-2_ik/))) then
      !
      a0(1,1) = 0.0_ark
      a0(1,2) = 0.0_ark
      a0(1,3) = -molec%req(1)*0.5_ark
      !
      a0(2,1) = 0.0_ark
      a0(2,2) = 0.0_ark
      a0(2,3) = molec%req(1)*0.5_ark
      !
      a0(3,1) = 0.0_ark
      a0(3,2) = molec%req(2)*sin(molec%alphaeq(1))
      a0(3,3) = molec%req(2)*cos(molec%alphaeq(1))-molec%req(1)*0.5_ark
      !
      a0(4,1) = a0(3,2)*sin(molec%taueq(1))+a0(3,1)*cos(molec%taueq(1))
      a0(4,2) = a0(3,2)*cos(molec%taueq(1))-a0(3,1)*sin(molec%taueq(1))
      a0(4,3) = molec%req(2)*cos(molec%alphaeq(1))-molec%req(1)*0.5_ark
      !
      a0(5,1) = -molec%req(2)*sin(molec%taueq(2))*sin(molec%alphaeq(1))
      a0(5,2) =  molec%req(2)*cos(molec%taueq(2))*sin(molec%alphaeq(1))
      a0(5,3) = -molec%req(2)*cos(molec%alphaeq(1))+molec%req(1)*0.5_ark
      !
      a0(6,1) =  a0(5,2)*sin(molec%taueq(3))+a0(5,1)*cos(molec%taueq(3))
      a0(6,2) =  a0(5,2)*cos(molec%taueq(3))-a0(5,1)*sin(molec%taueq(3))
      a0(6,3) = -molec%req(2)*cos(molec%alphaeq(1))+molec%req(1)*0.5_ark
      !
    elseif (all(molec%dihedtype(1:3)==(/2_ik,2_ik,2_ik/))) then
      !
      alpha0 = molec%alphaeq(1)
      !
      a0(1,1) = 0.0_ark
      a0(1,2) = 0.0_ark
      a0(1,3) = molec%req(1)*0.5_ark
      !
      a0(2,2) = 0.0_ark
      a0(2,1) = 0.0_ark
      a0(2,3) = -molec%req(1)*0.5_ark
      !
      a0(3,1) = molec%req(2)*sin(alpha0)
      a0(3,2) = 0
      a0(3,3) =-molec%req(2)*cos(alpha0)+molec%req(1)*0.5_ark
      !
      a0(4,1) = molec%req(2)*sin(alpha0)
      a0(4,2) = 0
      a0(4,3) = molec%req(2)*cos(alpha0)-molec%req(1)*0.5_ark
      !
      a0(5,1) =-molec%req(2)*sin(alpha0)
      a0(5,2) = 0
      a0(5,3) =-molec%req(2)*cos(alpha0)+molec%req(1)*0.5_ark
      !
      a0(6,1) =-molec%req(2)*sin(alpha0)
      a0(6,2) = 0
      a0(6,3) = molec%req(2)*cos(alpha0)-molec%req(1)*0.5_ark
      !
      !a0(5:6,:) = a0(6:5:-1,:)
      !
      a0  = cshift(a0,shift=-1,dim=2)
      !
      !a0(:,1:2) = a0(:,2:1:-1)
      !
      a0(:,3) =-a0(:,3)
      !
    else
      !
      write(out, '(/a,3(1x,i3),a)') &
      'ML_b0_C2H4 error: combination of dihedral angles of types = (', molec%dihedtype(1:3), ') is not permitted'
      stop 'ML_b0_C2H4 error: bad dihedral angles'
      !
    endif
    !
    do n=1, 3
      CM_shift = sum(a0(:,n)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
      a0(:,n) = a0(:,n) - CM_shift
    enddo
    !
    if (verbose>=3) then
      do iatom=1, Natoms
        write(out, '(1x,a,1x,3(1x,es16.8))') 'H', a0(iatom,1:3)
      enddo
    endif
    !
    b0(:,:,0) = a0(:,:)
    !
    if (Npoints/=0) then
      !
      if (.not.present(rho_borders).or..not.present(rho_ref)) then
        write(out, '(/a)') 'ML_b0_C2H4: rho_borders and rho_ref must be presented if Npoints > 0'
        stop 'ML_b0_C2H4: rho_borders or rho_ref is not specified'
      endif
      !
      rho_ref = 0
      !
      do i=0, npoints
        !
        tau = rho_i(i)
        !
        if (all(molec%dihedtype(1:3)==(/-2_ik,2_ik,-2_ik/))) then
          !
          b0(1,1,i) = 0.0_ark
          b0(1,2,i) = 0.0_ark
          b0(1,3,i) = -molec%req(1)*0.5_ark
          !
          b0(2,1,i) = 0.0_ark
          b0(2,2,i) = 0.0_ark
          b0(2,3,i) = molec%req(1)*0.5_ark
          !
          b0(3,1,i) = 0.0_ark
          b0(3,2,i) = molec%req(2)*sin(molec%alphaeq(1))
          b0(3,3,i) = molec%req(2)*cos(molec%alphaeq(1))-molec%req(1)*0.5_ark
          !
          b0(4,1,i) = a0(3,2)*sin(molec%taueq(1))+a0(3,1)*cos(molec%taueq(1))
          b0(4,2,i) = a0(3,2)*cos(molec%taueq(1))-a0(3,1)*sin(molec%taueq(1))
          b0(4,3,i) = molec%req(2)*cos(molec%alphaeq(1))-molec%req(1)*0.5_ark
          !
          b0(5,1,i) = -molec%req(2)*sin(tau)*sin(molec%alphaeq(1))
          b0(5,2,i) =  molec%req(2)*cos(tau)*sin(molec%alphaeq(1))
          b0(5,3,i) = -molec%req(2)*cos(molec%alphaeq(1))+molec%req(1)*0.5_ark
          !
          b0(6,1,i) =  a0(5,2)*sin(molec%taueq(3))+a0(5,1)*cos(molec%taueq(3))
          b0(6,2,i) =  a0(5,2)*cos(molec%taueq(3))-a0(5,1)*sin(molec%taueq(3))
          b0(6,3,i) = -molec%req(2)*cos(molec%alphaeq(1))+molec%req(1)*0.5_ark
          !
        else
          !
          write(out, '(/a,3(1x,i3),a)') &
          'ML_b0_C2H4 error: combination of dihedral angles of types = (', molec%dihedtype(1:3), ') is not permitted'
          stop 'ML_b0_C2H4 error: bad dihedral angles'
          !
        endif
        !
        call MLorienting_a0(molec%Natoms,molec%AtomMasses,b0(:,:,i))
        !
        do n = 1,3
          CM_shift = sum(b0(:,n,i)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
          b0(:,n,i) = b0(:,n,i) - CM_shift
        enddo
        !
        if (verbose>=5) then
          write(out, '(1x,a,1x,i6,100(1x,es16.8))') 'b0', i, b0(:,:,i)
        endif
        !
      enddo
      !
    endif
    !
    if (verbose>=5) write(out, '(/a)') 'ML_b0_C2H4/end'
    !
  end subroutine ML_b0_C2H4



  subroutine ML_symmetry_transformation_C2H4(ioper, nmodes, src, dst)
    !
    ! Symmetry transformation rules of coordinates
    !
    integer(ik), intent(in)  ::  ioper
    integer(ik), intent(in)  ::  nmodes
    real(ark), intent(in)    ::  src(1:nmodes)
    real(ark), intent(out)   ::  dst(1:nmodes)
    !
    if (verbose>=5) write(out, '(/a)') 'ML_symmetry_transformation_C2H4/start'
    !
    select case(trim(molec%coords_transform))
      !
    case default
      !
      write(out, '(/a,1x,a,1x,a)') &
      'ML_symmetry_transformation_C2H4 error: coordinate type =', trim(molec%coords_transform), 'is unknown'
      stop 'ML_symmetry_transformation_C2H4 error: bad coordinate type'
      !
    case('ZMAT_2BETA_1TAU')
      !
      select case(trim(molec%symmetry))
        !
      case default
        !
        write(out, '(/a,1x,a,1x,a)') &
        'ML_symmetry_transformation_C2H4 error: symmetry =', trim(molec%symmetry), 'is unknown'
        stop
        !
      case('D2H(M)')
        !
        select case(ioper)
          !
        case default
          !
          write(out, '(/a,1x,i3,1x,a)') &
          'ML_symmetry_transformation_C2H4 error: symmetry operation ', ioper, 'is unknown'
          stop
          !
        case (1) ! E
          !
          dst(1:12) = src(1:12)
          !
        case (2) ! C2a == (34)(56)
          !
          dst(1) = src(1)
          dst(2) = src(3)
          dst(3) = src(2)
          dst(4) = src(5)
          dst(5) = src(4)
          dst(6) = src(7)
          dst(7) = src(6)
          dst(8) = src(9)
          dst(9) = src(8)
          dst(10) = -src(10)
          dst(11) = -src(11)
          dst(12) = src(12)
          !
        case (3) ! C2b == (35)(46)(12)
          !
          dst(1) = src(1)
          dst(2) = src(4)
          dst(3) = src(5)
          dst(4) = src(2)
          dst(5) = src(3)
          dst(6) = src(8)
          dst(7) = src(9)
          dst(8) = src(6)
          dst(9) = src(7)
          dst(10) = src(11)
          dst(11) = src(10)
          dst(12) = src(12)
          !
        case (4) ! C2c == (36)(45)(12)
          !
          dst(1) = src(1)
          dst(2) = src(5)
          dst(3) = src(4)
          dst(4) = src(3)
          dst(5) = src(2)
          dst(6) = src(9)
          dst(7) = src(8)
          dst(8) = src(7)
          dst(9) = src(6)
          dst(10) = -src(11)
          dst(11) = -src(10)
          dst(12) = src(12)
          !
        case (5) ! sigma_ab == E*
          !
          dst(1) = src(1)
          dst(2) = src(2)
          dst(3) = src(3)
          dst(4) = src(4)
          dst(5) = src(5)
          dst(6) = src(6)
          dst(7) = src(7)
          dst(8) = src(8)
          dst(9) = src(9)
          dst(10) = -src(10)
          dst(11) = -src(11)
          dst(12) = -src(12)
          !
        case (6) ! sigma_ac == (34)(56)*
          !
          dst(1) = src(1)
          dst(2) = src(3)
          dst(3) = src(2)
          dst(4) = src(5)
          dst(5) = src(4)
          dst(6) = src(7)
          dst(7) = src(6)
          dst(8) = src(9)
          dst(9) = src(8)
          dst(10) = src(10)
          dst(11) = src(11)
          dst(12) = -src(12)
          !
        case (7) ! sigma_bc == (35)(46)(12)*
          !
          dst(1) = src(1)
          dst(2) = src(4)
          dst(3) = src(5)
          dst(4) = src(2)
          dst(5) = src(3)
          dst(6) = src(8)
          dst(7) = src(9)
          dst(8) = src(6)
          dst(9) = src(7)
          dst(10) = -src(11)
          dst(11) = -src(10)
          dst(12) = -src(12)
          !
        case (8) ! i == (36)(45)(12)*
          !
          dst(1) = src(1)
          dst(2) = src(5)
          dst(3) = src(4)
          dst(4) = src(3)
          dst(5) = src(2)
          dst(6) = src(9)
          dst(7) = src(8)
          dst(8) = src(7)
          dst(9) = src(6)
          dst(10) = src(11)
          dst(11) = src(10)
          dst(12) = -src(12)
          !
        end select
        !
       case('D2H')
        !
        select case(ioper)
          !
        case default
          !
          write(out, '(/a,1x,i3,1x,a)') &
          'ML_symmetry_transformation_C2H4 error: symmetry operation ', ioper, 'is unknown'
          stop
          !
        case (1) ! E
          !
          dst(1:12) = src(1:12)
          !
        case (4) ! C2a == (34)(56)
          !
          dst(1) = src(1)
          dst(2) = src(3)
          dst(3) = src(2)
          dst(4) = src(5)
          dst(5) = src(4)
          dst(6) = src(7)
          dst(7) = src(6)
          dst(8) = src(9)
          dst(9) = src(8)
          dst(10) = -src(10)
          dst(11) = -src(11)
          dst(12) = src(12)
          !
        case (3) ! C2b == (35)(46)(12)
          !
          dst(1) = src(1)
          dst(2) = src(4)
          dst(3) = src(5)
          dst(4) = src(2)
          dst(5) = src(3)
          dst(6) = src(8)
          dst(7) = src(9)
          dst(8) = src(6)
          dst(9) = src(7)
          dst(10) = src(11)
          dst(11) = src(10)
          dst(12) = src(12)
          !
        case (2) ! C2c == (36)(45)(12)
          !
          dst(1) = src(1)
          dst(2) = src(5)
          dst(3) = src(4)
          dst(4) = src(3)
          dst(5) = src(2)
          dst(6) = src(9)
          dst(7) = src(8)
          dst(8) = src(7)
          dst(9) = src(6)
          dst(10) = -src(11)
          dst(11) = -src(10)
          dst(12) = src(12)
          !
        case (6) ! sigma_ab == E*
          !
          dst(1) = src(1)
          dst(2) = src(2)
          dst(3) = src(3)
          dst(4) = src(4)
          dst(5) = src(5)
          dst(6) = src(6)
          dst(7) = src(7)
          dst(8) = src(8)
          dst(9) = src(9)
          dst(10) = -src(10)
          dst(11) = -src(11)
          dst(12) = -src(12)
          !
        case (7) ! sigma_ac == (34)(56)*
          !
          dst(1) = src(1)
          dst(2) = src(3)
          dst(3) = src(2)
          dst(4) = src(5)
          dst(5) = src(4)
          dst(6) = src(7)
          dst(7) = src(6)
          dst(8) = src(9)
          dst(9) = src(8)
          dst(10) = src(10)
          dst(11) = src(11)
          dst(12) = -src(12)
          !
        case (8) ! sigma_bc == (35)(46)(12)*
          !
          dst(1) = src(1)
          dst(2) = src(4)
          dst(3) = src(5)
          dst(4) = src(2)
          dst(5) = src(3)
          dst(6) = src(8)
          dst(7) = src(9)
          dst(8) = src(6)
          dst(9) = src(7)
          dst(10) = -src(11)
          dst(11) = -src(10)
          dst(12) = -src(12)
          !
        case (5) ! i == (36)(45)(12)*
          !
          dst(1) = src(1)
          dst(2) = src(5)
          dst(3) = src(4)
          dst(4) = src(3)
          dst(5) = src(2)
          dst(6) = src(9)
          dst(7) = src(8)
          dst(8) = src(7)
          dst(9) = src(6)
          dst(10) = src(11)
          dst(11) = src(10)
          dst(12) = -src(12)
          !
        end select
        !
      end select
      !
    case('C2H4_2BETA_1TAU')
      !
      select case(trim(molec%symmetry))
        !
      case default
        !
        write(out, '(/a,1x,a,1x,a)') &
        'ML_symmetry_transformation_C2H4 error: symmetry =', trim(molec%symmetry), 'is unknown'
        stop
        !
      case('D2H','D2H(M)')
        !
        select case(ioper)
          !
        case default
          !
          write(out, '(/a,1x,i3,1x,a)') &
          'ML_symmetry_transformation_C2H4 error: symmetry operation ', ioper, 'is unknown'
          stop
          !
        case (1) ! E
          !
          dst(1:12) = src(1:12)
          !
        case (2) ! C2a == (34)(56)
          !
          dst(1) = src(1)
          dst(2) = src(3)
          dst(3) = src(2)
          dst(4) = src(5)
          dst(5) = src(4)
          dst(6) = src(7)
          dst(7) = src(6)
          dst(8) = src(9)
          dst(9) = src(8)
          dst(10) =  -src(10)
          dst(11) =  -src(11)
          dst(12) = src(12)
          !
        case (3) ! C2b == (35)(46)(12)
          !
          dst(1) = src(1)
          dst(2) = src(4)
          dst(3) = src(5)
          dst(4) = src(2)
          dst(5) = src(3)
          dst(6) = src(8)
          dst(7) = src(9)
          dst(8) = src(6)
          dst(9) = src(7)
          dst(10) = -src(11)
          dst(11) = -src(10)
          dst(12) = src(12)
          !
        case (4) ! C2c == (36)(45)(12)
          !
          dst(1) = src(1)
          dst(2) = src(5)
          dst(3) = src(4)
          dst(4) = src(3)
          dst(5) = src(2)
          dst(6) = src(9)
          dst(7) = src(8)
          dst(8) = src(7)
          dst(9) = src(6)
          dst(10) = src(11)
          dst(11) = src(10)
          dst(12) = src(12)
          !
        case (5) ! sigma_ab == E*
          !
          dst(1) = src(1)
          dst(2) = src(2)
          dst(3) = src(3)
          dst(4) = src(4)
          dst(5) = src(5)
          dst(6) = src(6)
          dst(7) = src(7)
          dst(8) = src(8)
          dst(9) = src(9)
          dst(10) = -src(10)
          dst(11) = -src(11)
          dst(12) = -src(12)
          !
        case (6) ! sigma_ac == (34)(56)*
          !
          dst(1) = src(1)
          dst(2) = src(3)
          dst(3) = src(2)
          dst(4) = src(5)
          dst(5) = src(4)
          dst(6) = src(7)
          dst(7) = src(6)
          dst(8) = src(9)
          dst(9) = src(8)
          dst(10) = src(10)
          dst(11) = src(11)
          dst(12) = -src(12)
          !
        case (7) ! sigma_bc == (35)(46)(12)*
          !
          dst(1) = src(1)
          dst(2) = src(4)
          dst(3) = src(5)
          dst(4) = src(2)
          dst(5) = src(3)
          dst(6) = src(8)
          dst(7) = src(9)
          dst(8) = src(6)
          dst(9) = src(7)
          dst(10) = src(11)
          dst(11) = src(10)
          dst(12) = -src(12)
          !
        case (8) ! i == (36)(45)(12)*
          !
          dst(1) = src(1)
          dst(2) = src(5)
          dst(3) = src(4)
          dst(4) = src(3)
          dst(5) = src(2)
          dst(6) = src(9)
          dst(7) = src(8)
          dst(8) = src(7)
          dst(9) = src(6)
          dst(10) = -src(11)
          dst(11) = -src(10)
          dst(12) = -src(12)
          !
        end select
        !
      end select
      !
    case('R_ALPHA_4TAU')
      !
      select case(trim(molec%symmetry))
        !
      case default
        !
        write(out, '(/a,1x,a,1x,a)') &
        'ML_symmetry_transformation_C2H4 error: symmetry =', trim(molec%symmetry), 'is unknown'
        stop 'ML_symmetry_transformation_C2H4 error: bad symmetry type'
        !
      case('D2H','D2H(M)')
        !
        select case(ioper)
          !
        case default
          !
          write(out, '(/a,1x,i3,1x,a)') &
          'ML_symmetry_transformation_C2H4 error: symmetry operation ', ioper, 'is unknown'
          stop 'ML_symmetry_transformation_C2H4 error: bad symmetry operation type'
          !
        case (1) ! E
          !
          dst(1:12) = src(1:12)
          !
        case (2) ! C2a == (35)(46)
          !
          dst(1) = src(1)
          dst(2) = src(4)
          dst(3) = src(5)
          dst(4) = src(2)
          dst(5) = src(3)
          dst(6) = src(8)
          dst(7) = src(9)
          dst(8) = src(6)
          dst(9) = src(7)
          dst(10) =  src(11)
          dst(11) =  src(10)
          dst(12) = -src(12)
          !
        case (3) ! C2b == (34)56)(12)
          !
          dst(1) = src(1)
          dst(2) = src(3)
          dst(3) = src(2)
          dst(4) = src(5)
          dst(5) = src(4)
          dst(6) = src(7)
          dst(7) = src(6)
          dst(8) = src(9)
          dst(9) = src(8)
          dst(10) =-src(10)
          dst(11) =-src(11)
          dst(12) = src(12)
          !
        case (4) ! C2c == (36)(45)(12)
          !
          dst(1) = src(1)
          dst(2) = src(5)
          dst(3) = src(4)
          dst(4) = src(3)
          dst(5) = src(2)
          dst(6) = src(9)
          dst(7) = src(8)
          dst(8) = src(7)
          dst(9) = src(6)
          dst(10) =-src(11)
          dst(11) =-src(10)
          dst(12) =-src(12)
          !
        case (8) ! sigma_ab == E*
          !
          dst(1) = src(1)
          dst(2) = src(5)
          dst(3) = src(4)
          dst(4) = src(3)
          dst(5) = src(2)
          dst(6) = src(9)
          dst(7) = src(8)
          dst(8) = src(7)
          dst(9) = src(6)
          dst(10) = src(11)
          dst(11) = src(10)
          dst(12) = src(12)
          !
        case (7) ! sigma_ac == (35)(46)*
          !
          dst(1) = src(1)
          dst(2) = src(3)
          dst(3) = src(2)
          dst(4) = src(5)
          dst(5) = src(4)
          dst(6) = src(7)
          dst(7) = src(6)
          dst(8) = src(9)
          dst(9) = src(8)
          dst(10) = src(10)
          dst(11) = src(11)
          dst(12) =-src(12)
          !
        case (6) ! sigma_bc == (34)(56)(12)*
          !
          dst(1) = src(1)
          dst(2) = src(4)
          dst(3) = src(5)
          dst(4) = src(2)
          dst(5) = src(3)
          dst(6) = src(8)
          dst(7) = src(9)
          dst(8) = src(6)
          dst(9) = src(7)
          dst(10) = -src(11)
          dst(11) = -src(10)
          dst(12) =  src(12)
          !
        case (5) ! i == (36)(45)(12)*
          !
          dst(1) = src(1)
          dst(2) = src(2)
          dst(3) = src(3)
          dst(4) = src(4)
          dst(5) = src(5)
          dst(6) = src(6)
          dst(7) = src(7)
          dst(8) = src(8)
          dst(9) = src(9)
          dst(10) =-src(10)
          dst(11) =-src(11)
          dst(12) =-src(12)
          !
        end select
        !
      case('D2H(S)')
        !
        select case(ioper)
          !
        case default
          !
          write(out, '(/a,1x,i3,1x,a)') &
          'ML_symmetry_transformation_C2H4 error: symmetry operation ', ioper, 'is unknown'
          stop 'ML_symmetry_transformation_C2H4 error: bad symmetry operation type'
          !
        case (1) ! E
          !
          dst(1:12) = src(1:12)
          !
        case (4) ! C2x == (35)(46)
          !
          dst(1) = src(1)
          dst(2) = src(4)
          dst(3) = src(5)
          dst(4) = src(2)
          dst(5) = src(3)
          dst(6) = src(8)
          dst(7) = src(9)
          dst(8) = src(6)
          dst(9) = src(7)
          dst(10) =  src(11)
          dst(11) =  src(10)
          dst(12) = -src(12)
          !
        case (3) ! C2y == (34)56)(12)
          !
          dst(1) = src(1)
          dst(2) = src(3)
          dst(3) = src(2)
          dst(4) = src(5)
          dst(5) = src(4)
          dst(6) = src(7)
          dst(7) = src(6)
          dst(8) = src(9)
          dst(9) = src(8)
          dst(10) = src(10)
          dst(11) = src(11)
          dst(12) =-src(12)
          !
        case (2) ! C2z == (36)(45)(12)
          !
          dst(1) = src(1)
          dst(2) = src(5)
          dst(3) = src(4)
          dst(4) = src(3)
          dst(5) = src(2)
          dst(6) = src(9)
          dst(7) = src(8)
          dst(8) = src(7)
          dst(9) = src(6)
          dst(10) = src(11)
          dst(11) = src(10)
          dst(12) = src(12)
          !
        case (5) ! i == (36)(45)(12)*
          !
          dst(1) = src(1)
          dst(2) = src(5)
          dst(3) = src(4)
          dst(4) = src(3)
          dst(5) = src(2)
          dst(6) = src(9)
          dst(7) = src(8)
          dst(8) = src(7)
          dst(9) = src(6)
          dst(10) =-src(11)
          dst(11) =-src(10)
          dst(12) =-src(12)
          !
        case (8) ! sigma_yz == (35)(46)*
          !
          dst(1) = src(1)
          dst(2) = src(3)
          dst(3) = src(2)
          dst(4) = src(5)
          dst(5) = src(4)
          dst(6) = src(7)
          dst(7) = src(6)
          dst(8) = src(9)
          dst(9) = src(8)
          dst(10) =-src(10)
          dst(11) =-src(11)
          dst(12) = src(12)
          !
        case (6) ! sigma_xy == E*
          !
          dst(1) = src(1)
          dst(2) = src(2)
          dst(3) = src(3)
          dst(4) = src(4)
          dst(5) = src(5)
          dst(6) = src(6)
          dst(7) = src(7)
          dst(8) = src(8)
          dst(9) = src(9)
          dst(10) =-src(10)
          dst(11) =-src(11)
          dst(12) =-src(12)
          !
        case (7) ! sigma_xz == (34)(56)(12)*
          !
          dst(1) = src(1)
          dst(2) = src(4)
          dst(3) = src(5)
          dst(4) = src(2)
          dst(5) = src(3)
          dst(6) = src(8)
          dst(7) = src(9)
          dst(8) = src(6)
          dst(9) = src(7)
          dst(10) = -src(11)
          dst(11) = -src(10)
          dst(12) =  src(12)
          !
        end select
        !
      end select
      !
    end select
    !
    if (verbose>=5) write(out, '(/a)') 'ML_symmetry_transformation_C2H4/end'
    !
  end subroutine ML_symmetry_transformation_C2H4



  subroutine ML_rotsymmetry_C2H4(J,K,tau,gamma,ideg)
    !
    ! Symmetry transformation rules of TROVE rotational functions
    !
    integer(ik),intent(in) :: J,K,tau
    integer(ik),intent(out) :: gamma,ideg
    !
    if (verbose>=5) write(out, '(/a)') 'ML_rotsymmetry_C2H4/start'
    !
    select case(trim(molec%coords_transform))
      !
      !
    case default
      !
      write(out, '(/a,1x,a,1x,a)') &
      'ML_rotsymmetry_C2H4 error: coordinate type =', trim(molec%coords_transform), 'is unknown'
      stop 'ML_rotsymmetry_C2H4 error: bad coordinate type'
      !
      !
    case('ZMAT_2BETA_1TAU','C2H4_2BETA_1TAU')
      !
      select case(trim(molec%symmetry))
        !
      case default
        !
        write(out, '(/a,1x,a,1x,a)') &
        'ML_rotsymmetry_C2H4 error: symmetry =', trim(molec%symmetry), 'is unknown'
        stop 'ML_rotsymmetry_C2H4 error: bad symmetry type'
        !
      case('D2H(M)')
        !
        gamma = 0
        ideg = 1
        !
        if (mod(K+2,2)==0.and.tau==0) gamma = 1 ! Ag
        if (mod(K+2,2)==0.and.tau==1) gamma = 7 ! B1g
        if (mod(K+2,2)/=0.and.tau==1) gamma = 3 ! B3g
        if (mod(K+2,2)/=0.and.tau==0) gamma = 5 ! B2g
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
        'ML_rotsymmetry_C2H4 error: symmetry =', trim(molec%symmetry), 'is unknown'
        stop 'ML_rotsymmetry_C2H4 error: bad symmetry type'
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
      !
    end select
    !
    if (verbose>=5) write(out, '(/a)') 'ML_rotsymmetry_C2H4/end'
    !
  end subroutine ML_rotsymmetry_C2H4


end module mol_c2h4
