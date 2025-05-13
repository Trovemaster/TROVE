! This unit defines all specific routines for a eight-atomic molecule of ethane-type

module mol_c3h6
  use accuracy
  use moltype

  implicit none

  public ML_coordinate_transform_C3H6, ML_b0_C3H6, ML_symmetry_transformation_C3H6, ML_rotsymmetry_C3H6
  private

  integer(ik), parameter :: verbose = 4 ! Verbosity level

  !--------------------------
  !      ZMAT_prop1
  !--------------------------
  !
  !ZMAT
  !C   0  0  0  0  12.00000000   CH3 C
  !C   1  0  0  0  12.00000000   CH  C
  !C   2  1  0  0  12.00000000   CH2 C
  !H   1  2  3  2   1.00782505   CH3 H in plane
  !H   1  2  4  2   1.00782505   CH3 H
  !H   1  2  4  2   1.00782505   CH3 H
  !H   3  2  1  2   1.00782505   CH2 H   0 deg from CH3 group (trans)
  !H   3  2  7  2   1.00782505   CH2 H   0 deg from CH3 group (cis)
  !H   2  1  7  2   1.00782505   CH  H 
  !end
  ! 1  r_12        0         1.5022
  ! 2  r_23        0         1.3380
  ! 3  r_14        0         1.0914
  ! 4  r_15        0         1.09331
  ! 5  r_16        0         1.09331
  ! 6  r_37        0         1.0827
  ! 7  r_38        0         1.0847
  ! 8  r_26        0         1.0866
  ! 9  alp_321     0         124.4552
  !10  alp_214     0         110.9662
  !11  alp_215     0         110.9662
  !12  alp_216     0         110.9662
  !13  alp_238     0         121.1556
  !14  alp_239     0         121.4643
  !15  alp_127     0         118.8949
  !16  di_CCCH     0         0.0
  !17  di_HCCH     0         120.0000
  !18  di_HCCH     0         240.0000
  !19  di_CCCH     0         0.0
  !20  di_HCCH     0         180.0000
  !21  di_HCCHg    0         180.0000
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
    real(ark) :: thet78,thet79,thet89,tbar,A1,A2,theta12,theta13,theta23
    real(ark) :: T1,T2,T3,tol,T1_

    tol = 1.0d-11

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
    case('ZMAT_7ALF_5TAU','C3H6_7ALF_5TAU')
      !
      if (direct) then ! transform from Z-matrix coords to TROVE coords
        !
        !FOR STRETCHES AND 'ALPHA' BENDS JUST SUBTRACT EQUILIBRIUM COORDINATES
        dst(1:15) = src(1:15)-molec%local_eq(1:15)
        !
        !FOR CCCH 'C1 TO H5' DIHEDRAL ANGLE, ALSO RIGID, TRY JUST SUBTRACTING EQUILIBRIUM
          dst(16) = src(19) - molec%local_eq(19)
        !MAKE SURE THAT ANGLE IS AROUND 0
        if( abs( dst(16)-2.0_ark*pi).lt.1d-4 ) then
          dst(16) = dst(16) - 2.0_ark*pi
        end if
        !
        !FOR CH2 HCH 'BOOK' DIHEDRAL ANGLE, RIGID, TRY JUST SUBTRACTING EQUILIBRIUM
        dst(17) = src(20) - molec%local_eq(20)
        !
        !FOR CH HCCH DIHEDRAL ANGLE, TRY SUBTRACTING EQUILIBRIUM
        dst(18) = src(21) - molec%local_eq(21)
        !
        !T1 = src(16)-pi
        !
        T1 = src(16)
        !
        !       HERE MAKE 'T' ANGLES BY ADDING T1 TO THETA ANGLES
        !       ALWAYS WANT TO MEASURE ANGLES IN SAME WAY
        !       TO GET CONSISTENT TRANSFORMS, DEFINE RANGES FOR ANGLES ELSE ADD/SUB 2*PI
        !
        !IF(T1.LT.0.0.AND.ABS(T1).GT.TOL) THEN
        ! T1 = T1 + 2.0_ark*PI
        !END IF
        !IF(T1.GT.(2.0_ark*PI).AND.ABS(T1-2.0_ARK*PI).GT.TOL) THEN
        ! T1 = T1 - 2.0_ark*PI
        !END IF
        !
        T2 = src(17) + T1 
        T3 = src(18) + T1 
        !
        !IF(T2.LT.2.0_ark*Pi/3.0_ark.AND.abs(T2-2.0_ark*Pi/3.0_ark).GT.tol) THEN
        !T2 = T2 + 2.0_ark*PI
        !END IF
        !IF(T2.GT.8.0_ark*Pi/3.0_ark.AND.(T2-8.0_ark*Pi/3.0_ark).GT.tol) THEN
        !T2 = T2 - 2.0_ark*PI
        !END IF
        !
        !IF(T3.LT.(4.0_ark*PI/3.0_ark).AND.ABS(T3-4.0_ark*PI/3.0_ark).GT.TOL) THEN
        !T3 = T3 + 2.0_ark*PI
        !END IF
        !IF(T3.GT.(10.0_ark*PI/3.0_ark).AND.ABS(T3-10.0_ark*PI/3.0_ark).GT.TOL) THEN
        !T3 = T3 - 2.0_ark*PI
        !END IF
        !
        T1 = T1 - molec%local_eq(16)
        T2 = T2 - molec%local_eq(17)
        T3 = T3 - molec%local_eq(18)
        !
        !       SUBTRACT EQUILBRIUM THETA VALUES TO MAKE A1/A2 ZERO AT EQUILIBRIUM
        !       AND ENSURE CONSISTENT TRANSFROMS
        !
        A1  = (2.0_ark*T1 - T2 - T3)/sqrt(6.0_ark)   
        A2  = (             T2 - T3)/sqrt(2.0_ark)   
        tbar = (T1 + T2 + T3)/3.0_ark 
        !
        dst(19) = A1 
        dst(20) = A2 
        dst(21) = tbar 
        !
      else !  transform from TROVE coords to Z-matrix coords
        !
        dst(1:15) = src(1:15)+molec%local_eq(1:15)
        !
        dst(19) = src(16) + molec%local_eq(19)
        if( abs( dst(19)-2.0_ark*pi).lt.1d-4 ) then
        dst(19) =  dst(19) - 2.0_ark*pi
        end if
        dst(20) = src(17) + molec%local_eq(20)
        dst(21) = src(18) + molec%local_eq(21)
        !
        A1 = src(19) 
        A2 = src(20) 
        tbar = src(21) 
        !
        T1 = (3.0_ark*tbar + sqrt(6.0_ark)*A1)/3.0_ark
        !
        !T1 = T1 + pi 
        !
        IF(T1.LT.0.0.AND.ABS(T1).GT.TOL) THEN
         T1 = T1 + 2.0_ark*PI
        END IF
        IF(T1.GT.(2.0_ark*PI).AND.ABS(T1-2.0_ARK*PI).GT.TOL) THEN
         T1 = T1 - 2.0_ark*PI
        END IF
        !
        T2 = (-A1*sqrt(3.0_ark)+A2 )/sqrt(2.0_ark)
        T3 = (-A1*sqrt(3.0_ark)-A2 )/sqrt(2.0_ark)
        !
        T1 = tbar-(T2+T3)/3.0_ark
        !
        !IF(T1.LT.0.0.AND.ABS(T1).GT.TOL) THEN
        !  T1 = T1 + 2.0_ark*PI
        !END IF
        !IF(T1.GT.(2.0_ark*PI).AND.ABS(T1-2.0_ARK*PI).GT.TOL) THEN
        !  T1 = T1 - 2.0_ark*PI
        !END IF
        !
        !T2 = (sqrt(6.0_ark)*A1 - sqrt(2.0_ark)*A2 - 2.0_ark*T1)/(-2.0_ark)
        !T3 = T2 - sqrt(2.0_ark)*A2
        !
        dst(16) = T1  +  molec%local_eq(16) 
        dst(17) = T2  +  molec%local_eq(17)   
        dst(18) = T3  +  molec%local_eq(18)
        !
      endif
      !
    case('7ALF_TAU_1RHO')
      !
      if (direct) then ! transform from Z-matrix coords to TROVE coords

        !FOR STRETCHES AND 'ALPHA' BENDS JUST SUBTRACT EQUILIBRIUM COORDINATES
        dst(1:15) = src(1:15)-molec%local_eq(1:15)

        !FOR CCCH 'C1 TO H5' DIHEDRAL ANGLE, ALSO RIGID, TRY JUST SUBTRACTING EQUILIBRIUM
        dst(16) = src(19) - molec%local_eq(19)
        !MAKE SURE THAT ANGLE IS AROUND 0
        if( abs( dst(16)-2.0_ark*pi).lt.1d-4 ) then
          dst(16) = dst(16) - 2.0_ark*pi
        end if
        !
        !FOR CH2 HCH 'BOOK' DIHEDRAL ANGLE, RIGID, TRY JUST SUBTRACTING EQUILIBRIUM
        dst(17) = src(20) - molec%local_eq(20)
        !
        !FOR CH HCCH DIHEDRAL ANGLE, TRY SUBTRACTING EQUILIBRIUM
        dst(18) = src(21) - molec%local_eq(21)

        !T1 = src(16)-pi
        !
        T1 = src(16)

!       HERE MAKE 'T' ANGLES BY ADDING T1 TO THETA ANGLES

!       ALWAYS WANT TO MEASURE ANGLES IN SAME WAY
!       TO GET CONSISTENT TRANSFORMS, DEFINE RANGES FOR ANGLES ELSE ADD/SUB 2*PI
        !
        !IF(T1.LT.0.0.AND.ABS(T1).GT.TOL) THEN
        !  T1 = T1 + 2.0_ark*PI
        !END IF
        !IF(T1.GT.(2.0_ark*PI).AND.ABS(T1-2.0_ARK*PI).GT.TOL) THEN
        !  T1 = T1 - 2.0_ark*PI
        !END IF
        !
        !T2 = src(17) + T1 
        !T3 = src(18) + T1 
        !
!        IF(T2.LT.2.0_ark*Pi/3.0_ark.AND.abs(T2-2.0_ark*Pi/3.0_ark).GT.tol) THEN
!        T2 = T2 + 2.0_ark*PI
!        END IF
!        IF(T2.GT.8.0_ark*Pi/3.0_ark.AND.(T2-8.0_ark*Pi/3.0_ark).GT.tol) THEN
!        T2 = T2 - 2.0_ark*PI
!        END IF

!        IF(T3.LT.(4.0_ark*PI/3.0_ark).AND.ABS(T3-4.0_ark*PI/3.0_ark).GT.TOL) THEN
!        T3 = T3 + 2.0_ark*PI
!        END IF
!        IF(T3.GT.(10.0_ark*PI/3.0_ark).AND.ABS(T3-10.0_ark*PI/3.0_ark).GT.TOL) THEN
!        T3 = T3 - 2.0_ark*PI
!        END IF
        !
        !       SUBTRACT EQUILBRIUM THETA VALUES TO MAKE A1/A2 ZERO AT EQUILIBRIUM
        !       AND ENSURE CONSISTENT TRANSFROMS
        !
        T2 = src(17)
        T3 = src(18)
        !
        T1 = T1 - molec%local_eq(16)
        T2 = T2 - molec%local_eq(17)
        T3 = T3 - molec%local_eq(18)
        !
        A1  = (T2 + T3)*3.0_ark/sqrt(6.0_ark)   
        A2  = (T2 - T3)/sqrt(2.0_ark)   
        !
        !tbar = (src(16) + src(17) + src(18))/3.0_ark 
        !
        T2 = T2 + T1 
        T3 =-T3 + T1 
        !
        tbar = (T1 + T2 + T3)/3.0_ark 
        !
        !tbar = (T1 + T2 + T3)/3.0_ark 
        !
        dst(19) = A1 
        dst(20) = A2 
        dst(21) = tbar 
        !
      else !  transform from TROVE coords to Z-matrix coords
        !
        dst(1:15) = src(1:15)+molec%local_eq(1:15)
        !
        dst(19) = src(16) + molec%local_eq(19)
        if( abs( dst(19)-2.0_ark*pi).lt.1d-4 ) then
          dst(19) =  dst(19) - 2.0_ark*pi
        end if
        !
        dst(20) = src(17) + molec%local_eq(20)
        dst(21) = src(18) + molec%local_eq(21)
        !
        A1 = src(19) 
        A2 = src(20) 
        tbar = src(21) 
        !
        !T1 = (3.0_ark*tbar + sqrt(6.0_ark)*A1)/3.0_ark
        !
        !T1 = tbar-sqrt(2.0_ark)/3.0_ark*A2
        !
        !T1 = T1 + pi 
        !
        !T2 = (sqrt(6.0_ark)*A1 - sqrt(2.0_ark)*A2 - 2.0_ark*T1)/(-2.0_ark)
        !T3 = T2 - sqrt(2.0_ark)*A2
        !
        T2 = ( A1/sqrt(3.0_ark)+A2 )/sqrt(2.0_ark)
        T3 = ( A1/sqrt(3.0_ark)-A2 )/sqrt(2.0_ark)
        !
        T1 = tbar-(T2-T3)/3.0_ark
        !
        !IF(T1.LT.0.0.AND.ABS(T1).GT.TOL) THEN
        ! T1 = T1 + 2.0_ark*PI
        !END IF
        !IF(T1.GT.(2.0_ark*PI).AND.ABS(T1-2.0_ARK*PI).GT.TOL) THEN
        ! T1 = T1 - 2.0_ark*PI
        !END IF
        !
        dst(16) = T1  +  molec%local_eq(16) 
        dst(17) = T2  +  molec%local_eq(17)   
        dst(18) = T3  +  molec%local_eq(18)
        !
      endif
      !
    case('7ALF_5THETA_1TAU')
      !
      if (direct) then 
        !
        !for stretches and 'alpha' bends just subtract equilibrium coordinates
        dst(1:15) = src(1:15)-molec%local_eq(1:15)
        !
        !for ccch 'c1 to h5' dihedral angle, also rigid, try just subtracting equilibrium
        dst(16) = src(19) - molec%local_eq(19)
        !make sure that angle is around 0
        if( abs( dst(16)-2.0_ark*pi).lt.1d-4 ) then
          dst(16) = dst(16) - 2.0_ark*pi
        end if
        !
        !for ch2 hch 'book' dihedral angle, rigid, try just subtracting equilibrium
        dst(17) = src(20) - molec%local_eq(20)
        !
        !for ch hcch dihedral angle, try subtracting equilibrium
        dst(18) = src(21) - molec%local_eq(21)
        !
        t1 = src(16)
        !
        ! subtract equilbrium theta values to make a1/a2 zero at equilibrium
        ! and ensure consistent transfroms
        !
        t2 = src(17)
        t3 = src(18)
        !
        if (t2-t1<small_) t2 = t2 + 2.0_ark*pi
        if (t3-t2<small_) t3 = t3 + 2.0_ark*pi
        !
        theta12 = mod(t2-t1+2.0_ark*pi,2.0_ark*pi)
        theta23 = mod(t3-t2+2.0_ark*pi,2.0_ark*pi)
        theta13 = mod(t1-t3+2.0_ark*pi,2.0_ark*pi)
        !
        !t1 = t1 - molec%local_eq(16)
        !t2 = t2 - molec%local_eq(17)
        !t3 = t3 - molec%local_eq(18)
        !
        a1  = ( 2.0_ark*theta23 - theta13 - theta12 )/sqrt(6.0_ark)
        a2  = (                   theta13 - theta12 )/sqrt(2.0_ark)
        !
        tbar = (t1 + t2 + t3-2.0_ark*pi)/3.0_ark
        !
        dst(19) = a1
        dst(20) = a2
        dst(21) = tbar 
        !
      else !  transform from TROVE coords to Z-matrix coords
        !
        dst(1:15) = src(1:15)+molec%local_eq(1:15)
        !
        dst(19) = src(16) + molec%local_eq(19)
        if( abs( dst(19)-2.0_ark*pi).lt.1d-4 ) then
          dst(19) =  dst(19) - 2.0_ark*pi
        end if
        !
        dst(20) = src(17) + molec%local_eq(20)
        dst(21) = src(18) + molec%local_eq(21)
        !
        A1 = src(19) 
        A2 = src(20) 
        tbar = src(21) ! + 2.0_ark*pi/3.0_ark
        !
        !T2 = ( A1/sqrt(3.0_ark)+A2 )/sqrt(2.0_ark)
        !T3 = ( A1/sqrt(3.0_ark)-A2 )/sqrt(2.0_ark)
        !
        !T1 = tbar-(T2-T3)/3.0_ark
        !
        !t1 = -2.0_ark/3.0_ark*pi+1.0_ark/6.0_ark*sqrt(2.0_ark)*A2+1.0_ark/6.0_ark*sqrt(6.0_ark)*A1+tbar
        !t2 = -1.0_ark/3.0_ark*sqrt(2.0_ark)*A2+tbar
        !t3 = 1.0_ark/6.0_ark*sqrt(2.0_ark)*A2+2.0_ark/3.0_ark*pi-1.0_ark/6.0_ark*sqrt(6.0_ark)*A1+tbar 
        
        !t1 = tbar+1.0_ark/3*sqrt(2.0_ark)*A2
        !t2 = 2.0_ark/3.0_ark*Pi+tbar-1.0_ark/6.0_ark*sqrt(2.0_ark)*A2-1.0_ark/6.0_ark*sqrt(6.0_ark)*A1 
        !t3 = 4.0_ark/3.0_ark*Pi+tbar+1.0_ark/6.0_ark*sqrt(6.0_ark)*A1-1.0_ark/6.0_ark*sqrt(2.0_ark)*A2
        !
        t1 = tbar+1.0_ark/3.0_ark*sqrt(2.0_ark)*A2 
        t2 = 2.0_ark/3.0_ark*Pi+tbar-1.0_ark/6.0_ark*sqrt(2.0_ark)*A2-1.0_ark/6.0_ark*sqrt(6.0_ark)*A1
        t3 = 4.0_ark/3.0_ark*Pi+tbar+1.0_ark/6.0_ark*sqrt(6.0_ark)*A1-1.0_ark/6.0_ark*sqrt(2.0_ark)*A2
        !
        dst(16) =  mod(t1+4.0_ark*pi,4.0_ark*pi)
        dst(17) =  mod(t2+4.0_ark*pi,4.0_ark*pi)
        dst(18) =  mod(t3+4.0_ark*pi,4.0_ark*pi)
        !
      endif
      !
    end select
    !
    if (verbose>=5) write(out, '(/a)') 'ML_coordinate_transform_C2H6/end'
    !
  end function ML_coordinate_transform_C3H6



  subroutine ML_b0_C3H6(Npoints,Natoms,b0,rho_i,rho_ref,rho_borders)
    !
    ! Reference structure
    !
    integer(ik),intent(in) :: Npoints, Natoms
    real(ark),intent(out) :: b0(Natoms,3,0:Npoints)
    real(ark),intent(in),optional :: rho_i(0:Npoints)
    real(ark),intent(out),optional :: rho_ref
    real(ark),intent(in),optional :: rho_borders(2)  ! rhomim, rhomax - borders
    real(ark) :: rad,thet,r3
    !
    real(ark) :: a0(molec%Natoms,3),CM_shift,tau,alpha0,alpha,theta,r,r12
    real(ark) :: rC1e,rC2e,rH1e,rH2e,rH3e,rH4e,rH5e,rH6e,alpha1e,alpha2e,alpha3e
    real(ark) :: alpha4e,alpha5e,alpha6e,alpha7e,delta1e,delta2e,delta3e,delta4e
    real(ark) :: delta5e,delta6e,transform(3,3)
    integer(ik) :: i, n, iatom, ix, alloc
    real(ark),allocatable    :: db0(:,:,:)
    !
    if (verbose>=5) write(out, '(/a)') 'ML_b0_C3H6/start'
    !
    rad = pi/180.0_ark
    !
    select case(trim(molec%coords_transform))
      !
    case('ZMAT_7ALF_5TAU','C3H6_7ALF_5TAU')
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
      alpha1e    = molec%alphaeq(1)
      alpha2e    = molec%alphaeq(2)
      alpha3e    = molec%alphaeq(3)
      alpha4e    = molec%alphaeq(4)
      alpha5e    = molec%alphaeq(5)
      alpha6e    = molec%alphaeq(6)
      alpha7e    = molec%alphaeq(7)
      !
      delta1e    = molec%taueq(1)
      delta2e    = molec%taueq(2)
      delta3e    = molec%taueq(3)
      delta4e    = molec%taueq(4)
      delta5e    = molec%taueq(5)
      delta6e    = molec%taueq(6)
      a0 = 0
      !
      a0(1,1) = 0.0_ark
      a0(1,2) = 0.0_ark
      a0(1,3) = -rC1e
      !
      a0(2,1) = 0.0_ark
      a0(2,2) = 0.0_ark
      a0(2,3) = 0.0_ark
      !
      a0(3,1) = rC2e*sin(PI - alpha1e)
      a0(3,2) = 0.0_ark  
      a0(3,3) = rC2e*cos(PI - alpha1e)
      !
      a0(4,1) = rH1e*sin(PI - alpha2e)
      a0(4,2) = 0.0_ark  
      a0(4,3) = -rH1e*cos(PI - alpha2e) - rC1e
      !
      a0(5,1) = rH2e*sin(PI - alpha3e)*cos(delta2e)
      a0(5,2) = rH2e*sin(PI - alpha3e)*sin(delta2e)
      a0(5,3) = -rH2e*cos(PI - alpha3e) - rC1e
      !
      a0(6,1) = rH3e*sin(PI - alpha4e)*cos(delta3e)
      a0(6,2) = rH3e*sin(PI - alpha4e)*sin(delta3e)
      a0(6,3) = -rH3e*cos(PI - alpha4e) - rC1e
      !
      r3 = cosine1(rh4e,rc2e,alpha5e)
      thet = cosine2(rc2e,r3,rh4e)
      !
      a0(7,1) = r3*sin(PI + thet - alpha1e)
      a0(7,2) = 0.0_ark
      a0(7,3) = r3*cos(PI + thet - alpha1e)   
      !
      r3 = cosine1(rh5e,rc2e,alpha6e)
      thet = cosine2(rc2e,r3,rh5e)
      !
      a0(8,1) = r3*sin(PI - thet - alpha1e)
      a0(8,2) = 0.0_ark
      a0(8,3) = r3*cos(PI - thet - alpha1e)   
      !
      a0(9,1) = rh6e*SIN(PI - alpha7e)*COS(delta6e)
      a0(9,2) = rh6e*SIN(PI - alpha7e)*SIN(delta6e)
      a0(9,3) = rh6e*cos(pi - alpha7e) 
      !
    case('7ALF_TAU_1RHO')
      !
      !ZMAT
      !C   0  0  0  0  12.00000000   CH3 C
      !C   1  0  0  0  12.00000000   CH  C
      !C   2  1  0  0  12.00000000   CH2 C
      !H   1  2  3  2   1.00782505   CH3 H in plane
      !H   1  2  3  2   1.00782505   CH3 H
      !H   1  2  3  2   1.00782505   CH3 H
      !H   3  2  1  2   1.00782505   CH2 H   0 deg from CH3 group (trans)
      !H   3  2  7  2   1.00782505   CH2 H   0 deg from CH3 group (cis)
      !H   2  1  7  2   1.00782505   CH  H 
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
      alpha1e    = molec%alphaeq(1)
      alpha2e    = molec%alphaeq(2)
      alpha3e    = molec%alphaeq(3)
      alpha4e    = molec%alphaeq(4)
      alpha5e    = molec%alphaeq(5)
      alpha6e    = molec%alphaeq(6)
      alpha7e    = molec%alphaeq(7)
      !
      delta1e    = molec%taueq(1)
      delta2e    = molec%taueq(2)
      delta3e    = molec%taueq(3)
      delta4e    = molec%taueq(4)
      delta5e    = molec%taueq(5)
      delta6e    = molec%taueq(6)
      !
      a0 = 0
      !
      a0(1,1) = 0.0_ark
      a0(1,2) = 0.0_ark
      a0(1,3) = -rC1e
      !
      a0(2,1) = 0.0_ark
      a0(2,2) = 0.0_ark
      a0(2,3) = 0.0_ark
      !
      a0(3,1) = rC2e*sin(PI - alpha1e)
      a0(3,2) = 0.0_ark  
      a0(3,3) = rC2e*cos(PI - alpha1e)
      !
      a0(4,1) = rH1e*sin(PI - alpha2e)
      a0(4,2) = 0.0_ark  
      a0(4,3) =-rH1e*cos(PI - alpha2e) - rC1e
      !
      a0(5,1) = rH2e*sin(PI - alpha3e)*cos(delta2e)
      a0(5,2) = rH2e*sin(PI - alpha3e)*sin(delta2e)
      a0(5,3) =-rH2e*cos(PI - alpha3e) - rC1e
      !
      a0(6,1) = rH3e*sin(PI - alpha4e)*cos(delta3e)
      a0(6,2) =-rH3e*sin(PI - alpha4e)*sin(delta3e)
      a0(6,3) =-rH3e*cos(PI - alpha4e) - rC1e
      !
      r3 = cosine1(rh4e,rc2e,alpha5e)
      thet = cosine2(rc2e,r3,rh4e)
      !
      a0(7,1) = r3*sin(PI + thet - alpha1e)
      a0(7,2) = 0.0_ark
      a0(7,3) = r3*cos(PI + thet - alpha1e)   
      !
      r3 = cosine1(rh5e,rc2e,alpha6e)
      thet = cosine2(rc2e,r3,rh5e)
      !
      a0(8,1) = r3*sin(PI - thet - alpha1e)
      a0(8,2) = 0.0_ark
      a0(8,3) = r3*cos(PI - thet - alpha1e)   
      !
      a0(9,1) = rh6e*SIN(PI - alpha7e)*COS(delta6e)
      a0(9,2) = rh6e*SIN(PI - alpha7e)*SIN(delta6e)
      a0(9,3) = rh6e*cos(pi - alpha7e) 
      !
      !
    case('7ALF_5THETA_1TAU')
      !
      !ZMAT
      !C   0  0  0  0  12.00000000
      !C   1  0  0  0  12.00000000
      !C   2  1  0  0  12.00000000
      !H   1  2  3  202 1.00782505
      !H   1  2  3  202 1.00782505
      !H   1  2  3  202 1.00782505
      !H   3  2  1  202 1.00782505
      !H   3  2  7  202 1.00782505
      !H   2  1  7  202 1.00782505
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
      alpha1e    = molec%alphaeq(1)
      alpha2e    = molec%alphaeq(2)
      alpha3e    = molec%alphaeq(3)
      alpha4e    = molec%alphaeq(4)
      alpha5e    = molec%alphaeq(5)
      alpha6e    = molec%alphaeq(6)
      alpha7e    = molec%alphaeq(7)
      !
      delta1e    = molec%taueq(1)
      delta2e    = molec%taueq(2)
      delta3e    = molec%taueq(3)
      delta4e    = molec%taueq(4)
      delta5e    = molec%taueq(5)
      delta6e    = molec%taueq(6)
      !
      a0 = 0
      !
      a0(1,1) = 0.0_ark
      a0(1,2) = 0.0_ark
      a0(1,3) = -rC1e
      !
      a0(2,1) = 0.0_ark
      a0(2,2) = 0.0_ark
      a0(2,3) = 0.0_ark
      !
      a0(3,1) = rC2e*sin(PI - alpha1e)
      a0(3,2) = 0.0_ark  
      a0(3,3) = rC2e*cos(PI - alpha1e)
      !
      a0(4,1) = rH1e*sin(PI - alpha2e)
      a0(4,2) = 0.0_ark  
      a0(4,3) =-rH1e*cos(PI - alpha2e) - rC1e
      !
      a0(5,1) = rH2e*sin(PI - alpha3e)*cos(delta2e)
      a0(5,2) = rH2e*sin(PI - alpha3e)*sin(delta2e)
      a0(5,3) =-rH2e*cos(PI - alpha3e) - rC1e
      !
      a0(6,1) = rH3e*sin(PI - alpha4e)*cos(delta3e)
      a0(6,2) = rH3e*sin(PI - alpha4e)*sin(delta3e)
      a0(6,3) =-rH3e*cos(PI - alpha4e) - rC1e
      !
      r3 = cosine1(rh4e,rc2e,alpha5e)
      thet = cosine2(rc2e,r3,rh4e)
      !
      a0(7,1) = r3*sin(PI + thet - alpha1e)
      a0(7,2) = 0.0_ark
      a0(7,3) = r3*cos(PI + thet - alpha1e)   
      !
      r3 = cosine1(rh5e,rc2e,alpha6e)
      thet = cosine2(rc2e,r3,rh5e)
      !
      a0(8,1) = r3*sin(PI - thet - alpha1e)
      a0(8,2) = 0.0_ark
      a0(8,3) = r3*cos(PI - thet - alpha1e)   
      !
      a0(9,1) = rh6e*SIN(pi - alpha7e)*COS(delta6e)
      a0(9,2) = rh6e*SIN(pi - alpha7e)*SIN(delta6e)
      a0(9,3) = rh6e*cos(pi - alpha7e) 
      !
    end select  
    !
    do ix=1, 3
      CM_shift = sum(a0(:,ix)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
      a0(:,ix) = a0(:,ix) - CM_shift
    enddo
    !
    if (verbose>=4) then
      !
      write(out,*) Natoms
      write(out,*)
      write(out, '(1x,a,1x,3(1x,es16.8))') 'C', a0(1,1:3)
      write(out, '(1x,a,1x,3(1x,es16.8))') 'C', a0(2,1:3)
      write(out, '(1x,a,1x,3(1x,es16.8))') 'C', a0(3,1:3)
      do iatom=4, Natoms
        write(out, '(1x,a,1x,3(1x,es16.8))') 'H', a0(iatom,1:3)
      enddo
      !
    endif
    !
    b0(:,:,0) = a0(:,:)
    !
    if (Npoints/=0) then
      !
      if (.not.present(rho_borders).or..not.present(rho_ref)) then
        write(out, '(/a)') 'ML_b0_C3H6: rho_borders and rho_ref must be presented if Npoints > 0'
        stop 'ML_b0_C3H6: rho_borders or rho_ref is not specified'
      endif
      !
      !rho_ref = pi/3.0_ark
      !
      select case(trim(molec%coords_transform))
        !
      case('ZMAT_7ALF_5TAU','C3H6_7ALF_5TAU')
        !
        rho_ref = 0
        !
        do i=0, npoints
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
          alpha1e    = molec%alphaeq(1)
          alpha2e    = molec%alphaeq(2)
          alpha3e    = molec%alphaeq(3)
          alpha4e    = molec%alphaeq(4)
          alpha5e    = molec%alphaeq(5)
          alpha6e    = molec%alphaeq(6)
          alpha7e    = molec%alphaeq(7)
          !
          delta1e    = rho_i(i)
          delta2e    = molec%taueq(2) + delta1e
          delta3e    = molec%taueq(3) + delta1e
          delta4e    = molec%taueq(4)
          delta5e    = molec%taueq(5)
          delta6e    = molec%taueq(6)
          !
          b0(1,1,i) = 0.0_ark
          b0(1,2,i) = 0.0_ark
          b0(1,3,i) = -rC1e
          !     
          b0(2,1,i) = 0.0_ark
          b0(2,2,i) = 0.0_ark
          b0(2,3,i) = 0.0_ark
          !
          b0(3,1,i) = rC2e*sin(PI - alpha1e)
          b0(3,2,i) = 0.0_ark  
          b0(3,3,i) = rC2e*cos(PI - alpha1e)
          !      
          b0(4,1,i) = rh1e*sin(PI - alpha2e)*cos(delta1e)
          b0(4,2,i) = rh1e*sin(PI - alpha2e)*sin(delta1e)
          b0(4,3,i) = -rH1e*cos(PI - alpha2e) - rC1e
          !      
          b0(5,1,i) = rH2e*sin(PI - alpha3e)*cos(delta2e)
          b0(5,2,i) = rH2e*sin(PI - alpha3e)*sin(delta2e)
          b0(5,3,i) = -rH2e*cos(PI - alpha3e) - rC1e
          !     
          b0(6,1,i) = rH3e*sin(PI - alpha4e)*cos(delta3e)
          b0(6,2,i) = rH3e*sin(PI - alpha4e)*sin(delta3e)
          b0(6,3,i) = -rH3e*cos(PI - alpha4e) - rC1e
          !     
          r3 = cosine1(rh4e,rc2e,alpha5e)
          thet = cosine2(rc2e,r3,rh4e)
          !        
          b0(7,1,i) = r3*sin(pi + thet - alpha1e)
          b0(7,2,i) = 0.0_ark
          b0(7,3,i) = r3*cos(pi + thet - alpha1e)   
          !
          r3 = cosine1(rh5e,rc2e,alpha6e)
          thet = cosine2(rc2e,r3,rh5e)
          !
          b0(8,1,i) = r3*sin(pi - thet - alpha1e)
          b0(8,2,i) = 0.0_ark
          b0(8,3,i) = r3*cos(pi - thet - alpha1e)   
          !
          b0(9,1,i) = rh6e*sin(pi - alpha7e)*cos(delta6e)
          b0(9,2,i) = rh6e*sin(pi - alpha7e)*sin(delta6e)
          b0(9,3,i) = rh6e*cos(pi - alpha7e) 
          !
        enddo
        !
      case('7ALF_TAU_1RHO')
        !
        do i=0, npoints
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
          alpha1e    = molec%alphaeq(1)
          alpha2e    = molec%alphaeq(2)
          alpha3e    = molec%alphaeq(3)
          alpha4e    = molec%alphaeq(4)
          alpha5e    = molec%alphaeq(5)
          alpha6e    = molec%alphaeq(6)
          alpha7e    = molec%alphaeq(7)
          !
          delta1e    = rho_i(i)
          delta2e    = molec%taueq(2) + delta1e
          delta3e    =-molec%taueq(3) + delta1e
          delta4e    = molec%taueq(4)
          delta5e    = molec%taueq(5)
          delta6e    = molec%taueq(6)
          !
          b0(1,1,i) = 0.0_ark
          b0(1,2,i) = 0.0_ark
          b0(1,3,i) = -rC1e
          !     
          b0(2,1,i) = 0.0_ark
          b0(2,2,i) = 0.0_ark
          b0(2,3,i) = 0.0_ark
          !
          b0(3,1,i) = rC2e*sin(PI - alpha1e)
          b0(3,2,i) = 0.0_ark  
          b0(3,3,i) = rC2e*cos(PI - alpha1e)
          !      
          b0(4,1,i) = rh1e*sin(PI - alpha2e)*cos(delta1e)
          b0(4,2,i) = rh1e*sin(PI - alpha2e)*sin(delta1e)
          b0(4,3,i) = -rH1e*cos(PI - alpha2e) - rC1e
          !      
          b0(5,1,i) = rH2e*sin(PI - alpha3e)*cos(delta2e)
          b0(5,2,i) = rH2e*sin(PI - alpha3e)*sin(delta2e)
          b0(5,3,i) = -rH2e*cos(PI - alpha3e) - rC1e
          !     
          b0(6,1,i) = rH3e*sin(PI - alpha4e)*cos(delta3e)
          b0(6,2,i) = rH3e*sin(PI - alpha4e)*sin(delta3e)
          b0(6,3,i) = -rH3e*cos(PI - alpha4e) - rC1e
          !     
          r3 = cosine1(rh4e,rc2e,alpha5e)
          thet = cosine2(rc2e,r3,rh4e)
          !        
          b0(7,1,i) = r3*sin(pi + thet - alpha1e)
          b0(7,2,i) = 0.0_ark
          b0(7,3,i) = r3*cos(pi + thet - alpha1e)   
          !
          r3 = cosine1(rh5e,rc2e,alpha6e)
          thet = cosine2(rc2e,r3,rh5e)
          !
          b0(8,1,i) = r3*sin(pi - thet - alpha1e)
          b0(8,2,i) = 0.0_ark
          b0(8,3,i) = r3*cos(pi - thet - alpha1e)   
          !
          b0(9,1,i) = rh6e*sin(pi - alpha7e)*cos(delta6e)
          b0(9,2,i) = rh6e*sin(pi - alpha7e)*sin(delta6e)
          b0(9,3,i) = rh6e*cos(pi - alpha7e) 
          !
        enddo
        !
      case('7ALF_5THETA_1TAU')
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
        alpha1e    = molec%alphaeq(1)
        alpha2e    = molec%alphaeq(2)
        alpha3e    = molec%alphaeq(3)
        alpha4e    = molec%alphaeq(4)
        alpha5e    = molec%alphaeq(5)
        alpha6e    = molec%alphaeq(6)
        alpha7e    = molec%alphaeq(7)
        delta4e    = molec%taueq(4)
        delta5e    = molec%taueq(5)
        delta6e    = molec%taueq(6)
        !
        do i=0, npoints
          !
          delta1e    = rho_i(i)
          delta2e    = molec%taueq(2) + delta1e
          delta3e    = molec%taueq(3) + delta1e
          !
          b0(1,1,i) = 0.0_ark
          b0(1,2,i) = 0.0_ark
          b0(1,3,i) = -rC1e
          !     
          b0(2,1,i) = 0.0_ark
          b0(2,2,i) = 0.0_ark
          b0(2,3,i) = 0.0_ark
          !
          b0(3,1,i) = rC2e*sin(PI - alpha1e)
          b0(3,2,i) = 0.0_ark  
          b0(3,3,i) = rC2e*cos(PI - alpha1e)
          !      
          b0(4,1,i) = rh1e*sin(PI - alpha2e)*cos(delta1e)
          b0(4,2,i) = rh1e*sin(PI - alpha2e)*sin(delta1e)
          b0(4,3,i) = -rH1e*cos(PI - alpha2e) - rC1e
          !      
          b0(5,1,i) = rH2e*sin(PI - alpha3e)*cos(delta2e)
          b0(5,2,i) = rH2e*sin(PI - alpha3e)*sin(delta2e)
          b0(5,3,i) = -rH2e*cos(PI - alpha3e) - rC1e
          !     
          b0(6,1,i) = rH3e*sin(PI - alpha4e)*cos(delta3e)
          b0(6,2,i) = rH3e*sin(PI - alpha4e)*sin(delta3e)
          b0(6,3,i) = -rH3e*cos(PI - alpha4e) - rC1e
          !     
          r3 = cosine1(rh4e,rc2e,alpha5e)
          thet = cosine2(rc2e,r3,rh4e)
          !        
          b0(7,1,i) = r3*sin(pi + thet - alpha1e)
          b0(7,2,i) = 0.0_ark
          b0(7,3,i) = r3*cos(pi + thet - alpha1e)   
          !
          r3 = cosine1(rh5e,rc2e,alpha6e)
          thet = cosine2(rc2e,r3,rh5e)
          !
          b0(8,1,i) = r3*sin(pi - thet - alpha1e)
          b0(8,2,i) = 0.0_ark
          b0(8,3,i) = r3*cos(pi - thet - alpha1e)   
          !
          b0(9,1,i) = rh6e*sin(pi - alpha7e)*cos(delta6e)
          b0(9,2,i) = rh6e*sin(pi - alpha7e)*sin(delta6e)
          b0(9,3,i) = rh6e*cos(pi - alpha7e) 
          !
        enddo
        !
      end select 
      !
      select case(trim(molec%coords_transform))
        !
      case default
        !
        ! standart PAS 
        !
        do i=0, npoints
          !
          ! Saywitz-like (jensen's) reorientation
          !
          call MLorienting_a0(molec%Natoms,molec%AtomMasses,b0(:,:,i))
          !
        enddo
        !
      case('7ALF_5THETA_1TAU-XXX')
        !
        allocate(db0(Natoms,3,0:Npoints),stat=alloc)
        if (alloc/=0) then
            write (out,"(' Error ',i9,' trying to allocate db0-field')") alloc
            stop 'ML_b0_C3H6_db0, db0 - out of memory'
        end if
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
        alpha1e    = molec%alphaeq(1)
        alpha2e    = molec%alphaeq(2)
        alpha3e    = molec%alphaeq(3)
        alpha4e    = molec%alphaeq(4)
        alpha5e    = molec%alphaeq(5)
        alpha6e    = molec%alphaeq(6)
        alpha7e    = molec%alphaeq(7)
        delta4e    = molec%taueq(4)
        delta5e    = molec%taueq(5)
        delta6e    = molec%taueq(6)
        !
        db0 = 0
        !
        do i=0, npoints
          !
          delta1e    = rho_i(i)
          delta2e    = molec%taueq(2) + delta1e
          delta3e    = molec%taueq(3) + delta1e
          !      
          db0(4,1,i) =-rh1e*sin(PI - alpha2e)*sin(delta1e)
          db0(4,2,i) = rh1e*sin(PI - alpha2e)*cos(delta1e)
          !      
          db0(5,1,i) =-rH2e*sin(PI - alpha3e)*sin(delta2e)
          db0(5,2,i) = rH2e*sin(PI - alpha3e)*cos(delta2e)
          !     
          db0(6,1,i) =-rH3e*sin(PI - alpha4e)*sin(delta3e)
          db0(6,2,i) = rH3e*sin(PI - alpha4e)*cos(delta3e)
          !
        enddo
        !
        do i=0, npoints
          !
          ! Saywitz-like (jensen's) reorientation
          !
          !call MLorienting_a0(molec%Natoms,molec%AtomMasses,b0(:,:,i),transform=transform)
          !
          !do ix = 1,Natoms
          !   db0(ix,:,i) = matmul(transpose(transform),db0(ix,:,i))
          !enddo
          !
        enddo
        !
        Call MLorienting_a0_across_dadrho(Natoms,npoints,molec%AtomMasses,rho_borders,b0,db0,periodic=.true.)
        !
        deallocate(db0)
        !
      end select 
      !
      do i=0, npoints
        !
        do n = 1,3
          CM_shift = sum(b0(:,n,i)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
          b0(:,n,i) = b0(:,n,i) - CM_shift
        enddo
        !
        if (verbose>=4) then
         write(out,*) Natoms
         write(out,*)
          write(out, '(1x,a,1x,3(1x,es16.8))') 'C', b0(1,1:3,i)
          write(out, '(1x,a,1x,3(1x,es16.8))') 'C', b0(2,1:3,i)
          write(out, '(1x,a,1x,3(1x,es16.8))') 'C', b0(3,1:3,i)
          do iatom=4, Natoms
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



  subroutine ML_symmetry_transformation_C3H6(ioper, nmodes, src, dst)
    !
    ! Symmetry transformation rules of coordinates
    !
    integer(ik), intent(in)  ::  ioper
    integer(ik), intent(in)  ::  nmodes
    real(ark), intent(in)    ::  src(1:nmodes)
    real(ark), intent(out)   ::  dst(1:nmodes)
    real(ark) :: a,b,p
    !
    a = 0.5_ark
    b = 0.5_ark*sqrt(3.0_ark)
    p = 2.0_ark*pi/3.0_ark
    if (verbose>=5) write(out, '(/a)') 'ML_symmetry_transformation_C3H6/start'
    !
    select case(trim(molec%coords_transform))
      !
    case default
      !
      write(out, '(/a,1x,a,1x,a)') &
      'ML_symmetry_transformation_C3H6 error: coordinate type =', trim(molec%coords_transform), 'is unknown'
      stop 'ML_symmetry_transformation_C3H6 error: bad coordinate type'
      !
    case('ZMAT_7ALF_5TAU','C3H6_7ALF_5TAU')
      !
      select case(trim(molec%symmetry))
        !
      case default
        !
        write(out, '(/a,1x,a,1x,a)') &
        'ML_symmetry_transformation_C3H6 error: symmetry =', trim(molec%symmetry), 'is unknown'
        stop
        !
      case('C3V(M)')
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
          dst(1:21) = src(1:21)
          !
        case (2) ! (123)
          !
          dst(1:2) = src(1:2)

          dst(3) = src(4)
          dst(4) = src(5)
          dst(5) = src(3)

          dst(6:9) = src(6:9)

          dst(10) = src(11)
          dst(11) = src(12)
          dst(12) = src(10)

          dst(13:18) = src(13:18)

          dst(19) = -a*src(19) + b*src(20)
          dst(20) = -b*src(19) - a*src(20)
          !
          dst(21) = mod(src(21) + p,2.0_ark*pi)
          !
        case (3) !(132)
          !
          dst(1:2) = src(1:2)

          dst(3) = src(5)
          dst(4) = src(3)
          dst(5) = src(4)

          dst(6:9) = src(6:9)

          dst(10) = src(12)
          dst(11) = src(10)
          dst(12) = src(11)

          dst(13:18) = src(13:18)

          dst(19) = -a*src(19) - b*src(20)
          dst(20) =  b*src(19) - a*src(20)
          dst(21) = mod(src(21) + 2.0_ark*p,2.0_ark*pi)
          !
        case (4) ! (32)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(3)
          dst(4) = src(5)
          dst(5) = src(4)
          !
          dst(6:9) = src(6:9)
          !
          dst(10) = src(10)
          dst(11) = src(12)
          dst(12) = src(11)
          !
          dst(13:15) = src(13:15)
          dst(16:18) = -src(16:18)
          !
          dst(19) = src(19)
          dst(20) = -src(20)
          dst(21) = 2.0_ark*pi-src(21)
          !
        case (5) ! (12)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(4)
          dst(4) = src(3)
          dst(5) = src(5)
          !
          dst(6:9) = src(6:9)
          !
          dst(10) = src(11)
          dst(11) = src(10)
          dst(12) = src(12)
          !
          dst(13:15) = src(13:15)
          dst(16:18) = -src(16:18)
          !
          dst(19) = -a*src(19) + b*src(20)
          dst(20) =  b*src(19) + a*src(20)
          dst(21) = 2.0_ark*pi-mod(src(21)+p,2.0_ark*pi)
          !
        case (6) ! (13)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(5)
          dst(4) = src(4)
          dst(5) = src(3)
          !
          dst(6:9) = src(6:9)
          !
          dst(10) = src(12)
          dst(11) = src(11)
          dst(12) = src(10)
          !
          dst(13:15) = src(13:15)
          dst(16:18) = -src(16:18)
          !
          dst(19) = -a*src(19) - b*src(20)
          dst(20) = -b*src(19) + a*src(20)
          !
          dst(21) = 2.0_ark*pi-mod(src(21)+2.0_ark*p,2.0_ark*pi)
          !
        end select
          !
      case('C(M)')
          !
        select case(ioper)
          !
        case default
          !
          write(out, '(/a,1x,i3,1x,a)') &
          'ML_symmetry_transformation_C2H6 error: symmetry operation ', ioper, 'is unknown'
          stop
        case (1) ! E
          !
          dst(1:21) = src(1:21)
        end select
        !
      end select 
      !
    case('7ALF_TAU_1RHO','7ALF_5THETA_1TAU')
      !
      select case(trim(molec%symmetry))
        !
      case default
        !
        write(out, '(/a,1x,a,1x,a)') &
        'ML_symmetry_transformation_C3H6 error: symmetry =', trim(molec%symmetry), 'is unknown'
        stop
        !
      case('C3V(M)')
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
          dst(1:21) = src(1:21)
          !
        case (2) ! (123)
          !
          dst(1:2) = src(1:2)

          dst(3) = src(4)
          dst(4) = src(5)
          dst(5) = src(3)

          dst(6:9) = src(6:9)

          dst(10) = src(11)
          dst(11) = src(12)
          dst(12) = src(10)

          dst(13:18) = src(13:18)

          dst(19) = -a*src(19) + b*src(20)
          dst(20) = -b*src(19) - a*src(20)
          !
          dst(21) = mod(src(21) + p,2.0_ark*pi)
          !
        case (3) !(132)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(5)
          dst(4) = src(3)
          dst(5) = src(4)
          !
          dst(6:9) = src(6:9)
          !
          dst(10) = src(12)
          dst(11) = src(10)
          dst(12) = src(11)
          !
          dst(13:18) = src(13:18)
          !
          dst(19) = -a*src(19) - b*src(20)
          dst(20) = +b*src(19) - a*src(20)
          !
          dst(21) = mod(src(21) + 2.0_ark*p,2.0_ark*pi)
          !
        case (4) ! (32)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(3)
          dst(4) = src(5)
          dst(5) = src(4)
          !
          dst(6:9) = src(6:9)
          !
          dst(10) = src(10)
          dst(11) = src(12)
          dst(12) = src(11)
          !
          dst(13:15) = src(13:15)
          dst(16:18) = -src(16:18)
          !
          dst(19) = src(19)
          dst(20) = -src(20)
          dst(21) = 2.0_ark*pi-src(21)
          !
        case (5) ! (12)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(4)
          dst(4) = src(3)
          dst(5) = src(5)
          !
          dst(6:9) = src(6:9)
          !
          dst(10) = src(11)
          dst(11) = src(10)
          dst(12) = src(12)
          !
          dst(13:15) = src(13:15)
          dst(16:18) = -src(16:18)
          !
          dst(19) = -a*src(19) + b*src(20)
          dst(20) = +b*src(19) + a*src(20)
          dst(21) = 2.0_ark*pi-mod(src(21)+p,2.0_ark*pi)
          !
        case (6) ! (13)
          !
          dst(1:2) = src(1:2)
          !
          dst(3) = src(5)
          dst(4) = src(4)
          dst(5) = src(3)
          !
          dst(6:9) = src(6:9)
          !
          dst(10) = src(12)
          dst(11) = src(11)
          dst(12) = src(10)
          !
          dst(13:15) = src(13:15)
          dst(16:18) = -src(16:18)
          !
          dst(19) = -a*src(19) - b*src(20)
          dst(20) = -b*src(19) + a*src(20)
          !
          dst(21) = 2.0_ark*pi-mod(src(21)+2.0_ark*p,2.0_ark*pi)
          !
        end select
        !
      end select
      !
    end select
    !
    if (verbose>=5) write(out, '(/a)') 'ML_symmetry_transformation_C3H6/end'
    !
  end subroutine ML_symmetry_transformation_C3H6



  subroutine ML_rotsymmetry_C3H6(J,K,tau,gamma,ideg)
    !
    ! Symmetry transformation rules of TROVE rotational functions
    !
    integer(ik),intent(in) :: J,K,tau
    integer(ik),intent(out) :: gamma,ideg
    !
    if (verbose>=5) write(out, '(/a)') 'ML_rotsymmetry_C3H6/start'
    !
    select case(trim(molec%coords_transform))
      !
      !
    case default
      !
      write(out, '(/a,1x,a,1x,a)') &
      'ML_rotsymmetry_C3H6 error: coordinate type =', trim(molec%coords_transform), 'is unknown'
      stop 'ML_rotsymmetry_C3H6 error: bad coordinate type'
      !
      !
    case('ZMAT_7ALF_5TAU','C3H6_7ALF_5TAU','7ALF_TAU_1RHO','7ALF_5THETA_1TAU')
      !
      select case(trim(molec%symmetry))
        !
      case default
        !
        write(out, '(/a,1x,a,1x,a)') &
        'ML_rotsymmetry_C2H6 error: symmetry =', trim(molec%symmetry), 'is unknown'
        stop 'ML_rotsymmetry_C3H6 error: bad symmetry type'
        !
      case('C','C(M)')
        !
        gamma = 1
        ideg = 1
        !
      case('C3V','C3V(M)')
        !
        gamma = 0 
        ideg = 1 
        !
        if (mod(tau+2,2)==0) gamma = 1 !; return
        if (mod(tau+2,2)/=0) gamma = 2 !; return
        !
      end select
      !
    end select
    !
    if (verbose>=5) write(out, '(/a)') 'ML_rotsymmetry_C3H6/end'
    !
  end subroutine ML_rotsymmetry_C3H6


  function cosine1(r1,r2,thet) result(r3)
    !
    real(ark) :: r1,r2,thet,r3
    !
    r3 = sqrt(r1**2 + r2**2 - 2.0_ark*r1*r2*cos(thet) )
    !
  end function cosine1
  !
  function cosine2(r1,r2,r3) result(thet)
     !
     real(ark) :: r1,r2,r3,thet
     !
     thet = acos( (r1**2 + r2**2 - r3**2)/(2.0_ark*r1*r2) )
     !
  end function cosine2
  !
end module mol_c3h6
