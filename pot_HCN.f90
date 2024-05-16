!
!  This unit is for a user defined potential 
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
     write(out,"('a,a,a,a')") 'Wrong Potential ',trim(name),'; should be ',trim(poten_name)
   endif
   !
   write(out,"('a')") '  Using USER-tpye PES ',trim(poten_name)
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
   real(ark)              ::  f
   !
   f = MLpoten_xyz_Zobov(ncoords,natoms,local,xyz,force)
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
    ! xyz are undefined for the local case
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
    cos_theta = sum(n1*n2)
    !
    bigr  = r1/bohr
    smallr = r2/bohr
    !
    mu(1)=x_comp_dpl(cos_theta, bigr, smallr)
    mu(3)=z_comp_dpl(cos_theta, bigr, smallr)
    !
    mu(2) = 0
    !
    CN_CM = (x(2,:)*molec%AtomMasses(3))/(sum(molec%AtomMasses(:)))
    !
    u3 = x(1,:)-CN_CM
    !
    u3 = u3 / sqrt(sum(u3(:)**2))
    !
    u2 = (/0._ark,1.0_ark,0._ark/)
    !
    !u2 = n2
    !u2 = u2 / sqrt(sum(u2(:)**2))
    !
    u1 = MLvector_product(u2,u3)
    !
    tmat(1, :) = u1
    tmat(2, :) = u2
    tmat(3, :) = u3
    !
    f(1:3) = matmul(mu,tmat)
    !
    f = f/0.393430_ark ! to debye
    !
  end subroutine MLdipole



  !
  ! Defining potential energy function
  !
  function MLpoten_xyz_Zobov(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force( :)
   real(ark)              ::  f
   !
   real(ark)            :: r1,r2,alpha,xcos,v,v1,v2,v3,v4,v5,v6,v7,v8

   real(ark)            :: aa1,aa2,re1,re2,alphae,xst,xs1,xs2,v0
   integer(ik)          :: N
   real(ark)            :: ycos,v_t,q1,q2
   real(ark)            :: a1(3),a2(3),t1,t2,w(3),cosalpha
   character(len=cl)    :: txt
   !
   if (verbose>=6) write(out,"('MLpoten_xyz_tyuterev/start')")
   !
   r1  = local(1) ;  r2    = local(2) ;  alpha = local(3)
   !
   ! for the second atom = N swap r1(CN) and r2(CH)
   if (molec%AtomMasses(2)>molec%AtomMasses(3)) then
     !
     r1  = local(2) ;  r2    = local(1) ;  alpha = local(3)
     !
   endif
   !
   if (molec%Nangles>0) then
     alphae  = molec%alphaeq(1)
   else
     !
     a1(:) = xyz(2,:) - xyz(1,:)
     a2(:) = xyz(3,:) - xyz(1,:)
     !
     t1 =  sqrt(sum(a1(:)**2))
     t2 =  sqrt(sum(a2(:)**2))
     !
     a1 =  a1(:)/t1
     a2 =  a2(:)/t2
     !
     w(:) = MLvector_product(a1,a2)
     !
     cosalpha = sum(a1(:)*a2(:))
     !
     txt = "pot_xyz"
     !
     alpha = aacos(cosalpha,txt)
     !
     if (alpha<sqrt(small_)) alpha = pi-asin( sqrt( sin(local(3))**2+sin(local(4))**2 ))
     !
   endif 
   !
   re1    = force( 1)
   re2    = force( 2)
   alphae = force( 3)*pi/180.0_ark
   !
   ! calculate potential energy function values
   !
   xs1=(r1-re1)
   xs2=(r2-re2)
   !
   xst = cos(alpha)-cos(alphae)
   !
   N = size(force)
   !
  v0= force(  4) *xs1**0*xs2**0*xst**0
  v1= force(  5) *xs1**1*xs2**0*xst**0& 
     +force(  6) *xs1**0*xs2**1*xst**0& 
     +force(  7) *xs1**0*xs2**0*xst**1& 
     +force(  8) *xs1**2*xs2**0*xst**0& 
     +force(  9) *xs1**1*xs2**1*xst**0& 
     +force( 10) *xs1**1*xs2**0*xst**1& 
     +force( 11) *xs1**0*xs2**2*xst**0& 
     +force( 12) *xs1**0*xs2**1*xst**1& 
     +force( 13) *xs1**0*xs2**0*xst**2&   ! end of 2
     +force( 14) *xs1**3*xs2**0*xst**0& 
     +force( 15) *xs1**2*xs2**1*xst**0& 
     +force( 16) *xs1**2*xs2**0*xst**1& 
     +force( 17) *xs1**1*xs2**2*xst**0& 
     +force( 18) *xs1**1*xs2**1*xst**1& 
     +force( 19) *xs1**1*xs2**0*xst**2& 
     +force( 20) *xs1**0*xs2**3*xst**0& 
     +force( 21) *xs1**0*xs2**2*xst**1& 
     +force( 22) *xs1**0*xs2**1*xst**2& 
     +force( 23) *xs1**0*xs2**0*xst**3&   !  end of 3
     +force( 24) *xs1**4*xs2**0*xst**0& 
     +force( 25) *xs1**3*xs2**1*xst**0& 
     +force( 26) *xs1**3*xs2**0*xst**1& 
     +force( 27) *xs1**2*xs2**2*xst**0& 
     +force( 28) *xs1**2*xs2**1*xst**1& 
     +force( 29) *xs1**2*xs2**0*xst**2& 
     +force( 30) *xs1**1*xs2**3*xst**0& 
     +force( 31) *xs1**1*xs2**2*xst**1&  !
     +force( 32) *xs1**1*xs2**1*xst**2& 
     +force( 33) *xs1**1*xs2**0*xst**3& 
     +force( 34) *xs1**0*xs2**4*xst**0& 
     +force( 35) *xs1**0*xs2**3*xst**1& 
     +force( 36) *xs1**0*xs2**2*xst**2& 
     +force( 37) *xs1**0*xs2**1*xst**3& 
     +force( 38) *xs1**0*xs2**0*xst**4&  !  end of 4
     +force( 39) *xs1**5*xs2**0*xst**0& 
     +force( 40) *xs1**4*xs2**1*xst**0& 
     +force( 41) *xs1**4*xs2**0*xst**1& 
     +force( 42) *xs1**3*xs2**2*xst**0& 
     +force( 43) *xs1**3*xs2**1*xst**1& 
     +force( 44) *xs1**3*xs2**0*xst**2& 
     +force( 45) *xs1**2*xs2**3*xst**0& 
     +force( 46) *xs1**2*xs2**2*xst**1& 
     +force( 47) *xs1**2*xs2**1*xst**2& 
     +force( 48) *xs1**2*xs2**0*xst**3& 
     +force( 49) *xs1**1*xs2**4*xst**0& 
     +force( 50) *xs1**1*xs2**3*xst**1& 
     +force( 51) *xs1**1*xs2**2*xst**2& 
     +force( 52) *xs1**1*xs2**1*xst**3& 
     +force( 53) *xs1**0*xs2**5*xst**0& 
     +force( 54) *xs1**0*xs2**4*xst**1& 
     +force( 55) *xs1**0*xs2**3*xst**2& 
     +force( 56) *xs1**0*xs2**2*xst**3& 
     +force( 57) *xs1**0*xs2**1*xst**4& 
     +force( 58) *xs1**0*xs2**0*xst**5&  !  end of 5
     +force( 59) *xs1**6*xs2**0*xst**0& 
     +force( 60) *xs1**5*xs2**1*xst**0& 
     +force( 61) *xs1**5*xs2**0*xst**1& 
     +force( 62) *xs1**4*xs2**2*xst**0& 
     +force( 63) *xs1**4*xs2**1*xst**1& 
     +force( 64) *xs1**4*xs2**0*xst**2& 
     +force( 65) *xs1**3*xs2**3*xst**0& 
     +force( 66) *xs1**3*xs2**2*xst**1& 
     +force( 67) *xs1**3*xs2**1*xst**2& 
     +force( 68) *xs1**3*xs2**0*xst**3& 
     +force( 69) *xs1**2*xs2**4*xst**0& 
     +force( 70) *xs1**2*xs2**3*xst**1& 
     +force( 71) *xs1**2*xs2**2*xst**2& 
     +force( 72) *xs1**2*xs2**1*xst**3& 
     +force( 73) *xs1**2*xs2**0*xst**4& 
     +force( 74) *xs1**1*xs2**5*xst**0& 
     +force( 75) *xs1**1*xs2**4*xst**1& 
     +force( 76) *xs1**1*xs2**3*xst**2& 
     +force( 77) *xs1**1*xs2**2*xst**3& 
     +force( 78) *xs1**1*xs2**1*xst**4& 
     +force( 79) *xs1**1*xs2**0*xst**5& 
     +force( 80) *xs1**0*xs2**6*xst**0& 
     +force( 81) *xs1**0*xs2**1*xst**5& 
     +force( 82) *xs1**0*xs2**2*xst**4& 
     +force( 83) *xs1**0*xs2**3*xst**3& 
     +force( 84) *xs1**0*xs2**2*xst**4& 
     +force( 85) *xs1**0*xs2**1*xst**5& 
     +force( 86) *xs1**0*xs2**0*xst**6&   !  enf of 6
     +force( 87) *xs1**7*xs2**0*xst**0& 
     +force( 88) *xs1**6*xs2**1*xst**0& 
     +force( 89) *xs1**6*xs2**0*xst**1& 
     +force( 90) *xs1**5*xs2**2*xst**0& 
     +force( 91) *xs1**5*xs2**1*xst**1& 
     +force( 92) *xs1**5*xs2**0*xst**2& 
     +force( 93) *xs1**4*xs2**3*xst**0& 
     +force( 94) *xs1**4*xs2**2*xst**1& 
     +force( 95) *xs1**4*xs2**1*xst**2
  v2= force( 96) *xs1**4*xs2**0*xst**3& 
     +force( 97) *xs1**3*xs2**4*xst**0& 
     +force( 98) *xs1**3*xs2**3*xst**1& 
     +force( 99) *xs1**3*xs2**2*xst**2& 
     +force(100) *xs1**3*xs2**1*xst**3& 
     +force(101) *xs1**3*xs2**0*xst**4& 
     +force(102) *xs1**2*xs2**5*xst**0& 
     +force(103)*xs1**2*xs2**4*xst**1& 
     +force(104)*xs1**2*xs2**3*xst**2& 
     +force(105)*xs1**2*xs2**2*xst**3& 
     +force(106)*xs1**2*xs2**1*xst**4& 
     +force(107)*xs1**2*xs2**0*xst**5& 
     +force(108)*xs1**1*xs2**6*xst**0& 
     +force(109)*xs1**1*xs2**5*xst**1& 
     +force(110)*xs1**1*xs2**4*xst**2& 
     +force(111)*xs1**1*xs2**3*xst**3& 
     +force(112)*xs1**1*xs2**2*xst**4& 
     +force(113)*xs1**1*xs2**1*xst**5& 
     +force(114)*xs1**1*xs2**0*xst**6& 
     +force(115)*xs1**0*xs2**7*xst**0& 
     +force(116)*xs1**0*xs2**6*xst**1& 
     +force(117)*xs1**0*xs2**5*xst**2& 
     +force(118)*xs1**0*xs2**4*xst**3& 
     +force(119)*xs1**0*xs2**3*xst**4& 
     +force(120)*xs1**0*xs2**2*xst**5& 
     +force(121)*xs1**0*xs2**1*xst**6& 
     +force(122)*xs1**0*xs2**0*xst**7& ! end of 7
     +force(123)*xs1**0*xs2**0*xst**8& 
     +force(124)*xs1**6*xs2**2*xst**1& 
     +force(125)*xs1**7*xs2**0*xst**2& 
     +force(126)*xs1**7*xs2**2*xst**0& 
     +force(127)*xs1**8*xs2**0*xst**1& 
     +force(128)*xs1**9*xs2**0*xst**0& 
     +force(129)*xs1**0*xs2**2*xst**8& 
     +force(130)*xs1**0*xs2**4*xst**6& 
     +force(131)*xs1**0*xs2**6*xst**4& 
     +force(132)*xs1**0*xs2**8*xst**2& 
     +force(133)*xs1**1*xs2**0*xst**9& 
     +force(134)*xs1**1*xs2**2*xst**7& 
     +force(135)*xs1**1*xs2**4*xst**5& 
     +force(136)*xs1**1*xs2**6*xst**3& 
     +force(137)*xs1**1*xs2**8*xst**1& 
     +force(138)*xs1**2*xs2**0*xst**8& 
     +force(139)*xs1**2*xs2**2*xst**6& 
     +force(140)*xs1**2*xs2**4*xst**4& 
     +force(141)*xs1**2*xs2**6*xst**2& 
     +force(142)*xs1**3*xs2**0*xst**7& 
     +force(143)*xs1**3*xs2**2*xst**5& 
     +force(144)*xs1**3*xs2**4*xst**3& 
     +force(145)*xs1**3*xs2**6*xst**1& 
     +force(146)*xs1**4*xs2**0*xst**6& 
     +force(147)*xs1**4*xs2**2*xst**4& 
     +force(148)*xs1**4*xs2**4*xst**2& 
     +force(149)*xs1**5*xs2**0*xst**5& 
     +force(150)*xs1**5*xs2**2*xst**3& 
     +force(151)*xs1**5*xs2**4*xst**1& 
     +force(152)*xs1**6*xs2**0*xst**4& 
     +force(153)*xs1**6*xs2**2*xst**2& 
     +force(154)*xs1**6*xs2**4*xst**0& 
     +force(155)*xs1**7*xs2**0*xst**3& 
     +force(156)*xs1**7*xs2**2*xst**1& 
     +force(157)*xs1**8*xs2**0*xst**2& 
     +force(158)*xs1**8*xs2**2*xst**0& 
     +force(159)*xs1**9*xs2**0*xst**1& 
     +force(160)*xs1**10*xs2**0*xst**0& 
     +force(161)*xs1**0*xs2**2*xst**9& 
     +force(162)*xs1**0*xs2**4*xst**7& 
     +force(163)*xs1**0*xs2**6*xst**5& 
     +force(164)*xs1**0*xs2**8*xst**3& 
     +force(165)*xs1**0*xs2**10*xst**1& 
     +force(166)*xs1**1*xs2**0*xst**10& 
     +force(167)*xs1**1*xs2**2*xst**8& 
     +force(168)*xs1**1*xs2**4*xst**6& 
     +force(169)*xs1**1*xs2**6*xst**4& 
     +force(170)*xs1**1*xs2**8*xst**2& 
     +force(171)*xs1**1*xs2**10*xst**0& 
     +force(172)*xs1**2*xs2**0*xst**9& 
     +force(173)*xs1**2*xs2**2*xst**7& 
     +force(174)*xs1**2*xs2**4*xst**5& 
     +force(175)*xs1**2*xs2**6*xst**3& 
     +force(176)*xs1**2*xs2**8*xst**1& 
     +force(177)*xs1**3*xs2**0*xst**8 
  v3= force(178)*xs1**3*xs2**2*xst**6& 
     +force(179)*xs1**3*xs2**4*xst**4& 
     +force(180)*xs1**3*xs2**6*xst**2& 
     +force(181)*xs1**3*xs2**8*xst**0& 
     +force(182)*xs1**4*xs2**0*xst**7& 
     +force(183)*xs1**4*xs2**2*xst**5& 
     +force(184)*xs1**4*xs2**4*xst**3& 
     +force(185)*xs1**4*xs2**6*xst**1& 
     +force(186)*xs1**5*xs2**0*xst**6& 
     +force(187)*xs1**5*xs2**2*xst**4& 
     +force(188)*xs1**5*xs2**4*xst**2& 
     +force(189)*xs1**5*xs2**6*xst**0& 
     +force(190)*xs1**6*xs2**0*xst**5& 
     +force(191)*xs1**6*xs2**2*xst**3& 
     +force(192)*xs1**6*xs2**4*xst**1& 
     +force(193)*xs1**7*xs2**0*xst**4& 
     +force(194)*xs1**7*xs2**2*xst**2& 
     +force(195)*xs1**8*xs2**0*xst**3& 
     +force(196)*xs1**8*xs2**2*xst**1& 
     +force(197)*xs1**9*xs2**0*xst**2& 
     +force(198)*xs1**9*xs2**2*xst**0& 
     +force(199)*xs1**10*xs2**0*xst**1& 
     +force(200)*xs1**11*xs2**0*xst**0& 
     +force(201)*xs1**0*xs2**0*xst**12& 
     +force(202)*xs1**0*xs2**2*xst**10& 
     +force(203)*xs1**0*xs2**4*xst**8& 
     +force(204)*xs1**0*xs2**6*xst**6& 
     +force(205)*xs1**0*xs2**8*xst**4& 
     +force(206)*xs1**0*xs2**10*xst**2& 
     +force(207)*xs1**0*xs2**12*xst**0& 
     +force(208)*xs1**1*xs2**0*xst**11& 
     +force(209)*xs1**1*xs2**2*xst**9& 
     +force(210)*xs1**1*xs2**4*xst**7& 
     +force(211)*xs1**1*xs2**6*xst**5& 
     +force(212)*xs1**1*xs2**8*xst**3& 
     +force(213)*xs1**1*xs2**10*xst**1& 
     +force(214)*xs1**2*xs2**0*xst**10& 
     +force(215)*xs1**2*xs2**2*xst**8& 
     +force(216)*xs1**2*xs2**4*xst**6& 
     +force(217)*xs1**2*xs2**6*xst**4& 
     +force(218)*xs1**2*xs2**8*xst**2& 
     +force(219)*xs1**2*xs2**10*xst**0& 
     +force(220)*xs1**3*xs2**0*xst**9& 
     +force(221)*xs1**3*xs2**2*xst**7& 
     +force(222)*xs1**3*xs2**4*xst**5& 
     +force(223)*xs1**3*xs2**6*xst**3& 
     +force(224)*xs1**3*xs2**8*xst**1& 
     +force(225)*xs1**4*xs2**0*xst**8& 
     +force(226)*xs1**4*xs2**2*xst**6& 
     +force(227)*xs1**4*xs2**4*xst**4& 
     +force(228)*xs1**4*xs2**6*xst**2& 
     +force(229)*xs1**4*xs2**8*xst**0& 
     +force(230)*xs1**5*xs2**0*xst**7& 
     +force(231)*xs1**5*xs2**2*xst**5& 
     +force(232)*xs1**5*xs2**4*xst**3& 
     +force(233)*xs1**5*xs2**6*xst**1& 
     +force(234)*xs1**6*xs2**0*xst**6& 
     +force(235)*xs1**6*xs2**2*xst**4& 
     +force(236)*xs1**6*xs2**4*xst**2& 
     +force(237)*xs1**6*xs2**6*xst**0& 
     +force(238)*xs1**7*xs2**0*xst**5& 
     +force(239)*xs1**7*xs2**2*xst**3& 
     +force(240)*xs1**7*xs2**4*xst**1& 
     +force(241)*xs1**8*xs2**0*xst**4& 
     +force(242)*xs1**8*xs2**2*xst**2& 
     +force(243)*xs1**8*xs2**4*xst**0& 
     +force(244)*xs1**9*xs2**0*xst**3& 
     +force(245)*xs1**9*xs2**2*xst**1& 
     +force(246)*xs1**10*xs2**0*xst**2&  
     +force(247)*xs1**10*xs2**2*xst**0& 
     +force(248)*xs1**11*xs2**0*xst**1& 
     +force(249)*xs1**12*xs2**0*xst**0
    !
    f=v0+v1+v2+v3
    !
    if (verbose>=6) write(out,"('MLpoten_xyz_tyuterev/end')")

 end function MLpoten_xyz_Zobov



  
! This program calculates the dipole of:-
! T. van Mourik, G. J. Harris, O. L. Polyansky, J. Tennyson, 
! A. G. Csaszar and P. J. Knowles, J. Chem. Phys. 115, 3706 (2001).
! Written by by GJH March 16th 2001.
!
! 
! Fits are 3*3*5 and 3*3*7 fits to x and z respectively

!      A test main program
 
!      implicit double precision (a-h,o-z)
!      cosg=cos(1.3962634015952)
!      bigr= 2.16667
!      smallr=2.13333 
!      icomp=0

!      should then give dipd=0.499620862d0 for icomp=1
!      for icomp=0 dipd=-0.00166465088d0


!      smallr=2.2 
!      bigr=3.05172d0
!      cosg=-1.0

!      do 10 i=1, 201

!      icomp=0
!      call dipd(dipole0, smallr, bigr, cosg,  icomp )
!      icomp=1
!      call dipd(dipole1, smallr, bigr, cosg,  icomp )

!      write(6,20) cosg, bigr, smallr, dipole0*2.5417662, dipole1*2.5417662
!      cosg=cosg+0.009999d0

! 10    continue
! 20    format(5d12.4)

!      end


! Jacobi coordinates are used.
! Units of lengths in a0 and dipole in au.
!
! smallr is the C to N bond length. 

! bigr is the H to CN center of mass distance.

! cosg is the cosine of the angle between the bigr and smallr, 
!      an angle of 0 corresponds to a linear molecule of configuration
!      HCN and an angle of Pi is CNH

! icomp selects the x or z component of the dipole, 
!     if icomp=0 then returns z component of dipole
!     if icomp=1 then the x component of the dipole is returned. Dipole is the output 

! Z dipole component is along C to H bond.
! Xdipole component is perpendicular to z component, H bends in xz plane.



! *************************************************************************
      function z_comp_dpl(cosg, bigr, smallr) result (f)

      implicit none
      real(ark),intent(in) :: cosg, bigr, smallr
      real(ark) :: coeffz(63)
      real(ark) :: B(3), C(3) , f , temp,exprbig,exprsmall,apolynom
      integer(ik) :: icoefref,i,j,k


! code to calculate the z component of the dipole moment.
! Standard Deviation = 0.008 D
           

      coeffz(1)=-1039.93729474448763d0
      coeffz(2)=353.933029028493396d0
      coeffz(3)=-469.273752106125326d0
      coeffz(4)=5053.36372859932299d0
      coeffz(5)=-174.791252669703362d0
      coeffz(6)=-1790.40179998532812d0
      coeffz(7)=-1737.30128936150733d0
      coeffz(8)=1036.43863494317905d0
      coeffz(9)=-351.840317729231717d0
      coeffz(10)=466.654729016545031d0
      coeffz(11)=-5035.72026003557988d0
      coeffz(12)=174.888204603534292d0
      coeffz(13)=1783.68882790358092d0
      coeffz(14)=1730.16793139302162d0
      coeffz(15)=4.13996041131474449d0
      coeffz(16)=-1.77817321839038074d0
      coeffz(17)=2.46984600203225505d0
      coeffz(18)=-20.8620301949723831d0
      coeffz(19)=0.248557957243695296d0
      coeffz(20)=7.69363017734038415d0
      coeffz(21)=7.82369848101642729d0
      coeffz(22)=-170.934815933980384d0
      coeffz(23)=1038.00304354550415d0
      coeffz(24)=350.15335220496472d0
      coeffz(25)=-1945.5046107913207d0
      coeffz(26)=-192.883669686284673d0
      coeffz(27)=646.020865093223846d0
      coeffz(28)=464.968280321535981d0
      coeffz(29)=170.495969427259945d0
      coeffz(30)=-1035.59570071766628d0
      coeffz(31)=-348.571134507640448d0
      coeffz(32)=1938.59278642478573d0
      coeffz(33)=191.925410021956726d0
      coeffz(34)=-643.632109747355521d0
      coeffz(35)=-463.047969228713786d0
      coeffz(36)=0.60555784179134213d0
      coeffz(37)=-3.61463889815628707d0
      coeffz(38)=-1.65622524995037699d0
      coeffz(39)=8.06544948486997637d0
      coeffz(40)=0.985792239355321115d0
      coeffz(41)=-2.74472229985483116d0
      coeffz(42)=-2.10370116762671985d0
      coeffz(43)=11.6158010664658819d0
      coeffz(44)=-1.34518644018171973d0
      coeffz(45)=-30.6940205962599487d0
      coeffz(46)=12.6697949751099632d0
      coeffz(47)=20.7241897553344583d0
      coeffz(48)=-4.33083732019109644d0
      coeffz(49)=-9.04791261638706548d0
      coeffz(50)=-11.5726193461914656d0
      coeffz(51)=1.34352911201902514d0
      coeffz(52)=30.5762810875374454d0
      coeffz(53)=-12.6246369560126679d0
      coeffz(54)=-20.6427638916496952d0
      coeffz(55)=4.31545198368715627d0
      coeffz(56)=9.01190335420644576d0
      coeffz(57)=-0.0494259888418515414d0
      coeffz(58)=0.00383241424113694307d0
      coeffz(59)=0.132747818418275613d0
      coeffz(60)=-0.0526642790873265692d0
      coeffz(61)=-0.0909424933648184257d0
      coeffz(62)=0.0179744240449498035d0
      coeffz(63)=0.0400004758413293306d0

      B(1) = -0.33194847405555984d0
      B(2) = 0.26846609653383304d0
      B(3) = 1.4784268008920333d0

      C(1) = 0.41155209135581288d0
      C(2) = 0.41038551093500483d0
      C(3) = 0.58913645613028966d0

      temp=0
      icoefref=1
      
      do i=0, 2
         exprbig=exp(bigR *  B(i+1))
         do j=0, 2
            exprsmall=exp(smallr * C(j+1))
            do k=0, 6
               call plgndr(k, 0, cosg, apolynom)
                  temp=temp+(coeffz(icoefref)*exprbig*exprsmall*apolynom)
               icoefref=icoefref+1
            enddo
         enddo
      enddo
      f = temp
      return
      
      end function z_comp_dpl


! ************************************************************************

      function x_comp_dpl(cosg, bigr, smallr) result (f)


      implicit none
      !
      real(ark),intent(in) :: cosg, bigr, smallr
      real(ark) :: coeffx(45)
      real(ark) :: B(3), C(3) , f, temp,exprbig,exprsmall,apolynom
      integer(ik) :: icoefref,i,j,k


      coeffx(1)=3307.26744745775115d0
      coeffx(2)=31514.9042500608293d0
      coeffx(3)=-20148.594904970726d0
      coeffx(4)=546.081750420754701d0
      coeffx(5)=7723.26526967422282d0
      coeffx(6)=-584.587960975932364d0
      coeffx(7)=-32963.3776253884116d0
      coeffx(8)=21191.7974142975256d0
      coeffx(9)=-590.053047331642632d0
      coeffx(10)=-7990.27903779592029d0
      coeffx(11)=-0.031796628743699731d0
      coeffx(12)=-0.00940828117939172576d0
      coeffx(13)=0.00955414571501609397d0
      coeffx(14)=-0.00072965350671948823d0
      coeffx(15)=-0.00340116044617136415d0
      coeffx(16)=-793.754144943871578d0
      coeffx(17)=-3486.98023111115815d0
      coeffx(18)=2421.86358491171674d0
      coeffx(19)=-146.240552443951966d0
      coeffx(20)=-908.472444995532274d0
      coeffx(21)=521.168738604204326d0
      coeffx(22)=3646.28887952151323d0
      coeffx(23)=-2546.78588888870825d0
      coeffx(24)=154.568740785581287d0
      coeffx(25)=939.567725418272133d0
      coeffx(26)=0.00359417943193839101d0
      coeffx(27)=0.00108249483506098995d0
      coeffx(28)=-0.00108048813729769784d0
      coeffx(29)=0.000123617108271128898d0
      coeffx(30)=0.000402299503589065202d0
      coeffx(31)=-2556.9099117606792d0
      coeffx(32)=-28050.9552366028232d0
      coeffx(33)=17760.1943565608219d0
      coeffx(34)=-413.178836891090924d0
      coeffx(35)=-6829.82090935623486d0
      coeffx(36)=106.994087651491977d0
      coeffx(37)=29341.1779083554478d0
      coeffx(38)=-18680.3377835709173d0
      coeffx(39)=449.374983259753014d0
      coeffx(40)=7066.18946387678389d0
      coeffx(41)=0.0282356308831904958d0
      coeffx(42)=0.0083400647971135194d0
      coeffx(43)=-0.00848293720114831184d0
      coeffx(44)=0.000614212002770376261d0
      coeffx(45)=0.00300582095992624047d0




! Standard Deviation = 0.002 D

      temp=0
      icoefref=1

      B(1) = -0.47194808781916352d0
      B(2) = -0.54677773064856535d0
      B(3) = -0.46378831330917923d0

      C(1) = 0.0060793434403400061d0
      C(2) = -0.011738771795116187d0
      C(3) = 3.9847070586924296d0

      do i=0, 2
         exprbig=exp(bigR *  B(i+1))
         do j=0, 2
            exprsmall=exp(smallr * C(j+1))
            do k=1, 5
               call plgndr(k, 1, cosg, apolynom)
               temp=temp+(coeffx(icoefref)*exprbig*exprsmall*apolynom)
               icoefref=icoefref+1
            enddo
         enddo
      enddo
      f = temp

      return
      end function x_comp_dpl


! ************************************************************************

      subroutine plgndr(l, m, x, pmm)

! subroutine to calculate the asociated legenre polynomials, translated drectly from a C 
! version to fortran version. The C version comes from Numerical recipies in C,  Press et al. 
      implicit none 
      integer(ik),intent(in) :: l,m
      real(ark),intent(in) :: x
      real(ark),intent(out) :: pmm
      real(ark) ::  fact, pll, pmmp1, somx2
      integer(ik) ::  i, ll
      !
      if (m .LT. 0 .OR. m .GT. l .OR. (x*x) .GT. 1.0) then
         write(out,*)  "!ERROR! bad arguments in routine plgndr"
         stop '!ERROR! bad arguments in routine plgndr'
      endif

      pmm = 1.0_ark

      if (m .GT. 0) then
   
         somx2=sqrt((1.0_ark - x) * (1.0_ark + x))
         fact = 1.0_ark
         do i=1, m
            pmm = pmm*(-fact * somx2)
            fact = fact+ 2.0_ark
         enddo
      endif

      if (l .NE. m) then
         pmmp1 = x * real(2 * m + 1,ark) * pmm

         if (l .EQ. (m + 1)) then
            pmm=pmmp1
         else
       
            do ll = m+2, l
           
            pll=(x*real(2*ll-1,ark)*pmmp1-real(ll+m-1,ark)*pmm)/real(ll-m,ark)
            pmm=pmmp1
            pmmp1=pll
            enddo
            pmm=pll
         endif
      endif
      return 
      !
      end subroutine plgndr





end module pot_user
