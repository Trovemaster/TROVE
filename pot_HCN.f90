!
!  This unit is for a user defined potential 
!
module pot_user
  use accuracy
  use moltype

  implicit none

  public MLdipole,MLpoten,ML_MEP

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


 recursive subroutine MLdipole(rank,ncoords,natoms,local,xyz,f)
   !
   integer(ik),intent(in) ::  rank,ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
   real(ark),intent(out)  ::  f(rank)
       !
       f = 0
       !
  end subroutine MLdipole

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



end module pot_user
