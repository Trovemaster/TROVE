! This unit defines all specific routines for a six-atomic molecule of ethylene-type

module pot_c2h4
use accuracy
use moltype

implicit none

public MLpoten_c2h4_88, MLpoten_c2h4_lee

private

integer(ik), parameter :: verbose = 4 ! Verbosity level


contains

function MLpoten_c2h4_88(ncoords, natoms, local, xyz, force) result(f)

  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f

  real(ark) :: rad, y(12), a1, a2, a3, r1, r2, b1_zmat, b2_zmat, tau1_zmat, tau1, tau2, dtau, db1, db2
  real(ark) :: vshort, beta1, beta2, beta12, vlong, de_str2(5), de_str4(5), de_bnd2(7), de_bnd4(7)
  real(ark) :: vdamp, gamma1, gamma2, delta1, delta2, delta3, delta4, delta5

  rad = pi/180.0_ark

  ! expansion functions

  r1 = force(1)
  r2 = force(2)
  a1 = force(3)*rad
  a2 = force(4)*rad
  a3 = force(5)*rad
  beta1 = force(6)
  beta2 = force(7)

  ! long-range function parameters

  gamma1 = force(8)  ! cc stretch
  gamma2 = force(9)  ! ch stretch

  de_str2(1) = force(10)  ! cc stretch
  de_str4(1) = force(11)  ! cc stretch

  de_str2(2:5) = force(12)  ! ch stretch
  de_str4(2:5) = force(13)  ! ch stretch

  de_bnd2(1:4) = force(14)  ! hcc bend
  de_bnd4(1:4) = force(15)  ! hcc bend

  de_bnd2(5:6) = force(16)  ! beta bend
  de_bnd4(5:6) = force(17)  ! beta bend

  de_bnd2(7) = force(18)  ! tau bend
  de_bnd4(7) = force(19)  ! tau bend

  ! damping-function parameters

  delta1 = force(20)
  delta2 = force(21)
  delta3 = force(22)
  delta4 = force(23)
  delta5 = force(24)

  ! expansion functions

  y(1) = 1.0_ark-exp(-beta1*(local(1)-r1))

  y(2) = 1.0_ark-exp(-beta2*(local(2)-r2))
  y(3) = 1.0_ark-exp(-beta2*(local(3)-r2))
  y(4) = 1.0_ark-exp(-beta2*(local(4)-r2))
  y(5) = 1.0_ark-exp(-beta2*(local(5)-r2))

  y(6) = local(6) - a1
  y(7) = local(7) - a1
  y(8) = local(8) - a1
  y(9) = local(9) - a1

  b1_zmat   = local(10)
  tau1_zmat = local(11)
  b2_zmat   = local(12)

  db1 = b1_zmat - pi
  db2 = b2_zmat - pi

  if (tau1_zmat>pi) then
    tau1 = 2.0_ark*pi - tau1_zmat
  elseif (tau1_zmat<pi) then
    tau1 = -tau1_zmat
  endif

  tau2 = db1 + db2 + tau1
  dtau = tau1 + tau2

  y(10:12) = (/db1, db2, dtau/)

  ! short-range potential

  vshort = force(25) &!
         + c2h4_poten_n1_d8( y, force(26:57) ) &!
         + c2h4_poten_n2_d8( y, force(58:433) )

  ! long-range potential

  vlong = sum(de_str2(1:5)*y(1:5)**2 + de_str4(1:5)*y(1:5)**4) &!
        + exp( -gamma1*(local(1)-r1)**2 -gamma2*sum((local(2:5)-r2)**2) ) &!
                * sum( de_bnd2(1:7)*y(6:12)**2 + de_bnd4(1:7)*y(6:12)**4 )

  ! damping function

  vdamp = exp( -delta1*(local(1)-r1)**2        -2.0_ark*delta1*(local(1)-r1)**4 ) &!
        * exp( -delta2*sum((local(2:5)-r2)**2) -2.0_ark*delta2*sum((local(2:5)-r2)**4) ) &!
        * exp( -delta3*sum(y(6:9)**2)           -2.0_ark*delta3*sum(y(6:9)**4) ) &!
        * exp( -delta4*sum(y(10:11)**2)         -2.0_ark*delta4*sum(y(10:11)**4) ) &!
        * exp( -delta5*sum(y(12:12)**2)         -2.0_ark*delta5*sum(y(12:12)**4) )

  f = vshort*vdamp + vlong

end function MLpoten_c2h4_88


!---------------------------------------------------------------------


function c2h4_poten_n1_d8(coords, params) result(v)
real(ark), intent(in) :: coords(12)
real(ark), intent(in) :: params(32)
real(ark) :: v
real(ark) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
real(ark) :: f(32)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = 8*r0
f(2) = 8*r0**2
f(3) = 8*r0**3
f(4) = 8*r0**4
f(5) = 8*r0**5
f(6) = 8*r0**6
f(7) = 8*r0**7
f(8) = 8*r0**8
f(9) = 2*r1 + 2*r2 + 2*r3 + 2*r4
f(10) = 2*r1**2 + 2*r2**2 + 2*r3**2 + 2*r4**2
f(11) = 2*r1**3 + 2*r2**3 + 2*r3**3 + 2*r4**3
f(12) = 2*r1**4 + 2*r2**4 + 2*r3**4 + 2*r4**4
f(13) = 2*r1**5 + 2*r2**5 + 2*r3**5 + 2*r4**5
f(14) = 2*r1**6 + 2*r2**6 + 2*r3**6 + 2*r4**6
f(15) = 2*r1**7 + 2*r2**7 + 2*r3**7 + 2*r4**7
f(16) = 2*r1**8 + 2*r2**8 + 2*r3**8 + 2*r4**8
f(17) = 2*a1 + 2*a2 + 2*a3 + 2*a4
f(18) = 2*a1**2 + 2*a2**2 + 2*a3**2 + 2*a4**2
f(19) = 2*a1**3 + 2*a2**3 + 2*a3**3 + 2*a4**3
f(20) = 2*a1**4 + 2*a2**4 + 2*a3**4 + 2*a4**4
f(21) = 2*a1**5 + 2*a2**5 + 2*a3**5 + 2*a4**5
f(22) = 2*a1**6 + 2*a2**6 + 2*a3**6 + 2*a4**6
f(23) = 2*a1**7 + 2*a2**7 + 2*a3**7 + 2*a4**7
f(24) = 2*a1**8 + 2*a2**8 + 2*a3**8 + 2*a4**8
f(25) = 4*b1**2 + 4*b2**2
f(26) = 4*b1**4 + 4*b2**4
f(27) = 4*b1**6 + 4*b2**6
f(28) = 4*b1**8 + 4*b2**8
f(29) = 8*dtau**2
f(30) = 8*dtau**4
f(31) = 8*dtau**6
f(32) = 8*dtau**8
v = sum(f*params)
end function c2h4_poten_n1_d8


!---------------------------------------------------------------------


function c2h4_poten_n2_d8(coords, params) result(v)
real(ark), intent(in) :: coords(12)
real(ark), intent(in) :: params(376)
real(ark) :: v
real(ark) :: r0,r1,r2,r3,r4,a1,a2,a3,a4,b1,b2,dtau
real(ark) :: f(376)
r0 = coords(1)
r1 = coords(2)
r2 = coords(3)
r3 = coords(4)
r4 = coords(5)
a1 = coords(6)
a2 = coords(7)
a3 = coords(8)
a4 = coords(9)
b1 = coords(10)
b2 = coords(11)
dtau = coords(12)
f(1) = 2*r0*(r1 + r2 + r3 + r4)
f(2) = 2*r0*(r1**2 + r2**2 + r3**2 + r4**2)
f(3) = 2*r0*(r1**3 + r2**3 + r3**3 + r4**3)
f(4) = 2*r0*(r1**4 + r2**4 + r3**4 + r4**4)
f(5) = 2*r0*(r1**5 + r2**5 + r3**5 + r4**5)
f(6) = 2*r0*(r1**6 + r2**6 + r3**6 + r4**6)
f(7) = 2*r0*(r1**7 + r2**7 + r3**7 + r4**7)
f(8) = 2*r0**2*(r1 + r2 + r3 + r4)
f(9) = 2*r0**2*(r1**2 + r2**2 + r3**2 + r4**2)
f(10) = 2*r0**2*(r1**3 + r2**3 + r3**3 + r4**3)
f(11) = 2*r0**2*(r1**4 + r2**4 + r3**4 + r4**4)
f(12) = 2*r0**2*(r1**5 + r2**5 + r3**5 + r4**5)
f(13) = 2*r0**2*(r1**6 + r2**6 + r3**6 + r4**6)
f(14) = 2*r0**3*(r1 + r2 + r3 + r4)
f(15) = 2*r0**3*(r1**2 + r2**2 + r3**2 + r4**2)
f(16) = 2*r0**3*(r1**3 + r2**3 + r3**3 + r4**3)
f(17) = 2*r0**3*(r1**4 + r2**4 + r3**4 + r4**4)
f(18) = 2*r0**3*(r1**5 + r2**5 + r3**5 + r4**5)
f(19) = 2*r0**4*(r1 + r2 + r3 + r4)
f(20) = 2*r0**4*(r1**2 + r2**2 + r3**2 + r4**2)
f(21) = 2*r0**4*(r1**3 + r2**3 + r3**3 + r4**3)
f(22) = 2*r0**4*(r1**4 + r2**4 + r3**4 + r4**4)
f(23) = 2*r0**5*(r1 + r2 + r3 + r4)
f(24) = 2*r0**5*(r1**2 + r2**2 + r3**2 + r4**2)
f(25) = 2*r0**5*(r1**3 + r2**3 + r3**3 + r4**3)
f(26) = 2*r0**6*(r1 + r2 + r3 + r4)
f(27) = 2*r0**6*(r1**2 + r2**2 + r3**2 + r4**2)
f(28) = 2*r0**7*(r1 + r2 + r3 + r4)
f(29) = 2*r0*(a1 + a2 + a3 + a4)
f(30) = 2*r0*(a1**2 + a2**2 + a3**2 + a4**2)
f(31) = 2*r0*(a1**3 + a2**3 + a3**3 + a4**3)
f(32) = 2*r0*(a1**4 + a2**4 + a3**4 + a4**4)
f(33) = 2*r0*(a1**5 + a2**5 + a3**5 + a4**5)
f(34) = 2*r0*(a1**6 + a2**6 + a3**6 + a4**6)
f(35) = 2*r0*(a1**7 + a2**7 + a3**7 + a4**7)
f(36) = 2*r0**2*(a1 + a2 + a3 + a4)
f(37) = 2*r0**2*(a1**2 + a2**2 + a3**2 + a4**2)
f(38) = 2*r0**2*(a1**3 + a2**3 + a3**3 + a4**3)
f(39) = 2*r0**2*(a1**4 + a2**4 + a3**4 + a4**4)
f(40) = 2*r0**2*(a1**5 + a2**5 + a3**5 + a4**5)
f(41) = 2*r0**2*(a1**6 + a2**6 + a3**6 + a4**6)
f(42) = 2*r0**3*(a1 + a2 + a3 + a4)
f(43) = 2*r0**3*(a1**2 + a2**2 + a3**2 + a4**2)
f(44) = 2*r0**3*(a1**3 + a2**3 + a3**3 + a4**3)
f(45) = 2*r0**3*(a1**4 + a2**4 + a3**4 + a4**4)
f(46) = 2*r0**3*(a1**5 + a2**5 + a3**5 + a4**5)
f(47) = 2*r0**4*(a1 + a2 + a3 + a4)
f(48) = 2*r0**4*(a1**2 + a2**2 + a3**2 + a4**2)
f(49) = 2*r0**4*(a1**3 + a2**3 + a3**3 + a4**3)
f(50) = 2*r0**4*(a1**4 + a2**4 + a3**4 + a4**4)
f(51) = 2*r0**5*(a1 + a2 + a3 + a4)
f(52) = 2*r0**5*(a1**2 + a2**2 + a3**2 + a4**2)
f(53) = 2*r0**5*(a1**3 + a2**3 + a3**3 + a4**3)
f(54) = 2*r0**6*(a1 + a2 + a3 + a4)
f(55) = 2*r0**6*(a1**2 + a2**2 + a3**2 + a4**2)
f(56) = 2*r0**7*(a1 + a2 + a3 + a4)
f(57) = 4*r0*(b1**2 + b2**2)
f(58) = 4*r0*(b1**4 + b2**4)
f(59) = 4*r0*(b1**6 + b2**6)
f(60) = 4*r0**2*(b1**2 + b2**2)
f(61) = 4*r0**2*(b1**4 + b2**4)
f(62) = 4*r0**2*(b1**6 + b2**6)
f(63) = 4*r0**3*(b1**2 + b2**2)
f(64) = 4*r0**3*(b1**4 + b2**4)
f(65) = 4*r0**4*(b1**2 + b2**2)
f(66) = 4*r0**4*(b1**4 + b2**4)
f(67) = 4*r0**5*(b1**2 + b2**2)
f(68) = 4*r0**6*(b1**2 + b2**2)
f(69) = 8*dtau**2*r0
f(70) = 8*dtau**4*r0
f(71) = 8*dtau**6*r0
f(72) = 8*dtau**2*r0**2
f(73) = 8*dtau**4*r0**2
f(74) = 8*dtau**6*r0**2
f(75) = 8*dtau**2*r0**3
f(76) = 8*dtau**4*r0**3
f(77) = 8*dtau**2*r0**4
f(78) = 8*dtau**4*r0**4
f(79) = 8*dtau**2*r0**5
f(80) = 8*dtau**2*r0**6
f(81) = 4*r1*r2 + 4*r3*r4
f(82) = 2*r1**2*r2 + 2*r1*r2**2 + 2*r3**2*r4 + 2*r3*r4**2
f(83) = 2*r1**3*r2 + 2*r1*r2**3 + 2*r3**3*r4 + 2*r3*r4**3
f(84) = 2*r1**4*r2 + 2*r1*r2**4 + 2*r3**4*r4 + 2*r3*r4**4
f(85) = 2*r1**5*r2 + 2*r1*r2**5 + 2*r3**5*r4 + 2*r3*r4**5
f(86) = 2*r1**6*r2 + 2*r1*r2**6 + 2*r3**6*r4 + 2*r3*r4**6
f(87) = 2*r1**7*r2 + 2*r1*r2**7 + 2*r3**7*r4 + 2*r3*r4**7
f(88) = 4*r1**2*r2**2 + 4*r3**2*r4**2
f(89) = 2*r1**3*r2**2 + 2*r1**2*r2**3 + 2*r3**3*r4**2 + 2*r3**2*r4**3
f(90) = 2*r1**4*r2**2 + 2*r1**2*r2**4 + 2*r3**4*r4**2 + 2*r3**2*r4**4
f(91) = 2*r1**5*r2**2 + 2*r1**2*r2**5 + 2*r3**5*r4**2 + 2*r3**2*r4**5
f(92) = 2*r1**6*r2**2 + 2*r1**2*r2**6 + 2*r3**6*r4**2 + 2*r3**2*r4**6
f(93) = 4*r1**3*r2**3 + 4*r3**3*r4**3
f(94) = 2*r1**4*r2**3 + 2*r1**3*r2**4 + 2*r3**4*r4**3 + 2*r3**3*r4**4
f(95) = 2*r1**5*r2**3 + 2*r1**3*r2**5 + 2*r3**5*r4**3 + 2*r3**3*r4**5
f(96) = 4*r1**4*r2**4 + 4*r3**4*r4**4
f(97) = 4*r1*r3 + 4*r2*r4
f(98) = 2*r1**2*r3 + 2*r1*r3**2 + 2*r2**2*r4 + 2*r2*r4**2
f(99) = 2*r1**3*r3 + 2*r1*r3**3 + 2*r2**3*r4 + 2*r2*r4**3
f(100) = 2*r1**4*r3 + 2*r1*r3**4 + 2*r2**4*r4 + 2*r2*r4**4
f(101) = 2*r1**5*r3 + 2*r1*r3**5 + 2*r2**5*r4 + 2*r2*r4**5
f(102) = 2*r1**6*r3 + 2*r1*r3**6 + 2*r2**6*r4 + 2*r2*r4**6
f(103) = 2*r1**7*r3 + 2*r1*r3**7 + 2*r2**7*r4 + 2*r2*r4**7
f(104) = 4*r1**2*r3**2 + 4*r2**2*r4**2
f(105) = 2*r1**3*r3**2 + 2*r1**2*r3**3 + 2*r2**3*r4**2 + 2*r2**2*r4**3
f(106) = 2*r1**4*r3**2 + 2*r1**2*r3**4 + 2*r2**4*r4**2 + 2*r2**2*r4**4
f(107) = 2*r1**5*r3**2 + 2*r1**2*r3**5 + 2*r2**5*r4**2 + 2*r2**2*r4**5
f(108) = 2*r1**6*r3**2 + 2*r1**2*r3**6 + 2*r2**6*r4**2 + 2*r2**2*r4**6
f(109) = 4*r1**3*r3**3 + 4*r2**3*r4**3
f(110) = 2*r1**4*r3**3 + 2*r1**3*r3**4 + 2*r2**4*r4**3 + 2*r2**3*r4**4
f(111) = 2*r1**5*r3**3 + 2*r1**3*r3**5 + 2*r2**5*r4**3 + 2*r2**3*r4**5
f(112) = 4*r1**4*r3**4 + 4*r2**4*r4**4
f(113) = 4*r1*r4 + 4*r2*r3
f(114) = 2*r1**2*r4 + 2*r1*r4**2 + 2*r2**2*r3 + 2*r2*r3**2
f(115) = 2*r1**3*r4 + 2*r1*r4**3 + 2*r2**3*r3 + 2*r2*r3**3
f(116) = 2*r1**4*r4 + 2*r1*r4**4 + 2*r2**4*r3 + 2*r2*r3**4
f(117) = 2*r1**5*r4 + 2*r1*r4**5 + 2*r2**5*r3 + 2*r2*r3**5
f(118) = 2*r1**6*r4 + 2*r1*r4**6 + 2*r2**6*r3 + 2*r2*r3**6
f(119) = 2*r1**7*r4 + 2*r1*r4**7 + 2*r2**7*r3 + 2*r2*r3**7
f(120) = 4*r1**2*r4**2 + 4*r2**2*r3**2
f(121) = 2*r1**3*r4**2 + 2*r1**2*r4**3 + 2*r2**3*r3**2 + 2*r2**2*r3**3
f(122) = 2*r1**4*r4**2 + 2*r1**2*r4**4 + 2*r2**4*r3**2 + 2*r2**2*r3**4
f(123) = 2*r1**5*r4**2 + 2*r1**2*r4**5 + 2*r2**5*r3**2 + 2*r2**2*r3**5
f(124) = 2*r1**6*r4**2 + 2*r1**2*r4**6 + 2*r2**6*r3**2 + 2*r2**2*r3**6
f(125) = 4*r1**3*r4**3 + 4*r2**3*r3**3
f(126) = 2*r1**4*r4**3 + 2*r1**3*r4**4 + 2*r2**4*r3**3 + 2*r2**3*r3**4
f(127) = 2*r1**5*r4**3 + 2*r1**3*r4**5 + 2*r2**5*r3**3 + 2*r2**3*r3**5
f(128) = 4*r1**4*r4**4 + 4*r2**4*r3**4
f(129) = 2*a1*r1 + 2*a2*r2 + 2*a3*r3 + 2*a4*r4
f(130) = 2*a1**2*r1 + 2*a2**2*r2 + 2*a3**2*r3 + 2*a4**2*r4
f(131) = 2*a1**3*r1 + 2*a2**3*r2 + 2*a3**3*r3 + 2*a4**3*r4
f(132) = 2*a1**4*r1 + 2*a2**4*r2 + 2*a3**4*r3 + 2*a4**4*r4
f(133) = 2*a1**5*r1 + 2*a2**5*r2 + 2*a3**5*r3 + 2*a4**5*r4
f(134) = 2*a1**6*r1 + 2*a2**6*r2 + 2*a3**6*r3 + 2*a4**6*r4
f(135) = 2*a1**7*r1 + 2*a2**7*r2 + 2*a3**7*r3 + 2*a4**7*r4
f(136) = 2*a1*r1**2 + 2*a2*r2**2 + 2*a3*r3**2 + 2*a4*r4**2
f(137) = 2*a1**2*r1**2 + 2*a2**2*r2**2 + 2*a3**2*r3**2 + 2*a4**2*r4**2
f(138) = 2*a1**3*r1**2 + 2*a2**3*r2**2 + 2*a3**3*r3**2 + 2*a4**3*r4**2
f(139) = 2*a1**4*r1**2 + 2*a2**4*r2**2 + 2*a3**4*r3**2 + 2*a4**4*r4**2
f(140) = 2*a1**5*r1**2 + 2*a2**5*r2**2 + 2*a3**5*r3**2 + 2*a4**5*r4**2
f(141) = 2*a1**6*r1**2 + 2*a2**6*r2**2 + 2*a3**6*r3**2 + 2*a4**6*r4**2
f(142) = 2*a1*r1**3 + 2*a2*r2**3 + 2*a3*r3**3 + 2*a4*r4**3
f(143) = 2*a1**2*r1**3 + 2*a2**2*r2**3 + 2*a3**2*r3**3 + 2*a4**2*r4**3
f(144) = 2*a1**3*r1**3 + 2*a2**3*r2**3 + 2*a3**3*r3**3 + 2*a4**3*r4**3
f(145) = 2*a1**4*r1**3 + 2*a2**4*r2**3 + 2*a3**4*r3**3 + 2*a4**4*r4**3
f(146) = 2*a1**5*r1**3 + 2*a2**5*r2**3 + 2*a3**5*r3**3 + 2*a4**5*r4**3
f(147) = 2*a1*r1**4 + 2*a2*r2**4 + 2*a3*r3**4 + 2*a4*r4**4
f(148) = 2*a1**2*r1**4 + 2*a2**2*r2**4 + 2*a3**2*r3**4 + 2*a4**2*r4**4
f(149) = 2*a1**3*r1**4 + 2*a2**3*r2**4 + 2*a3**3*r3**4 + 2*a4**3*r4**4
f(150) = 2*a1**4*r1**4 + 2*a2**4*r2**4 + 2*a3**4*r3**4 + 2*a4**4*r4**4
f(151) = 2*a1*r1**5 + 2*a2*r2**5 + 2*a3*r3**5 + 2*a4*r4**5
f(152) = 2*a1**2*r1**5 + 2*a2**2*r2**5 + 2*a3**2*r3**5 + 2*a4**2*r4**5
f(153) = 2*a1**3*r1**5 + 2*a2**3*r2**5 + 2*a3**3*r3**5 + 2*a4**3*r4**5
f(154) = 2*a1*r1**6 + 2*a2*r2**6 + 2*a3*r3**6 + 2*a4*r4**6
f(155) = 2*a1**2*r1**6 + 2*a2**2*r2**6 + 2*a3**2*r3**6 + 2*a4**2*r4**6
f(156) = 2*a1*r1**7 + 2*a2*r2**7 + 2*a3*r3**7 + 2*a4*r4**7
f(157) = 2*a1*r2 + 2*a2*r1 + 2*a3*r4 + 2*a4*r3
f(158) = 2*a1**2*r2 + 2*a2**2*r1 + 2*a3**2*r4 + 2*a4**2*r3
f(159) = 2*a1**3*r2 + 2*a2**3*r1 + 2*a3**3*r4 + 2*a4**3*r3
f(160) = 2*a1**4*r2 + 2*a2**4*r1 + 2*a3**4*r4 + 2*a4**4*r3
f(161) = 2*a1**5*r2 + 2*a2**5*r1 + 2*a3**5*r4 + 2*a4**5*r3
f(162) = 2*a1**6*r2 + 2*a2**6*r1 + 2*a3**6*r4 + 2*a4**6*r3
f(163) = 2*a1**7*r2 + 2*a2**7*r1 + 2*a3**7*r4 + 2*a4**7*r3
f(164) = 2*a1*r2**2 + 2*a2*r1**2 + 2*a3*r4**2 + 2*a4*r3**2
f(165) = 2*a1**2*r2**2 + 2*a2**2*r1**2 + 2*a3**2*r4**2 + 2*a4**2*r3**2
f(166) = 2*a1**3*r2**2 + 2*a2**3*r1**2 + 2*a3**3*r4**2 + 2*a4**3*r3**2
f(167) = 2*a1**4*r2**2 + 2*a2**4*r1**2 + 2*a3**4*r4**2 + 2*a4**4*r3**2
f(168) = 2*a1**5*r2**2 + 2*a2**5*r1**2 + 2*a3**5*r4**2 + 2*a4**5*r3**2
f(169) = 2*a1**6*r2**2 + 2*a2**6*r1**2 + 2*a3**6*r4**2 + 2*a4**6*r3**2
f(170) = 2*a1*r2**3 + 2*a2*r1**3 + 2*a3*r4**3 + 2*a4*r3**3
f(171) = 2*a1**2*r2**3 + 2*a2**2*r1**3 + 2*a3**2*r4**3 + 2*a4**2*r3**3
f(172) = 2*a1**3*r2**3 + 2*a2**3*r1**3 + 2*a3**3*r4**3 + 2*a4**3*r3**3
f(173) = 2*a1**4*r2**3 + 2*a2**4*r1**3 + 2*a3**4*r4**3 + 2*a4**4*r3**3
f(174) = 2*a1**5*r2**3 + 2*a2**5*r1**3 + 2*a3**5*r4**3 + 2*a4**5*r3**3
f(175) = 2*a1*r2**4 + 2*a2*r1**4 + 2*a3*r4**4 + 2*a4*r3**4
f(176) = 2*a1**2*r2**4 + 2*a2**2*r1**4 + 2*a3**2*r4**4 + 2*a4**2*r3**4
f(177) = 2*a1**3*r2**4 + 2*a2**3*r1**4 + 2*a3**3*r4**4 + 2*a4**3*r3**4
f(178) = 2*a1**4*r2**4 + 2*a2**4*r1**4 + 2*a3**4*r4**4 + 2*a4**4*r3**4
f(179) = 2*a1*r2**5 + 2*a2*r1**5 + 2*a3*r4**5 + 2*a4*r3**5
f(180) = 2*a1**2*r2**5 + 2*a2**2*r1**5 + 2*a3**2*r4**5 + 2*a4**2*r3**5
f(181) = 2*a1**3*r2**5 + 2*a2**3*r1**5 + 2*a3**3*r4**5 + 2*a4**3*r3**5
f(182) = 2*a1*r2**6 + 2*a2*r1**6 + 2*a3*r4**6 + 2*a4*r3**6
f(183) = 2*a1**2*r2**6 + 2*a2**2*r1**6 + 2*a3**2*r4**6 + 2*a4**2*r3**6
f(184) = 2*a1*r2**7 + 2*a2*r1**7 + 2*a3*r4**7 + 2*a4*r3**7
f(185) = 2*a1*r3 + 2*a2*r4 + 2*a3*r1 + 2*a4*r2
f(186) = 2*a1**2*r3 + 2*a2**2*r4 + 2*a3**2*r1 + 2*a4**2*r2
f(187) = 2*a1**3*r3 + 2*a2**3*r4 + 2*a3**3*r1 + 2*a4**3*r2
f(188) = 2*a1**4*r3 + 2*a2**4*r4 + 2*a3**4*r1 + 2*a4**4*r2
f(189) = 2*a1**5*r3 + 2*a2**5*r4 + 2*a3**5*r1 + 2*a4**5*r2
f(190) = 2*a1**6*r3 + 2*a2**6*r4 + 2*a3**6*r1 + 2*a4**6*r2
f(191) = 2*a1**7*r3 + 2*a2**7*r4 + 2*a3**7*r1 + 2*a4**7*r2
f(192) = 2*a1*r3**2 + 2*a2*r4**2 + 2*a3*r1**2 + 2*a4*r2**2
f(193) = 2*a1**2*r3**2 + 2*a2**2*r4**2 + 2*a3**2*r1**2 + 2*a4**2*r2**2
f(194) = 2*a1**3*r3**2 + 2*a2**3*r4**2 + 2*a3**3*r1**2 + 2*a4**3*r2**2
f(195) = 2*a1**4*r3**2 + 2*a2**4*r4**2 + 2*a3**4*r1**2 + 2*a4**4*r2**2
f(196) = 2*a1**5*r3**2 + 2*a2**5*r4**2 + 2*a3**5*r1**2 + 2*a4**5*r2**2
f(197) = 2*a1**6*r3**2 + 2*a2**6*r4**2 + 2*a3**6*r1**2 + 2*a4**6*r2**2
f(198) = 2*a1*r3**3 + 2*a2*r4**3 + 2*a3*r1**3 + 2*a4*r2**3
f(199) = 2*a1**2*r3**3 + 2*a2**2*r4**3 + 2*a3**2*r1**3 + 2*a4**2*r2**3
f(200) = 2*a1**3*r3**3 + 2*a2**3*r4**3 + 2*a3**3*r1**3 + 2*a4**3*r2**3
f(201) = 2*a1**4*r3**3 + 2*a2**4*r4**3 + 2*a3**4*r1**3 + 2*a4**4*r2**3
f(202) = 2*a1**5*r3**3 + 2*a2**5*r4**3 + 2*a3**5*r1**3 + 2*a4**5*r2**3
f(203) = 2*a1*r3**4 + 2*a2*r4**4 + 2*a3*r1**4 + 2*a4*r2**4
f(204) = 2*a1**2*r3**4 + 2*a2**2*r4**4 + 2*a3**2*r1**4 + 2*a4**2*r2**4
f(205) = 2*a1**3*r3**4 + 2*a2**3*r4**4 + 2*a3**3*r1**4 + 2*a4**3*r2**4
f(206) = 2*a1**4*r3**4 + 2*a2**4*r4**4 + 2*a3**4*r1**4 + 2*a4**4*r2**4
f(207) = 2*a1*r3**5 + 2*a2*r4**5 + 2*a3*r1**5 + 2*a4*r2**5
f(208) = 2*a1**2*r3**5 + 2*a2**2*r4**5 + 2*a3**2*r1**5 + 2*a4**2*r2**5
f(209) = 2*a1**3*r3**5 + 2*a2**3*r4**5 + 2*a3**3*r1**5 + 2*a4**3*r2**5
f(210) = 2*a1*r3**6 + 2*a2*r4**6 + 2*a3*r1**6 + 2*a4*r2**6
f(211) = 2*a1**2*r3**6 + 2*a2**2*r4**6 + 2*a3**2*r1**6 + 2*a4**2*r2**6
f(212) = 2*a1*r3**7 + 2*a2*r4**7 + 2*a3*r1**7 + 2*a4*r2**7
f(213) = 2*a1*r4 + 2*a2*r3 + 2*a3*r2 + 2*a4*r1
f(214) = 2*a1**2*r4 + 2*a2**2*r3 + 2*a3**2*r2 + 2*a4**2*r1
f(215) = 2*a1**3*r4 + 2*a2**3*r3 + 2*a3**3*r2 + 2*a4**3*r1
f(216) = 2*a1**4*r4 + 2*a2**4*r3 + 2*a3**4*r2 + 2*a4**4*r1
f(217) = 2*a1**5*r4 + 2*a2**5*r3 + 2*a3**5*r2 + 2*a4**5*r1
f(218) = 2*a1**6*r4 + 2*a2**6*r3 + 2*a3**6*r2 + 2*a4**6*r1
f(219) = 2*a1**7*r4 + 2*a2**7*r3 + 2*a3**7*r2 + 2*a4**7*r1
f(220) = 2*a1*r4**2 + 2*a2*r3**2 + 2*a3*r2**2 + 2*a4*r1**2
f(221) = 2*a1**2*r4**2 + 2*a2**2*r3**2 + 2*a3**2*r2**2 + 2*a4**2*r1**2
f(222) = 2*a1**3*r4**2 + 2*a2**3*r3**2 + 2*a3**3*r2**2 + 2*a4**3*r1**2
f(223) = 2*a1**4*r4**2 + 2*a2**4*r3**2 + 2*a3**4*r2**2 + 2*a4**4*r1**2
f(224) = 2*a1**5*r4**2 + 2*a2**5*r3**2 + 2*a3**5*r2**2 + 2*a4**5*r1**2
f(225) = 2*a1**6*r4**2 + 2*a2**6*r3**2 + 2*a3**6*r2**2 + 2*a4**6*r1**2
f(226) = 2*a1*r4**3 + 2*a2*r3**3 + 2*a3*r2**3 + 2*a4*r1**3
f(227) = 2*a1**2*r4**3 + 2*a2**2*r3**3 + 2*a3**2*r2**3 + 2*a4**2*r1**3
f(228) = 2*a1**3*r4**3 + 2*a2**3*r3**3 + 2*a3**3*r2**3 + 2*a4**3*r1**3
f(229) = 2*a1**4*r4**3 + 2*a2**4*r3**3 + 2*a3**4*r2**3 + 2*a4**4*r1**3
f(230) = 2*a1**5*r4**3 + 2*a2**5*r3**3 + 2*a3**5*r2**3 + 2*a4**5*r1**3
f(231) = 2*a1*r4**4 + 2*a2*r3**4 + 2*a3*r2**4 + 2*a4*r1**4
f(232) = 2*a1**2*r4**4 + 2*a2**2*r3**4 + 2*a3**2*r2**4 + 2*a4**2*r1**4
f(233) = 2*a1**3*r4**4 + 2*a2**3*r3**4 + 2*a3**3*r2**4 + 2*a4**3*r1**4
f(234) = 2*a1**4*r4**4 + 2*a2**4*r3**4 + 2*a3**4*r2**4 + 2*a4**4*r1**4
f(235) = 2*a1*r4**5 + 2*a2*r3**5 + 2*a3*r2**5 + 2*a4*r1**5
f(236) = 2*a1**2*r4**5 + 2*a2**2*r3**5 + 2*a3**2*r2**5 + 2*a4**2*r1**5
f(237) = 2*a1**3*r4**5 + 2*a2**3*r3**5 + 2*a3**3*r2**5 + 2*a4**3*r1**5
f(238) = 2*a1*r4**6 + 2*a2*r3**6 + 2*a3*r2**6 + 2*a4*r1**6
f(239) = 2*a1**2*r4**6 + 2*a2**2*r3**6 + 2*a3**2*r2**6 + 2*a4**2*r1**6
f(240) = 2*a1*r4**7 + 2*a2*r3**7 + 2*a3*r2**7 + 2*a4*r1**7
f(241) = 2*b1**2*r1 + 2*b1**2*r2 + 2*b2**2*r3 + 2*b2**2*r4
f(242) = 2*b1**4*r1 + 2*b1**4*r2 + 2*b2**4*r3 + 2*b2**4*r4
f(243) = 2*b1**6*r1 + 2*b1**6*r2 + 2*b2**6*r3 + 2*b2**6*r4
f(244) = 2*b1**2*r1**2 + 2*b1**2*r2**2 + 2*b2**2*r3**2 + 2*b2**2*r4**2
f(245) = 2*b1**4*r1**2 + 2*b1**4*r2**2 + 2*b2**4*r3**2 + 2*b2**4*r4**2
f(246) = 2*b1**6*r1**2 + 2*b1**6*r2**2 + 2*b2**6*r3**2 + 2*b2**6*r4**2
f(247) = 2*b1**2*r1**3 + 2*b1**2*r2**3 + 2*b2**2*r3**3 + 2*b2**2*r4**3
f(248) = 2*b1**4*r1**3 + 2*b1**4*r2**3 + 2*b2**4*r3**3 + 2*b2**4*r4**3
f(249) = 2*b1**2*r1**4 + 2*b1**2*r2**4 + 2*b2**2*r3**4 + 2*b2**2*r4**4
f(250) = 2*b1**4*r1**4 + 2*b1**4*r2**4 + 2*b2**4*r3**4 + 2*b2**4*r4**4
f(251) = 2*b1**2*r1**5 + 2*b1**2*r2**5 + 2*b2**2*r3**5 + 2*b2**2*r4**5
f(252) = 2*b1**2*r1**6 + 2*b1**2*r2**6 + 2*b2**2*r3**6 + 2*b2**2*r4**6
f(253) = 2*b1**2*r3 + 2*b1**2*r4 + 2*b2**2*r1 + 2*b2**2*r2
f(254) = 2*b1**4*r3 + 2*b1**4*r4 + 2*b2**4*r1 + 2*b2**4*r2
f(255) = 2*b1**6*r3 + 2*b1**6*r4 + 2*b2**6*r1 + 2*b2**6*r2
f(256) = 2*b1**2*r3**2 + 2*b1**2*r4**2 + 2*b2**2*r1**2 + 2*b2**2*r2**2
f(257) = 2*b1**4*r3**2 + 2*b1**4*r4**2 + 2*b2**4*r1**2 + 2*b2**4*r2**2
f(258) = 2*b1**6*r3**2 + 2*b1**6*r4**2 + 2*b2**6*r1**2 + 2*b2**6*r2**2
f(259) = 2*b1**2*r3**3 + 2*b1**2*r4**3 + 2*b2**2*r1**3 + 2*b2**2*r2**3
f(260) = 2*b1**4*r3**3 + 2*b1**4*r4**3 + 2*b2**4*r1**3 + 2*b2**4*r2**3
f(261) = 2*b1**2*r3**4 + 2*b1**2*r4**4 + 2*b2**2*r1**4 + 2*b2**2*r2**4
f(262) = 2*b1**4*r3**4 + 2*b1**4*r4**4 + 2*b2**4*r1**4 + 2*b2**4*r2**4
f(263) = 2*b1**2*r3**5 + 2*b1**2*r4**5 + 2*b2**2*r1**5 + 2*b2**2*r2**5
f(264) = 2*b1**2*r3**6 + 2*b1**2*r4**6 + 2*b2**2*r1**6 + 2*b2**2*r2**6
f(265) = 2*dtau**2*(r1 + r2 + r3 + r4)
f(266) = 2*dtau**4*(r1 + r2 + r3 + r4)
f(267) = 2*dtau**6*(r1 + r2 + r3 + r4)
f(268) = 2*dtau**2*(r1**2 + r2**2 + r3**2 + r4**2)
f(269) = 2*dtau**4*(r1**2 + r2**2 + r3**2 + r4**2)
f(270) = 2*dtau**6*(r1**2 + r2**2 + r3**2 + r4**2)
f(271) = 2*dtau**2*(r1**3 + r2**3 + r3**3 + r4**3)
f(272) = 2*dtau**4*(r1**3 + r2**3 + r3**3 + r4**3)
f(273) = 2*dtau**2*(r1**4 + r2**4 + r3**4 + r4**4)
f(274) = 2*dtau**4*(r1**4 + r2**4 + r3**4 + r4**4)
f(275) = 2*dtau**2*(r1**5 + r2**5 + r3**5 + r4**5)
f(276) = 2*dtau**2*(r1**6 + r2**6 + r3**6 + r4**6)
f(277) = 4*a1*a2 + 4*a3*a4
f(278) = 2*a1**2*a2 + 2*a1*a2**2 + 2*a3**2*a4 + 2*a3*a4**2
f(279) = 2*a1**3*a2 + 2*a1*a2**3 + 2*a3**3*a4 + 2*a3*a4**3
f(280) = 2*a1**4*a2 + 2*a1*a2**4 + 2*a3**4*a4 + 2*a3*a4**4
f(281) = 2*a1**5*a2 + 2*a1*a2**5 + 2*a3**5*a4 + 2*a3*a4**5
f(282) = 2*a1**6*a2 + 2*a1*a2**6 + 2*a3**6*a4 + 2*a3*a4**6
f(283) = 2*a1**7*a2 + 2*a1*a2**7 + 2*a3**7*a4 + 2*a3*a4**7
f(284) = 4*a1**2*a2**2 + 4*a3**2*a4**2
f(285) = 2*a1**3*a2**2 + 2*a1**2*a2**3 + 2*a3**3*a4**2 + 2*a3**2*a4**3
f(286) = 2*a1**4*a2**2 + 2*a1**2*a2**4 + 2*a3**4*a4**2 + 2*a3**2*a4**4
f(287) = 2*a1**5*a2**2 + 2*a1**2*a2**5 + 2*a3**5*a4**2 + 2*a3**2*a4**5
f(288) = 2*a1**6*a2**2 + 2*a1**2*a2**6 + 2*a3**6*a4**2 + 2*a3**2*a4**6
f(289) = 4*a1**3*a2**3 + 4*a3**3*a4**3
f(290) = 2*a1**4*a2**3 + 2*a1**3*a2**4 + 2*a3**4*a4**3 + 2*a3**3*a4**4
f(291) = 2*a1**5*a2**3 + 2*a1**3*a2**5 + 2*a3**5*a4**3 + 2*a3**3*a4**5
f(292) = 4*a1**4*a2**4 + 4*a3**4*a4**4
f(293) = 4*a1*a3 + 4*a2*a4
f(294) = 2*a1**2*a3 + 2*a1*a3**2 + 2*a2**2*a4 + 2*a2*a4**2
f(295) = 2*a1**3*a3 + 2*a1*a3**3 + 2*a2**3*a4 + 2*a2*a4**3
f(296) = 2*a1**4*a3 + 2*a1*a3**4 + 2*a2**4*a4 + 2*a2*a4**4
f(297) = 2*a1**5*a3 + 2*a1*a3**5 + 2*a2**5*a4 + 2*a2*a4**5
f(298) = 2*a1**6*a3 + 2*a1*a3**6 + 2*a2**6*a4 + 2*a2*a4**6
f(299) = 2*a1**7*a3 + 2*a1*a3**7 + 2*a2**7*a4 + 2*a2*a4**7
f(300) = 4*a1**2*a3**2 + 4*a2**2*a4**2
f(301) = 2*a1**3*a3**2 + 2*a1**2*a3**3 + 2*a2**3*a4**2 + 2*a2**2*a4**3
f(302) = 2*a1**4*a3**2 + 2*a1**2*a3**4 + 2*a2**4*a4**2 + 2*a2**2*a4**4
f(303) = 2*a1**5*a3**2 + 2*a1**2*a3**5 + 2*a2**5*a4**2 + 2*a2**2*a4**5
f(304) = 2*a1**6*a3**2 + 2*a1**2*a3**6 + 2*a2**6*a4**2 + 2*a2**2*a4**6
f(305) = 4*a1**3*a3**3 + 4*a2**3*a4**3
f(306) = 2*a1**4*a3**3 + 2*a1**3*a3**4 + 2*a2**4*a4**3 + 2*a2**3*a4**4
f(307) = 2*a1**5*a3**3 + 2*a1**3*a3**5 + 2*a2**5*a4**3 + 2*a2**3*a4**5
f(308) = 4*a1**4*a3**4 + 4*a2**4*a4**4
f(309) = 4*a1*a4 + 4*a2*a3
f(310) = 2*a1**2*a4 + 2*a1*a4**2 + 2*a2**2*a3 + 2*a2*a3**2
f(311) = 2*a1**3*a4 + 2*a1*a4**3 + 2*a2**3*a3 + 2*a2*a3**3
f(312) = 2*a1**4*a4 + 2*a1*a4**4 + 2*a2**4*a3 + 2*a2*a3**4
f(313) = 2*a1**5*a4 + 2*a1*a4**5 + 2*a2**5*a3 + 2*a2*a3**5
f(314) = 2*a1**6*a4 + 2*a1*a4**6 + 2*a2**6*a3 + 2*a2*a3**6
f(315) = 2*a1**7*a4 + 2*a1*a4**7 + 2*a2**7*a3 + 2*a2*a3**7
f(316) = 4*a1**2*a4**2 + 4*a2**2*a3**2
f(317) = 2*a1**3*a4**2 + 2*a1**2*a4**3 + 2*a2**3*a3**2 + 2*a2**2*a3**3
f(318) = 2*a1**4*a4**2 + 2*a1**2*a4**4 + 2*a2**4*a3**2 + 2*a2**2*a3**4
f(319) = 2*a1**5*a4**2 + 2*a1**2*a4**5 + 2*a2**5*a3**2 + 2*a2**2*a3**5
f(320) = 2*a1**6*a4**2 + 2*a1**2*a4**6 + 2*a2**6*a3**2 + 2*a2**2*a3**6
f(321) = 4*a1**3*a4**3 + 4*a2**3*a3**3
f(322) = 2*a1**4*a4**3 + 2*a1**3*a4**4 + 2*a2**4*a3**3 + 2*a2**3*a3**4
f(323) = 2*a1**5*a4**3 + 2*a1**3*a4**5 + 2*a2**5*a3**3 + 2*a2**3*a3**5
f(324) = 4*a1**4*a4**4 + 4*a2**4*a3**4
f(325) = 2*a1*b1**2 + 2*a2*b1**2 + 2*a3*b2**2 + 2*a4*b2**2
f(326) = 2*a1*b1**4 + 2*a2*b1**4 + 2*a3*b2**4 + 2*a4*b2**4
f(327) = 2*a1*b1**6 + 2*a2*b1**6 + 2*a3*b2**6 + 2*a4*b2**6
f(328) = 2*a1**2*b1**2 + 2*a2**2*b1**2 + 2*a3**2*b2**2 + 2*a4**2*b2**2
f(329) = 2*a1**2*b1**4 + 2*a2**2*b1**4 + 2*a3**2*b2**4 + 2*a4**2*b2**4
f(330) = 2*a1**2*b1**6 + 2*a2**2*b1**6 + 2*a3**2*b2**6 + 2*a4**2*b2**6
f(331) = 2*a1**3*b1**2 + 2*a2**3*b1**2 + 2*a3**3*b2**2 + 2*a4**3*b2**2
f(332) = 2*a1**3*b1**4 + 2*a2**3*b1**4 + 2*a3**3*b2**4 + 2*a4**3*b2**4
f(333) = 2*a1**4*b1**2 + 2*a2**4*b1**2 + 2*a3**4*b2**2 + 2*a4**4*b2**2
f(334) = 2*a1**4*b1**4 + 2*a2**4*b1**4 + 2*a3**4*b2**4 + 2*a4**4*b2**4
f(335) = 2*a1**5*b1**2 + 2*a2**5*b1**2 + 2*a3**5*b2**2 + 2*a4**5*b2**2
f(336) = 2*a1**6*b1**2 + 2*a2**6*b1**2 + 2*a3**6*b2**2 + 2*a4**6*b2**2
f(337) = 2*a1*b2**2 + 2*a2*b2**2 + 2*a3*b1**2 + 2*a4*b1**2
f(338) = 2*a1*b2**4 + 2*a2*b2**4 + 2*a3*b1**4 + 2*a4*b1**4
f(339) = 2*a1*b2**6 + 2*a2*b2**6 + 2*a3*b1**6 + 2*a4*b1**6
f(340) = 2*a1**2*b2**2 + 2*a2**2*b2**2 + 2*a3**2*b1**2 + 2*a4**2*b1**2
f(341) = 2*a1**2*b2**4 + 2*a2**2*b2**4 + 2*a3**2*b1**4 + 2*a4**2*b1**4
f(342) = 2*a1**2*b2**6 + 2*a2**2*b2**6 + 2*a3**2*b1**6 + 2*a4**2*b1**6
f(343) = 2*a1**3*b2**2 + 2*a2**3*b2**2 + 2*a3**3*b1**2 + 2*a4**3*b1**2
f(344) = 2*a1**3*b2**4 + 2*a2**3*b2**4 + 2*a3**3*b1**4 + 2*a4**3*b1**4
f(345) = 2*a1**4*b2**2 + 2*a2**4*b2**2 + 2*a3**4*b1**2 + 2*a4**4*b1**2
f(346) = 2*a1**4*b2**4 + 2*a2**4*b2**4 + 2*a3**4*b1**4 + 2*a4**4*b1**4
f(347) = 2*a1**5*b2**2 + 2*a2**5*b2**2 + 2*a3**5*b1**2 + 2*a4**5*b1**2
f(348) = 2*a1**6*b2**2 + 2*a2**6*b2**2 + 2*a3**6*b1**2 + 2*a4**6*b1**2
f(349) = 2*dtau**2*(a1 + a2 + a3 + a4)
f(350) = 2*dtau**4*(a1 + a2 + a3 + a4)
f(351) = 2*dtau**6*(a1 + a2 + a3 + a4)
f(352) = 2*dtau**2*(a1**2 + a2**2 + a3**2 + a4**2)
f(353) = 2*dtau**4*(a1**2 + a2**2 + a3**2 + a4**2)
f(354) = 2*dtau**6*(a1**2 + a2**2 + a3**2 + a4**2)
f(355) = 2*dtau**2*(a1**3 + a2**3 + a3**3 + a4**3)
f(356) = 2*dtau**4*(a1**3 + a2**3 + a3**3 + a4**3)
f(357) = 2*dtau**2*(a1**4 + a2**4 + a3**4 + a4**4)
f(358) = 2*dtau**4*(a1**4 + a2**4 + a3**4 + a4**4)
f(359) = 2*dtau**2*(a1**5 + a2**5 + a3**5 + a4**5)
f(360) = 2*dtau**2*(a1**6 + a2**6 + a3**6 + a4**6)
f(361) = 8*b1*b2
f(362) = 4*b1*b2*(b1**2 + b2**2)
f(363) = 4*b1*b2*(b1**4 + b2**4)
f(364) = 4*b1*b2*(b1**6 + b2**6)
f(365) = 8*b1**2*b2**2
f(366) = 4*b1**2*b2**2*(b1**2 + b2**2)
f(367) = 4*b1**2*b2**2*(b1**4 + b2**4)
f(368) = 8*b1**3*b2**3
f(369) = 4*b1**3*b2**3*(b1**2 + b2**2)
f(370) = 8*b1**4*b2**4
f(371) = 4*dtau**2*(b1**2 + b2**2)
f(372) = 4*dtau**4*(b1**2 + b2**2)
f(373) = 4*dtau**6*(b1**2 + b2**2)
f(374) = 4*dtau**2*(b1**4 + b2**4)
f(375) = 4*dtau**4*(b1**4 + b2**4)
f(376) = 4*dtau**2*(b1**6 + b2**6)
v = sum(f*params)
end function c2h4_poten_n2_d8



function MLpoten_c2h4_lee(ncoords, natoms, local, xyz, force) result(f)

  integer(ik),intent(in) :: ncoords, natoms
  real(ark),intent(in)   :: local(ncoords)
  real(ark),intent(in)   :: xyz(natoms,3)
  real(ark),intent(in)   :: force(:)
  real(ark)              :: f,W(12,12)

  real(ark) :: y(12)
  real(ark) :: tau4213,tau5124,tau5126,tau6213
  integer(ik) :: k1,k2


  tau4213 =-local(10)
  tau4213 = mod(tau4213+2.0_ark*pi,2.0_ark*pi)
  if (tau4213>pi) tau4213 = tau4213 - 2.0_ark*pi 
  !
  tau5124 = local(11)
  tau5124 = mod(tau5124+2.0_ark*pi,2.0_ark*pi)
  !
  tau5126 =-local(12)
  tau5126 = mod(tau5126+2.0_ark*pi,2.0_ark*pi)
  if (tau5126>pi) tau5126 = tau5126 - 2.0_ark*pi 
  !
  tau6213 = 2.0_ark*pi-tau5124-tau4213-tau5126
  tau6213 = mod(tau6213+2.0_ark*pi,2.0_ark*pi)

  y(1) = 0.5_ark*sum(local(2:5))-molec%req(2)*2.0_ark
  y(2) = local(1)-molec%req(1)
  y(3) = 0.5_ark*sum(local(6:9))-molec%alphaeq(1)*2.0_ark
  y(4) = (tau5126+tau4213-molec%taueq(1)*2.0_ark)/sqrt(2.0_ark)
  y(5) = 0.5_ark*( local(2)-local(3)-local(4)+local(5) )
  y(6) = 0.5_ark*( local(6)-local(7)-local(8)+local(9) )
  y(7) = (tau6213-tau5124)/sqrt(2.0_ark)
  y(8) = (tau5126-tau4213)/sqrt(2.0_ark)
  y(9) = 0.5_ark*(-local(2)-local(3)+local(4)+local(5) )
  y(10)= 0.5_ark*(-local(6)-local(7)+local(8)+local(9) )
  y(11)= 0.5_ark*(-local(2)+local(3)-local(4)+local(5) )
  y(12)= 0.5_ark*(-local(6)+local(7)-local(8)+local(9) )
  !
  W = 0
  !
  W( 1, 1)= force( 1)
  W( 2, 1)= force( 2)
  W( 3, 1)= force( 3)
  W( 2, 2)= force( 4)
  W( 3, 2)= force( 5)
  W( 3, 3)= force( 6)
  W( 4, 4)= force( 7)
  W( 5, 5)= force( 8)
  W( 6, 5)= force( 9)
  W( 6, 6)= force(10)
  W( 7, 7)= force(11)
  W( 8, 8)= force(12)
  W( 9, 9)= force(13)
  W(10, 9)= force(14)
  W(10,10)= force(15)
  W(11,11)= force(16)
  W(12,11)= force(17)
  W(12,12)= force(18)
  !
  f = 0
  !
  do k1 = 1,12
      f = f + w(k1,k1)*y(k1)**2*0.5_ark
  enddo
  !
  do k1 = 1,12
    do k2 = k1+1,12
      !
      W(k1,k2) = W(k2,k1)
      !
      f = f + w(k2,k1)*y(k1)*y(k2)
      !
    enddo
  enddo
  !
  f = f*1.0e-11/planck/vellgt 
  !
  !f = dot_product(y(:),matmul(W,y(:)))*1.0e-11/planck/vellgt 
  !
end function MLpoten_c2h4_lee

end module pot_c2h4