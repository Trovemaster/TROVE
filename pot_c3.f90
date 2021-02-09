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
  contains
  !
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
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   f = MLpoten_c3_tyuterev_cotan(ncoords,natoms,local,xyz,force)
   !
 end function MLpoten
 !

 !
 function MLpoten_c3_tyuterev_cotan(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   real(ark)            :: r12,r32,alpha,xcos,v,v1,v2,v3,v4,v5,v6,v7,v8

   real(ark)            :: aa1,re12,alphae,xst,y1,y2,y3,xs1,xs2,v0,vp1,vp2,vp3
   real(ark)            :: g1,g2,b1,b2,rhh,vhh
   integer(ik)          :: N
   real(ark)             :: ycos,v_t,q1,q2
   !
   if (verbose>=6) write(out,"('MLpoten_xy2_tyuterev/start')")
   !
   r12  = local(1) ;  r32    = local(2) ;  alpha = local(3)
   !
   re12    = molec%req(1)
   aa1     = molec%specparam(1)
   !
   if (molec%Nangles>0) then
     alphae  = molec%alphaeq(1)
   else
     alphae = pi+(-molec%taueq(1)+molec%taueq(2))
     alpha = pi - local(3) 
   endif 
   !
   if  (molec%Ndihedrals>1) then
     !
     !alpha = pi-  asin( sqrt( sin(local(3))**2+sin(local(4))**2 ))
     !
     alpha = pi - asin( sqrt( (local(3))**2+(local(4))**2 ))
     !
   endif 
   !
   b1   = force(1)
   b2   = force(2)
   g1   = force(3)
   g2   = force(4)
   !
   rhh=sqrt(r12**2+r32**2-2.d0*r12*r32*cos(alpha))
   vhh=b1*exp(-g1*rhh)+b2*exp(-g2*rhh**2)
   !
   ! calculate potential energy function values
   !
   y1=1.0_ark-exp(-aa1*(r12-re12))
   y2=1.0_ark-exp(-aa1*(r32-re12))
   !
   y3=cos(alpha)-cos(alphae)
   y3=cotan(alpha*0.5d0)
   !
   N = size(force)
   !
   v4 = 0 ; v5 = 0 ; v6 = 0 ; v7 = 0 ; v8 = 0
   !
 v0 = force(5)*y1**0*y2**0*y3**0
 v1 = force(6)*y1**0*y2**0*y3**1&
    + force(7)*y1**1*y2**0*y3**0&
    + force(7)*y1**0*y2**1*y3**0
 v2 = force(8)*y1**0*y2**0*y3**2&
    + force(9)*y1**1*y2**0*y3**1&
    + force(9)*y1**0*y2**1*y3**1&
    + force(10)*y1**1*y2**1*y3**0&
    + force(11)*y1**2*y2**0*y3**0&
    + force(11)*y1**0*y2**2*y3**0
    !
    y3=(cos(alpha)-cos(alphae)) 
    !
 v3 = force(12)*y1**0*y2**0*y3**3&
    + force(13)*y1**1*y2**0*y3**2&
    + force(13)*y1**0*y2**1*y3**2&
    + force(14)*y1**1*y2**1*y3**1&
    + force(15)*y1**2*y2**0*y3**1&
    + force(15)*y1**0*y2**2*y3**1&
    + force(16)*y1**2*y2**1*y3**0&
    + force(16)*y1**1*y2**2*y3**0&
    + force(17)*y1**3*y2**0*y3**0&
    + force(17)*y1**0*y2**3*y3**0

 if (N>18) then
  v4 = force(18)*y1**0*y2**0*y3**4&
    + force(19)*y1**1*y2**0*y3**3&
    + force(19)*y1**0*y2**1*y3**3&
    + force(20)*y1**1*y2**1*y3**2&
    + force(21)*y1**2*y2**0*y3**2&
    + force(21)*y1**0*y2**2*y3**2&
    + force(22)*y1**2*y2**1*y3**1&
    + force(22)*y1**1*y2**2*y3**1&
    + force(23)*y1**2*y2**2*y3**0&
    + force(24)*y1**3*y2**0*y3**1&
    + force(24)*y1**0*y2**3*y3**1&
    + force(25)*y1**3*y2**1*y3**0&
    + force(25)*y1**1*y2**3*y3**0&
    + force(26)*y1**4*y2**0*y3**0&
    + force(26)*y1**0*y2**4*y3**0
endif

 if (N>26) then
  v5 = force(27)*y1**0*y2**0*y3**5&
    + force(28)*y1**1*y2**0*y3**4&
    + force(28)*y1**0*y2**1*y3**4&
    + force(29)*y1**1*y2**1*y3**3&
    + force(30)*y1**2*y2**0*y3**3&
    + force(30)*y1**0*y2**2*y3**3&
    + force(31)*y1**2*y2**1*y3**2&
    + force(31)*y1**1*y2**2*y3**2&
    + force(32)*y1**2*y2**2*y3**1&
    + force(33)*y1**3*y2**0*y3**2&
    + force(33)*y1**0*y2**3*y3**2&
    + force(34)*y1**3*y2**1*y3**1&
    + force(34)*y1**1*y2**3*y3**1&
    + force(35)*y1**3*y2**2*y3**0&
    + force(35)*y1**2*y2**3*y3**0&
    + force(36)*y1**4*y2**0*y3**1&
    + force(36)*y1**0*y2**4*y3**1&
    + force(37)*y1**4*y2**1*y3**0&
    + force(37)*y1**1*y2**4*y3**0&
    + force(38)*y1**5*y2**0*y3**0&
    + force(38)*y1**0*y2**5*y3**0
endif

 if (N>38) then
  v6 = force(39)*y1**0*y2**0*y3**6&
    + force(40)*y1**1*y2**0*y3**5&
    + force(40)*y1**0*y2**1*y3**5&
    + force(41)*y1**1*y2**1*y3**4&
    + force(42)*y1**2*y2**0*y3**4&
    + force(42)*y1**0*y2**2*y3**4&
    + force(43)*y1**2*y2**1*y3**3&
    + force(43)*y1**1*y2**2*y3**3&
    + force(44)*y1**2*y2**2*y3**2&
    + force(45)*y1**3*y2**0*y3**3&
    + force(45)*y1**0*y2**3*y3**3&
    + force(46)*y1**3*y2**1*y3**2&
    + force(46)*y1**1*y2**3*y3**2&
    + force(47)*y1**3*y2**2*y3**1&
    + force(47)*y1**2*y2**3*y3**1&
    + force(48)*y1**3*y2**3*y3**0&
    + force(49)*y1**4*y2**0*y3**2&
    + force(49)*y1**0*y2**4*y3**2&
    + force(50)*y1**4*y2**1*y3**1&
    + force(50)*y1**1*y2**4*y3**1&
    + force(51)*y1**4*y2**2*y3**0&
    + force(51)*y1**2*y2**4*y3**0&
    + force(52)*y1**5*y2**0*y3**1&
    + force(52)*y1**0*y2**5*y3**1&
    + force(53)*y1**5*y2**1*y3**0&
    + force(53)*y1**1*y2**5*y3**0&
    + force(54)*y1**6*y2**0*y3**0&
    + force(54)*y1**0*y2**6*y3**0
 endif

 if (N>54) then
 v7 = force(55)*y1**0*y2**0*y3**7&
    + force(56)*y1**1*y2**0*y3**6&
    + force(56)*y1**0*y2**1*y3**6&
    + force(57)*y1**1*y2**1*y3**5&
    + force(58)*y1**2*y2**0*y3**5&
    + force(58)*y1**0*y2**2*y3**5&
    + force(59)*y1**2*y2**1*y3**4&
    + force(59)*y1**1*y2**2*y3**4&
    + force(60)*y1**2*y2**2*y3**3&
    + force(61)*y1**3*y2**0*y3**4&
    + force(61)*y1**0*y2**3*y3**4&
    + force(62)*y1**3*y2**1*y3**3&
    + force(62)*y1**1*y2**3*y3**3&
    + force(63)*y1**3*y2**2*y3**2&
    + force(63)*y1**2*y2**3*y3**2&
    + force(64)*y1**3*y2**3*y3**1&
    + force(65)*y1**4*y2**0*y3**3&
    + force(65)*y1**0*y2**4*y3**3&
    + force(66)*y1**4*y2**1*y3**2&
    + force(66)*y1**1*y2**4*y3**2&
    + force(67)*y1**4*y2**2*y3**1&
    + force(67)*y1**2*y2**4*y3**1&
    + force(68)*y1**4*y2**3*y3**0&
    + force(68)*y1**3*y2**4*y3**0&
    + force(69)*y1**5*y2**0*y3**2&
    + force(69)*y1**0*y2**5*y3**2&
    + force(70)*y1**5*y2**1*y3**1&
    + force(70)*y1**1*y2**5*y3**1&
    + force(71)*y1**5*y2**2*y3**0&
    + force(71)*y1**2*y2**5*y3**0&
    + force(72)*y1**6*y2**0*y3**1&
    + force(72)*y1**0*y2**6*y3**1&
    + force(73)*y1**6*y2**1*y3**0&
    + force(73)*y1**1*y2**6*y3**0&
    + force(74)*y1**7*y2**0*y3**0&
    + force(74)*y1**0*y2**7*y3**0
 endif

 if (N>74) then
 v8 = force(75)*y1**0*y2**0*y3**8&
    + force(76)*y1**1*y2**0*y3**7&
    + force(76)*y1**0*y2**1*y3**7&
    + force(77)*y1**1*y2**1*y3**6&
    + force(78)*y1**2*y2**0*y3**6&
    + force(78)*y1**0*y2**2*y3**6&
    + force(79)*y1**2*y2**1*y3**5&
    + force(79)*y1**1*y2**2*y3**5&
    + force(80)*y1**2*y2**2*y3**4&
    + force(81)*y1**3*y2**0*y3**5&
    + force(81)*y1**0*y2**3*y3**5&
    + force(82)*y1**3*y2**1*y3**4&
    + force(82)*y1**1*y2**3*y3**4&
    + force(83)*y1**3*y2**2*y3**3&
    + force(83)*y1**2*y2**3*y3**3&
    + force(84)*y1**3*y2**3*y3**2&
    + force(85)*y1**4*y2**0*y3**4&
    + force(85)*y1**0*y2**4*y3**4&
    + force(86)*y1**4*y2**1*y3**3&
    + force(86)*y1**1*y2**4*y3**3&
    + force(87)*y1**4*y2**2*y3**2&
    + force(87)*y1**2*y2**4*y3**2&
    + force(88)*y1**4*y2**3*y3**1&
    + force(88)*y1**3*y2**4*y3**1&
    + force(89)*y1**4*y2**4*y3**0&
    + force(90)*y1**5*y2**0*y3**3&
    + force(90)*y1**0*y2**5*y3**3&
    + force(91)*y1**5*y2**1*y3**2&
    + force(91)*y1**1*y2**5*y3**2&
    + force(92)*y1**5*y2**2*y3**1&
    + force(92)*y1**2*y2**5*y3**1&
    + force(93)*y1**5*y2**3*y3**0&
    + force(93)*y1**3*y2**5*y3**0&
    + force(94)*y1**6*y2**0*y3**2&
    + force(94)*y1**0*y2**6*y3**2&
    + force(95)*y1**6*y2**1*y3**1&
    + force(95)*y1**1*y2**6*y3**1&
    + force(96)*y1**6*y2**2*y3**0&
    + force(96)*y1**2*y2**6*y3**0&
    + force(97)*y1**7*y2**0*y3**1&
    + force(97)*y1**0*y2**7*y3**1&
    + force(98)*y1**7*y2**1*y3**0&
    + force(98)*y1**1*y2**7*y3**0&
    + force(99)*y1**8*y2**0*y3**0&
    + force(99)*y1**0*y2**8*y3**0
endif
    !
    !th1 = 0.5d0*( 1.0d0-tanh( 0.0001_ark*( v0+v1+v2-50000_ark ) ) )
    !
    !f=(v0+v1+v2)+(v3+v4+v5+v6+v7+v8)*th1+vhh
    !
    !vhh = (v0+v1+v2)*(1.0_ark/r12**12+1.0_ark/r32**12+1.0_ark/rhh**12)*b1
    !
    f=v0+v1+v2+v3+v4+v5+v6+v7+v8+vhh
    !
    q1 = local(1)/bohr ; q2 = local(2)/bohr ;  ycos = cos(alpha)
    !
    v_t = 0
    !
    !call potv(v_t,q1,q2,ycos)
    !
    f = f + v_t*219474.630670_ark
    !
    if (verbose>=6) write(out,"('MLpoten_xy2_tyuterev/end')")
      !
  end function MLpoten_c3_tyuterev_cotan

end module pot_user
