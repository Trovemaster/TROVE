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
 !
 !
 recursive subroutine MLdipole(rank,ncoords,natoms,local,xyz,f)
   !
   integer(ik),intent(in) ::  rank,ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
   real(ark),intent(out)  ::  f(rank)
   !
   f = MLdipole_xy4_dF(xyz)
   !
 end subroutine MLdipole
 !
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
   f = MLpoten_xy4_qz_f12(ncoords,natoms,local,xyz,force)
   !
 end function MLpoten
 !
 !
 ! CH4 PES ccsd(t)-f12/cc-vpQZ 
 ! Yurchenko, June 2011
 ! 
 function MLpoten_xy4_qz_f12(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
    !
    real(ark) :: dx(5,3)
    integer(ik),parameter :: n = 289
    integer(ik) :: i,k,nmax
    real(ark) :: dF(n)
     !
     call potch4_diff_V(n,local,xyz,dF)
     !
     f = 0
     !
     nmax = min(size(force),molec%parmax)
     !
     do i = 3,nmax
       !
       k = molec%pot_ind(1,i)
       !
       f = f + force(i)*dF(k)
       !
     enddo
     !
 end function MLpoten_xy4_qz_f12



  subroutine potch4_diff_V(n,local,xyz,dF)
    !
    integer(ik),intent(in)  :: n
    real(ark),intent(in)  :: local(:),xyz(5,3)
    real(ark),intent(out) :: dF(n)
    !
    real(ark) :: r1,r2,r3,r4,cosa12,alpha12,cosa13,alpha13,cosa23,alpha23,cosa14,alpha14,cosa24,alpha24,cosa34,alpha34,beta312,beta412,cosbeta
    !
    real(ark) :: s1,s2,s3,s4,dx(4,3)
    real(ark) ::   betac, betaa,rc2,r02,a12,a13,a14,a23,a24,a34,t12,t13,t14,t23,t24,t34,vdump,cosae
    !
    real(ark) :: v
    real(ark) :: re,alphae,a,y1,y2,y3,y4,y5,y6,y7,y8,y9
    character(len=cl)         :: txt
    !
    real(ark) :: ae= 1.9106332362490185563277142050315_ark

      if (size(local)==9.and.molec%NDihedrals==0) then
        !
        r1  = local(1) 
        r2  = local(2) 
        r3  = local(3) 
        r4  = local(4) 
        !
        alpha12 = local(5)
        alpha13 = local(6)
        alpha23 = local(7)
        alpha14 = local(8)
        alpha24 = local(9)
        !
        cosbeta = (cos(alpha23)-cos(alpha12)*cos(alpha13) )/(sin(alpha12)*sin(alpha13))
        beta312 = aacos(cosbeta,txt)
        !
        cosbeta = (cos(alpha24)-cos(alpha12)*cos(alpha14) )/(sin(alpha12)*sin(alpha14))
        beta412 = aacos(cosbeta,txt)
        !
        cosa34 = cos(alpha13)*cos(alpha14)+cos(beta312+beta412)*sin(alpha13)*sin(alpha14)
        alpha34 = aacos(cosa34,txt)
        ! 
      else
        !
        !  set up for ch4
        !  c is atom 1 and then come the four h atoms
        !  get the four ch distances
        !
        dx(1,:)=xyz(2,:)-xyz(1,:)
        dx(2,:)=xyz(3,:)-xyz(1,:)
        dx(3,:)=xyz(4,:)-xyz(1,:)
        dx(4,:)=xyz(5,:)-xyz(1,:)
        !
        r1=sqrt(sum(dx(1,:)**2))
        r2=sqrt(sum(dx(2,:)**2))
        r3=sqrt(sum(dx(3,:)**2))
        r4=sqrt(sum(dx(4,:)**2))
        !
        !  now get the angles
        !
        cosa12=(sum(dx(1,:)*dx(2,:)))/(r1*r2)
        cosa13=(sum(dx(1,:)*dx(3,:)))/(r1*r3)
        cosa14=(sum(dx(1,:)*dx(4,:)))/(r1*r4)
        cosa23=(sum(dx(2,:)*dx(3,:)))/(r2*r3)
        cosa24=(sum(dx(2,:)*dx(4,:)))/(r2*r4)
        cosa34=(sum(dx(3,:)*dx(4,:)))/(r3*r4)
        !
        alpha12=acos(cosa12)
        alpha13=acos(cosa13)
        alpha14=acos(cosa14)
        alpha23=acos(cosa23)
        alpha24=acos(cosa24)
        alpha34=acos(cosa34)
        !
      endif
      !
      re      = molec%force(1)
      alphae  = ae
      !
      cosae = -1.0_ark/3.0_ark
      !
      !ae = acos(cosae)
      !
      a       = molec%force(2)
      !
      !betaa   = molec%force(3)
      !betac   = molec%force(4)
      !
      y1=1.0_ark-exp(-a*(r1-re))
      y2=1.0_ark-exp(-a*(r2-re))
      y3=1.0_ark-exp(-a*(r3-re))
      y4=1.0_ark-exp(-a*(r4-re))
      !
      y5=(2.0_ark*alpha12-alpha13-alpha14-alpha23-alpha24+2.0_ark*alpha34)/sqrt(12.0_ark)
      y6=(alpha13-alpha14-alpha23+alpha24)*0.5_ark
      y7=(alpha24-alpha13)/sqrt(2.0_ark)
      y8=(alpha23-alpha14)/sqrt(2.0_ark)
      y9=(alpha34-alpha12)/sqrt(2.0_ark)
      !
      rc2 = (r1-re)**2+(r2-re)**2+(r3-re)**2+(r4-re)**2
      !vdump = exp(-betac*rc2)
      !
      vdump = 1.0_ark
      
      !t12 = exp(-betaa*((r1-re)**2+(r2-re)**2))
      !t13 = exp(-betaa*((r1-re)**2+(r3-re)**2))
      !t14 = exp(-betaa*((r1-re)**2+(r4-re)**2))
      !t23 = exp(-betaa*((r2-re)**2+(r3-re)**2))
      !t24 = exp(-betaa*((r2-re)**2+(r4-re)**2))
      !t34 = exp(-betaa*((r3-re)**2+(r4-re)**2))
      
      a12 = cos(alpha12)-cosae
      a13 = cos(alpha13)-cosae
      a14 = cos(alpha14)-cosae
      a23 = cos(alpha23)-cosae
      a24 = cos(alpha24)-cosae
      a34 = cos(alpha34)-cosae
      !
      dF(1) = 0._ark
      dF(2) = 0._ark
      dF(3) = 1.0_ark
      dF(4) = y2+y3+y4+y1
      dF(5) = y8**2+y7**2+y9**2
      dF(6) = y6**2+y5**2
      dF(7) = (-y7-y8-y9)*y1+(y7-y9+y8)*y2+(y8+y9-y7)*y3+(y9+y7-y8)*y4
      dF(8) = (y4+y3+y2)*y1+(y4+y3)*y2+y3*y4
      dF(9) = y2**2+y3**2+y4**2+y1**2
      dF(10) = y7*y8*y9
      dF(11) = (-sqrt(3._ark)*y7**2/3._ark-sqrt(3._ark)*y8**2/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**2)*y5+(y7**2-y8**2)*y6
      dF(12) = y5**3-3._ark*y5*y6**2
      dF(13) = ((y8-2._ark*y9+y7)*y5+(sqrt(3._ark)*y8-sqrt(3._ark)*y7)*y6)*y1+((-y8-2._ark*y9- &
          y7)*y5+(sqrt(3._ark)*y7-sqrt(3._ark)*y8)*y6)*y2+((2._ark*y9+y7-y8)*y5+(-sqrt(3._ark)*y8- &
          sqrt(3._ark)*y7)*y6)*y3+((2._ark*y9+y8-y7)*y5+(sqrt(3._ark)*y7+sqrt(3._ark)*y8)*y6)*y4
      dF(14) = ((y9+y8)*y7+y8*y9)*y1+((y8-y9)*y7-y8*y9)*y2+((-y9-y8)*y7+y8*y9)*y3+((- &
          y8+y9)*y7-y8*y9)*y4
      dF(15) = (y8**2+y7**2+y9**2)*y1+(y8**2+y7**2+y9**2)*y2+(y8**2+y7**2+y9**2)*y3+ &
          (y8**2+y7**2+y9**2)*y4
      dF(16) = (y6**2+y5**2)*y1+(y6**2+y5**2)*y2+(y6**2+y5**2)*y3+(y6**2+y5**2)*y4
      dF(17) = (y3*y7+y4*y8+y2*y9)*y1+(-y3*y8-y4*y7)*y2-y3*y4*y9
      dF(18) = (y2*y5+(-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y3+(-y5/2._ark-sqrt(3._ark)*y6/ &
          2._ark)*y4)*y1+((-y5/2._ark-sqrt(3._ark)*y6/2._ark)*y3+(-y5/2._ark+sqrt(3._ark)*y6/ &
          2._ark)*y4)*y2+y3*y4*y5
      dF(19) = ((y4+y3)*y2+y3*y4)*y1+y2*y3*y4
      dF(20) = (y9+y7+y8)*y1**2+(-y7+y9-y8)*y2**2+(-y8-y9+y7)*y3**2+(-y7-y9+y8)*y4**2
      dF(21) = (y4+y3+y2)*y1**2+(y3**2+y2**2+y4**2)*y1+(y4+y3)*y2**2+(y3**2+y4**2)*y2+ &
          y3**2*y4+y3*y4**2
      dF(22) = y3**3+y1**3+y4**3+y2**3
      dF(23) = (y9**2+y8**2)*y7**2+y8**2*y9**2
      dF(24) = y9**4+y8**4+y7**4
      dF(25) = -sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/ &
          6._ark-sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2
      dF(26) = (y8**2+y7**2+y9**2)*y5**2+(y8**2+y7**2+y9**2)*y6**2
      dF(27) = y6**4+y5**4+2._ark*y5**2*y6**2
      dF(28) = (y3**2+y2**2+y4**2)*y1**2+(y3**2+y4**2)*y2**2+y3**2*y4**2
      dF(29) = (-y7-y8-y9)*y1**3+(y7-y9+y8)*y2**3+(y8+y9-y7)*y3**3+(y9+y7-y8)*y4**3
      dF(30) = y4**4+y3**4+y2**4+y1**4
      dF(31) = y3*y7*y8*y9+y2*y7*y8*y9+y4*y7*y8*y9+y1*y7*y8*y9
      dF(32) = ((y9+y8)*y7**2+(y9**2+y8**2)*y7+y8*y9**2+y8**2*y9)*y1+((-y8+y9)*y7**2+ &
          (-y9**2-y8**2)*y7+y8**2*y9-y8*y9**2)*y2+((-y9-y8)*y7**2+(y9**2+y8**2)*y7- &
          y8**2*y9-y8*y9**2)*y3+((y8-y9)*y7**2+(-y9**2-y8**2)*y7+y8*y9**2-y8**2*y9)*y4
      dF(33) = (y9**3+y8**3+y7**3)*y1+(y9**3-y8**3-y7**3)*y2+(-y8**3-y9**3+y7**3)*y3+ &
          (-y7**3-y9**3+y8**3)*y4
      dF(34) = ((-sqrt(3._ark)*y7**2/3._ark-sqrt(3._ark)*y8**2/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**2)*y5+(y7**2-y8**2)*y6)*y1+((-sqrt(3._ark)*y7**2/3._ark- &
          sqrt(3._ark)*y8**2/3._ark+2._ark/3._ark*sqrt(3._ark)*y9**2)*y5+(y7**2-y8**2)*y6)*y2+((- &
          sqrt(3._ark)*y7**2/3._ark-sqrt(3._ark)*y8**2/3._ark+2._ark/3._ark*sqrt(3._ark)*y9**2)*y5+ &
          (y7**2-y8**2)*y6)*y3+((-sqrt(3._ark)*y7**2/3._ark-sqrt(3._ark)*y8**2/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**2)*y5+(y7**2-y8**2)*y6)*y4
      dF(35) = (((y8-y9/2._ark)*y7-y8*y9/2._ark)*y5+(-sqrt(3._ark)*y7*y9/2._ark+ &
          sqrt(3._ark)*y8*y9/2._ark)*y6)*y1+(((y8+y9/2._ark)*y7+y8*y9/2._ark)*y5+(- &
          sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y2+(((y9/2._ark-y8)*y7-y8*y9/ &
          2._ark)*y5+(sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y3+(((-y9/2._ark- &
          y8)*y7+y8*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y6)*y4
      dF(36) = (-sqrt(3._ark)*y5**2*y9/2._ark+(y7-y8)*y6*y5+(sqrt(3._ark)*y9/6._ark- &
          sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark)*y6**2)*y1+(-sqrt(3._ark)*y5**2*y9/2._ark+(y8- &
          y7)*y6*y5+(sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/3._ark)*y6**2)*y2+ &
          (sqrt(3._ark)*y5**2*y9/2._ark+(y7+y8)*y6*y5+(sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y7/3._ark- &
          sqrt(3._ark)*y9/6._ark)*y6**2)*y3+(sqrt(3._ark)*y5**2*y9/2._ark+(-y8-y7)*y6*y5+(- &
          sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark)*y6**2)*y4
      dF(37) = ((y9+y7+y8)*y5**2+(y9+y7+y8)*y6**2)*y1+((-y7+y9-y8)*y5**2+(-y7+y9- &
          y8)*y6**2)*y2+((-y8-y9+y7)*y5**2+(-y8-y9+y7)*y6**2)*y3+((-y7-y9+y8)*y5**2+(-y7- &
          y9+y8)*y6**2)*y4
      dF(38) = (y5**3-3._ark*y5*y6**2)*y1+(y5**3-3._ark*y5*y6**2)*y2+(y5**3- &
          3._ark*y5*y6**2)*y3+(y5**3-3._ark*y5*y6**2)*y4
      dF(39) = (y3*y7**2+y2*y9**2+y4*y8**2)*y1+(y4*y7**2+y3*y8**2)*y2+y3*y4*y9**2
      dF(40) = (y4*y7*y9+y3*y8*y9+y2*y7*y8)*y1+(-y4*y8*y9-y3*y7*y9)*y2-y3*y4*y7*y8
      dF(41) = ((y8**2+y7**2)*y2+(y9**2+y8**2)*y3+(y7**2+y9**2)*y4)*y1+((y7**2+ &
          y9**2)*y3+(y9**2+y8**2)*y4)*y2+(y8**2+y7**2)*y4*y3
      dF(42) = ((y6**2+y5**2)*y2+(y6**2+y5**2)*y3+(y6**2+y5**2)*y4)*y1+((y6**2+ &
          y5**2)*y3+(y6**2+y5**2)*y4)*y2+(y6**2+y5**2)*y4*y3
      dF(43) = ((sqrt(3._ark)*y6**2/2._ark-sqrt(3._ark)*y5**2/6._ark)*y2+(y5*y6+ &
          sqrt(3._ark)*y5**2/3._ark)*y3+(sqrt(3._ark)*y5**2/3._ark-y5*y6)*y4)*y1+ &
          ((sqrt(3._ark)*y5**2/3._ark-y5*y6)*y3+(y5*y6+sqrt(3._ark)*y5**2/3._ark)*y4)*y2+ &
          (sqrt(3._ark)*y6**2/2._ark-sqrt(3._ark)*y5**2/6._ark)*y4*y3
      dF(44) = (y2*y5*y9+(-y5*y7/2._ark+sqrt(3._ark)*y6*y7/2._ark)*y3+(-y5*y8/2._ark- &
          sqrt(3._ark)*y6*y8/2._ark)*y4)*y1+((y5*y8/2._ark+sqrt(3._ark)*y6*y8/2._ark)*y3+(y5*y7/ &
          2._ark-sqrt(3._ark)*y6*y7/2._ark)*y4)*y2-y3*y4*y5*y9
      dF(45) = (((y9+y7-y8)*y3+(y8+y9-y7)*y4)*y2+(y7-y9+y8)*y4*y3)*y1+(-y7-y8- &
          y9)*y4*y3*y2
      dF(46) = y1*y2*y3*y4
      dF(47) = (y3*y7+y4*y8+y2*y9)*y1**2+(y4**2*y8+y3**2*y7+y2**2*y9)*y1+(-y3*y8- &
          y4*y7)*y2**2+(-y3**2*y8-y4**2*y7)*y2-y3*y4**2*y9-y3**2*y4*y9
      dF(48) = ((-y8-y7)*y2+(-y9-y8)*y3+(-y9-y7)*y4)*y1**2+((y7+y8)*y2**2+(y9+ &
          y8)*y3**2+(y7+y9)*y4**2)*y1+((y7-y9)*y3+(y8-y9)*y4)*y2**2+((y9-y7)*y3**2+(-y8+ &
          y9)*y4**2)*y2+(y8-y7)*y4*y3**2+(y7-y8)*y4**2*y3
      dF(49) = (y2*y5+(-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y3+(-y5/2._ark-sqrt(3._ark)*y6/ &
          2._ark)*y4)*y1**2+(y2**2*y5+(-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y3**2+(-y5/2._ark- &
          sqrt(3._ark)*y6/2._ark)*y4**2)*y1+((-y5/2._ark-sqrt(3._ark)*y6/2._ark)*y3+(-y5/2._ark+ &
          sqrt(3._ark)*y6/2._ark)*y4)*y2**2+((-y5/2._ark-sqrt(3._ark)*y6/2._ark)*y3**2+(-y5/2._ark+ &
          sqrt(3._ark)*y6/2._ark)*y4**2)*y2+y3*y4**2*y5+y3**2*y4*y5
      dF(50) = ((y4+y3)*y2+y3*y4)*y1**2+((y4+y3)*y2**2+(y3**2+y4**2)*y2+y3*y4**2+ &
          y3**2*y4)*y1+y2**2*y3*y4+(y3*y4**2+y3**2*y4)*y2
      dF(51) = (y4+y3+y2)*y1**3+(y4**3+y3**3+y2**3)*y1+(y4+y3)*y2**3+(y4**3+y3**3)*y2+ &
          y3*y4**3+y3**3*y4
      dF(52) = ((y9+y8)*y7+y8*y9)*y1**2+((y8-y9)*y7-y8*y9)*y2**2+((-y9-y8)*y7+ &
          y8*y9)*y3**2+((-y8+y9)*y7-y8*y9)*y4**2
      dF(53) = (y8**2+y7**2+y9**2)*y1**2+(y8**2+y7**2+y9**2)*y2**2+(y8**2+y7**2+ &
          y9**2)*y3**2+(y8**2+y7**2+y9**2)*y4**2
      dF(54) = (y6**2+y5**2)*y1**2+(y6**2+y5**2)*y2**2+(y6**2+y5**2)*y3**2+(y6**2+ &
          y5**2)*y4**2
      dF(55) = ((y9-y8/2._ark-y7/2._ark)*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
          2._ark)*y6)*y1**2+((y9+y7/2._ark+y8/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/ &
          2._ark)*y6)*y2**2+((-y7/2._ark-y9+y8/2._ark)*y5+(sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/ &
          2._ark)*y6)*y3**2+((y7/2._ark-y8/2._ark-y9)*y5+(-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
          2._ark)*y6)*y4**2
      dF(56) = y7**3*y8*y9+(y8**3*y9+y8*y9**3)*y7
      dF(57) = (2._ark/3._ark*sqrt(3._ark)*y9**4-sqrt(3._ark)*y8**4/3._ark-sqrt(3._ark)*y7**4/ &
          3._ark)*y5+(-y8**4+y7**4)*y6
      dF(58) = sqrt(3._ark)*y5**3*y9**2+(y7**2-y8**2)*y6*y5**2+(-sqrt(3._ark)*y9**2/3._ark- &
         4._ark/3._ark*sqrt(3._ark)*y7**2-4._ark/3._ark*sqrt(3._ark)*y8**2)*y6**2*y5+(y7**2- &
         y8**2)*y6**3
      dF(59) = ((-y9**2/2._ark+y8**2)*y7**2-y8**2*y9**2/2._ark)*y5+ &
          (sqrt(3._ark)*y8**2*y9**2/2._ark-sqrt(3._ark)*y7**2*y9**2/2._ark)*y6
      dF(60) = y5**2*y7*y8*y9+y6**2*y7*y8*y9
      dF(61) = (y8**2+y7**2+y9**2)*y5**3+(-3._ark*y8**2-3._ark*y7**2-3._ark*y9**2)*y6**2*y5
      dF(62) = -3._ark*y5*y6**4+y5**5-2._ark*y5**3*y6**2
      s1 = (((y9/2._ark-y8)*y7**2+(y9**2/2._ark-y8**2)*y7+y8*y9**2/2._ark+y8**2*y9/2._ark)*y5+ &
          (-sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y8**2*y9/2._ark+sqrt(3._ark)*y7*y9**2/2._ark+ &
          sqrt(3._ark)*y7**2*y9/2._ark)*y6)*y1+(((y8+y9/2._ark)*y7**2+(-y9**2/2._ark+y8**2)*y7- &
          y8*y9**2/2._ark+y8**2*y9/2._ark)*y5+(-sqrt(3._ark)*y7*y9**2/2._ark+sqrt(3._ark)*y8*y9**2/ &
          2._ark+sqrt(3._ark)*y7**2*y9/2._ark-sqrt(3._ark)*y8**2*y9/2._ark)*y6)*y2
      dF(63) = s1+(((y8-y9/2._ark)*y7**2+(y9**2/2._ark-y8**2)*y7-y8**2*y9/2._ark-y8*y9**2/ &
          2._ark)*y5+(-sqrt(3._ark)*y7**2*y9/2._ark+sqrt(3._ark)*y8*y9**2/2._ark+ &
          sqrt(3._ark)*y7*y9**2/2._ark+sqrt(3._ark)*y8**2*y9/2._ark)*y6)*y3+(((-y9/2._ark-y8)*y7**2+ &
          (-y9**2/2._ark+y8**2)*y7-y8**2*y9/2._ark+y8*y9**2/2._ark)*y5+(sqrt(3._ark)*y8**2*y9/ &
          2._ark-sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7*y9**2/2._ark-sqrt(3._ark)*y7**2*y9/ &
          2._ark)*y6)*y4
      dF(64) = ((sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5**2+(y8*y9- &
          y7*y9)*y6*y5+((sqrt(3._ark)*y9/6._ark+2._ark/3._ark*sqrt(3._ark)*y8)*y7+sqrt(3._ark)*y8*y9/ &
          6._ark)*y6**2)*y1+((-sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y5**2+(-y8*y9+ &
          y7*y9)*y6*y5+((-sqrt(3._ark)*y9/6._ark+2._ark/3._ark*sqrt(3._ark)*y8)*y7-sqrt(3._ark)*y8*y9/ &
          6._ark)*y6**2)*y2+((-sqrt(3._ark)*y7*y9/2._ark+sqrt(3._ark)*y8*y9/2._ark)*y5**2+(y8*y9+ &
          y7*y9)*y6*y5+((-sqrt(3._ark)*y9/6._ark-2._ark/3._ark*sqrt(3._ark)*y8)*y7+sqrt(3._ark)*y8*y9/ &
          6._ark)*y6**2)*y3+((-sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5**2+(-y8*y9- &
          y7*y9)*y6*y5+((-2._ark/3._ark*sqrt(3._ark)*y8+sqrt(3._ark)*y9/6._ark)*y7-sqrt(3._ark)*y8*y9/ &
          6._ark)*y6**2)*y4
      dF(65) = (-sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/ &
          6._ark-sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2)*y1+(- &
          sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark- &
          sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2)*y2+(-sqrt(3._ark)*y5**2*y9**2/ &
          2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark-sqrt(3._ark)*y8**2/3._ark- &
          sqrt(3._ark)*y7**2/3._ark)*y6**2)*y3+(-sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2- &
          y8**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark-sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/ &
          3._ark)*y6**2)*y4
      dF(66) = (((y9+y8)*y7+y8*y9)*y5**2+((y9+y8)*y7+y8*y9)*y6**2)*y1+(((y8-y9)*y7- &
          y8*y9)*y5**2+((y8-y9)*y7-y8*y9)*y6**2)*y2+(((-y9-y8)*y7+y8*y9)*y5**2+((-y9- &
          y8)*y7+y8*y9)*y6**2)*y3+(((-y8+y9)*y7-y8*y9)*y5**2+((-y8+y9)*y7-y8*y9)*y6**2)*y4
      dF(67) = ((y8**2+y7**2+y9**2)*y5**2+(y8**2+y7**2+y9**2)*y6**2)*y1+((y8**2+y7**2+ &
          y9**2)*y5**2+(y8**2+y7**2+y9**2)*y6**2)*y2+((y8**2+y7**2+y9**2)*y5**2+(y8**2+ &
          y7**2+y9**2)*y6**2)*y3+((y8**2+y7**2+y9**2)*y5**2+(y8**2+y7**2+y9**2)*y6**2)*y4
      dF(68) = ((sqrt(3._ark)*y7+sqrt(3._ark)*y8)*y5**3+(y8-y7)*y6*y5**2+(-5._ark/ &
          3._ark*sqrt(3._ark)*y7-5._ark/3._ark*sqrt(3._ark)*y8-8._ark/3._ark*sqrt(3._ark)*y9)*y6**2*y5+ &
          (y8-y7)*y6**3)*y1+((-sqrt(3._ark)*y8-sqrt(3._ark)*y7)*y5**3+(y7-y8)*y6*y5**2+(-8._ark/ &
          3._ark*sqrt(3._ark)*y9+5._ark/3._ark*sqrt(3._ark)*y7+5._ark/3._ark*sqrt(3._ark)*y8)*y6**2*y5+ &
          (y7-y8)*y6**3)*y2+((sqrt(3._ark)*y7-sqrt(3._ark)*y8)*y5**3+(-y8-y7)*y6*y5**2+(8._ark/ &
          3._ark*sqrt(3._ark)*y9+5._ark/3._ark*sqrt(3._ark)*y8-5._ark/3._ark*sqrt(3._ark)*y7)*y6**2*y5+(- &
          y8-y7)*y6**3)*y3+((sqrt(3._ark)*y8-sqrt(3._ark)*y7)*y5**3+(y7+y8)*y6*y5**2+(5._ark/ &
          3._ark*sqrt(3._ark)*y7-5._ark/3._ark*sqrt(3._ark)*y8+8._ark/3._ark*sqrt(3._ark)*y9)*y6**2*y5+ &
          (y7+y8)*y6**3)*y4
      dF(69) = ((y9+y7+y8)*y5**3+(-3._ark*y8-3._ark*y9-3._ark*y7)*y6**2*y5)*y1+((-y7+y9- &
          y8)*y5**3+(3._ark*y7+3._ark*y8-3._ark*y9)*y6**2*y5)*y2+((-y8-y9+y7)*y5**3+(-3._ark*y7+ &
          3._ark*y8+3._ark*y9)*y6**2*y5)*y3+((-y7-y9+y8)*y5**3+(3._ark*y7+3._ark*y9- &
          3._ark*y8)*y6**2*y5)*y4
      dF(70) = (y6**4+y5**4+2._ark*y5**2*y6**2)*y1+(y6**4+y5**4+2._ark*y5**2*y6**2)*y2+ &
          (y6**4+y5**4+2._ark*y5**2*y6**2)*y3+(y6**4+y5**4+2._ark*y5**2*y6**2)*y4
      dF(71) = (y4*y7*y8*y9+y3*y7*y8*y9+y2*y7*y8*y9)*y1+(y4*y7*y8*y9+y3*y7*y8*y9)*y2+ &
          y3*y4*y7*y8*y9
      dF(72) = ((-y8**2*y9-y7**2*y9)*y2+(-y9**2-y8**2)*y7*y3+(-y8*y9**2- &
          y7**2*y8)*y4)*y1+((y7**2*y8+y8*y9**2)*y3+(y9**2+y8**2)*y7*y4)*y2+(y7**2*y9+ &
          y8**2*y9)*y4*y3
      dF(73) = (-y2*y9**3-y4*y8**3-y3*y7**3)*y1+(y4*y7**3+y3*y8**3)*y2+y3*y4*y9**3
      dF(74) = (((sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y7**2/2._ark)*y5+(y8**2/2._ark-y7**2/ &
          2._ark)*y6)*y2+(-sqrt(3._ark)*y5*y9**2/2._ark+(y9**2/2._ark+y8**2)*y6)*y3+(- &
          sqrt(3._ark)*y5*y9**2/2._ark+(-y7**2-y9**2/2._ark)*y6)*y4)*y1+((-sqrt(3._ark)*y5*y9**2/ &
          2._ark+(-y7**2-y9**2/2._ark)*y6)*y3+(-sqrt(3._ark)*y5*y9**2/2._ark+(y9**2/2._ark+ &
          y8**2)*y6)*y4)*y2+((sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y7**2/2._ark)*y5+(y8**2/2._ark- &
          y7**2/2._ark)*y6)*y4*y3
      dF(75) = (2._ark*y2*y5*y7*y8+(-y5*y8*y9+sqrt(3._ark)*y6*y8*y9)*y3+(- &
          sqrt(3._ark)*y6*y7*y9-y5*y7*y9)*y4)*y1+((y5*y7*y9+sqrt(3._ark)*y6*y7*y9)*y3+(- &
          sqrt(3._ark)*y6*y8*y9+y5*y8*y9)*y4)*y2-2._ark*y3*y4*y5*y7*y8
      dF(76) = (((-y8**2/2._ark-y7**2/2._ark)*y5+(sqrt(3._ark)*y8**2/2._ark-sqrt(3._ark)*y7**2/ &
          2._ark)*y6)*y2+((-y9**2/2._ark+y8**2)*y5-sqrt(3._ark)*y6*y9**2/2._ark)*y3+((y7**2-y9**2/ &
          2._ark)*y5+sqrt(3._ark)*y6*y9**2/2._ark)*y4)*y1+(((y7**2-y9**2/2._ark)*y5+ &
          sqrt(3._ark)*y6*y9**2/2._ark)*y3+((-y9**2/2._ark+y8**2)*y5-sqrt(3._ark)*y6*y9**2/ &
          2._ark)*y4)*y2+((-y8**2/2._ark-y7**2/2._ark)*y5+(sqrt(3._ark)*y8**2/2._ark- &
          sqrt(3._ark)*y7**2/2._ark)*y6)*y4*y3
      dF(77) = (-2._ark*y2*y5*y9**2+(-sqrt(3._ark)*y6*y7**2+y5*y7**2)*y3+(y5*y8**2+ &
          sqrt(3._ark)*y6*y8**2)*y4)*y1+((y5*y8**2+sqrt(3._ark)*y6*y8**2)*y3+(- &
          sqrt(3._ark)*y6*y7**2+y5*y7**2)*y4)*y2-2._ark*y3*y4*y5*y9**2
      dF(78) = ((-sqrt(3._ark)*y6**2*y9/6._ark+sqrt(3._ark)*y5**2*y9/2._ark)*y2+(-y5*y6*y7+ &
          sqrt(3._ark)*y6**2*y7/3._ark)*y3+(sqrt(3._ark)*y6**2*y8/3._ark+y5*y6*y8)*y4)*y1+((- &
          sqrt(3._ark)*y6**2*y8/3._ark-y5*y6*y8)*y3+(-sqrt(3._ark)*y6**2*y7/3._ark+ &
          y5*y6*y7)*y4)*y2+(sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y4*y3
      dF(79) = ((-y5**3/3._ark+y5*y6**2)*y2+(-y5**3/3._ark+y5*y6**2)*y3+(-y5**3/3._ark+ &
          y5*y6**2)*y4)*y1+((-y5**3/3._ark+y5*y6**2)*y3+(-y5**3/3._ark+y5*y6**2)*y4)*y2+(- &
          y5**3/3._ark+y5*y6**2)*y4*y3
      dF(80) = ((-y9*y6**2-y9*y5**2)*y2+(-y5**2*y7-y7*y6**2)*y3+(-y5**2*y8- &
          y8*y6**2)*y4)*y1+((y8*y6**2+y5**2*y8)*y3+(y7*y6**2+y5**2*y7)*y4)*y2+(y9*y6**2+ &
          y9*y5**2)*y4*y3
      dF(81) = ((5._ark/9._ark*sqrt(3._ark)*y5**3+sqrt(3._ark)*y5*y6**2)*y2+(y6**3+y5**2*y6- &
          4._ark/9._ark*sqrt(3._ark)*y5**3)*y3+(-y6**3-y5**2*y6-4._ark/ &
          9._ark*sqrt(3._ark)*y5**3)*y4)*y1+((-y6**3-y5**2*y6-4._ark/9._ark*sqrt(3._ark)*y5**3)*y3+ &
          (y6**3+y5**2*y6-4._ark/9._ark*sqrt(3._ark)*y5**3)*y4)*y2+(5._ark/9._ark*sqrt(3._ark)*y5**3+ &
          sqrt(3._ark)*y5*y6**2)*y4*y3
      dF(82) = (-y3*y8*y9-y2*y7*y8-y4*y7*y9)*y1**2+(-y4**2*y7*y9-y2**2*y7*y8- &
          y3**2*y8*y9)*y1+(y3*y7*y9+y4*y8*y9)*y2**2+(y4**2*y8*y9+y3**2*y7*y9)*y2+ &
          y3**2*y4*y7*y8+y3*y4**2*y7*y8
      dF(83) = ((y8**2+y7**2)*y2+(y9**2+y8**2)*y3+(y7**2+y9**2)*y4)*y1**2+((y8**2+ &
          y7**2)*y2**2+(y9**2+y8**2)*y3**2+(y7**2+y9**2)*y4**2)*y1+((y7**2+y9**2)*y3+ &
          (y9**2+y8**2)*y4)*y2**2+((y7**2+y9**2)*y3**2+(y9**2+y8**2)*y4**2)*y2+(y8**2+ &
          y7**2)*y4*y3**2+(y8**2+y7**2)*y4**2*y3
      dF(84) = ((-y8*y9-y7*y9)*y2+(-y9-y8)*y7*y3+(-y7*y8-y8*y9)*y4)*y1**2+((y8*y9+ &
          y7*y9)*y2**2+(y9+y8)*y7*y3**2+(y8*y9+y7*y8)*y4**2)*y1+((-y7*y8+y8*y9)*y3+(-y8+ &
          y9)*y7*y4)*y2**2+((-y8*y9+y7*y8)*y3**2+(y8-y9)*y7*y4**2)*y2+(-y8*y9+ &
          y7*y9)*y4*y3**2+(y8*y9-y7*y9)*y4**2*y3
      dF(85) = (y3*y7**2+y2*y9**2+y4*y8**2)*y1**2+(y4**2*y8**2+y3**2*y7**2+ &
          y2**2*y9**2)*y1+(y4*y7**2+y3*y8**2)*y2**2+(y4**2*y7**2+y3**2*y8**2)*y2+ &
          y3*y4**2*y9**2+y3**2*y4*y9**2
      s1 = (((sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y8/2._ark-y7/2._ark)*y6)*y2+(- &
          sqrt(3._ark)*y5*y9/2._ark+(y8+y9/2._ark)*y6)*y3+(-sqrt(3._ark)*y5*y9/2._ark+(-y9/2._ark- &
          y7)*y6)*y4)*y1**2+(((-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y5+(y7/2._ark-y8/ &
          2._ark)*y6)*y2**2+(sqrt(3._ark)*y5*y9/2._ark+(-y9/2._ark-y8)*y6)*y3**2+ &
          (sqrt(3._ark)*y5*y9/2._ark+(y7+y9/2._ark)*y6)*y4**2)*y1+((-sqrt(3._ark)*y5*y9/2._ark+(-y9/ &
         2._ark+y7)*y6)*y3+(-sqrt(3._ark)*y5*y9/2._ark+(y9/2._ark-y8)*y6)*y4)*y2**2
      dF(86) = s1+((sqrt(3._ark)*y5*y9/2._ark+(-y7+y9/2._ark)*y6)*y3**2+(sqrt(3._ark)*y5*y9/ &
          2._ark+(y8-y9/2._ark)*y6)*y4**2)*y2+((sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y5+(- &
          y8/2._ark-y7/2._ark)*y6)*y4*y3**2+((-sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y7/ &
          2._ark+y8/2._ark)*y6)*y4**2*y3
      dF(87) = ((y6**2+y5**2)*y2+(y6**2+y5**2)*y3+(y6**2+y5**2)*y4)*y1**2+((y6**2+ &
          y5**2)*y2**2+(y6**2+y5**2)*y3**2+(y6**2+y5**2)*y4**2)*y1+((y6**2+y5**2)*y3+ &
          (y6**2+y5**2)*y4)*y2**2+((y6**2+y5**2)*y3**2+(y6**2+y5**2)*y4**2)*y2+(y6**2+ &
          y5**2)*y4*y3**2+(y6**2+y5**2)*y4**2*y3
      s1 = (((-y8/2._ark-y7/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y6)*y2+ &
          ((y8-y9/2._ark)*y5-sqrt(3._ark)*y6*y9/2._ark)*y3+((-y9/2._ark+y7)*y5+sqrt(3._ark)*y6*y9/ &
          2._ark)*y4)*y1**2+(((y7/2._ark+y8/2._ark)*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
          2._ark)*y6)*y2**2+((y9/2._ark-y8)*y5+sqrt(3._ark)*y6*y9/2._ark)*y3**2+((-y7+y9/2._ark)*y5- &
          sqrt(3._ark)*y6*y9/2._ark)*y4**2)*y1+(((-y9/2._ark-y7)*y5+sqrt(3._ark)*y6*y9/2._ark)*y3+ &
          ((-y9/2._ark-y8)*y5-sqrt(3._ark)*y6*y9/2._ark)*y4)*y2**2
      dF(88) = s1+(((y7+y9/2._ark)*y5-sqrt(3._ark)*y6*y9/2._ark)*y3**2+((y8+y9/2._ark)*y5+ &
          sqrt(3._ark)*y6*y9/2._ark)*y4**2)*y2+((y8/2._ark-y7/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark- &
          sqrt(3._ark)*y8/2._ark)*y6)*y4*y3**2+((y7/2._ark-y8/2._ark)*y5+(sqrt(3._ark)*y7/2._ark+ &
          sqrt(3._ark)*y8/2._ark)*y6)*y4**2*y3
      dF(89) = ((((-y8+y9)*y7-y8*y9)*y3+((-y9-y8)*y7+y8*y9)*y4)*y2+((y8-y9)*y7- &
          y8*y9)*y4*y3)*y1+((y9+y8)*y7+y8*y9)*y4*y3*y2
      dF(90) = (((y8**2+y7**2+y9**2)*y3+(y8**2+y7**2+y9**2)*y4)*y2+(y8**2+y7**2+ &
          y9**2)*y4*y3)*y1+(y8**2+y7**2+y9**2)*y4*y3*y2
      dF(91) = (((y6**2+y5**2)*y3+(y6**2+y5**2)*y4)*y2+(y6**2+y5**2)*y4*y3)*y1+(y6**2+ &
          y5**2)*y4*y3*y2
      dF(92) = ((-sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y2+(-y5*y6- &
          sqrt(3._ark)*y5**2/3._ark)*y3+(y5*y6-sqrt(3._ark)*y5**2/3._ark)*y4)*y1**2+((- &
          sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y2**2+(-y5*y6-sqrt(3._ark)*y5**2/ &
          3._ark)*y3**2+(y5*y6-sqrt(3._ark)*y5**2/3._ark)*y4**2)*y1+((y5*y6-sqrt(3._ark)*y5**2/ &
          3._ark)*y3+(-y5*y6-sqrt(3._ark)*y5**2/3._ark)*y4)*y2**2+((y5*y6-sqrt(3._ark)*y5**2/ &
          3._ark)*y3**2+(-y5*y6-sqrt(3._ark)*y5**2/3._ark)*y4**2)*y2+(-sqrt(3._ark)*y6**2/2._ark+ &
          sqrt(3._ark)*y5**2/6._ark)*y4*y3**2+(-sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/ &
          6._ark)*y4**2*y3
      dF(93) = ((y3*y8+y4*y7)*y2+y3*y4*y9)*y1**2+((-y3*y7-y4*y8)*y2**2+(-y4**2*y9- &
          y3**2*y9)*y2-y3*y4**2*y7-y3**2*y4*y8)*y1+y2**2*y3*y4*y9+(y3*y4**2*y8+ &
          y3**2*y4*y7)*y2
      dF(94) = (((-sqrt(3._ark)*y5/3._ark-y6)*y3+(-sqrt(3._ark)*y5/3._ark+y6)*y4)*y2+2._ark/ &
          3._ark*sqrt(3._ark)*y3*y4*y5)*y1**2+(((-sqrt(3._ark)*y5/3._ark+y6)*y3+(-sqrt(3._ark)*y5/ &
          3._ark-y6)*y4)*y2**2+(2._ark/3._ark*sqrt(3._ark)*y4**2*y5+2._ark/ &
          3._ark*sqrt(3._ark)*y3**2*y5)*y2+(-sqrt(3._ark)*y5/3._ark-y6)*y4*y3**2+(-sqrt(3._ark)*y5/ &
          3._ark+y6)*y4**2*y3)*y1+2._ark/3._ark*sqrt(3._ark)*y2**2*y3*y4*y5+((-sqrt(3._ark)*y5/3._ark+ &
          y6)*y4*y3**2+(-sqrt(3._ark)*y5/3._ark-y6)*y4**2*y3)*y2
      dF(95) = ((y4+y3)*y2**2+(y3**2+y4**2)*y2+y3*y4**2+y3**2*y4)*y1**2+((y3**2+ &
          y4**2)*y2**2+y3**2*y4**2)*y1+(y3*y4**2+y3**2*y4)*y2**2+y2*y3**2*y4**2
      dF(96) = (-y4*y8-y3*y7-y2*y9)*y1**3+(-y4**3*y8-y3**3*y7-y2**3*y9)*y1+(y3*y8+ &
          y4*y7)*y2**3+(y3**3*y8+y4**3*y7)*y2+y3**3*y4*y9+y3*y4**3*y9
      dF(97) = ((y7+y8)*y2+(y9+y8)*y3+(y7+y9)*y4)*y1**3+((-y8-y7)*y2**3+(-y9- &
          y8)*y3**3+(-y9-y7)*y4**3)*y1+((y9-y7)*y3+(-y8+y9)*y4)*y2**3+((y7-y9)*y3**3+(y8- &
          y9)*y4**3)*y2+(y7-y8)*y4*y3**3+(y8-y7)*y4**3*y3
      dF(98) = (-2._ark/3._ark*sqrt(3._ark)*y2*y5+(-y6+sqrt(3._ark)*y5/3._ark)*y3+(y6+ &
          sqrt(3._ark)*y5/3._ark)*y4)*y1**3+(-2._ark/3._ark*sqrt(3._ark)*y2**3*y5+(-y6+ &
          sqrt(3._ark)*y5/3._ark)*y3**3+(y6+sqrt(3._ark)*y5/3._ark)*y4**3)*y1+((y6+sqrt(3._ark)*y5/ &
          3._ark)*y3+(-y6+sqrt(3._ark)*y5/3._ark)*y4)*y2**3+((y6+sqrt(3._ark)*y5/3._ark)*y3**3+(-y6+ &
          sqrt(3._ark)*y5/3._ark)*y4**3)*y2-2._ark/3._ark*sqrt(3._ark)*y3**3*y4*y5-2._ark/ &
          3._ark*sqrt(3._ark)*y3*y4**3*y5
      dF(99) = ((y4+y3)*y2+y3*y4)*y1**3+((y4+y3)*y2**3+(y4**3+y3**3)*y2+y3**3*y4+ &
          y3*y4**3)*y1+y2**3*y3*y4+(y3**3*y4+y3*y4**3)*y2
      dF(100) = (y4+y3+y2)*y1**4+(y4**4+y2**4+y3**4)*y1+(y4+y3)*y2**4+(y4**4+ &
          y3**4)*y2+y3*y4**4+y3**4*y4
      dF(101) = y2**2*y7*y8*y9+y3**2*y7*y8*y9+y1**2*y7*y8*y9+y4**2*y7*y8*y9
      dF(102) = ((-y9-y8)*y7**2+(-y9**2-y8**2)*y7-y8*y9**2-y8**2*y9)*y1**2+((y8- &
          y9)*y7**2+(y9**2+y8**2)*y7-y8**2*y9+y8*y9**2)*y2**2+((y9+y8)*y7**2+(-y9**2- &
          y8**2)*y7+y8**2*y9+y8*y9**2)*y3**2+((-y8+y9)*y7**2+(y9**2+y8**2)*y7-y8*y9**2+ &
          y8**2*y9)*y4**2
      dF(103) = (-y8**3-y7**3-y9**3)*y1**2+(-y9**3+y8**3+y7**3)*y2**2+(y8**3+y9**3- &
          y7**3)*y3**2+(y7**3+y9**3-y8**3)*y4**2
      dF(104) = (((y8-y9/2._ark)*y7-y8*y9/2._ark)*y5+(-sqrt(3._ark)*y7*y9/2._ark+ &
          sqrt(3._ark)*y8*y9/2._ark)*y6)*y1**2+(((y8+y9/2._ark)*y7+y8*y9/2._ark)*y5+(- &
          sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y2**2+(((y9/2._ark-y8)*y7-y8*y9/ &
          2._ark)*y5+(sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y3**2+(((-y9/2._ark- &
          y8)*y7+y8*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y6)*y4**2
      dF(105) = ((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y1**2+((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y2**2+((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y3**2+((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y4**2
      dF(106) = ((-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y5**2+(y8-y7)*y6*y5+(- &
          sqrt(3._ark)*y8/6._ark-2._ark/3._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y7/6._ark)*y6**2)*y1**2+ &
          ((sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5**2+(y7-y8)*y6*y5+(-2._ark/ &
          3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y8/6._ark+sqrt(3._ark)*y7/6._ark)*y6**2)*y2**2+((- &
          sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5**2+(-y8-y7)*y6*y5+(-sqrt(3._ark)*y7/ &
          6._ark+2._ark/3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y8/6._ark)*y6**2)*y3**2+((sqrt(3._ark)*y7/ &
          2._ark-sqrt(3._ark)*y8/2._ark)*y5**2+(y7+y8)*y6*y5+(2._ark/3._ark*sqrt(3._ark)*y9+ &
          sqrt(3._ark)*y7/6._ark-sqrt(3._ark)*y8/6._ark)*y6**2)*y4**2
      dF(107) = ((y9+y7+y8)*y5**2+(y9+y7+y8)*y6**2)*y1**2+((-y7+y9-y8)*y5**2+(-y7+y9- &
          y8)*y6**2)*y2**2+((-y8-y9+y7)*y5**2+(-y8-y9+y7)*y6**2)*y3**2+((-y7-y9+y8)*y5**2+ &
          (-y7-y9+y8)*y6**2)*y4**2
      dF(108) = (y5**3-3._ark*y5*y6**2)*y1**2+(y5**3-3._ark*y5*y6**2)*y2**2+(y5**3- &
          3._ark*y5*y6**2)*y3**2+(y5**3-3._ark*y5*y6**2)*y4**2
      dF(109) = (2._ark/3._ark*sqrt(3._ark)*y2*y5*y9+(-sqrt(3._ark)*y5*y7/3._ark+y6*y7)*y3+(- &
          sqrt(3._ark)*y5*y8/3._ark-y6*y8)*y4)*y1**2+(2._ark/3._ark*sqrt(3._ark)*y2**2*y5*y9+(- &
          sqrt(3._ark)*y5*y7/3._ark+y6*y7)*y3**2+(-sqrt(3._ark)*y5*y8/3._ark-y6*y8)*y4**2)*y1+ &
          ((sqrt(3._ark)*y5*y8/3._ark+y6*y8)*y3+(-y6*y7+sqrt(3._ark)*y5*y7/3._ark)*y4)*y2**2+ &
          ((sqrt(3._ark)*y5*y8/3._ark+y6*y8)*y3**2+(-y6*y7+sqrt(3._ark)*y5*y7/3._ark)*y4**2)*y2- &
          2._ark/3._ark*sqrt(3._ark)*y3*y4**2*y5*y9-2._ark/3._ark*sqrt(3._ark)*y3**2*y4*y5*y9
      dF(110) = (-y3**2*y7-y4**2*y8-y2**2*y9)*y1**2+(y4**2*y7+y3**2*y8)*y2**2+ &
          y3**2*y4**2*y9
      dF(111) = (-2._ark/3._ark*sqrt(3._ark)*y2**2*y5+(-y6+sqrt(3._ark)*y5/3._ark)*y3**2+(y6+ &
          sqrt(3._ark)*y5/3._ark)*y4**2)*y1**2+((y6+sqrt(3._ark)*y5/3._ark)*y3**2+(-y6+ &
          sqrt(3._ark)*y5/3._ark)*y4**2)*y2**2-2._ark/3._ark*sqrt(3._ark)*y3**2*y4**2*y5
      dF(112) = ((y9+y8)*y7+y8*y9)*y1**3+((y8-y9)*y7-y8*y9)*y2**3+((-y9-y8)*y7+ &
          y8*y9)*y3**3+((-y8+y9)*y7-y8*y9)*y4**3
      dF(113) = (y8**2+y7**2+y9**2)*y1**3+(y8**2+y7**2+y9**2)*y2**3+(y8**2+y7**2+ &
          y9**2)*y3**3+(y8**2+y7**2+y9**2)*y4**3
      dF(114) = ((-2._ark/3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/3._ark)*y5+ &
          (y8-y7)*y6)*y1**3+((-sqrt(3._ark)*y7/3._ark-2._ark/3._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y8/ &
          3._ark)*y5+(y7-y8)*y6)*y2**3+((sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9)*y5+(-y8-y7)*y6)*y3**3+((2._ark/3._ark*sqrt(3._ark)*y9+ &
          sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y7/3._ark)*y5+(y7+y8)*y6)*y4**3
      dF(115) = (y6**2+y5**2)*y1**3+(y6**2+y5**2)*y2**3+(y6**2+y5**2)*y3**3+(y6**2+ &
          y5**2)*y4**3
      dF(116) = (y3**2+y2**2+y4**2)*y1**3+(y4**3+y3**3+y2**3)*y1**2+(y3**2+ &
          y4**2)*y2**3+(y4**3+y3**3)*y2**2+y3**2*y4**3+y3**3*y4**2
      dF(117) = (-y7-y8-y9)*y1**4+(y7-y9+y8)*y2**4+(y8+y9-y7)*y3**4+(y9+y7-y8)*y4**4
      dF(118) = y4**5+y3**5+y2**5+y1**5
      dF(119) = (y7**2*y8*y9+(y8*y9**2+y8**2*y9)*y7)*y1+(-y7**2*y8*y9+(y8*y9**2- &
          y8**2*y9)*y7)*y2+(y7**2*y8*y9+(-y8**2*y9-y8*y9**2)*y7)*y3+(-y7**2*y8*y9+(- &
          y8*y9**2+y8**2*y9)*y7)*y4
      dF(120) = ((y9**2+y8**2)*y7**2+y8**2*y9**2)*y1+((y9**2+y8**2)*y7**2+ &
          y8**2*y9**2)*y2+((y9**2+y8**2)*y7**2+y8**2*y9**2)*y3+((y9**2+y8**2)*y7**2+ &
          y8**2*y9**2)*y4
      dF(121) = ((y9+y8)*y7**3+(y9**3+y8**3)*y7+y8**3*y9+y8*y9**3)*y1+((y8-y9)*y7**3+ &
          (y8**3-y9**3)*y7-y8*y9**3-y8**3*y9)*y2+((-y9-y8)*y7**3+(-y9**3-y8**3)*y7+ &
          y8**3*y9+y8*y9**3)*y3+((-y8+y9)*y7**3+(y9**3-y8**3)*y7-y8*y9**3-y8**3*y9)*y4
      dF(122) = (y9**4+y8**4+y7**4)*y1+(y9**4+y8**4+y7**4)*y2+(y9**4+y8**4+y7**4)*y3+ &
          (y9**4+y8**4+y7**4)*y4
      s1 = ((sqrt(3._ark)*y7*y9**2/2._ark-sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y7**2*y9/ &
          2._ark+sqrt(3._ark)*y8*y9**2/2._ark)*y5+((y8+y9/2._ark)*y7**2+(-y9**2/2._ark-y8**2)*y7- &
          y8**2*y9/2._ark+y8*y9**2/2._ark)*y6)*y1+((-sqrt(3._ark)*y7*y9**2/2._ark- &
          sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7**2*y9/2._ark-sqrt(3._ark)*y8**2*y9/2._ark)*y5+ &
          ((y9/2._ark-y8)*y7**2+(y9**2/2._ark+y8**2)*y7-y8*y9**2/2._ark-y8**2*y9/2._ark)*y6)*y2
      dF(123) = s1+((sqrt(3._ark)*y7**2*y9/2._ark+sqrt(3._ark)*y7*y9**2/2._ark+ &
          sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y8*y9**2/2._ark)*y5+((-y9/2._ark-y8)*y7**2+(- &
          y9**2/2._ark-y8**2)*y7+y8**2*y9/2._ark-y8*y9**2/2._ark)*y6)*y3+((sqrt(3._ark)*y8**2*y9/ &
          2._ark+sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7*y9**2/2._ark+sqrt(3._ark)*y7**2*y9/ &
          2._ark)*y5+((y8-y9/2._ark)*y7**2+(y9**2/2._ark+y8**2)*y7+y8**2*y9/2._ark+y8*y9**2/ &
          2._ark)*y6)*y4
      dF(124) = ((-sqrt(3._ark)*y7**3/3._ark-sqrt(3._ark)*y8**3/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**3)*y5+(-y8**3+y7**3)*y6)*y1+((sqrt(3._ark)*y7**3/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**3+sqrt(3._ark)*y8**3/3._ark)*y5+(-y7**3+y8**3)*y6)*y2+((-2._ark/ &
          3._ark*sqrt(3._ark)*y9**3-sqrt(3._ark)*y7**3/3._ark+sqrt(3._ark)*y8**3/3._ark)*y5+(y8**3+ &
          y7**3)*y6)*y3+((-2._ark/3._ark*sqrt(3._ark)*y9**3+sqrt(3._ark)*y7**3/3._ark- &
          sqrt(3._ark)*y8**3/3._ark)*y5+(-y7**3-y8**3)*y6)*y4
      dF(125) = ((((-sqrt(3._ark)*y8/3._ark-2._ark/3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y7/ &
          3._ark)*y5+(-y8-y7)*y6)*y3+((-sqrt(3._ark)*y7/3._ark+sqrt(3._ark)*y8/3._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y9)*y5+(y7+y8)*y6)*y4)*y2+((sqrt(3._ark)*y8/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y7/3._ark)*y5+(y8-y7)*y6)*y4*y3)*y1+((2._ark/ &
          3._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark)*y5+(y7- &
          y8)*y6)*y4*y3*y2
      dF(126) = (((y7+y9)*y3+(y9+y8)*y4)*y2+(y7+y8)*y4*y3)*y1**2+(((-y8+y9)*y3+(y9- &
          y7)*y4)*y2**2+((y7-y8)*y3**2+(y8-y7)*y4**2)*y2+(y7-y9)*y4*y3**2+(y8- &
          y9)*y4**2*y3)*y1+(-y8-y7)*y4*y3*y2**2+((-y9-y8)*y4*y3**2+(-y9-y7)*y4**2*y3)*y2
      dF(127) = y1**2*y2*y3*y4+(y2**2*y3*y4+(y3*y4**2+y3**2*y4)*y2)*y1
      dF(128) = (y9**2+y8**2)*y7**4+(y8**4+y9**4)*y7**2+y8**4*y9**2+y8**2*y9**4
      dF(129) = y7**6+y9**6+y8**6
      dF(130) = y7**2*y8**2*y9**2
      dF(131) = ((y9**2+y8**2)*y7**2+y8**2*y9**2)*y5**2+((y9**2+y8**2)*y7**2+ &
          y8**2*y9**2)*y6**2
      dF(132) = -6._ark*y5**2*y6**4+9._ark*y5**4*y6**2+y6**6
      dF(133) = (y7**3*y8*y9+(-2._ark*y8*y9**3+y8**3*y9)*y7)*y5+(sqrt(3._ark)*y7*y8**3*y9- &
          sqrt(3._ark)*y7**3*y8*y9)*y6
      dF(134) = ((2._ark/3._ark*sqrt(3._ark)*y8**2+sqrt(3._ark)*y9**2/6._ark)*y7**2+ &
          sqrt(3._ark)*y8**2*y9**2/6._ark)*y5**2+(-y8**2*y9**2+y7**2*y9**2)*y6*y5+ &
          (sqrt(3._ark)*y7**2*y9**2/2._ark+sqrt(3._ark)*y8**2*y9**2/2._ark)*y6**2
      dF(135) = -sqrt(3._ark)*y5**2*y9**4/2._ark+(-y8**4+y7**4)*y6*y5+(-sqrt(3._ark)*y7**4/ &
          3._ark-sqrt(3._ark)*y8**4/3._ark+sqrt(3._ark)*y9**4/6._ark)*y6**2
      dF(136) = -y5**3*y7*y8*y9/3._ark+y5*y6**2*y7*y8*y9
      dF(137) = (y9**4+y8**4+y7**4)*y5**2+(y9**4+y8**4+y7**4)*y6**2
      dF(138) = 3._ark/4._ark*y5**4*y9**2+(-y9**2/2._ark+y7**2+y8**2)*y6**2*y5**2+(-2._ark/ &
          3._ark*sqrt(3._ark)*y7**2+2._ark/3._ark*sqrt(3._ark)*y8**2)*y6**3*y5+(y8**2/3._ark+y9**2/ &
          12._ark+y7**2/3._ark)*y6**4
      dF(139) = -sqrt(3._ark)*y5**4*y9**2/4._ark+(y7**2-y8**2)*y6*y5**3- &
          sqrt(3._ark)*y5**2*y6**2*y9**2/2._ark+(-y8**2/3._ark+y7**2/3._ark)*y6**3*y5+(-2._ark/ &
          9._ark*sqrt(3._ark)*y8**2+7._ark/36._ark*sqrt(3._ark)*y9**2-2._ark/ &
          9._ark*sqrt(3._ark)*y7**2)*y6**4
      dF(140) = (-y9**2/2._ark+y7**2+y8**2)*y5**4+3._ark*y5**2*y6**2*y9**2+(4._ark/ &
          3._ark*sqrt(3._ark)*y7**2-4._ark/3._ark*sqrt(3._ark)*y8**2)*y6**3*y5+(y8**2/3._ark+5._ark/ &
          6._ark*y9**2+y7**2/3._ark)*y6**4
      dF(141) = 9._ark*y5**2*y6**4+y5**6-6._ark*y5**4*y6**2
      dF(142) = ((y9**2+y8**2)*y7**3+(y9**3+y8**3)*y7**2+y8**3*y9**2+y8**2*y9**3)*y1+ &
          ((-y9**2-y8**2)*y7**3+(y9**3-y8**3)*y7**2+y8**2*y9**3-y8**3*y9**2)*y2+((y9**2+ &
          y8**2)*y7**3+(-y9**3-y8**3)*y7**2-y8**2*y9**3-y8**3*y9**2)*y3+((-y9**2- &
          y8**2)*y7**3+(y8**3-y9**3)*y7**2+y8**3*y9**2-y8**2*y9**3)*y4
      dF(143) = ((y8*y9**2+y8**2*y9)*y7**2+y7*y8**2*y9**2)*y1+((-y8*y9**2+ &
          y8**2*y9)*y7**2-y7*y8**2*y9**2)*y2+((-y8**2*y9-y8*y9**2)*y7**2+ &
          y7*y8**2*y9**2)*y3+((y8*y9**2-y8**2*y9)*y7**2-y7*y8**2*y9**2)*y4
      dF(144) = (y7**3*y8*y9+(y8**3*y9+y8*y9**3)*y7)*y1+(y7**3*y8*y9+(y8**3*y9+ &
          y8*y9**3)*y7)*y2+(y7**3*y8*y9+(y8**3*y9+y8*y9**3)*y7)*y3+(y7**3*y8*y9+(y8**3*y9+ &
          y8*y9**3)*y7)*y4
      dF(145) = ((y9+y8)*y7**4+(y8**4+y9**4)*y7+y8**4*y9+y8*y9**4)*y1+((-y8+y9)*y7**4+ &
          (-y8**4-y9**4)*y7+y8**4*y9-y8*y9**4)*y2+((-y9-y8)*y7**4+(y8**4+y9**4)*y7- &
          y8**4*y9-y8*y9**4)*y3+((y8-y9)*y7**4+(-y8**4-y9**4)*y7+y8*y9**4-y8**4*y9)*y4
      dF(146) = (-y7**5-y8**5-y9**5)*y1+(y7**5-y9**5+y8**5)*y2+(-y7**5+y8**5+ &
          y9**5)*y3+(y7**5+y9**5-y8**5)*y4
      s1 = ((-sqrt(3._ark)*y7**3*y9/2._ark+sqrt(3._ark)*y7*y9**3/2._ark+sqrt(3._ark)*y8*y9**3/ &
          2._ark-sqrt(3._ark)*y8**3*y9/2._ark)*y5+((y8+y9/2._ark)*y7**3+(-y9**3/2._ark-y8**3)*y7- &
          y8**3*y9/2._ark+y8*y9**3/2._ark)*y6)*y1+((sqrt(3._ark)*y8**3*y9/2._ark+ &
          sqrt(3._ark)*y7**3*y9/2._ark-sqrt(3._ark)*y7*y9**3/2._ark-sqrt(3._ark)*y8*y9**3/2._ark)*y5+ &
          ((y8-y9/2._ark)*y7**3+(y9**3/2._ark-y8**3)*y7+y8**3*y9/2._ark-y8*y9**3/2._ark)*y6)*y2
      dF(147) = s1+((-sqrt(3._ark)*y7*y9**3/2._ark-sqrt(3._ark)*y8**3*y9/2._ark+ &
          sqrt(3._ark)*y7**3*y9/2._ark+sqrt(3._ark)*y8*y9**3/2._ark)*y5+((-y9/2._ark-y8)*y7**3+ &
          (y8**3+y9**3/2._ark)*y7-y8**3*y9/2._ark+y8*y9**3/2._ark)*y6)*y3+((sqrt(3._ark)*y8**3*y9/ &
          2._ark-sqrt(3._ark)*y7**3*y9/2._ark+sqrt(3._ark)*y7*y9**3/2._ark-sqrt(3._ark)*y8*y9**3/ &
          2._ark)*y5+((y9/2._ark-y8)*y7**3+(-y9**3/2._ark+y8**3)*y7+y8**3*y9/2._ark-y8*y9**3/ &
          2._ark)*y6)*y4
      dF(148) = ((y7**2*y8*y9/2._ark+(y8**2*y9/2._ark-y8*y9**2)*y7)*y5+ &
          (sqrt(3._ark)*y7*y8**2*y9/2._ark-sqrt(3._ark)*y7**2*y8*y9/2._ark)*y6)*y1+((-y7**2*y8*y9/ &
          2._ark+(-y8**2*y9/2._ark-y8*y9**2)*y7)*y5+(-sqrt(3._ark)*y7*y8**2*y9/2._ark+ &
          sqrt(3._ark)*y7**2*y8*y9/2._ark)*y6)*y2+((y7**2*y8*y9/2._ark+(-y8**2*y9/2._ark+ &
          y8*y9**2)*y7)*y5+(-sqrt(3._ark)*y7*y8**2*y9/2._ark-sqrt(3._ark)*y7**2*y8*y9/ &
          2._ark)*y6)*y3+((-y7**2*y8*y9/2._ark+(y8*y9**2+y8**2*y9/2._ark)*y7)*y5+ &
          (sqrt(3._ark)*y7*y8**2*y9/2._ark+sqrt(3._ark)*y7**2*y8*y9/2._ark)*y6)*y4
      dF(149) = ((-2._ark*y9**4+y7**4+y8**4)*y5+(sqrt(3._ark)*y8**4- &
          sqrt(3._ark)*y7**4)*y6)*y1+((-2._ark*y9**4+y7**4+y8**4)*y5+(sqrt(3._ark)*y8**4- &
          sqrt(3._ark)*y7**4)*y6)*y2+((-2._ark*y9**4+y7**4+y8**4)*y5+(sqrt(3._ark)*y8**4- &
          sqrt(3._ark)*y7**4)*y6)*y3+((-2._ark*y9**4+y7**4+y8**4)*y5+(sqrt(3._ark)*y8**4- &
          sqrt(3._ark)*y7**4)*y6)*y4
      dF(150) = ((y9+y8)*y7+y8*y9)*y1**4+((y8-y9)*y7-y8*y9)*y2**4+((-y9-y8)*y7+ &
          y8*y9)*y3**4+((-y8+y9)*y7-y8*y9)*y4**4
      dF(151) = (((-y9**2/2._ark+y8**2)*y7**2-y8**2*y9**2/2._ark)*y5+ &
          (sqrt(3._ark)*y8**2*y9**2/2._ark-sqrt(3._ark)*y7**2*y9**2/2._ark)*y6)*y1+(((-y9**2/2._ark+ &
          y8**2)*y7**2-y8**2*y9**2/2._ark)*y5+(sqrt(3._ark)*y8**2*y9**2/2._ark- &
          sqrt(3._ark)*y7**2*y9**2/2._ark)*y6)*y2+(((-y9**2/2._ark+y8**2)*y7**2-y8**2*y9**2/ &
          2._ark)*y5+(sqrt(3._ark)*y8**2*y9**2/2._ark-sqrt(3._ark)*y7**2*y9**2/2._ark)*y6)*y3+(((- &
          y9**2/2._ark+y8**2)*y7**2-y8**2*y9**2/2._ark)*y5+(sqrt(3._ark)*y8**2*y9**2/2._ark- &
          sqrt(3._ark)*y7**2*y9**2/2._ark)*y6)*y4
      s1 = (((y9/2._ark-y8)*y7**3+(y9**3/2._ark-y8**3)*y7+y8**3*y9/2._ark+y8*y9**3/2._ark)*y5+ &
          (-sqrt(3._ark)*y8*y9**3/2._ark+sqrt(3._ark)*y7*y9**3/2._ark+sqrt(3._ark)*y7**3*y9/2._ark- &
          sqrt(3._ark)*y8**3*y9/2._ark)*y6)*y1+(((-y9/2._ark-y8)*y7**3+(-y9**3/2._ark-y8**3)*y7- &
          y8**3*y9/2._ark-y8*y9**3/2._ark)*y5+(-sqrt(3._ark)*y7*y9**3/2._ark+sqrt(3._ark)*y8**3*y9/ &
          2._ark+sqrt(3._ark)*y8*y9**3/2._ark-sqrt(3._ark)*y7**3*y9/2._ark)*y6)*y2
      dF(152) = s1+(((y8-y9/2._ark)*y7**3+(-y9**3/2._ark+y8**3)*y7+y8**3*y9/2._ark+y8*y9**3/ &
          2._ark)*y5+(-sqrt(3._ark)*y8*y9**3/2._ark-sqrt(3._ark)*y7*y9**3/2._ark- &
          sqrt(3._ark)*y8**3*y9/2._ark-sqrt(3._ark)*y7**3*y9/2._ark)*y6)*y3+(((y8+y9/2._ark)*y7**3+ &
          (y8**3+y9**3/2._ark)*y7-y8**3*y9/2._ark-y8*y9**3/2._ark)*y5+(sqrt(3._ark)*y7*y9**3/2._ark+ &
          sqrt(3._ark)*y8**3*y9/2._ark+sqrt(3._ark)*y7**3*y9/2._ark+sqrt(3._ark)*y8*y9**3/ &
          2._ark)*y6)*y4
      s1 = (((-sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y8/3._ark)*y7**2+(sqrt(3._ark)*y8**2/3._ark+ &
          sqrt(3._ark)*y9**2/3._ark)*y7+sqrt(3._ark)*y8*y9**2/3._ark-sqrt(3._ark)*y8**2*y9/ &
          6._ark)*y5**2+(-y7**2*y8+(y9**2+y8**2)*y7-y8*y9**2)*y6*y5+(sqrt(3._ark)*y7**2*y9/ &
          2._ark+sqrt(3._ark)*y8**2*y9/2._ark)*y6**2)*y1+(((-sqrt(3._ark)*y9/6._ark-sqrt(3._ark)*y8/ &
          3._ark)*y7**2+(-sqrt(3._ark)*y9**2/3._ark-sqrt(3._ark)*y8**2/3._ark)*y7- &
          sqrt(3._ark)*y8*y9**2/3._ark-sqrt(3._ark)*y8**2*y9/6._ark)*y5**2+(y7**2*y8+(-y9**2- &
          y8**2)*y7+y8*y9**2)*y6*y5+(sqrt(3._ark)*y7**2*y9/2._ark+sqrt(3._ark)*y8**2*y9/ &
          2._ark)*y6**2)*y2
      dF(153) = s1+(((sqrt(3._ark)*y9/6._ark-sqrt(3._ark)*y8/3._ark)*y7**2+(sqrt(3._ark)*y8**2/ &
          3._ark+sqrt(3._ark)*y9**2/3._ark)*y7+sqrt(3._ark)*y8**2*y9/6._ark-sqrt(3._ark)*y8*y9**2/ &
          3._ark)*y5**2+(y7**2*y8+(y9**2+y8**2)*y7+y8*y9**2)*y6*y5+(-sqrt(3._ark)*y8**2*y9/ &
          2._ark-sqrt(3._ark)*y7**2*y9/2._ark)*y6**2)*y3+(((sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y9/ &
          6._ark)*y7**2+(-sqrt(3._ark)*y9**2/3._ark-sqrt(3._ark)*y8**2/3._ark)*y7+ &
          sqrt(3._ark)*y8*y9**2/3._ark+sqrt(3._ark)*y8**2*y9/6._ark)*y5**2+(-y7**2*y8+(-y9**2- &
          y8**2)*y7-y8*y9**2)*y6*y5+(-sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y7**2*y9/ &
          2._ark)*y6**2)*y4
      s1 = ((sqrt(3._ark)*y8/6._ark-5._ark/24._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y7/6._ark)*y5**4+(- &
          sqrt(3._ark)*y7/6._ark-sqrt(3._ark)*y8/6._ark+7._ark/12._ark*sqrt(3._ark)*y9)*y6**2*y5**2+(y7- &
          y8)*y6**3*y5+sqrt(3._ark)*y6**4*y9/8._ark)*y1+((-sqrt(3._ark)*y8/6._ark-5._ark/ &
          24._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y7/6._ark)*y5**4+(7._ark/12._ark*sqrt(3._ark)*y9+ &
          sqrt(3._ark)*y7/6._ark+sqrt(3._ark)*y8/6._ark)*y6**2*y5**2+(y8-y7)*y6**3*y5+ &
          sqrt(3._ark)*y6**4*y9/8._ark)*y2
      dF(154) = s1+((-sqrt(3._ark)*y8/6._ark+5._ark/24._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y7/ &
          6._ark)*y5**4+(-7._ark/12._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y7/6._ark+sqrt(3._ark)*y8/ &
          6._ark)*y6**2*y5**2+(y7+y8)*y6**3*y5-sqrt(3._ark)*y6**4*y9/8._ark)*y3+((sqrt(3._ark)*y8/ &
          6._ark-sqrt(3._ark)*y7/6._ark+5._ark/24._ark*sqrt(3._ark)*y9)*y5**4+(sqrt(3._ark)*y7/6._ark- &
          sqrt(3._ark)*y8/6._ark-7._ark/12._ark*sqrt(3._ark)*y9)*y6**2*y5**2+(-y8-y7)*y6**3*y5- &
          sqrt(3._ark)*y6**4*y9/8._ark)*y4
      dF(155) = (((-y9-y8)*y7-y8*y9)*y5**2+((-y9-y8)*y7-y8*y9)*y6**2)*y1**2+(((-y8+ &
          y9)*y7+y8*y9)*y5**2+((-y8+y9)*y7+y8*y9)*y6**2)*y2**2+(((y9+y8)*y7-y8*y9)*y5**2+ &
          ((y9+y8)*y7-y8*y9)*y6**2)*y3**2+(((y8-y9)*y7+y8*y9)*y5**2+((y8-y9)*y7+ &
          y8*y9)*y6**2)*y4**2
      dF(156) = (y6**4+y5**4+2._ark*y5**2*y6**2)*y1**2+(y6**4+y5**4+ &
          2._ark*y5**2*y6**2)*y2**2+(y6**4+y5**4+2._ark*y5**2*y6**2)*y3**2+(y6**4+y5**4+ &
          2._ark*y5**2*y6**2)*y4**2
      dF(157) = ((-y7**3-y8**3)*y2+(-y9**3-y8**3)*y3+(-y7**3-y9**3)*y4)*y1**2+((y8**3+ &
          y7**3)*y2**2+(y9**3+y8**3)*y3**2+(y9**3+y7**3)*y4**2)*y1+((-y9**3+y7**3)*y3+ &
          (y8**3-y9**3)*y4)*y2**2+((-y7**3+y9**3)*y3**2+(y9**3-y8**3)*y4**2)*y2+(-y7**3+ &
          y8**3)*y4*y3**2+(-y8**3+y7**3)*y4**2*y3
      dF(158) = ((y8*y9**2+y7*y9**2)*y2+(y9+y8)*y7**2*y3+(y8**2*y9+ &
          y7*y8**2)*y4)*y1**2+((-y8*y9**2-y7*y9**2)*y2**2+(-y9-y8)*y7**2*y3**2+(-y8**2*y9- &
          y7*y8**2)*y4**2)*y1+((-y7*y8**2+y8**2*y9)*y3+(-y8+y9)*y7**2*y4)*y2**2+((- &
          y8**2*y9+y7*y8**2)*y3**2+(y8-y9)*y7**2*y4**2)*y2+(-y8*y9**2+y7*y9**2)*y4*y3**2+ &
          (-y7*y9**2+y8*y9**2)*y4**2*y3
      dF(159) = (y4*y7*y8*y9+y3*y7*y8*y9+y2*y7*y8*y9)*y1**2+(y4**2*y7*y8*y9+ &
          y3**2*y7*y8*y9+y2**2*y7*y8*y9)*y1+(y4*y7*y8*y9+y3*y7*y8*y9)*y2**2+ &
          (y4**2*y7*y8*y9+y3**2*y7*y8*y9)*y2+y3**2*y4*y7*y8*y9+y3*y4**2*y7*y8*y9
      dF(160) = ((y7*y8**2+y7**2*y8)*y2+(y8*y9**2+y8**2*y9)*y3+(y7*y9**2+ &
          y7**2*y9)*y4)*y1**2+((-y7**2*y8-y7*y8**2)*y2**2+(-y8**2*y9-y8*y9**2)*y3**2+(- &
          y7*y9**2-y7**2*y9)*y4**2)*y1+((-y7*y9**2+y7**2*y9)*y3+(-y8*y9**2+ &
          y8**2*y9)*y4)*y2**2+((y7*y9**2-y7**2*y9)*y3**2+(y8*y9**2-y8**2*y9)*y4**2)*y2+ &
          (y7*y8**2-y7**2*y8)*y4*y3**2+(-y7*y8**2+y7**2*y8)*y4**2*y3
      dF(161) = ((-y8**2*y9-y7**2*y9)*y2+(-y9**2-y8**2)*y7*y3+(-y8*y9**2- &
          y7**2*y8)*y4)*y1**2+((-y8**2*y9-y7**2*y9)*y2**2+(-y9**2-y8**2)*y7*y3**2+(- &
          y8*y9**2-y7**2*y8)*y4**2)*y1+((y7**2*y8+y8*y9**2)*y3+(y9**2+y8**2)*y7*y4)*y2**2+ &
          ((y7**2*y8+y8*y9**2)*y3**2+(y9**2+y8**2)*y7*y4**2)*y2+(y7**2*y9+ &
          y8**2*y9)*y4*y3**2+(y7**2*y9+y8**2*y9)*y4**2*y3
      dF(162) = (y4**2*y8**2+y3**2*y7**2+y2**2*y9**2)*y1**2+(y4**2*y7**2+ &
          y3**2*y8**2)*y2**2+y3**2*y4**2*y9**2
      dF(163) = ((y8**2+y7**2)*y2**2+(y9**2+y8**2)*y3**2+(y7**2+y9**2)*y4**2)*y1**2+ &
          ((y7**2+y9**2)*y3**2+(y9**2+y8**2)*y4**2)*y2**2+(y8**2+y7**2)*y4**2*y3**2
      dF(164) = ((-y9-y8)*y7**2+(-y9**2-y8**2)*y7-y8*y9**2-y8**2*y9)*y1**3+((y8- &
          y9)*y7**2+(y9**2+y8**2)*y7-y8**2*y9+y8*y9**2)*y2**3+((y9+y8)*y7**2+(-y9**2- &
          y8**2)*y7+y8**2*y9+y8*y9**2)*y3**3+((-y8+y9)*y7**2+(y9**2+y8**2)*y7-y8*y9**2+ &
          y8**2*y9)*y4**3
      dF(165) = (-y8**3-y7**3-y9**3)*y1**3+(-y9**3+y8**3+y7**3)*y2**3+(y8**3+y9**3- &
          y7**3)*y3**3+(y7**3+y9**3-y8**3)*y4**3
      dF(166) = y4**3*y7*y8*y9+y2**3*y7*y8*y9+y1**3*y7*y8*y9+y3**3*y7*y8*y9
      dF(167) = ((sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y7**2/3._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y9**2)*y5+(-y7**2+y8**2)*y6)*y1**3+((sqrt(3._ark)*y8**2/3._ark+ &
          sqrt(3._ark)*y7**2/3._ark-2._ark/3._ark*sqrt(3._ark)*y9**2)*y5+(-y7**2+y8**2)*y6)*y2**3+ &
          ((sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y7**2/3._ark-2._ark/3._ark*sqrt(3._ark)*y9**2)*y5+(- &
          y7**2+y8**2)*y6)*y3**3+((sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y7**2/3._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y9**2)*y5+(-y7**2+y8**2)*y6)*y4**3
      dF(168) = ((2._ark/3._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark)*y5+ &
          (y7-y8)*y6)*y1**4+((sqrt(3._ark)*y8/3._ark+2._ark/3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y7/ &
          3._ark)*y5+(y8-y7)*y6)*y2**4+((-sqrt(3._ark)*y7/3._ark+sqrt(3._ark)*y8/3._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y9)*y5+(y7+y8)*y6)*y3**4+((-sqrt(3._ark)*y8/3._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y7/3._ark)*y5+(-y8-y7)*y6)*y4**4
      s1 = (((-sqrt(3._ark)*y9/3._ark-sqrt(3._ark)*y8/3._ark)*y7**2+(-sqrt(3._ark)*y8**2/3._ark+ &
          sqrt(3._ark)*y9**2/6._ark)*y7+sqrt(3._ark)*y8*y9**2/6._ark-sqrt(3._ark)*y8**2*y9/ &
          3._ark)*y5**2+((-y9-y8)*y7**2+y7*y8**2+y8**2*y9)*y6*y5+(-sqrt(3._ark)*y8*y9**2/2._ark- &
          sqrt(3._ark)*y7*y9**2/2._ark)*y6**2)*y1+(((-sqrt(3._ark)*y9/3._ark+sqrt(3._ark)*y8/ &
          3._ark)*y7**2+(sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y9**2/6._ark)*y7- &
          sqrt(3._ark)*y8**2*y9/3._ark-sqrt(3._ark)*y8*y9**2/6._ark)*y5**2+((y8-y9)*y7**2- &
          y7*y8**2+y8**2*y9)*y6*y5+(sqrt(3._ark)*y7*y9**2/2._ark+sqrt(3._ark)*y8*y9**2/ &
          2._ark)*y6**2)*y2
      dF(169) = s1+(((sqrt(3._ark)*y9/3._ark+sqrt(3._ark)*y8/3._ark)*y7**2+(-sqrt(3._ark)*y8**2/ &
          3._ark+sqrt(3._ark)*y9**2/6._ark)*y7-sqrt(3._ark)*y8*y9**2/6._ark+sqrt(3._ark)*y8**2*y9/ &
          3._ark)*y5**2+((y9+y8)*y7**2+y7*y8**2-y8**2*y9)*y6*y5+(sqrt(3._ark)*y8*y9**2/2._ark- &
          sqrt(3._ark)*y7*y9**2/2._ark)*y6**2)*y3+(((-sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y9/ &
          3._ark)*y7**2+(sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y9**2/6._ark)*y7+ &
          sqrt(3._ark)*y8*y9**2/6._ark+sqrt(3._ark)*y8**2*y9/3._ark)*y5**2+((-y8+y9)*y7**2- &
          y7*y8**2-y8**2*y9)*y6*y5+(-sqrt(3._ark)*y8*y9**2/2._ark+sqrt(3._ark)*y7*y9**2/ &
          2._ark)*y6**2)*y4
      dF(170) = ((-sqrt(3._ark)*y9**3/6._ark+sqrt(3._ark)*y7**3/3._ark+sqrt(3._ark)*y8**3/ &
          3._ark)*y5**2+(-y8**3+y7**3)*y6*y5+sqrt(3._ark)*y6**2*y9**3/2._ark)*y1+((- &
          sqrt(3._ark)*y9**3/6._ark-sqrt(3._ark)*y8**3/3._ark-sqrt(3._ark)*y7**3/3._ark)*y5**2+(- &
          y7**3+y8**3)*y6*y5+sqrt(3._ark)*y6**2*y9**3/2._ark)*y2+((sqrt(3._ark)*y7**3/3._ark- &
          sqrt(3._ark)*y8**3/3._ark+sqrt(3._ark)*y9**3/6._ark)*y5**2+(y8**3+y7**3)*y6*y5- &
          sqrt(3._ark)*y6**2*y9**3/2._ark)*y3+((-sqrt(3._ark)*y7**3/3._ark+sqrt(3._ark)*y8**3/3._ark+ &
          sqrt(3._ark)*y9**3/6._ark)*y5**2+(-y7**3-y8**3)*y6*y5-sqrt(3._ark)*y6**2*y9**3/ &
          2._ark)*y4
      dF(171) = ((-y8**2/3._ark-y7**2/3._ark-y9**2/3._ark)*y5**3+(y8**2+y7**2+ &
          y9**2)*y6**2*y5)*y1+((-y8**2/3._ark-y7**2/3._ark-y9**2/3._ark)*y5**3+(y8**2+y7**2+ &
          y9**2)*y6**2*y5)*y2+((-y8**2/3._ark-y7**2/3._ark-y9**2/3._ark)*y5**3+(y8**2+y7**2+ &
          y9**2)*y6**2*y5)*y3+((-y8**2/3._ark-y7**2/3._ark-y9**2/3._ark)*y5**3+(y8**2+y7**2+ &
          y9**2)*y6**2*y5)*y4
      dF(172) = (((y8/3._ark+y9/3._ark)*y7+y8*y9/3._ark)*y5**3+((-y9-y8)*y7- &
          y8*y9)*y6**2*y5)*y1+(((y8/3._ark-y9/3._ark)*y7-y8*y9/3._ark)*y5**3+((-y8+y9)*y7+ &
          y8*y9)*y6**2*y5)*y2+(((-y8/3._ark-y9/3._ark)*y7+y8*y9/3._ark)*y5**3+((y9+y8)*y7- &
          y8*y9)*y6**2*y5)*y3+(((y9/3._ark-y8/3._ark)*y7-y8*y9/3._ark)*y5**3+((y8-y9)*y7+ &
          y8*y9)*y6**2*y5)*y4
      dF(173) = ((y7**2*y9**2+y8**2*y9**2)*y2+(y9**2+y8**2)*y7**2*y3+(y8**2*y9**2+ &
          y7**2*y8**2)*y4)*y1+((y8**2*y9**2+y7**2*y8**2)*y3+(y9**2+y8**2)*y7**2*y4)*y2+ &
          (y7**2*y9**2+y8**2*y9**2)*y4*y3
      dF(174) = (((sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y7**2/2._ark)*y5**2+(- &
          sqrt(3._ark)*y8**2/6._ark-sqrt(3._ark)*y7**2/6._ark)*y6**2)*y2+((-y9**2-y8**2)*y6*y5+ &
          (sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y9**2/3._ark)*y6**2)*y3+((y7**2+y9**2)*y6*y5+ &
          (sqrt(3._ark)*y7**2/3._ark+sqrt(3._ark)*y9**2/3._ark)*y6**2)*y4)*y1+(((y7**2+ &
          y9**2)*y6*y5+(sqrt(3._ark)*y7**2/3._ark+sqrt(3._ark)*y9**2/3._ark)*y6**2)*y3+((-y9**2- &
          y8**2)*y6*y5+(sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y9**2/3._ark)*y6**2)*y4)*y2+ &
          ((sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y7**2/2._ark)*y5**2+(-sqrt(3._ark)*y8**2/6._ark- &
          sqrt(3._ark)*y7**2/6._ark)*y6**2)*y4*y3
      dF(175) = (((y8**2+y7**2)*y5**2+(y8**2+y7**2)*y6**2)*y2+((y9**2+y8**2)*y5**2+ &
          (y9**2+y8**2)*y6**2)*y3+((y7**2+y9**2)*y5**2+(y7**2+y9**2)*y6**2)*y4)*y1+ &
          (((y7**2+y9**2)*y5**2+(y7**2+y9**2)*y6**2)*y3+((y9**2+y8**2)*y5**2+(y9**2+ &
          y8**2)*y6**2)*y4)*y2+((y8**2+y7**2)*y5**2+(y8**2+y7**2)*y6**2)*y4*y3
      s1 = (((-sqrt(3._ark)*y8**2-sqrt(3._ark)*y7**2)*y5+(-y7**2+y8**2)*y6)*y2+ &
          (sqrt(3._ark)*y5*y8**2+(-y8**2-2._ark*y9**2)*y6)*y3+(sqrt(3._ark)*y5*y7**2+ &
          (2._ark*y9**2+y7**2)*y6)*y4)*y1**2+(((-sqrt(3._ark)*y8**2-sqrt(3._ark)*y7**2)*y5+(- &
          y7**2+y8**2)*y6)*y2**2+(sqrt(3._ark)*y5*y8**2+(-y8**2-2._ark*y9**2)*y6)*y3**2+ &
          (sqrt(3._ark)*y5*y7**2+(2._ark*y9**2+y7**2)*y6)*y4**2)*y1+((sqrt(3._ark)*y5*y7**2+ &
          (2._ark*y9**2+y7**2)*y6)*y3+(sqrt(3._ark)*y5*y8**2+(-y8**2-2._ark*y9**2)*y6)*y4)*y2**2
      dF(176) = s1+((sqrt(3._ark)*y5*y7**2+(2._ark*y9**2+y7**2)*y6)*y3**2+ &
          (sqrt(3._ark)*y5*y8**2+(-y8**2-2._ark*y9**2)*y6)*y4**2)*y2+((-sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y5+(-y7**2+y8**2)*y6)*y4*y3**2+((-sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y5+(-y7**2+y8**2)*y6)*y4**2*y3
      dF(177) = (y4+y3+y2)*y1**5+(y2**5+y4**5+y3**5)*y1+(y4+y3)*y2**5+(y4**5+ &
          y3**5)*y2+y3**5*y4+y3*y4**5
      s2 = (((-sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark)*y5**2+(y7-y8)*y6*y5)*y2+ &
          ((sqrt(3._ark)*y8/6._ark-sqrt(3._ark)*y9/3._ark)*y5**2+y5*y6*y9-sqrt(3._ark)*y6**2*y8/ &
          2._ark)*y3+((sqrt(3._ark)*y7/6._ark-sqrt(3._ark)*y9/3._ark)*y5**2-y5*y6*y9- &
          sqrt(3._ark)*y6**2*y7/2._ark)*y4)*y1**2
      s3 = (((sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/3._ark)*y5**2+(y8-y7)*y6*y5)*y2**2+((- &
          sqrt(3._ark)*y8/6._ark+sqrt(3._ark)*y9/3._ark)*y5**2-y5*y6*y9+sqrt(3._ark)*y6**2*y8/ &
          2._ark)*y3**2+((-sqrt(3._ark)*y7/6._ark+sqrt(3._ark)*y9/3._ark)*y5**2+y5*y6*y9+ &
          sqrt(3._ark)*y6**2*y7/2._ark)*y4**2)*y1+(((-sqrt(3._ark)*y9/3._ark-sqrt(3._ark)*y7/ &
          6._ark)*y5**2-y5*y6*y9+sqrt(3._ark)*y6**2*y7/2._ark)*y3+((-sqrt(3._ark)*y9/3._ark- &
          sqrt(3._ark)*y8/6._ark)*y5**2+y5*y6*y9+sqrt(3._ark)*y6**2*y8/2._ark)*y4)*y2**2
      s1 = s2+s3
      dF(178) = s1+(((sqrt(3._ark)*y7/6._ark+sqrt(3._ark)*y9/3._ark)*y5**2+y5*y6*y9- &
          sqrt(3._ark)*y6**2*y7/2._ark)*y3**2+((sqrt(3._ark)*y8/6._ark+sqrt(3._ark)*y9/3._ark)*y5**2- &
          y5*y6*y9-sqrt(3._ark)*y6**2*y8/2._ark)*y4**2)*y2+((sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y7/ &
          3._ark)*y5**2+(y7+y8)*y6*y5)*y4*y3**2+((sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/ &
          3._ark)*y5**2+(-y8-y7)*y6*y5)*y4**2*y3
      s2 = (((-sqrt(3._ark)*y7/6._ark-sqrt(3._ark)*y8/6._ark)*y5**2+(y7-y8)*y6*y5+(- &
          sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y6**2)*y2+((-2._ark/3._ark*sqrt(3._ark)*y9- &
          sqrt(3._ark)*y8/6._ark)*y5**2-y5*y6*y8-sqrt(3._ark)*y6**2*y8/2._ark)*y3+((-2._ark/ &
          3._ark*sqrt(3._ark)*y9-sqrt(3._ark)*y7/6._ark)*y5**2+y5*y6*y7-sqrt(3._ark)*y6**2*y7/ &
          2._ark)*y4)*y1**2
      s3 = (((sqrt(3._ark)*y8/6._ark+sqrt(3._ark)*y7/6._ark)*y5**2+(y8-y7)*y6*y5+ &
          (sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y6**2)*y2**2+((2._ark/3._ark*sqrt(3._ark)*y9+ &
          sqrt(3._ark)*y8/6._ark)*y5**2+y5*y6*y8+sqrt(3._ark)*y6**2*y8/2._ark)*y3**2+ &
          ((sqrt(3._ark)*y7/6._ark+2._ark/3._ark*sqrt(3._ark)*y9)*y5**2-y5*y6*y7+ &
          sqrt(3._ark)*y6**2*y7/2._ark)*y4**2)*y1+(((sqrt(3._ark)*y7/6._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y9)*y5**2-y5*y6*y7+sqrt(3._ark)*y6**2*y7/2._ark)*y3+((sqrt(3._ark)*y8/ &
          6._ark-2._ark/3._ark*sqrt(3._ark)*y9)*y5**2+y5*y6*y8+sqrt(3._ark)*y6**2*y8/2._ark)*y4)*y2**2
      s1 = s2+s3
      dF(179) = s1+(((-sqrt(3._ark)*y7/6._ark+2._ark/3._ark*sqrt(3._ark)*y9)*y5**2+y5*y6*y7- &
          sqrt(3._ark)*y6**2*y7/2._ark)*y3**2+((-sqrt(3._ark)*y8/6._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9)*y5**2-y5*y6*y8-sqrt(3._ark)*y6**2*y8/2._ark)*y4**2)*y2+((- &
          sqrt(3._ark)*y7/6._ark+sqrt(3._ark)*y8/6._ark)*y5**2+(y7+y8)*y6*y5+(-sqrt(3._ark)*y7/2._ark+ &
          sqrt(3._ark)*y8/2._ark)*y6**2)*y4*y3**2+((sqrt(3._ark)*y7/6._ark-sqrt(3._ark)*y8/ &
          6._ark)*y5**2+(-y8-y7)*y6*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
          2._ark)*y6**2)*y4**2*y3
      dF(180) = ((y4*y7**2+y3*y8**2)*y2+y3*y4*y9**2)*y1**2+((y4*y8**2+y3*y7**2)*y2**2+ &
          (y4**2*y9**2+y3**2*y9**2)*y2+y3**2*y4*y8**2+y3*y4**2*y7**2)*y1+ &
          y2**2*y3*y4*y9**2+(y3*y4**2*y8**2+y3**2*y4*y7**2)*y2
      dF(181) = (((-y7*y8-y8*y9)*y3+(-y9-y8)*y7*y4)*y2+(-y8*y9-y7*y9)*y4*y3)*y1**2+ &
          (((-y8+y9)*y7*y3+(-y7*y8+y8*y9)*y4)*y2**2+((-y8*y9+y7*y9)*y3**2+(y8*y9- &
          y7*y9)*y4**2)*y2+(-y8*y9+y7*y8)*y4*y3**2+(y8-y9)*y7*y4**2*y3)*y1+(y8*y9+ &
          y7*y9)*y4*y3*y2**2+((y9+y8)*y7*y4*y3**2+(y8*y9+y7*y8)*y4**2*y3)*y2
      dF(182) = (y5**2*y7*y8*y9+y6**2*y7*y8*y9)*y1+(y5**2*y7*y8*y9+y6**2*y7*y8*y9)*y2+ &
          (y5**2*y7*y8*y9+y6**2*y7*y8*y9)*y3+(y5**2*y7*y8*y9+y6**2*y7*y8*y9)*y4
      dF(183) = (((y9+y8)*y7**2+(y9**2+y8**2)*y7+y8*y9**2+y8**2*y9)*y5**2+((y9+ &
          y8)*y7**2+(y9**2+y8**2)*y7+y8*y9**2+y8**2*y9)*y6**2)*y1+(((-y8+y9)*y7**2+(- &
          y9**2-y8**2)*y7+y8**2*y9-y8*y9**2)*y5**2+((-y8+y9)*y7**2+(-y9**2-y8**2)*y7+ &
          y8**2*y9-y8*y9**2)*y6**2)*y2+(((-y9-y8)*y7**2+(y9**2+y8**2)*y7-y8**2*y9- &
          y8*y9**2)*y5**2+((-y9-y8)*y7**2+(y9**2+y8**2)*y7-y8**2*y9-y8*y9**2)*y6**2)*y3+ &
          (((y8-y9)*y7**2+(-y9**2-y8**2)*y7+y8*y9**2-y8**2*y9)*y5**2+((y8-y9)*y7**2+(- &
          y9**2-y8**2)*y7+y8*y9**2-y8**2*y9)*y6**2)*y4
      dF(184) = ((y9**3+y8**3+y7**3)*y5**2+(y9**3+y8**3+y7**3)*y6**2)*y1+((y9**3- &
          y8**3-y7**3)*y5**2+(y9**3-y8**3-y7**3)*y6**2)*y2+((-y8**3-y9**3+y7**3)*y5**2+(- &
          y8**3-y9**3+y7**3)*y6**2)*y3+((-y7**3-y9**3+y8**3)*y5**2+(-y7**3-y9**3+ &
          y8**3)*y6**2)*y4
      dF(185) = (((4._ark/9._ark*sqrt(3._ark)*y9-5._ark/9._ark*sqrt(3._ark)*y8)*y7+4._ark/ &
          9._ark*sqrt(3._ark)*y8*y9)*y5**3+(-y8*y9+y7*y9)*y6*y5**2-sqrt(3._ark)*y5*y6**2*y7*y8+ &
          (-y8*y9+y7*y9)*y6**3)*y1+(((-5._ark/9._ark*sqrt(3._ark)*y8-4._ark/ &
          9._ark*sqrt(3._ark)*y9)*y7-4._ark/9._ark*sqrt(3._ark)*y8*y9)*y5**3+(y8*y9-y7*y9)*y6*y5**2- &
          sqrt(3._ark)*y5*y6**2*y7*y8+(y8*y9-y7*y9)*y6**3)*y2+(((-4._ark/9._ark*sqrt(3._ark)*y9+ &
          5._ark/9._ark*sqrt(3._ark)*y8)*y7+4._ark/9._ark*sqrt(3._ark)*y8*y9)*y5**3+(-y8*y9- &
          y7*y9)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y7*y8+(-y8*y9-y7*y9)*y6**3)*y3+(((4._ark/ &
          9._ark*sqrt(3._ark)*y9+5._ark/9._ark*sqrt(3._ark)*y8)*y7-4._ark/ &
          9._ark*sqrt(3._ark)*y8*y9)*y5**3+(y8*y9+y7*y9)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y7*y8+ &
          (y8*y9+y7*y9)*y6**3)*y4
      dF(186) = ((-4._ark/9._ark*sqrt(3._ark)*y8**2+5._ark/9._ark*sqrt(3._ark)*y9**2-4._ark/ &
          9._ark*sqrt(3._ark)*y7**2)*y5**3+(y7**2-y8**2)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y9**2+ &
          (y7**2-y8**2)*y6**3)*y1+((-4._ark/9._ark*sqrt(3._ark)*y8**2+5._ark/ &
          9._ark*sqrt(3._ark)*y9**2-4._ark/9._ark*sqrt(3._ark)*y7**2)*y5**3+(y7**2-y8**2)*y6*y5**2+ &
          sqrt(3._ark)*y5*y6**2*y9**2+(y7**2-y8**2)*y6**3)*y2+((-4._ark/9._ark*sqrt(3._ark)*y8**2+ &
          5._ark/9._ark*sqrt(3._ark)*y9**2-4._ark/9._ark*sqrt(3._ark)*y7**2)*y5**3+(y7**2- &
          y8**2)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y9**2+(y7**2-y8**2)*y6**3)*y3+((-4._ark/ &
          9._ark*sqrt(3._ark)*y8**2+5._ark/9._ark*sqrt(3._ark)*y9**2-4._ark/ &
          9._ark*sqrt(3._ark)*y7**2)*y5**3+(y7**2-y8**2)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y9**2+ &
          (y7**2-y8**2)*y6**3)*y4
      dF(187) = ((y9+y7+y8)*y5**4+(2._ark*y8+2._ark*y9+2._ark*y7)*y6**2*y5**2+(y9+y7+ &
          y8)*y6**4)*y1+((-y7+y9-y8)*y5**4+(-2._ark*y7+2._ark*y9-2._ark*y8)*y6**2*y5**2+(-y7+y9- &
          y8)*y6**4)*y2+((-y8-y9+y7)*y5**4+(-2._ark*y8+2._ark*y7-2._ark*y9)*y6**2*y5**2+(-y8-y9+ &
          y7)*y6**4)*y3+((-y7-y9+y8)*y5**4+(-2._ark*y9+2._ark*y8-2._ark*y7)*y6**2*y5**2+(-y7-y9+ &
          y8)*y6**4)*y4
      s1 = ((sqrt(3._ark)*y9/24._ark+sqrt(3._ark)*y8/6._ark+sqrt(3._ark)*y7/6._ark)*y5**4+(y7- &
          y8)*y6*y5**3+(sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y9/4._ark+sqrt(3._ark)*y7/ &
          2._ark)*y6**2*y5**2+3._ark/8._ark*sqrt(3._ark)*y6**4*y9)*y1+((-sqrt(3._ark)*y8/6._ark- &
          sqrt(3._ark)*y7/6._ark+sqrt(3._ark)*y9/24._ark)*y5**4+(y8-y7)*y6*y5**3+(-sqrt(3._ark)*y8/ &
          2._ark-sqrt(3._ark)*y9/4._ark-sqrt(3._ark)*y7/2._ark)*y6**2*y5**2+3._ark/ &
          8._ark*sqrt(3._ark)*y6**4*y9)*y2
      dF(188) = s1+((-sqrt(3._ark)*y8/6._ark-sqrt(3._ark)*y9/24._ark+sqrt(3._ark)*y7/ &
          6._ark)*y5**4+(y7+y8)*y6*y5**3+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark+ &
          sqrt(3._ark)*y9/4._ark)*y6**2*y5**2-3._ark/8._ark*sqrt(3._ark)*y6**4*y9)*y3+ &
          ((sqrt(3._ark)*y8/6._ark-sqrt(3._ark)*y9/24._ark-sqrt(3._ark)*y7/6._ark)*y5**4+(-y8- &
          y7)*y6*y5**3+(sqrt(3._ark)*y9/4._ark-sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/ &
          2._ark)*y6**2*y5**2-3._ark/8._ark*sqrt(3._ark)*y6**4*y9)*y4
      dF(189) = (-3._ark*y5*y6**4+y5**5-2._ark*y5**3*y6**2)*y1+(-3._ark*y5*y6**4+y5**5- &
          2._ark*y5**3*y6**2)*y2+(-3._ark*y5*y6**4+y5**5-2._ark*y5**3*y6**2)*y3+(-3._ark*y5*y6**4+ &
          y5**5-2._ark*y5**3*y6**2)*y4
      dF(190) = ((sqrt(3._ark)*y5**2*y9**2/2._ark-sqrt(3._ark)*y6**2*y9**2/6._ark)*y2+(- &
          y5*y6*y7**2+sqrt(3._ark)*y6**2*y7**2/3._ark)*y3+(y5*y6*y8**2+sqrt(3._ark)*y6**2*y8**2/ &
          3._ark)*y4)*y1+((y5*y6*y8**2+sqrt(3._ark)*y6**2*y8**2/3._ark)*y3+(-y5*y6*y7**2+ &
          sqrt(3._ark)*y6**2*y7**2/3._ark)*y4)*y2+(sqrt(3._ark)*y5**2*y9**2/2._ark- &
          sqrt(3._ark)*y6**2*y9**2/6._ark)*y4*y3
      dF(191) = ((sqrt(3._ark)*y5**2*y7*y8/2._ark-sqrt(3._ark)*y6**2*y7*y8/6._ark)*y2+ &
          (sqrt(3._ark)*y6**2*y8*y9/3._ark-y5*y6*y8*y9)*y3+(sqrt(3._ark)*y6**2*y7*y9/3._ark+ &
          y5*y6*y7*y9)*y4)*y1+((-sqrt(3._ark)*y6**2*y7*y9/3._ark-y5*y6*y7*y9)*y3+(- &
          sqrt(3._ark)*y6**2*y8*y9/3._ark+y5*y6*y8*y9)*y4)*y2+(sqrt(3._ark)*y6**2*y7*y8/6._ark- &
          sqrt(3._ark)*y5**2*y7*y8/2._ark)*y4*y3
      dF(192) = ((y8*y7*y5**2+y8*y7*y6**2)*y2+(y5**2*y8*y9+y9*y8*y6**2)*y3+ &
          (y5**2*y7*y9+y9*y7*y6**2)*y4)*y1+((-y9*y7*y6**2-y5**2*y7*y9)*y3+(-y5**2*y8*y9- &
          y9*y8*y6**2)*y4)*y2+(-y8*y7*y5**2-y8*y7*y6**2)*y4*y3
      dF(193) = ((2._ark/9._ark*sqrt(3._ark)*y6**4+2._ark*sqrt(3._ark)*y5**2*y6**2)*y2+(5._ark/ &
          3._ark*y5*y6**3+7._ark/18._ark*sqrt(3._ark)*y6**4-y5**3*y6+sqrt(3._ark)*y5**4/2._ark)*y3+ &
          (7._ark/18._ark*sqrt(3._ark)*y6**4+sqrt(3._ark)*y5**4/2._ark+y5**3*y6-5._ark/ &
          3._ark*y5*y6**3)*y4)*y1+((7._ark/18._ark*sqrt(3._ark)*y6**4+sqrt(3._ark)*y5**4/2._ark+ &
          y5**3*y6-5._ark/3._ark*y5*y6**3)*y3+(5._ark/3._ark*y5*y6**3+7._ark/18._ark*sqrt(3._ark)*y6**4- &
          y5**3*y6+sqrt(3._ark)*y5**4/2._ark)*y4)*y2+(2._ark/9._ark*sqrt(3._ark)*y6**4+ &
          2._ark*sqrt(3._ark)*y5**2*y6**2)*y4*y3
      dF(194) = (y4*y8**3+y2*y9**3+y3*y7**3)*y1**2+(y2**2*y9**3+y4**2*y8**3+ &
          y3**2*y7**3)*y1+(-y4*y7**3-y3*y8**3)*y2**2+(-y3**2*y8**3-y4**2*y7**3)*y2- &
          y3**2*y4*y9**3-y3*y4**2*y9**3
      s1 = (((sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5+(-y7*y9/2._ark+y8*y9/ &
          2._ark)*y6)*y2+(-sqrt(3._ark)*y5*y7*y9/2._ark+(y8+y9/2._ark)*y7*y6)*y3+(- &
          sqrt(3._ark)*y5*y8*y9/2._ark+(-y8*y9/2._ark-y7*y8)*y6)*y4)*y1**2+(((-sqrt(3._ark)*y8*y9/ &
          2._ark-sqrt(3._ark)*y7*y9/2._ark)*y5+(y7*y9/2._ark-y8*y9/2._ark)*y6)*y2**2+ &
          (sqrt(3._ark)*y5*y7*y9/2._ark+(-y9/2._ark-y8)*y7*y6)*y3**2+(sqrt(3._ark)*y5*y8*y9/2._ark+ &
          (y7*y8+y8*y9/2._ark)*y6)*y4**2)*y1+((sqrt(3._ark)*y5*y8*y9/2._ark+(y8*y9/2._ark- &
          y7*y8)*y6)*y3+(sqrt(3._ark)*y5*y7*y9/2._ark+(y8-y9/2._ark)*y7*y6)*y4)*y2**2
      dF(195) = s1+((-sqrt(3._ark)*y5*y8*y9/2._ark+(y7*y8-y8*y9/2._ark)*y6)*y3**2+(- &
          sqrt(3._ark)*y5*y7*y9/2._ark+(y9/2._ark-y8)*y7*y6)*y4**2)*y2+((-sqrt(3._ark)*y7*y9/2._ark+ &
          sqrt(3._ark)*y8*y9/2._ark)*y5+(y8*y9/2._ark+y7*y9/2._ark)*y6)*y4*y3**2+((- &
          sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5+(-y8*y9/2._ark-y7*y9/ &
          2._ark)*y6)*y4**2*y3
      dF(196) = (((-y8-y7)*y5**2+(-y8-y7)*y6**2)*y2+((-y9-y8)*y5**2+(-y9- &
          y8)*y6**2)*y3+((-y9-y7)*y5**2+(-y9-y7)*y6**2)*y4)*y1**2+(((y7+y8)*y5**2+(y7+ &
          y8)*y6**2)*y2**2+((y9+y8)*y5**2+(y9+y8)*y6**2)*y3**2+((y7+y9)*y5**2+(y7+ &
          y9)*y6**2)*y4**2)*y1+(((y7-y9)*y5**2+(y7-y9)*y6**2)*y3+((y8-y9)*y5**2+(y8- &
          y9)*y6**2)*y4)*y2**2+(((y9-y7)*y5**2+(y9-y7)*y6**2)*y3**2+((-y8+y9)*y5**2+(-y8+ &
          y9)*y6**2)*y4**2)*y2+((y8-y7)*y5**2+(y8-y7)*y6**2)*y4*y3**2+((y7-y8)*y5**2+(y7- &
          y8)*y6**2)*y4**2*y3
      s1 = (((y8**2+y7**2)*y5+(sqrt(3._ark)*y7**2-sqrt(3._ark)*y8**2)*y6)*y2+((y9**2- &
          2._ark*y8**2)*y5+sqrt(3._ark)*y6*y9**2)*y3+((-2._ark*y7**2+y9**2)*y5- &
          sqrt(3._ark)*y6*y9**2)*y4)*y1**2+(((y8**2+y7**2)*y5+(sqrt(3._ark)*y7**2- &
          sqrt(3._ark)*y8**2)*y6)*y2**2+((y9**2-2._ark*y8**2)*y5+sqrt(3._ark)*y6*y9**2)*y3**2+ &
          ((-2._ark*y7**2+y9**2)*y5-sqrt(3._ark)*y6*y9**2)*y4**2)*y1+(((-2._ark*y7**2+y9**2)*y5- &
          sqrt(3._ark)*y6*y9**2)*y3+((y9**2-2._ark*y8**2)*y5+sqrt(3._ark)*y6*y9**2)*y4)*y2**2
      dF(197) = s1+(((-2._ark*y7**2+y9**2)*y5-sqrt(3._ark)*y6*y9**2)*y3**2+((y9**2- &
          2._ark*y8**2)*y5+sqrt(3._ark)*y6*y9**2)*y4**2)*y2+((y8**2+y7**2)*y5+ &
          (sqrt(3._ark)*y7**2-sqrt(3._ark)*y8**2)*y6)*y4*y3**2+((y8**2+y7**2)*y5+ &
          (sqrt(3._ark)*y7**2-sqrt(3._ark)*y8**2)*y6)*y4**2*y3
      dF(198) = (-2._ark*y2*y5*y9**2+(-sqrt(3._ark)*y6*y7**2+y5*y7**2)*y3+(y5*y8**2+ &
          sqrt(3._ark)*y6*y8**2)*y4)*y1**2+(-2._ark*y2**2*y5*y9**2+(-sqrt(3._ark)*y6*y7**2+ &
          y5*y7**2)*y3**2+(y5*y8**2+sqrt(3._ark)*y6*y8**2)*y4**2)*y1+((y5*y8**2+ &
          sqrt(3._ark)*y6*y8**2)*y3+(-sqrt(3._ark)*y6*y7**2+y5*y7**2)*y4)*y2**2+((y5*y8**2+ &
          sqrt(3._ark)*y6*y8**2)*y3**2+(-sqrt(3._ark)*y6*y7**2+y5*y7**2)*y4**2)*y2- &
          2._ark*y3*y4**2*y5*y9**2-2._ark*y3**2*y4*y5*y9**2
      s1 = (((y8*y9/2._ark+y7*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/ &
          2._ark)*y6)*y2+((y9/2._ark-y8)*y7*y5+sqrt(3._ark)*y6*y7*y9/2._ark)*y3+((y8*y9/2._ark- &
          y7*y8)*y5-sqrt(3._ark)*y6*y8*y9/2._ark)*y4)*y1**2+(((-y8*y9/2._ark-y7*y9/2._ark)*y5+(- &
          sqrt(3._ark)*y7*y9/2._ark+sqrt(3._ark)*y8*y9/2._ark)*y6)*y2**2+((y8-y9/2._ark)*y7*y5- &
          sqrt(3._ark)*y6*y7*y9/2._ark)*y3**2+((y7*y8-y8*y9/2._ark)*y5+sqrt(3._ark)*y6*y8*y9/ &
          2._ark)*y4**2)*y1+(((-y8*y9/2._ark-y7*y8)*y5+sqrt(3._ark)*y6*y8*y9/2._ark)*y3+((-y9/ &
          2._ark-y8)*y7*y5-sqrt(3._ark)*y6*y7*y9/2._ark)*y4)*y2**2
      dF(199) = s1+(((y7*y8+y8*y9/2._ark)*y5-sqrt(3._ark)*y6*y8*y9/2._ark)*y3**2+((y8+y9/ &
          2._ark)*y7*y5+sqrt(3._ark)*y6*y7*y9/2._ark)*y4**2)*y2+((-y7*y9/2._ark+y8*y9/2._ark)*y5+(- &
          sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y6)*y4*y3**2+((y7*y9/2._ark-y8*y9/ &
          2._ark)*y5+(sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y4**2*y3
      dF(200) = ((y9*y6**2+y9*y5**2)*y2+(y7*y6**2+y5**2*y7)*y3+(y8*y6**2+ &
          y5**2*y8)*y4)*y1**2+((y9*y6**2+y9*y5**2)*y2**2+(y7*y6**2+y5**2*y7)*y3**2+ &
          (y8*y6**2+y5**2*y8)*y4**2)*y1+((-y5**2*y8-y8*y6**2)*y3+(-y5**2*y7- &
          y7*y6**2)*y4)*y2**2+((-y5**2*y8-y8*y6**2)*y3**2+(-y5**2*y7-y7*y6**2)*y4**2)*y2+ &
          (-y9*y6**2-y9*y5**2)*y4*y3**2+(-y9*y6**2-y9*y5**2)*y4**2*y3
      dF(201) = ((y5**3-3._ark*y5*y6**2)*y2+(y5**3-3._ark*y5*y6**2)*y3+(y5**3- &
          3._ark*y5*y6**2)*y4)*y1**2+((y5**3-3._ark*y5*y6**2)*y2**2+(y5**3- &
          3._ark*y5*y6**2)*y3**2+(y5**3-3._ark*y5*y6**2)*y4**2)*y1+((y5**3-3._ark*y5*y6**2)*y3+ &
          (y5**3-3._ark*y5*y6**2)*y4)*y2**2+((y5**3-3._ark*y5*y6**2)*y3**2+(y5**3- &
          3._ark*y5*y6**2)*y4**2)*y2+(y5**3-3._ark*y5*y6**2)*y4*y3**2+(y5**3- &
          3._ark*y5*y6**2)*y4**2*y3
      dF(202) = ((y8**2+y7**2)*y2+(y9**2+y8**2)*y3+(y7**2+y9**2)*y4)*y1**3+((y8**2+ &
          y7**2)*y2**3+(y9**2+y8**2)*y3**3+(y7**2+y9**2)*y4**3)*y1+((y7**2+y9**2)*y3+ &
          (y9**2+y8**2)*y4)*y2**3+((y7**2+y9**2)*y3**3+(y9**2+y8**2)*y4**3)*y2+(y8**2+ &
          y7**2)*y4*y3**3+(y8**2+y7**2)*y4**3*y3
      dF(203) = ((-y8*y9-y7*y9)*y2+(-y9-y8)*y7*y3+(-y7*y8-y8*y9)*y4)*y1**3+((y8*y9+ &
          y7*y9)*y2**3+(y9+y8)*y7*y3**3+(y8*y9+y7*y8)*y4**3)*y1+((-y7*y8+y8*y9)*y3+(-y8+ &
          y9)*y7*y4)*y2**3+((-y8*y9+y7*y8)*y3**3+(y8-y9)*y7*y4**3)*y2+(-y8*y9+ &
          y7*y9)*y4*y3**3+(y8*y9-y7*y9)*y4**3*y3
      dF(204) = (y3*y7**2+y2*y9**2+y4*y8**2)*y1**3+(y3**3*y7**2+y4**3*y8**2+ &
          y2**3*y9**2)*y1+(y4*y7**2+y3*y8**2)*y2**3+(y3**3*y8**2+y4**3*y7**2)*y2+ &
          y3**3*y4*y9**2+y3*y4**3*y9**2
      dF(205) = (y4*y7*y9+y3*y8*y9+y2*y7*y8)*y1**3+(y3**3*y8*y9+y4**3*y7*y9+ &
          y2**3*y7*y8)*y1+(-y4*y8*y9-y3*y7*y9)*y2**3+(-y3**3*y7*y9-y4**3*y8*y9)*y2- &
          y3**3*y4*y7*y8-y3*y4**3*y7*y8
      s1 = (((sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y7/2._ark-y8/2._ark)*y6)*y2+(- &
          sqrt(3._ark)*y5*y8/2._ark+(y9+y8/2._ark)*y6)*y3+(-sqrt(3._ark)*y5*y7/2._ark+(-y9-y7/ &
          2._ark)*y6)*y4)*y1**3+(((-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y5+(y8/2._ark-y7/ &
          2._ark)*y6)*y2**3+(sqrt(3._ark)*y5*y8/2._ark+(-y8/2._ark-y9)*y6)*y3**3+ &
          (sqrt(3._ark)*y5*y7/2._ark+(y9+y7/2._ark)*y6)*y4**3)*y1+((sqrt(3._ark)*y5*y7/2._ark+(y7/ &
          2._ark-y9)*y6)*y3+(sqrt(3._ark)*y5*y8/2._ark+(-y8/2._ark+y9)*y6)*y4)*y2**3
      dF(206) = s1+((-sqrt(3._ark)*y5*y7/2._ark+(-y7/2._ark+y9)*y6)*y3**3+(- &
          sqrt(3._ark)*y5*y8/2._ark+(-y9+y8/2._ark)*y6)*y4**3)*y2+((sqrt(3._ark)*y7/2._ark- &
          sqrt(3._ark)*y8/2._ark)*y5+(y7/2._ark+y8/2._ark)*y6)*y4*y3**3+((-sqrt(3._ark)*y7/2._ark+ &
          sqrt(3._ark)*y8/2._ark)*y5+(-y8/2._ark-y7/2._ark)*y6)*y4**3*y3
      s1 = (((y7/2._ark+y8/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y6)*y2+((- &
          y9+y8/2._ark)*y5+sqrt(3._ark)*y6*y8/2._ark)*y3+((y7/2._ark-y9)*y5-sqrt(3._ark)*y6*y7/ &
          2._ark)*y4)*y1**3+(((-y8/2._ark-y7/2._ark)*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
          2._ark)*y6)*y2**3+((-y8/2._ark+y9)*y5-sqrt(3._ark)*y6*y8/2._ark)*y3**3+((-y7/2._ark+ &
          y9)*y5+sqrt(3._ark)*y6*y7/2._ark)*y4**3)*y1+(((-y9-y7/2._ark)*y5+sqrt(3._ark)*y6*y7/ &
          2._ark)*y3+((-y8/2._ark-y9)*y5-sqrt(3._ark)*y6*y8/2._ark)*y4)*y2**3
      dF(207) = s1+(((y9+y7/2._ark)*y5-sqrt(3._ark)*y6*y7/2._ark)*y3**3+((y9+y8/2._ark)*y5+ &
          sqrt(3._ark)*y6*y8/2._ark)*y4**3)*y2+((y7/2._ark-y8/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark- &
          sqrt(3._ark)*y8/2._ark)*y6)*y4*y3**3+((y8/2._ark-y7/2._ark)*y5+(sqrt(3._ark)*y7/2._ark+ &
          sqrt(3._ark)*y8/2._ark)*y6)*y4**3*y3
      dF(208) = ((-sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/2._ark)*y2+(-y5*y6+ &
          sqrt(3._ark)*y6**2/3._ark)*y3+(y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4)*y1**3+((- &
          sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/2._ark)*y2**3+(-y5*y6+sqrt(3._ark)*y6**2/ &
          3._ark)*y3**3+(y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4**3)*y1+((y5*y6+sqrt(3._ark)*y6**2/ &
          3._ark)*y3+(-y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4)*y2**3+((y5*y6+sqrt(3._ark)*y6**2/ &
          3._ark)*y3**3+(-y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4**3)*y2+(-sqrt(3._ark)*y6**2/6._ark+ &
          sqrt(3._ark)*y5**2/2._ark)*y4*y3**3+(-sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/ &
          2._ark)*y4**3*y3
      dF(209) = ((y6**2+y5**2)*y2+(y6**2+y5**2)*y3+(y6**2+y5**2)*y4)*y1**3+((y6**2+ &
          y5**2)*y2**3+(y6**2+y5**2)*y3**3+(y6**2+y5**2)*y4**3)*y1+((y6**2+y5**2)*y3+ &
          (y6**2+y5**2)*y4)*y2**3+((y6**2+y5**2)*y3**3+(y6**2+y5**2)*y4**3)*y2+(y6**2+ &
          y5**2)*y4*y3**3+(y6**2+y5**2)*y4**3*y3
      dF(210) = (y3*y7+y4*y8+y2*y9)*y1**4+(y2**4*y9+y4**4*y8+y3**4*y7)*y1+(-y3*y8- &
          y4*y7)*y2**4+(-y4**4*y7-y3**4*y8)*y2-y3**4*y4*y9-y3*y4**4*y9
      dF(211) = ((-y8-y7)*y2+(-y9-y8)*y3+(-y9-y7)*y4)*y1**4+((y7+y8)*y2**4+(y9+ &
          y8)*y3**4+(y7+y9)*y4**4)*y1+((y7-y9)*y3+(y8-y9)*y4)*y2**4+((y9-y7)*y3**4+(-y8+ &
          y9)*y4**4)*y2+(y8-y7)*y4*y3**4+(y7-y8)*y4**4*y3
      dF(212) = (-2._ark/3._ark*sqrt(3._ark)*y2*y5+(-y6+sqrt(3._ark)*y5/3._ark)*y3+(y6+ &
          sqrt(3._ark)*y5/3._ark)*y4)*y1**4+(-2._ark/3._ark*sqrt(3._ark)*y2**4*y5+(-y6+ &
          sqrt(3._ark)*y5/3._ark)*y3**4+(y6+sqrt(3._ark)*y5/3._ark)*y4**4)*y1+((y6+sqrt(3._ark)*y5/ &
          3._ark)*y3+(-y6+sqrt(3._ark)*y5/3._ark)*y4)*y2**4+((y6+sqrt(3._ark)*y5/3._ark)*y3**4+(-y6+ &
          sqrt(3._ark)*y5/3._ark)*y4**4)*y2-2._ark/3._ark*sqrt(3._ark)*y3*y4**4*y5-2._ark/ &
          3._ark*sqrt(3._ark)*y3**4*y4*y5
      dF(213) = ((y8**4+y7**4)*y2+(y8**4+y9**4)*y3+(y9**4+y7**4)*y4)*y1+((y9**4+ &
          y7**4)*y3+(y8**4+y9**4)*y4)*y2+(y8**4+y7**4)*y4*y3
      dF(214) = ((y7**3*y8+y7*y8**3)*y2+(y8**3*y9+y8*y9**3)*y3+(y7*y9**3+ &
          y7**3*y9)*y4)*y1+((-y7*y9**3-y7**3*y9)*y3+(-y8**3*y9-y8*y9**3)*y4)*y2+(- &
          y7**3*y8-y7*y8**3)*y4*y3
      dF(215) = (y4*y7**2*y9**2+y2*y7**2*y8**2+y3*y8**2*y9**2)*y1+(y3*y7**2*y9**2+ &
          y4*y8**2*y9**2)*y2+y3*y4*y7**2*y8**2
      dF(216) = (y3*y7**2*y8*y9+y4*y7*y8**2*y9+y2*y7*y8*y9**2)*y1+(-y4*y7**2*y8*y9- &
          y3*y7*y8**2*y9)*y2-y3*y4*y7*y8*y9**2
      dF(217) = (y4*y8**4+y3*y7**4+y2*y9**4)*y1+(y3*y8**4+y4*y7**4)*y2+y3*y4*y9**4
      dF(218) = (8._ark/3._ark*sqrt(3._ark)*y2*y5*y6**2*y9+(-sqrt(3._ark)*y5**3*y7+ &
          y5**2*y6*y7+5._ark/3._ark*sqrt(3._ark)*y5*y6**2*y7+y6**3*y7)*y3+(-y5**2*y6*y8+5._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**2*y8-sqrt(3._ark)*y5**3*y8-y6**3*y8)*y4)*y1+((y6**3*y8+ &
          y5**2*y6*y8-5._ark/3._ark*sqrt(3._ark)*y5*y6**2*y8+sqrt(3._ark)*y5**3*y8)*y3+(- &
          y5**2*y6*y7-y6**3*y7+sqrt(3._ark)*y5**3*y7-5._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**2*y7)*y4)*y2-8._ark/3._ark*sqrt(3._ark)*y3*y4*y5*y6**2*y9
      dF(219) = ((y5**2*y9**2+y6**2*y9**2)*y2+(y6**2*y7**2+y5**2*y7**2)*y3+ &
          (y6**2*y8**2+y5**2*y8**2)*y4)*y1+((y6**2*y8**2+y5**2*y8**2)*y3+(y6**2*y7**2+ &
          y5**2*y7**2)*y4)*y2+(y5**2*y9**2+y6**2*y9**2)*y4*y3
      dF(220) = ((4._ark/3._ark*y6**4+4._ark*y5**2*y6**2)*y2+(4._ark/3._ark*sqrt(3._ark)*y5*y6**3+ &
          y5**2*y6**2+5._ark/6._ark*y6**4+3._ark/2._ark*y5**4)*y3+(3._ark/2._ark*y5**4-4._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**3+y5**2*y6**2+5._ark/6._ark*y6**4)*y4)*y1+((3._ark/2._ark*y5**4- &
          4._ark/3._ark*sqrt(3._ark)*y5*y6**3+y5**2*y6**2+5._ark/6._ark*y6**4)*y3+(4._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**3+y5**2*y6**2+5._ark/6._ark*y6**4+3._ark/2._ark*y5**4)*y4)*y2+ &
          (4._ark/3._ark*y6**4+4._ark*y5**2*y6**2)*y4*y3
      dF(221) = ((y4*y7*y8*y9+y3*y7*y8*y9)*y2+y3*y4*y7*y8*y9)*y1+y2*y3*y4*y7*y8*y9
      dF(222) = (((((-y9/2._ark-y8)*y7+y8*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark- &
          sqrt(3._ark)*y7*y9/2._ark)*y6)*y3+(((y9/2._ark-y8)*y7-y8*y9/2._ark)*y5+ &
          (sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y4)*y2+(((y8+y9/2._ark)*y7+ &
          y8*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y4*y3)*y1+ &
          (((y8-y9/2._ark)*y7-y8*y9/2._ark)*y5+(-sqrt(3._ark)*y7*y9/2._ark+sqrt(3._ark)*y8*y9/ &
          2._ark)*y6)*y4*y3*y2
      dF(223) = ((((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y3+((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y4)*y2+((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y4*y3)*y1+((y8**2+y7**2-2._ark*y9**2)*y5+(sqrt(3._ark)*y8**2- &
          sqrt(3._ark)*y7**2)*y6)*y4*y3*y2
      dF(224) = ((((y9+y7-y8)*y5**2+(y9+y7-y8)*y6**2)*y3+((y8+y9-y7)*y5**2+(y8+y9- &
          y7)*y6**2)*y4)*y2+((y7-y9+y8)*y5**2+(y7-y9+y8)*y6**2)*y4*y3)*y1+((-y7-y8- &
          y9)*y5**2+(-y7-y8-y9)*y6**2)*y4*y3*y2
      dF(225) = (((y5**3-3._ark*y5*y6**2)*y3+(y5**3-3._ark*y5*y6**2)*y4)*y2+(y5**3- &
          3._ark*y5*y6**2)*y4*y3)*y1+(y5**3-3._ark*y5*y6**2)*y4*y3*y2
      dF(226) = (((-sqrt(3._ark)*y6*y8-y5*y8)*y3+(sqrt(3._ark)*y6*y7-y5*y7)*y4)*y2+ &
          2._ark*y3*y4*y5*y9)*y1**2+(((-sqrt(3._ark)*y6*y7+y5*y7)*y3+(sqrt(3._ark)*y6*y8+ &
          y5*y8)*y4)*y2**2+(-2._ark*y3**2*y5*y9-2._ark*y4**2*y5*y9)*y2+(sqrt(3._ark)*y6*y8+ &
          y5*y8)*y4*y3**2+(-sqrt(3._ark)*y6*y7+y5*y7)*y4**2*y3)*y1+2._ark*y2**2*y3*y4*y5*y9+ &
          ((sqrt(3._ark)*y6*y7-y5*y7)*y4*y3**2+(-sqrt(3._ark)*y6*y8-y5*y8)*y4**2*y3)*y2
      dF(227) = (((y6**2+y5**2)*y3+(y6**2+y5**2)*y4)*y2+(y6**2+y5**2)*y4*y3)*y1**2+ &
          (((y6**2+y5**2)*y3+(y6**2+y5**2)*y4)*y2**2+((y6**2+y5**2)*y3**2+(y6**2+ &
          y5**2)*y4**2)*y2+(y6**2+y5**2)*y4*y3**2+(y6**2+y5**2)*y4**2*y3)*y1+(y6**2+ &
          y5**2)*y4*y3*y2**2+((y6**2+y5**2)*y4*y3**2+(y6**2+y5**2)*y4**2*y3)*y2
      dF(228) = ((sqrt(3._ark)*y5**3-sqrt(3._ark)*y5*y6**2/3._ark)*y2+(y6**3-4._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**2+y5**2*y6)*y3+(-y6**3-y5**2*y6-4._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**2)*y4)*y1**2+((sqrt(3._ark)*y5**3-sqrt(3._ark)*y5*y6**2/ &
          3._ark)*y2**2+(y6**3-4._ark/3._ark*sqrt(3._ark)*y5*y6**2+y5**2*y6)*y3**2+(-y6**3- &
          y5**2*y6-4._ark/3._ark*sqrt(3._ark)*y5*y6**2)*y4**2)*y1+((-y6**3-y5**2*y6-4._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**2)*y3+(y6**3-4._ark/3._ark*sqrt(3._ark)*y5*y6**2+ &
          y5**2*y6)*y4)*y2**2+((-y6**3-y5**2*y6-4._ark/3._ark*sqrt(3._ark)*y5*y6**2)*y3**2+ &
          (y6**3-4._ark/3._ark*sqrt(3._ark)*y5*y6**2+y5**2*y6)*y4**2)*y2+(sqrt(3._ark)*y5**3- &
          sqrt(3._ark)*y5*y6**2/3._ark)*y4*y3**2+(sqrt(3._ark)*y5**3-sqrt(3._ark)*y5*y6**2/ &
          3._ark)*y4**2*y3
      dF(229) = ((-y8**2*y9+y7**2*y9)*y6*y2+((-sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y9**2/ &
          2._ark)*y7*y5+(y9**2/2._ark-y8**2/2._ark)*y7*y6)*y3+((-sqrt(3._ark)*y7**2*y8/2._ark+ &
          sqrt(3._ark)*y8*y9**2/2._ark)*y5+(-y8*y9**2/2._ark+y7**2*y8/2._ark)*y6)*y4)*y1+(((- &
          sqrt(3._ark)*y8*y9**2/2._ark+sqrt(3._ark)*y7**2*y8/2._ark)*y5+(-y7**2*y8/2._ark+y8*y9**2/ &
          2._ark)*y6)*y3+((sqrt(3._ark)*y8**2/2._ark-sqrt(3._ark)*y9**2/2._ark)*y7*y5+(y8**2/2._ark- &
          y9**2/2._ark)*y7*y6)*y4)*y2+(-y7**2*y9+y8**2*y9)*y6*y4*y3
      dF(230) = (y2*y5*y9**3+(sqrt(3._ark)*y6*y7**3/2._ark-y5*y7**3/2._ark)*y3+(- &
          sqrt(3._ark)*y6*y8**3/2._ark-y5*y8**3/2._ark)*y4)*y1+((y5*y8**3/2._ark+ &
          sqrt(3._ark)*y6*y8**3/2._ark)*y3+(-sqrt(3._ark)*y6*y7**3/2._ark+y5*y7**3/2._ark)*y4)*y2- &
          y3*y4*y5*y9**3
      dF(231) = (y2*y5*y7*y8*y9+(sqrt(3._ark)*y6*y7*y8*y9/2._ark-y5*y7*y8*y9/2._ark)*y3+(- &
          y5*y7*y8*y9/2._ark-sqrt(3._ark)*y6*y7*y8*y9/2._ark)*y4)*y1+((-y5*y7*y8*y9/2._ark- &
          sqrt(3._ark)*y6*y7*y8*y9/2._ark)*y3+(sqrt(3._ark)*y6*y7*y8*y9/2._ark-y5*y7*y8*y9/ &
          2._ark)*y4)*y2+y3*y4*y5*y7*y8*y9
      dF(232) = ((y7**2*y9+y8**2*y9)*y5*y2+((-y8**2/2._ark-y9**2/2._ark)*y7*y5+ &
          (sqrt(3._ark)*y9**2/2._ark+sqrt(3._ark)*y8**2/2._ark)*y7*y6)*y3+((-y8*y9**2/2._ark- &
          y7**2*y8/2._ark)*y5+(-sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7**2*y8/ &
          2._ark)*y6)*y4)*y1+(((y8*y9**2/2._ark+y7**2*y8/2._ark)*y5+(sqrt(3._ark)*y8*y9**2/2._ark+ &
          sqrt(3._ark)*y7**2*y8/2._ark)*y6)*y3+((y8**2/2._ark+y9**2/2._ark)*y7*y5+(- &
          sqrt(3._ark)*y9**2/2._ark-sqrt(3._ark)*y8**2/2._ark)*y7*y6)*y4)*y2+(-y8**2*y9- &
          y7**2*y9)*y5*y4*y3
      dF(233) = (((-sqrt(3._ark)*y7**2/2._ark-sqrt(3._ark)*y8**2/2._ark)*y5**2+(y7**2- &
          y8**2)*y6*y5+(-sqrt(3._ark)*y8**2/6._ark-sqrt(3._ark)*y7**2/6._ark)*y6**2)*y2+(- &
          sqrt(3._ark)*y5**2*y9**2/2._ark+y5*y6*y9**2+(-sqrt(3._ark)*y9**2/6._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y8**2)*y6**2)*y3+(-sqrt(3._ark)*y5**2*y9**2/2._ark-y5*y6*y9**2+(- &
          sqrt(3._ark)*y9**2/6._ark-2._ark/3._ark*sqrt(3._ark)*y7**2)*y6**2)*y4)*y1+((- &
          sqrt(3._ark)*y5**2*y9**2/2._ark-y5*y6*y9**2+(-sqrt(3._ark)*y9**2/6._ark-2._ark/ &
          3._ark*sqrt(3._ark)*y7**2)*y6**2)*y3+(-sqrt(3._ark)*y5**2*y9**2/2._ark+y5*y6*y9**2+(- &
          sqrt(3._ark)*y9**2/6._ark-2._ark/3._ark*sqrt(3._ark)*y8**2)*y6**2)*y4)*y2+((- &
          sqrt(3._ark)*y7**2/2._ark-sqrt(3._ark)*y8**2/2._ark)*y5**2+(y7**2-y8**2)*y6*y5+(- &
          sqrt(3._ark)*y8**2/6._ark-sqrt(3._ark)*y7**2/6._ark)*y6**2)*y4*y3
      dF(234) = ((-3._ark*y5*y6**2*y9+y5**3*y9)*y2+(y5**3*y7-3._ark*y5*y6**2*y7)*y3+ &
          (y5**3*y8-3._ark*y5*y6**2*y8)*y4)*y1+((-y5**3*y8+3._ark*y5*y6**2*y8)*y3+(-y5**3*y7+ &
          3._ark*y5*y6**2*y7)*y4)*y2+(-y5**3*y9+3._ark*y5*y6**2*y9)*y4*y3
      dF(235) = ((-5._ark/3._ark*y6**4-6._ark*y5**2*y6**2+y5**4)*y2+(-8._ark/ &
          3._ark*sqrt(3._ark)*y5*y6**3-2._ark/3._ark*y6**4-2._ark*y5**4)*y3+(-2._ark*y5**4-2._ark/ &
          3._ark*y6**4+8._ark/3._ark*sqrt(3._ark)*y5*y6**3)*y4)*y1+((-2._ark*y5**4-2._ark/3._ark*y6**4+ &
          8._ark/3._ark*sqrt(3._ark)*y5*y6**3)*y3+(-8._ark/3._ark*sqrt(3._ark)*y5*y6**3-2._ark/ &
          3._ark*y6**4-2._ark*y5**4)*y4)*y2+(-5._ark/3._ark*y6**4-6._ark*y5**2*y6**2+y5**4)*y4*y3
      dF(236) = (((y7**3+y9**3-y8**3)*y3+(y8**3+y9**3-y7**3)*y4)*y2+(-y9**3+y8**3+ &
          y7**3)*y4*y3)*y1+(-y8**3-y7**3-y9**3)*y4*y3*y2
      dF(237) = ((((-y8+y9)*y7**2+(y9**2+y8**2)*y7-y8*y9**2+y8**2*y9)*y3+((y9+ &
          y8)*y7**2+(-y9**2-y8**2)*y7+y8**2*y9+y8*y9**2)*y4)*y2+((y8-y9)*y7**2+(y9**2+ &
          y8**2)*y7-y8**2*y9+y8*y9**2)*y4*y3)*y1+((-y9-y8)*y7**2+(-y9**2-y8**2)*y7- &
          y8*y9**2-y8**2*y9)*y4*y3*y2
      dF(238) = (((-sqrt(3._ark)*y5**2*y9/2._ark+(y7+y8)*y6*y5+(sqrt(3._ark)*y9/6._ark+ &
          sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y7/3._ark)*y6**2)*y3+(-sqrt(3._ark)*y5**2*y9/2._ark+(- &
          y8-y7)*y6*y5+(sqrt(3._ark)*y9/6._ark-sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/ &
          3._ark)*y6**2)*y4)*y2+(sqrt(3._ark)*y5**2*y9/2._ark+(y7-y8)*y6*y5+(-sqrt(3._ark)*y7/ &
          3._ark-sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y9/6._ark)*y6**2)*y4*y3)*y1+ &
          (sqrt(3._ark)*y5**2*y9/2._ark+(y8-y7)*y6*y5+(sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/3._ark- &
          sqrt(3._ark)*y9/6._ark)*y6**2)*y4*y3*y2
      dF(239) = (y8**2+y7**2+y9**2)*y4*y3*y2*y1
      dF(240) = (y6**2+y5**2)*y4*y3*y2*y1
      dF(241) = (y9**4+y8**4+y7**4)*y1**2+(y9**4+y8**4+y7**4)*y2**2+(y9**4+y8**4+ &
          y7**4)*y3**2+(y9**4+y8**4+y7**4)*y4**2
      dF(242) = (y7**2*y8*y9+(y8*y9**2+y8**2*y9)*y7)*y1**2+(-y7**2*y8*y9+(y8*y9**2- &
          y8**2*y9)*y7)*y2**2+(y7**2*y8*y9+(-y8**2*y9-y8*y9**2)*y7)*y3**2+(-y7**2*y8*y9+(- &
          y8*y9**2+y8**2*y9)*y7)*y4**2
      dF(243) = ((y9**2+y8**2)*y7**2+y8**2*y9**2)*y1**2+((y9**2+y8**2)*y7**2+ &
          y8**2*y9**2)*y2**2+((y9**2+y8**2)*y7**2+y8**2*y9**2)*y3**2+((y9**2+y8**2)*y7**2+ &
          y8**2*y9**2)*y4**2
      dF(244) = ((y9+y8)*y7**3+(y9**3+y8**3)*y7+y8**3*y9+y8*y9**3)*y1**2+((y8- &
          y9)*y7**3+(y8**3-y9**3)*y7-y8*y9**3-y8**3*y9)*y2**2+((-y9-y8)*y7**3+(-y9**3- &
          y8**3)*y7+y8**3*y9+y8*y9**3)*y3**2+((-y8+y9)*y7**3+(y9**3-y8**3)*y7-y8*y9**3- &
          y8**3*y9)*y4**2
      s1 = ((sqrt(3._ark)*y7*y9**2/2._ark-sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y7**2*y9/ &
          2._ark+sqrt(3._ark)*y8*y9**2/2._ark)*y5+((y8+y9/2._ark)*y7**2+(-y9**2/2._ark-y8**2)*y7- &
          y8**2*y9/2._ark+y8*y9**2/2._ark)*y6)*y1**2+((-sqrt(3._ark)*y7*y9**2/2._ark- &
          sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7**2*y9/2._ark-sqrt(3._ark)*y8**2*y9/2._ark)*y5+ &
          ((y9/2._ark-y8)*y7**2+(y9**2/2._ark+y8**2)*y7-y8*y9**2/2._ark-y8**2*y9/2._ark)*y6)*y2**2
      dF(245) = s1+((sqrt(3._ark)*y7**2*y9/2._ark+sqrt(3._ark)*y7*y9**2/2._ark+ &
          sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y8*y9**2/2._ark)*y5+((-y9/2._ark-y8)*y7**2+(- &
          y9**2/2._ark-y8**2)*y7+y8**2*y9/2._ark-y8*y9**2/2._ark)*y6)*y3**2+ &
          ((sqrt(3._ark)*y8**2*y9/2._ark+sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7*y9**2/2._ark+ &
          sqrt(3._ark)*y7**2*y9/2._ark)*y5+((y8-y9/2._ark)*y7**2+(y9**2/2._ark+y8**2)*y7+y8**2*y9/ &
          2._ark+y8*y9**2/2._ark)*y6)*y4**2
      dF(246) = ((-sqrt(3._ark)*y7**3/3._ark-sqrt(3._ark)*y8**3/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**3)*y5+(-y8**3+y7**3)*y6)*y1**2+((sqrt(3._ark)*y7**3/3._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y9**3+sqrt(3._ark)*y8**3/3._ark)*y5+(-y7**3+y8**3)*y6)*y2**2+((- &
          2._ark/3._ark*sqrt(3._ark)*y9**3-sqrt(3._ark)*y7**3/3._ark+sqrt(3._ark)*y8**3/3._ark)*y5+ &
          (y8**3+y7**3)*y6)*y3**2+((-2._ark/3._ark*sqrt(3._ark)*y9**3+sqrt(3._ark)*y7**3/3._ark- &
          sqrt(3._ark)*y8**3/3._ark)*y5+(-y7**3-y8**3)*y6)*y4**2
      s1 = (((y8-y9/2._ark)*y7**2+(-y9**2/2._ark+y8**2)*y7-y8*y9**2/2._ark-y8**2*y9/ &
          2._ark)*y5+(-sqrt(3._ark)*y7**2*y9/2._ark+sqrt(3._ark)*y8*y9**2/2._ark+ &
          sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y7*y9**2/2._ark)*y6)*y1**2+(((-y9/2._ark- &
          y8)*y7**2+(y9**2/2._ark-y8**2)*y7+y8*y9**2/2._ark-y8**2*y9/2._ark)*y5+(- &
          sqrt(3._ark)*y7**2*y9/2._ark+sqrt(3._ark)*y7*y9**2/2._ark+sqrt(3._ark)*y8**2*y9/2._ark- &
          sqrt(3._ark)*y8*y9**2/2._ark)*y6)*y2**2
      dF(247) = s1+(((y9/2._ark-y8)*y7**2+(-y9**2/2._ark+y8**2)*y7+y8**2*y9/2._ark+y8*y9**2/ &
          2._ark)*y5+(sqrt(3._ark)*y7**2*y9/2._ark-sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y8*y9**2/ &
          2._ark-sqrt(3._ark)*y7*y9**2/2._ark)*y6)*y3**2+(((y8+y9/2._ark)*y7**2+(y9**2/2._ark- &
          y8**2)*y7+y8**2*y9/2._ark-y8*y9**2/2._ark)*y5+(sqrt(3._ark)*y7*y9**2/2._ark- &
          sqrt(3._ark)*y8**2*y9/2._ark+sqrt(3._ark)*y8*y9**2/2._ark+sqrt(3._ark)*y7**2*y9/ &
          2._ark)*y6)*y4**2
      dF(248) = ((-sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y5**2+(-y8*y9+ &
          y7*y9)*y6*y5+((-sqrt(3._ark)*y9/6._ark-2._ark/3._ark*sqrt(3._ark)*y8)*y7-sqrt(3._ark)*y8*y9/ &
          6._ark)*y6**2)*y1**2+((sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5**2+(y8*y9- &
          y7*y9)*y6*y5+((-2._ark/3._ark*sqrt(3._ark)*y8+sqrt(3._ark)*y9/6._ark)*y7+sqrt(3._ark)*y8*y9/ &
          6._ark)*y6**2)*y2**2+((-sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5**2+(- &
          y8*y9-y7*y9)*y6*y5+((sqrt(3._ark)*y9/6._ark+2._ark/3._ark*sqrt(3._ark)*y8)*y7- &
          sqrt(3._ark)*y8*y9/6._ark)*y6**2)*y3**2+((-sqrt(3._ark)*y7*y9/2._ark+sqrt(3._ark)*y8*y9/ &
          2._ark)*y5**2+(y8*y9+y7*y9)*y6*y5+((-sqrt(3._ark)*y9/6._ark+2._ark/ &
          3._ark*sqrt(3._ark)*y8)*y7+sqrt(3._ark)*y8*y9/6._ark)*y6**2)*y4**2
      dF(249) = (-sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/ &
          6._ark-sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2)*y1**2+(- &
          sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark- &
          sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2)*y2**2+(- &
          sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark- &
          sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2)*y3**2+(- &
          sqrt(3._ark)*y5**2*y9**2/2._ark+(y7**2-y8**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark- &
          sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y6**2)*y4**2
      dF(250) = ((-y9/3._ark-y8/3._ark-y7/3._ark)*y5**3+(y9+y7+y8)*y6**2*y5)*y1**2+((y8/ &
          3._ark-y9/3._ark+y7/3._ark)*y5**3+(-y7+y9-y8)*y6**2*y5)*y2**2+((y9/3._ark+y8/3._ark-y7/ &
          3._ark)*y5**3+(-y8-y9+y7)*y6**2*y5)*y3**2+((y7/3._ark+y9/3._ark-y8/3._ark)*y5**3+(-y7- &
          y9+y8)*y6**2*y5)*y4**2
      dF(251) = ((y8**2+y7**2+y9**2)*y5**2+(y8**2+y7**2+y9**2)*y6**2)*y1**2+((y8**2+ &
          y7**2+y9**2)*y5**2+(y8**2+y7**2+y9**2)*y6**2)*y2**2+((y8**2+y7**2+y9**2)*y5**2+ &
          (y8**2+y7**2+y9**2)*y6**2)*y3**2+((y8**2+y7**2+y9**2)*y5**2+(y8**2+y7**2+ &
          y9**2)*y6**2)*y4**2
      dF(252) = ((-4._ark/9._ark*sqrt(3._ark)*y8-4._ark/9._ark*sqrt(3._ark)*y7+5._ark/ &
          9._ark*sqrt(3._ark)*y9)*y5**3+(y7-y8)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y9+(y7- &
          y8)*y6**3)*y1**2+((4._ark/9._ark*sqrt(3._ark)*y7+5._ark/9._ark*sqrt(3._ark)*y9+4._ark/ &
          9._ark*sqrt(3._ark)*y8)*y5**3+(y8-y7)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y9+(y8- &
          y7)*y6**3)*y2**2+((4._ark/9._ark*sqrt(3._ark)*y8-5._ark/9._ark*sqrt(3._ark)*y9-4._ark/ &
          9._ark*sqrt(3._ark)*y7)*y5**3+(y7+y8)*y6*y5**2-sqrt(3._ark)*y5*y6**2*y9+(y7+ &
          y8)*y6**3)*y3**2+((-5._ark/9._ark*sqrt(3._ark)*y9+4._ark/9._ark*sqrt(3._ark)*y7-4._ark/ &
          9._ark*sqrt(3._ark)*y8)*y5**3+(-y8-y7)*y6*y5**2-sqrt(3._ark)*y5*y6**2*y9+(-y8- &
          y7)*y6**3)*y4**2
      dF(253) = ((sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y2+(- &
          sqrt(3._ark)*y6**2*y7/3._ark+y5*y6*y7)*y3+(-sqrt(3._ark)*y6**2*y8/3._ark- &
          y5*y6*y8)*y4)*y1**2+((sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y2**2+ &
          (-sqrt(3._ark)*y6**2*y7/3._ark+y5*y6*y7)*y3**2+(-sqrt(3._ark)*y6**2*y8/3._ark- &
          y5*y6*y8)*y4**2)*y1+((sqrt(3._ark)*y6**2*y8/3._ark+y5*y6*y8)*y3+(-y5*y6*y7+ &
          sqrt(3._ark)*y6**2*y7/3._ark)*y4)*y2**2+((sqrt(3._ark)*y6**2*y8/3._ark+y5*y6*y8)*y3**2+ &
          (-y5*y6*y7+sqrt(3._ark)*y6**2*y7/3._ark)*y4**2)*y2+(-sqrt(3._ark)*y6**2*y9/6._ark+ &
          sqrt(3._ark)*y5**2*y9/2._ark)*y4*y3**2+(-sqrt(3._ark)*y6**2*y9/6._ark+ &
          sqrt(3._ark)*y5**2*y9/2._ark)*y4**2*y3
      dF(254) = (2._ark/3._ark*sqrt(3._ark)*y2**2*y6**2+(y5*y6+sqrt(3._ark)*y5**2/2._ark+ &
          sqrt(3._ark)*y6**2/6._ark)*y3**2+(-y5*y6+sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/ &
          2._ark)*y4**2)*y1**2+((-y5*y6+sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/2._ark)*y3**2+ &
          (y5*y6+sqrt(3._ark)*y5**2/2._ark+sqrt(3._ark)*y6**2/6._ark)*y4**2)*y2**2+2._ark/ &
          3._ark*sqrt(3._ark)*y3**2*y4**2*y6**2
      dF(255) = (y2*y5*y7*y8+(sqrt(3._ark)*y6*y8*y9/2._ark-y5*y8*y9/2._ark)*y3+(- &
          sqrt(3._ark)*y6*y7*y9/2._ark-y5*y7*y9/2._ark)*y4)*y1**2+(y2**2*y5*y7*y8+ &
          (sqrt(3._ark)*y6*y8*y9/2._ark-y5*y8*y9/2._ark)*y3**2+(-sqrt(3._ark)*y6*y7*y9/2._ark- &
          y5*y7*y9/2._ark)*y4**2)*y1+((sqrt(3._ark)*y6*y7*y9/2._ark+y5*y7*y9/2._ark)*y3+(- &
          sqrt(3._ark)*y6*y8*y9/2._ark+y5*y8*y9/2._ark)*y4)*y2**2+((sqrt(3._ark)*y6*y7*y9/2._ark+ &
          y5*y7*y9/2._ark)*y3**2+(-sqrt(3._ark)*y6*y8*y9/2._ark+y5*y8*y9/2._ark)*y4**2)*y2- &
          y3**2*y4*y5*y7*y8-y3*y4**2*y5*y7*y8
      dF(256) = ((y3*y7*y9+y4*y8*y9)*y2+y3*y4*y7*y8)*y1**2+((-y4*y7*y9- &
          y3*y8*y9)*y2**2+(-y3**2*y7*y8-y4**2*y7*y8)*y2-y3*y4**2*y8*y9-y3**2*y4*y7*y9)*y1+ &
          y2**2*y3*y4*y7*y8+(y3**2*y4*y8*y9+y3*y4**2*y7*y9)*y2
      dF(257) = (((y7**2+y9**2)*y3+(y9**2+y8**2)*y4)*y2+(y8**2+y7**2)*y4*y3)*y1**2+ &
          (((y9**2+y8**2)*y3+(y7**2+y9**2)*y4)*y2**2+((y8**2+y7**2)*y3**2+(y8**2+ &
          y7**2)*y4**2)*y2+(y7**2+y9**2)*y4*y3**2+(y9**2+y8**2)*y4**2*y3)*y1+(y8**2+ &
          y7**2)*y4*y3*y2**2+((y9**2+y8**2)*y4*y3**2+(y7**2+y9**2)*y4**2*y3)*y2
      s1 = (((sqrt(3._ark)*y5*y9/2._ark+(y7+y9/2._ark)*y6)*y3+(sqrt(3._ark)*y5*y9/2._ark+(-y9/ &
          2._ark-y8)*y6)*y4)*y2+((-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y5+(y7/2._ark-y8/ &
          2._ark)*y6)*y4*y3)*y1**2+(((sqrt(3._ark)*y5*y9/2._ark+(y8-y9/2._ark)*y6)*y3+ &
          (sqrt(3._ark)*y5*y9/2._ark+(-y7+y9/2._ark)*y6)*y4)*y2**2+(((-sqrt(3._ark)*y7/2._ark+ &
          sqrt(3._ark)*y8/2._ark)*y5+(y7/2._ark+y8/2._ark)*y6)*y3**2+((sqrt(3._ark)*y7/2._ark- &
          sqrt(3._ark)*y8/2._ark)*y5+(-y8/2._ark-y7/2._ark)*y6)*y4**2)*y2+(-sqrt(3._ark)*y5*y9/2._ark+ &
          (-y9/2._ark+y7)*y6)*y4*y3**2+(-sqrt(3._ark)*y5*y9/2._ark+(y9/2._ark-y8)*y6)*y4**2*y3)*y1
      dF(258) = s1+((sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y8/2._ark-y7/ &
          2._ark)*y6)*y4*y3*y2**2+((-sqrt(3._ark)*y5*y9/2._ark+(y8+y9/2._ark)*y6)*y4*y3**2+(- &
          sqrt(3._ark)*y5*y9/2._ark+(-y9/2._ark-y7)*y6)*y4**2*y3)*y2
      s1 = ((((-y9/2._ark+y7)*y5+sqrt(3._ark)*y6*y9/2._ark)*y3+((y8-y9/2._ark)*y5- &
          sqrt(3._ark)*y6*y9/2._ark)*y4)*y2+((-y8/2._ark-y7/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark+ &
          sqrt(3._ark)*y8/2._ark)*y6)*y4*y3)*y1**2+((((-y9/2._ark-y8)*y5-sqrt(3._ark)*y6*y9/ &
          2._ark)*y3+((-y9/2._ark-y7)*y5+sqrt(3._ark)*y6*y9/2._ark)*y4)*y2**2+(((y8/2._ark-y7/ &
          2._ark)*y5+(-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y6)*y3**2+((y7/2._ark-y8/ &
          2._ark)*y5+(sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y6)*y4**2)*y2+((y7+y9/2._ark)*y5- &
          sqrt(3._ark)*y6*y9/2._ark)*y4*y3**2+((y8+y9/2._ark)*y5+sqrt(3._ark)*y6*y9/ &
          2._ark)*y4**2*y3)*y1
      dF(259) = s1+((y7/2._ark+y8/2._ark)*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
          2._ark)*y6)*y4*y3*y2**2+(((y9/2._ark-y8)*y5+sqrt(3._ark)*y6*y9/2._ark)*y4*y3**2+((-y7+ &
          y9/2._ark)*y5-sqrt(3._ark)*y6*y9/2._ark)*y4**2*y3)*y2
      dF(260) = (((y5*y6+sqrt(3._ark)*y6**2/3._ark)*y3+(-y5*y6+sqrt(3._ark)*y6**2/ &
          3._ark)*y4)*y2+(-sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/2._ark)*y4*y3)*y1**2+(((- &
          y5*y6+sqrt(3._ark)*y6**2/3._ark)*y3+(y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4)*y2**2+((- &
          sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/2._ark)*y3**2+(-sqrt(3._ark)*y6**2/6._ark+ &
          sqrt(3._ark)*y5**2/2._ark)*y4**2)*y2+(y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4*y3**2+(-y5*y6+ &
          sqrt(3._ark)*y6**2/3._ark)*y4**2*y3)*y1+(-sqrt(3._ark)*y6**2/6._ark+sqrt(3._ark)*y5**2/ &
          2._ark)*y4*y3*y2**2+((-y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4*y3**2+(y5*y6+ &
          sqrt(3._ark)*y6**2/3._ark)*y4**2*y3)*y2
      dF(261) = (y9+y7+y8)*y4*y3*y2*y1**2+((-y7+y9-y8)*y4*y3*y2**2+((-y8-y9+ &
          y7)*y4*y3**2+(-y7-y9+y8)*y4**2*y3)*y2)*y1
      dF(262) = (y4**2*y7*y9+y2**2*y7*y8+y3**2*y8*y9)*y1**2+(-y4**2*y8*y9- &
          y3**2*y7*y9)*y2**2-y3**2*y4**2*y7*y8
      dF(263) = (y2**2*y5*y9+(-y5*y7/2._ark+sqrt(3._ark)*y6*y7/2._ark)*y3**2+(-y5*y8/2._ark- &
          sqrt(3._ark)*y6*y8/2._ark)*y4**2)*y1**2+((y5*y8/2._ark+sqrt(3._ark)*y6*y8/2._ark)*y3**2+ &
          (y5*y7/2._ark-sqrt(3._ark)*y6*y7/2._ark)*y4**2)*y2**2-y3**2*y4**2*y5*y9
      dF(264) = ((y6**2+y5**2)*y2**2+(y6**2+y5**2)*y3**2+(y6**2+y5**2)*y4**2)*y1**2+ &
          ((y6**2+y5**2)*y3**2+(y6**2+y5**2)*y4**2)*y2**2+(y6**2+y5**2)*y4**2*y3**2
      dF(265) = ((y3*y9+y4*y9)*y2**2+(y4**2*y8+y3**2*y7)*y2+y3*y4**2*y8+ &
          y3**2*y4*y7)*y1**2+((-y3**2*y8-y4**2*y7)*y2**2-y3**2*y4**2*y9)*y1+(-y3*y4**2*y7- &
          y3**2*y4*y8)*y2**2-y2*y3**2*y4**2*y9
      dF(266) = (((y7-y8)*y3+(y8-y7)*y4)*y2**2+((-y8+y9)*y3**2+(y9-y7)*y4**2)*y2+(y8- &
          y9)*y4*y3**2+(y7-y9)*y4**2*y3)*y1**2+(((y7+y9)*y3**2+(y9+y8)*y4**2)*y2**2+(y7+ &
          y8)*y4**2*y3**2)*y1+((-y9-y7)*y4*y3**2+(-y9-y8)*y4**2*y3)*y2**2+(-y8- &
          y7)*y4**2*y3**2*y2
      dF(267) = ((y4*y5+y3*y5)*y2**2+((-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y3**2+(-y5/2._ark- &
          sqrt(3._ark)*y6/2._ark)*y4**2)*y2+(-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y4*y3**2+(-y5/2._ark- &
          sqrt(3._ark)*y6/2._ark)*y4**2*y3)*y1**2+(((-y5/2._ark-sqrt(3._ark)*y6/2._ark)*y3**2+(-y5/ &
          2._ark+sqrt(3._ark)*y6/2._ark)*y4**2)*y2**2+y3**2*y4**2*y5)*y1+((-y5/2._ark- &
          sqrt(3._ark)*y6/2._ark)*y4*y3**2+(-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y4**2*y3)*y2**2+ &
          y2*y3**2*y4**2*y5
      dF(268) = (y2**2*y3*y4+(y3*y4**2+y3**2*y4)*y2)*y1**2+((y3*y4**2+y3**2*y4)*y2**2+ &
          y2*y3**2*y4**2)*y1
      dF(269) = ((y3**2+y4**2)*y2**2+y3**2*y4**2)*y1**2+y2**2*y3**2*y4**2
      dF(270) = (((y8-y9/2._ark)*y7-y8*y9/2._ark)*y5+(-sqrt(3._ark)*y7*y9/2._ark+ &
          sqrt(3._ark)*y8*y9/2._ark)*y6)*y1**3+(((y8+y9/2._ark)*y7+y8*y9/2._ark)*y5+(- &
          sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y2**3+(((y9/2._ark-y8)*y7-y8*y9/ &
          2._ark)*y5+(sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y3**3+(((-y9/2._ark- &
          y8)*y7+y8*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y6)*y4**3
      dF(271) = (-sqrt(3._ark)*y5**2*y9/2._ark+(y7-y8)*y6*y5+(sqrt(3._ark)*y9/6._ark- &
          sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/3._ark)*y6**2)*y1**3+(-sqrt(3._ark)*y5**2*y9/2._ark+ &
          (y8-y7)*y6*y5+(sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/ &
          3._ark)*y6**2)*y2**3+(sqrt(3._ark)*y5**2*y9/2._ark+(y7+y8)*y6*y5+(sqrt(3._ark)*y8/3._ark- &
          sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y9/6._ark)*y6**2)*y3**3+(sqrt(3._ark)*y5**2*y9/2._ark+(- &
          y8-y7)*y6*y5+(-sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y8/ &
          3._ark)*y6**2)*y4**3
      dF(272) = ((y9+y7+y8)*y5**2+(y9+y7+y8)*y6**2)*y1**3+((-y7+y9-y8)*y5**2+(-y7+y9- &
          y8)*y6**2)*y2**3+((-y8-y9+y7)*y5**2+(-y8-y9+y7)*y6**2)*y3**3+((-y7-y9+y8)*y5**2+ &
          (-y7-y9+y8)*y6**2)*y4**3
      dF(273) = (y5**3-3._ark*y5*y6**2)*y1**3+(y5**3-3._ark*y5*y6**2)*y2**3+(y5**3- &
          3._ark*y5*y6**2)*y3**3+(y5**3-3._ark*y5*y6**2)*y4**3
      dF(274) = (y4**3+y3**3+y2**3)*y1**3+(y4**3+y3**3)*y2**3+y3**3*y4**3
      dF(275) = (2._ark/3._ark*sqrt(3._ark)*y2*y5*y9+(-sqrt(3._ark)*y5*y7/3._ark+y6*y7)*y3+(- &
          sqrt(3._ark)*y5*y8/3._ark-y6*y8)*y4)*y1**3+(2._ark/3._ark*sqrt(3._ark)*y2**3*y5*y9+(- &
          sqrt(3._ark)*y5*y7/3._ark+y6*y7)*y3**3+(-sqrt(3._ark)*y5*y8/3._ark-y6*y8)*y4**3)*y1+ &
          ((sqrt(3._ark)*y5*y8/3._ark+y6*y8)*y3+(-y6*y7+sqrt(3._ark)*y5*y7/3._ark)*y4)*y2**3+ &
          ((sqrt(3._ark)*y5*y8/3._ark+y6*y8)*y3**3+(-y6*y7+sqrt(3._ark)*y5*y7/3._ark)*y4**3)*y2- &
          2._ark/3._ark*sqrt(3._ark)*y3*y4**3*y5*y9-2._ark/3._ark*sqrt(3._ark)*y3**3*y4*y5*y9
      dF(276) = ((y3*y8+y4*y7)*y2+y3*y4*y9)*y1**3+((-y3*y7-y4*y8)*y2**3+(-y3**3*y9- &
          y4**3*y9)*y2-y3**3*y4*y8-y3*y4**3*y7)*y1+y2**3*y3*y4*y9+(y3*y4**3*y8+ &
          y3**3*y4*y7)*y2
      dF(277) = (((y7+y9)*y3+(y9+y8)*y4)*y2+(y7+y8)*y4*y3)*y1**3+(((-y8+y9)*y3+(y9- &
          y7)*y4)*y2**3+((y7-y8)*y3**3+(y8-y7)*y4**3)*y2+(y7-y9)*y4*y3**3+(y8- &
          y9)*y4**3*y3)*y1+(-y8-y7)*y4*y3*y2**3+((-y9-y8)*y4*y3**3+(-y9-y7)*y4**3*y3)*y2
      dF(278) = (((-y5/2._ark-sqrt(3._ark)*y6/2._ark)*y3+(-y5/2._ark+sqrt(3._ark)*y6/ &
          2._ark)*y4)*y2+y3*y4*y5)*y1**3+(((-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y3+(-y5/2._ark- &
          sqrt(3._ark)*y6/2._ark)*y4)*y2**3+(y3**3*y5+y4**3*y5)*y2+(-y5/2._ark-sqrt(3._ark)*y6/ &
          2._ark)*y4*y3**3+(-y5/2._ark+sqrt(3._ark)*y6/2._ark)*y4**3*y3)*y1+y2**3*y3*y4*y5+((-y5/ &
          2._ark+sqrt(3._ark)*y6/2._ark)*y4*y3**3+(-y5/2._ark-sqrt(3._ark)*y6/2._ark)*y4**3*y3)*y2
      dF(279) = ((y7+y8)*y2**2+(y9+y8)*y3**2+(y7+y9)*y4**2)*y1**3+((-y8-y7)*y2**3+(- &
          y9-y8)*y3**3+(-y9-y7)*y4**3)*y1**2+((y9-y7)*y3**2+(-y8+y9)*y4**2)*y2**3+((y7- &
          y9)*y3**3+(y8-y9)*y4**3)*y2**2+(y7-y8)*y4**2*y3**3+(y8-y7)*y4**3*y3**2
      dF(280) = (y4**2*y8+y3**2*y7+y2**2*y9)*y1**3+(y3**3*y7+y2**3*y9+y4**3*y8)*y1**2+ &
          (-y3**2*y8-y4**2*y7)*y2**3+(-y4**3*y7-y3**3*y8)*y2**2-y3**3*y4**2*y9- &
          y3**2*y4**3*y9
      dF(281) = (2._ark/3._ark*sqrt(3._ark)*y2**2*y5+(-sqrt(3._ark)*y5/3._ark+y6)*y3**2+(- &
          sqrt(3._ark)*y5/3._ark-y6)*y4**2)*y1**3+(2._ark/3._ark*sqrt(3._ark)*y2**3*y5+(- &
          sqrt(3._ark)*y5/3._ark+y6)*y3**3+(-sqrt(3._ark)*y5/3._ark-y6)*y4**3)*y1**2+((- &
          sqrt(3._ark)*y5/3._ark-y6)*y3**2+(-sqrt(3._ark)*y5/3._ark+y6)*y4**2)*y2**3+((- &
          sqrt(3._ark)*y5/3._ark-y6)*y3**3+(-sqrt(3._ark)*y5/3._ark+y6)*y4**3)*y2**2+2._ark/ &
          3._ark*sqrt(3._ark)*y3**3*y4**2*y5+2._ark/3._ark*sqrt(3._ark)*y3**2*y4**3*y5
      dF(282) = ((y4+y3)*y2**2+(y3**2+y4**2)*y2+y3*y4**2+y3**2*y4)*y1**3+((y4+ &
          y3)*y2**3+(y4**3+y3**3)*y2+y3**3*y4+y3*y4**3)*y1**2+((y3**2+y4**2)*y2**3+(y4**3+ &
          y3**3)*y2**2+y3**3*y4**2+y3**2*y4**3)*y1+(y3*y4**2+y3**2*y4)*y2**3+(y3**3*y4+ &
          y3*y4**3)*y2**2+(y3**3*y4**2+y3**2*y4**3)*y2
      dF(283) = y1**3*y2*y3*y4+(y2**3*y3*y4+(y3**3*y4+y3*y4**3)*y2)*y1
      dF(284) = (y8**2+y7**2+y9**2)*y1**4+(y8**2+y7**2+y9**2)*y2**4+(y8**2+y7**2+ &
          y9**2)*y3**4+(y8**2+y7**2+y9**2)*y4**4
      dF(285) = (y6**2+y5**2)*y1**4+(y6**2+y5**2)*y2**4+(y6**2+y5**2)*y3**4+(y6**2+ &
          y5**2)*y4**4
      dF(286) = (y3**2+y2**2+y4**2)*y1**4+(y4**4+y2**4+y3**4)*y1**2+(y3**2+ &
          y4**2)*y2**4+(y4**4+y3**4)*y2**2+y3**4*y4**2+y3**2*y4**4
      dF(287) = ((y4+y3)*y2+y3*y4)*y1**4+((y4+y3)*y2**4+(y4**4+y3**4)*y2+y3*y4**4+ &
          y3**4*y4)*y1+y2**4*y3*y4+(y3*y4**4+y3**4*y4)*y2
      dF(288) = (y9+y7+y8)*y1**5+(-y7+y9-y8)*y2**5+(-y8-y9+y7)*y3**5+(-y7-y9+y8)*y4**5
      dF(289) = y3**6+y4**6+y1**6+y2**6
      !
      !dF(290:292) = 0
      !
      !dF(290)  = (a12**2*t12+a13**2*t13+a14**2*t14+a23**2*t23+a24**2*t24+a34**2*t34)
      !dF(291)  = (a12**3*t12+a13**3*t13+a14**3*t14+a23**3*t23+a24**3*t24+a34**3*t34)
      !dF(292)  = (a12**4*t12+a13**4*t13+a14**4*t14+a23**4*t23+a24**4*t24+a34**4*t34)
      
      !dF(10:289) = vdump*dF(10:289)
      !
      !
 end subroutine potch4_diff_V


!
! CH4 DMS ccsd(t)-f12/cc-vpTZ 
! Yurchenko, November 2011
! 
  !
  function MLdipole_xy4_dF(xyz) result(f)
    !
    real(ark), intent(in) :: xyz(molec%natoms, 3)
    real(ark)             :: f(3)
    !
    integer(ik),parameter :: n = 681, lspace = 150
    integer(ik) :: i,k,nparams,ierror,rank0
    real(ark)   :: dF(3,n),local(10),x(4,3),tmat(4,3),rmat(3,3),dipin(3),dipout(3)
    real(rk)    :: dip_rk(3, 1), rmat_rk(3, 3), tsing(3), wspace(lspace),tol = -1.0d-12
     !
     forall(i=1:4) x(i,:) = xyz(i+1,:)-xyz(1,:)
     !  
     local(1) = sqrt(sum(x(1,:)**2))
     local(2) = sqrt(sum(x(2,:)**2))
     local(3) = sqrt(sum(x(3,:)**2))
     local(4) = sqrt(sum(x(4,:)**2))
     !
     local(5) = acos(sum(x(1,:)*x(2,:))/(local(1)*local(2)))
     local(6) = acos(sum(x(1,:)*x(3,:))/(local(1)*local(3)))
     local(7) = acos(sum(x(1,:)*x(4,:))/(local(1)*local(4)))
     local(8) = acos(sum(x(2,:)*x(3,:))/(local(2)*local(3)))
     local(9) = acos(sum(x(2,:)*x(4,:))/(local(2)*local(4)))
     local(10)= acos(sum(x(3,:)*x(4,:))/(local(3)*local(4)))
     !
     tmat(1, :) = x(1,:) / local(1)
     tmat(2, :) = x(2,:) / local(2)
     tmat(3, :) = x(3,:) / local(3)
     tmat(4, :) = x(4,:) / local(4)
     !
     rmat(1,:) = 1.d0/2.d0*(tmat(1,:)-tmat(2,:)+tmat(3,:)-tmat(4,:))
     rmat(2,:) = 1.d0/2.d0*(tmat(1,:)-tmat(2,:)-tmat(3,:)+tmat(4,:))
     rmat(3,:) = 1.d0/2.d0*(tmat(1,:)+tmat(2,:)-tmat(3,:)-tmat(4,:))
     !
     forall(i=1:3) rmat(i,:) = rmat(i,:)/sqrt(sum(rmat(i,:)**2))
     !
     call dipch4_diff_mu(n,local,dF)
     !
     f = 0
     !
     nparams = extF%nterms(1)
     !
     do i = 2,nparams
       !
       k = extF%ifit(i,1)
       !
       dipin(:) = dF(:,k)
       !
       call MLlinurark(3,rmat,dipin,dipout,ierror)
       !
       if (ierror>0) then
         !
         rmat_rk = rmat
         dip_rk(:,1) = dipin(:)
         !
         call dgelss(3,3,1,rmat_rk,3,dip_rk,3,tsing,tol,rank0,wspace,lspace,ierror)
         !
         dipout = real(dip_rk(:,1),ark)
         !
         if (ierror>0) then
           !
           print *,i,k,ierror,rmat,dipin
           write(out,"('MLdipole_xy4_dF: dgelss error = ',i)") ierror
           stop 'MLdipole_xy4_dF: dgelss error'
           !
         endif
         !
       endif
       !
       f(:) = f(:) + extF%coef(i,1)*dipout(:)
       !
     enddo
     !
  end function MLdipole_xy4_dF


  subroutine dipch4_diff_mu(n,local,dF)
    !
    integer(ik),intent(in)  :: n
    real(ark),intent(in)  :: local(10)
    real(ark),intent(out) :: dF(3,n)
    !
    real(ark) :: r1,r2,r3,r4,alpha12,alpha13,alpha23,alpha14,alpha24,alpha34
    !
    real(ark) :: re,beta,y1,y2,y3,y4,y5,y6,y7,y8,y9
    !
    beta = 1.0d0
    !
    re = extF%coef(1,1)
    !
    r1      = local(1)     
    r2      = local(2)
    r3      = local(3)
    r4      = local(4)
    alpha12 = local(5)
    alpha13 = local(6)
    alpha14 = local(7)
    alpha23 = local(8)
    alpha24 = local(9)
    alpha34 = local(10)
    !
    y1=1.0d0*(r1-re) *exp(-beta*(r1-re)**2)
    y2=1.0d0*(r2-re) *exp(-beta*(r2-re)**2)
    y3=1.0d0*(r3-re) *exp(-beta*(r3-re)**2)
    y4=1.0d0*(r4-re) *exp(-beta*(r4-re)**2)
    !
    !y1=(r1-re)
    !y2=(r2-re)
    !y3=(r3-re)
    !y4=(r4-re)
    !
    y5=(2.0d0*alpha12-alpha13-alpha14-alpha23-alpha24+2.0d0*alpha34)/sqrt(12.0d0)
    y6=(alpha13-alpha14-alpha23+alpha24)*0.5d0
    y7=(alpha24-alpha13)/sqrt(2.0d0)
    y8=(alpha23-alpha14)/sqrt(2.0d0)
    y9=(alpha34-alpha12)/sqrt(2.0d0)
    !
    dF(1,1) = 0._ark
    dF(1,2) = y7
    dF(1,3) = -y4-y2+y1+y3
    dF(1,4) = y8*y9
    dF(1,5) = -sqrt(3._ark)*y6*y7+y5*y7
    dF(1,6) = y1**2-y2**2-y4**2+y3**2
    dF(1,7) = y2*y7+y3*y7+y1*y7+y4*y7
    dF(1,8) = (y8+y9)*y1+(y8-y9)*y2+(-y8-y9)*y3+(-y8+y9)*y4
    dF(1,9) = (sqrt(3._ark)*y6/2._ark-y5/2._ark)*y1+(-sqrt(3._ark)*y6/2._ark+y5/2._ark)*y2+ &
      (sqrt(3._ark)*y6/2._ark-y5/2._ark)*y3+(-sqrt(3._ark)*y6/2._ark+y5/2._ark)*y4
    dF(1,10) = -y2*y4+y1*y3
    dF(1,11) = (y8**2+y9**2)*y7
    dF(1,12) = y7**3
    dF(1,13) = -y5*y8*y9/2._ark+sqrt(3._ark)*y6*y8*y9/2._ark
    dF(1,14) = -sqrt(3.D0)*y6**2*y7/3._ark+y5*y6*y7
    dF(1,15) = y6**2*y7+y5**2*y7
    dF(1,16) = (-y8-y9)*y7*y1+(y8-y9)*y7*y2+(y8+y9)*y7*y3+(-y8+y9)*y7*y4
    dF(1,17) = y4*y8*y9+y1*y8*y9+y3*y8*y9+y2*y8*y9
    dF(1,18) = y4*y7**2-y1*y7**2-y3*y7**2+y2*y7**2
    dF(1,19) = (-y9**2-y8**2)*y1+(y8**2+y9**2)*y2+(-y9**2-y8**2)*y3+(y8**2+y9**2)*y4
    dF(1,20) = (y6*y7-sqrt(3._ark)*y5*y7/3._ark)*y1+(y6*y7-sqrt(3._ark)*y5*y7/3._ark)*y2+ &
      (y6*y7-sqrt(3._ark)*y5*y7/3._ark)*y3+(y6*y7-sqrt(3._ark)*y5*y7/3._ark)*y4
    dF(1,21) = ((-sqrt(3._ark)*y8/2._ark+sqrt(3._ark)*y9/2._ark)*y5+(y9/2._ark-y8/ &
      2._ark)*y6)*y1+((-sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y9/2._ark)*y5+(-y8/2._ark-y9/ &
      2._ark)*y6)*y2+((-sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y8/2._ark-y9/ &
      2._ark)*y6)*y3+((sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y8/2._ark+y9/ &
      2._ark)*y6)*y4
    dF(1,22) = ((-y8/2._ark-y9/2._ark)*y5+(sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/ &
      2._ark)*y6)*y1+((y9/2._ark-y8/2._ark)*y5+(-sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/ &
      2._ark)*y6)*y2+((y8/2._ark+y9/2._ark)*y5+(-sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y9/ &
      2._ark)*y6)*y3+((y8/2._ark-y9/2._ark)*y5+(-sqrt(3._ark)*y8/2._ark+sqrt(3._ark)*y9/ &
      2._ark)*y6)*y4
    dF(1,23) = (-y5*y6-sqrt(3._ark)*y5**2/3._ark)*y1+(y5*y6+sqrt(3._ark)*y5**2/3._ark)*y2+(- &
      y5*y6-sqrt(3._ark)*y5**2/3._ark)*y3+(y5*y6+sqrt(3._ark)*y5**2/3._ark)*y4
    dF(1,24) = (-sqrt(3._ark)*y6/2._ark+y5/2._ark)*y1**2+(sqrt(3._ark)*y6/2._ark-y5/ &
      2._ark)*y2**2+(-sqrt(3._ark)*y6/2._ark+y5/2._ark)*y3**2+(sqrt(3._ark)*y6/2._ark-y5/ &
      2._ark)*y4**2
    dF(1,25) = y2**3-y3**3-y1**3+y4**3
    dF(1,26) = y2*y4*y7+y1*y3*y7
    dF(1,27) = (y2*y7+y4*y7)*y1+y2*y3*y7+y3*y4*y7
    dF(1,28) = (-y2*y8-y4*y9)*y1+y3*y4*y8+y2*y3*y9
    dF(1,29) = (-sqrt(3._ark)*y6/2._ark+y5/2._ark)*y3*y1+(sqrt(3._ark)*y6/2._ark-y5/ &
      2._ark)*y4*y2
    dF(1,30) = (y4+y2)*y1**2+(-y4**2-y2**2)*y1-y2**2*y3+y3**2*y4-y3*y4**2+y2*y3**2
    dF(1,31) = y2**2*y4+y2*y4**2-y1**2*y3-y1*y3**2
    dF(1,32) = (-y8-y9)*y1**2+(-y8+y9)*y2**2+(y8+y9)*y3**2+(y8-y9)*y4**2
    dF(1,33) = y2**2*y7+y3**2*y7+y4**2*y7+y1**2*y7
    dF(1,34) = (-y6**2-y5**2)*y1+(y6**2+y5**2)*y2+(-y6**2-y5**2)*y3+(y6**2+y5**2)*y4
    dF(1,35) = ((y4-y3)*y2-y3*y4)*y1+y2*y3*y4
    dF(1,36) = y7**2*y8*y9
    dF(1,37) = y8*y9**3+y8**3*y9
    dF(1,38) = (sqrt(3._ark)*y8**2-sqrt(3._ark)*y9**2)*y7*y5+(y8**2-y9**2)*y7*y6
    dF(1,39) = (-2._ark*y8**2+y9**2)*y7*y5+sqrt(3._ark)*y6*y7*y9**2
    dF(1,40) = y5*y7**3-sqrt(3._ark)*y6*y7**3
    dF(1,41) = -y5*y6*y8*y9-sqrt(3._ark)*y6**2*y8*y9/6._ark-sqrt(3._ark)*y5**2*y8*y9/2._ark
    dF(1,42) = y6**2*y8*y9+y5**2*y8*y9
    dF(1,43) = y5**2*y6*y7-sqrt(3._ark)*y5**3*y7+y6**3*y7+5._ark/ &
      3._ark*sqrt(3._ark)*y5*y6**2*y7
    dF(1,44) = y5**3*y7-3._ark*y5*y6**2*y7
    dF(1,45) = y4*y7**3+y3*y7**3+y2*y7**3+y1*y7**3
    dF(1,46) = ((sqrt(3._ark)*y9/6._ark-sqrt(3._ark)*y8/3._ark)*y5**2+y5*y6*y8- &
      sqrt(3._ark)*y6**2*y9/2._ark)*y1+((-sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y9/6._ark)*y5**2+ &
      y5*y6*y8+sqrt(3._ark)*y6**2*y9/2._ark)*y2+((-sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y8/ &
      3._ark)*y5**2-y5*y6*y8+sqrt(3._ark)*y6**2*y9/2._ark)*y3+((sqrt(3._ark)*y9/6._ark+ &
      sqrt(3._ark)*y8/3._ark)*y5**2-y5*y6*y8-sqrt(3._ark)*y6**2*y9/2._ark)*y4
    dF(1,47) = ((-sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y8/2._ark-y9/ &
      2._ark)*y6)*y1**2+((sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y8/2._ark+y9/ &
      2._ark)*y6)*y2**2+((-sqrt(3._ark)*y8/2._ark+sqrt(3._ark)*y9/2._ark)*y5+(y9/2._ark-y8/ &
      2._ark)*y6)*y3**2+((-sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y9/2._ark)*y5+(-y8/2._ark-y9/ &
      2._ark)*y6)*y4**2
    dF(1,48) = (-y9**2-y8**2)*y1**2+(y8**2+y9**2)*y2**2+(-y9**2-y8**2)*y3**2+(y8**2+ &
      y9**2)*y4**2
    dF(1,49) = (sqrt(3._ark)*y6/2._ark-y5/2._ark)*y1**3+(-sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y2**3+(sqrt(3._ark)*y6/2._ark-y5/2._ark)*y3**3+(-sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y4**3
    dF(1,50) = y2**4-y3**4+y4**4-y1**4
    dF(1,51) = (y8**2+y9**2)*y7*y1+(y8**2+y9**2)*y7*y2+(y8**2+y9**2)*y7*y3+(y8**2+ &
      y9**2)*y7*y4
    dF(1,52) = (y8**2*y9+y8*y9**2)*y1+(y8*y9**2-y8**2*y9)*y2+(-y8*y9**2- &
      y8**2*y9)*y3+(-y8*y9**2+y8**2*y9)*y4
    dF(1,53) = (y8+y9)*y7**2*y1+(y8-y9)*y7**2*y2+(-y8-y9)*y7**2*y3+(-y8+y9)*y7**2*y4
    dF(1,54) = (y8**3+y9**3)*y1+(-y9**3+y8**3)*y2+(-y9**3-y8**3)*y3+(y9**3-y8**3)*y4
    dF(1,55) = y1*y7*y8*y9-y2*y7*y8*y9-y4*y7*y8*y9+y3*y7*y8*y9
    dF(1,56) = (sqrt(3._ark)*y5*y8**2/2._ark+(-y9**2-y8**2/2._ark)*y6)*y1+(- &
      sqrt(3._ark)*y5*y8**2/2._ark+(y9**2+y8**2/2._ark)*y6)*y2+(sqrt(3._ark)*y5*y8**2/2._ark+(- &
      y9**2-y8**2/2._ark)*y6)*y3+(-sqrt(3._ark)*y5*y8**2/2._ark+(y9**2+y8**2/2._ark)*y6)*y4
    dF(1,57) = (y6*y8*y9-sqrt(3._ark)*y5*y8*y9/3._ark)*y1+(y6*y8*y9-sqrt(3._ark)*y5*y8*y9/ &
      3._ark)*y2+(y6*y8*y9-sqrt(3._ark)*y5*y8*y9/3._ark)*y3+(y6*y8*y9-sqrt(3._ark)*y5*y8*y9/ &
      3._ark)*y4
    dF(1,58) = (sqrt(3._ark)*y5*y7*y8+(-y8-2._ark*y9)*y7*y6)*y1+(-sqrt(3._ark)*y5*y7*y8+ &
      (y8-2._ark*y9)*y7*y6)*y2+(-sqrt(3._ark)*y5*y7*y8+(y8+2._ark*y9)*y7*y6)*y3+ &
      (sqrt(3._ark)*y5*y7*y8+(-y8+2._ark*y9)*y7*y6)*y4
    dF(1,59) = ((-y8/3._ark+2._ark/3._ark*y9)*y5**2-2._ark/3._ark*sqrt(3._ark)*y5*y6*y9+ &
      y6**2*y8)*y1+((-2._ark/3._ark*y9-y8/3._ark)*y5**2+2._ark/3._ark*sqrt(3._ark)*y5*y6*y9+ &
      y6**2*y8)*y2+((-2._ark/3._ark*y9+y8/3._ark)*y5**2+2._ark/3._ark*sqrt(3._ark)*y5*y6*y9- &
      y6**2*y8)*y3+((2._ark/3._ark*y9+y8/3._ark)*y5**2-2._ark/3._ark*sqrt(3._ark)*y5*y6*y9- &
      y6**2*y8)*y4
    dF(1,60) = ((y9/3._ark+4._ark/3._ark*y8)*y5**2+2._ark/3._ark*sqrt(3._ark)*y5*y6*y9+ &
      y6**2*y9)*y1+((4._ark/3._ark*y8-y9/3._ark)*y5**2-2._ark/3._ark*sqrt(3._ark)*y5*y6*y9- &
      y6**2*y9)*y2+((-4._ark/3._ark*y8-y9/3._ark)*y5**2-2._ark/3._ark*sqrt(3._ark)*y5*y6*y9- &
      y6**2*y9)*y3+((y9/3._ark-4._ark/3._ark*y8)*y5**2+2._ark/3._ark*sqrt(3._ark)*y5*y6*y9+ &
      y6**2*y9)*y4
    dF(1,61) = ((-y8**2/2._ark+y9**2)*y5-sqrt(3._ark)*y6*y8**2/2._ark)*y1+((y8**2/2._ark- &
      y9**2)*y5+sqrt(3._ark)*y6*y8**2/2._ark)*y2+((-y8**2/2._ark+y9**2)*y5- &
      sqrt(3._ark)*y6*y8**2/2._ark)*y3+((y8**2/2._ark-y9**2)*y5+sqrt(3._ark)*y6*y8**2/2._ark)*y4
    dF(1,62) = (-y5*y7**2/2._ark+sqrt(3._ark)*y6*y7**2/2._ark)*y1+(y5*y7**2/2._ark- &
      sqrt(3._ark)*y6*y7**2/2._ark)*y2+(-y5*y7**2/2._ark+sqrt(3._ark)*y6*y7**2/2._ark)*y3+ &
      (y5*y7**2/2._ark-sqrt(3._ark)*y6*y7**2/2._ark)*y4
    dF(1,63) = ((-2._ark*y8+y9)*y7*y5+sqrt(3._ark)*y6*y7*y9)*y1+((y9+2._ark*y8)*y7*y5+ &
      sqrt(3._ark)*y6*y7*y9)*y2+((2._ark*y8-y9)*y7*y5-sqrt(3._ark)*y6*y7*y9)*y3+((-y9- &
      2._ark*y8)*y7*y5-sqrt(3._ark)*y6*y7*y9)*y4
    dF(1,64) = (-sqrt(3._ark)*y6**2*y7/3._ark+y5*y6*y7)*y1+(-sqrt(3._ark)*y6**2*y7/3._ark+ &
      y5*y6*y7)*y2+(-sqrt(3._ark)*y6**2*y7/3._ark+y5*y6*y7)*y3+(-sqrt(3._ark)*y6**2*y7/3._ark+ &
      y5*y6*y7)*y4
    dF(1,65) = (y5*y6**2-y5**3/3._ark)*y1+(y5**3/3._ark-y5*y6**2)*y2+(y5*y6**2-y5**3/ &
      3._ark)*y3+(y5**3/3._ark-y5*y6**2)*y4
    dF(1,66) = (y6**2*y7+y5**2*y7)*y1+(y6**2*y7+y5**2*y7)*y2+(y6**2*y7+y5**2*y7)*y3+ &
      (y6**2*y7+y5**2*y7)*y4
    dF(1,67) = (sqrt(3._ark)*y5*y6**2+sqrt(3._ark)*y5**3/9._ark-y6**3-y5**2*y6)*y1+ &
      (y5**2*y6-sqrt(3._ark)*y5**3/9._ark-sqrt(3._ark)*y5*y6**2+y6**3)*y2+ &
      (sqrt(3._ark)*y5*y6**2+sqrt(3._ark)*y5**3/9._ark-y6**3-y5**2*y6)*y3+(y5**2*y6- &
      sqrt(3._ark)*y5**3/9._ark-sqrt(3._ark)*y5*y6**2+y6**3)*y4
    dF(1,68) = (-y5*y6+sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y3*y1+(- &
      sqrt(3._ark)*y6**2/2._ark-sqrt(3._ark)*y5**2/6._ark+y5*y6)*y4*y2
    dF(1,69) = (sqrt(3._ark)*y5/3._ark-y6)*y3*y1**2+(sqrt(3._ark)*y5/3._ark-y6)*y3**2*y1+ &
      (y6-sqrt(3._ark)*y5/3._ark)*y4*y2**2+(y6-sqrt(3._ark)*y5/3._ark)*y4**2*y2
    dF(1,70) = y1**3*y3+y1*y3**3-y2**3*y4-y2*y4**3
    dF(1,71) = (-y2-y4)*y1**3+(y2**3+y4**3)*y1-y2*y3**3-y3**3*y4+y2**3*y3+y3*y4**3
    dF(1,72) = y1*y3*y7**2-y2*y4*y7**2
    dF(1,73) = (y2*y7*y9+y4*y7*y8)*y1-y3*y4*y7*y9-y2*y3*y7*y8
    dF(1,74) = (y8**2+y9**2)*y3*y1+(-y9**2-y8**2)*y4*y2
    dF(1,75) = (y2*y8*y9+y4*y8*y9)*y1+y3*y4*y8*y9+y2*y3*y8*y9
    dF(1,76) = y2*y4*y8*y9+y1*y3*y8*y9
    dF(1,77) = (y2*y6*y7+(y6*y7/2._ark-sqrt(3._ark)*y5*y7/2._ark)*y4)*y1+(y6*y7/2._ark- &
      sqrt(3._ark)*y5*y7/2._ark)*y3*y2+y3*y4*y6*y7
    dF(1,78) = (-y2*y6*y8+(-y6*y9/2._ark+sqrt(3._ark)*y5*y9/2._ark)*y4)*y1+(y6*y9/2._ark- &
      sqrt(3._ark)*y5*y9/2._ark)*y3*y2+y3*y4*y6*y8
    dF(1,79) = (y6**2+y5**2)*y3*y1+(-y6**2-y5**2)*y4*y2
    dF(1,80) = (-y5*y7/2._ark+sqrt(3._ark)*y6*y7/2._ark)*y3*y1+(-y5*y7/2._ark+ &
      sqrt(3._ark)*y6*y7/2._ark)*y4*y2
    dF(1,81) = (y2*y5*y7+(-sqrt(3._ark)*y6*y7/2._ark-y5*y7/2._ark)*y4)*y1+(- &
      sqrt(3._ark)*y6*y7/2._ark-y5*y7/2._ark)*y3*y2+y3*y4*y5*y7
    dF(1,82) = (y2*y5*y8+(-y5*y9/2._ark-sqrt(3._ark)*y6*y9/2._ark)*y4)*y1+ &
      (sqrt(3._ark)*y6*y9/2._ark+y5*y9/2._ark)*y3*y2-y3*y4*y5*y8
    dF(1,83) = ((y4*y7+y3*y7)*y2+y3*y4*y7)*y1+y2*y3*y4*y7
    dF(1,84) = (((sqrt(3._ark)*y6/2._ark-y5/2._ark)*y3+(-sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y4)*y2+(sqrt(3._ark)*y6/2._ark-y5/2._ark)*y4*y3)*y1+(-sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y4*y3*y2
    dF(1,85) = (((y8-y9)*y3+(y8+y9)*y4)*y2+(-y8+y9)*y4*y3)*y1+(-y8-y9)*y4*y3*y2
    dF(1,86) = (-y2*y5+(sqrt(3._ark)*y6/2._ark+y5/2._ark)*y4)*y1**2+(y2**2*y5+(- &
      sqrt(3._ark)*y6/2._ark-y5/2._ark)*y4**2)*y1+(-sqrt(3._ark)*y6/2._ark-y5/2._ark)*y3*y2**2+ &
      (sqrt(3._ark)*y6/2._ark+y5/2._ark)*y3**2*y2+y3*y4**2*y5-y3**2*y4*y5
    dF(1,87) = (y8+y9)*y7*y1**2+(-y8+y9)*y7*y2**2+(-y8-y9)*y7*y3**2+(y8-y9)*y7*y4**2
    dF(1,88) = y4**2*y8*y9+y3**2*y8*y9+y2**2*y8*y9+y1**2*y8*y9
    dF(1,89) = y1**2*y7**2+y3**2*y7**2-y4**2*y7**2-y2**2*y7**2
    dF(1,90) = (y6*y7-sqrt(3._ark)*y5*y7/3._ark)*y1**2+(y6*y7-sqrt(3._ark)*y5*y7/ &
      3._ark)*y2**2+(y6*y7-sqrt(3._ark)*y5*y7/3._ark)*y3**2+(y6*y7-sqrt(3._ark)*y5*y7/ &
      3._ark)*y4**2
    dF(1,91) = (y6**2+y5**2)*y1**2+(-y6**2-y5**2)*y2**2+(y6**2+y5**2)*y3**2+(-y6**2- &
      y5**2)*y4**2
    dF(1,92) = ((-y8/2._ark-y9/2._ark)*y5+(sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/ &
      2._ark)*y6)*y1**2+((y9/2._ark-y8/2._ark)*y5+(-sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/ &
      2._ark)*y6)*y2**2+((y8/2._ark+y9/2._ark)*y5+(-sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y9/ &
      2._ark)*y6)*y3**2+((y8/2._ark-y9/2._ark)*y5+(-sqrt(3._ark)*y8/2._ark+sqrt(3._ark)*y9/ &
      2._ark)*y6)*y4**2
    dF(1,93) = (-y5*y6-sqrt(3._ark)*y5**2/3._ark)*y1**2+(y5*y6+sqrt(3._ark)*y5**2/ &
      3._ark)*y2**2+(-y5*y6-sqrt(3._ark)*y5**2/3._ark)*y3**2+(y5*y6+sqrt(3._ark)*y5**2/ &
      3._ark)*y4**2
    dF(1,94) = (y4*y9+y2*y8)*y1**2+(y2**2*y8+y4**2*y9)*y1-y3*y4**2*y8-y3**2*y4*y8- &
      y2**2*y3*y9-y2*y3**2*y9
    dF(1,95) = (y8+y9)*y3*y1**2+(-y8-y9)*y3**2*y1+(y8-y9)*y4*y2**2+(-y8+y9)*y4**2*y2
    dF(1,96) = (y2*y7+y4*y7)*y1**2+(y2**2*y7+y4**2*y7)*y1+y3**2*y4*y7+y2*y3**2*y7+ &
      y2**2*y3*y7+y3*y4**2*y7
    dF(1,97) = y1**2*y3*y7+y2*y4**2*y7+y2**2*y4*y7+y1*y3**2*y7
    dF(1,98) = (y4*y8+y2*y9)*y1**2+(-y2**2*y9-y4**2*y8)*y1+y3*y4**2*y9-y3**2*y4*y9+ &
      y2**2*y3*y8-y2*y3**2*y8
    dF(1,99) = y1**2*y3**2-y2**2*y4**2
    dF(1,100) = (y2*y6+(-sqrt(3._ark)*y5/2._ark+y6/2._ark)*y4)*y1**2+(-y2**2*y6+(-y6/2._ark+ &
      sqrt(3._ark)*y5/2._ark)*y4**2)*y1+(-y6/2._ark+sqrt(3._ark)*y5/2._ark)*y3*y2**2+(- &
      sqrt(3._ark)*y5/2._ark+y6/2._ark)*y3**2*y2+y3**2*y4*y6-y3*y4**2*y6
    dF(1,101) = (y2*y3+y3*y4)*y1**2+(-y2**2*y4+(y3**2-y4**2)*y2+y3**2*y4)*y1- &
      y2**2*y3*y4-y2*y3*y4**2
    dF(1,102) = y1**2*y2*y4+(-y2**2*y3-y3*y4**2)*y1+y2*y3**2*y4
    dF(1,103) = (y8+y9)*y1**3+(y8-y9)*y2**3+(-y8-y9)*y3**3+(-y8+y9)*y4**3
    dF(1,104) = y4**3*y7+y3**3*y7+y2**3*y7+y1**3*y7
    dF(1,105) = (y8**2+y9**2)*y7**3
    dF(1,106) = y7*y8**2*y9**2
    dF(1,107) = (y9**4+y8**4)*y7
    dF(1,108) = y7**5
    dF(1,109) = -sqrt(3._ark)*y5*y8**3*y9+(2._ark*y8*y9**3+y8**3*y9)*y6
    dF(1,110) = 3._ark*y5**2*y7*y9**2-2._ark*sqrt(3._ark)*y5*y6*y7*y9**2+(y9**2+ &
      4._ark*y8**2)*y7*y6**2
    dF(1,111) = y6**4*y7+y5**4*y7/3._ark-4._ark/3._ark*sqrt(3._ark)*y5**3*y6*y7
    dF(1,112) = (y8*y9**3+y8**3*y9)*y5+(-sqrt(3._ark)*y8*y9**3-sqrt(3._ark)*y8**3*y9)*y6
    dF(1,113) = y5*y7**2*y8*y9-sqrt(3._ark)*y6*y7**2*y8*y9
    dF(1,114) = -sqrt(3._ark)*y5**2*y7*y9**2+(y9**2-y8**2)*y7*y6*y5- &
      sqrt(3._ark)*y6**2*y7*y8**2
    dF(1,115) = y5*y6*y7**3-sqrt(3._ark)*y6**2*y7**3/3._ark
    dF(1,116) = y5*y6**3*y7+2._ark/9._ark*sqrt(3._ark)*y5**4*y7+y5**3*y6*y7/3._ark
    dF(1,117) = (-2._ark*y9**2+y8**2)*y7*y5**2+2._ark*sqrt(3._ark)*y5*y6*y7*y9**2- &
      3._ark*y6**2*y7*y8**2
    dF(1,118) = y5**2*y7**3+y6**2*y7**3
    dF(1,119) = -5._ark/3._ark*sqrt(3._ark)*y5*y6**2*y8*y9-y5**2*y6*y8*y9+ &
      sqrt(3._ark)*y5**3*y8*y9-y6**3*y8*y9
    dF(1,120) = y5**2*y6**2*y7+2._ark/3._ark*sqrt(3._ark)*y5**3*y6*y7+y5**4*y7/3._ark
    dF(1,121) = -3._ark*y5*y6**2*y8*y9+y5**3*y8*y9
    dF(1,122) = (y9**4+y8**4)*y1+(-y8**4-y9**4)*y2+(y9**4+y8**4)*y3+(-y8**4- &
      y9**4)*y4
    dF(1,123) = (-y9**3-y8**3)*y7*y1+(-y9**3+y8**3)*y7*y2+(y8**3+y9**3)*y7*y3+ &
      (y9**3-y8**3)*y7*y4
    dF(1,124) = (sqrt(3._ark)*y5*y9**3/2._ark+(-y9**3/2._ark-y8**3)*y6)*y1+(- &
      sqrt(3._ark)*y5*y9**3/2._ark+(y9**3/2._ark-y8**3)*y6)*y2+(-sqrt(3._ark)*y5*y9**3/2._ark+ &
      (y9**3/2._ark+y8**3)*y6)*y3+(sqrt(3._ark)*y5*y9**3/2._ark+(-y9**3/2._ark+y8**3)*y6)*y4
    dF(1,125) = (sqrt(3._ark)*y5**3*y9/2._ark+(-y8-7._ark/2._ark*y9)*y6*y5**2+(5._ark/ &
      2._ark*sqrt(3._ark)*y9+2._ark*sqrt(3._ark)*y8)*y6**2*y5+(-3._ark/2._ark*y9- &
      3._ark*y8)*y6**3)*y1+(-sqrt(3._ark)*y5**3*y9/2._ark+(7._ark/2._ark*y9-y8)*y6*y5**2+ &
      (2._ark*sqrt(3._ark)*y8-5._ark/2._ark*sqrt(3._ark)*y9)*y6**2*y5+(3._ark/2._ark*y9- &
      3._ark*y8)*y6**3)*y2+(-sqrt(3._ark)*y5**3*y9/2._ark+(7._ark/2._ark*y9+y8)*y6*y5**2+(- &
      2._ark*sqrt(3._ark)*y8-5._ark/2._ark*sqrt(3._ark)*y9)*y6**2*y5+(3._ark*y8+3._ark/ &
      2._ark*y9)*y6**3)*y3+(sqrt(3._ark)*y5**3*y9/2._ark+(y8-7._ark/2._ark*y9)*y6*y5**2+(5._ark/ &
      2._ark*sqrt(3._ark)*y9-2._ark*sqrt(3._ark)*y8)*y6**2*y5+(-3._ark/2._ark*y9+ &
      3._ark*y8)*y6**3)*y4
    dF(1,126) = y1**3*y8*y9+y4**3*y8*y9+y3**3*y8*y9+y2**3*y8*y9
    dF(1,127) = ((y9**3/2._ark-y8**3)*y5+sqrt(3._ark)*y6*y9**3/2._ark)*y1+((-y9**3/2._ark- &
      y8**3)*y5-sqrt(3._ark)*y6*y9**3/2._ark)*y2+((-y9**3/2._ark+y8**3)*y5- &
      sqrt(3._ark)*y6*y9**3/2._ark)*y3+((y9**3/2._ark+y8**3)*y5+sqrt(3._ark)*y6*y9**3/2._ark)*y4
    dF(1,128) = ((-y9**2-y8**2)*y5**2+(-y9**2-y8**2)*y6**2)*y1+((y8**2+y9**2)*y5**2+ &
      (y8**2+y9**2)*y6**2)*y2+((-y9**2-y8**2)*y5**2+(-y9**2-y8**2)*y6**2)*y3+((y8**2+ &
      y9**2)*y5**2+(y8**2+y9**2)*y6**2)*y4
    dF(1,129) = (y5**3*y7-3._ark*y5*y6**2*y7)*y1+(y5**3*y7-3._ark*y5*y6**2*y7)*y2+ &
      (y5**3*y7-3._ark*y5*y6**2*y7)*y3+(y5**3*y7-3._ark*y5*y6**2*y7)*y4
    dF(1,130) = ((7._ark/2._ark*y9-y8)*y5**3-15._ark/2._ark*sqrt(3._ark)*y5**2*y6*y9+(9._ark*y8+ &
      27._ark/2._ark*y9)*y6**2*y5+(-3._ark/2._ark*sqrt(3._ark)*y9-6._ark*sqrt(3._ark)*y8)*y6**3)*y1+ &
      ((-y8-7._ark/2._ark*y9)*y5**3+15._ark/2._ark*sqrt(3._ark)*y5**2*y6*y9+(9._ark*y8-27._ark/ &
      2._ark*y9)*y6**2*y5+(-6._ark*sqrt(3._ark)*y8+3._ark/2._ark*sqrt(3._ark)*y9)*y6**3)*y2+((y8- &
      7._ark/2._ark*y9)*y5**3+15._ark/2._ark*sqrt(3._ark)*y5**2*y6*y9+(-9._ark*y8-27._ark/ &
      2._ark*y9)*y6**2*y5+(6._ark*sqrt(3._ark)*y8+3._ark/2._ark*sqrt(3._ark)*y9)*y6**3)*y3+((7._ark/ &
      2._ark*y9+y8)*y5**3-15._ark/2._ark*sqrt(3._ark)*y5**2*y6*y9+(27._ark/2._ark*y9- &
      9._ark*y8)*y6**2*y5+(6._ark*sqrt(3._ark)*y8-3._ark/2._ark*sqrt(3._ark)*y9)*y6**3)*y4
    dF(1,131) = ((-y8**2+2._ark*y9**2)*y5-sqrt(3._ark)*y6*y8**2)*y3*y1+((-2._ark*y9**2+ &
      y8**2)*y5+sqrt(3._ark)*y6*y8**2)*y4*y2
    dF(1,132) = (y2*y5+(-sqrt(3._ark)*y6/2._ark-y5/2._ark)*y4)*y1**3+(-y2**3*y5+ &
      (sqrt(3._ark)*y6/2._ark+y5/2._ark)*y4**3)*y1+(sqrt(3._ark)*y6/2._ark+y5/2._ark)*y3*y2**3+(- &
      sqrt(3._ark)*y6/2._ark-y5/2._ark)*y3**3*y2+y3**3*y4*y5-y3*y4**3*y5
    dF(1,133) = (-sqrt(3._ark)*y6/2._ark+y5/2._ark)*y3*y1**3+(-sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y3**3*y1+(sqrt(3._ark)*y6/2._ark-y5/2._ark)*y4*y2**3+(sqrt(3._ark)*y6/2._ark-y5/ &
      2._ark)*y4**3*y2
    dF(1,134) = ((y8+y9)*y5**2+(y8+y9)*y6**2)*y1**2+((y8-y9)*y5**2+(y8- &
      y9)*y6**2)*y2**2+((-y8-y9)*y5**2+(-y8-y9)*y6**2)*y3**2+((-y8+y9)*y5**2+(-y8+ &
      y9)*y6**2)*y4**2
    dF(1,135) = ((-sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y8/2._ark-y9/ &
      2._ark)*y6)*y1**3+((sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y8/2._ark+y9/ &
      2._ark)*y6)*y2**3+((-sqrt(3._ark)*y8/2._ark+sqrt(3._ark)*y9/2._ark)*y5+(y9/2._ark-y8/ &
      2._ark)*y6)*y3**3+((-sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y9/2._ark)*y5+(-y8/2._ark-y9/ &
      2._ark)*y6)*y4**3
    dF(1,136) = (y6**2+y5**2)*y1**3+(-y6**2-y5**2)*y2**3+(y6**2+y5**2)*y3**3+(- &
      y6**2-y5**2)*y4**3
    dF(1,137) = (sqrt(3._ark)*y5*y8**2*y9/2._ark+(-y8*y9**2-y8**2*y9/2._ark)*y6)*y1+(- &
      sqrt(3._ark)*y5*y8**2*y9/2._ark+(-y8*y9**2+y8**2*y9/2._ark)*y6)*y2+(- &
      sqrt(3._ark)*y5*y8**2*y9/2._ark+(y8**2*y9/2._ark+y8*y9**2)*y6)*y3+ &
      (sqrt(3._ark)*y5*y8**2*y9/2._ark+(-y8**2*y9/2._ark+y8*y9**2)*y6)*y4
    dF(1,138) = (y6*y7*y8*y9-sqrt(3._ark)*y5*y7*y8*y9/3._ark)*y1+(-y6*y7*y8*y9+ &
      sqrt(3._ark)*y5*y7*y8*y9/3._ark)*y2+(y6*y7*y8*y9-sqrt(3._ark)*y5*y7*y8*y9/3._ark)*y3+(- &
      y6*y7*y8*y9+sqrt(3._ark)*y5*y7*y8*y9/3._ark)*y4
    dF(1,139) = (-sqrt(3._ark)*y5**2*y9**2/2._ark-y5*y6*y8**2+(sqrt(3._ark)*y9**2/6._ark- &
      sqrt(3._ark)*y8**2/3._ark)*y6**2)*y1+(sqrt(3._ark)*y5**2*y9**2/2._ark+y5*y6*y8**2+ &
      (sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y9**2/6._ark)*y6**2)*y2+(- &
      sqrt(3._ark)*y5**2*y9**2/2._ark-y5*y6*y8**2+(sqrt(3._ark)*y9**2/6._ark-sqrt(3._ark)*y8**2/ &
      3._ark)*y6**2)*y3+(sqrt(3._ark)*y5**2*y9**2/2._ark+y5*y6*y8**2+(sqrt(3._ark)*y8**2/3._ark- &
      sqrt(3._ark)*y9**2/6._ark)*y6**2)*y4
    dF(1,140) = (y4*y7*y9+y2*y7*y8)*y1**2+(-y4**2*y7*y9-y2**2*y7*y8)*y1+ &
      y2**2*y3*y7*y9+y3*y4**2*y7*y8-y3**2*y4*y7*y8-y2*y3**2*y7*y9
    dF(1,141) = (-y9**2-y8**2)*y3*y1**2+(-y9**2-y8**2)*y3**2*y1+(y8**2+ &
      y9**2)*y4*y2**2+(y8**2+y9**2)*y4**2*y2
    dF(1,142) = (y2*y7**2+y4*y7**2)*y1**2+(-y2**2*y7**2-y4**2*y7**2)*y1- &
      y3*y4**2*y7**2+y3**2*y4*y7**2+y2*y3**2*y7**2-y2**2*y3*y7**2
    dF(1,143) = (-y2*y8**2-y4*y9**2)*y1**2+(y4**2*y9**2+y2**2*y8**2)*y1- &
      y2*y3**2*y9**2+y3*y4**2*y8**2+y2**2*y3*y9**2-y3**2*y4*y8**2
    dF(1,144) = (y8+y9)*y7*y3*y1**2+(-y8-y9)*y7*y3**2*y1+(-y8+y9)*y7*y4*y2**2+(y8- &
      y9)*y7*y4**2*y2
    dF(1,145) = (-y4*y7*y8-y2*y7*y9)*y1**2+(-y2**2*y7*y9-y4**2*y7*y8)*y1+ &
      y2*y3**2*y7*y8+y3**2*y4*y7*y9+y2**2*y3*y7*y8+y3*y4**2*y7*y9
    dF(1,146) = (y2*y8*y9+y4*y8*y9)*y1**2+(y4**2*y8*y9+y2**2*y8*y9)*y1+ &
      y3*y4**2*y8*y9+y3**2*y4*y8*y9+y2*y3**2*y8*y9+y2**2*y3*y8*y9
    dF(1,147) = (2._ark*y2*y6*y9+(-sqrt(3._ark)*y5*y8+y6*y8)*y4)*y1**2+(- &
      2._ark*y2**2*y6*y9+(sqrt(3._ark)*y5*y8-y6*y8)*y4**2)*y1+(-sqrt(3._ark)*y5*y8+ &
      y6*y8)*y3*y2**2+(sqrt(3._ark)*y5*y8-y6*y8)*y3**2*y2-2._ark*y3**2*y4*y6*y9+ &
      2._ark*y3*y4**2*y6*y9
    dF(1,148) = ((-sqrt(3._ark)*y6*y7/2._ark-y5*y7/2._ark)*y2+y4*y5*y7)*y1**2+((- &
      sqrt(3._ark)*y6*y7/2._ark-y5*y7/2._ark)*y2**2+y4**2*y5*y7)*y1+y2**2*y3*y5*y7+ &
      y2*y3**2*y5*y7+(-sqrt(3._ark)*y6*y7/2._ark-y5*y7/2._ark)*y4*y3**2+(-sqrt(3._ark)*y6*y7/ &
      2._ark-y5*y7/2._ark)*y4**2*y3
    dF(1,149) = ((-sqrt(3._ark)*y6**2/6._ark-y5*y6-sqrt(3._ark)*y5**2/2._ark)*y2+(- &
      sqrt(3._ark)*y6**2/6._ark-y5*y6-sqrt(3._ark)*y5**2/2._ark)*y4)*y1**2+((y5*y6+ &
      sqrt(3._ark)*y5**2/2._ark+sqrt(3._ark)*y6**2/6._ark)*y2**2+(y5*y6+sqrt(3._ark)*y5**2/2._ark+ &
      sqrt(3._ark)*y6**2/6._ark)*y4**2)*y1+(y5*y6+sqrt(3._ark)*y5**2/2._ark+sqrt(3._ark)*y6**2/ &
      6._ark)*y3*y2**2+(-sqrt(3._ark)*y6**2/6._ark-y5*y6-sqrt(3._ark)*y5**2/2._ark)*y3**2*y2+(- &
      sqrt(3._ark)*y6**2/6._ark-y5*y6-sqrt(3._ark)*y5**2/2._ark)*y4*y3**2+(y5*y6+ &
      sqrt(3._ark)*y5**2/2._ark+sqrt(3._ark)*y6**2/6._ark)*y4**2*y3
    dF(1,150) = ((-sqrt(3._ark)*y6**2/3._ark-y5*y6)*y2+(sqrt(3._ark)*y6**2/6._ark- &
      sqrt(3._ark)*y5**2/2._ark)*y4)*y1**2+((y5*y6+sqrt(3._ark)*y6**2/3._ark)*y2**2+ &
      (sqrt(3._ark)*y5**2/2._ark-sqrt(3._ark)*y6**2/6._ark)*y4**2)*y1+(sqrt(3._ark)*y5**2/2._ark- &
      sqrt(3._ark)*y6**2/6._ark)*y3*y2**2+(sqrt(3._ark)*y6**2/6._ark-sqrt(3._ark)*y5**2/ &
      2._ark)*y3**2*y2+(-sqrt(3._ark)*y6**2/3._ark-y5*y6)*y4*y3**2+(y5*y6+sqrt(3._ark)*y6**2/ &
      3._ark)*y4**2*y3
    dF(1,151) = (-y6**2-y5**2)*y3*y1**2+(-y6**2-y5**2)*y3**2*y1+(y6**2+ &
      y5**2)*y4*y2**2+(y6**2+y5**2)*y4**2*y2
    dF(1,152) = (y8+y9)*y3*y1**3+(-y8-y9)*y3**3*y1+(y8-y9)*y4*y2**3+(-y8+ &
      y9)*y4**3*y2
    dF(1,153) = (-y2*y8-y4*y9)*y1**3+(-y4**3*y9-y2**3*y8)*y1+y2**3*y3*y9+ &
      y3*y4**3*y8+y3**3*y4*y8+y2*y3**3*y9
    dF(1,154) = y2*y4**3*y7+y2**3*y4*y7+y1**3*y3*y7+y1*y3**3*y7
    dF(1,155) = (y4*y8+y2*y9)*y1**3+(-y4**3*y8-y2**3*y9)*y1+y2**3*y3*y8+y3*y4**3*y9- &
      y2*y3**3*y8-y3**3*y4*y9
    dF(1,156) = (-y2-y4)*y1**4+(y4**4+y2**4)*y1+y3*y4**4-y3**4*y4-y2*y3**4+y2**4*y3
    dF(1,157) = (y4*y7**3+y2*y7**3)*y1+y3*y4*y7**3+y2*y3*y7**3
    dF(1,158) = (y8**2+y9**2)*y7*y3*y1+(y8**2+y9**2)*y7*y4*y2
    dF(1,159) = (-y2*y8*y9**2-y4*y8**2*y9)*y1+y2*y3*y8**2*y9+y3*y4*y8*y9**2
    dF(1,160) = y2*y4*y7**3+y1*y3*y7**3
    dF(1,161) = -y1*y3*y7*y8*y9+y2*y4*y7*y8*y9
    dF(1,162) = (y4*y7*y8**2+y2*y7*y9**2)*y1+y2*y3*y7*y8**2+y3*y4*y7*y9**2
    dF(1,163) = (-y2*y7**2*y8-y4*y7**2*y9)*y1+y2*y3*y7**2*y9+y3*y4*y7**2*y8
    dF(1,164) = (y4*y7*y9**2+y2*y7*y8**2)*y1+y3*y4*y7*y8**2+y2*y3*y7*y9**2
    dF(1,165) = (-y2*y8**3-y4*y9**3)*y1+y3*y4*y8**3+y2*y3*y9**3
    dF(1,166) = ((-sqrt(3._ark)*y8**2+sqrt(3._ark)*y9**2)*y5+(y9**2-y8**2)*y6)*y3*y1+ &
      ((sqrt(3._ark)*y8**2-sqrt(3._ark)*y9**2)*y5+(y8**2-y9**2)*y6)*y4*y2
    dF(1,167) = ((sqrt(3._ark)*y5*y7*y9-y6*y7*y9)*y2-2._ark*y4*y6*y7*y8)*y1+ &
      2._ark*y2*y3*y6*y7*y8+(y6*y7*y9-sqrt(3._ark)*y5*y7*y9)*y4*y3
    dF(1,168) = ((-y5*y7*y9+sqrt(3._ark)*y6*y7*y9)*y2+(sqrt(3._ark)*y6*y7*y8- &
      y5*y7*y8)*y4)*y1+(y5*y7*y8-sqrt(3._ark)*y6*y7*y8)*y3*y2+(-sqrt(3._ark)*y6*y7*y9+ &
      y5*y7*y9)*y4*y3
    dF(1,169) = ((sqrt(3._ark)*y5**2*y7/4._ark-sqrt(3._ark)*y6**2*y7/4._ark+y5*y6*y7/ &
      2._ark)*y2+y4*y5*y6*y7)*y1+y2*y3*y5*y6*y7+(sqrt(3._ark)*y5**2*y7/4._ark- &
      sqrt(3._ark)*y6**2*y7/4._ark+y5*y6*y7/2._ark)*y4*y3
    dF(1,170) = ((y5*y6*y8-sqrt(3._ark)*y5**2*y8/3._ark)*y2+(-sqrt(3._ark)*y6**2*y9/2._ark+ &
      sqrt(3._ark)*y5**2*y9/6._ark)*y4)*y1+(sqrt(3._ark)*y6**2*y9/2._ark-sqrt(3._ark)*y5**2*y9/ &
      6._ark)*y3*y2+(sqrt(3._ark)*y5**2*y8/3._ark-y5*y6*y8)*y4*y3
    dF(1,171) = (sqrt(3._ark)*y5/3._ark-y6)*y4*y2*y1**2+((y6-sqrt(3._ark)*y5/ &
      3._ark)*y3*y2**2+(y6-sqrt(3._ark)*y5/3._ark)*y4**2*y3)*y1+(sqrt(3._ark)*y5/3._ark- &
      y6)*y4*y3**2*y2
    dF(1,172) = (y3*y4*y7+y2*y3*y7)*y1**2+(y2**2*y4*y7+(y3**2*y7+y4**2*y7)*y2+ &
      y3**2*y4*y7)*y1+y2*y3*y4**2*y7+y2**2*y3*y4*y7
    dF(1,173) = (y3*y4*y8+y2*y3*y9)*y1**2+(-y2**2*y4*y9+(-y4**2*y8-y3**2*y8)*y2- &
      y3**2*y4*y9)*y1+y2**2*y3*y4*y8+y2*y3*y4**2*y9
    dF(1,174) = ((y6+sqrt(3._ark)*y5)*y3*y2+(-y6-sqrt(3._ark)*y5)*y4*y3)*y1**2+((-y6- &
      sqrt(3._ark)*y5)*y4*y2**2+((-y6-sqrt(3._ark)*y5)*y3**2+(y6+sqrt(3._ark)*y5)*y4**2)*y2+ &
      (y6+sqrt(3._ark)*y5)*y4*y3**2)*y1+(y6+sqrt(3._ark)*y5)*y4*y3*y2**2+(-y6- &
      sqrt(3._ark)*y5)*y4**2*y3*y2
    dF(1,175) = y1**3*y2*y4+(-y2**3*y3-y3*y4**3)*y1+y2*y3**3*y4
    dF(1,176) = (-y2*y3-y3*y4)*y1**3+(y2**3*y4+(-y3**3+y4**3)*y2-y3**3*y4)*y1+ &
      y2*y3*y4**3+y2**3*y3*y4
    dF(1,177) = (y8**2*y9+y8*y9**2)*y1**2+(y8*y9**2-y8**2*y9)*y2**2+(-y8*y9**2- &
      y8**2*y9)*y3**2+(-y8*y9**2+y8**2*y9)*y4**2
    dF(1,178) = (y8**2+y9**2)*y7*y1**2+(y8**2+y9**2)*y7*y2**2+(y8**2+ &
      y9**2)*y7*y3**2+(y8**2+y9**2)*y7*y4**2
    dF(1,179) = (-y8-y9)*y7**2*y1**2+(-y8+y9)*y7**2*y2**2+(y8+y9)*y7**2*y3**2+(y8- &
      y9)*y7**2*y4**2
    dF(1,180) = (sqrt(3._ark)*y5*y7*y8/2._ark+(-y8/2._ark-y9)*y7*y6)*y1**2+(- &
      sqrt(3._ark)*y5*y7*y8/2._ark+(y8/2._ark-y9)*y7*y6)*y2**2+(-sqrt(3._ark)*y5*y7*y8/2._ark+ &
      (y8/2._ark+y9)*y7*y6)*y3**2+(sqrt(3._ark)*y5*y7*y8/2._ark+(-y8/2._ark+y9)*y7*y6)*y4**2
    dF(1,181) = ((-2._ark*y9**2+y8**2)*y5+sqrt(3._ark)*y6*y8**2)*y1**2+((-y8**2+ &
      2._ark*y9**2)*y5-sqrt(3._ark)*y6*y8**2)*y2**2+((-2._ark*y9**2+y8**2)*y5+ &
      sqrt(3._ark)*y6*y8**2)*y3**2+((-y8**2+2._ark*y9**2)*y5-sqrt(3._ark)*y6*y8**2)*y4**2
    dF(1,182) = ((-y8/2._ark+y9)*y7*y5-sqrt(3._ark)*y6*y7*y8/2._ark)*y1**2+((y8/2._ark+ &
      y9)*y7*y5+sqrt(3._ark)*y6*y7*y8/2._ark)*y2**2+((y8/2._ark-y9)*y7*y5+ &
      sqrt(3._ark)*y6*y7*y8/2._ark)*y3**2+((-y8/2._ark-y9)*y7*y5-sqrt(3._ark)*y6*y7*y8/ &
      2._ark)*y4**2
    dF(1,183) = y1**2*y3*y8*y9+y1*y3**2*y8*y9+y2**2*y4*y8*y9+y2*y4**2*y8*y9
    dF(1,184) = (2._ark*y2*y6*y8+(-sqrt(3._ark)*y5*y9+y6*y9)*y4)*y1**2+ &
      (2._ark*y2**2*y6*y8+(-sqrt(3._ark)*y5*y9+y6*y9)*y4**2)*y1+(-y6*y9+ &
      sqrt(3._ark)*y5*y9)*y3*y2**2+(-y6*y9+sqrt(3._ark)*y5*y9)*y3**2*y2- &
      2._ark*y3**2*y4*y6*y8-2._ark*y3*y4**2*y6*y8
    dF(1,185) = -y1**3*y7**2-y3**3*y7**2+y4**3*y7**2+y2**3*y7**2
    dF(1,186) = (y8**2+y9**2)*y1**3+(-y9**2-y8**2)*y2**3+(y8**2+y9**2)*y3**3+(- &
      y9**2-y8**2)*y4**3
    dF(1,187) = (-y8-y9)*y7*y1**3+(y8-y9)*y7*y2**3+(y8+y9)*y7*y3**3+(-y8+ &
      y9)*y7*y4**3
    dF(1,188) = ((y8/2._ark+y9/2._ark)*y5+(-sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y9/ &
      2._ark)*y6)*y1**3+((y8/2._ark-y9/2._ark)*y5+(-sqrt(3._ark)*y8/2._ark+sqrt(3._ark)*y9/ &
      2._ark)*y6)*y2**3+((-y8/2._ark-y9/2._ark)*y5+(sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/ &
      2._ark)*y6)*y3**3+((y9/2._ark-y8/2._ark)*y5+(-sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/ &
      2._ark)*y6)*y4**3
    dF(1,189) = -y4*y8**2*y9**2-y2*y8**2*y9**2+y1*y8**2*y9**2+y3*y8**2*y9**2
    dF(1,190) = (y8**2+y9**2)*y7**2*y1+(-y9**2-y8**2)*y7**2*y2+(y8**2+ &
      y9**2)*y7**2*y3+(-y9**2-y8**2)*y7**2*y4
    dF(1,191) = (y8+y9)*y7**3*y1+(-y8+y9)*y7**3*y2+(-y8-y9)*y7**3*y3+(y8- &
      y9)*y7**3*y4
    dF(1,192) = -y4*y7**4-y2*y7**4+y3*y7**4+y1*y7**4
    dF(1,193) = (y8*y9**3+y8**3*y9)*y1+(y8*y9**3+y8**3*y9)*y2+(y8*y9**3+ &
      y8**3*y9)*y3+(y8*y9**3+y8**3*y9)*y4
    dF(1,194) = y1*y7**2*y8*y9+y4*y7**2*y8*y9+y3*y7**2*y8*y9+y2*y7**2*y8*y9
    dF(1,195) = (y8**2*y9+y8*y9**2)*y7*y1+(-y8*y9**2+y8**2*y9)*y7*y2+(-y8*y9**2- &
      y8**2*y9)*y7*y3+(y8*y9**2-y8**2*y9)*y7*y4
    dF(1,196) = (sqrt(3._ark)*y5*y7*y8**2+(-2._ark*y9**2-y8**2)*y7*y6)*y1+ &
      (sqrt(3._ark)*y5*y7*y8**2+(-2._ark*y9**2-y8**2)*y7*y6)*y2+(sqrt(3._ark)*y5*y7*y8**2+(- &
      2._ark*y9**2-y8**2)*y7*y6)*y3+(sqrt(3._ark)*y5*y7*y8**2+(-2._ark*y9**2- &
      y8**2)*y7*y6)*y4
    dF(1,197) = (-sqrt(3._ark)*y5*y7**2*y8+(y8+2._ark*y9)*y7**2*y6)*y1+(- &
      sqrt(3._ark)*y5*y7**2*y8+(y8-2._ark*y9)*y7**2*y6)*y2+(sqrt(3._ark)*y5*y7**2*y8+(-y8- &
      2._ark*y9)*y7**2*y6)*y3+(sqrt(3._ark)*y5*y7**2*y8+(-y8+2._ark*y9)*y7**2*y6)*y4
    dF(1,198) = (y5**2*y7**2+y6**2*y7**2)*y1+(-y5**2*y7**2-y6**2*y7**2)*y2+ &
      (y5**2*y7**2+y6**2*y7**2)*y3+(-y5**2*y7**2-y6**2*y7**2)*y4
    dF(1,199) = (y6**2*y8*y9+y5**2*y8*y9)*y1+(y6**2*y8*y9+y5**2*y8*y9)*y2+ &
      (y6**2*y8*y9+y5**2*y8*y9)*y3+(y6**2*y8*y9+y5**2*y8*y9)*y4
    dF(1,200) = (3._ark/2._ark*y5**2*y7*y8+sqrt(3._ark)*y5*y6*y7*y9+(-y8/2._ark+ &
      y9)*y7*y6**2)*y1+(-3._ark/2._ark*y5**2*y7*y8+sqrt(3._ark)*y5*y6*y7*y9+(y8/2._ark+ &
      y9)*y7*y6**2)*y2+(-3._ark/2._ark*y5**2*y7*y8-sqrt(3._ark)*y5*y6*y7*y9+(y8/2._ark- &
      y9)*y7*y6**2)*y3+(3._ark/2._ark*y5**2*y7*y8-sqrt(3._ark)*y5*y6*y7*y9+(-y8/2._ark- &
      y9)*y7*y6**2)*y4
    dF(1,201) = (-4._ark/3._ark*sqrt(3._ark)*y5*y6**2*y7+y5**2*y6*y7+y6**3*y7)*y1+(-4._ark/ &
      3._ark*sqrt(3._ark)*y5*y6**2*y7+y5**2*y6*y7+y6**3*y7)*y2+(-4._ark/ &
      3._ark*sqrt(3._ark)*y5*y6**2*y7+y5**2*y6*y7+y6**3*y7)*y3+(-4._ark/ &
      3._ark*sqrt(3._ark)*y5*y6**2*y7+y5**2*y6*y7+y6**3*y7)*y4
    dF(1,202) = (-3._ark/2._ark*sqrt(3._ark)*y5**3*y9+27._ark/2._ark*y5**2*y6*y9+(- &
      6._ark*sqrt(3._ark)*y8-15._ark/2._ark*sqrt(3._ark)*y9)*y6**2*y5+(7._ark/2._ark*y9+ &
      10._ark*y8)*y6**3)*y1+(3._ark/2._ark*sqrt(3._ark)*y5**3*y9-27._ark/2._ark*y5**2*y6*y9+(- &
      6._ark*sqrt(3._ark)*y8+15._ark/2._ark*sqrt(3._ark)*y9)*y6**2*y5+(10._ark*y8-7._ark/ &
      2._ark*y9)*y6**3)*y2+(3._ark/2._ark*sqrt(3._ark)*y5**3*y9-27._ark/2._ark*y5**2*y6*y9+(15._ark/ &
      2._ark*sqrt(3._ark)*y9+6._ark*sqrt(3._ark)*y8)*y6**2*y5+(-10._ark*y8-7._ark/ &
      2._ark*y9)*y6**3)*y3+(-3._ark/2._ark*sqrt(3._ark)*y5**3*y9+27._ark/2._ark*y5**2*y6*y9+ &
      (6._ark*sqrt(3._ark)*y8-15._ark/2._ark*sqrt(3._ark)*y9)*y6**2*y5+(-10._ark*y8+7._ark/ &
      2._ark*y9)*y6**3)*y4
    dF(1,203) = (3._ark*y5**2*y6**2-2._ark*sqrt(3._ark)*y5*y6**3+y6**4)*y1+(-y6**4- &
      3._ark*y5**2*y6**2+2._ark*sqrt(3._ark)*y5*y6**3)*y2+(3._ark*y5**2*y6**2- &
      2._ark*sqrt(3._ark)*y5*y6**3+y6**4)*y3+(-y6**4-3._ark*y5**2*y6**2+ &
      2._ark*sqrt(3._ark)*y5*y6**3)*y4
    dF(1,204) = (sqrt(3._ark)*y6*y7**3/2._ark-y5*y7**3/2._ark)*y1+(sqrt(3._ark)*y6*y7**3/ &
      2._ark-y5*y7**3/2._ark)*y2+(sqrt(3._ark)*y6*y7**3/2._ark-y5*y7**3/2._ark)*y3+ &
      (sqrt(3._ark)*y6*y7**3/2._ark-y5*y7**3/2._ark)*y4
    dF(1,205) = ((-2._ark*y8+y9)*y7**2*y5+sqrt(3._ark)*y6*y7**2*y9)*y1+((-y9- &
      2._ark*y8)*y7**2*y5-sqrt(3._ark)*y6*y7**2*y9)*y2+((2._ark*y8-y9)*y7**2*y5- &
      sqrt(3._ark)*y6*y7**2*y9)*y3+((y9+2._ark*y8)*y7**2*y5+sqrt(3._ark)*y6*y7**2*y9)*y4
    dF(1,206) = ((y8**2+y9**2)*y7*y5+(-sqrt(3._ark)*y8**2-sqrt(3._ark)*y9**2)*y7*y6)*y1+ &
      ((y8**2+y9**2)*y7*y5+(-sqrt(3._ark)*y8**2-sqrt(3._ark)*y9**2)*y7*y6)*y2+((y8**2+ &
      y9**2)*y7*y5+(-sqrt(3._ark)*y8**2-sqrt(3._ark)*y9**2)*y7*y6)*y3+((y8**2+ &
      y9**2)*y7*y5+(-sqrt(3._ark)*y8**2-sqrt(3._ark)*y9**2)*y7*y6)*y4
    dF(1,207) = ((-y8**2*y9/2._ark+y8*y9**2)*y5-sqrt(3._ark)*y6*y8**2*y9/2._ark)*y1+ &
      ((y8**2*y9/2._ark+y8*y9**2)*y5+sqrt(3._ark)*y6*y8**2*y9/2._ark)*y2+((-y8*y9**2+ &
      y8**2*y9/2._ark)*y5+sqrt(3._ark)*y6*y8**2*y9/2._ark)*y3+((-y8*y9**2-y8**2*y9/2._ark)*y5- &
      sqrt(3._ark)*y6*y8**2*y9/2._ark)*y4
    dF(1,208) = ((-sqrt(3._ark)*y9**2/2._ark-sqrt(3._ark)*y8**2/2._ark)*y5**2+(-y9**2- &
      y8**2)*y6*y5+(-sqrt(3._ark)*y9**2/6._ark-sqrt(3._ark)*y8**2/6._ark)*y6**2)*y1+ &
      ((sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y9**2/2._ark)*y5**2+(y8**2+y9**2)*y6*y5+ &
      (sqrt(3._ark)*y9**2/6._ark+sqrt(3._ark)*y8**2/6._ark)*y6**2)*y2+((-sqrt(3._ark)*y9**2/ &
      2._ark-sqrt(3._ark)*y8**2/2._ark)*y5**2+(-y9**2-y8**2)*y6*y5+(-sqrt(3._ark)*y9**2/6._ark- &
      sqrt(3._ark)*y8**2/6._ark)*y6**2)*y3+((sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y9**2/ &
      2._ark)*y5**2+(y8**2+y9**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark+sqrt(3._ark)*y8**2/ &
      6._ark)*y6**2)*y4
    dF(1,209) = (-sqrt(3._ark)*y5**2*y8*y9/3._ark-y5*y6*y8*y9)*y1+(- &
      sqrt(3._ark)*y5**2*y8*y9/3._ark-y5*y6*y8*y9)*y2+(-sqrt(3._ark)*y5**2*y8*y9/3._ark- &
      y5*y6*y8*y9)*y3+(-sqrt(3._ark)*y5**2*y8*y9/3._ark-y5*y6*y8*y9)*y4
    dF(1,210) = ((-sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/2._ark)*y7*y5**2+(-y8+ &
      y9)*y7*y6*y5+(-sqrt(3._ark)*y8/2._ark+sqrt(3._ark)*y9/2._ark)*y7*y6**2)*y1+((- &
      sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y9/2._ark)*y7*y5**2+(y8+y9)*y7*y6*y5+(sqrt(3._ark)*y9/ &
      2._ark+sqrt(3._ark)*y8/2._ark)*y7*y6**2)*y2+((-sqrt(3._ark)*y8/2._ark+sqrt(3._ark)*y9/ &
      2._ark)*y7*y5**2+(y8-y9)*y7*y6*y5+(-sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/ &
      2._ark)*y7*y6**2)*y3+((sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/2._ark)*y7*y5**2+(-y8- &
      y9)*y7*y6*y5+(-sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y9/2._ark)*y7*y6**2)*y4
    dF(1,211) = (y5*y6*y7**2+sqrt(3._ark)*y5**2*y7**2/3._ark)*y1+(-y5*y6*y7**2- &
      sqrt(3._ark)*y5**2*y7**2/3._ark)*y2+(y5*y6*y7**2+sqrt(3._ark)*y5**2*y7**2/3._ark)*y3+(- &
      y5*y6*y7**2-sqrt(3._ark)*y5**2*y7**2/3._ark)*y4
    dF(1,212) = (-3._ark/2._ark*y5**3*y9+5._ark/2._ark*sqrt(3._ark)*y5**2*y6*y9+(-2._ark*y8- &
      7._ark/2._ark*y9)*y6**2*y5+(2._ark*sqrt(3._ark)*y8+sqrt(3._ark)*y9/2._ark)*y6**3)*y1+(3._ark/ &
      2._ark*y5**3*y9-5._ark/2._ark*sqrt(3._ark)*y5**2*y6*y9+(7._ark/2._ark*y9-2._ark*y8)*y6**2*y5+ &
      (2._ark*sqrt(3._ark)*y8-sqrt(3._ark)*y9/2._ark)*y6**3)*y2+(3._ark/2._ark*y5**3*y9-5._ark/ &
      2._ark*sqrt(3._ark)*y5**2*y6*y9+(2._ark*y8+7._ark/2._ark*y9)*y6**2*y5+(- &
      2._ark*sqrt(3._ark)*y8-sqrt(3._ark)*y9/2._ark)*y6**3)*y3+(-3._ark/2._ark*y5**3*y9+5._ark/ &
      2._ark*sqrt(3._ark)*y5**2*y6*y9+(2._ark*y8-7._ark/2._ark*y9)*y6**2*y5+(- &
      2._ark*sqrt(3._ark)*y8+sqrt(3._ark)*y9/2._ark)*y6**3)*y4
    dF(1,213) = ((-y8/2._ark+y9)*y7*y5**2-sqrt(3._ark)*y5*y6*y7*y9+3._ark/ &
      2._ark*y6**2*y7*y8)*y1+((y8/2._ark+y9)*y7*y5**2-sqrt(3._ark)*y5*y6*y7*y9-3._ark/ &
      2._ark*y6**2*y7*y8)*y2+((y8/2._ark-y9)*y7*y5**2+sqrt(3._ark)*y5*y6*y7*y9-3._ark/ &
      2._ark*y6**2*y7*y8)*y3+((-y8/2._ark-y9)*y7*y5**2+sqrt(3._ark)*y5*y6*y7*y9+3._ark/ &
      2._ark*y6**2*y7*y8)*y4
    dF(1,214) = (-y5**3*y6-2._ark/3._ark*sqrt(3._ark)*y5**2*y6**2+y5*y6**3)*y1+(y5**3*y6+ &
      2._ark/3._ark*sqrt(3._ark)*y5**2*y6**2-y5*y6**3)*y2+(-y5**3*y6-2._ark/ &
      3._ark*sqrt(3._ark)*y5**2*y6**2+y5*y6**3)*y3+(y5**3*y6+2._ark/ &
      3._ark*sqrt(3._ark)*y5**2*y6**2-y5*y6**3)*y4
    dF(1,215) = (-y5**2*y6**2+2._ark*sqrt(3._ark)*y5*y6**3+y5**4)*y1+(-y5**4- &
      2._ark*sqrt(3._ark)*y5*y6**3+y5**2*y6**2)*y2+(-y5**2*y6**2+2._ark*sqrt(3._ark)*y5*y6**3+ &
      y5**4)*y3+(-y5**4-2._ark*sqrt(3._ark)*y5*y6**3+y5**2*y6**2)*y4
    dF(1,216) = (y6**2*y7+y5**2*y7)*y3*y1+(y6**2*y7+y5**2*y7)*y4*y2
    dF(1,217) = (-y5*y6*y7-sqrt(3._ark)*y5**2*y7/3._ark)*y3*y1+(-y5*y6*y7- &
      sqrt(3._ark)*y5**2*y7/3._ark)*y4*y2
    dF(1,218) = ((3._ark/4._ark*y6**2*y7+y5**2*y7/4._ark+sqrt(3._ark)*y5*y6*y7/2._ark)*y2+ &
      y4*y5**2*y7)*y1+y2*y3*y5**2*y7+(3._ark/4._ark*y6**2*y7+y5**2*y7/4._ark+ &
      sqrt(3._ark)*y5*y6*y7/2._ark)*y4*y3
    dF(1,219) = (y5**3-3._ark*y5*y6**2)*y3*y1+(-y5**3+3._ark*y5*y6**2)*y4*y2
    dF(1,220) = y1*y3**2*y7**2+y1**2*y3*y7**2-y2**2*y4*y7**2-y2*y4**2*y7**2
    dF(1,221) = (-y4*y8**2-y2*y9**2)*y1**2+(y2**2*y9**2+y4**2*y8**2)*y1+ &
      y3*y4**2*y9**2-y2*y3**2*y8**2+y2**2*y3*y8**2-y3**2*y4*y9**2
    dF(1,222) = ((-sqrt(3._ark)*y5*y6*y7/2._ark+y6**2*y7/4._ark+3._ark/4._ark*y5**2*y7)*y2+ &
      y4*y6**2*y7)*y1+y2*y3*y6**2*y7+(-sqrt(3._ark)*y5*y6*y7/2._ark+y6**2*y7/4._ark+3._ark/ &
      4._ark*y5**2*y7)*y4*y3
    dF(1,223) = ((y5**2*y8+y6**2*y8)*y2+(y6**2*y9+y5**2*y9)*y4)*y1+(-y6**2*y9- &
      y5**2*y9)*y3*y2+(-y6**2*y8-y5**2*y8)*y4*y3
    dF(1,224) = (y5**2*y6+y6**3-4._ark/3._ark*sqrt(3._ark)*y5*y6**2)*y3*y1+(4._ark/ &
      3._ark*sqrt(3._ark)*y5*y6**2-y6**3-y5**2*y6)*y4*y2
    dF(1,225) = (-2._ark*y2*y5*y8*y9+(sqrt(3._ark)*y6*y8*y9+y5*y8*y9)*y4)*y1+ &
      (sqrt(3._ark)*y6*y8*y9+y5*y8*y9)*y3*y2-2._ark*y3*y4*y5*y8*y9
    dF(1,226) = (-sqrt(3._ark)*y6*y7**2+y5*y7**2)*y3*y1+(-y5*y7**2+ &
      sqrt(3._ark)*y6*y7**2)*y4*y2
    dF(1,227) = ((-y5*y6*y8+sqrt(3._ark)*y6**2*y8/2._ark+sqrt(3._ark)*y5**2*y8/6._ark)*y2+(- &
      y5*y6*y9+sqrt(3._ark)*y6**2*y9/2._ark+sqrt(3._ark)*y5**2*y9/6._ark)*y4)*y1+(- &
      sqrt(3._ark)*y6**2*y9/2._ark+y5*y6*y9-sqrt(3._ark)*y5**2*y9/6._ark)*y3*y2+(y5*y6*y8- &
      sqrt(3._ark)*y6**2*y8/2._ark-sqrt(3._ark)*y5**2*y8/6._ark)*y4*y3
    dF(1,228) = (((-y8+y9)*y7*y3+(y8+y9)*y7*y4)*y2+(y8-y9)*y7*y4*y3)*y1+(-y8- &
      y9)*y7*y4*y3*y2
    dF(1,229) = ((-y4*y7**2+y3*y7**2)*y2+y3*y4*y7**2)*y1-y2*y3*y4*y7**2
    dF(1,230) = (((-y6**2-y5**2)*y3+(y6**2+y5**2)*y4)*y2+(-y6**2-y5**2)*y4*y3)*y1+ &
      (y6**2+y5**2)*y4*y3*y2
    dF(1,231) = (y8+y9)*y4*y2*y1**2+((y8-y9)*y3*y2**2+(-y8+y9)*y4**2*y3)*y1+(-y8- &
      y9)*y4*y3**2*y2
    dF(1,232) = ((-y6**2-y5**2)*y2+(-y6**2-y5**2)*y4)*y1**2+((y6**2+y5**2)*y2**2+ &
      (y6**2+y5**2)*y4**2)*y1+(y6**2+y5**2)*y3*y2**2+(-y6**2-y5**2)*y3**2*y2+(-y6**2- &
      y5**2)*y4*y3**2+(y6**2+y5**2)*y4**2*y3
    dF(1,233) = (y2*y7+y4*y7)*y1**3+(y4**3*y7+y2**3*y7)*y1+y2*y3**3*y7+y2**3*y3*y7+ &
      y3*y4**3*y7+y3**3*y4*y7
    dF(1,234) = ((sqrt(3._ark)*y5*y8*y9+y6*y8*y9)*y2+(-y6*y8*y9- &
      sqrt(3._ark)*y5*y8*y9)*y4)*y1+(-y6*y8*y9-sqrt(3._ark)*y5*y8*y9)*y3*y2+ &
      (sqrt(3._ark)*y5*y8*y9+y6*y8*y9)*y4*y3
    dF(1,235) = (-y5*y8*y9/2._ark+sqrt(3._ark)*y6*y8*y9/2._ark)*y3*y1+(-y5*y8*y9/2._ark+ &
      sqrt(3._ark)*y6*y8*y9/2._ark)*y4*y2
    dF(1,236) = (((y8**2+y9**2)*y3+(-y9**2-y8**2)*y4)*y2+(y8**2+y9**2)*y4*y3)*y1+(- &
      y9**2-y8**2)*y4*y3*y2
    dF(1,237) = ((y4*y8*y9+y3*y8*y9)*y2+y3*y4*y8*y9)*y1+y2*y3*y4*y8*y9
    dF(1,238) = ((((-sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y9/2._ark)*y5+(-y8/2._ark-y9/ &
      2._ark)*y6)*y3+((-sqrt(3._ark)*y8/2._ark+sqrt(3._ark)*y9/2._ark)*y5+(y9/2._ark-y8/ &
      2._ark)*y6)*y4)*y2+((sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y8/2._ark+y9/ &
      2._ark)*y6)*y4*y3)*y1+((-sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y8/2._ark-y9/ &
      2._ark)*y6)*y4*y3*y2
    dF(1,239) = (((-sqrt(3._ark)*y6*y7+y5*y7)*y3+(-sqrt(3._ark)*y6*y7+y5*y7)*y4)*y2+(- &
      sqrt(3._ark)*y6*y7+y5*y7)*y4*y3)*y1+(-sqrt(3._ark)*y6*y7+y5*y7)*y4*y3*y2
    dF(1,240) = ((((y8/2._ark-y9/2._ark)*y5+(-sqrt(3._ark)*y8/2._ark+sqrt(3._ark)*y9/ &
      2._ark)*y6)*y3+((y8/2._ark+y9/2._ark)*y5+(-sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y9/ &
      2._ark)*y6)*y4)*y2+((y9/2._ark-y8/2._ark)*y5+(-sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/ &
      2._ark)*y6)*y4*y3)*y1+((-y8/2._ark-y9/2._ark)*y5+(sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/ &
      2._ark)*y6)*y4*y3*y2
    dF(1,241) = (((y5*y6+sqrt(3._ark)*y5**2/2._ark+sqrt(3._ark)*y6**2/6._ark)*y3+(- &
      sqrt(3._ark)*y6**2/6._ark-y5*y6-sqrt(3._ark)*y5**2/2._ark)*y4)*y2+(y5*y6+ &
      sqrt(3._ark)*y5**2/2._ark+sqrt(3._ark)*y6**2/6._ark)*y4*y3)*y1+(-sqrt(3._ark)*y6**2/6._ark- &
      y5*y6-sqrt(3._ark)*y5**2/2._ark)*y4*y3*y2
    dF(1,242) = y1*y2*y3*y4*y7
    dF(1,243) = (-y3*y4*y9-y2*y3*y8)*y1**2+(-y2**2*y4*y8+(y3**2*y9-y4**2*y9)*y2+ &
      y3**2*y4*y8)*y1+y2**2*y3*y4*y9+y2*y3*y4**2*y8
    dF(1,244) = y3**2*y7**3+y4**2*y7**3+y2**2*y7**3+y1**2*y7**3
    dF(1,245) = (y8**3+y9**3)*y1**2+(-y9**3+y8**3)*y2**2+(-y9**3-y8**3)*y3**2+ &
      (y9**3-y8**3)*y4**2
    dF(1,246) = -y4**2*y7*y8*y9+y1**2*y7*y8*y9+y3**2*y7*y8*y9-y2**2*y7*y8*y9
    dF(1,247) = ((-sqrt(3._ark)*y8**2+sqrt(3._ark)*y9**2)*y5+(y9**2-y8**2)*y6)*y1**2+ &
      ((sqrt(3._ark)*y8**2-sqrt(3._ark)*y9**2)*y5+(y8**2-y9**2)*y6)*y2**2+((- &
      sqrt(3._ark)*y8**2+sqrt(3._ark)*y9**2)*y5+(y9**2-y8**2)*y6)*y3**2+ &
      ((sqrt(3._ark)*y8**2-sqrt(3._ark)*y9**2)*y5+(y8**2-y9**2)*y6)*y4**2
    dF(1,248) = (y6*y7**2-sqrt(3._ark)*y5*y7**2/3._ark)*y1**2+(-y6*y7**2+ &
      sqrt(3._ark)*y5*y7**2/3._ark)*y2**2+(y6*y7**2-sqrt(3._ark)*y5*y7**2/3._ark)*y3**2+(- &
      y6*y7**2+sqrt(3._ark)*y5*y7**2/3._ark)*y4**2
    dF(1,249) = (-y5**2*y6+sqrt(3._ark)*y5**3-5._ark/3._ark*sqrt(3._ark)*y5*y6**2- &
      y6**3)*y1**2+(-sqrt(3._ark)*y5**3+y6**3+5._ark/3._ark*sqrt(3._ark)*y5*y6**2+ &
      y5**2*y6)*y2**2+(-y5**2*y6+sqrt(3._ark)*y5**3-5._ark/3._ark*sqrt(3._ark)*y5*y6**2- &
      y6**3)*y3**2+(-sqrt(3._ark)*y5**3+y6**3+5._ark/3._ark*sqrt(3._ark)*y5*y6**2+ &
      y5**2*y6)*y4**2
    dF(1,250) = (-y5*y8*y9/2._ark+sqrt(3._ark)*y6*y8*y9/2._ark)*y1**2+(-y5*y8*y9/2._ark+ &
      sqrt(3._ark)*y6*y8*y9/2._ark)*y2**2+(-y5*y8*y9/2._ark+sqrt(3._ark)*y6*y8*y9/2._ark)*y3**2+ &
      (-y5*y8*y9/2._ark+sqrt(3._ark)*y6*y8*y9/2._ark)*y4**2
    dF(1,251) = ((sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y8/6._ark)*y5**2+(-y8-y9)*y6*y5+ &
      (sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/2._ark)*y6**2)*y1**2+((-sqrt(3._ark)*y9/6._ark+ &
      sqrt(3._ark)*y8/6._ark)*y5**2+(-y8+y9)*y6*y5+(-sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/ &
      2._ark)*y6**2)*y2**2+((-sqrt(3._ark)*y8/6._ark-sqrt(3._ark)*y9/6._ark)*y5**2+(y8+ &
      y9)*y6*y5+(-sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y9/2._ark)*y6**2)*y3**2+((sqrt(3._ark)*y9/ &
      6._ark-sqrt(3._ark)*y8/6._ark)*y5**2+(y8-y9)*y6*y5+(-sqrt(3._ark)*y8/2._ark+sqrt(3._ark)*y9/ &
      2._ark)*y6**2)*y4**2
    dF(1,252) = (-sqrt(3._ark)*y6**2*y7/3._ark+y5*y6*y7)*y1**2+(-sqrt(3._ark)*y6**2*y7/ &
      3._ark+y5*y6*y7)*y2**2+(-sqrt(3._ark)*y6**2*y7/3._ark+y5*y6*y7)*y3**2+(- &
      sqrt(3._ark)*y6**2*y7/3._ark+y5*y6*y7)*y4**2
    dF(1,253) = ((-sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y8/3._ark)*y5**2-y5*y6*y8+ &
      sqrt(3._ark)*y6**2*y9/2._ark)*y1**2+((sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y8/3._ark)*y5**2- &
      y5*y6*y8-sqrt(3._ark)*y6**2*y9/2._ark)*y2**2+((sqrt(3._ark)*y9/6._ark-sqrt(3._ark)*y8/ &
      3._ark)*y5**2+y5*y6*y8-sqrt(3._ark)*y6**2*y9/2._ark)*y3**2+((-sqrt(3._ark)*y8/3._ark- &
      sqrt(3._ark)*y9/6._ark)*y5**2+y5*y6*y8+sqrt(3._ark)*y6**2*y9/2._ark)*y4**2
    dF(1,254) = (y6**2*y7+y5**2*y7)*y1**2+(y6**2*y7+y5**2*y7)*y2**2+(y6**2*y7+ &
      y5**2*y7)*y3**2+(y6**2*y7+y5**2*y7)*y4**2
    dF(1,255) = (y5**3-3._ark*y5*y6**2)*y1**2+(-y5**3+3._ark*y5*y6**2)*y2**2+(y5**3- &
      3._ark*y5*y6**2)*y3**2+(-y5**3+3._ark*y5*y6**2)*y4**2
    dF(1,256) = y1**2*y3**2*y7+y2**2*y4**2*y7
    dF(1,257) = (y2**2*y7+y4**2*y7)*y1**2+y3**2*y4**2*y7+y2**2*y3**2*y7
    dF(1,258) = (y2**2*y8+y4**2*y9)*y1**2-y3**2*y4**2*y8-y2**2*y3**2*y9
    dF(1,259) = (sqrt(3._ark)*y5/3._ark-y6)*y3**2*y1**2+(y6-sqrt(3._ark)*y5/ &
      3._ark)*y4**2*y2**2
    dF(1,260) = ((-sqrt(3._ark)*y9+sqrt(3._ark)*y8)*y5+(y8-y9)*y6)*y3*y1**2+ &
      ((sqrt(3._ark)*y9-sqrt(3._ark)*y8)*y5+(-y8+y9)*y6)*y3**2*y1+((sqrt(3._ark)*y9+ &
      sqrt(3._ark)*y8)*y5+(y8+y9)*y6)*y4*y2**2+((-sqrt(3._ark)*y9-sqrt(3._ark)*y8)*y5+(-y8- &
      y9)*y6)*y4**2*y2
    dF(1,261) = ((-y6*y7/2._ark+sqrt(3._ark)*y5*y7/2._ark)*y2-y4*y6*y7)*y1**2+((-y6*y7/ &
      2._ark+sqrt(3._ark)*y5*y7/2._ark)*y2**2-y4**2*y6*y7)*y1-y2**2*y3*y6*y7-y2*y3**2*y6*y7+ &
      (-y6*y7/2._ark+sqrt(3._ark)*y5*y7/2._ark)*y4*y3**2+(-y6*y7/2._ark+sqrt(3._ark)*y5*y7/ &
      2._ark)*y4**2*y3
    dF(1,262) = (y5*y6-sqrt(3._ark)*y6**2/3._ark)*y3*y1**2+(y5*y6-sqrt(3._ark)*y6**2/ &
      3._ark)*y3**2*y1+(sqrt(3._ark)*y6**2/3._ark-y5*y6)*y4*y2**2+(sqrt(3._ark)*y6**2/3._ark- &
      y5*y6)*y4**2*y2
    dF(1,263) = y1**2*y2*y4*y7+(y3*y4**2*y7+y2**2*y3*y7)*y1+y2*y3**2*y4*y7
    dF(1,264) = (y3**2*y4+y2*y3**2)*y1**2-y1*y2**2*y4**2-y2**2*y3*y4**2
    dF(1,265) = ((y4-y3)*y2**2+y2*y4**2-y3*y4**2)*y1**2+(-y3**2*y4**2- &
      y2**2*y3**2)*y1+y2**2*y3**2*y4+y2*y3**2*y4**2
    dF(1,266) = ((sqrt(3._ark)*y6*y9+y5*y9)*y2-2._ark*y4*y5*y8)*y1**2+((- &
      sqrt(3._ark)*y6*y9-y5*y9)*y2**2+2._ark*y4**2*y5*y8)*y1-2._ark*y2**2*y3*y5*y8+ &
      2._ark*y2*y3**2*y5*y8+(-sqrt(3._ark)*y6*y9-y5*y9)*y4*y3**2+(sqrt(3._ark)*y6*y9+ &
      y5*y9)*y4**2*y3
    dF(1,267) = (-y5*y7/2._ark+sqrt(3._ark)*y6*y7/2._ark)*y3*y1**2+(-y5*y7/2._ark+ &
      sqrt(3._ark)*y6*y7/2._ark)*y3**2*y1+(-y5*y7/2._ark+sqrt(3._ark)*y6*y7/2._ark)*y4*y2**2+(- &
      y5*y7/2._ark+sqrt(3._ark)*y6*y7/2._ark)*y4**2*y2
    dF(1,268) = ((y5*y8+sqrt(3._ark)*y6*y8)*y2-2._ark*y4*y5*y9)*y1**2+((y5*y8+ &
      sqrt(3._ark)*y6*y8)*y2**2-2._ark*y4**2*y5*y9)*y1+2._ark*y2**2*y3*y5*y9+ &
      2._ark*y2*y3**2*y5*y9+(-y5*y8-sqrt(3._ark)*y6*y8)*y4*y3**2+(-y5*y8- &
      sqrt(3._ark)*y6*y8)*y4**2*y3
    dF(1,269) = ((-2._ark*y8+y9)*y5+sqrt(3._ark)*y6*y9)*y3*y1**2+((2._ark*y8-y9)*y5- &
      sqrt(3._ark)*y6*y9)*y3**2*y1+((-y9-2._ark*y8)*y5-sqrt(3._ark)*y6*y9)*y4*y2**2+((y9+ &
      2._ark*y8)*y5+sqrt(3._ark)*y6*y9)*y4**2*y2
    dF(1,270) = (-2._ark*y2*y3*y5+(sqrt(3._ark)*y6+y5)*y4*y3)*y1**2+(2._ark*y2**2*y4*y5+ &
      ((sqrt(3._ark)*y6+y5)*y3**2+(-y5-sqrt(3._ark)*y6)*y4**2)*y2-2._ark*y3**2*y4*y5)*y1+(- &
      y5-sqrt(3._ark)*y6)*y4*y3*y2**2+2._ark*y2*y3*y4**2*y5
    dF(1,271) = y1**2*y2*y3*y4+(-y2**2*y3*y4+(-y3*y4**2+y3**2*y4)*y2)*y1
    dF(1,272) = (y6*y7-sqrt(3._ark)*y5*y7/3._ark)*y1**3+(y6*y7-sqrt(3._ark)*y5*y7/ &
      3._ark)*y2**3+(y6*y7-sqrt(3._ark)*y5*y7/3._ark)*y3**3+(y6*y7-sqrt(3._ark)*y5*y7/ &
      3._ark)*y4**3
    dF(1,273) = (-y5*y6-sqrt(3._ark)*y5**2/3._ark)*y1**3+(y5*y6+sqrt(3._ark)*y5**2/ &
      3._ark)*y2**3+(-y5*y6-sqrt(3._ark)*y5**2/3._ark)*y3**3+(y5*y6+sqrt(3._ark)*y5**2/ &
      3._ark)*y4**3
    dF(1,274) = y1**2*y3**3+y1**3*y3**2-y2**2*y4**3-y2**3*y4**2
    dF(1,275) = (y4**2+y2**2)*y1**3+(-y4**3-y2**3)*y1**2+y3**3*y4**2+y2**2*y3**3- &
      y2**3*y3**2-y3**2*y4**3
    dF(1,276) = (y2*y6+(-sqrt(3._ark)*y5/2._ark+y6/2._ark)*y4)*y1**3+(-y2**3*y6+(-y6/2._ark+ &
      sqrt(3._ark)*y5/2._ark)*y4**3)*y1+(-y6/2._ark+sqrt(3._ark)*y5/2._ark)*y3*y2**3+(- &
      sqrt(3._ark)*y5/2._ark+y6/2._ark)*y3**3*y2-y3*y4**3*y6+y3**3*y4*y6
    dF(1,277) = y1**4*y7+y4**4*y7+y3**4*y7+y2**4*y7
    dF(1,278) = (y8+y9)*y1**4+(y8-y9)*y2**4+(-y8-y9)*y3**4+(-y8+y9)*y4**4
    dF(1,279) = (sqrt(3._ark)*y6/2._ark-y5/2._ark)*y1**4+(-sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y2**4+(sqrt(3._ark)*y6/2._ark-y5/2._ark)*y3**4+(-sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y4**4
    dF(1,280) = y1**4*y3+y1*y3**4-y2*y4**4-y2**4*y4
    dF(1,281) = y1**5+y3**5-y2**5-y4**5
    dF(1,282) = y8*y9**5+y8**5*y9
    dF(1,283) = y8**3*y9**3
    dF(1,284) = (y8*y9**3+y8**3*y9)*y7**2
    dF(1,285) = y7**4*y8*y9
    dF(1,286) = -sqrt(3._ark)*y5*y7**3*y8**2+(y8**2+2._ark*y9**2)*y7**3*y6
    dF(1,287) = (sqrt(3._ark)*y9**4-sqrt(3._ark)*y8**4)*y7*y5+(y9**4-y8**4)*y7*y6
    dF(1,288) = -4._ark/3._ark*sqrt(3._ark)*y5*y6**2*y7**3+y5**2*y6*y7**3+y6**3*y7**3
    dF(1,289) = -y5*y7*y8**2*y9**2/2._ark+sqrt(3._ark)*y6*y7*y8**2*y9**2/2._ark
    dF(1,290) = (-2._ark*y8**2+y9**2)*y7**3*y5+sqrt(3._ark)*y6*y7**3*y9**2
    dF(1,291) = (y9**4-2._ark*y8**4)*y7*y5+sqrt(3._ark)*y6*y7*y9**4
    dF(1,292) = y5*y7**5-sqrt(3._ark)*y6*y7**5
    dF(1,293) = (-sqrt(3._ark)*y8**3*y9/2._ark-sqrt(3._ark)*y8*y9**3/2._ark)*y5**2+(- &
      y8**3*y9-y8*y9**3)*y6*y5+(-sqrt(3._ark)*y8*y9**3/6._ark-sqrt(3._ark)*y8**3*y9/ &
      6._ark)*y6**2
    dF(1,294) = y5*y6*y7**2*y8*y9-sqrt(3._ark)*y6**2*y7**2*y8*y9/3._ark
    dF(1,295) = -sqrt(3._ark)*y5**2*y8*y9**3/2._ark-y5*y6*y8**3*y9+(- &
      sqrt(3._ark)*y8**3*y9/3._ark+sqrt(3._ark)*y8*y9**3/6._ark)*y6**2
    dF(1,296) = (-y9**2/3._ark-4._ark/3._ark*y8**2)*y7*y5**3+sqrt(3._ark)*y5**2*y6*y7*y9**2+ &
      y5*y6**2*y7*y9**2+(5._ark/9._ark*sqrt(3._ark)*y9**2+4._ark/ &
      9._ark*sqrt(3._ark)*y8**2)*y7*y6**3
    dF(1,297) = y5**3*y7*y8**2-sqrt(3._ark)*y5**2*y6*y7*y9**2+y5*y6**2*y7*y8**2+(- &
      5._ark/9._ark*sqrt(3._ark)*y9**2-4._ark/9._ark*sqrt(3._ark)*y8**2)*y7*y6**3
    dF(1,298) = y6**2*y7**2*y8*y9+y5**2*y7**2*y8*y9
    dF(1,299) = (y8*y9**3+y8**3*y9)*y5**2+(y8*y9**3+y8**3*y9)*y6**2
    dF(1,300) = (4._ark/9._ark*sqrt(3._ark)*y8**2+4._ark/9._ark*sqrt(3._ark)*y9**2)*y7*y5**3- &
      2._ark*y5**2*y6*y7*y9**2+(-2._ark/3._ark*y9**2-4._ark/3._ark*y8**2)*y7*y6**3
    dF(1,301) = (4._ark/9._ark*sqrt(3._ark)*y8**2+4._ark/9._ark*sqrt(3._ark)*y9**2)*y7*y5**3+(- &
      y9**2-y8**2)*y7*y6*y5**2+(-y9**2-y8**2)*y7*y6**3
    dF(1,302) = 5._ark/24._ark*y6**4*y8*y9+3._ark/8._ark*y5**4*y8*y9+ &
      sqrt(3._ark)*y5*y6**3*y8*y9/3._ark+y5**2*y6**2*y8*y9/4._ark
    dF(1,303) = y5**2*y6**3*y7+17._ark/6._ark*sqrt(3._ark)*y5*y6**4*y7+2._ark*y6**5*y7+ &
      3._ark*y5**4*y6*y7-3._ark/2._ark*sqrt(3._ark)*y5**5*y7
    dF(1,304) = y5**3*y7**3-3._ark*y5*y6**2*y7**3
    dF(1,305) = -sqrt(3._ark)*y6**4*y8*y9/36._ark-sqrt(3._ark)*y5**4*y8*y9/4._ark- &
      y5**3*y6*y8*y9-y5*y6**3*y8*y9/3._ark-sqrt(3._ark)*y5**2*y6**2*y8*y9/2._ark
    dF(1,306) = -6._ark*y5*y6**4*y7+y5**3*y6**2*y7-3._ark/2._ark*sqrt(3._ark)*y6**5*y7-5._ark/ &
      2._ark*sqrt(3._ark)*y5**4*y6*y7+3._ark*y5**5*y7
    dF(1,307) = 7._ark/12._ark*y6**4*y8*y9+y5**4*y8*y9/4._ark-2._ark/ &
      3._ark*sqrt(3._ark)*y5*y6**3*y8*y9+3._ark/2._ark*y5**2*y6**2*y8*y9
    dF(1,308) = -15._ark*y5*y6**4*y7-3._ark*sqrt(3._ark)*y6**5*y7- &
      5._ark*sqrt(3._ark)*y5**4*y6*y7+7._ark*y5**5*y7
    dF(1,309) = (-y8-y9)*y7**4*y1+(-y8+y9)*y7**4*y2+(y8+y9)*y7**4*y3+(y8- &
      y9)*y7**4*y4
    dF(1,310) = ((sqrt(3._ark)*y8**2-sqrt(3._ark)*y9**2)*y7**2*y5+(y8**2- &
      y9**2)*y7**2*y6)*y1+((-sqrt(3._ark)*y8**2+sqrt(3._ark)*y9**2)*y7**2*y5+(y9**2- &
      y8**2)*y7**2*y6)*y2+((sqrt(3._ark)*y8**2-sqrt(3._ark)*y9**2)*y7**2*y5+(y8**2- &
      y9**2)*y7**2*y6)*y3+((-sqrt(3._ark)*y8**2+sqrt(3._ark)*y9**2)*y7**2*y5+(y9**2- &
      y8**2)*y7**2*y6)*y4
    dF(1,311) = (sqrt(3._ark)*y5*y8*y9**3/2._ark+(-y8*y9**3/2._ark-y8**3*y9)*y6)*y1+ &
      (sqrt(3._ark)*y5*y8*y9**3/2._ark+(-y8*y9**3/2._ark-y8**3*y9)*y6)*y2+ &
      (sqrt(3._ark)*y5*y8*y9**3/2._ark+(-y8*y9**3/2._ark-y8**3*y9)*y6)*y3+ &
      (sqrt(3._ark)*y5*y8*y9**3/2._ark+(-y8*y9**3/2._ark-y8**3*y9)*y6)*y4
    dF(1,312) = ((-y8/2._ark-y9/2._ark)*y7**3*y5+(sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/ &
      2._ark)*y7**3*y6)*y1+((y8/2._ark-y9/2._ark)*y7**3*y5+(-sqrt(3._ark)*y8/2._ark+ &
      sqrt(3._ark)*y9/2._ark)*y7**3*y6)*y2+((y8/2._ark+y9/2._ark)*y7**3*y5+(-sqrt(3._ark)*y8/ &
      2._ark-sqrt(3._ark)*y9/2._ark)*y7**3*y6)*y3+((y9/2._ark-y8/2._ark)*y7**3*y5+(- &
      sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/2._ark)*y7**3*y6)*y4
    dF(1,313) = ((y8**3*y9-y8*y9**3/2._ark)*y5-sqrt(3._ark)*y6*y8*y9**3/2._ark)*y1+ &
      ((y8**3*y9-y8*y9**3/2._ark)*y5-sqrt(3._ark)*y6*y8*y9**3/2._ark)*y2+((y8**3*y9- &
      y8*y9**3/2._ark)*y5-sqrt(3._ark)*y6*y8*y9**3/2._ark)*y3+((y8**3*y9-y8*y9**3/2._ark)*y5- &
      sqrt(3._ark)*y6*y8*y9**3/2._ark)*y4
    dF(1,314) = ((-sqrt(3._ark)*y8*y9**2/6._ark-2._ark/3._ark*sqrt(3._ark)*y8**2*y9)*y5**2- &
      y5*y6*y8*y9**2-sqrt(3._ark)*y6**2*y8*y9**2/2._ark)*y1+((-sqrt(3._ark)*y8*y9**2/6._ark+ &
      2._ark/3._ark*sqrt(3._ark)*y8**2*y9)*y5**2-y5*y6*y8*y9**2-sqrt(3._ark)*y6**2*y8*y9**2/ &
      2._ark)*y2+((2._ark/3._ark*sqrt(3._ark)*y8**2*y9+sqrt(3._ark)*y8*y9**2/6._ark)*y5**2+ &
      y5*y6*y8*y9**2+sqrt(3._ark)*y6**2*y8*y9**2/2._ark)*y3+((-2._ark/ &
      3._ark*sqrt(3._ark)*y8**2*y9+sqrt(3._ark)*y8*y9**2/6._ark)*y5**2+y5*y6*y8*y9**2+ &
      sqrt(3._ark)*y6**2*y8*y9**2/2._ark)*y4
    dF(1,315) = ((-y8-y9)*y7**2*y5**2+(-y8-y9)*y7**2*y6**2)*y1+((-y8+ &
      y9)*y7**2*y5**2+(-y8+y9)*y7**2*y6**2)*y2+((y8+y9)*y7**2*y5**2+(y8+ &
      y9)*y7**2*y6**2)*y3+((y8-y9)*y7**2*y5**2+(y8-y9)*y7**2*y6**2)*y4
    dF(1,316) = ((y8**2+y9**2)*y7*y5**2+(y8**2+y9**2)*y7*y6**2)*y1+((y8**2+ &
      y9**2)*y7*y5**2+(y8**2+y9**2)*y7*y6**2)*y2+((y8**2+y9**2)*y7*y5**2+(y8**2+ &
      y9**2)*y7*y6**2)*y3+((y8**2+y9**2)*y7*y5**2+(y8**2+y9**2)*y7*y6**2)*y4
    dF(1,317) = (6._ark*y5**4*y8+4._ark*sqrt(3._ark)*y5**3*y6*y9+(10._ark*y8+ &
      10._ark*y9)*y6**2*y5**2-4._ark*sqrt(3._ark)*y5*y6**3*y8+6._ark*y6**4*y9)*y1+ &
      (6._ark*y5**4*y8-4._ark*sqrt(3._ark)*y5**3*y6*y9+(10._ark*y8-10._ark*y9)*y6**2*y5**2- &
      4._ark*sqrt(3._ark)*y5*y6**3*y8-6._ark*y6**4*y9)*y2+(-6._ark*y5**4*y8- &
      4._ark*sqrt(3._ark)*y5**3*y6*y9+(-10._ark*y8-10._ark*y9)*y6**2*y5**2+ &
      4._ark*sqrt(3._ark)*y5*y6**3*y8-6._ark*y6**4*y9)*y3+(-6._ark*y5**4*y8+ &
      4._ark*sqrt(3._ark)*y5**3*y6*y9+(-10._ark*y8+10._ark*y9)*y6**2*y5**2+ &
      4._ark*sqrt(3._ark)*y5*y6**3*y8+6._ark*y6**4*y9)*y4
    dF(1,318) = y4**2*y7**2*y8*y9+y3**2*y7**2*y8*y9+y2**2*y7**2*y8*y9+ &
      y1**2*y7**2*y8*y9
    dF(1,319) = ((sqrt(3._ark)*y8**2*y9-sqrt(3._ark)*y8*y9**2)*y5+(-y8*y9**2+ &
      y8**2*y9)*y6)*y1**2+((-sqrt(3._ark)*y8**2*y9-sqrt(3._ark)*y8*y9**2)*y5+(-y8*y9**2- &
      y8**2*y9)*y6)*y2**2+((sqrt(3._ark)*y8*y9**2-sqrt(3._ark)*y8**2*y9)*y5+(y8*y9**2- &
      y8**2*y9)*y6)*y3**2+((sqrt(3._ark)*y8*y9**2+sqrt(3._ark)*y8**2*y9)*y5+(y8**2*y9+ &
      y8*y9**2)*y6)*y4**2
    dF(1,320) = (4._ark/3._ark*sqrt(3._ark)*y5*y6**2*y7-y5**2*y6*y7-y6**3*y7)*y1**2+(4._ark/ &
      3._ark*sqrt(3._ark)*y5*y6**2*y7-y5**2*y6*y7-y6**3*y7)*y2**2+(4._ark/ &
      3._ark*sqrt(3._ark)*y5*y6**2*y7-y5**2*y6*y7-y6**3*y7)*y3**2+(4._ark/ &
      3._ark*sqrt(3._ark)*y5*y6**2*y7-y5**2*y6*y7-y6**3*y7)*y4**2
    dF(1,321) = (y5**3*y7-3._ark*y5*y6**2*y7)*y1**2+(y5**3*y7-3._ark*y5*y6**2*y7)*y2**2+ &
      (y5**3*y7-3._ark*y5*y6**2*y7)*y3**2+(y5**3*y7-3._ark*y5*y6**2*y7)*y4**2
    dF(1,322) = (y6**2*y7+y5**2*y7)*y1**3+(y6**2*y7+y5**2*y7)*y2**3+(y6**2*y7+ &
      y5**2*y7)*y3**3+(y6**2*y7+y5**2*y7)*y4**3
    dF(1,323) = ((2._ark*y8+2._ark*y9)*y5**2+(2._ark*sqrt(3._ark)*y8+ &
      2._ark*sqrt(3._ark)*y9)*y6*y5)*y1**3+((-2._ark*y9+2._ark*y8)*y5**2+(2._ark*sqrt(3._ark)*y8- &
      2._ark*sqrt(3._ark)*y9)*y6*y5)*y2**3+((-2._ark*y8-2._ark*y9)*y5**2+(-2._ark*sqrt(3._ark)*y9- &
      2._ark*sqrt(3._ark)*y8)*y6*y5)*y3**3+((-2._ark*y8+2._ark*y9)*y5**2+(-2._ark*sqrt(3._ark)*y8+ &
      2._ark*sqrt(3._ark)*y9)*y6*y5)*y4**3
    dF(1,324) = ((y9/2._ark-y8)*y5+sqrt(3._ark)*y6*y9/2._ark)*y1**4+((-y9/2._ark-y8)*y5- &
      sqrt(3._ark)*y6*y9/2._ark)*y2**4+((y8-y9/2._ark)*y5-sqrt(3._ark)*y6*y9/2._ark)*y3**4+((y8+ &
      y9/2._ark)*y5+sqrt(3._ark)*y6*y9/2._ark)*y4**4
    dF(1,325) = (-sqrt(3._ark)*y6**2/6._ark-y5*y6-sqrt(3._ark)*y5**2/2._ark)*y1**4+(y5*y6+ &
      sqrt(3._ark)*y5**2/2._ark+sqrt(3._ark)*y6**2/6._ark)*y2**4+(-sqrt(3._ark)*y6**2/6._ark- &
      y5*y6-sqrt(3._ark)*y5**2/2._ark)*y3**4+(y5*y6+sqrt(3._ark)*y5**2/2._ark+ &
      sqrt(3._ark)*y6**2/6._ark)*y4**4
    dF(1,326) = (-y6**2-y5**2)*y1**4+(y6**2+y5**2)*y2**4+(-y6**2-y5**2)*y3**4+ &
      (y6**2+y5**2)*y4**4
    dF(1,327) = y4**5*y7+y3**5*y7+y1**5*y7+y2**5*y7
    dF(1,328) = (-y8-y9)*y1**5+(-y8+y9)*y2**5+(y8+y9)*y3**5+(y8-y9)*y4**5
    dF(1,329) = (-y8**5-y9**5)*y1+(y9**5-y8**5)*y2+(y8**5+y9**5)*y3+(y8**5-y9**5)*y4
    dF(1,330) = (-y9**3-y8**3)*y7**2*y1+(y9**3-y8**3)*y7**2*y2+(y8**3+ &
      y9**3)*y7**2*y3+(-y9**3+y8**3)*y7**2*y4
    dF(1,331) = (y8**2*y9+y8*y9**2)*y7**2*y1+(y8*y9**2-y8**2*y9)*y7**2*y2+(- &
      y8*y9**2-y8**2*y9)*y7**2*y3+(-y8*y9**2+y8**2*y9)*y7**2*y4
    dF(1,332) = (-y6*y7**2*y8*y9+sqrt(3._ark)*y5*y7**2*y8*y9/3._ark)*y1+(- &
      y6*y7**2*y8*y9+sqrt(3._ark)*y5*y7**2*y8*y9/3._ark)*y2+(-y6*y7**2*y8*y9+ &
      sqrt(3._ark)*y5*y7**2*y8*y9/3._ark)*y3+(-y6*y7**2*y8*y9+sqrt(3._ark)*y5*y7**2*y8*y9/ &
      3._ark)*y4
    dF(1,333) = (-3._ark*sqrt(3._ark)*y5**3*y7*y8+9._ark*y5**2*y6*y7*y8- &
      3._ark*sqrt(3._ark)*y5*y6**2*y7*y8+(8._ark*y9+y8)*y7*y6**3)*y1+ &
      (3._ark*sqrt(3._ark)*y5**3*y7*y8-9._ark*y5**2*y6*y7*y8+3._ark*sqrt(3._ark)*y5*y6**2*y7*y8+ &
      (-y8+8._ark*y9)*y7*y6**3)*y2+(3._ark*sqrt(3._ark)*y5**3*y7*y8-9._ark*y5**2*y6*y7*y8+ &
      3._ark*sqrt(3._ark)*y5*y6**2*y7*y8+(-8._ark*y9-y8)*y7*y6**3)*y3+(- &
      3._ark*sqrt(3._ark)*y5**3*y7*y8+9._ark*y5**2*y6*y7*y8-3._ark*sqrt(3._ark)*y5*y6**2*y7*y8+ &
      (-8._ark*y9+y8)*y7*y6**3)*y4
    dF(1,334) = ((y8**4/2._ark+y9**4/2._ark)*y5+(-sqrt(3._ark)*y9**4/2._ark- &
      sqrt(3._ark)*y8**4/2._ark)*y6)*y1+((-y8**4/2._ark-y9**4/2._ark)*y5+(sqrt(3._ark)*y8**4/ &
      2._ark+sqrt(3._ark)*y9**4/2._ark)*y6)*y2+((y8**4/2._ark+y9**4/2._ark)*y5+(- &
      sqrt(3._ark)*y9**4/2._ark-sqrt(3._ark)*y8**4/2._ark)*y6)*y3+((-y8**4/2._ark-y9**4/ &
      2._ark)*y5+(sqrt(3._ark)*y8**4/2._ark+sqrt(3._ark)*y9**4/2._ark)*y6)*y4
    dF(1,335) = (-sqrt(3._ark)*y5**2*y8**3/2._ark-y5*y6*y9**3+(-sqrt(3._ark)*y9**3/3._ark+ &
      sqrt(3._ark)*y8**3/6._ark)*y6**2)*y1+(-sqrt(3._ark)*y5**2*y8**3/2._ark+y5*y6*y9**3+ &
      (sqrt(3._ark)*y8**3/6._ark+sqrt(3._ark)*y9**3/3._ark)*y6**2)*y2+(sqrt(3._ark)*y5**2*y8**3/ &
      2._ark+y5*y6*y9**3+(-sqrt(3._ark)*y8**3/6._ark+sqrt(3._ark)*y9**3/3._ark)*y6**2)*y3+ &
      (sqrt(3._ark)*y5**2*y8**3/2._ark-y5*y6*y9**3+(-sqrt(3._ark)*y9**3/3._ark- &
      sqrt(3._ark)*y8**3/6._ark)*y6**2)*y4
    dF(1,336) = (sqrt(3._ark)*y5**2*y7**2*y8/2._ark-y5*y6*y7**2*y8+(sqrt(3._ark)*y8/6._ark+ &
      2._ark/3._ark*sqrt(3._ark)*y9)*y7**2*y6**2)*y1+(sqrt(3._ark)*y5**2*y7**2*y8/2._ark- &
      y5*y6*y7**2*y8+(-2._ark/3._ark*sqrt(3._ark)*y9+sqrt(3._ark)*y8/6._ark)*y7**2*y6**2)*y2+(- &
      sqrt(3._ark)*y5**2*y7**2*y8/2._ark+y5*y6*y7**2*y8+(-sqrt(3._ark)*y8/6._ark-2._ark/ &
      3._ark*sqrt(3._ark)*y9)*y7**2*y6**2)*y3+(-sqrt(3._ark)*y5**2*y7**2*y8/2._ark+ &
      y5*y6*y7**2*y8+(-sqrt(3._ark)*y8/6._ark+2._ark/3._ark*sqrt(3._ark)*y9)*y7**2*y6**2)*y4
    dF(1,337) = (-y5*y6**2*y7**2+y5**3*y7**2/3._ark)*y1+(y5*y6**2*y7**2-y5**3*y7**2/ &
      3._ark)*y2+(-y5*y6**2*y7**2+y5**3*y7**2/3._ark)*y3+(y5*y6**2*y7**2-y5**3*y7**2/ &
      3._ark)*y4
    dF(1,338) = (y5*y6**2*y8*y9-y5**3*y8*y9/3._ark)*y1+(y5*y6**2*y8*y9-y5**3*y8*y9/ &
      3._ark)*y2+(y5*y6**2*y8*y9-y5**3*y8*y9/3._ark)*y3+(y5*y6**2*y8*y9-y5**3*y8*y9/ &
      3._ark)*y4
    dF(1,339) = (12._ark/5._ark*sqrt(3._ark)*y5**4*y8+27._ark/5._ark*y5**3*y6*y9+(24._ark/ &
      5._ark*sqrt(3._ark)*y8+21._ark/5._ark*sqrt(3._ark)*y9)*y6**2*y5**2+(-28._ark/5._ark*y8- &
      y9)*y6**3*y5+13._ark/5._ark*sqrt(3._ark)*y6**4*y9)*y1+(12._ark/5._ark*sqrt(3._ark)*y5**4*y8- &
      27._ark/5._ark*y5**3*y6*y9+(-21._ark/5._ark*sqrt(3._ark)*y9+24._ark/ &
      5._ark*sqrt(3._ark)*y8)*y6**2*y5**2+(y9-28._ark/5._ark*y8)*y6**3*y5-13._ark/ &
      5._ark*sqrt(3._ark)*y6**4*y9)*y2+(-12._ark/5._ark*sqrt(3._ark)*y5**4*y8-27._ark/ &
      5._ark*y5**3*y6*y9+(-24._ark/5._ark*sqrt(3._ark)*y8-21._ark/ &
      5._ark*sqrt(3._ark)*y9)*y6**2*y5**2+(y9+28._ark/5._ark*y8)*y6**3*y5-13._ark/ &
      5._ark*sqrt(3._ark)*y6**4*y9)*y3+(-12._ark/5._ark*sqrt(3._ark)*y5**4*y8+27._ark/ &
      5._ark*y5**3*y6*y9+(21._ark/5._ark*sqrt(3._ark)*y9-24._ark/ &
      5._ark*sqrt(3._ark)*y8)*y6**2*y5**2+(-y9+28._ark/5._ark*y8)*y6**3*y5+13._ark/ &
      5._ark*sqrt(3._ark)*y6**4*y9)*y4
    dF(1,340) = (-sqrt(3._ark)*y5**4*y7/6._ark-y5*y6**3*y7+sqrt(3._ark)*y5**2*y6**2*y7/ &
      6._ark)*y1+(-sqrt(3._ark)*y5**4*y7/6._ark-y5*y6**3*y7+sqrt(3._ark)*y5**2*y6**2*y7/ &
      6._ark)*y2+(-sqrt(3._ark)*y5**4*y7/6._ark-y5*y6**3*y7+sqrt(3._ark)*y5**2*y6**2*y7/ &
      6._ark)*y3+(-sqrt(3._ark)*y5**4*y7/6._ark-y5*y6**3*y7+sqrt(3._ark)*y5**2*y6**2*y7/ &
      6._ark)*y4
    dF(1,341) = (3._ark/2._ark*sqrt(3._ark)*y6**5+6._ark*y5*y6**4+5._ark/ &
      2._ark*sqrt(3._ark)*y5**4*y6-3._ark*y5**5-y5**3*y6**2)*y1+(-6._ark*y5*y6**4+3._ark*y5**5- &
      3._ark/2._ark*sqrt(3._ark)*y6**5-5._ark/2._ark*sqrt(3._ark)*y5**4*y6+y5**3*y6**2)*y2+(3._ark/ &
      2._ark*sqrt(3._ark)*y6**5+6._ark*y5*y6**4+5._ark/2._ark*sqrt(3._ark)*y5**4*y6-3._ark*y5**5- &
      y5**3*y6**2)*y3+(-6._ark*y5*y6**4+3._ark*y5**5-3._ark/2._ark*sqrt(3._ark)*y6**5-5._ark/ &
      2._ark*sqrt(3._ark)*y5**4*y6+y5**3*y6**2)*y4
    dF(1,342) = y2*y4*y7**4-y1*y3*y7**4
    dF(1,343) = ((-sqrt(3._ark)*y5*y7*y9**2+y6*y7*y9**2)*y2+2._ark*y4*y6*y7*y8**2)*y1+ &
      2._ark*y2*y3*y6*y7*y8**2+(-sqrt(3._ark)*y5*y7*y9**2+y6*y7*y9**2)*y4*y3
    dF(1,344) = (y2*y6*y8*y9**2+(-sqrt(3._ark)*y5*y8**2*y9/2._ark+y6*y8**2*y9/ &
      2._ark)*y4)*y1+(-y6*y8**2*y9/2._ark+sqrt(3._ark)*y5*y8**2*y9/2._ark)*y3*y2- &
      y3*y4*y6*y8*y9**2
    dF(1,345) = (y2*y6*y7**3+(-sqrt(3._ark)*y5*y7**3/2._ark+y6*y7**3/2._ark)*y4)*y1+(- &
      sqrt(3._ark)*y5*y7**3/2._ark+y6*y7**3/2._ark)*y3*y2+y3*y4*y6*y7**3
    dF(1,346) = (y2*y6*y8**3+(-sqrt(3._ark)*y5*y9**3/2._ark+y6*y9**3/2._ark)*y4)*y1+ &
      (sqrt(3._ark)*y5*y9**3/2._ark-y6*y9**3/2._ark)*y3*y2-y3*y4*y6*y8**3
    dF(1,347) = ((-y6**2*y7*y9-y5**2*y7*y9)*y2+(-y6**2*y7*y8-y5**2*y7*y8)*y4)*y1+ &
      (y5**2*y7*y8+y6**2*y7*y8)*y3*y2+(y6**2*y7*y9+y5**2*y7*y9)*y4*y3
    dF(1,348) = (y6**2*y8*y9+y5**2*y8*y9)*y3*y1+(y6**2*y8*y9+y5**2*y8*y9)*y4*y2
    dF(1,349) = (-y2*y5*y7**2*y8+(y5*y7**2*y9/2._ark+sqrt(3._ark)*y6*y7**2*y9/ &
      2._ark)*y4)*y1+(-sqrt(3._ark)*y6*y7**2*y9/2._ark-y5*y7**2*y9/2._ark)*y3*y2+ &
      y3*y4*y5*y7**2*y8
    dF(1,350) = ((y5*y7*y8**2-sqrt(3._ark)*y6*y7*y8**2)*y2+(y5*y7*y9**2- &
      sqrt(3._ark)*y6*y7*y9**2)*y4)*y1+(y5*y7*y9**2-sqrt(3._ark)*y6*y7*y9**2)*y3*y2+ &
      (y5*y7*y8**2-sqrt(3._ark)*y6*y7*y8**2)*y4*y3
    dF(1,351) = (y2*y5*y7**3+(-y5*y7**3/2._ark-sqrt(3._ark)*y6*y7**3/2._ark)*y4)*y1+(- &
      y5*y7**3/2._ark-sqrt(3._ark)*y6*y7**3/2._ark)*y3*y2+y3*y4*y5*y7**3
    dF(1,352) = (-y2*y5*y8**3+(y5*y9**3/2._ark+sqrt(3._ark)*y6*y9**3/2._ark)*y4)*y1+(- &
      sqrt(3._ark)*y6*y9**3/2._ark-y5*y9**3/2._ark)*y3*y2+y3*y4*y5*y8**3
    dF(1,353) = (sqrt(3._ark)*y5**2*y8**2/2._ark+y5*y6*y9**2+(-sqrt(3._ark)*y8**2/6._ark+ &
      sqrt(3._ark)*y9**2/3._ark)*y6**2)*y3*y1+(-sqrt(3._ark)*y5**2*y8**2/2._ark-y5*y6*y9**2+ &
      (sqrt(3._ark)*y8**2/6._ark-sqrt(3._ark)*y9**2/3._ark)*y6**2)*y4*y2
    dF(1,354) = ((sqrt(3._ark)*y5**2*y7*y9/3._ark-y5*y6*y7*y9)*y2+ &
      (sqrt(3._ark)*y6**2*y7*y8/2._ark-sqrt(3._ark)*y5**2*y7*y8/6._ark)*y4)*y1+(- &
      sqrt(3._ark)*y6**2*y7*y8/2._ark+sqrt(3._ark)*y5**2*y7*y8/6._ark)*y3*y2+(- &
      sqrt(3._ark)*y5**2*y7*y9/3._ark+y5*y6*y7*y9)*y4*y3
    dF(1,355) = (3._ark/8._ark*y6**4-y5**4/8._ark-sqrt(3._ark)*y5*y6**3+5._ark/ &
      4._ark*y5**2*y6**2)*y3*y1+(y5**4/8._ark-5._ark/4._ark*y5**2*y6**2-3._ark/8._ark*y6**4+ &
      sqrt(3._ark)*y5*y6**3)*y4*y2
    dF(1,356) = ((-y6*y8*y9-sqrt(3._ark)*y5*y8*y9)*y2+(sqrt(3._ark)*y5*y8*y9+ &
      y6*y8*y9)*y4)*y1**2+((-y6*y8*y9-sqrt(3._ark)*y5*y8*y9)*y2**2+(sqrt(3._ark)*y5*y8*y9+ &
      y6*y8*y9)*y4**2)*y1+(sqrt(3._ark)*y5*y8*y9+y6*y8*y9)*y3*y2**2+ &
      (sqrt(3._ark)*y5*y8*y9+y6*y8*y9)*y3**2*y2+(-y6*y8*y9- &
      sqrt(3._ark)*y5*y8*y9)*y4*y3**2+(-y6*y8*y9-sqrt(3._ark)*y5*y8*y9)*y4**2*y3
    dF(1,357) = ((-5._ark/9._ark*sqrt(3._ark)*y5*y6**2-y6**3-sqrt(3._ark)*y5**3/9._ark+11._ark/ &
      9._ark*y5**2*y6)*y2+(4._ark/9._ark*sqrt(3._ark)*y5**3-14._ark/9._ark*y5**2*y6+2._ark/ &
      3._ark*y6**3)*y4)*y1**2+((y6**3-11._ark/9._ark*y5**2*y6+5._ark/9._ark*sqrt(3._ark)*y5*y6**2+ &
      sqrt(3._ark)*y5**3/9._ark)*y2**2+(14._ark/9._ark*y5**2*y6-4._ark/9._ark*sqrt(3._ark)*y5**3- &
      2._ark/3._ark*y6**3)*y4**2)*y1+(14._ark/9._ark*y5**2*y6-4._ark/9._ark*sqrt(3._ark)*y5**3-2._ark/ &
      3._ark*y6**3)*y3*y2**2+(4._ark/9._ark*sqrt(3._ark)*y5**3-14._ark/9._ark*y5**2*y6+2._ark/ &
      3._ark*y6**3)*y3**2*y2+(-5._ark/9._ark*sqrt(3._ark)*y5*y6**2-y6**3-sqrt(3._ark)*y5**3/ &
      9._ark+11._ark/9._ark*y5**2*y6)*y4*y3**2+(y6**3-11._ark/9._ark*y5**2*y6+5._ark/ &
      9._ark*sqrt(3._ark)*y5*y6**2+sqrt(3._ark)*y5**3/9._ark)*y4**2*y3
    dF(1,358) = y1**2*y8**2*y9**2-y4**2*y8**2*y9**2-y2**2*y8**2*y9**2+ &
      y3**2*y8**2*y9**2
    dF(1,359) = (-y9**2-y8**2)*y7**2*y1**2+(y8**2+y9**2)*y7**2*y2**2+(-y9**2- &
      y8**2)*y7**2*y3**2+(y8**2+y9**2)*y7**2*y4**2
    dF(1,360) = (-y8*y9**2-y8**2*y9)*y7*y1**2+(y8*y9**2-y8**2*y9)*y7*y2**2+ &
      (y8**2*y9+y8*y9**2)*y7*y3**2+(-y8*y9**2+y8**2*y9)*y7*y4**2
    dF(1,361) = (y8*y9**3+y8**3*y9)*y1**2+(y8*y9**3+y8**3*y9)*y2**2+(y8*y9**3+ &
      y8**3*y9)*y3**2+(y8*y9**3+y8**3*y9)*y4**2
    dF(1,362) = (y8**3+y9**3)*y7*y1**2+(y9**3-y8**3)*y7*y2**2+(-y9**3- &
      y8**3)*y7*y3**2+(-y9**3+y8**3)*y7*y4**2
    dF(1,363) = (-y5**2*y6**2-y6**4/2._ark-y5**4/2._ark)*y1**2+(y5**4/2._ark+y5**2*y6**2+ &
      y6**4/2._ark)*y2**2+(-y5**2*y6**2-y6**4/2._ark-y5**4/2._ark)*y3**2+(y5**4/2._ark+ &
      y5**2*y6**2+y6**4/2._ark)*y4**2
    dF(1,364) = ((-y8-y9)*y5**3+(3._ark*sqrt(3._ark)*y8+3._ark*sqrt(3._ark)*y9)*y6*y5**2+(- &
      9._ark*y9-9._ark*y8)*y6**2*y5+(3._ark*sqrt(3._ark)*y8+3._ark*sqrt(3._ark)*y9)*y6**3)*y1**2+ &
      ((-y8+y9)*y5**3+(3._ark*sqrt(3._ark)*y8-3._ark*sqrt(3._ark)*y9)*y6*y5**2+(-9._ark*y8+ &
      9._ark*y9)*y6**2*y5+(3._ark*sqrt(3._ark)*y8-3._ark*sqrt(3._ark)*y9)*y6**3)*y2**2+((y8+ &
      y9)*y5**3+(-3._ark*sqrt(3._ark)*y8-3._ark*sqrt(3._ark)*y9)*y6*y5**2+(9._ark*y8+ &
      9._ark*y9)*y6**2*y5+(-3._ark*sqrt(3._ark)*y8-3._ark*sqrt(3._ark)*y9)*y6**3)*y3**2+((y8- &
      y9)*y5**3+(3._ark*sqrt(3._ark)*y9-3._ark*sqrt(3._ark)*y8)*y6*y5**2+(-9._ark*y9+ &
      9._ark*y8)*y6**2*y5+(3._ark*sqrt(3._ark)*y9-3._ark*sqrt(3._ark)*y8)*y6**3)*y4**2
    dF(1,365) = (-y4*y7**2*y8-y2*y7**2*y9)*y1**2+(y4**2*y7**2*y8+y2**2*y7**2*y9)*y1+ &
      y2*y3**2*y7**2*y8-y3*y4**2*y7**2*y9-y2**2*y3*y7**2*y8+y3**2*y4*y7**2*y9
    dF(1,366) = (y4*y8*y9**2+y2*y8**2*y9)*y1**2+(-y2**2*y8**2*y9-y4**2*y8*y9**2)*y1+ &
      y2**2*y3*y8*y9**2-y3**2*y4*y8**2*y9+y3*y4**2*y8**2*y9-y2*y3**2*y8*y9**2
    dF(1,367) = (-y2**2*y7*y9-y4**2*y7*y8)*y1**2+y3**2*y4**2*y7*y9+y2**2*y3**2*y7*y8
    dF(1,368) = y2**2*y4**2*y8*y9+y1**2*y3**2*y8*y9
    dF(1,369) = (-y9**2-y8**2)*y3**2*y1**2+(y8**2+y9**2)*y4**2*y2**2
    dF(1,370) = (y2**2*y6*y7+(y6*y7/2._ark-sqrt(3._ark)*y5*y7/2._ark)*y4**2)*y1**2+(y6*y7/ &
      2._ark-sqrt(3._ark)*y5*y7/2._ark)*y3**2*y2**2+y3**2*y4**2*y6*y7
    dF(1,371) = ((-y5*y8-sqrt(3._ark)*y6*y8)*y2**2+2._ark*y4**2*y5*y9)*y1**2- &
      2._ark*y9*y5*y3**2*y2**2+(y5*y8+sqrt(3._ark)*y6*y8)*y4**2*y3**2
    dF(1,372) = (y2**2*y5*y7+(-sqrt(3._ark)*y6*y7/2._ark-y5*y7/2._ark)*y4**2)*y1**2+(- &
      sqrt(3._ark)*y6*y7/2._ark-y5*y7/2._ark)*y3**2*y2**2+y3**2*y4**2*y5*y7
    dF(1,373) = (-sqrt(3._ark)*y6/2._ark+y5/2._ark)*y3**2*y1**3+(-sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y3**3*y1**2+(sqrt(3._ark)*y6/2._ark-y5/2._ark)*y4**2*y2**3+(sqrt(3._ark)*y6/2._ark- &
      y5/2._ark)*y4**3*y2**2
    dF(1,374) = (-y4**2-y2**2)*y1**4+(y4**4+y2**4)*y1**2+y2**4*y3**2+y3**2*y4**4- &
      y3**4*y4**2-y2**2*y3**4
    dF(1,375) = ((-3._ark*y8-3._ark*y9)*y5**2+(-2._ark*sqrt(3._ark)*y9- &
      2._ark*sqrt(3._ark)*y8)*y6*y5+(-y8-y9)*y6**2)*y1**3+((-3._ark*y8+3._ark*y9)*y5**2+(- &
      2._ark*sqrt(3._ark)*y8+2._ark*sqrt(3._ark)*y9)*y6*y5+(-y8+y9)*y6**2)*y2**3+((3._ark*y8+ &
      3._ark*y9)*y5**2+(2._ark*sqrt(3._ark)*y8+2._ark*sqrt(3._ark)*y9)*y6*y5+(y8+ &
      y9)*y6**2)*y3**3+((3._ark*y8-3._ark*y9)*y5**2+(2._ark*sqrt(3._ark)*y8- &
      2._ark*sqrt(3._ark)*y9)*y6*y5+(y8-y9)*y6**2)*y4**3
    dF(1,376) = (-y5**3+3._ark*y5*y6**2)*y1**3+(y5**3-3._ark*y5*y6**2)*y2**3+(-y5**3+ &
      3._ark*y5*y6**2)*y3**3+(y5**3-3._ark*y5*y6**2)*y4**3
    dF(1,377) = (y2*y8*y9+y4*y8*y9)*y1**3+(y2**3*y8*y9+y4**3*y8*y9)*y1+ &
      y3**3*y4*y8*y9+y2**3*y3*y8*y9+y2*y3**3*y8*y9+y3*y4**3*y8*y9
    dF(1,378) = (-y6**2-y5**2)*y3*y1**3+(-y6**2-y5**2)*y3**3*y1+(y6**2+ &
      y5**2)*y4*y2**3+(y6**2+y5**2)*y4**3*y2
    dF(1,379) = ((y5*y9-sqrt(3._ark)*y6*y9)*y2+(-sqrt(3._ark)*y6*y8+y5*y8)*y4)*y1**3+ &
      ((sqrt(3._ark)*y6*y9-y5*y9)*y2**3+(sqrt(3._ark)*y6*y8-y5*y8)*y4**3)*y1+(- &
      sqrt(3._ark)*y6*y8+y5*y8)*y3*y2**3+(sqrt(3._ark)*y6*y8-y5*y8)*y3**3*y2+ &
      (sqrt(3._ark)*y6*y9-y5*y9)*y4*y3**3+(y5*y9-sqrt(3._ark)*y6*y9)*y4**3*y3
    dF(1,380) = (-y5*y7/2._ark+sqrt(3._ark)*y6*y7/2._ark)*y3*y1**3+(-y5*y7/2._ark+ &
      sqrt(3._ark)*y6*y7/2._ark)*y3**3*y1+(-y5*y7/2._ark+sqrt(3._ark)*y6*y7/2._ark)*y4*y2**3+(- &
      y5*y7/2._ark+sqrt(3._ark)*y6*y7/2._ark)*y4**3*y2
    dF(1,381) = ((y5*y6-sqrt(3._ark)*y6**2/3._ark)*y2+(y5*y6-sqrt(3._ark)*y6**2/ &
      3._ark)*y4)*y1**3+((sqrt(3._ark)*y6**2/3._ark-y5*y6)*y2**3+(sqrt(3._ark)*y6**2/3._ark- &
      y5*y6)*y4**3)*y1+(sqrt(3._ark)*y6**2/3._ark-y5*y6)*y3*y2**3+(y5*y6-sqrt(3._ark)*y6**2/ &
      3._ark)*y3**3*y2+(y5*y6-sqrt(3._ark)*y6**2/3._ark)*y4*y3**3+(sqrt(3._ark)*y6**2/3._ark- &
      y5*y6)*y4**3*y3
    dF(1,382) = ((-y6**2-y5**2)*y2+(-y6**2-y5**2)*y4)*y1**3+((y6**2+y5**2)*y2**3+ &
      (y6**2+y5**2)*y4**3)*y1+(y6**2+y5**2)*y3*y2**3+(-y6**2-y5**2)*y3**3*y2+(-y6**2- &
      y5**2)*y4*y3**3+(y6**2+y5**2)*y4**3*y3
    dF(1,383) = y2**3*y4**2*y7+y1**3*y3**2*y7+y2**2*y4**3*y7+y1**2*y3**3*y7
    dF(1,384) = (-y4**2*y9-y2**2*y8)*y1**3+(-y4**3*y9-y2**3*y8)*y1**2+ &
      y3**2*y4**3*y8+y2**3*y3**2*y9+y2**2*y3**3*y9+y3**3*y4**2*y8
    dF(1,385) = (y8+y9)*y3**2*y1**3+(-y8-y9)*y3**3*y1**2+(y8-y9)*y4**2*y2**3+(-y8+ &
      y9)*y4**3*y2**2
    dF(1,386) = ((-sqrt(3._ark)*y5+y6)*y2**2+2._ark*y4**2*y6)*y1**3+((-y6+ &
      sqrt(3._ark)*y5)*y2**3-2._ark*y4**3*y6)*y1**2-2._ark*y2**3*y3**2*y6+ &
      2._ark*y2**2*y3**3*y6+(-sqrt(3._ark)*y5+y6)*y4**2*y3**3+(-y6+ &
      sqrt(3._ark)*y5)*y4**3*y3**2
    dF(1,387) = -y1**3*y3**3+y2**3*y4**3
    dF(1,388) = (y8**2+y9**2)*y1**4+(-y9**2-y8**2)*y2**4+(y8**2+y9**2)*y3**4+(- &
      y9**2-y8**2)*y4**4
    dF(1,389) = (y6-sqrt(3._ark)*y5/3._ark)*y1**5+(sqrt(3._ark)*y5/3._ark-y6)*y2**5+(y6- &
      sqrt(3._ark)*y5/3._ark)*y3**5+(sqrt(3._ark)*y5/3._ark-y6)*y4**5
    dF(1,390) = (-y8**2*y9**3-y8**3*y9**2)*y1+(y8**2*y9**3-y8**3*y9**2)*y2+ &
      (y8**2*y9**3+y8**3*y9**2)*y3+(-y8**2*y9**3+y8**3*y9**2)*y4
    dF(1,391) = (y8**2+y9**2)*y7**3*y1+(y8**2+y9**2)*y7**3*y2+(y8**2+ &
      y9**2)*y7**3*y3+(y8**2+y9**2)*y7**3*y4
    dF(1,392) = (y9**4+y8**4)*y7*y1+(y9**4+y8**4)*y7*y2+(y9**4+y8**4)*y7*y3+(y9**4+ &
      y8**4)*y7*y4
    dF(1,393) = (y8**4*y9+y8*y9**4)*y1+(y8*y9**4-y8**4*y9)*y2+(-y8*y9**4- &
      y8**4*y9)*y3+(-y8*y9**4+y8**4*y9)*y4
    dF(1,394) = y3*y7**3*y8*y9+y1*y7**3*y8*y9-y4*y7**3*y8*y9-y2*y7**3*y8*y9
    dF(1,395) = ((-sqrt(3._ark)*y9**4/2._ark+sqrt(3._ark)*y8**4/2._ark)*y5+(y8**4/2._ark- &
      y9**4/2._ark)*y6)*y1+((sqrt(3._ark)*y9**4/2._ark-sqrt(3._ark)*y8**4/2._ark)*y5+(y9**4/ &
      2._ark-y8**4/2._ark)*y6)*y2+((-sqrt(3._ark)*y9**4/2._ark+sqrt(3._ark)*y8**4/2._ark)*y5+ &
      (y8**4/2._ark-y9**4/2._ark)*y6)*y3+((sqrt(3._ark)*y9**4/2._ark-sqrt(3._ark)*y8**4/ &
      2._ark)*y5+(y9**4/2._ark-y8**4/2._ark)*y6)*y4
    dF(1,396) = ((-sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/2._ark)*y7**3*y5+(y8/2._ark-y9/ &
      2._ark)*y7**3*y6)*y1+((-sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y9/2._ark)*y7**3*y5+(-y8/2._ark- &
      y9/2._ark)*y7**3*y6)*y2+((-sqrt(3._ark)*y8/2._ark+sqrt(3._ark)*y9/2._ark)*y7**3*y5+(y9/ &
      2._ark-y8/2._ark)*y7**3*y6)*y3+((sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/2._ark)*y7**3*y5+ &
      (y8/2._ark+y9/2._ark)*y7**3*y6)*y4
    dF(1,397) = (sqrt(3._ark)*y5*y7*y8*y9**2+(-y8*y9**2-2._ark*y8**2*y9)*y7*y6)*y1+(- &
      sqrt(3._ark)*y5*y7*y8*y9**2+(-2._ark*y8**2*y9+y8*y9**2)*y7*y6)*y2+(- &
      sqrt(3._ark)*y5*y7*y8*y9**2+(y8*y9**2+2._ark*y8**2*y9)*y7*y6)*y3+ &
      (sqrt(3._ark)*y5*y7*y8*y9**2+(-y8*y9**2+2._ark*y8**2*y9)*y7*y6)*y4
    dF(1,398) = (y5**4*y7+2._ark*y5**2*y6**2*y7+y6**4*y7)*y1+(y5**4*y7+ &
      2._ark*y5**2*y6**2*y7+y6**4*y7)*y2+(y5**4*y7+2._ark*y5**2*y6**2*y7+y6**4*y7)*y3+ &
      (y5**4*y7+2._ark*y5**2*y6**2*y7+y6**4*y7)*y4
    dF(1,399) = (-27._ark/5._ark*y5**4*y8-24._ark/5._ark*sqrt(3._ark)*y5**3*y6*y9+(-36._ark/ &
      5._ark*y9-54._ark/5._ark*y8)*y6**2*y5**2+16._ark/5._ark*sqrt(3._ark)*y5*y6**3*y8+(-28._ark/ &
      5._ark*y9+y8)*y6**4)*y1+(-27._ark/5._ark*y5**4*y8+24._ark/5._ark*sqrt(3._ark)*y5**3*y6*y9+ &
      (36._ark/5._ark*y9-54._ark/5._ark*y8)*y6**2*y5**2+16._ark/5._ark*sqrt(3._ark)*y5*y6**3*y8+(y8+ &
      28._ark/5._ark*y9)*y6**4)*y2+(27._ark/5._ark*y5**4*y8+24._ark/5._ark*sqrt(3._ark)*y5**3*y6*y9+ &
      (36._ark/5._ark*y9+54._ark/5._ark*y8)*y6**2*y5**2-16._ark/5._ark*sqrt(3._ark)*y5*y6**3*y8+ &
      (28._ark/5._ark*y9-y8)*y6**4)*y3+(27._ark/5._ark*y5**4*y8-24._ark/ &
      5._ark*sqrt(3._ark)*y5**3*y6*y9+(54._ark/5._ark*y8-36._ark/5._ark*y9)*y6**2*y5**2-16._ark/ &
      5._ark*sqrt(3._ark)*y5*y6**3*y8+(-28._ark/5._ark*y9-y8)*y6**4)*y4
    dF(1,400) = (sqrt(3._ark)*y6*y7**4/2._ark-y5*y7**4/2._ark)*y1+(-sqrt(3._ark)*y6*y7**4/ &
      2._ark+y5*y7**4/2._ark)*y2+(sqrt(3._ark)*y6*y7**4/2._ark-y5*y7**4/2._ark)*y3+(- &
      sqrt(3._ark)*y6*y7**4/2._ark+y5*y7**4/2._ark)*y4
    dF(1,401) = ((-sqrt(3._ark)*y9**2/2._ark-sqrt(3._ark)*y8**2/2._ark)*y7*y5**2+(-y9**2- &
      y8**2)*y7*y6*y5+(-sqrt(3._ark)*y9**2/6._ark-sqrt(3._ark)*y8**2/6._ark)*y7*y6**2)*y1+((- &
      sqrt(3._ark)*y9**2/2._ark-sqrt(3._ark)*y8**2/2._ark)*y7*y5**2+(-y9**2-y8**2)*y7*y6*y5+(- &
      sqrt(3._ark)*y9**2/6._ark-sqrt(3._ark)*y8**2/6._ark)*y7*y6**2)*y2+((-sqrt(3._ark)*y9**2/ &
      2._ark-sqrt(3._ark)*y8**2/2._ark)*y7*y5**2+(-y9**2-y8**2)*y7*y6*y5+(-sqrt(3._ark)*y9**2/ &
      6._ark-sqrt(3._ark)*y8**2/6._ark)*y7*y6**2)*y3+((-sqrt(3._ark)*y9**2/2._ark- &
      sqrt(3._ark)*y8**2/2._ark)*y7*y5**2+(-y9**2-y8**2)*y7*y6*y5+(-sqrt(3._ark)*y9**2/6._ark- &
      sqrt(3._ark)*y8**2/6._ark)*y7*y6**2)*y4
    dF(1,402) = ((-sqrt(3._ark)*y8*y9**2/3._ark-sqrt(3._ark)*y8**2*y9/3._ark)*y5**2+(- &
      y8*y9**2-y8**2*y9)*y6*y5)*y1+((-sqrt(3._ark)*y8*y9**2/3._ark+sqrt(3._ark)*y8**2*y9/ &
      3._ark)*y5**2+(-y8*y9**2+y8**2*y9)*y6*y5)*y2+((sqrt(3._ark)*y8**2*y9/3._ark+ &
      sqrt(3._ark)*y8*y9**2/3._ark)*y5**2+(y8**2*y9+y8*y9**2)*y6*y5)*y3+ &
      ((sqrt(3._ark)*y8*y9**2/3._ark-sqrt(3._ark)*y8**2*y9/3._ark)*y5**2+(y8*y9**2- &
      y8**2*y9)*y6*y5)*y4
    dF(1,403) = (sqrt(3._ark)*y5**2*y7**2*y8/2._ark+y5*y6*y7**2*y9+(-sqrt(3._ark)*y8/6._ark+ &
      sqrt(3._ark)*y9/3._ark)*y7**2*y6**2)*y1+(sqrt(3._ark)*y5**2*y7**2*y8/2._ark- &
      y5*y6*y7**2*y9+(-sqrt(3._ark)*y8/6._ark-sqrt(3._ark)*y9/3._ark)*y7**2*y6**2)*y2+(- &
      sqrt(3._ark)*y5**2*y7**2*y8/2._ark-y5*y6*y7**2*y9+(sqrt(3._ark)*y8/6._ark-sqrt(3._ark)*y9/ &
      3._ark)*y7**2*y6**2)*y3+(-sqrt(3._ark)*y5**2*y7**2*y8/2._ark+y5*y6*y7**2*y9+ &
      (sqrt(3._ark)*y9/3._ark+sqrt(3._ark)*y8/6._ark)*y7**2*y6**2)*y4
    dF(1,404) = (sqrt(3._ark)*y5**2*y7*y9**2/2._ark+y5*y6*y7*y8**2+(sqrt(3._ark)*y8**2/ &
      3._ark-sqrt(3._ark)*y9**2/6._ark)*y7*y6**2)*y1+(sqrt(3._ark)*y5**2*y7*y9**2/2._ark+ &
      y5*y6*y7*y8**2+(sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y9**2/6._ark)*y7*y6**2)*y2+ &
      (sqrt(3._ark)*y5**2*y7*y9**2/2._ark+y5*y6*y7*y8**2+(sqrt(3._ark)*y8**2/3._ark- &
      sqrt(3._ark)*y9**2/6._ark)*y7*y6**2)*y3+(sqrt(3._ark)*y5**2*y7*y9**2/2._ark+ &
      y5*y6*y7*y8**2+(sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y9**2/6._ark)*y7*y6**2)*y4
    dF(1,405) = (sqrt(3._ark)*y5**2*y9**3/2._ark+y5*y6*y8**3+(-sqrt(3._ark)*y9**3/6._ark+ &
      sqrt(3._ark)*y8**3/3._ark)*y6**2)*y1+(-sqrt(3._ark)*y5**2*y9**3/2._ark+y5*y6*y8**3+ &
      (sqrt(3._ark)*y9**3/6._ark+sqrt(3._ark)*y8**3/3._ark)*y6**2)*y2+(- &
      sqrt(3._ark)*y5**2*y9**3/2._ark-y5*y6*y8**3+(-sqrt(3._ark)*y8**3/3._ark+ &
      sqrt(3._ark)*y9**3/6._ark)*y6**2)*y3+(sqrt(3._ark)*y5**2*y9**3/2._ark-y5*y6*y8**3+(- &
      sqrt(3._ark)*y8**3/3._ark-sqrt(3._ark)*y9**3/6._ark)*y6**2)*y4
    dF(1,406) = (-y5*y6*y7**3+sqrt(3._ark)*y6**2*y7**3/2._ark+sqrt(3._ark)*y5**2*y7**3/ &
      6._ark)*y1+(-y5*y6*y7**3+sqrt(3._ark)*y6**2*y7**3/2._ark+sqrt(3._ark)*y5**2*y7**3/ &
      6._ark)*y2+(-y5*y6*y7**3+sqrt(3._ark)*y6**2*y7**3/2._ark+sqrt(3._ark)*y5**2*y7**3/ &
      6._ark)*y3+(-y5*y6*y7**3+sqrt(3._ark)*y6**2*y7**3/2._ark+sqrt(3._ark)*y5**2*y7**3/ &
      6._ark)*y4
    dF(1,407) = (y5*y6*y7*y8*y9-sqrt(3._ark)*y5**2*y7*y8*y9/6._ark- &
      sqrt(3._ark)*y6**2*y7*y8*y9/2._ark)*y1+(sqrt(3._ark)*y5**2*y7*y8*y9/6._ark+ &
      sqrt(3._ark)*y6**2*y7*y8*y9/2._ark-y5*y6*y7*y8*y9)*y2+(y5*y6*y7*y8*y9- &
      sqrt(3._ark)*y5**2*y7*y8*y9/6._ark-sqrt(3._ark)*y6**2*y7*y8*y9/2._ark)*y3+ &
      (sqrt(3._ark)*y5**2*y7*y8*y9/6._ark+sqrt(3._ark)*y6**2*y7*y8*y9/2._ark- &
      y5*y6*y7*y8*y9)*y4
    dF(1,408) = (3._ark/16._ark*y5**3*y8**2+5._ark/16._ark*sqrt(3._ark)*y5**2*y6*y9**2+(- &
      13._ark/16._ark*y8**2-y9**2)*y6**2*y5+(sqrt(3._ark)*y8**2/4._ark+sqrt(3._ark)*y9**2/ &
      16._ark)*y6**3)*y1+(-3._ark/16._ark*y5**3*y8**2-5._ark/16._ark*sqrt(3._ark)*y5**2*y6*y9**2+ &
      (y9**2+13._ark/16._ark*y8**2)*y6**2*y5+(-sqrt(3._ark)*y9**2/16._ark-sqrt(3._ark)*y8**2/ &
      4._ark)*y6**3)*y2+(3._ark/16._ark*y5**3*y8**2+5._ark/16._ark*sqrt(3._ark)*y5**2*y6*y9**2+(- &
      13._ark/16._ark*y8**2-y9**2)*y6**2*y5+(sqrt(3._ark)*y8**2/4._ark+sqrt(3._ark)*y9**2/ &
      16._ark)*y6**3)*y3+(-3._ark/16._ark*y5**3*y8**2-5._ark/16._ark*sqrt(3._ark)*y5**2*y6*y9**2+ &
      (y9**2+13._ark/16._ark*y8**2)*y6**2*y5+(-sqrt(3._ark)*y9**2/16._ark-sqrt(3._ark)*y8**2/ &
      4._ark)*y6**3)*y4
    dF(1,409) = ((-y9**3-y8**3)*y5**2+(-y9**3-y8**3)*y6**2)*y1+((y9**3-y8**3)*y5**2+ &
      (y9**3-y8**3)*y6**2)*y2+((y8**3+y9**3)*y5**2+(y8**3+y9**3)*y6**2)*y3+((-y9**3+ &
      y8**3)*y5**2+(-y9**3+y8**3)*y6**2)*y4
    dF(1,410) = (sqrt(3._ark)*y5**3*y8**2/4._ark+(-y8**2+y9**2/4._ark)*y6*y5**2+ &
      sqrt(3._ark)*y5*y6**2*y8**2/4._ark-3._ark/4._ark*y6**3*y9**2)*y1+(- &
      sqrt(3._ark)*y5**3*y8**2/4._ark+(-y9**2/4._ark+y8**2)*y6*y5**2- &
      sqrt(3._ark)*y5*y6**2*y8**2/4._ark+3._ark/4._ark*y6**3*y9**2)*y2+ &
      (sqrt(3._ark)*y5**3*y8**2/4._ark+(-y8**2+y9**2/4._ark)*y6*y5**2+ &
      sqrt(3._ark)*y5*y6**2*y8**2/4._ark-3._ark/4._ark*y6**3*y9**2)*y3+(- &
      sqrt(3._ark)*y5**3*y8**2/4._ark+(-y9**2/4._ark+y8**2)*y6*y5**2- &
      sqrt(3._ark)*y5*y6**2*y8**2/4._ark+3._ark/4._ark*y6**3*y9**2)*y4
    dF(1,411) = (-sqrt(3._ark)*y5*y6**2*y7**2-sqrt(3._ark)*y5**3*y7**2/9._ark+y6**3*y7**2+ &
      y5**2*y6*y7**2)*y1+(-y6**3*y7**2+sqrt(3._ark)*y5**3*y7**2/9._ark-y5**2*y6*y7**2+ &
      sqrt(3._ark)*y5*y6**2*y7**2)*y2+(-sqrt(3._ark)*y5*y6**2*y7**2- &
      sqrt(3._ark)*y5**3*y7**2/9._ark+y6**3*y7**2+y5**2*y6*y7**2)*y3+(-y6**3*y7**2+ &
      sqrt(3._ark)*y5**3*y7**2/9._ark-y5**2*y6*y7**2+sqrt(3._ark)*y5*y6**2*y7**2)*y4
    dF(1,412) = (sqrt(3._ark)*y5**3*y7*y8+(-4._ark*y8+y9)*y7*y6*y5**2+ &
      sqrt(3._ark)*y5*y6**2*y7*y8-3._ark*y6**3*y7*y9)*y1+(-sqrt(3._ark)*y5**3*y7*y8+(y9+ &
      4._ark*y8)*y7*y6*y5**2-sqrt(3._ark)*y5*y6**2*y7*y8-3._ark*y6**3*y7*y9)*y2+(- &
      sqrt(3._ark)*y5**3*y7*y8+(4._ark*y8-y9)*y7*y6*y5**2-sqrt(3._ark)*y5*y6**2*y7*y8+ &
      3._ark*y6**3*y7*y9)*y3+(sqrt(3._ark)*y5**3*y7*y8+(-y9-4._ark*y8)*y7*y6*y5**2+ &
      sqrt(3._ark)*y5*y6**2*y7*y8+3._ark*y6**3*y7*y9)*y4
    dF(1,413) = ((-y9**2-7._ark/16._ark*y8**2)*y5**3+15._ark/ &
      16._ark*sqrt(3._ark)*y5**2*y6*y9**2+9._ark/16._ark*y5*y6**2*y8**2+(3._ark/ &
      4._ark*sqrt(3._ark)*y8**2+3._ark/16._ark*sqrt(3._ark)*y9**2)*y6**3)*y1+((y9**2+7._ark/ &
      16._ark*y8**2)*y5**3-15._ark/16._ark*sqrt(3._ark)*y5**2*y6*y9**2-9._ark/ &
      16._ark*y5*y6**2*y8**2+(-3._ark/16._ark*sqrt(3._ark)*y9**2-3._ark/ &
      4._ark*sqrt(3._ark)*y8**2)*y6**3)*y2+((-y9**2-7._ark/16._ark*y8**2)*y5**3+15._ark/ &
      16._ark*sqrt(3._ark)*y5**2*y6*y9**2+9._ark/16._ark*y5*y6**2*y8**2+(3._ark/ &
      4._ark*sqrt(3._ark)*y8**2+3._ark/16._ark*sqrt(3._ark)*y9**2)*y6**3)*y3+((y9**2+7._ark/ &
      16._ark*y8**2)*y5**3-15._ark/16._ark*sqrt(3._ark)*y5**2*y6*y9**2-9._ark/ &
      16._ark*y5*y6**2*y8**2+(-3._ark/16._ark*sqrt(3._ark)*y9**2-3._ark/ &
      4._ark*sqrt(3._ark)*y8**2)*y6**3)*y4
    dF(1,414) = (13._ark/5._ark*sqrt(3._ark)*y5**4*y8+(y8+28._ark/5._ark*y9)*y6*y5**3+(21._ark/ &
      5._ark*sqrt(3._ark)*y8+24._ark/5._ark*sqrt(3._ark)*y9)*y6**2*y5**2-27._ark/5._ark*y5*y6**3*y8+ &
      12._ark/5._ark*sqrt(3._ark)*y6**4*y9)*y1+(13._ark/5._ark*sqrt(3._ark)*y5**4*y8+(-28._ark/ &
      5._ark*y9+y8)*y6*y5**3+(-24._ark/5._ark*sqrt(3._ark)*y9+21._ark/ &
      5._ark*sqrt(3._ark)*y8)*y6**2*y5**2-27._ark/5._ark*y5*y6**3*y8-12._ark/ &
      5._ark*sqrt(3._ark)*y6**4*y9)*y2+(-13._ark/5._ark*sqrt(3._ark)*y5**4*y8+(-28._ark/5._ark*y9- &
      y8)*y6*y5**3+(-24._ark/5._ark*sqrt(3._ark)*y9-21._ark/5._ark*sqrt(3._ark)*y8)*y6**2*y5**2+ &
      27._ark/5._ark*y5*y6**3*y8-12._ark/5._ark*sqrt(3._ark)*y6**4*y9)*y3+(-13._ark/ &
      5._ark*sqrt(3._ark)*y5**4*y8+(28._ark/5._ark*y9-y8)*y6*y5**3+(24._ark/5._ark*sqrt(3._ark)*y9- &
      21._ark/5._ark*sqrt(3._ark)*y8)*y6**2*y5**2+27._ark/5._ark*y5*y6**3*y8+12._ark/ &
      5._ark*sqrt(3._ark)*y6**4*y9)*y4
    dF(1,415) = (-sqrt(3._ark)*y5**2*y6**2*y7/2._ark-sqrt(3._ark)*y5**4*y7/6._ark- &
      y5**3*y6*y7)*y1+(-sqrt(3._ark)*y5**2*y6**2*y7/2._ark-sqrt(3._ark)*y5**4*y7/6._ark- &
      y5**3*y6*y7)*y2+(-sqrt(3._ark)*y5**2*y6**2*y7/2._ark-sqrt(3._ark)*y5**4*y7/6._ark- &
      y5**3*y6*y7)*y3+(-sqrt(3._ark)*y5**2*y6**2*y7/2._ark-sqrt(3._ark)*y5**4*y7/6._ark- &
      y5**3*y6*y7)*y4
    dF(1,416) = ((-y9+28._ark/5._ark*y8)*y5**4+16._ark/5._ark*sqrt(3._ark)*y5**3*y6*y9+(36._ark/ &
      5._ark*y8+54._ark/5._ark*y9)*y6**2*y5**2-24._ark/5._ark*sqrt(3._ark)*y5*y6**3*y8+27._ark/ &
      5._ark*y6**4*y9)*y1+((y9+28._ark/5._ark*y8)*y5**4-16._ark/5._ark*sqrt(3._ark)*y5**3*y6*y9+(- &
      54._ark/5._ark*y9+36._ark/5._ark*y8)*y6**2*y5**2-24._ark/5._ark*sqrt(3._ark)*y5*y6**3*y8- &
      27._ark/5._ark*y6**4*y9)*y2+((y9-28._ark/5._ark*y8)*y5**4-16._ark/ &
      5._ark*sqrt(3._ark)*y5**3*y6*y9+(-54._ark/5._ark*y9-36._ark/5._ark*y8)*y6**2*y5**2+24._ark/ &
      5._ark*sqrt(3._ark)*y5*y6**3*y8-27._ark/5._ark*y6**4*y9)*y3+((-28._ark/5._ark*y8-y9)*y5**4+ &
      16._ark/5._ark*sqrt(3._ark)*y5**3*y6*y9+(54._ark/5._ark*y9-36._ark/5._ark*y8)*y6**2*y5**2+ &
      24._ark/5._ark*sqrt(3._ark)*y5*y6**3*y8+27._ark/5._ark*y6**4*y9)*y4
    dF(1,417) = (-3._ark*sqrt(3._ark)*y6**5+7._ark*y5**5-15._ark*y5*y6**4- &
      5._ark*sqrt(3._ark)*y5**4*y6)*y1+(5._ark*sqrt(3._ark)*y5**4*y6+3._ark*sqrt(3._ark)*y6**5- &
      7._ark*y5**5+15._ark*y5*y6**4)*y2+(-3._ark*sqrt(3._ark)*y6**5+7._ark*y5**5-15._ark*y5*y6**4- &
      5._ark*sqrt(3._ark)*y5**4*y6)*y3+(5._ark*sqrt(3._ark)*y5**4*y6+3._ark*sqrt(3._ark)*y6**5- &
      7._ark*y5**5+15._ark*y5*y6**4)*y4
    dF(1,418) = (2._ark*y2*y6*y7*y8**2+(-sqrt(3._ark)*y5*y7*y9**2+y6*y7*y9**2)*y4)*y1+(- &
      sqrt(3._ark)*y5*y7*y9**2+y6*y7*y9**2)*y3*y2+2._ark*y3*y4*y6*y7*y8**2
    dF(1,419) = ((3._ark/5._ark*sqrt(3._ark)*y5**2*y6*y7-4._ark/5._ark*y5**3*y7+3._ark/ &
      5._ark*sqrt(3._ark)*y6**3*y7)*y2+(3._ark/5._ark*sqrt(3._ark)*y5**2*y6*y7-4._ark/ &
      5._ark*y5**3*y7+3._ark/5._ark*sqrt(3._ark)*y6**3*y7)*y4)*y1+(3._ark/ &
      5._ark*sqrt(3._ark)*y5**2*y6*y7-4._ark/5._ark*y5**3*y7+3._ark/ &
      5._ark*sqrt(3._ark)*y6**3*y7)*y3*y2+(3._ark/5._ark*sqrt(3._ark)*y5**2*y6*y7-4._ark/ &
      5._ark*y5**3*y7+3._ark/5._ark*sqrt(3._ark)*y6**3*y7)*y4*y3
    dF(1,420) = ((9._ark/4._ark*y5*y6**2*y8+5._ark/4._ark*y5**3*y8)*y2+(-3._ark/ &
      4._ark*sqrt(3._ark)*y5**2*y6*y9-y5**3*y9-3._ark/4._ark*sqrt(3._ark)*y6**3*y9)*y4)*y1+ &
      (3._ark/4._ark*sqrt(3._ark)*y6**3*y9+3._ark/4._ark*sqrt(3._ark)*y5**2*y6*y9+y5**3*y9)*y3*y2+ &
      (-5._ark/4._ark*y5**3*y8-9._ark/4._ark*y5*y6**2*y8)*y4*y3
    dF(1,421) = (y8**3+y9**3)*y3*y1**2+(-y9**3-y8**3)*y3**2*y1+(-y9**3+ &
      y8**3)*y4*y2**2+(y9**3-y8**3)*y4**2*y2
    dF(1,422) = (-y2*y8**3-y4*y9**3)*y1**2+(-y4**2*y9**3-y2**2*y8**3)*y1+ &
      y3*y4**2*y8**3+y3**2*y4*y8**3+y2**2*y3*y9**3+y2*y3**2*y9**3
    dF(1,423) = (y4*y7**3+y2*y7**3)*y1**2+(y4**2*y7**3+y2**2*y7**3)*y1+ &
      y3**2*y4*y7**3+y3*y4**2*y7**3+y2*y3**2*y7**3+y2**2*y3*y7**3
    dF(1,424) = (-y8*y9**2-y8**2*y9)*y3*y1**2+(y8**2*y9+y8*y9**2)*y3**2*y1+(- &
      y8*y9**2+y8**2*y9)*y4*y2**2+(y8*y9**2-y8**2*y9)*y4**2*y2
    dF(1,425) = ((-y5*y7**2/2._ark-sqrt(3._ark)*y6*y7**2/2._ark)*y2+y4*y5*y7**2)*y1**2+ &
      ((sqrt(3._ark)*y6*y7**2/2._ark+y5*y7**2/2._ark)*y2**2-y4**2*y5*y7**2)*y1- &
      y2**2*y3*y5*y7**2+y2*y3**2*y5*y7**2+(-y5*y7**2/2._ark-sqrt(3._ark)*y6*y7**2/ &
      2._ark)*y4*y3**2+(sqrt(3._ark)*y6*y7**2/2._ark+y5*y7**2/2._ark)*y4**2*y3
    dF(1,426) = ((-sqrt(3._ark)*y6*y8**2-y5*y8**2)*y2+2._ark*y4*y5*y9**2)*y1**2+ &
      ((sqrt(3._ark)*y6*y8**2+y5*y8**2)*y2**2-2._ark*y4**2*y5*y9**2)*y1- &
      2._ark*y2**2*y3*y5*y9**2+2._ark*y2*y3**2*y5*y9**2+(-sqrt(3._ark)*y6*y8**2- &
      y5*y8**2)*y4*y3**2+(sqrt(3._ark)*y6*y8**2+y5*y8**2)*y4**2*y3
    dF(1,427) = ((-sqrt(3._ark)*y5+y6)*y2+2._ark*y4*y6)*y1**4+((-y6+ &
      sqrt(3._ark)*y5)*y2**4-2._ark*y4**4*y6)*y1-2._ark*y2**4*y3*y6+2._ark*y2*y3**4*y6+(- &
      sqrt(3._ark)*y5+y6)*y4*y3**4+(-y6+sqrt(3._ark)*y5)*y4**4*y3
    dF(1,428) = ((y5-sqrt(3._ark)*y6)*y2+(y5-sqrt(3._ark)*y6)*y4)*y1**4+((sqrt(3._ark)*y6- &
      y5)*y2**4+(sqrt(3._ark)*y6-y5)*y4**4)*y1+(sqrt(3._ark)*y6-y5)*y3*y2**4+(y5- &
      sqrt(3._ark)*y6)*y3**4*y2+(y5-sqrt(3._ark)*y6)*y4*y3**4+(sqrt(3._ark)*y6-y5)*y4**4*y3
    dF(1,429) = (y4*y8**3*y9+y2*y8*y9**3)*y1+y3*y4*y8*y9**3+y2*y3*y8**3*y9
    dF(1,430) = (-y4*y7*y8**3-y2*y7*y9**3)*y1+y3*y4*y7*y9**3+y2*y3*y7*y8**3
    dF(1,431) = y2*y4*y7**2*y8*y9+y1*y3*y7**2*y8*y9
    dF(1,432) = ((y8**2+y9**2)*y7*y5+(-sqrt(3._ark)*y8**2- &
      sqrt(3._ark)*y9**2)*y7*y6)*y3*y1+((y8**2+y9**2)*y7*y5+(-sqrt(3._ark)*y8**2- &
      sqrt(3._ark)*y9**2)*y7*y6)*y4*y2
    dF(1,433) = (-2._ark*y2*y5*y7*y9**2+(y5*y7*y8**2+sqrt(3._ark)*y6*y7*y8**2)*y4)*y1+ &
      (y5*y7*y8**2+sqrt(3._ark)*y6*y7*y8**2)*y3*y2-2._ark*y3*y4*y5*y7*y9**2
    dF(1,434) = (sqrt(3._ark)*y5**2*y7**2/2._ark+sqrt(3._ark)*y6**2*y7**2/6._ark+ &
      y5*y6*y7**2)*y3*y1+(-sqrt(3._ark)*y6**2*y7**2/6._ark-y5*y6*y7**2- &
      sqrt(3._ark)*y5**2*y7**2/2._ark)*y4*y2
    dF(1,435) = (((-y8+y9)*y7**2*y3+(-y8-y9)*y7**2*y4)*y2+(y8-y9)*y7**2*y4*y3)*y1+ &
      (y8+y9)*y7**2*y4*y3*y2
    dF(1,436) = (((-y8*y9**2+y8**2*y9)*y3+(-y8*y9**2-y8**2*y9)*y4)*y2+(y8*y9**2- &
      y8**2*y9)*y4*y3)*y1+(y8**2*y9+y8*y9**2)*y4*y3*y2
    dF(1,437) = (((y8**2+y9**2)*y7*y3+(y8**2+y9**2)*y7*y4)*y2+(y8**2+ &
      y9**2)*y7*y4*y3)*y1+(y8**2+y9**2)*y7*y4*y3*y2
    dF(1,438) = ((((-sqrt(3._ark)*y8**2+sqrt(3._ark)*y9**2)*y5+(y9**2-y8**2)*y6)*y3+ &
      ((sqrt(3._ark)*y8**2-sqrt(3._ark)*y9**2)*y5+(y8**2-y9**2)*y6)*y4)*y2+((- &
      sqrt(3._ark)*y8**2+sqrt(3._ark)*y9**2)*y5+(y9**2-y8**2)*y6)*y4*y3)*y1+ &
      ((sqrt(3._ark)*y8**2-sqrt(3._ark)*y9**2)*y5+(y8**2-y9**2)*y6)*y4*y3*y2
    dF(1,439) = y1**3*y2*y4*y7+(y3*y4**3*y7+y2**3*y3*y7)*y1+y2*y3**3*y4*y7
    dF(1,440) = (y3*y4*y8+y2*y3*y9)*y1**3+(-y2**3*y4*y9+(-y3**3*y8-y4**3*y8)*y2- &
      y3**3*y4*y9)*y1+y2*y3*y4**3*y9+y2**3*y3*y4*y8
    dF(1,441) = ((-y6-sqrt(3._ark)*y5)*y3*y2+(y6+sqrt(3._ark)*y5)*y4*y3)*y1**3+((y6+ &
      sqrt(3._ark)*y5)*y4*y2**3+((y6+sqrt(3._ark)*y5)*y3**3+(-y6-sqrt(3._ark)*y5)*y4**3)*y2+ &
      (-y6-sqrt(3._ark)*y5)*y4*y3**3)*y1+(-y6-sqrt(3._ark)*y5)*y4*y3*y2**3+(y6+ &
      sqrt(3._ark)*y5)*y4**3*y3*y2
    dF(1,442) = (2._ark*y2*y3*y5+(-y5-sqrt(3._ark)*y6)*y4*y3)*y1**3+(-2._ark*y2**3*y4*y5+ &
      ((-y5-sqrt(3._ark)*y6)*y3**3+(sqrt(3._ark)*y6+y5)*y4**3)*y2+2._ark*y3**3*y4*y5)*y1+ &
      (sqrt(3._ark)*y6+y5)*y4*y3*y2**3-2._ark*y2*y3*y4**3*y5
    dF(1,443) = ((-y6*y8**2-sqrt(3._ark)*y5*y8**2)*y2+(sqrt(3._ark)*y5*y9**2+ &
      y6*y9**2)*y4)*y1**2+((sqrt(3._ark)*y5*y8**2+y6*y8**2)*y2**2+(-sqrt(3._ark)*y5*y9**2- &
      y6*y9**2)*y4**2)*y1+(-sqrt(3._ark)*y5*y9**2-y6*y9**2)*y3*y2**2+ &
      (sqrt(3._ark)*y5*y9**2+y6*y9**2)*y3**2*y2+(-y6*y8**2- &
      sqrt(3._ark)*y5*y8**2)*y4*y3**2+(sqrt(3._ark)*y5*y8**2+y6*y8**2)*y4**2*y3
    dF(1,444) = (-sqrt(3._ark)*y5*y9**2+(2._ark*y8**2+y9**2)*y6)*y3*y1**2+(- &
      sqrt(3._ark)*y5*y9**2+(2._ark*y8**2+y9**2)*y6)*y3**2*y1+(sqrt(3._ark)*y5*y9**2+(- &
      2._ark*y8**2-y9**2)*y6)*y4*y2**2+(sqrt(3._ark)*y5*y9**2+(-2._ark*y8**2- &
      y9**2)*y6)*y4**2*y2
    dF(1,445) = ((-sqrt(3._ark)*y5*y7**2/2._ark+y6*y7**2/2._ark)*y2+y4*y6*y7**2)*y1**2+ &
      ((sqrt(3._ark)*y5*y7**2/2._ark-y6*y7**2/2._ark)*y2**2-y4**2*y6*y7**2)*y1- &
      y2**2*y3*y6*y7**2+y2*y3**2*y6*y7**2+(-sqrt(3._ark)*y5*y7**2/2._ark+y6*y7**2/ &
      2._ark)*y4*y3**2+(sqrt(3._ark)*y5*y7**2/2._ark-y6*y7**2/2._ark)*y4**2*y3
    dF(1,446) = ((-y6**2*y8-y5**2*y8)*y2+(-y6**2*y9-y5**2*y9)*y4)*y1**2+((-y6**2*y8- &
      y5**2*y8)*y2**2+(-y6**2*y9-y5**2*y9)*y4**2)*y1+(y6**2*y9+y5**2*y9)*y3*y2**2+ &
      (y6**2*y9+y5**2*y9)*y3**2*y2+(y5**2*y8+y6**2*y8)*y4*y3**2+(y5**2*y8+ &
      y6**2*y8)*y4**2*y3
    dF(1,447) = (y6**3-4._ark/9._ark*sqrt(3._ark)*y5**3+y5**2*y6)*y3*y1**2+(y6**3-4._ark/ &
      9._ark*sqrt(3._ark)*y5**3+y5**2*y6)*y3**2*y1+(-y5**2*y6+4._ark/9._ark*sqrt(3._ark)*y5**3- &
      y6**3)*y4*y2**2+(-y5**2*y6+4._ark/9._ark*sqrt(3._ark)*y5**3-y6**3)*y4**2*y2
    dF(1,448) = ((y5**3+y5*y6**2-4._ark/3._ark*sqrt(3._ark)*y5**2*y6)*y2+ &
      (sqrt(3._ark)*y5**2*y6/3._ark-sqrt(3._ark)*y6**3)*y4)*y1**2+((4._ark/ &
      3._ark*sqrt(3._ark)*y5**2*y6-y5*y6**2-y5**3)*y2**2+(-sqrt(3._ark)*y5**2*y6/3._ark+ &
      sqrt(3._ark)*y6**3)*y4**2)*y1+(-sqrt(3._ark)*y5**2*y6/3._ark+ &
      sqrt(3._ark)*y6**3)*y3*y2**2+(sqrt(3._ark)*y5**2*y6/3._ark-sqrt(3._ark)*y6**3)*y3**2*y2+ &
      (y5**3+y5*y6**2-4._ark/3._ark*sqrt(3._ark)*y5**2*y6)*y4*y3**2+(4._ark/ &
      3._ark*sqrt(3._ark)*y5**2*y6-y5*y6**2-y5**3)*y4**2*y3
    dF(1,449) = (-y2*y3*y7**2-y3*y4*y7**2)*y1**2+(y2**2*y4*y7**2+(y4**2*y7**2- &
      y3**2*y7**2)*y2-y3**2*y4*y7**2)*y1+y2*y3*y4**2*y7**2+y2**2*y3*y4*y7**2
    dF(1,450) = (y3*y4*y7*y9+y2*y3*y7*y8)*y1**2+(-y2**2*y4*y7*y8+(-y4**2*y7*y9- &
      y3**2*y7*y9)*y2-y3**2*y4*y7*y8)*y1+y2*y3*y4**2*y7*y8+y2**2*y3*y4*y7*y9
    dF(1,451) = ((-y3*y8-y4*y8)*y2**2-y2*y4**2*y9-y3*y4**2*y9)*y1**2+ &
      (y3**2*y4**2*y8+y2**2*y3**2*y9)*y1+y2*y3**2*y4**2*y8+y2**2*y3**2*y4*y9
    dF(1,452) = ((-y8+y9)*y3**2*y2+(y8-y9)*y4*y3**2)*y1**2+(-y8-y9)*y4**2*y2**2*y1+ &
      (y8+y9)*y4**2*y3*y2**2
    dF(1,453) = ((-sqrt(3._ark)*y6/2._ark+y5/2._ark)*y3**2*y2+(-sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y4*y3**2)*y1**2+(sqrt(3._ark)*y6/2._ark-y5/2._ark)*y4**2*y2**2*y1+ &
      (sqrt(3._ark)*y6/2._ark-y5/2._ark)*y4**2*y3*y2**2
    dF(1,454) = (y2*y7**2+y4*y7**2)*y1**3+(-y4**3*y7**2-y2**3*y7**2)*y1- &
      y3*y4**3*y7**2-y2**3*y3*y7**2+y2*y3**3*y7**2+y3**3*y4*y7**2
    dF(1,455) = (y4+y2)*y1**5+(-y2**5-y4**5)*y1-y2**5*y3+y3**5*y4-y3*y4**5+y2*y3**5
    dF(1,456) = (-sqrt(3._ark)*y5*y7**2*y8/2._ark+(y8/2._ark+y9)*y7**2*y6)*y1**2+(- &
      sqrt(3._ark)*y5*y7**2*y8/2._ark+(y8/2._ark-y9)*y7**2*y6)*y2**2+ &
      (sqrt(3._ark)*y5*y7**2*y8/2._ark+(-y8/2._ark-y9)*y7**2*y6)*y3**2+ &
      (sqrt(3._ark)*y5*y7**2*y8/2._ark+(-y8/2._ark+y9)*y7**2*y6)*y4**2
    dF(1,457) = (sqrt(3._ark)*y5*y7*y8**2+(-2._ark*y9**2-y8**2)*y7*y6)*y1**2+ &
      (sqrt(3._ark)*y5*y7*y8**2+(-2._ark*y9**2-y8**2)*y7*y6)*y2**2+ &
      (sqrt(3._ark)*y5*y7*y8**2+(-2._ark*y9**2-y8**2)*y7*y6)*y3**2+ &
      (sqrt(3._ark)*y5*y7*y8**2+(-2._ark*y9**2-y8**2)*y7*y6)*y4**2
    dF(1,458) = ((y8**2+y9**2)*y5**2+(y8**2+y9**2)*y6**2)*y1**2+((-y9**2- &
      y8**2)*y5**2+(-y9**2-y8**2)*y6**2)*y2**2+((y8**2+y9**2)*y5**2+(y8**2+ &
      y9**2)*y6**2)*y3**2+((-y9**2-y8**2)*y5**2+(-y9**2-y8**2)*y6**2)*y4**2
    dF(1,459) = (y6**2*y8*y9+y5**2*y8*y9)*y1**2+(y6**2*y8*y9+y5**2*y8*y9)*y2**2+ &
      (y6**2*y8*y9+y5**2*y8*y9)*y3**2+(y6**2*y8*y9+y5**2*y8*y9)*y4**2
    dF(1,460) = ((y8**2+y9**2)*y7*y5+(-sqrt(3._ark)*y8**2- &
      sqrt(3._ark)*y9**2)*y7*y6)*y1**2+((y8**2+y9**2)*y7*y5+(-sqrt(3._ark)*y8**2- &
      sqrt(3._ark)*y9**2)*y7*y6)*y2**2+((y8**2+y9**2)*y7*y5+(-sqrt(3._ark)*y8**2- &
      sqrt(3._ark)*y9**2)*y7*y6)*y3**2+((y8**2+y9**2)*y7*y5+(-sqrt(3._ark)*y8**2- &
      sqrt(3._ark)*y9**2)*y7*y6)*y4**2
    dF(1,461) = ((-y8*y9**2+2._ark*y8**2*y9)*y5-sqrt(3._ark)*y6*y8*y9**2)*y1**2+((- &
      y8*y9**2-2._ark*y8**2*y9)*y5-sqrt(3._ark)*y6*y8*y9**2)*y2**2+((-2._ark*y8**2*y9+ &
      y8*y9**2)*y5+sqrt(3._ark)*y6*y8*y9**2)*y3**2+((y8*y9**2+2._ark*y8**2*y9)*y5+ &
      sqrt(3._ark)*y6*y8*y9**2)*y4**2
    dF(1,462) = ((y8/2._ark-y9)*y7**2*y5+sqrt(3._ark)*y6*y7**2*y8/2._ark)*y1**2+((y8/2._ark+ &
      y9)*y7**2*y5+sqrt(3._ark)*y6*y7**2*y8/2._ark)*y2**2+((-y8/2._ark+y9)*y7**2*y5- &
      sqrt(3._ark)*y6*y7**2*y8/2._ark)*y3**2+((-y8/2._ark-y9)*y7**2*y5- &
      sqrt(3._ark)*y6*y7**2*y8/2._ark)*y4**2
    dF(1,463) = (sqrt(3._ark)*y5**2*y7*y9/2._ark+y5*y6*y7*y8+(-sqrt(3._ark)*y9/6._ark+ &
      sqrt(3._ark)*y8/3._ark)*y7*y6**2)*y1**2+(sqrt(3._ark)*y5**2*y7*y9/2._ark-y5*y6*y7*y8+(- &
      sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y9/6._ark)*y7*y6**2)*y2**2+(-sqrt(3._ark)*y5**2*y7*y9/ &
      2._ark-y5*y6*y7*y8+(sqrt(3._ark)*y9/6._ark-sqrt(3._ark)*y8/3._ark)*y7*y6**2)*y3**2+(- &
      sqrt(3._ark)*y5**2*y7*y9/2._ark+y5*y6*y7*y8+(sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y8/ &
      3._ark)*y7*y6**2)*y4**2
    dF(1,464) = ((-sqrt(3._ark)*y9**2/6._ark-sqrt(3._ark)*y8**2/6._ark)*y5**2+(y8**2+ &
      y9**2)*y6*y5+(-sqrt(3._ark)*y9**2/2._ark-sqrt(3._ark)*y8**2/2._ark)*y6**2)*y1**2+ &
      ((sqrt(3._ark)*y9**2/6._ark+sqrt(3._ark)*y8**2/6._ark)*y5**2+(-y9**2-y8**2)*y6*y5+ &
      (sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y9**2/2._ark)*y6**2)*y2**2+((-sqrt(3._ark)*y9**2/ &
      6._ark-sqrt(3._ark)*y8**2/6._ark)*y5**2+(y8**2+y9**2)*y6*y5+(-sqrt(3._ark)*y9**2/2._ark- &
      sqrt(3._ark)*y8**2/2._ark)*y6**2)*y3**2+((sqrt(3._ark)*y9**2/6._ark+sqrt(3._ark)*y8**2/ &
      6._ark)*y5**2+(-y9**2-y8**2)*y6*y5+(sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y9**2/ &
      2._ark)*y6**2)*y4**2
    dF(1,465) = ((sqrt(3._ark)*y8**2/6._ark-sqrt(3._ark)*y9**2/3._ark)*y5**2+y5*y6*y9**2- &
      sqrt(3._ark)*y6**2*y8**2/2._ark)*y1**2+((-sqrt(3._ark)*y8**2/6._ark+sqrt(3._ark)*y9**2/ &
      3._ark)*y5**2-y5*y6*y9**2+sqrt(3._ark)*y6**2*y8**2/2._ark)*y2**2+((sqrt(3._ark)*y8**2/ &
      6._ark-sqrt(3._ark)*y9**2/3._ark)*y5**2+y5*y6*y9**2-sqrt(3._ark)*y6**2*y8**2/ &
      2._ark)*y3**2+((-sqrt(3._ark)*y8**2/6._ark+sqrt(3._ark)*y9**2/3._ark)*y5**2-y5*y6*y9**2+ &
      sqrt(3._ark)*y6**2*y8**2/2._ark)*y4**2
    dF(1,466) = (y8+y9)*y7**2*y3*y1**2+(-y8-y9)*y7**2*y3**2*y1+(y8- &
      y9)*y7**2*y4*y2**2+(-y8+y9)*y7**2*y4**2*y2
    dF(1,467) = (-y2*y7**2*y8-y4*y7**2*y9)*y1**2+(-y2**2*y7**2*y8- &
      y4**2*y7**2*y9)*y1+y3**2*y4*y7**2*y8+y2*y3**2*y7**2*y9+y3*y4**2*y7**2*y8+ &
      y2**2*y3*y7**2*y9
    dF(1,468) = (y4*y7*y9**2+y2*y7*y8**2)*y1**2+(y4**2*y7*y9**2+y2**2*y7*y8**2)*y1+ &
      y3**2*y4*y7*y8**2+y2*y3**2*y7*y9**2+y2**2*y3*y7*y9**2+y3*y4**2*y7*y8**2
    dF(1,469) = (y6*y7**2-sqrt(3._ark)*y5*y7**2/3._ark)*y3*y1**2+(y6*y7**2- &
      sqrt(3._ark)*y5*y7**2/3._ark)*y3**2*y1+(-y6*y7**2+sqrt(3._ark)*y5*y7**2/ &
      3._ark)*y4*y2**2+(-y6*y7**2+sqrt(3._ark)*y5*y7**2/3._ark)*y4**2*y2
    dF(1,470) = ((y6*y9**2-sqrt(3._ark)*y5*y9**2)*y2+2._ark*y4*y6*y8**2)*y1**2+ &
      ((sqrt(3._ark)*y5*y9**2-y6*y9**2)*y2**2-2._ark*y4**2*y6*y8**2)*y1- &
      2._ark*y2**2*y3*y6*y8**2+2._ark*y2*y3**2*y6*y8**2+(y6*y9**2- &
      sqrt(3._ark)*y5*y9**2)*y4*y3**2+(sqrt(3._ark)*y5*y9**2-y6*y9**2)*y4**2*y3
    dF(1,471) = ((y6*y7*y8/2._ark+sqrt(3._ark)*y5*y7*y8/2._ark)*y2+(-y6*y7*y9/2._ark- &
      sqrt(3._ark)*y5*y7*y9/2._ark)*y4)*y1**2+((-y6*y7*y8/2._ark-sqrt(3._ark)*y5*y7*y8/ &
      2._ark)*y2**2+(y6*y7*y9/2._ark+sqrt(3._ark)*y5*y7*y9/2._ark)*y4**2)*y1+(-y6*y7*y9/2._ark- &
      sqrt(3._ark)*y5*y7*y9/2._ark)*y3*y2**2+(y6*y7*y9/2._ark+sqrt(3._ark)*y5*y7*y9/ &
      2._ark)*y3**2*y2+(-y6*y7*y8/2._ark-sqrt(3._ark)*y5*y7*y8/2._ark)*y4*y3**2+(y6*y7*y8/ &
      2._ark+sqrt(3._ark)*y5*y7*y8/2._ark)*y4**2*y3
    dF(1,472) = ((-y6**2*y9-sqrt(3._ark)*y5*y6*y9)*y2+(-3._ark/2._ark*y5**2*y8+y6**2*y8/ &
      2._ark)*y4)*y1**2+((y6**2*y9+sqrt(3._ark)*y5*y6*y9)*y2**2+(-y6**2*y8/2._ark+3._ark/ &
      2._ark*y5**2*y8)*y4**2)*y1+(-3._ark/2._ark*y5**2*y8+y6**2*y8/2._ark)*y3*y2**2+(- &
      y6**2*y8/2._ark+3._ark/2._ark*y5**2*y8)*y3**2*y2+(y6**2*y9+ &
      sqrt(3._ark)*y5*y6*y9)*y4*y3**2+(-y6**2*y9-sqrt(3._ark)*y5*y6*y9)*y4**2*y3
    dF(1,473) = ((2._ark*sqrt(3._ark)*y5*y6*y7+y6**2*y7+3._ark*y5**2*y7)*y2+ &
      (2._ark*sqrt(3._ark)*y5*y6*y7+y6**2*y7+3._ark*y5**2*y7)*y4)*y1**2+ &
      ((2._ark*sqrt(3._ark)*y5*y6*y7+y6**2*y7+3._ark*y5**2*y7)*y2**2+ &
      (2._ark*sqrt(3._ark)*y5*y6*y7+y6**2*y7+3._ark*y5**2*y7)*y4**2)*y1+ &
      (2._ark*sqrt(3._ark)*y5*y6*y7+y6**2*y7+3._ark*y5**2*y7)*y3*y2**2+ &
      (2._ark*sqrt(3._ark)*y5*y6*y7+y6**2*y7+3._ark*y5**2*y7)*y3**2*y2+ &
      (2._ark*sqrt(3._ark)*y5*y6*y7+y6**2*y7+3._ark*y5**2*y7)*y4*y3**2+ &
      (2._ark*sqrt(3._ark)*y5*y6*y7+y6**2*y7+3._ark*y5**2*y7)*y4**2*y3
    dF(1,474) = ((-y5*y9**2+sqrt(3._ark)*y6*y9**2)*y2+(-y5*y8**2+ &
      sqrt(3._ark)*y6*y8**2)*y4)*y1**2+((y5*y9**2-sqrt(3._ark)*y6*y9**2)*y2**2+(y5*y8**2- &
      sqrt(3._ark)*y6*y8**2)*y4**2)*y1+(y5*y8**2-sqrt(3._ark)*y6*y8**2)*y3*y2**2+(- &
      y5*y8**2+sqrt(3._ark)*y6*y8**2)*y3**2*y2+(-y5*y9**2+sqrt(3._ark)*y6*y9**2)*y4*y3**2+ &
      (y5*y9**2-sqrt(3._ark)*y6*y9**2)*y4**2*y3
    dF(1,475) = (-2._ark*y2*y5*y8*y9+(sqrt(3._ark)*y6*y8*y9+y5*y8*y9)*y4)*y1**2+(- &
      2._ark*y2**2*y5*y8*y9+(sqrt(3._ark)*y6*y8*y9+y5*y8*y9)*y4**2)*y1+ &
      (sqrt(3._ark)*y6*y8*y9+y5*y8*y9)*y3*y2**2+(sqrt(3._ark)*y6*y8*y9+y5*y8*y9)*y3**2*y2- &
      2._ark*y3*y4**2*y5*y8*y9-2._ark*y3**2*y4*y5*y8*y9
    dF(1,476) = (-y3*y4*y8**2-y2*y3*y9**2)*y1**2+(y2**2*y4*y9**2+(-y3**2*y8**2+ &
      y4**2*y8**2)*y2-y3**2*y4*y9**2)*y1+y2**2*y3*y4*y8**2+y2*y3*y4**2*y9**2
    dF(1,477) = (y8**2+y9**2)*y4*y2*y1**2+((-y9**2-y8**2)*y3*y2**2+(-y9**2- &
      y8**2)*y4**2*y3)*y1+(y8**2+y9**2)*y4*y3**2*y2
    dF(1,478) = (2._ark*y2*y3*y6*y9+(-sqrt(3._ark)*y5*y8+y6*y8)*y4*y3)*y1**2+(- &
      2._ark*y2**2*y4*y6*y9+((sqrt(3._ark)*y5*y8-y6*y8)*y3**2+(sqrt(3._ark)*y5*y8- &
      y6*y8)*y4**2)*y2-2._ark*y3**2*y4*y6*y9)*y1+(-sqrt(3._ark)*y5*y8+y6*y8)*y4*y3*y2**2+ &
      2._ark*y2*y3*y4**2*y6*y9
    dF(1,479) = (-2._ark*y2**2*y6*y8+(-y6*y9+sqrt(3._ark)*y5*y9)*y4**2)*y1**2+(- &
      sqrt(3._ark)*y5*y9+y6*y9)*y3**2*y2**2+2._ark*y3**2*y4**2*y6*y8
    dF(1,480) = (-y6**2-y5**2)*y3**2*y1**2+(y6**2+y5**2)*y4**2*y2**2
    dF(1,481) = (-sqrt(3._ark)*y6*y7+y5*y7)*y3**2*y1**2+(-sqrt(3._ark)*y6*y7+ &
      y5*y7)*y4**2*y2**2
    dF(1,482) = ((y4*y7+y3*y7)*y2**2+y2*y4**2*y7+y3*y4**2*y7)*y1**2+(y2**2*y3**2*y7+ &
      y3**2*y4**2*y7)*y1+y2*y3**2*y4**2*y7+y2**2*y3**2*y4*y7
    dF(1,483) = (-2._ark*y2**2*y5+(sqrt(3._ark)*y6+y5)*y4**2)*y1**3+(2._ark*y2**3*y5+(-y5- &
      sqrt(3._ark)*y6)*y4**3)*y1**2+(-y5-sqrt(3._ark)*y6)*y3**2*y2**3+(sqrt(3._ark)*y6+ &
      y5)*y3**3*y2**2-2._ark*y3**3*y4**2*y5+2._ark*y3**2*y4**3*y5
    dF(1,484) = (-y9**3-y8**3)*y1**3+(y9**3-y8**3)*y2**3+(y8**3+y9**3)*y3**3+(- &
      y9**3+y8**3)*y4**3
    dF(1,485) = y2**3*y7**3+y4**3*y7**3+y3**3*y7**3+y1**3*y7**3
    dF(1,486) = (y8**2+y9**2)*y7*y1**3+(y8**2+y9**2)*y7*y2**3+(y8**2+ &
      y9**2)*y7*y3**3+(y8**2+y9**2)*y7*y4**3
    dF(1,487) = (y8**2*y9+y8*y9**2)*y1**3+(y8*y9**2-y8**2*y9)*y2**3+(-y8*y9**2- &
      y8**2*y9)*y3**3+(-y8*y9**2+y8**2*y9)*y4**3
    dF(1,488) = (-y8-y9)*y7**2*y1**3+(-y8+y9)*y7**2*y2**3+(y8+y9)*y7**2*y3**3+(y8- &
      y9)*y7**2*y4**3
    dF(1,489) = y4**3*y7*y8*y9-y3**3*y7*y8*y9-y1**3*y7*y8*y9+y2**3*y7*y8*y9
    dF(1,490) = (-sqrt(3._ark)*y5*y9**2+(2._ark*y8**2+y9**2)*y6)*y1**3+ &
      (sqrt(3._ark)*y5*y9**2+(-2._ark*y8**2-y9**2)*y6)*y2**3+(-sqrt(3._ark)*y5*y9**2+ &
      (2._ark*y8**2+y9**2)*y6)*y3**3+(sqrt(3._ark)*y5*y9**2+(-2._ark*y8**2-y9**2)*y6)*y4**3
    dF(1,491) = ((-sqrt(3._ark)*y8/2._ark+sqrt(3._ark)*y9/2._ark)*y7*y5+(y9/2._ark-y8/ &
      2._ark)*y7*y6)*y1**3+((sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/2._ark)*y7*y5+(y8/2._ark+y9/ &
      2._ark)*y7*y6)*y2**3+((-sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/2._ark)*y7*y5+(y8/2._ark-y9/ &
      2._ark)*y7*y6)*y3**3+((-sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y9/2._ark)*y7*y5+(-y8/2._ark-y9/ &
      2._ark)*y7*y6)*y4**3
    dF(1,492) = ((y8/2._ark+y9/2._ark)*y7*y5+(-sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y9/ &
      2._ark)*y7*y6)*y1**3+((y9/2._ark-y8/2._ark)*y7*y5+(-sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/ &
      2._ark)*y7*y6)*y2**3+((-y8/2._ark-y9/2._ark)*y7*y5+(sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/ &
      2._ark)*y7*y6)*y3**3+((y8/2._ark-y9/2._ark)*y7*y5+(-sqrt(3._ark)*y8/2._ark+sqrt(3._ark)*y9/ &
      2._ark)*y7*y6)*y4**3
    dF(1,493) = ((-y8**2+2._ark*y9**2)*y5-sqrt(3._ark)*y6*y8**2)*y1**3+((-2._ark*y9**2+ &
      y8**2)*y5+sqrt(3._ark)*y6*y8**2)*y2**3+((-y8**2+2._ark*y9**2)*y5- &
      sqrt(3._ark)*y6*y8**2)*y3**3+((-2._ark*y9**2+y8**2)*y5+sqrt(3._ark)*y6*y8**2)*y4**3
    dF(1,494) = (-y4*y7*y8-y2*y7*y9)*y1**3+(-y4**3*y7*y8-y2**3*y7*y9)*y1+ &
      y3*y4**3*y7*y9+y3**3*y4*y7*y9+y2*y3**3*y7*y8+y2**3*y3*y7*y8
    dF(1,495) = (y4*y8**2+y2*y9**2)*y1**3+(-y4**3*y8**2-y2**3*y9**2)*y1- &
      y2**3*y3*y8**2-y3*y4**3*y9**2+y2*y3**3*y8**2+y3**3*y4*y9**2
    dF(1,496) = (y2**2*y7+y4**2*y7)*y1**3+(y4**3*y7+y2**3*y7)*y1**2+y2**3*y3**2*y7+ &
      y3**3*y4**2*y7+y3**2*y4**3*y7+y2**2*y3**3*y7
    dF(1,497) = (y2**2*y9+y4**2*y8)*y1**3+(-y4**3*y8-y2**3*y9)*y1**2+y2**3*y3**2*y8+ &
      y3**2*y4**3*y9-y2**2*y3**3*y8-y3**3*y4**2*y9
    dF(1,498) = y1*y7**5+y2*y7**5+y4*y7**5+y3*y7**5
    dF(1,499) = (y8*y9**3+y8**3*y9)*y7*y1+(-y8**3*y9-y8*y9**3)*y7*y2+(y8*y9**3+ &
      y8**3*y9)*y7*y3+(-y8**3*y9-y8*y9**3)*y7*y4
    dF(1,500) = y2*y7*y8**2*y9**2+y1*y7*y8**2*y9**2+y3*y7*y8**2*y9**2+ &
      y4*y7*y8**2*y9**2
    dF(1,501) = (sqrt(3._ark)*y5*y7*y8**3+(-y8**3-2._ark*y9**3)*y7*y6)*y1+(- &
      sqrt(3._ark)*y5*y7*y8**3+(y8**3-2._ark*y9**3)*y7*y6)*y2+(-sqrt(3._ark)*y5*y7*y8**3+ &
      (y8**3+2._ark*y9**3)*y7*y6)*y3+(sqrt(3._ark)*y5*y7*y8**3+(2._ark*y9**3- &
      y8**3)*y7*y6)*y4
    dF(1,502) = (y5**2*y7**3+y6**2*y7**3)*y1+(y5**2*y7**3+y6**2*y7**3)*y2+ &
      (y5**2*y7**3+y6**2*y7**3)*y3+(y5**2*y7**3+y6**2*y7**3)*y4
    dF(1,503) = ((y8**2*y9+y8*y9**2)*y5**2+(y8**2*y9+y8*y9**2)*y6**2)*y1+((y8*y9**2- &
      y8**2*y9)*y5**2+(y8*y9**2-y8**2*y9)*y6**2)*y2+((-y8*y9**2-y8**2*y9)*y5**2+(- &
      y8*y9**2-y8**2*y9)*y6**2)*y3+((-y8*y9**2+y8**2*y9)*y5**2+(-y8*y9**2+ &
      y8**2*y9)*y6**2)*y4
    dF(1,504) = (y5**2*y7*y8*y9+y6**2*y7*y8*y9)*y1+(-y5**2*y7*y8*y9- &
      y6**2*y7*y8*y9)*y2+(y5**2*y7*y8*y9+y6**2*y7*y8*y9)*y3+(-y5**2*y7*y8*y9- &
      y6**2*y7*y8*y9)*y4
    dF(1,505) = (4._ark/9._ark*sqrt(3._ark)*y5**3*y8*y9-y6**3*y8*y9-y5**2*y6*y8*y9)*y1+ &
      (4._ark/9._ark*sqrt(3._ark)*y5**3*y8*y9-y6**3*y8*y9-y5**2*y6*y8*y9)*y2+(4._ark/ &
      9._ark*sqrt(3._ark)*y5**3*y8*y9-y6**3*y8*y9-y5**2*y6*y8*y9)*y3+(4._ark/ &
      9._ark*sqrt(3._ark)*y5**3*y8*y9-y6**3*y8*y9-y5**2*y6*y8*y9)*y4
    dF(1,506) = (-3._ark/4._ark*sqrt(3._ark)*y5**3*y8**2+9._ark/4._ark*y5**2*y6*y9**2-3._ark/ &
      4._ark*sqrt(3._ark)*y5*y6**2*y8**2+(y8**2+5._ark/4._ark*y9**2)*y6**3)*y1+(3._ark/ &
      4._ark*sqrt(3._ark)*y5**3*y8**2-9._ark/4._ark*y5**2*y6*y9**2+3._ark/ &
      4._ark*sqrt(3._ark)*y5*y6**2*y8**2+(-5._ark/4._ark*y9**2-y8**2)*y6**3)*y2+(-3._ark/ &
      4._ark*sqrt(3._ark)*y5**3*y8**2+9._ark/4._ark*y5**2*y6*y9**2-3._ark/ &
      4._ark*sqrt(3._ark)*y5*y6**2*y8**2+(y8**2+5._ark/4._ark*y9**2)*y6**3)*y3+(3._ark/ &
      4._ark*sqrt(3._ark)*y5**3*y8**2-9._ark/4._ark*y5**2*y6*y9**2+3._ark/ &
      4._ark*sqrt(3._ark)*y5*y6**2*y8**2+(-5._ark/4._ark*y9**2-y8**2)*y6**3)*y4
    dF(1,507) = ((-2._ark*y9**2+y8**2)*y7**2*y5+sqrt(3._ark)*y6*y7**2*y8**2)*y1+((- &
      y8**2+2._ark*y9**2)*y7**2*y5-sqrt(3._ark)*y6*y7**2*y8**2)*y2+((-2._ark*y9**2+ &
      y8**2)*y7**2*y5+sqrt(3._ark)*y6*y7**2*y8**2)*y3+((-y8**2+2._ark*y9**2)*y7**2*y5- &
      sqrt(3._ark)*y6*y7**2*y8**2)*y4
    dF(1,508) = ((-2._ark*y8**3+y9**3)*y7*y5+sqrt(3._ark)*y6*y7*y9**3)*y1+((2._ark*y8**3+ &
      y9**3)*y7*y5+sqrt(3._ark)*y6*y7*y9**3)*y2+((2._ark*y8**3-y9**3)*y7*y5- &
      sqrt(3._ark)*y6*y7*y9**3)*y3+((-y9**3-2._ark*y8**3)*y7*y5-sqrt(3._ark)*y6*y7*y9**3)*y4
    dF(1,509) = (-sqrt(3._ark)*y6*y8**2*y9**2+y5*y8**2*y9**2)*y1+ &
      (sqrt(3._ark)*y6*y8**2*y9**2-y5*y8**2*y9**2)*y2+(-sqrt(3._ark)*y6*y8**2*y9**2+ &
      y5*y8**2*y9**2)*y3+(sqrt(3._ark)*y6*y8**2*y9**2-y5*y8**2*y9**2)*y4
    dF(1,510) = ((y8**2*y9-2._ark*y8*y9**2)*y7*y5+sqrt(3._ark)*y6*y7*y8**2*y9)*y1+ &
      ((y8**2*y9+2._ark*y8*y9**2)*y7*y5+sqrt(3._ark)*y6*y7*y8**2*y9)*y2+((-y8**2*y9+ &
      2._ark*y8*y9**2)*y7*y5-sqrt(3._ark)*y6*y7*y8**2*y9)*y3+((-y8**2*y9- &
      2._ark*y8*y9**2)*y7*y5-sqrt(3._ark)*y6*y7*y8**2*y9)*y4
    dF(1,511) = (3._ark/4._ark*y5**3*y7*y8-5._ark/4._ark*sqrt(3._ark)*y5**2*y6*y7*y8+(y9+7._ark/ &
      4._ark*y8)*y7*y6**2*y5+(-sqrt(3._ark)*y9-sqrt(3._ark)*y8/4._ark)*y7*y6**3)*y1+(-3._ark/ &
      4._ark*y5**3*y7*y8+5._ark/4._ark*sqrt(3._ark)*y5**2*y6*y7*y8+(-7._ark/4._ark*y8+ &
      y9)*y7*y6**2*y5+(-sqrt(3._ark)*y9+sqrt(3._ark)*y8/4._ark)*y7*y6**3)*y2+(-3._ark/ &
      4._ark*y5**3*y7*y8+5._ark/4._ark*sqrt(3._ark)*y5**2*y6*y7*y8+(-7._ark/4._ark*y8- &
      y9)*y7*y6**2*y5+(sqrt(3._ark)*y8/4._ark+sqrt(3._ark)*y9)*y7*y6**3)*y3+(3._ark/ &
      4._ark*y5**3*y7*y8-5._ark/4._ark*sqrt(3._ark)*y5**2*y6*y7*y8+(-y9+7._ark/ &
      4._ark*y8)*y7*y6**2*y5+(sqrt(3._ark)*y9-sqrt(3._ark)*y8/4._ark)*y7*y6**3)*y4
    dF(1,512) = (-3._ark/2._ark*sqrt(3._ark)*y5**5+17._ark/6._ark*sqrt(3._ark)*y5*y6**4+ &
      2._ark*y6**5+3._ark*y5**4*y6+y5**2*y6**3)*y1+(-y5**2*y6**3-17._ark/ &
      6._ark*sqrt(3._ark)*y5*y6**4-3._ark*y5**4*y6+3._ark/2._ark*sqrt(3._ark)*y5**5- &
      2._ark*y6**5)*y2+(-3._ark/2._ark*sqrt(3._ark)*y5**5+17._ark/6._ark*sqrt(3._ark)*y5*y6**4+ &
      2._ark*y6**5+3._ark*y5**4*y6+y5**2*y6**3)*y3+(-y5**2*y6**3-17._ark/ &
      6._ark*sqrt(3._ark)*y5*y6**4-3._ark*y5**4*y6+3._ark/2._ark*sqrt(3._ark)*y5**5-2._ark*y6**5)*y4
    dF(1,513) = ((13._ark/4._ark*y8+y9)*y7*y5**3-15._ark/4._ark*sqrt(3._ark)*y5**2*y6*y7*y8+ &
      9._ark/4._ark*y5*y6**2*y7*y8+(-3._ark/4._ark*sqrt(3._ark)*y8- &
      3._ark*sqrt(3._ark)*y9)*y7*y6**3)*y1+((-13._ark/4._ark*y8+y9)*y7*y5**3+15._ark/ &
      4._ark*sqrt(3._ark)*y5**2*y6*y7*y8-9._ark/4._ark*y5*y6**2*y7*y8+(-3._ark*sqrt(3._ark)*y9+ &
      3._ark/4._ark*sqrt(3._ark)*y8)*y7*y6**3)*y2+((-13._ark/4._ark*y8-y9)*y7*y5**3+15._ark/ &
      4._ark*sqrt(3._ark)*y5**2*y6*y7*y8-9._ark/4._ark*y5*y6**2*y7*y8+(3._ark*sqrt(3._ark)*y9+ &
      3._ark/4._ark*sqrt(3._ark)*y8)*y7*y6**3)*y3+((-y9+13._ark/4._ark*y8)*y7*y5**3-15._ark/ &
      4._ark*sqrt(3._ark)*y5**2*y6*y7*y8+9._ark/4._ark*y5*y6**2*y7*y8+(-3._ark/ &
      4._ark*sqrt(3._ark)*y8+3._ark*sqrt(3._ark)*y9)*y7*y6**3)*y4
    dF(1,514) = (y4*y7**2*y8*y9+y2*y7**2*y8*y9)*y1+y3*y4*y7**2*y8*y9+ &
      y2*y3*y7**2*y8*y9
    dF(1,515) = ((-3._ark/4._ark*sqrt(3._ark)*y5*y6**2*y7-27._ark/20._ark*y5**2*y6*y7-3._ark/ &
      20._ark*sqrt(3._ark)*y5**3*y7-7._ark/20._ark*y6**3*y7)*y2+(9._ark/10._ark*y5**2*y6*y7- &
      y6**3*y7/10._ark+3._ark/5._ark*sqrt(3._ark)*y5**3*y7)*y4)*y1+(9._ark/10._ark*y5**2*y6*y7- &
      y6**3*y7/10._ark+3._ark/5._ark*sqrt(3._ark)*y5**3*y7)*y3*y2+(-3._ark/ &
      4._ark*sqrt(3._ark)*y5*y6**2*y7-27._ark/20._ark*y5**2*y6*y7-3._ark/ &
      20._ark*sqrt(3._ark)*y5**3*y7-7._ark/20._ark*y6**3*y7)*y4*y3
    dF(1,516) = ((sqrt(3._ark)*y5**2*y6*y7/5._ark-3._ark/5._ark*y5**3*y7+ &
      sqrt(3._ark)*y6**3*y7/5._ark+y5*y6**2*y7)*y2+(sqrt(3._ark)*y5**2*y6*y7/5._ark-3._ark/ &
      5._ark*y5**3*y7+sqrt(3._ark)*y6**3*y7/5._ark+y5*y6**2*y7)*y4)*y1+ &
      (sqrt(3._ark)*y5**2*y6*y7/5._ark-3._ark/5._ark*y5**3*y7+sqrt(3._ark)*y6**3*y7/5._ark+ &
      y5*y6**2*y7)*y3*y2+(sqrt(3._ark)*y5**2*y6*y7/5._ark-3._ark/5._ark*y5**3*y7+ &
      sqrt(3._ark)*y6**3*y7/5._ark+y5*y6**2*y7)*y4*y3
    dF(1,517) = ((-3._ark/4._ark*y5**3*y8+y5*y6**2*y8/4._ark)*y2+(sqrt(3._ark)*y5**2*y6*y9/ &
      4._ark+sqrt(3._ark)*y6**3*y9/4._ark+y9*y6**2*y5)*y4)*y1+(-y9*y6**2*y5- &
      sqrt(3._ark)*y6**3*y9/4._ark-sqrt(3._ark)*y5**2*y6*y9/4._ark)*y3*y2+(3._ark/4._ark*y5**3*y8- &
      y5*y6**2*y8/4._ark)*y4*y3
    dF(1,518) = ((y6**2*y8*y9+y5**2*y8*y9)*y2+(y6**2*y8*y9+y5**2*y8*y9)*y4)*y1+ &
      (y6**2*y8*y9+y5**2*y8*y9)*y3*y2+(y6**2*y8*y9+y5**2*y8*y9)*y4*y3
    dF(1,519) = -y2**5*y4+y1*y3**5-y2*y4**5+y1**5*y3
    dF(1,520) = -y2*y4*y8**2*y9**2+y1*y3*y8**2*y9**2
    dF(1,521) = (y4*y8*y9**3+y2*y8**3*y9)*y1+y3*y4*y8**3*y9+y2*y3*y8*y9**3
    dF(1,522) = (y4*y7**3*y8+y2*y7**3*y9)*y1-y3*y4*y7**3*y9-y2*y3*y7**3*y8
    dF(1,523) = (y8**2+y9**2)*y7**2*y3*y1+(-y9**2-y8**2)*y7**2*y4*y2
    dF(1,524) = (-sqrt(3._ark)*y5**3*y7/9._ark+y5**2*y6*y7-sqrt(3._ark)*y5*y6**2*y7+ &
      y6**3*y7)*y3*y1+(-sqrt(3._ark)*y5**3*y7/9._ark+y5**2*y6*y7-sqrt(3._ark)*y5*y6**2*y7+ &
      y6**3*y7)*y4*y2
    dF(1,525) = (-sqrt(3._ark)*y6**2*y8*y9/2._ark+y5*y6*y8*y9-sqrt(3._ark)*y5**2*y8*y9/ &
      6._ark)*y3*y1+(-sqrt(3._ark)*y6**2*y8*y9/2._ark+y5*y6*y8*y9-sqrt(3._ark)*y5**2*y8*y9/ &
      6._ark)*y4*y2
    dF(1,526) = ((sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y9**2/2._ark)*y5**2+(y8**2+ &
      y9**2)*y6*y5+(sqrt(3._ark)*y9**2/6._ark+sqrt(3._ark)*y8**2/6._ark)*y6**2)*y3*y1+((- &
      sqrt(3._ark)*y9**2/2._ark-sqrt(3._ark)*y8**2/2._ark)*y5**2+(-y9**2-y8**2)*y6*y5+(- &
      sqrt(3._ark)*y9**2/6._ark-sqrt(3._ark)*y8**2/6._ark)*y6**2)*y4*y2
    dF(1,527) = (-2._ark/3._ark*sqrt(3._ark)*y2*y5**2*y7*y9+(-sqrt(3._ark)*y5**2*y7*y8/6._ark- &
      sqrt(3._ark)*y6**2*y7*y8/2._ark-y5*y6*y7*y8)*y4)*y1+(sqrt(3._ark)*y5**2*y7*y8/6._ark+ &
      y5*y6*y7*y8+sqrt(3._ark)*y6**2*y7*y8/2._ark)*y3*y2+2._ark/ &
      3._ark*sqrt(3._ark)*y3*y4*y5**2*y7*y9
    dF(1,528) = ((sqrt(3._ark)*y6**2*y8*y9/6._ark-sqrt(3._ark)*y5**2*y8*y9/2._ark)*y2+(- &
      sqrt(3._ark)*y6**2*y8*y9/3._ark-y5*y6*y8*y9)*y4)*y1+(-sqrt(3._ark)*y6**2*y8*y9/3._ark- &
      y5*y6*y8*y9)*y3*y2+(sqrt(3._ark)*y6**2*y8*y9/6._ark-sqrt(3._ark)*y5**2*y8*y9/ &
      2._ark)*y4*y3
    dF(1,529) = ((sqrt(3._ark)*y5**3*y7/20._ark+sqrt(3._ark)*y5*y6**2*y7/4._ark+9._ark/ &
      20._ark*y6**3*y7-11._ark/20._ark*y5**2*y6*y7)*y2+(-sqrt(3._ark)*y5**3*y7/5._ark-3._ark/ &
      10._ark*y6**3*y7+7._ark/10._ark*y5**2*y6*y7)*y4)*y1+(-sqrt(3._ark)*y5**3*y7/5._ark-3._ark/ &
      10._ark*y6**3*y7+7._ark/10._ark*y5**2*y6*y7)*y3*y2+(sqrt(3._ark)*y5**3*y7/20._ark+ &
      sqrt(3._ark)*y5*y6**2*y7/4._ark+9._ark/20._ark*y6**3*y7-11._ark/20._ark*y5**2*y6*y7)*y4*y3
    dF(1,530) = (3._ark*y5*y6**3-3._ark/4._ark*sqrt(3._ark)*y6**4-sqrt(3._ark)*y5**4/12._ark- &
      3._ark/2._ark*sqrt(3._ark)*y5**2*y6**2+y5**3*y6)*y3*y1+(-y5**3*y6+3._ark/ &
      4._ark*sqrt(3._ark)*y6**4+3._ark/2._ark*sqrt(3._ark)*y5**2*y6**2-3._ark*y5*y6**3+ &
      sqrt(3._ark)*y5**4/12._ark)*y4*y2
    dF(1,531) = (sqrt(3._ark)*y5*y7*y8+(-y8-2._ark*y9)*y7*y6)*y3*y1**2+(- &
      sqrt(3._ark)*y5*y7*y8+(y8+2._ark*y9)*y7*y6)*y3**2*y1+(-sqrt(3._ark)*y5*y7*y8+(y8- &
      2._ark*y9)*y7*y6)*y4*y2**2+(sqrt(3._ark)*y5*y7*y8+(-y8+2._ark*y9)*y7*y6)*y4**2*y2
    dF(1,532) = ((sqrt(3._ark)*y5*y7*y9+y6*y7*y9)*y2+(-sqrt(3._ark)*y5*y7*y8- &
      y6*y7*y8)*y4)*y1**2+((sqrt(3._ark)*y5*y7*y9+y6*y7*y9)*y2**2+(-sqrt(3._ark)*y5*y7*y8- &
      y6*y7*y8)*y4**2)*y1+(y6*y7*y8+sqrt(3._ark)*y5*y7*y8)*y3*y2**2+(y6*y7*y8+ &
      sqrt(3._ark)*y5*y7*y8)*y3**2*y2+(-sqrt(3._ark)*y5*y7*y9-y6*y7*y9)*y4*y3**2+(- &
      sqrt(3._ark)*y5*y7*y9-y6*y7*y9)*y4**2*y3
    dF(1,533) = ((-5._ark/9._ark*sqrt(3._ark)*y5*y6**2-16._ark/9._ark*y5**2*y6- &
      sqrt(3._ark)*y5**3/9._ark)*y2+(13._ark/9._ark*y5**2*y6-y6**3/3._ark+4._ark/ &
      9._ark*sqrt(3._ark)*y5**3)*y4)*y1**2+((16._ark/9._ark*y5**2*y6+sqrt(3._ark)*y5**3/9._ark+ &
      5._ark/9._ark*sqrt(3._ark)*y5*y6**2)*y2**2+(-13._ark/9._ark*y5**2*y6-4._ark/ &
      9._ark*sqrt(3._ark)*y5**3+y6**3/3._ark)*y4**2)*y1+(-13._ark/9._ark*y5**2*y6-4._ark/ &
      9._ark*sqrt(3._ark)*y5**3+y6**3/3._ark)*y3*y2**2+(13._ark/9._ark*y5**2*y6-y6**3/3._ark+4._ark/ &
      9._ark*sqrt(3._ark)*y5**3)*y3**2*y2+(-5._ark/9._ark*sqrt(3._ark)*y5*y6**2-16._ark/ &
      9._ark*y5**2*y6-sqrt(3._ark)*y5**3/9._ark)*y4*y3**2+(16._ark/9._ark*y5**2*y6+ &
      sqrt(3._ark)*y5**3/9._ark+5._ark/9._ark*sqrt(3._ark)*y5*y6**2)*y4**2*y3
    dF(1,534) = ((2._ark*y8-y9)*y7*y5-sqrt(3._ark)*y6*y7*y9)*y3*y1**2+((-2._ark*y8+ &
      y9)*y7*y5+sqrt(3._ark)*y6*y7*y9)*y3**2*y1+((-y9-2._ark*y8)*y7*y5- &
      sqrt(3._ark)*y6*y7*y9)*y4*y2**2+((y9+2._ark*y8)*y7*y5+sqrt(3._ark)*y6*y7*y9)*y4**2*y2
    dF(1,535) = (sqrt(3._ark)*y6**2*y7/6._ark+sqrt(3._ark)*y5**2*y7/2._ark+ &
      y5*y6*y7)*y3*y1**2+(sqrt(3._ark)*y6**2*y7/6._ark+sqrt(3._ark)*y5**2*y7/2._ark+ &
      y5*y6*y7)*y3**2*y1+(sqrt(3._ark)*y6**2*y7/6._ark+sqrt(3._ark)*y5**2*y7/2._ark+ &
      y5*y6*y7)*y4*y2**2+(sqrt(3._ark)*y6**2*y7/6._ark+sqrt(3._ark)*y5**2*y7/2._ark+ &
      y5*y6*y7)*y4**2*y2
    dF(1,536) = ((y5*y6*y9/2._ark-sqrt(3._ark)*y6**2*y9/2._ark)*y2+(-sqrt(3._ark)*y6**2*y8/ &
      4._ark+y5*y6*y8-sqrt(3._ark)*y5**2*y8/4._ark)*y4)*y1**2+((sqrt(3._ark)*y6**2*y9/2._ark- &
      y5*y6*y9/2._ark)*y2**2+(sqrt(3._ark)*y5**2*y8/4._ark-y5*y6*y8+sqrt(3._ark)*y6**2*y8/ &
      4._ark)*y4**2)*y1+(-sqrt(3._ark)*y6**2*y8/4._ark+y5*y6*y8-sqrt(3._ark)*y5**2*y8/ &
      4._ark)*y3*y2**2+(sqrt(3._ark)*y5**2*y8/4._ark-y5*y6*y8+sqrt(3._ark)*y6**2*y8/ &
      4._ark)*y3**2*y2+(sqrt(3._ark)*y6**2*y9/2._ark-y5*y6*y9/2._ark)*y4*y3**2+(y5*y6*y9/2._ark- &
      sqrt(3._ark)*y6**2*y9/2._ark)*y4**2*y3
    dF(1,537) = (y3*y4*y8*y9+y2*y3*y8*y9)*y1**2+(y2**2*y4*y8*y9+(y4**2*y8*y9+ &
      y3**2*y8*y9)*y2+y3**2*y4*y8*y9)*y1+y2*y3*y4**2*y8*y9+y2**2*y3*y4*y8*y9
    dF(1,538) = (-y8-y9)*y7*y4*y2*y1**2+((y8-y9)*y7*y3*y2**2+(-y8+ &
      y9)*y7*y4**2*y3)*y1+(y8+y9)*y7*y4*y3**2*y2
    dF(1,539) = -y1**2*y2*y4*y7**2+(y2**2*y3*y7**2+y3*y4**2*y7**2)*y1- &
      y2*y3**2*y4*y7**2
    dF(1,540) = (-y3*y4*y9**2-y2*y3*y8**2)*y1**2+(y2**2*y4*y8**2+(-y3**2*y9**2+ &
      y4**2*y9**2)*y2-y3**2*y4*y8**2)*y1+y2*y3*y4**2*y8**2+y2**2*y3*y4*y9**2
    dF(1,541) = (y4*y7*y8*y9**2+y2*y7*y8**2*y9)*y1-y3*y4*y7*y8**2*y9- &
      y2*y3*y7*y8*y9**2
    dF(1,542) = (y8*y9**3+y8**3*y9)*y3*y1+(y8*y9**3+y8**3*y9)*y4*y2
    dF(1,543) = (y9**4+y8**4)*y3*y1+(-y8**4-y9**4)*y4*y2
    dF(1,544) = (-sqrt(3._ark)*y5*y7*y8**2+(y8**2+2._ark*y9**2)*y7*y6)*y3*y1+(- &
      sqrt(3._ark)*y5*y7*y8**2+(y8**2+2._ark*y9**2)*y7*y6)*y4*y2
    dF(1,545) = (y2*y6*y7**2*y8+(-sqrt(3._ark)*y5*y7**2*y9/2._ark+y6*y7**2*y9/ &
      2._ark)*y4)*y1+(-y6*y7**2*y9/2._ark+sqrt(3._ark)*y5*y7**2*y9/2._ark)*y3*y2- &
      y3*y4*y6*y7**2*y8
    dF(1,546) = ((-3._ark/4._ark*sqrt(3._ark)*y5**3*y8-y6**3*y8-3._ark/ &
      4._ark*sqrt(3._ark)*y5*y6**2*y8)*y2+(3._ark/4._ark*sqrt(3._ark)*y5*y6**2*y9+y6**3*y9+3._ark/ &
      4._ark*sqrt(3._ark)*y5**3*y9)*y4)*y1+(-3._ark/4._ark*sqrt(3._ark)*y5**3*y9-y6**3*y9-3._ark/ &
      4._ark*sqrt(3._ark)*y5*y6**2*y9)*y3*y2+(3._ark/4._ark*sqrt(3._ark)*y5*y6**2*y8+3._ark/ &
      4._ark*sqrt(3._ark)*y5**3*y8+y6**3*y8)*y4*y3
    dF(1,547) = (-2._ark*sqrt(3._ark)*y5*y6**3+3._ark/4._ark*y5**4+9._ark/2._ark*y5**2*y6**2+ &
      7._ark/4._ark*y6**4)*y3*y1+(-3._ark/4._ark*y5**4-9._ark/2._ark*y5**2*y6**2-7._ark/4._ark*y6**4+ &
      2._ark*sqrt(3._ark)*y5*y6**3)*y4*y2
    dF(1,548) = (sqrt(3._ark)*y6*y7**3/2._ark-y5*y7**3/2._ark)*y3*y1+(sqrt(3._ark)*y6*y7**3/ &
      2._ark-y5*y7**3/2._ark)*y4*y2
    dF(1,549) = (y2*y5*y8*y9**2+(-sqrt(3._ark)*y6*y8**2*y9/2._ark-y5*y8**2*y9/ &
      2._ark)*y4)*y1+(y5*y8**2*y9/2._ark+sqrt(3._ark)*y6*y8**2*y9/2._ark)*y3*y2- &
      y3*y4*y5*y8*y9**2
    dF(1,550) = (-y5*y7*y8*y9/2._ark+sqrt(3._ark)*y6*y7*y8*y9/2._ark)*y3*y1+(y5*y7*y8*y9/ &
      2._ark-sqrt(3._ark)*y6*y7*y8*y9/2._ark)*y4*y2
    dF(1,551) = ((-y5*y6*y8*y9-sqrt(3._ark)*y6**2*y8*y9/6._ark-sqrt(3._ark)*y5**2*y8*y9/ &
      2._ark)*y2+(-y5*y6*y8*y9-sqrt(3._ark)*y6**2*y8*y9/6._ark-sqrt(3._ark)*y5**2*y8*y9/ &
      2._ark)*y4)*y1+(-y5*y6*y8*y9-sqrt(3._ark)*y6**2*y8*y9/6._ark-sqrt(3._ark)*y5**2*y8*y9/ &
      2._ark)*y3*y2+(-y5*y6*y8*y9-sqrt(3._ark)*y6**2*y8*y9/6._ark-sqrt(3._ark)*y5**2*y8*y9/ &
      2._ark)*y4*y3
    dF(1,552) = (-y5**3*y7/3._ark+y5*y6**2*y7)*y3*y1+(-y5**3*y7/3._ark+ &
      y5*y6**2*y7)*y4*y2
    dF(1,553) = (y5**2*y7**2+y6**2*y7**2)*y3*y1+(-y5**2*y7**2-y6**2*y7**2)*y4*y2
    dF(1,554) = ((y8**2+y9**2)*y5**2+(y8**2+y9**2)*y6**2)*y3*y1+((-y9**2- &
      y8**2)*y5**2+(-y9**2-y8**2)*y6**2)*y4*y2
    dF(1,555) = ((-sqrt(3._ark)*y5**3*y8/4._ark-y5**2*y6*y8-sqrt(3._ark)*y5*y6**2*y8/ &
      4._ark)*y2+(y5**2*y6*y9+sqrt(3._ark)*y5*y6**2*y9/4._ark+sqrt(3._ark)*y5**3*y9/ &
      4._ark)*y4)*y1+(-y5**2*y6*y9-sqrt(3._ark)*y5**3*y9/4._ark-sqrt(3._ark)*y5*y6**2*y9/ &
      4._ark)*y3*y2+(sqrt(3._ark)*y5*y6**2*y8/4._ark+sqrt(3._ark)*y5**3*y8/4._ark+ &
      y5**2*y6*y8)*y4*y3
    dF(1,556) = ((y4*y7**3+y3*y7**3)*y2+y3*y4*y7**3)*y1+y2*y3*y4*y7**3
    dF(1,557) = (((-sqrt(3._ark)*y5*y7*y9/2._ark+(y9/2._ark-y8)*y7*y6)*y3+(- &
      sqrt(3._ark)*y5*y7*y9/2._ark+(y8+y9/2._ark)*y7*y6)*y4)*y2+(sqrt(3._ark)*y5*y7*y9/2._ark+ &
      (y8-y9/2._ark)*y7*y6)*y4*y3)*y1+(sqrt(3._ark)*y5*y7*y9/2._ark+(-y9/2._ark- &
      y8)*y7*y6)*y4*y3*y2
    dF(1,558) = ((((-2._ark*y8**2+y9**2)*y5+sqrt(3._ark)*y6*y9**2)*y3+((-y9**2+ &
      2._ark*y8**2)*y5-sqrt(3._ark)*y6*y9**2)*y4)*y2+((-2._ark*y8**2+y9**2)*y5+ &
      sqrt(3._ark)*y6*y9**2)*y4*y3)*y1+((-y9**2+2._ark*y8**2)*y5- &
      sqrt(3._ark)*y6*y9**2)*y4*y3*y2
    dF(1,559) = ((((-y9/2._ark-y8)*y7*y5-sqrt(3._ark)*y6*y7*y9/2._ark)*y3+((y8-y9/ &
      2._ark)*y7*y5-sqrt(3._ark)*y6*y7*y9/2._ark)*y4)*y2+((y8+y9/2._ark)*y7*y5+ &
      sqrt(3._ark)*y6*y7*y9/2._ark)*y4*y3)*y1+((y9/2._ark-y8)*y7*y5+sqrt(3._ark)*y6*y7*y9/ &
      2._ark)*y4*y3*y2
    dF(1,560) = (((sqrt(3._ark)*y5**2*y9/2._ark+(y8-2._ark*y9)*y6*y5+(-sqrt(3._ark)*y8+ &
      sqrt(3._ark)*y9/2._ark)*y6**2)*y3+(-sqrt(3._ark)*y5**2*y9/2._ark+(y8+2._ark*y9)*y6*y5+(- &
      sqrt(3._ark)*y8-sqrt(3._ark)*y9/2._ark)*y6**2)*y4)*y2+(-sqrt(3._ark)*y5**2*y9/2._ark+(-y8+ &
      2._ark*y9)*y6*y5+(sqrt(3._ark)*y8-sqrt(3._ark)*y9/2._ark)*y6**2)*y4*y3)*y1+ &
      (sqrt(3._ark)*y5**2*y9/2._ark+(-y8-2._ark*y9)*y6*y5+(sqrt(3._ark)*y9/2._ark+ &
      sqrt(3._ark)*y8)*y6**2)*y4*y3*y2
    dF(1,561) = (((y5*y6**2-y5**3/3._ark)*y3+(y5**3/3._ark-y5*y6**2)*y4)*y2+(y5*y6**2- &
      y5**3/3._ark)*y4*y3)*y1+(y5**3/3._ark-y5*y6**2)*y4*y3*y2
    dF(1,562) = (-y3*y4*y9-y2*y3*y8)*y1**3+(-y2**3*y4*y8+(-y4**3*y9+y3**3*y9)*y2+ &
      y3**3*y4*y8)*y1+y2**3*y3*y4*y9+y2*y3*y4**3*y8
    dF(1,563) = (-sqrt(3._ark)*y6/2._ark+y5/2._ark)*y4*y2*y1**3+((sqrt(3._ark)*y6/2._ark-y5/ &
      2._ark)*y3*y2**3+(sqrt(3._ark)*y6/2._ark-y5/2._ark)*y4**3*y3)*y1+(-sqrt(3._ark)*y6/2._ark+ &
      y5/2._ark)*y4*y3**3*y2
    dF(1,564) = (((-y9**3+y8**3)*y3+(y8**3+y9**3)*y4)*y2+(y9**3-y8**3)*y4*y3)*y1+(- &
      y9**3-y8**3)*y4*y3*y2
    dF(1,565) = ((-y4*y7*y8*y9+y3*y7*y8*y9)*y2+y3*y4*y7*y8*y9)*y1-y2*y3*y4*y7*y8*y9
    dF(1,566) = (((y6**2*y7+y5**2*y7)*y3+(y6**2*y7+y5**2*y7)*y4)*y2+(y6**2*y7+ &
      y5**2*y7)*y4*y3)*y1+(y6**2*y7+y5**2*y7)*y4*y3*y2
    dF(1,567) = ((((-2._ark*sqrt(3._ark)*y8+2._ark*sqrt(3._ark)*y9)*y6*y5+(-2._ark*y9+ &
      2._ark*y8)*y6**2)*y3+((-2._ark*sqrt(3._ark)*y9-2._ark*sqrt(3._ark)*y8)*y6*y5+(2._ark*y8+ &
      2._ark*y9)*y6**2)*y4)*y2+((2._ark*sqrt(3._ark)*y8-2._ark*sqrt(3._ark)*y9)*y6*y5+(-2._ark*y8+ &
      2._ark*y9)*y6**2)*y4*y3)*y1+((2._ark*sqrt(3._ark)*y8+2._ark*sqrt(3._ark)*y9)*y6*y5+(- &
      2._ark*y8-2._ark*y9)*y6**2)*y4*y3*y2
    dF(1,568) = (((y5**2*y6-sqrt(3._ark)*y5**3/9._ark-sqrt(3._ark)*y5*y6**2+y6**3)*y3+ &
      (sqrt(3._ark)*y5*y6**2+sqrt(3._ark)*y5**3/9._ark-y6**3-y5**2*y6)*y4)*y2+(y5**2*y6- &
      sqrt(3._ark)*y5**3/9._ark-sqrt(3._ark)*y5*y6**2+y6**3)*y4*y3)*y1+(sqrt(3._ark)*y5*y6**2+ &
      sqrt(3._ark)*y5**3/9._ark-y6**3-y5**2*y6)*y4*y3*y2
    dF(1,569) = (((-y5*y8*y9/2._ark+sqrt(3._ark)*y6*y8*y9/2._ark)*y3+(-y5*y8*y9/2._ark+ &
      sqrt(3._ark)*y6*y8*y9/2._ark)*y4)*y2+(-y5*y8*y9/2._ark+sqrt(3._ark)*y6*y8*y9/ &
      2._ark)*y4*y3)*y1+(-y5*y8*y9/2._ark+sqrt(3._ark)*y6*y8*y9/2._ark)*y4*y3*y2
    dF(1,570) = (((-sqrt(3._ark)*y6*y7**2+y5*y7**2)*y3+(-y5*y7**2+ &
      sqrt(3._ark)*y6*y7**2)*y4)*y2+(-sqrt(3._ark)*y6*y7**2+y5*y7**2)*y4*y3)*y1+(- &
      y5*y7**2+sqrt(3._ark)*y6*y7**2)*y4*y3*y2
    dF(1,571) = (((-sqrt(3._ark)*y5**2*y7/6._ark-sqrt(3._ark)*y6**2*y7/2._ark+y5*y6*y7)*y3+ &
      (-sqrt(3._ark)*y5**2*y7/6._ark-sqrt(3._ark)*y6**2*y7/2._ark+y5*y6*y7)*y4)*y2+(- &
      sqrt(3._ark)*y5**2*y7/6._ark-sqrt(3._ark)*y6**2*y7/2._ark+y5*y6*y7)*y4*y3)*y1+(- &
      sqrt(3._ark)*y5**2*y7/6._ark-sqrt(3._ark)*y6**2*y7/2._ark+y5*y6*y7)*y4*y3*y2
    dF(1,572) = ((((y8-y9)*y5**2+(-2._ark*sqrt(3._ark)*y8+2._ark*sqrt(3._ark)*y9)*y6*y5+ &
      (3._ark*y8-3._ark*y9)*y6**2)*y3+((y8+y9)*y5**2+(-2._ark*sqrt(3._ark)*y9- &
      2._ark*sqrt(3._ark)*y8)*y6*y5+(3._ark*y8+3._ark*y9)*y6**2)*y4)*y2+((-y8+y9)*y5**2+ &
      (2._ark*sqrt(3._ark)*y8-2._ark*sqrt(3._ark)*y9)*y6*y5+(-3._ark*y8+ &
      3._ark*y9)*y6**2)*y4*y3)*y1+((-y8-y9)*y5**2+(2._ark*sqrt(3._ark)*y8+ &
      2._ark*sqrt(3._ark)*y9)*y6*y5+(-3._ark*y8-3._ark*y9)*y6**2)*y4*y3*y2
    dF(1,573) = y1*y2*y3*y4*y8*y9
    dF(1,574) = (y6*y7-sqrt(3._ark)*y5*y7/3._ark)*y4*y3*y2*y1
    dF(1,575) = (y3*y4*y7+y2*y3*y7)*y1**3+(y2**3*y4*y7+(y3**3*y7+y4**3*y7)*y2+ &
      y3**3*y4*y7)*y1+y2*y3*y4**3*y7+y2**3*y3*y4*y7
    dF(1,576) = (-y2*y3-y3*y4)*y1**4+(y2**4*y4+(y4**4-y3**4)*y2-y3**4*y4)*y1+ &
      y2*y3*y4**4+y2**4*y3*y4
    dF(1,577) = -y1**4*y2*y4+(y3*y4**4+y2**4*y3)*y1-y2*y3**4*y4
    dF(1,578) = (y4*y7*y8**2+y2*y7*y9**2)*y1**2+(y4**2*y7*y8**2+y2**2*y7*y9**2)*y1+ &
      y2*y3**2*y7*y8**2+y2**2*y3*y7*y8**2+y3**2*y4*y7*y9**2+y3*y4**2*y7*y9**2
    dF(1,579) = ((-y5*y6+sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y3*y2+(-y5*y6+ &
      sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y4*y3)*y1**2+((-sqrt(3._ark)*y6**2/ &
      2._ark-sqrt(3._ark)*y5**2/6._ark+y5*y6)*y4*y2**2+((-y5*y6+sqrt(3._ark)*y6**2/2._ark+ &
      sqrt(3._ark)*y5**2/6._ark)*y3**2+(-sqrt(3._ark)*y6**2/2._ark-sqrt(3._ark)*y5**2/6._ark+ &
      y5*y6)*y4**2)*y2+(-y5*y6+sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/ &
      6._ark)*y4*y3**2)*y1+(-sqrt(3._ark)*y6**2/2._ark-sqrt(3._ark)*y5**2/6._ark+ &
      y5*y6)*y4*y3*y2**2+(-sqrt(3._ark)*y6**2/2._ark-sqrt(3._ark)*y5**2/6._ark+ &
      y5*y6)*y4**2*y3*y2
    dF(1,580) = (y5*y6-sqrt(3._ark)*y6**2/3._ark)*y4*y2*y1**2+((sqrt(3._ark)*y6**2/3._ark- &
      y5*y6)*y3*y2**2+(sqrt(3._ark)*y6**2/3._ark-y5*y6)*y4**2*y3)*y1+(y5*y6- &
      sqrt(3._ark)*y6**2/3._ark)*y4*y3**2*y2
    dF(1,581) = (y2*y3**2*y7+y3**2*y4*y7)*y1**2+y1*y2**2*y4**2*y7+y2**2*y3*y4**2*y7
    dF(1,582) = ((y4*y9-y3*y9)*y2**2+y2*y4**2*y8-y3*y4**2*y8)*y1**2+(y2**2*y3**2*y8+ &
      y3**2*y4**2*y9)*y1-y2**2*y3**2*y4*y8-y2*y3**2*y4**2*y9
    dF(1,583) = (((sqrt(3._ark)*y6+y5)*y3+(-y5-sqrt(3._ark)*y6)*y4)*y2**2+ &
      2._ark*y2*y4**2*y5-2._ark*y3*y4**2*y5)*y1**2+(-2._ark*y2**2*y3**2*y5+(sqrt(3._ark)*y6+ &
      y5)*y4**2*y3**2)*y1+2._ark*y2**2*y3**2*y4*y5+(-y5-sqrt(3._ark)*y6)*y4**2*y3**2*y2
    dF(1,584) = (y3*y4**2+y2**2*y3)*y1**3+(-y2**3*y4-y2*y4**3)*y1**2+(y2**2*y3**3+ &
      y3**3*y4**2)*y1-y2**3*y3**2*y4-y2*y3**2*y4**3
    dF(1,585) = y1**2*y2*y4*y8*y9+(y3*y4**2*y8*y9+y2**2*y3*y8*y9)*y1+ &
      y2*y3**2*y4*y8*y9
    dF(1,586) = (sqrt(3._ark)*y5*y9+(-y9-2._ark*y8)*y6)*y4*y2*y1**2+((-sqrt(3._ark)*y5*y9+ &
      (-2._ark*y8+y9)*y6)*y3*y2**2+(sqrt(3._ark)*y5*y9+(2._ark*y8-y9)*y6)*y4**2*y3)*y1+(- &
      sqrt(3._ark)*y5*y9+(y9+2._ark*y8)*y6)*y4*y3**2*y2
    dF(1,587) = (-y8-y9)*y4*y3*y2*y1**2+((-y8+y9)*y4*y3*y2**2+((y8+y9)*y4*y3**2+(y8- &
      y9)*y4**2*y3)*y2)*y1
    dF(1,588) = y1**2*y2*y3*y4*y7+(y2**2*y3*y4*y7+(y3**2*y4*y7+y3*y4**2*y7)*y2)*y1
    dF(1,589) = (sqrt(3._ark)*y6-y5)*y4*y3*y2*y1**2+((y5-sqrt(3._ark)*y6)*y4*y3*y2**2+ &
      ((sqrt(3._ark)*y6-y5)*y4*y3**2+(y5-sqrt(3._ark)*y6)*y4**2*y3)*y2)*y1
    dF(1,590) = y1**2*y7**4-y2**2*y7**4+y3**2*y7**4-y4**2*y7**4
    dF(1,591) = (y8+y9)*y7**3*y1**2+(-y8+y9)*y7**3*y2**2+(-y8-y9)*y7**3*y3**2+(y8- &
      y9)*y7**3*y4**2
    dF(1,592) = (y9**4+y8**4)*y1**2+(-y8**4-y9**4)*y2**2+(y9**4+y8**4)*y3**2+(- &
      y8**4-y9**4)*y4**2
    dF(1,593) = ((sqrt(3._ark)*y9**3-sqrt(3._ark)*y8**3)*y5+(y9**3-y8**3)*y6)*y1**2+((- &
      sqrt(3._ark)*y8**3-sqrt(3._ark)*y9**3)*y5+(-y9**3-y8**3)*y6)*y2**2+ &
      ((sqrt(3._ark)*y8**3-sqrt(3._ark)*y9**3)*y5+(-y9**3+y8**3)*y6)*y3**2+ &
      ((sqrt(3._ark)*y8**3+sqrt(3._ark)*y9**3)*y5+(y8**3+y9**3)*y6)*y4**2
    dF(1,594) = (-3._ark/4._ark*sqrt(3._ark)*y5**3*y8+(9._ark/2._ark*y9+9._ark/ &
      4._ark*y8)*y6*y5**2+(-3._ark*sqrt(3._ark)*y9-15._ark/4._ark*sqrt(3._ark)*y8)*y6**2*y5+ &
      (13._ark/4._ark*y8+7._ark/2._ark*y9)*y6**3)*y1**2+(-3._ark/4._ark*sqrt(3._ark)*y5**3*y8+(- &
      9._ark/2._ark*y9+9._ark/4._ark*y8)*y6*y5**2+(3._ark*sqrt(3._ark)*y9-15._ark/ &
      4._ark*sqrt(3._ark)*y8)*y6**2*y5+(13._ark/4._ark*y8-7._ark/2._ark*y9)*y6**3)*y2**2+(3._ark/ &
      4._ark*sqrt(3._ark)*y5**3*y8+(-9._ark/4._ark*y8-9._ark/2._ark*y9)*y6*y5**2+(15._ark/ &
      4._ark*sqrt(3._ark)*y8+3._ark*sqrt(3._ark)*y9)*y6**2*y5+(-7._ark/2._ark*y9-13._ark/ &
      4._ark*y8)*y6**3)*y3**2+(3._ark/4._ark*sqrt(3._ark)*y5**3*y8+(-9._ark/4._ark*y8+9._ark/ &
      2._ark*y9)*y6*y5**2+(-3._ark*sqrt(3._ark)*y9+15._ark/4._ark*sqrt(3._ark)*y8)*y6**2*y5+(7._ark/ &
      2._ark*y9-13._ark/4._ark*y8)*y6**3)*y4**2
    dF(1,595) = (-y5*y7*y8*y9/2._ark+sqrt(3._ark)*y6*y7*y8*y9/2._ark)*y1**2+(y5*y7*y8*y9/ &
      2._ark-sqrt(3._ark)*y6*y7*y8*y9/2._ark)*y2**2+(-y5*y7*y8*y9/2._ark+ &
      sqrt(3._ark)*y6*y7*y8*y9/2._ark)*y3**2+(y5*y7*y8*y9/2._ark-sqrt(3._ark)*y6*y7*y8*y9/ &
      2._ark)*y4**2
    dF(1,596) = (y5*y7**3-sqrt(3._ark)*y6*y7**3)*y1**2+(y5*y7**3- &
      sqrt(3._ark)*y6*y7**3)*y2**2+(y5*y7**3-sqrt(3._ark)*y6*y7**3)*y3**2+(y5*y7**3- &
      sqrt(3._ark)*y6*y7**3)*y4**2
    dF(1,597) = ((y8**3-2._ark*y9**3)*y5+sqrt(3._ark)*y6*y8**3)*y1**2+((y8**3+ &
      2._ark*y9**3)*y5+sqrt(3._ark)*y6*y8**3)*y2**2+((2._ark*y9**3-y8**3)*y5- &
      sqrt(3._ark)*y6*y8**3)*y3**2+((-y8**3-2._ark*y9**3)*y5-sqrt(3._ark)*y6*y8**3)*y4**2
    dF(1,598) = (-sqrt(3._ark)*y6**2*y7**2/6._ark-y5*y6*y7**2-sqrt(3._ark)*y5**2*y7**2/ &
      2._ark)*y1**2+(sqrt(3._ark)*y5**2*y7**2/2._ark+sqrt(3._ark)*y6**2*y7**2/6._ark+ &
      y5*y6*y7**2)*y2**2+(-sqrt(3._ark)*y6**2*y7**2/6._ark-y5*y6*y7**2- &
      sqrt(3._ark)*y5**2*y7**2/2._ark)*y3**2+(sqrt(3._ark)*y5**2*y7**2/2._ark+ &
      sqrt(3._ark)*y6**2*y7**2/6._ark+y5*y6*y7**2)*y4**2
    dF(1,599) = ((sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/2._ark)*y7*y5**2+(y8+y9)*y7*y6*y5+ &
      (sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y8/6._ark)*y7*y6**2)*y1**2+((-sqrt(3._ark)*y8/2._ark+ &
      sqrt(3._ark)*y9/2._ark)*y7*y5**2+(-y8+y9)*y7*y6*y5+(sqrt(3._ark)*y9/6._ark- &
      sqrt(3._ark)*y8/6._ark)*y7*y6**2)*y2**2+((-sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y9/ &
      2._ark)*y7*y5**2+(-y8-y9)*y7*y6*y5+(-sqrt(3._ark)*y8/6._ark-sqrt(3._ark)*y9/ &
      6._ark)*y7*y6**2)*y3**2+((-sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8/2._ark)*y7*y5**2+(y8- &
      y9)*y7*y6*y5+(-sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y8/6._ark)*y7*y6**2)*y4**2
    dF(1,600) = (sqrt(3._ark)*y6**2*y8*y9/2._ark+sqrt(3._ark)*y5**2*y8*y9/6._ark- &
      y5*y6*y8*y9)*y1**2+(sqrt(3._ark)*y6**2*y8*y9/2._ark+sqrt(3._ark)*y5**2*y8*y9/6._ark- &
      y5*y6*y8*y9)*y2**2+(sqrt(3._ark)*y6**2*y8*y9/2._ark+sqrt(3._ark)*y5**2*y8*y9/6._ark- &
      y5*y6*y8*y9)*y3**2+(sqrt(3._ark)*y6**2*y8*y9/2._ark+sqrt(3._ark)*y5**2*y8*y9/6._ark- &
      y5*y6*y8*y9)*y4**2
    dF(1,601) = ((-sqrt(3._ark)*y9-sqrt(3._ark)*y8)*y6*y5**2+(4._ark*y9+4._ark*y8)*y6**2*y5+ &
      (-sqrt(3._ark)*y9-sqrt(3._ark)*y8)*y6**3)*y1**2+((sqrt(3._ark)*y9- &
      sqrt(3._ark)*y8)*y6*y5**2+(4._ark*y8-4._ark*y9)*y6**2*y5+(sqrt(3._ark)*y9- &
      sqrt(3._ark)*y8)*y6**3)*y2**2+((sqrt(3._ark)*y9+sqrt(3._ark)*y8)*y6*y5**2+(-4._ark*y8- &
      4._ark*y9)*y6**2*y5+(sqrt(3._ark)*y9+sqrt(3._ark)*y8)*y6**3)*y3**2+((-sqrt(3._ark)*y9+ &
      sqrt(3._ark)*y8)*y6*y5**2+(4._ark*y9-4._ark*y8)*y6**2*y5+(-sqrt(3._ark)*y9+ &
      sqrt(3._ark)*y8)*y6**3)*y4**2
    dF(1,602) = (-7._ark/24._ark*sqrt(3._ark)*y6**4-3._ark/4._ark*sqrt(3._ark)*y5**2*y6**2+ &
      y5*y6**3-sqrt(3._ark)*y5**4/8._ark)*y1**2+(-y5*y6**3+7._ark/24._ark*sqrt(3._ark)*y6**4+ &
      sqrt(3._ark)*y5**4/8._ark+3._ark/4._ark*sqrt(3._ark)*y5**2*y6**2)*y2**2+(-7._ark/ &
      24._ark*sqrt(3._ark)*y6**4-3._ark/4._ark*sqrt(3._ark)*y5**2*y6**2+y5*y6**3- &
      sqrt(3._ark)*y5**4/8._ark)*y3**2+(-y5*y6**3+7._ark/24._ark*sqrt(3._ark)*y6**4+ &
      sqrt(3._ark)*y5**4/8._ark+3._ark/4._ark*sqrt(3._ark)*y5**2*y6**2)*y4**2
    dF(1,603) = (y5**2*y7**2+y6**2*y7**2)*y1**2+(-y5**2*y7**2-y6**2*y7**2)*y2**2+ &
      (y5**2*y7**2+y6**2*y7**2)*y3**2+(-y5**2*y7**2-y6**2*y7**2)*y4**2
    dF(1,604) = ((y8+y9)*y7*y5**2+(y8+y9)*y7*y6**2)*y1**2+((-y8+y9)*y7*y5**2+(-y8+ &
      y9)*y7*y6**2)*y2**2+((-y8-y9)*y7*y5**2+(-y8-y9)*y7*y6**2)*y3**2+((y8- &
      y9)*y7*y5**2+(y8-y9)*y7*y6**2)*y4**2
    dF(1,605) = (-sqrt(3._ark)*y5**3*y8/4._ark+(7._ark/4._ark*y8+y9/2._ark)*y6*y5**2+(- &
      sqrt(3._ark)*y9-5._ark/4._ark*sqrt(3._ark)*y8)*y6**2*y5+(3._ark/4._ark*y8+3._ark/ &
      2._ark*y9)*y6**3)*y1**2+(-sqrt(3._ark)*y5**3*y8/4._ark+(-y9/2._ark+7._ark/ &
      4._ark*y8)*y6*y5**2+(sqrt(3._ark)*y9-5._ark/4._ark*sqrt(3._ark)*y8)*y6**2*y5+(-3._ark/ &
      2._ark*y9+3._ark/4._ark*y8)*y6**3)*y2**2+(sqrt(3._ark)*y5**3*y8/4._ark+(-y9/2._ark-7._ark/ &
      4._ark*y8)*y6*y5**2+(sqrt(3._ark)*y9+5._ark/4._ark*sqrt(3._ark)*y8)*y6**2*y5+(-3._ark/ &
      2._ark*y9-3._ark/4._ark*y8)*y6**3)*y3**2+(sqrt(3._ark)*y5**3*y8/4._ark+(y9/2._ark-7._ark/ &
      4._ark*y8)*y6*y5**2+(-sqrt(3._ark)*y9+5._ark/4._ark*sqrt(3._ark)*y8)*y6**2*y5+(3._ark/ &
      2._ark*y9-3._ark/4._ark*y8)*y6**3)*y4**2
    dF(1,606) = (-7._ark/24._ark*sqrt(3._ark)*y5**4-sqrt(3._ark)*y6**4/8._ark-y5**3*y6-3._ark/ &
      4._ark*sqrt(3._ark)*y5**2*y6**2)*y1**2+(7._ark/24._ark*sqrt(3._ark)*y5**4+ &
      sqrt(3._ark)*y6**4/8._ark+y5**3*y6+3._ark/4._ark*sqrt(3._ark)*y5**2*y6**2)*y2**2+(-7._ark/ &
      24._ark*sqrt(3._ark)*y5**4-sqrt(3._ark)*y6**4/8._ark-y5**3*y6-3._ark/ &
      4._ark*sqrt(3._ark)*y5**2*y6**2)*y3**2+(7._ark/24._ark*sqrt(3._ark)*y5**4+ &
      sqrt(3._ark)*y6**4/8._ark+y5**3*y6+3._ark/4._ark*sqrt(3._ark)*y5**2*y6**2)*y4**2
    dF(1,607) = y2*y4**2*y7**3+y1**2*y3*y7**3+y2**2*y4*y7**3+y1*y3**2*y7**3
    dF(1,608) = (y2*y9**3+y4*y8**3)*y1**2+(-y4**2*y8**3-y2**2*y9**3)*y1+ &
      y2**2*y3*y8**3+y3*y4**2*y9**3-y2*y3**2*y8**3-y3**2*y4*y9**3
    dF(1,609) = y1*y3**2*y7*y8*y9+y1**2*y3*y7*y8*y9-y2*y4**2*y7*y8*y9- &
      y2**2*y4*y7*y8*y9
    dF(1,610) = (y2*y7*y8*y9+y4*y7*y8*y9)*y1**2+(-y2**2*y7*y8*y9-y4**2*y7*y8*y9)*y1- &
      y2**2*y3*y7*y8*y9-y3*y4**2*y7*y8*y9+y3**2*y4*y7*y8*y9+y2*y3**2*y7*y8*y9
    dF(1,611) = (y2*y8*y9**2+y4*y8**2*y9)*y1**2+(y4**2*y8**2*y9+y2**2*y8*y9**2)*y1- &
      y2*y3**2*y8**2*y9-y3**2*y4*y8*y9**2-y2**2*y3*y8**2*y9-y3*y4**2*y8*y9**2
    dF(1,612) = (-y6*y8*y9+sqrt(3._ark)*y5*y8*y9/3._ark)*y3*y1**2+(-y6*y8*y9+ &
      sqrt(3._ark)*y5*y8*y9/3._ark)*y3**2*y1+(-y6*y8*y9+sqrt(3._ark)*y5*y8*y9/ &
      3._ark)*y4*y2**2+(-y6*y8*y9+sqrt(3._ark)*y5*y8*y9/3._ark)*y4**2*y2
    dF(1,613) = (3._ark/4._ark*y5**2*y9-sqrt(3._ark)*y5*y6*y9/2._ark+(y9/4._ark+ &
      y8)*y6**2)*y3*y1**2+(-3._ark/4._ark*y5**2*y9+sqrt(3._ark)*y5*y6*y9/2._ark+(-y9/4._ark- &
      y8)*y6**2)*y3**2*y1+(-3._ark/4._ark*y5**2*y9+sqrt(3._ark)*y5*y6*y9/2._ark+(y8-y9/ &
      4._ark)*y6**2)*y4*y2**2+(3._ark/4._ark*y5**2*y9-sqrt(3._ark)*y5*y6*y9/2._ark+(-y8+y9/ &
      4._ark)*y6**2)*y4**2*y2
    dF(1,614) = ((-y5*y7*y8/2._ark+sqrt(3._ark)*y6*y7*y8/2._ark)*y2+(-y5*y7*y9/2._ark+ &
      sqrt(3._ark)*y6*y7*y9/2._ark)*y4)*y1**2+((y5*y7*y8/2._ark-sqrt(3._ark)*y6*y7*y8/ &
      2._ark)*y2**2+(y5*y7*y9/2._ark-sqrt(3._ark)*y6*y7*y9/2._ark)*y4**2)*y1+(-y5*y7*y9/2._ark+ &
      sqrt(3._ark)*y6*y7*y9/2._ark)*y3*y2**2+(y5*y7*y9/2._ark-sqrt(3._ark)*y6*y7*y9/ &
      2._ark)*y3**2*y2+(y5*y7*y8/2._ark-sqrt(3._ark)*y6*y7*y8/2._ark)*y4*y3**2+(-y5*y7*y8/ &
      2._ark+sqrt(3._ark)*y6*y7*y8/2._ark)*y4**2*y3
    dF(1,615) = (y5*y6**2-y5**3/3._ark)*y3*y1**2+(y5*y6**2-y5**3/3._ark)*y3**2*y1+ &
      (y5**3/3._ark-y5*y6**2)*y4*y2**2+(y5**3/3._ark-y5*y6**2)*y4**2*y2
    dF(1,616) = ((-2._ark*y9**2+y8**2)*y5+sqrt(3._ark)*y6*y8**2)*y3*y1**2+((-2._ark*y9**2+ &
      y8**2)*y5+sqrt(3._ark)*y6*y8**2)*y3**2*y1+((-y8**2+2._ark*y9**2)*y5- &
      sqrt(3._ark)*y6*y8**2)*y4*y2**2+((-y8**2+2._ark*y9**2)*y5- &
      sqrt(3._ark)*y6*y8**2)*y4**2*y2
    dF(1,617) = ((y5*y7*y9+sqrt(3._ark)*y6*y7*y9)*y2-2._ark*y4*y5*y7*y8)*y1**2+ &
      ((y5*y7*y9+sqrt(3._ark)*y6*y7*y9)*y2**2-2._ark*y4**2*y5*y7*y8)*y1+ &
      2._ark*y2**2*y3*y5*y7*y8+2._ark*y2*y3**2*y5*y7*y8+(-sqrt(3._ark)*y6*y7*y9- &
      y5*y7*y9)*y4*y3**2+(-sqrt(3._ark)*y6*y7*y9-y5*y7*y9)*y4**2*y3
    dF(1,618) = ((-y5*y6*y8-sqrt(3._ark)*y5**2*y8/2._ark-sqrt(3._ark)*y6**2*y8/6._ark)*y2+(- &
      y5*y6*y9-sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y4)*y1**2+((- &
      y5*y6*y8-sqrt(3._ark)*y5**2*y8/2._ark-sqrt(3._ark)*y6**2*y8/6._ark)*y2**2+(-y5*y6*y9- &
      sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y4**2)*y1+ &
      (sqrt(3._ark)*y6**2*y9/6._ark+y5*y6*y9+sqrt(3._ark)*y5**2*y9/2._ark)*y3*y2**2+ &
      (sqrt(3._ark)*y6**2*y9/6._ark+y5*y6*y9+sqrt(3._ark)*y5**2*y9/2._ark)*y3**2*y2+ &
      (sqrt(3._ark)*y6**2*y8/6._ark+sqrt(3._ark)*y5**2*y8/2._ark+y5*y6*y8)*y4*y3**2+ &
      (sqrt(3._ark)*y6**2*y8/6._ark+sqrt(3._ark)*y5**2*y8/2._ark+y5*y6*y8)*y4**2*y3
    dF(1,619) = ((2._ark*y5*y6*y7+sqrt(3._ark)*y5**2*y7/2._ark+sqrt(3._ark)*y6**2*y7/ &
      2._ark)*y2+(y5*y6*y7+sqrt(3._ark)*y5**2*y7)*y4)*y1**2+((2._ark*y5*y6*y7+ &
      sqrt(3._ark)*y5**2*y7/2._ark+sqrt(3._ark)*y6**2*y7/2._ark)*y2**2+(y5*y6*y7+ &
      sqrt(3._ark)*y5**2*y7)*y4**2)*y1+(y5*y6*y7+sqrt(3._ark)*y5**2*y7)*y3*y2**2+ &
      (y5*y6*y7+sqrt(3._ark)*y5**2*y7)*y3**2*y2+(2._ark*y5*y6*y7+sqrt(3._ark)*y5**2*y7/2._ark+ &
      sqrt(3._ark)*y6**2*y7/2._ark)*y4*y3**2+(2._ark*y5*y6*y7+sqrt(3._ark)*y5**2*y7/2._ark+ &
      sqrt(3._ark)*y6**2*y7/2._ark)*y4**2*y3
    dF(1,620) = (sqrt(3._ark)*y5**2*y9/4._ark+(y8+y9/2._ark)*y6*y5-sqrt(3._ark)*y6**2*y9/ &
      4._ark)*y3*y1**2+(-sqrt(3._ark)*y5**2*y9/4._ark+(-y9/2._ark-y8)*y6*y5+ &
      sqrt(3._ark)*y6**2*y9/4._ark)*y3**2*y1+(-sqrt(3._ark)*y5**2*y9/4._ark+(y8-y9/ &
      2._ark)*y6*y5+sqrt(3._ark)*y6**2*y9/4._ark)*y4*y2**2+(sqrt(3._ark)*y5**2*y9/4._ark+(y9/ &
      2._ark-y8)*y6*y5-sqrt(3._ark)*y6**2*y9/4._ark)*y4**2*y2
    dF(1,621) = ((y5*y6*y8+sqrt(3._ark)*y6**2*y8/3._ark)*y2+(-sqrt(3._ark)*y6**2*y9/6._ark+ &
      sqrt(3._ark)*y5**2*y9/2._ark)*y4)*y1**2+((y5*y6*y8+sqrt(3._ark)*y6**2*y8/3._ark)*y2**2+ &
      (-sqrt(3._ark)*y6**2*y9/6._ark+sqrt(3._ark)*y5**2*y9/2._ark)*y4**2)*y1+ &
      (sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y3*y2**2+ &
      (sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y3**2*y2+(- &
      sqrt(3._ark)*y6**2*y8/3._ark-y5*y6*y8)*y4*y3**2+(-sqrt(3._ark)*y6**2*y8/3._ark- &
      y5*y6*y8)*y4**2*y3
    dF(1,622) = ((-2._ark/3._ark*y5**3+4._ark/9._ark*sqrt(3._ark)*y5**2*y6+2._ark/ &
      3._ark*y5*y6**2)*y2+(y5*y6**2-y5**3/3._ark+sqrt(3._ark)*y6**3/3._ark- &
      sqrt(3._ark)*y5**2*y6/9._ark)*y4)*y1**2+((2._ark/3._ark*y5**3-2._ark/3._ark*y5*y6**2-4._ark/ &
      9._ark*sqrt(3._ark)*y5**2*y6)*y2**2+(-y5*y6**2-sqrt(3._ark)*y6**3/3._ark+ &
      sqrt(3._ark)*y5**2*y6/9._ark+y5**3/3._ark)*y4**2)*y1+(-y5*y6**2-sqrt(3._ark)*y6**3/3._ark+ &
      sqrt(3._ark)*y5**2*y6/9._ark+y5**3/3._ark)*y3*y2**2+(y5*y6**2-y5**3/3._ark+ &
      sqrt(3._ark)*y6**3/3._ark-sqrt(3._ark)*y5**2*y6/9._ark)*y3**2*y2+(-2._ark/3._ark*y5**3+4._ark/ &
      9._ark*sqrt(3._ark)*y5**2*y6+2._ark/3._ark*y5*y6**2)*y4*y3**2+(2._ark/3._ark*y5**3-2._ark/ &
      3._ark*y5*y6**2-4._ark/9._ark*sqrt(3._ark)*y5**2*y6)*y4**2*y3
    dF(1,623) = ((-2._ark*y5**2*y7-2._ark*sqrt(3._ark)*y5*y6*y7)*y2+(-2._ark*y5**2*y7- &
      2._ark*sqrt(3._ark)*y5*y6*y7)*y4)*y1**2+((-2._ark*y5**2*y7- &
      2._ark*sqrt(3._ark)*y5*y6*y7)*y2**2+(-2._ark*y5**2*y7- &
      2._ark*sqrt(3._ark)*y5*y6*y7)*y4**2)*y1+(-2._ark*y5**2*y7- &
      2._ark*sqrt(3._ark)*y5*y6*y7)*y3*y2**2+(-2._ark*y5**2*y7- &
      2._ark*sqrt(3._ark)*y5*y6*y7)*y3**2*y2+(-2._ark*y5**2*y7- &
      2._ark*sqrt(3._ark)*y5*y6*y7)*y4*y3**2+(-2._ark*y5**2*y7- &
      2._ark*sqrt(3._ark)*y5*y6*y7)*y4**2*y3
    dF(1,624) = ((y9/4._ark+y8)*y5**2+sqrt(3._ark)*y5*y6*y9/2._ark+3._ark/ &
      4._ark*y6**2*y9)*y3*y1**2+((-y9/4._ark-y8)*y5**2-sqrt(3._ark)*y5*y6*y9/2._ark-3._ark/ &
      4._ark*y6**2*y9)*y3**2*y1+((y8-y9/4._ark)*y5**2-sqrt(3._ark)*y5*y6*y9/2._ark-3._ark/ &
      4._ark*y6**2*y9)*y4*y2**2+((-y8+y9/4._ark)*y5**2+sqrt(3._ark)*y5*y6*y9/2._ark+3._ark/ &
      4._ark*y6**2*y9)*y4**2*y2
    dF(1,625) = (y4**2*y8*y9+y2**2*y8*y9)*y1**2+y3**2*y4**2*y8*y9+y2**2*y3**2*y8*y9
    dF(1,626) = (y5*y6+sqrt(3._ark)*y5**2/3._ark)*y3**2*y1**2+(-y5*y6-sqrt(3._ark)*y5**2/ &
      3._ark)*y4**2*y2**2
    dF(1,627) = (y8**2+y9**2)*y7*y3*y1**2+(y8**2+y9**2)*y7*y3**2*y1+(y8**2+ &
      y9**2)*y7*y4*y2**2+(y8**2+y9**2)*y7*y4**2*y2
    dF(1,628) = ((-sqrt(3._ark)*y5*y6*y9+y5**2*y9)*y2+(-y5**2*y8/2._ark+3._ark/ &
      2._ark*y6**2*y8)*y4)*y1**2+((sqrt(3._ark)*y5*y6*y9-y5**2*y9)*y2**2+(y5**2*y8/2._ark- &
      3._ark/2._ark*y6**2*y8)*y4**2)*y1+(-y5**2*y8/2._ark+3._ark/2._ark*y6**2*y8)*y3*y2**2+ &
      (y5**2*y8/2._ark-3._ark/2._ark*y6**2*y8)*y3**2*y2+(sqrt(3._ark)*y5*y6*y9- &
      y5**2*y9)*y4*y3**2+(-sqrt(3._ark)*y5*y6*y9+y5**2*y9)*y4**2*y3
    dF(1,629) = (y6**2*y7+y5**2*y7)*y3*y1**2+(y6**2*y7+y5**2*y7)*y3**2*y1+(y6**2*y7+ &
      y5**2*y7)*y4*y2**2+(y6**2*y7+y5**2*y7)*y4**2*y2
    dF(1,630) = (y2*y3*y7*y9+y3*y4*y7*y8)*y1**2+(y2**2*y4*y7*y9+(y4**2*y7*y8- &
      y3**2*y7*y8)*y2-y3**2*y4*y7*y9)*y1-y2*y3*y4**2*y7*y9-y2**2*y3*y4*y7*y8
    dF(1,631) = ((-y6*y8/2._ark-sqrt(3._ark)*y5*y8/2._ark)*y3*y2+(y6*y9/2._ark+ &
      sqrt(3._ark)*y5*y9/2._ark)*y4*y3)*y1**2+((-y6*y8/2._ark-sqrt(3._ark)*y5*y8/ &
      2._ark)*y4*y2**2+((-sqrt(3._ark)*y5*y9/2._ark-y6*y9/2._ark)*y3**2+(y6*y9/2._ark+ &
      sqrt(3._ark)*y5*y9/2._ark)*y4**2)*y2+(sqrt(3._ark)*y5*y8/2._ark+y6*y8/ &
      2._ark)*y4*y3**2)*y1+(-sqrt(3._ark)*y5*y9/2._ark-y6*y9/2._ark)*y4*y3*y2**2+ &
      (sqrt(3._ark)*y5*y8/2._ark+y6*y8/2._ark)*y4**2*y3*y2
    dF(1,632) = ((y8-2._ark*y9)*y5+sqrt(3._ark)*y6*y8)*y4*y2*y1**2+(((y8+2._ark*y9)*y5+ &
      sqrt(3._ark)*y6*y8)*y3*y2**2+((-y8-2._ark*y9)*y5-sqrt(3._ark)*y6*y8)*y4**2*y3)*y1+((- &
      y8+2._ark*y9)*y5-sqrt(3._ark)*y6*y8)*y4*y3**2*y2
    dF(1,633) = ((-sqrt(3._ark)*y6*y7+y5*y7)*y3*y2+(-sqrt(3._ark)*y6*y7+ &
      y5*y7)*y4*y3)*y1**2+((-sqrt(3._ark)*y6*y7+y5*y7)*y4*y2**2+((-sqrt(3._ark)*y6*y7+ &
      y5*y7)*y3**2+(-sqrt(3._ark)*y6*y7+y5*y7)*y4**2)*y2+(-sqrt(3._ark)*y6*y7+ &
      y5*y7)*y4*y3**2)*y1+(-sqrt(3._ark)*y6*y7+y5*y7)*y4*y3*y2**2+(-sqrt(3._ark)*y6*y7+ &
      y5*y7)*y4**2*y3*y2
    dF(1,634) = (-sqrt(3._ark)*y6*y7+y5*y7)*y4*y2*y1**2+((-sqrt(3._ark)*y6*y7+ &
      y5*y7)*y3*y2**2+(-sqrt(3._ark)*y6*y7+y5*y7)*y4**2*y3)*y1+(-sqrt(3._ark)*y6*y7+ &
      y5*y7)*y4*y3**2*y2
    dF(1,635) = ((sqrt(3._ark)*y6*y8/2._ark-y5*y8/2._ark)*y3*y2+(sqrt(3._ark)*y6*y9/2._ark- &
      y5*y9/2._ark)*y4*y3)*y1**2+((sqrt(3._ark)*y6*y8/2._ark-y5*y8/2._ark)*y4*y2**2+((- &
      sqrt(3._ark)*y6*y9/2._ark+y5*y9/2._ark)*y3**2+(sqrt(3._ark)*y6*y9/2._ark-y5*y9/ &
      2._ark)*y4**2)*y2+(-sqrt(3._ark)*y6*y8/2._ark+y5*y8/2._ark)*y4*y3**2)*y1+(- &
      sqrt(3._ark)*y6*y9/2._ark+y5*y9/2._ark)*y4*y3*y2**2+(-sqrt(3._ark)*y6*y8/2._ark+y5*y8/ &
      2._ark)*y4**2*y3*y2
    dF(1,636) = ((-sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y3*y2+(y5*y6- &
      sqrt(3._ark)*y5**2/3._ark)*y4*y3)*y1**2+((-sqrt(3._ark)*y5**2/6._ark+sqrt(3._ark)*y6**2/ &
      2._ark)*y4*y2**2+((y5*y6-sqrt(3._ark)*y5**2/3._ark)*y3**2+(sqrt(3._ark)*y5**2/3._ark- &
      y5*y6)*y4**2)*y2+(-sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y4*y3**2)*y1+ &
      (sqrt(3._ark)*y5**2/3._ark-y5*y6)*y4*y3*y2**2+(-sqrt(3._ark)*y5**2/6._ark+ &
      sqrt(3._ark)*y6**2/2._ark)*y4**2*y3*y2
    dF(1,637) = (y6**2+y5**2)*y4*y2*y1**2+((-y6**2-y5**2)*y3*y2**2+(-y6**2- &
      y5**2)*y4**2*y3)*y1+(y6**2+y5**2)*y4*y3**2*y2
    dF(1,638) = (((y6+sqrt(3._ark)*y5)*y3+(-y6-sqrt(3._ark)*y5)*y4)*y2**2+(y6+ &
      sqrt(3._ark)*y5)*y4**2*y2+(-y6-sqrt(3._ark)*y5)*y4**2*y3)*y1**2+((-y6- &
      sqrt(3._ark)*y5)*y3**2*y2**2+(y6+sqrt(3._ark)*y5)*y4**2*y3**2)*y1+(y6+ &
      sqrt(3._ark)*y5)*y4*y3**2*y2**2+(-y6-sqrt(3._ark)*y5)*y4**2*y3**2*y2
    dF(1,639) = (y3**2*y4+y2*y3**2)*y1**3+(y3**3*y4+y2*y3**3)*y1**2+(-y2**2*y4**3- &
      y2**3*y4**2)*y1-y2**2*y3*y4**3-y2**3*y3*y4**2
    dF(1,640) = (-y2*y4**2-y2**2*y4)*y1**3+(y3*y4**3+y2**3*y3)*y1**2+(y3**2*y4**3+ &
      y2**3*y3**2)*y1-y2*y3**3*y4**2-y2**2*y3**3*y4
    dF(1,641) = (2._ark*y2*y3*y6*y7+(y6*y7-sqrt(3._ark)*y5*y7)*y4*y3)*y1**2+ &
      (2._ark*y2**2*y4*y6*y7+((y6*y7-sqrt(3._ark)*y5*y7)*y3**2+(y6*y7- &
      sqrt(3._ark)*y5*y7)*y4**2)*y2+2._ark*y3**2*y4*y6*y7)*y1+(y6*y7- &
      sqrt(3._ark)*y5*y7)*y4*y3*y2**2+2._ark*y2*y3*y4**2*y6*y7
    dF(1,642) = ((y6**2+y5**2)*y3*y2+(y6**2+y5**2)*y4*y3)*y1**2+((-y6**2- &
      y5**2)*y4*y2**2+((y6**2+y5**2)*y3**2+(-y6**2-y5**2)*y4**2)*y2+(y6**2+ &
      y5**2)*y4*y3**2)*y1+(-y6**2-y5**2)*y4*y3*y2**2+(-y6**2-y5**2)*y4**2*y3*y2
    dF(1,643) = ((sqrt(3._ark)*y6*y9+y5*y9)*y3*y2-2._ark*y3*y4*y5*y8)*y1**2+((- &
      sqrt(3._ark)*y6*y9-y5*y9)*y4*y2**2+(2._ark*y4**2*y5*y8+2._ark*y3**2*y5*y8)*y2+(- &
      sqrt(3._ark)*y6*y9-y5*y9)*y4*y3**2)*y1-2._ark*y2**2*y3*y4*y5*y8+(sqrt(3._ark)*y6*y9+ &
      y5*y9)*y4**2*y3*y2
    dF(1,644) = y1**2*y2*y3**2*y4-y1*y2**2*y3*y4**2
    dF(1,645) = y1**2*y3**2*y7**2-y2**2*y4**2*y7**2
    dF(1,646) = ((y3**2-y4**2)*y2**2+y3**2*y4**2)*y1**2-y2**2*y3**2*y4**2
    dF(1,647) = (y6*y7**2-sqrt(3._ark)*y5*y7**2/3._ark)*y1**3+(-y6*y7**2+ &
      sqrt(3._ark)*y5*y7**2/3._ark)*y2**3+(y6*y7**2-sqrt(3._ark)*y5*y7**2/3._ark)*y3**3+(- &
      y6*y7**2+sqrt(3._ark)*y5*y7**2/3._ark)*y4**3
    dF(1,648) = (-y5*y8*y9/2._ark+sqrt(3._ark)*y6*y8*y9/2._ark)*y1**3+(-y5*y8*y9/2._ark+ &
      sqrt(3._ark)*y6*y8*y9/2._ark)*y2**3+(-y5*y8*y9/2._ark+sqrt(3._ark)*y6*y8*y9/2._ark)*y3**3+ &
      (-y5*y8*y9/2._ark+sqrt(3._ark)*y6*y8*y9/2._ark)*y4**3
    dF(1,649) = (-sqrt(3._ark)*y6**2*y7/3._ark+y5*y6*y7)*y1**3+(-sqrt(3._ark)*y6**2*y7/ &
      3._ark+y5*y6*y7)*y2**3+(-sqrt(3._ark)*y6**2*y7/3._ark+y5*y6*y7)*y3**3+(- &
      sqrt(3._ark)*y6**2*y7/3._ark+y5*y6*y7)*y4**3
    dF(1,650) = ((sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y8)*y5**2+(y8+2._ark*y9)*y6*y5+ &
      sqrt(3._ark)*y6**2*y9/2._ark)*y1**3+((sqrt(3._ark)*y8-sqrt(3._ark)*y9/2._ark)*y5**2+(y8- &
      2._ark*y9)*y6*y5-sqrt(3._ark)*y6**2*y9/2._ark)*y2**3+((-sqrt(3._ark)*y8-sqrt(3._ark)*y9/ &
      2._ark)*y5**2+(-y8-2._ark*y9)*y6*y5-sqrt(3._ark)*y6**2*y9/2._ark)*y3**3+((- &
      sqrt(3._ark)*y8+sqrt(3._ark)*y9/2._ark)*y5**2+(-y8+2._ark*y9)*y6*y5+sqrt(3._ark)*y6**2*y9/ &
      2._ark)*y4**3
    dF(1,651) = (-y5**2*y6+sqrt(3._ark)*y5**3-5._ark/3._ark*sqrt(3._ark)*y5*y6**2- &
      y6**3)*y1**3+(-sqrt(3._ark)*y5**3+y6**3+5._ark/3._ark*sqrt(3._ark)*y5*y6**2+ &
      y5**2*y6)*y2**3+(-y5**2*y6+sqrt(3._ark)*y5**3-5._ark/3._ark*sqrt(3._ark)*y5*y6**2- &
      y6**3)*y3**3+(-sqrt(3._ark)*y5**3+y6**3+5._ark/3._ark*sqrt(3._ark)*y5*y6**2+ &
      y5**2*y6)*y4**3
    dF(1,652) = y2**3*y4*y8*y9+y1**3*y3*y8*y9+y1*y3**3*y8*y9+y2*y4**3*y8*y9
    dF(1,653) = (y4*y7*y9+y2*y7*y8)*y1**3+(-y4**3*y7*y9-y2**3*y7*y8)*y1- &
      y2*y3**3*y7*y9-y3**3*y4*y7*y8+y3*y4**3*y7*y8+y2**3*y3*y7*y9
    dF(1,654) = (y8+y9)*y7*y3*y1**3+(-y8-y9)*y7*y3**3*y1+(-y8+y9)*y7*y4*y2**3+(y8- &
      y9)*y7*y4**3*y2
    dF(1,655) = (y4*y9**2+y2*y8**2)*y1**3+(-y4**3*y9**2-y2**3*y8**2)*y1- &
      y2**3*y3*y9**2+y2*y3**3*y9**2-y3*y4**3*y8**2+y3**3*y4*y8**2
    dF(1,656) = (y8**2+y9**2)*y3*y1**3+(y8**2+y9**2)*y3**3*y1+(-y9**2- &
      y8**2)*y4*y2**3+(-y9**2-y8**2)*y4**3*y2
    dF(1,657) = (-sqrt(3._ark)*y5*y9/2._ark+(y8+y9/2._ark)*y6)*y3*y1**3+(sqrt(3._ark)*y5*y9/ &
      2._ark+(-y9/2._ark-y8)*y6)*y3**3*y1+(sqrt(3._ark)*y5*y9/2._ark+(y8-y9/ &
      2._ark)*y6)*y4*y2**3+(-sqrt(3._ark)*y5*y9/2._ark+(y9/2._ark-y8)*y6)*y4**3*y2
    dF(1,658) = ((-y6*y7/2._ark+sqrt(3._ark)*y5*y7/2._ark)*y2-y4*y6*y7)*y1**3+((-y6*y7/ &
      2._ark+sqrt(3._ark)*y5*y7/2._ark)*y2**3-y4**3*y6*y7)*y1-y2**3*y3*y6*y7-y2*y3**3*y6*y7+ &
      (-y6*y7/2._ark+sqrt(3._ark)*y5*y7/2._ark)*y4*y3**3+(-y6*y7/2._ark+sqrt(3._ark)*y5*y7/ &
      2._ark)*y4**3*y3
    dF(1,659) = ((sqrt(3._ark)*y5*y8/2._ark+y6*y8/2._ark)*y2+(-sqrt(3._ark)*y5*y9/2._ark- &
      y6*y9/2._ark)*y4)*y1**3+((sqrt(3._ark)*y5*y8/2._ark+y6*y8/2._ark)*y2**3+(- &
      sqrt(3._ark)*y5*y9/2._ark-y6*y9/2._ark)*y4**3)*y1+(y6*y9/2._ark+sqrt(3._ark)*y5*y9/ &
      2._ark)*y3*y2**3+(y6*y9/2._ark+sqrt(3._ark)*y5*y9/2._ark)*y3**3*y2+(-y6*y8/2._ark- &
      sqrt(3._ark)*y5*y8/2._ark)*y4*y3**3+(-y6*y8/2._ark-sqrt(3._ark)*y5*y8/2._ark)*y4**3*y3
    dF(1,660) = (2._ark*y2*y6*y9+(-sqrt(3._ark)*y5*y8+y6*y8)*y4)*y1**3+(- &
      2._ark*y2**3*y6*y9+(sqrt(3._ark)*y5*y8-y6*y8)*y4**3)*y1+(-sqrt(3._ark)*y5*y8+ &
      y6*y8)*y3*y2**3+(sqrt(3._ark)*y5*y8-y6*y8)*y3**3*y2+2._ark*y3*y4**3*y6*y9- &
      2._ark*y3**3*y4*y6*y9
    dF(1,661) = ((y8-y9/2._ark)*y5-sqrt(3._ark)*y6*y9/2._ark)*y3*y1**3+((y9/2._ark-y8)*y5+ &
      sqrt(3._ark)*y6*y9/2._ark)*y3**3*y1+((y8+y9/2._ark)*y5+sqrt(3._ark)*y6*y9/ &
      2._ark)*y4*y2**3+((-y9/2._ark-y8)*y5-sqrt(3._ark)*y6*y9/2._ark)*y4**3*y2
    dF(1,662) = ((-sqrt(3._ark)*y6*y7/2._ark-y5*y7/2._ark)*y2+y4*y5*y7)*y1**3+((- &
      sqrt(3._ark)*y6*y7/2._ark-y5*y7/2._ark)*y2**3+y4**3*y5*y7)*y1+y2**3*y3*y5*y7+ &
      y2*y3**3*y5*y7+(-sqrt(3._ark)*y6*y7/2._ark-y5*y7/2._ark)*y4*y3**3+(-sqrt(3._ark)*y6*y7/ &
      2._ark-y5*y7/2._ark)*y4**3*y3
    dF(1,663) = ((sqrt(3._ark)*y6*y8/2._ark-y5*y8/2._ark)*y2+(sqrt(3._ark)*y6*y9/2._ark-y5*y9/ &
      2._ark)*y4)*y1**3+((sqrt(3._ark)*y6*y8/2._ark-y5*y8/2._ark)*y2**3+(sqrt(3._ark)*y6*y9/ &
      2._ark-y5*y9/2._ark)*y4**3)*y1+(-sqrt(3._ark)*y6*y9/2._ark+y5*y9/2._ark)*y3*y2**3+(- &
      sqrt(3._ark)*y6*y9/2._ark+y5*y9/2._ark)*y3**3*y2+(-sqrt(3._ark)*y6*y8/2._ark+y5*y8/ &
      2._ark)*y4*y3**3+(-sqrt(3._ark)*y6*y8/2._ark+y5*y8/2._ark)*y4**3*y3
    dF(1,664) = (-sqrt(3._ark)*y6**2/2._ark-sqrt(3._ark)*y5**2/6._ark+y5*y6)*y3*y1**3+(- &
      sqrt(3._ark)*y6**2/2._ark-sqrt(3._ark)*y5**2/6._ark+y5*y6)*y3**3*y1+(-y5*y6+ &
      sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y4*y2**3+(-y5*y6+sqrt(3._ark)*y6**2/ &
      2._ark+sqrt(3._ark)*y5**2/6._ark)*y4**3*y2
    dF(1,665) = (2._ark/3._ark*sqrt(3._ark)*y2*y6**2+(sqrt(3._ark)*y5**2/2._ark-y5*y6+ &
      sqrt(3._ark)*y6**2/6._ark)*y4)*y1**3+(-2._ark/3._ark*sqrt(3._ark)*y2**3*y6**2+(- &
      sqrt(3._ark)*y5**2/2._ark+y5*y6-sqrt(3._ark)*y6**2/6._ark)*y4**3)*y1+(-sqrt(3._ark)*y5**2/ &
      2._ark+y5*y6-sqrt(3._ark)*y6**2/6._ark)*y3*y2**3+(sqrt(3._ark)*y5**2/2._ark-y5*y6+ &
      sqrt(3._ark)*y6**2/6._ark)*y3**3*y2+2._ark/3._ark*sqrt(3._ark)*y3**3*y4*y6**2-2._ark/ &
      3._ark*sqrt(3._ark)*y3*y4**3*y6**2
    dF(1,666) = (y8+y9)*y4*y2*y1**3+((y8-y9)*y3*y2**3+(-y8+y9)*y4**3*y3)*y1+(-y8- &
      y9)*y4*y3**3*y2
    dF(1,667) = -y2*y4**3*y7**2-y2**3*y4*y7**2+y1**3*y3*y7**2+y1*y3**3*y7**2
    dF(1,668) = y1**3*y2*y3*y4+(-y2**3*y3*y4+(-y3*y4**3+y3**3*y4)*y2)*y1
    dF(1,669) = (y8+y9)*y7*y1**4+(-y8+y9)*y7*y2**4+(-y8-y9)*y7*y3**4+(y8- &
      y9)*y7*y4**4
    dF(1,670) = y2**4*y8*y9+y1**4*y8*y9+y4**4*y8*y9+y3**4*y8*y9
    dF(1,671) = -y2**4*y7**2-y4**4*y7**2+y3**4*y7**2+y1**4*y7**2
    dF(1,672) = (y6*y7-sqrt(3._ark)*y5*y7/3._ark)*y1**4+(y6*y7-sqrt(3._ark)*y5*y7/ &
      3._ark)*y2**4+(y6*y7-sqrt(3._ark)*y5*y7/3._ark)*y3**4+(y6*y7-sqrt(3._ark)*y5*y7/ &
      3._ark)*y4**4
    dF(1,673) = (sqrt(3._ark)*y5*y9/2._ark+(-y9/2._ark-y8)*y6)*y1**4+(-sqrt(3._ark)*y5*y9/ &
      2._ark+(y9/2._ark-y8)*y6)*y2**4+(-sqrt(3._ark)*y5*y9/2._ark+(y8+y9/2._ark)*y6)*y3**4+ &
      (sqrt(3._ark)*y5*y9/2._ark+(y8-y9/2._ark)*y6)*y4**4
    dF(1,674) = (y4*y8+y2*y9)*y1**4+(-y4**4*y8-y2**4*y9)*y1-y3**4*y4*y9-y2*y3**4*y8+ &
      y2**4*y3*y8+y3*y4**4*y9
    dF(1,675) = y2*y4**4*y7+y2**4*y4*y7+y1*y3**4*y7+y1**4*y3*y7
    dF(1,676) = (y2*y7+y4*y7)*y1**4+(y4**4*y7+y2**4*y7)*y1+y3*y4**4*y7+y3**4*y4*y7+ &
      y2*y3**4*y7+y2**4*y3*y7
    dF(1,677) = (y4*y9+y2*y8)*y1**4+(y4**4*y9+y2**4*y8)*y1-y3*y4**4*y8-y2*y3**4*y9- &
      y2**4*y3*y9-y3**4*y4*y8
    dF(1,678) = (y8+y9)*y3*y1**4+(-y8-y9)*y3**4*y1+(y8-y9)*y4*y2**4+(-y8+ &
      y9)*y4**4*y2
    dF(1,679) = (sqrt(3._ark)*y6/2._ark-y5/2._ark)*y3*y1**4+(sqrt(3._ark)*y6/2._ark-y5/ &
      2._ark)*y3**4*y1+(-sqrt(3._ark)*y6/2._ark+y5/2._ark)*y4*y2**4+(-sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y4**4*y2
    dF(1,680) = -y2**4*y4**2+y1**2*y3**4+y1**4*y3**2-y2**2*y4**4
    dF(1,681) = -y4**6+y1**6+y3**6-y2**6
    dF(2,1) = 0._ark
    dF(2,2) = y8
    dF(2,3) = y4-y2+y1-y3
    dF(2,4) = y7*y9
    dF(2,5) = y5*y8+sqrt(3._ark)*y6*y8
    dF(2,6) = y1**2-y2**2+y4**2-y3**2
    dF(2,7) = y1*y8+y4*y8+y3*y8+y2*y8
    dF(2,8) = (y9+y7)*y1+(-y9+y7)*y2+(-y7+y9)*y3+(-y9-y7)*y4
    dF(2,9) = (-sqrt(3._ark)*y6/2._ark-y5/2._ark)*y1+(sqrt(3._ark)*y6/2._ark+y5/2._ark)*y2+ &
      (sqrt(3._ark)*y6/2._ark+y5/2._ark)*y3+(-sqrt(3._ark)*y6/2._ark-y5/2._ark)*y4
    dF(2,10) = y1*y4-y2*y3
    dF(2,11) = y8*y9**2+y7**2*y8
    dF(2,12) = y8**3
    dF(2,13) = -sqrt(3._ark)*y6*y7*y9/2._ark-y5*y7*y9/2._ark
    dF(2,14) = -sqrt(3._ark)*y6**2*y8/3._ark-y5*y6*y8
    dF(2,15) = y5**2*y8+y6**2*y8
    dF(2,16) = (-y7*y8-y8*y9)*y1+(y7*y8-y8*y9)*y2+(y8*y9-y7*y8)*y3+(y8*y9+y7*y8)*y4
    dF(2,17) = y3*y7*y9+y4*y7*y9+y2*y7*y9+y1*y7*y9
    dF(2,18) = y3*y8**2+y2*y8**2-y1*y8**2-y4*y8**2
    dF(2,19) = (-y9**2-y7**2)*y1+(y9**2+y7**2)*y2+(y9**2+y7**2)*y3+(-y9**2-y7**2)*y4
    dF(2,20) = (-y6*y8-sqrt(3._ark)*y5*y8/3._ark)*y1+(-y6*y8-sqrt(3._ark)*y5*y8/3._ark)*y2+ &
      (-y6*y8-sqrt(3._ark)*y5*y8/3._ark)*y3+(-y6*y8-sqrt(3._ark)*y5*y8/3._ark)*y4
    dF(2,21) = ((sqrt(3._ark)*y9/2._ark-sqrt(3._ark)*y7/2._ark)*y5+(y7/2._ark-y9/2._ark)*y6)*y1+ &
      ((-sqrt(3._ark)*y9/2._ark-sqrt(3._ark)*y7/2._ark)*y5+(y9/2._ark+y7/2._ark)*y6)*y2+ &
      ((sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y9/2._ark)*y5+(-y9/2._ark-y7/2._ark)*y6)*y3+ &
      ((sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y9/2._ark)*y5+(y9/2._ark-y7/2._ark)*y6)*y4
    dF(2,22) = ((-y9/2._ark-y7/2._ark)*y5+(-sqrt(3._ark)*y9/2._ark-sqrt(3._ark)*y7/ &
      2._ark)*y6)*y1+((y9/2._ark-y7/2._ark)*y5+(sqrt(3._ark)*y9/2._ark-sqrt(3._ark)*y7/ &
      2._ark)*y6)*y2+((y7/2._ark-y9/2._ark)*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y9/ &
      2._ark)*y6)*y3+((y9/2._ark+y7/2._ark)*y5+(sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y9/ &
      2._ark)*y6)*y4
    dF(2,23) = (y5*y6-sqrt(3._ark)*y5**2/3._ark)*y1+(sqrt(3._ark)*y5**2/3._ark-y5*y6)*y2+ &
      (sqrt(3._ark)*y5**2/3._ark-y5*y6)*y3+(y5*y6-sqrt(3._ark)*y5**2/3._ark)*y4
    dF(2,24) = (sqrt(3._ark)*y6/2._ark+y5/2._ark)*y1**2+(-sqrt(3._ark)*y6/2._ark-y5/ &
      2._ark)*y2**2+(-sqrt(3._ark)*y6/2._ark-y5/2._ark)*y3**2+(sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y4**2
    dF(2,25) = y2**3+y3**3-y1**3-y4**3
    dF(2,26) = y2*y3*y8+y1*y4*y8
    dF(2,27) = (y2*y8+y3*y8)*y1+y3*y4*y8+y2*y4*y8
    dF(2,28) = (-y2*y7-y3*y9)*y1+y3*y4*y7+y2*y4*y9
    dF(2,29) = (sqrt(3._ark)*y6/2._ark+y5/2._ark)*y4*y1+(-sqrt(3._ark)*y6/2._ark-y5/ &
      2._ark)*y3*y2
    dF(2,30) = (y3+y2)*y1**2+(-y3**2-y2**2)*y1-y3**2*y4-y2**2*y4+y2*y4**2+y3*y4**2
    dF(2,31) = y2**2*y3+y2*y3**2-y1*y4**2-y1**2*y4
    dF(2,32) = (-y9-y7)*y1**2+(-y7+y9)*y2**2+(-y9+y7)*y3**2+(y9+y7)*y4**2
    dF(2,33) = y3**2*y8+y4**2*y8+y2**2*y8+y1**2*y8
    dF(2,34) = (-y6**2-y5**2)*y1+(y6**2+y5**2)*y2+(y6**2+y5**2)*y3+(-y6**2-y5**2)*y4
    dF(2,35) = ((-y4+y3)*y2-y3*y4)*y1+y2*y3*y4
    dF(2,36) = y7*y8**2*y9
    dF(2,37) = y7*y9**3+y7**3*y9
    dF(2,38) = (sqrt(3._ark)*y7**2*y8-sqrt(3._ark)*y8*y9**2)*y5+(y8*y9**2-y7**2*y8)*y6
    dF(2,39) = (y8*y9**2-2._ark*y7**2*y8)*y5-sqrt(3._ark)*y6*y8*y9**2
    dF(2,40) = y5*y8**3+sqrt(3._ark)*y6*y8**3
    dF(2,41) = y5*y6*y7*y9-sqrt(3._ark)*y5**2*y7*y9/2._ark-sqrt(3._ark)*y6**2*y7*y9/6._ark
    dF(2,42) = y6**2*y7*y9+y5**2*y7*y9
    dF(2,43) = 5._ark/3._ark*sqrt(3._ark)*y5*y6**2*y8-sqrt(3._ark)*y5**3*y8-y5**2*y6*y8- &
      y6**3*y8
    dF(2,44) = y5**3*y8-3._ark*y5*y6**2*y8
    dF(2,45) = y4*y8**3+y2*y8**3+y3*y8**3+y1*y8**3
    dF(2,46) = ((sqrt(3._ark)*y9/6._ark-sqrt(3._ark)*y7/3._ark)*y5**2-y5*y6*y7- &
      sqrt(3._ark)*y6**2*y9/2._ark)*y1+((-sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y9/6._ark)*y5**2- &
      y5*y6*y7+sqrt(3._ark)*y6**2*y9/2._ark)*y2+((sqrt(3._ark)*y7/3._ark+sqrt(3._ark)*y9/ &
      6._ark)*y5**2+y5*y6*y7-sqrt(3._ark)*y6**2*y9/2._ark)*y3+((-sqrt(3._ark)*y9/6._ark+ &
      sqrt(3._ark)*y7/3._ark)*y5**2+y5*y6*y7+sqrt(3._ark)*y6**2*y9/2._ark)*y4
    dF(2,47) = ((sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y9/2._ark)*y5+(y9/2._ark-y7/ &
      2._ark)*y6)*y1**2+((sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y9/2._ark)*y5+(-y9/2._ark-y7/ &
      2._ark)*y6)*y2**2+((-sqrt(3._ark)*y9/2._ark-sqrt(3._ark)*y7/2._ark)*y5+(y9/2._ark+y7/ &
      2._ark)*y6)*y3**2+((sqrt(3._ark)*y9/2._ark-sqrt(3._ark)*y7/2._ark)*y5+(y7/2._ark-y9/ &
      2._ark)*y6)*y4**2
    dF(2,48) = (-y9**2-y7**2)*y1**2+(y9**2+y7**2)*y2**2+(y9**2+y7**2)*y3**2+(-y9**2- &
      y7**2)*y4**2
    dF(2,49) = (-sqrt(3._ark)*y6/2._ark-y5/2._ark)*y1**3+(sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y2**3+(sqrt(3._ark)*y6/2._ark+y5/2._ark)*y3**3+(-sqrt(3._ark)*y6/2._ark-y5/ &
      2._ark)*y4**3
    dF(2,50) = y2**4+y3**4-y4**4-y1**4
    dF(2,51) = (y8*y9**2+y7**2*y8)*y1+(y8*y9**2+y7**2*y8)*y2+(y8*y9**2+y7**2*y8)*y3+ &
      (y8*y9**2+y7**2*y8)*y4
    dF(2,52) = (y7*y9**2+y7**2*y9)*y1+(-y7**2*y9+y7*y9**2)*y2+(y7**2*y9- &
      y7*y9**2)*y3+(-y7**2*y9-y7*y9**2)*y4
    dF(2,53) = (y8**2*y9+y7*y8**2)*y1+(-y8**2*y9+y7*y8**2)*y2+(-y7*y8**2+ &
      y8**2*y9)*y3+(-y7*y8**2-y8**2*y9)*y4
    dF(2,54) = (y7**3+y9**3)*y1+(y7**3-y9**3)*y2+(y9**3-y7**3)*y3+(-y7**3-y9**3)*y4
    dF(2,55) = y1*y7*y8*y9-y2*y7*y8*y9+y4*y7*y8*y9-y3*y7*y8*y9
    dF(2,56) = (sqrt(3._ark)*y5*y7**2/2._ark+(y9**2+y7**2/2._ark)*y6)*y1+(- &
      sqrt(3._ark)*y5*y7**2/2._ark+(-y9**2-y7**2/2._ark)*y6)*y2+(-sqrt(3._ark)*y5*y7**2/2._ark+ &
      (-y9**2-y7**2/2._ark)*y6)*y3+(sqrt(3._ark)*y5*y7**2/2._ark+(y9**2+y7**2/2._ark)*y6)*y4
    dF(2,57) = (-sqrt(3._ark)*y5*y7*y9/3._ark-y6*y7*y9)*y1+(-sqrt(3._ark)*y5*y7*y9/3._ark- &
      y6*y7*y9)*y2+(-sqrt(3._ark)*y5*y7*y9/3._ark-y6*y7*y9)*y3+(-sqrt(3._ark)*y5*y7*y9/3._ark- &
      y6*y7*y9)*y4
    dF(2,58) = (sqrt(3._ark)*y5*y7*y8+(y7*y8+2._ark*y8*y9)*y6)*y1+(-sqrt(3._ark)*y5*y7*y8+ &
      (-y7*y8+2._ark*y8*y9)*y6)*y2+(sqrt(3._ark)*y5*y7*y8+(y7*y8-2._ark*y8*y9)*y6)*y3+(- &
      sqrt(3._ark)*y5*y7*y8+(-2._ark*y8*y9-y7*y8)*y6)*y4
    dF(2,59) = ((2._ark/3._ark*y9-y7/3._ark)*y5**2+2._ark/3._ark*sqrt(3._ark)*y5*y6*y9+ &
      y6**2*y7)*y1+((-y7/3._ark-2._ark/3._ark*y9)*y5**2-2._ark/3._ark*sqrt(3._ark)*y5*y6*y9+ &
      y6**2*y7)*y2+((2._ark/3._ark*y9+y7/3._ark)*y5**2+2._ark/3._ark*sqrt(3._ark)*y5*y6*y9- &
      y6**2*y7)*y3+((y7/3._ark-2._ark/3._ark*y9)*y5**2-2._ark/3._ark*sqrt(3._ark)*y5*y6*y9- &
      y6**2*y7)*y4
    dF(2,60) = ((y9/3._ark+4._ark/3._ark*y7)*y5**2-2._ark/3._ark*sqrt(3._ark)*y5*y6*y9+ &
      y6**2*y9)*y1+((-y9/3._ark+4._ark/3._ark*y7)*y5**2+2._ark/3._ark*sqrt(3._ark)*y5*y6*y9- &
      y6**2*y9)*y2+((-4._ark/3._ark*y7+y9/3._ark)*y5**2-2._ark/3._ark*sqrt(3._ark)*y5*y6*y9+ &
      y6**2*y9)*y3+((-4._ark/3._ark*y7-y9/3._ark)*y5**2+2._ark/3._ark*sqrt(3._ark)*y5*y6*y9- &
      y6**2*y9)*y4
    dF(2,61) = ((y9**2-y7**2/2._ark)*y5+sqrt(3._ark)*y6*y7**2/2._ark)*y1+((y7**2/2._ark- &
      y9**2)*y5-sqrt(3._ark)*y6*y7**2/2._ark)*y2+((y7**2/2._ark-y9**2)*y5- &
      sqrt(3._ark)*y6*y7**2/2._ark)*y3+((y9**2-y7**2/2._ark)*y5+sqrt(3._ark)*y6*y7**2/2._ark)*y4
    dF(2,62) = (-y5*y8**2/2._ark-sqrt(3._ark)*y6*y8**2/2._ark)*y1+(y5*y8**2/2._ark+ &
      sqrt(3._ark)*y6*y8**2/2._ark)*y2+(y5*y8**2/2._ark+sqrt(3._ark)*y6*y8**2/2._ark)*y3+(- &
      y5*y8**2/2._ark-sqrt(3._ark)*y6*y8**2/2._ark)*y4
    dF(2,63) = ((y8*y9-2._ark*y7*y8)*y5-sqrt(3._ark)*y6*y8*y9)*y1+((2._ark*y7*y8+ &
      y8*y9)*y5-sqrt(3._ark)*y6*y8*y9)*y2+((-2._ark*y7*y8-y8*y9)*y5+ &
      sqrt(3._ark)*y6*y8*y9)*y3+((2._ark*y7*y8-y8*y9)*y5+sqrt(3._ark)*y6*y8*y9)*y4
    dF(2,64) = (-sqrt(3._ark)*y6**2*y8/3._ark-y5*y6*y8)*y1+(-sqrt(3._ark)*y6**2*y8/3._ark- &
      y5*y6*y8)*y2+(-sqrt(3._ark)*y6**2*y8/3._ark-y5*y6*y8)*y3+(-sqrt(3._ark)*y6**2*y8/3._ark- &
      y5*y6*y8)*y4
    dF(2,65) = (y5*y6**2-y5**3/3._ark)*y1+(y5**3/3._ark-y5*y6**2)*y2+(y5**3/3._ark- &
      y5*y6**2)*y3+(y5*y6**2-y5**3/3._ark)*y4
    dF(2,66) = (y5**2*y8+y6**2*y8)*y1+(y5**2*y8+y6**2*y8)*y2+(y5**2*y8+y6**2*y8)*y3+ &
      (y5**2*y8+y6**2*y8)*y4
    dF(2,67) = (sqrt(3._ark)*y5*y6**2+sqrt(3._ark)*y5**3/9._ark+y6**3+y5**2*y6)*y1+(- &
      y5**2*y6-sqrt(3._ark)*y5*y6**2-y6**3-sqrt(3._ark)*y5**3/9._ark)*y2+(-y5**2*y6- &
      sqrt(3._ark)*y5*y6**2-y6**3-sqrt(3._ark)*y5**3/9._ark)*y3+(sqrt(3._ark)*y5*y6**2+ &
      sqrt(3._ark)*y5**3/9._ark+y6**3+y5**2*y6)*y4
    dF(2,68) = (y5*y6+sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y4*y1+(- &
      sqrt(3._ark)*y6**2/2._ark-sqrt(3._ark)*y5**2/6._ark-y5*y6)*y3*y2
    dF(2,69) = (y6+sqrt(3._ark)*y5/3._ark)*y4*y1**2+(y6+sqrt(3._ark)*y5/3._ark)*y4**2*y1+(- &
      sqrt(3._ark)*y5/3._ark-y6)*y3*y2**2+(-sqrt(3._ark)*y5/3._ark-y6)*y3**2*y2
    dF(2,70) = y1**3*y4+y1*y4**3-y2*y3**3-y2**3*y3
    dF(2,71) = (-y3-y2)*y1**3+(y3**3+y2**3)*y1+y2**3*y4-y2*y4**3-y3*y4**3+y3**3*y4
    dF(2,72) = y1*y4*y8**2-y2*y3*y8**2
    dF(2,73) = (y3*y7*y8+y2*y8*y9)*y1-y3*y4*y8*y9-y2*y4*y7*y8
    dF(2,74) = (y9**2+y7**2)*y4*y1+(-y9**2-y7**2)*y3*y2
    dF(2,75) = (y3*y7*y9+y2*y7*y9)*y1+y3*y4*y7*y9+y2*y4*y7*y9
    dF(2,76) = y2*y3*y7*y9+y1*y4*y7*y9
    dF(2,77) = (-y2*y6*y8+(-y6*y8/2._ark-sqrt(3._ark)*y5*y8/2._ark)*y3)*y1+(-y6*y8/2._ark- &
      sqrt(3._ark)*y5*y8/2._ark)*y4*y2-y3*y4*y6*y8
    dF(2,78) = (y2*y6*y7+(y6*y9/2._ark+sqrt(3._ark)*y5*y9/2._ark)*y3)*y1+(- &
      sqrt(3._ark)*y5*y9/2._ark-y6*y9/2._ark)*y4*y2-y3*y4*y6*y7
    dF(2,79) = (y6**2+y5**2)*y4*y1+(-y6**2-y5**2)*y3*y2
    dF(2,80) = (-y5*y8/2._ark-sqrt(3._ark)*y6*y8/2._ark)*y4*y1+(-y5*y8/2._ark- &
      sqrt(3._ark)*y6*y8/2._ark)*y3*y2
    dF(2,81) = (y2*y5*y8+(sqrt(3._ark)*y6*y8/2._ark-y5*y8/2._ark)*y3)*y1+ &
      (sqrt(3._ark)*y6*y8/2._ark-y5*y8/2._ark)*y4*y2+y3*y4*y5*y8
    dF(2,82) = (y2*y5*y7+(sqrt(3._ark)*y6*y9/2._ark-y5*y9/2._ark)*y3)*y1+(- &
      sqrt(3._ark)*y6*y9/2._ark+y5*y9/2._ark)*y4*y2-y3*y4*y5*y7
    dF(2,83) = ((y4*y8+y3*y8)*y2+y3*y4*y8)*y1+y2*y3*y4*y8
    dF(2,84) = (((sqrt(3._ark)*y6/2._ark+y5/2._ark)*y3+(-sqrt(3._ark)*y6/2._ark-y5/ &
      2._ark)*y4)*y2+(-sqrt(3._ark)*y6/2._ark-y5/2._ark)*y4*y3)*y1+(sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y4*y3*y2
    dF(2,85) = (((y9+y7)*y3+(-y9+y7)*y4)*y2+(-y7+y9)*y4*y3)*y1+(-y9-y7)*y4*y3*y2
    dF(2,86) = (-y2*y5+(-sqrt(3._ark)*y6/2._ark+y5/2._ark)*y3)*y1**2+(y2**2*y5+ &
      (sqrt(3._ark)*y6/2._ark-y5/2._ark)*y3**2)*y1+(sqrt(3._ark)*y6/2._ark-y5/2._ark)*y4*y2**2+(- &
      sqrt(3._ark)*y6/2._ark+y5/2._ark)*y4**2*y2-y3*y4**2*y5+y3**2*y4*y5
    dF(2,87) = (y8*y9+y7*y8)*y1**2+(y8*y9-y7*y8)*y2**2+(y7*y8-y8*y9)*y3**2+(-y7*y8- &
      y8*y9)*y4**2
    dF(2,88) = y1**2*y7*y9+y4**2*y7*y9+y3**2*y7*y9+y2**2*y7*y9
    dF(2,89) = y4**2*y8**2+y1**2*y8**2-y2**2*y8**2-y3**2*y8**2
    dF(2,90) = (-y6*y8-sqrt(3._ark)*y5*y8/3._ark)*y1**2+(-y6*y8-sqrt(3._ark)*y5*y8/ &
      3._ark)*y2**2+(-y6*y8-sqrt(3._ark)*y5*y8/3._ark)*y3**2+(-y6*y8-sqrt(3._ark)*y5*y8/ &
      3._ark)*y4**2
    dF(2,91) = (y6**2+y5**2)*y1**2+(-y6**2-y5**2)*y2**2+(-y6**2-y5**2)*y3**2+(y6**2+ &
      y5**2)*y4**2
    dF(2,92) = ((-y9/2._ark-y7/2._ark)*y5+(-sqrt(3._ark)*y9/2._ark-sqrt(3._ark)*y7/ &
      2._ark)*y6)*y1**2+((y9/2._ark-y7/2._ark)*y5+(sqrt(3._ark)*y9/2._ark-sqrt(3._ark)*y7/ &
      2._ark)*y6)*y2**2+((y7/2._ark-y9/2._ark)*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y9/ &
      2._ark)*y6)*y3**2+((y9/2._ark+y7/2._ark)*y5+(sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y9/ &
      2._ark)*y6)*y4**2
    dF(2,93) = (y5*y6-sqrt(3._ark)*y5**2/3._ark)*y1**2+(sqrt(3._ark)*y5**2/3._ark- &
      y5*y6)*y2**2+(sqrt(3._ark)*y5**2/3._ark-y5*y6)*y3**2+(y5*y6-sqrt(3._ark)*y5**2/ &
      3._ark)*y4**2
    dF(2,94) = (y3*y9+y2*y7)*y1**2+(y3**2*y9+y2**2*y7)*y1-y2*y4**2*y9-y3**2*y4*y7- &
      y3*y4**2*y7-y2**2*y4*y9
    dF(2,95) = (y9+y7)*y4*y1**2+(-y9-y7)*y4**2*y1+(-y9+y7)*y3*y2**2+(-y7+ &
      y9)*y3**2*y2
    dF(2,96) = (y2*y8+y3*y8)*y1**2+(y3**2*y8+y2**2*y8)*y1+y3*y4**2*y8+y3**2*y4*y8+ &
      y2*y4**2*y8+y2**2*y4*y8
    dF(2,97) = y2*y3**2*y8+y2**2*y3*y8+y1*y4**2*y8+y1**2*y4*y8
    dF(2,98) = (y3*y7+y2*y9)*y1**2+(-y3**2*y7-y2**2*y9)*y1+y2**2*y4*y7-y3*y4**2*y9- &
      y2*y4**2*y7+y3**2*y4*y9
    dF(2,99) = y1**2*y4**2-y2**2*y3**2
    dF(2,100) = (-y2*y6+(-y6/2._ark-sqrt(3._ark)*y5/2._ark)*y3)*y1**2+(y2**2*y6+ &
      (sqrt(3._ark)*y5/2._ark+y6/2._ark)*y3**2)*y1+(sqrt(3._ark)*y5/2._ark+y6/2._ark)*y4*y2**2+(- &
      y6/2._ark-sqrt(3._ark)*y5/2._ark)*y4**2*y2+y3**2*y4*y6-y3*y4**2*y6
    dF(2,101) = (y2*y4+y3*y4)*y1**2+(-y2**2*y3+(-y3**2+y4**2)*y2+y3*y4**2)*y1- &
      y2*y3**2*y4-y2**2*y3*y4
    dF(2,102) = y1**2*y2*y3+(-y3**2*y4-y2**2*y4)*y1+y2*y3*y4**2
    dF(2,103) = (y9+y7)*y1**3+(-y9+y7)*y2**3+(-y7+y9)*y3**3+(-y9-y7)*y4**3
    dF(2,104) = y4**3*y8+y3**3*y8+y2**3*y8+y1**3*y8
    dF(2,105) = y7**2*y8**3+y8**3*y9**2
    dF(2,106) = y7**2*y8*y9**2
    dF(2,107) = y7**4*y8+y8*y9**4
    dF(2,108) = y8**5
    dF(2,109) = -sqrt(3._ark)*y5*y7**3*y9+(-2._ark*y7*y9**3-y7**3*y9)*y6
    dF(2,110) = 3._ark*y5**2*y8*y9**2+2._ark*sqrt(3._ark)*y5*y6*y8*y9**2+(y8*y9**2+ &
      4._ark*y7**2*y8)*y6**2
    dF(2,111) = 4._ark/3._ark*sqrt(3._ark)*y5**3*y6*y8+y6**4*y8+y5**4*y8/3._ark
    dF(2,112) = (y7*y9**3+y7**3*y9)*y5+(sqrt(3._ark)*y7*y9**3+sqrt(3._ark)*y7**3*y9)*y6
    dF(2,113) = y5*y7*y8**2*y9+sqrt(3._ark)*y6*y7*y8**2*y9
    dF(2,114) = -sqrt(3._ark)*y5**2*y8*y9**2+(y7**2*y8-y8*y9**2)*y6*y5- &
      sqrt(3._ark)*y6**2*y7**2*y8
    dF(2,115) = -sqrt(3._ark)*y6**2*y8**3/3._ark-y5*y6*y8**3
    dF(2,116) = -y5*y6**3*y8+2._ark/9._ark*sqrt(3._ark)*y5**4*y8-y5**3*y6*y8/3._ark
    dF(2,117) = (-2._ark*y8*y9**2+y7**2*y8)*y5**2-2._ark*sqrt(3._ark)*y5*y6*y8*y9**2- &
      3._ark*y6**2*y7**2*y8
    dF(2,118) = y6**2*y8**3+y5**2*y8**3
    dF(2,119) = -5._ark/3._ark*sqrt(3._ark)*y5*y6**2*y7*y9+sqrt(3._ark)*y5**3*y7*y9+ &
      y5**2*y6*y7*y9+y6**3*y7*y9
    dF(2,120) = -2._ark/3._ark*sqrt(3._ark)*y5**3*y6*y8+y5**2*y6**2*y8+y5**4*y8/3._ark
    dF(2,121) = -3._ark*y5*y6**2*y7*y9+y5**3*y7*y9
    dF(2,122) = (y9**4+y7**4)*y1+(-y9**4-y7**4)*y2+(-y9**4-y7**4)*y3+(y9**4+ &
      y7**4)*y4
    dF(2,123) = (-y7**3*y8-y8*y9**3)*y1+(y7**3*y8-y8*y9**3)*y2+(y8*y9**3- &
      y7**3*y8)*y3+(y8*y9**3+y7**3*y8)*y4
    dF(2,124) = (sqrt(3._ark)*y5*y9**3/2._ark+(y7**3+y9**3/2._ark)*y6)*y1+(- &
      sqrt(3._ark)*y5*y9**3/2._ark+(y7**3-y9**3/2._ark)*y6)*y2+(sqrt(3._ark)*y5*y9**3/2._ark+ &
      (y9**3/2._ark-y7**3)*y6)*y3+(-sqrt(3._ark)*y5*y9**3/2._ark+(-y9**3/2._ark-y7**3)*y6)*y4
    dF(2,125) = (sqrt(3._ark)*y5**3*y9/2._ark+(y7+7._ark/2._ark*y9)*y6*y5**2+ &
      (2._ark*sqrt(3._ark)*y7+5._ark/2._ark*sqrt(3._ark)*y9)*y6**2*y5+(3._ark*y7+3._ark/ &
      2._ark*y9)*y6**3)*y1+(-sqrt(3._ark)*y5**3*y9/2._ark+(-7._ark/2._ark*y9+y7)*y6*y5**2+(- &
      5._ark/2._ark*sqrt(3._ark)*y9+2._ark*sqrt(3._ark)*y7)*y6**2*y5+(-3._ark/2._ark*y9+ &
      3._ark*y7)*y6**3)*y2+(sqrt(3._ark)*y5**3*y9/2._ark+(-y7+7._ark/2._ark*y9)*y6*y5**2+(5._ark/ &
      2._ark*sqrt(3._ark)*y9-2._ark*sqrt(3._ark)*y7)*y6**2*y5+(3._ark/2._ark*y9- &
      3._ark*y7)*y6**3)*y3+(-sqrt(3._ark)*y5**3*y9/2._ark+(-7._ark/2._ark*y9-y7)*y6*y5**2+(- &
      2._ark*sqrt(3._ark)*y7-5._ark/2._ark*sqrt(3._ark)*y9)*y6**2*y5+(-3._ark*y7-3._ark/ &
      2._ark*y9)*y6**3)*y4
    dF(2,126) = y4**3*y7*y9+y3**3*y7*y9+y2**3*y7*y9+y1**3*y7*y9
    dF(2,127) = ((y9**3/2._ark-y7**3)*y5-sqrt(3._ark)*y6*y9**3/2._ark)*y1+((-y9**3/2._ark- &
      y7**3)*y5+sqrt(3._ark)*y6*y9**3/2._ark)*y2+((y7**3+y9**3/2._ark)*y5- &
      sqrt(3._ark)*y6*y9**3/2._ark)*y3+((y7**3-y9**3/2._ark)*y5+sqrt(3._ark)*y6*y9**3/2._ark)*y4
    dF(2,128) = ((-y9**2-y7**2)*y5**2+(-y9**2-y7**2)*y6**2)*y1+((y9**2+y7**2)*y5**2+ &
      (y9**2+y7**2)*y6**2)*y2+((y9**2+y7**2)*y5**2+(y9**2+y7**2)*y6**2)*y3+((-y9**2- &
      y7**2)*y5**2+(-y9**2-y7**2)*y6**2)*y4
    dF(2,129) = (y5**3*y8-3._ark*y5*y6**2*y8)*y1+(y5**3*y8-3._ark*y5*y6**2*y8)*y2+ &
      (y5**3*y8-3._ark*y5*y6**2*y8)*y3+(y5**3*y8-3._ark*y5*y6**2*y8)*y4
    dF(2,130) = ((-y7+7._ark/2._ark*y9)*y5**3+15._ark/2._ark*sqrt(3._ark)*y5**2*y6*y9+(27._ark/ &
      2._ark*y9+9._ark*y7)*y6**2*y5+(6._ark*sqrt(3._ark)*y7+3._ark/ &
      2._ark*sqrt(3._ark)*y9)*y6**3)*y1+((-7._ark/2._ark*y9-y7)*y5**3-15._ark/ &
      2._ark*sqrt(3._ark)*y5**2*y6*y9+(9._ark*y7-27._ark/2._ark*y9)*y6**2*y5+ &
      (6._ark*sqrt(3._ark)*y7-3._ark/2._ark*sqrt(3._ark)*y9)*y6**3)*y2+((y7+7._ark/2._ark*y9)*y5**3+ &
      15._ark/2._ark*sqrt(3._ark)*y5**2*y6*y9+(27._ark/2._ark*y9-9._ark*y7)*y6**2*y5+(3._ark/ &
      2._ark*sqrt(3._ark)*y9-6._ark*sqrt(3._ark)*y7)*y6**3)*y3+((-7._ark/2._ark*y9+y7)*y5**3- &
      15._ark/2._ark*sqrt(3._ark)*y5**2*y6*y9+(-9._ark*y7-27._ark/2._ark*y9)*y6**2*y5+(-3._ark/ &
      2._ark*sqrt(3._ark)*y9-6._ark*sqrt(3._ark)*y7)*y6**3)*y4
    dF(2,131) = ((-y7**2+2._ark*y9**2)*y5+sqrt(3._ark)*y6*y7**2)*y4*y1+((y7**2- &
      2._ark*y9**2)*y5-sqrt(3._ark)*y6*y7**2)*y3*y2
    dF(2,132) = (y2*y5+(sqrt(3._ark)*y6/2._ark-y5/2._ark)*y3)*y1**3+(-y2**3*y5+(- &
      sqrt(3._ark)*y6/2._ark+y5/2._ark)*y3**3)*y1+(-sqrt(3._ark)*y6/2._ark+y5/2._ark)*y4*y2**3+ &
      (sqrt(3._ark)*y6/2._ark-y5/2._ark)*y4**3*y2+y3*y4**3*y5-y3**3*y4*y5
    dF(2,133) = (sqrt(3._ark)*y6/2._ark+y5/2._ark)*y4*y1**3+(sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y4**3*y1+(-sqrt(3._ark)*y6/2._ark-y5/2._ark)*y3*y2**3+(-sqrt(3._ark)*y6/2._ark-y5/ &
      2._ark)*y3**3*y2
    dF(2,134) = ((y9+y7)*y5**2+(y9+y7)*y6**2)*y1**2+((-y9+y7)*y5**2+(-y9+ &
      y7)*y6**2)*y2**2+((-y7+y9)*y5**2+(-y7+y9)*y6**2)*y3**2+((-y9-y7)*y5**2+(-y9- &
      y7)*y6**2)*y4**2
    dF(2,135) = ((sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y9/2._ark)*y5+(y9/2._ark-y7/ &
      2._ark)*y6)*y1**3+((sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y9/2._ark)*y5+(-y9/2._ark-y7/ &
      2._ark)*y6)*y2**3+((-sqrt(3._ark)*y9/2._ark-sqrt(3._ark)*y7/2._ark)*y5+(y9/2._ark+y7/ &
      2._ark)*y6)*y3**3+((sqrt(3._ark)*y9/2._ark-sqrt(3._ark)*y7/2._ark)*y5+(y7/2._ark-y9/ &
      2._ark)*y6)*y4**3
    dF(2,136) = (y6**2+y5**2)*y1**3+(-y6**2-y5**2)*y2**3+(-y6**2-y5**2)*y3**3+ &
      (y6**2+y5**2)*y4**3
    dF(2,137) = (sqrt(3._ark)*y5*y7**2*y9/2._ark+(y7**2*y9/2._ark+y7*y9**2)*y6)*y1+(- &
      sqrt(3._ark)*y5*y7**2*y9/2._ark+(y7*y9**2-y7**2*y9/2._ark)*y6)*y2+ &
      (sqrt(3._ark)*y5*y7**2*y9/2._ark+(y7**2*y9/2._ark-y7*y9**2)*y6)*y3+(- &
      sqrt(3._ark)*y5*y7**2*y9/2._ark+(-y7*y9**2-y7**2*y9/2._ark)*y6)*y4
    dF(2,138) = (-y6*y7*y8*y9-sqrt(3._ark)*y5*y7*y8*y9/3._ark)*y1+(y6*y7*y8*y9+ &
      sqrt(3._ark)*y5*y7*y8*y9/3._ark)*y2+(y6*y7*y8*y9+sqrt(3._ark)*y5*y7*y8*y9/3._ark)*y3+(- &
      y6*y7*y8*y9-sqrt(3._ark)*y5*y7*y8*y9/3._ark)*y4
    dF(2,139) = (-sqrt(3._ark)*y5**2*y9**2/2._ark+y5*y6*y7**2+(sqrt(3._ark)*y9**2/6._ark- &
      sqrt(3._ark)*y7**2/3._ark)*y6**2)*y1+(sqrt(3._ark)*y5**2*y9**2/2._ark-y5*y6*y7**2+(- &
      sqrt(3._ark)*y9**2/6._ark+sqrt(3._ark)*y7**2/3._ark)*y6**2)*y2+(sqrt(3._ark)*y5**2*y9**2/ &
      2._ark-y5*y6*y7**2+(-sqrt(3._ark)*y9**2/6._ark+sqrt(3._ark)*y7**2/3._ark)*y6**2)*y3+(- &
      sqrt(3._ark)*y5**2*y9**2/2._ark+y5*y6*y7**2+(sqrt(3._ark)*y9**2/6._ark-sqrt(3._ark)*y7**2/ &
      3._ark)*y6**2)*y4
    dF(2,140) = (y3*y8*y9+y2*y7*y8)*y1**2+(-y2**2*y7*y8-y3**2*y8*y9)*y1- &
      y3*y4**2*y7*y8+y2**2*y4*y8*y9-y2*y4**2*y8*y9+y3**2*y4*y7*y8
    dF(2,141) = (-y9**2-y7**2)*y4*y1**2+(-y9**2-y7**2)*y4**2*y1+(y9**2+ &
      y7**2)*y3*y2**2+(y9**2+y7**2)*y3**2*y2
    dF(2,142) = (y3*y8**2+y2*y8**2)*y1**2+(-y2**2*y8**2-y3**2*y8**2)*y1+ &
      y2*y4**2*y8**2-y2**2*y4*y8**2-y3**2*y4*y8**2+y3*y4**2*y8**2
    dF(2,143) = (-y2*y7**2-y3*y9**2)*y1**2+(y3**2*y9**2+y2**2*y7**2)*y1+ &
      y2**2*y4*y9**2-y3*y4**2*y7**2-y2*y4**2*y9**2+y3**2*y4*y7**2
    dF(2,144) = (y8*y9+y7*y8)*y4*y1**2+(-y7*y8-y8*y9)*y4**2*y1+(y8*y9- &
      y7*y8)*y3*y2**2+(y7*y8-y8*y9)*y3**2*y2
    dF(2,145) = (-y3*y7*y8-y2*y8*y9)*y1**2+(-y2**2*y8*y9-y3**2*y7*y8)*y1+ &
      y2*y4**2*y7*y8+y3**2*y4*y8*y9+y2**2*y4*y7*y8+y3*y4**2*y8*y9
    dF(2,146) = (y3*y7*y9+y2*y7*y9)*y1**2+(y2**2*y7*y9+y3**2*y7*y9)*y1+ &
      y2*y4**2*y7*y9+y3**2*y4*y7*y9+y3*y4**2*y7*y9+y2**2*y4*y7*y9
    dF(2,147) = (-2._ark*y2*y6*y9+(-y6*y7-sqrt(3._ark)*y5*y7)*y3)*y1**2+ &
      (2._ark*y2**2*y6*y9+(sqrt(3._ark)*y5*y7+y6*y7)*y3**2)*y1+(-y6*y7- &
      sqrt(3._ark)*y5*y7)*y4*y2**2+(sqrt(3._ark)*y5*y7+y6*y7)*y4**2*y2- &
      2._ark*y3**2*y4*y6*y9+2._ark*y3*y4**2*y6*y9
    dF(2,148) = ((sqrt(3._ark)*y6*y8/2._ark-y5*y8/2._ark)*y2+y3*y5*y8)*y1**2+ &
      ((sqrt(3._ark)*y6*y8/2._ark-y5*y8/2._ark)*y2**2+y3**2*y5*y8)*y1+y2**2*y4*y5*y8+ &
      y2*y4**2*y5*y8+(sqrt(3._ark)*y6*y8/2._ark-y5*y8/2._ark)*y4*y3**2+(sqrt(3._ark)*y6*y8/ &
      2._ark-y5*y8/2._ark)*y4**2*y3
    dF(2,149) = ((-sqrt(3._ark)*y5**2/2._ark+y5*y6-sqrt(3._ark)*y6**2/6._ark)*y2+(- &
      sqrt(3._ark)*y5**2/2._ark+y5*y6-sqrt(3._ark)*y6**2/6._ark)*y3)*y1**2+((sqrt(3._ark)*y5**2/ &
      2._ark-y5*y6+sqrt(3._ark)*y6**2/6._ark)*y2**2+(sqrt(3._ark)*y5**2/2._ark-y5*y6+ &
      sqrt(3._ark)*y6**2/6._ark)*y3**2)*y1+(sqrt(3._ark)*y5**2/2._ark-y5*y6+sqrt(3._ark)*y6**2/ &
      6._ark)*y4*y2**2+(-sqrt(3._ark)*y5**2/2._ark+y5*y6-sqrt(3._ark)*y6**2/6._ark)*y4**2*y2+ &
      (sqrt(3._ark)*y5**2/2._ark-y5*y6+sqrt(3._ark)*y6**2/6._ark)*y4*y3**2+(-sqrt(3._ark)*y5**2/ &
      2._ark+y5*y6-sqrt(3._ark)*y6**2/6._ark)*y4**2*y3
    dF(2,150) = ((y5*y6-sqrt(3._ark)*y6**2/3._ark)*y2+(sqrt(3._ark)*y6**2/6._ark- &
      sqrt(3._ark)*y5**2/2._ark)*y3)*y1**2+((sqrt(3._ark)*y6**2/3._ark-y5*y6)*y2**2+ &
      (sqrt(3._ark)*y5**2/2._ark-sqrt(3._ark)*y6**2/6._ark)*y3**2)*y1+(sqrt(3._ark)*y5**2/2._ark- &
      sqrt(3._ark)*y6**2/6._ark)*y4*y2**2+(sqrt(3._ark)*y6**2/6._ark-sqrt(3._ark)*y5**2/ &
      2._ark)*y4**2*y2+(sqrt(3._ark)*y6**2/3._ark-y5*y6)*y4*y3**2+(y5*y6-sqrt(3._ark)*y6**2/ &
      3._ark)*y4**2*y3
    dF(2,151) = (-y6**2-y5**2)*y4*y1**2+(-y6**2-y5**2)*y4**2*y1+(y6**2+ &
      y5**2)*y3*y2**2+(y6**2+y5**2)*y3**2*y2
    dF(2,152) = (y9+y7)*y4*y1**3+(-y9-y7)*y4**3*y1+(-y9+y7)*y3*y2**3+(-y7+ &
      y9)*y3**3*y2
    dF(2,153) = (-y2*y7-y3*y9)*y1**3+(-y3**3*y9-y2**3*y7)*y1+y3*y4**3*y7+ &
      y3**3*y4*y7+y2*y4**3*y9+y2**3*y4*y9
    dF(2,154) = y2*y3**3*y8+y2**3*y3*y8+y1*y4**3*y8+y1**3*y4*y8
    dF(2,155) = (y3*y7+y2*y9)*y1**3+(-y2**3*y9-y3**3*y7)*y1-y2*y4**3*y7+y3**3*y4*y9- &
      y3*y4**3*y9+y2**3*y4*y7
    dF(2,156) = (-y3-y2)*y1**4+(y3**4+y2**4)*y1+y2**4*y4-y3*y4**4+y3**4*y4-y2*y4**4
    dF(2,157) = (y2*y8**3+y3*y8**3)*y1+y3*y4*y8**3+y2*y4*y8**3
    dF(2,158) = (y8*y9**2+y7**2*y8)*y4*y1+(y8*y9**2+y7**2*y8)*y3*y2
    dF(2,159) = (-y2*y7*y9**2-y3*y7**2*y9)*y1+y3*y4*y7*y9**2+y2*y4*y7**2*y9
    dF(2,160) = y2*y3*y8**3+y1*y4*y8**3
    dF(2,161) = y2*y3*y7*y8*y9-y1*y4*y7*y8*y9
    dF(2,162) = (y3*y7**2*y8+y2*y8*y9**2)*y1+y3*y4*y8*y9**2+y2*y4*y7**2*y8
    dF(2,163) = (-y2*y7*y8**2-y3*y8**2*y9)*y1+y3*y4*y7*y8**2+y2*y4*y8**2*y9
    dF(2,164) = (y2*y7**2*y8+y3*y8*y9**2)*y1+y3*y4*y7**2*y8+y2*y4*y8*y9**2
    dF(2,165) = (-y2*y7**3-y3*y9**3)*y1+y3*y4*y7**3+y2*y4*y9**3
    dF(2,166) = ((sqrt(3._ark)*y9**2-sqrt(3._ark)*y7**2)*y5+(y7**2-y9**2)*y6)*y4*y1+ &
      ((sqrt(3._ark)*y7**2-sqrt(3._ark)*y9**2)*y5+(-y7**2+y9**2)*y6)*y3*y2
    dF(2,167) = ((sqrt(3._ark)*y5*y8*y9+y6*y8*y9)*y2+2._ark*y3*y6*y7*y8)*y1- &
      2._ark*y2*y4*y6*y7*y8+(-y6*y8*y9-sqrt(3._ark)*y5*y8*y9)*y4*y3
    dF(2,168) = ((-y5*y8*y9-sqrt(3._ark)*y6*y8*y9)*y2+(-sqrt(3._ark)*y6*y7*y8- &
      y5*y7*y8)*y3)*y1+(sqrt(3._ark)*y6*y7*y8+y5*y7*y8)*y4*y2+(sqrt(3._ark)*y6*y8*y9+ &
      y5*y8*y9)*y4*y3
    dF(2,169) = ((-y5*y6*y8/2._ark-sqrt(3._ark)*y6**2*y8/4._ark+sqrt(3._ark)*y5**2*y8/ &
      4._ark)*y2-y3*y5*y6*y8)*y1-y2*y4*y5*y6*y8+(-y5*y6*y8/2._ark-sqrt(3._ark)*y6**2*y8/ &
      4._ark+sqrt(3._ark)*y5**2*y8/4._ark)*y4*y3
    dF(2,170) = ((-y5*y6*y7-sqrt(3._ark)*y5**2*y7/3._ark)*y2+(-sqrt(3._ark)*y6**2*y9/2._ark+ &
      sqrt(3._ark)*y5**2*y9/6._ark)*y3)*y1+(sqrt(3._ark)*y6**2*y9/2._ark-sqrt(3._ark)*y5**2*y9/ &
      6._ark)*y4*y2+(sqrt(3._ark)*y5**2*y7/3._ark+y5*y6*y7)*y4*y3
    dF(2,171) = (y6+sqrt(3._ark)*y5/3._ark)*y3*y2*y1**2+((-sqrt(3._ark)*y5/3._ark- &
      y6)*y4*y2**2+(-sqrt(3._ark)*y5/3._ark-y6)*y4*y3**2)*y1+(y6+sqrt(3._ark)*y5/ &
      3._ark)*y4**2*y3*y2
    dF(2,172) = (y3*y4*y8+y2*y4*y8)*y1**2+(y2**2*y3*y8+(y4**2*y8+y3**2*y8)*y2+ &
      y3*y4**2*y8)*y1+y2*y3**2*y4*y8+y2**2*y3*y4*y8
    dF(2,173) = (y3*y4*y7+y2*y4*y9)*y1**2+(-y2**2*y3*y9+(-y3**2*y7-y4**2*y7)*y2- &
      y3*y4**2*y9)*y1+y2*y3**2*y4*y9+y2**2*y3*y4*y7
    dF(2,174) = ((-y6+sqrt(3._ark)*y5)*y4*y2+(-sqrt(3._ark)*y5+y6)*y4*y3)*y1**2+((- &
      sqrt(3._ark)*y5+y6)*y3*y2**2+((-y6+sqrt(3._ark)*y5)*y3**2+(-sqrt(3._ark)*y5+ &
      y6)*y4**2)*y2+(-y6+sqrt(3._ark)*y5)*y4**2*y3)*y1+(-y6+sqrt(3._ark)*y5)*y4*y3*y2**2+ &
      (-sqrt(3._ark)*y5+y6)*y4*y3**2*y2
    dF(2,175) = y1**3*y2*y3+(-y3**3*y4-y2**3*y4)*y1+y2*y3*y4**3
    dF(2,176) = (-y3*y4-y2*y4)*y1**3+(y2**3*y3+(y3**3-y4**3)*y2-y3*y4**3)*y1+ &
      y2**3*y3*y4+y2*y3**3*y4
    dF(2,177) = (y7*y9**2+y7**2*y9)*y1**2+(-y7**2*y9+y7*y9**2)*y2**2+(y7**2*y9- &
      y7*y9**2)*y3**2+(-y7**2*y9-y7*y9**2)*y4**2
    dF(2,178) = (y8*y9**2+y7**2*y8)*y1**2+(y8*y9**2+y7**2*y8)*y2**2+(y8*y9**2+ &
      y7**2*y8)*y3**2+(y8*y9**2+y7**2*y8)*y4**2
    dF(2,179) = (-y7*y8**2-y8**2*y9)*y1**2+(-y7*y8**2+y8**2*y9)*y2**2+(-y8**2*y9+ &
      y7*y8**2)*y3**2+(y8**2*y9+y7*y8**2)*y4**2
    dF(2,180) = (sqrt(3._ark)*y5*y7*y8/2._ark+(y8*y9+y7*y8/2._ark)*y6)*y1**2+(- &
      sqrt(3._ark)*y5*y7*y8/2._ark+(y8*y9-y7*y8/2._ark)*y6)*y2**2+(sqrt(3._ark)*y5*y7*y8/2._ark+ &
      (y7*y8/2._ark-y8*y9)*y6)*y3**2+(-sqrt(3._ark)*y5*y7*y8/2._ark+(-y7*y8/2._ark- &
      y8*y9)*y6)*y4**2
    dF(2,181) = ((y7**2-2._ark*y9**2)*y5-sqrt(3._ark)*y6*y7**2)*y1**2+((-y7**2+ &
      2._ark*y9**2)*y5+sqrt(3._ark)*y6*y7**2)*y2**2+((-y7**2+2._ark*y9**2)*y5+ &
      sqrt(3._ark)*y6*y7**2)*y3**2+((y7**2-2._ark*y9**2)*y5-sqrt(3._ark)*y6*y7**2)*y4**2
    dF(2,182) = ((y8*y9-y7*y8/2._ark)*y5+sqrt(3._ark)*y6*y7*y8/2._ark)*y1**2+((y8*y9+ &
      y7*y8/2._ark)*y5-sqrt(3._ark)*y6*y7*y8/2._ark)*y2**2+((-y7*y8/2._ark-y8*y9)*y5+ &
      sqrt(3._ark)*y6*y7*y8/2._ark)*y3**2+((y7*y8/2._ark-y8*y9)*y5-sqrt(3._ark)*y6*y7*y8/ &
      2._ark)*y4**2
    dF(2,183) = y2**2*y3*y7*y9+y2*y3**2*y7*y9+y1*y4**2*y7*y9+y1**2*y4*y7*y9
    dF(2,184) = (-2._ark*y2*y6*y7+(-sqrt(3._ark)*y5*y9-y6*y9)*y3)*y1**2+(- &
      2._ark*y2**2*y6*y7+(-sqrt(3._ark)*y5*y9-y6*y9)*y3**2)*y1+(y6*y9+ &
      sqrt(3._ark)*y5*y9)*y4*y2**2+(y6*y9+sqrt(3._ark)*y5*y9)*y4**2*y2+ &
      2._ark*y3**2*y4*y6*y7+2._ark*y3*y4**2*y6*y7
    dF(2,185) = y2**3*y8**2+y3**3*y8**2-y1**3*y8**2-y4**3*y8**2
    dF(2,186) = (y9**2+y7**2)*y1**3+(-y9**2-y7**2)*y2**3+(-y9**2-y7**2)*y3**3+ &
      (y9**2+y7**2)*y4**3
    dF(2,187) = (-y7*y8-y8*y9)*y1**3+(y7*y8-y8*y9)*y2**3+(y8*y9-y7*y8)*y3**3+(y8*y9+ &
      y7*y8)*y4**3
    dF(2,188) = ((y9/2._ark+y7/2._ark)*y5+(sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y9/ &
      2._ark)*y6)*y1**3+((y7/2._ark-y9/2._ark)*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y9/ &
      2._ark)*y6)*y2**3+((y9/2._ark-y7/2._ark)*y5+(sqrt(3._ark)*y9/2._ark-sqrt(3._ark)*y7/ &
      2._ark)*y6)*y3**3+((-y9/2._ark-y7/2._ark)*y5+(-sqrt(3._ark)*y9/2._ark-sqrt(3._ark)*y7/ &
      2._ark)*y6)*y4**3
    dF(2,189) = y4*y7**2*y9**2+y1*y7**2*y9**2-y2*y7**2*y9**2-y3*y7**2*y9**2
    dF(2,190) = (y7**2*y8**2+y8**2*y9**2)*y1+(-y7**2*y8**2-y8**2*y9**2)*y2+(- &
      y7**2*y8**2-y8**2*y9**2)*y3+(y7**2*y8**2+y8**2*y9**2)*y4
    dF(2,191) = (y8**3*y9+y7*y8**3)*y1+(y8**3*y9-y7*y8**3)*y2+(y7*y8**3- &
      y8**3*y9)*y3+(-y7*y8**3-y8**3*y9)*y4
    dF(2,192) = y1*y8**4-y2*y8**4+y4*y8**4-y3*y8**4
    dF(2,193) = (y7*y9**3+y7**3*y9)*y1+(y7*y9**3+y7**3*y9)*y2+(y7*y9**3+ &
      y7**3*y9)*y3+(y7*y9**3+y7**3*y9)*y4
    dF(2,194) = y4*y7*y8**2*y9+y3*y7*y8**2*y9+y2*y7*y8**2*y9+y1*y7*y8**2*y9
    dF(2,195) = (y7*y8*y9**2+y7**2*y8*y9)*y1+(y7**2*y8*y9-y7*y8*y9**2)*y2+(- &
      y7**2*y8*y9+y7*y8*y9**2)*y3+(-y7*y8*y9**2-y7**2*y8*y9)*y4
    dF(2,196) = (sqrt(3._ark)*y5*y7**2*y8+(y7**2*y8+2._ark*y8*y9**2)*y6)*y1+ &
      (sqrt(3._ark)*y5*y7**2*y8+(y7**2*y8+2._ark*y8*y9**2)*y6)*y2+(sqrt(3._ark)*y5*y7**2*y8+ &
      (y7**2*y8+2._ark*y8*y9**2)*y6)*y3+(sqrt(3._ark)*y5*y7**2*y8+(y7**2*y8+ &
      2._ark*y8*y9**2)*y6)*y4
    dF(2,197) = (-sqrt(3._ark)*y5*y7*y8**2+(-y7*y8**2-2._ark*y8**2*y9)*y6)*y1+(- &
      sqrt(3._ark)*y5*y7*y8**2+(-y7*y8**2+2._ark*y8**2*y9)*y6)*y2+(sqrt(3._ark)*y5*y7*y8**2+ &
      (-2._ark*y8**2*y9+y7*y8**2)*y6)*y3+(sqrt(3._ark)*y5*y7*y8**2+(2._ark*y8**2*y9+ &
      y7*y8**2)*y6)*y4
    dF(2,198) = (y6**2*y8**2+y5**2*y8**2)*y1+(-y6**2*y8**2-y5**2*y8**2)*y2+(- &
      y6**2*y8**2-y5**2*y8**2)*y3+(y6**2*y8**2+y5**2*y8**2)*y4
    dF(2,199) = (y6**2*y7*y9+y5**2*y7*y9)*y1+(y6**2*y7*y9+y5**2*y7*y9)*y2+ &
      (y6**2*y7*y9+y5**2*y7*y9)*y3+(y6**2*y7*y9+y5**2*y7*y9)*y4
    dF(2,200) = (3._ark/2._ark*y5**2*y7*y8-sqrt(3._ark)*y5*y6*y8*y9+(y8*y9-y7*y8/ &
      2._ark)*y6**2)*y1+(-3._ark/2._ark*y5**2*y7*y8-sqrt(3._ark)*y5*y6*y8*y9+(y8*y9+y7*y8/ &
      2._ark)*y6**2)*y2+(3._ark/2._ark*y5**2*y7*y8+sqrt(3._ark)*y5*y6*y8*y9+(-y7*y8/2._ark- &
      y8*y9)*y6**2)*y3+(-3._ark/2._ark*y5**2*y7*y8+sqrt(3._ark)*y5*y6*y8*y9+(y7*y8/2._ark- &
      y8*y9)*y6**2)*y4
    dF(2,201) = (-4._ark/3._ark*sqrt(3._ark)*y5*y6**2*y8-y6**3*y8-y5**2*y6*y8)*y1+(-4._ark/ &
      3._ark*sqrt(3._ark)*y5*y6**2*y8-y6**3*y8-y5**2*y6*y8)*y2+(-4._ark/ &
      3._ark*sqrt(3._ark)*y5*y6**2*y8-y6**3*y8-y5**2*y6*y8)*y3+(-4._ark/ &
      3._ark*sqrt(3._ark)*y5*y6**2*y8-y6**3*y8-y5**2*y6*y8)*y4
    dF(2,202) = (-3._ark/2._ark*sqrt(3._ark)*y5**3*y9-27._ark/2._ark*y5**2*y6*y9+(- &
      6._ark*sqrt(3._ark)*y7-15._ark/2._ark*sqrt(3._ark)*y9)*y6**2*y5+(-7._ark/2._ark*y9- &
      10._ark*y7)*y6**3)*y1+(3._ark/2._ark*sqrt(3._ark)*y5**3*y9+27._ark/2._ark*y5**2*y6*y9+ &
      (15._ark/2._ark*sqrt(3._ark)*y9-6._ark*sqrt(3._ark)*y7)*y6**2*y5+(7._ark/2._ark*y9- &
      10._ark*y7)*y6**3)*y2+(-3._ark/2._ark*sqrt(3._ark)*y5**3*y9-27._ark/2._ark*y5**2*y6*y9+(- &
      15._ark/2._ark*sqrt(3._ark)*y9+6._ark*sqrt(3._ark)*y7)*y6**2*y5+(10._ark*y7-7._ark/ &
      2._ark*y9)*y6**3)*y3+(3._ark/2._ark*sqrt(3._ark)*y5**3*y9+27._ark/2._ark*y5**2*y6*y9+(15._ark/ &
      2._ark*sqrt(3._ark)*y9+6._ark*sqrt(3._ark)*y7)*y6**2*y5+(10._ark*y7+7._ark/ &
      2._ark*y9)*y6**3)*y4
    dF(2,203) = (y6**4+3._ark*y5**2*y6**2+2._ark*sqrt(3._ark)*y5*y6**3)*y1+(-y6**4- &
      3._ark*y5**2*y6**2-2._ark*sqrt(3._ark)*y5*y6**3)*y2+(-y6**4-3._ark*y5**2*y6**2- &
      2._ark*sqrt(3._ark)*y5*y6**3)*y3+(y6**4+3._ark*y5**2*y6**2+ &
      2._ark*sqrt(3._ark)*y5*y6**3)*y4
    dF(2,204) = (-y5*y8**3/2._ark-sqrt(3._ark)*y6*y8**3/2._ark)*y1+(-y5*y8**3/2._ark- &
      sqrt(3._ark)*y6*y8**3/2._ark)*y2+(-y5*y8**3/2._ark-sqrt(3._ark)*y6*y8**3/2._ark)*y3+(- &
      y5*y8**3/2._ark-sqrt(3._ark)*y6*y8**3/2._ark)*y4
    dF(2,205) = ((-2._ark*y7*y8**2+y8**2*y9)*y5-sqrt(3._ark)*y6*y8**2*y9)*y1+((- &
      2._ark*y7*y8**2-y8**2*y9)*y5+sqrt(3._ark)*y6*y8**2*y9)*y2+((2._ark*y7*y8**2+ &
      y8**2*y9)*y5-sqrt(3._ark)*y6*y8**2*y9)*y3+((2._ark*y7*y8**2-y8**2*y9)*y5+ &
      sqrt(3._ark)*y6*y8**2*y9)*y4
    dF(2,206) = ((y8*y9**2+y7**2*y8)*y5+(sqrt(3._ark)*y8*y9**2+ &
      sqrt(3._ark)*y7**2*y8)*y6)*y1+((y8*y9**2+y7**2*y8)*y5+(sqrt(3._ark)*y8*y9**2+ &
      sqrt(3._ark)*y7**2*y8)*y6)*y2+((y8*y9**2+y7**2*y8)*y5+(sqrt(3._ark)*y8*y9**2+ &
      sqrt(3._ark)*y7**2*y8)*y6)*y3+((y8*y9**2+y7**2*y8)*y5+(sqrt(3._ark)*y8*y9**2+ &
      sqrt(3._ark)*y7**2*y8)*y6)*y4
    dF(2,207) = ((y7*y9**2-y7**2*y9/2._ark)*y5+sqrt(3._ark)*y6*y7**2*y9/2._ark)*y1+ &
      ((y7**2*y9/2._ark+y7*y9**2)*y5-sqrt(3._ark)*y6*y7**2*y9/2._ark)*y2+((-y7*y9**2- &
      y7**2*y9/2._ark)*y5+sqrt(3._ark)*y6*y7**2*y9/2._ark)*y3+((y7**2*y9/2._ark-y7*y9**2)*y5- &
      sqrt(3._ark)*y6*y7**2*y9/2._ark)*y4
    dF(2,208) = ((-sqrt(3._ark)*y7**2/2._ark-sqrt(3._ark)*y9**2/2._ark)*y5**2+(y9**2+ &
      y7**2)*y6*y5+(-sqrt(3._ark)*y9**2/6._ark-sqrt(3._ark)*y7**2/6._ark)*y6**2)*y1+ &
      ((sqrt(3._ark)*y9**2/2._ark+sqrt(3._ark)*y7**2/2._ark)*y5**2+(-y9**2-y7**2)*y6*y5+ &
      (sqrt(3._ark)*y7**2/6._ark+sqrt(3._ark)*y9**2/6._ark)*y6**2)*y2+((sqrt(3._ark)*y9**2/2._ark+ &
      sqrt(3._ark)*y7**2/2._ark)*y5**2+(-y9**2-y7**2)*y6*y5+(sqrt(3._ark)*y7**2/6._ark+ &
      sqrt(3._ark)*y9**2/6._ark)*y6**2)*y3+((-sqrt(3._ark)*y7**2/2._ark-sqrt(3._ark)*y9**2/ &
      2._ark)*y5**2+(y9**2+y7**2)*y6*y5+(-sqrt(3._ark)*y9**2/6._ark-sqrt(3._ark)*y7**2/ &
      6._ark)*y6**2)*y4
    dF(2,209) = (-sqrt(3._ark)*y5**2*y7*y9/3._ark+y5*y6*y7*y9)*y1+(- &
      sqrt(3._ark)*y5**2*y7*y9/3._ark+y5*y6*y7*y9)*y2+(-sqrt(3._ark)*y5**2*y7*y9/3._ark+ &
      y5*y6*y7*y9)*y3+(-sqrt(3._ark)*y5**2*y7*y9/3._ark+y5*y6*y7*y9)*y4
    dF(2,210) = ((-sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y8/2._ark)*y5**2+(y7*y8- &
      y8*y9)*y6*y5+(sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y8/2._ark)*y6**2)*y1+((- &
      sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y8/2._ark)*y5**2+(-y7*y8-y8*y9)*y6*y5+ &
      (sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y8/2._ark)*y6**2)*y2+((sqrt(3._ark)*y8*y9/2._ark+ &
      sqrt(3._ark)*y7*y8/2._ark)*y5**2+(y8*y9+y7*y8)*y6*y5+(-sqrt(3._ark)*y8*y9/2._ark- &
      sqrt(3._ark)*y7*y8/2._ark)*y6**2)*y3+((sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y8/ &
      2._ark)*y5**2+(y8*y9-y7*y8)*y6*y5+(-sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y8/ &
      2._ark)*y6**2)*y4
    dF(2,211) = (-y5*y6*y8**2+sqrt(3._ark)*y5**2*y8**2/3._ark)*y1+(- &
      sqrt(3._ark)*y5**2*y8**2/3._ark+y5*y6*y8**2)*y2+(-sqrt(3._ark)*y5**2*y8**2/3._ark+ &
      y5*y6*y8**2)*y3+(-y5*y6*y8**2+sqrt(3._ark)*y5**2*y8**2/3._ark)*y4
    dF(2,212) = (-3._ark/2._ark*y5**3*y9-5._ark/2._ark*sqrt(3._ark)*y5**2*y6*y9+(-7._ark/ &
      2._ark*y9-2._ark*y7)*y6**2*y5+(-sqrt(3._ark)*y9/2._ark-2._ark*sqrt(3._ark)*y7)*y6**3)*y1+ &
      (3._ark/2._ark*y5**3*y9+5._ark/2._ark*sqrt(3._ark)*y5**2*y6*y9+(7._ark/2._ark*y9- &
      2._ark*y7)*y6**2*y5+(sqrt(3._ark)*y9/2._ark-2._ark*sqrt(3._ark)*y7)*y6**3)*y2+(-3._ark/ &
      2._ark*y5**3*y9-5._ark/2._ark*sqrt(3._ark)*y5**2*y6*y9+(2._ark*y7-7._ark/2._ark*y9)*y6**2*y5+ &
      (-sqrt(3._ark)*y9/2._ark+2._ark*sqrt(3._ark)*y7)*y6**3)*y3+(3._ark/2._ark*y5**3*y9+5._ark/ &
      2._ark*sqrt(3._ark)*y5**2*y6*y9+(7._ark/2._ark*y9+2._ark*y7)*y6**2*y5+(2._ark*sqrt(3._ark)*y7+ &
      sqrt(3._ark)*y9/2._ark)*y6**3)*y4
    dF(2,213) = ((y8*y9-y7*y8/2._ark)*y5**2+sqrt(3._ark)*y5*y6*y8*y9+3._ark/ &
      2._ark*y6**2*y7*y8)*y1+((y8*y9+y7*y8/2._ark)*y5**2+sqrt(3._ark)*y5*y6*y8*y9-3._ark/ &
      2._ark*y6**2*y7*y8)*y2+((-y7*y8/2._ark-y8*y9)*y5**2-sqrt(3._ark)*y5*y6*y8*y9+3._ark/ &
      2._ark*y6**2*y7*y8)*y3+((y7*y8/2._ark-y8*y9)*y5**2-sqrt(3._ark)*y5*y6*y8*y9-3._ark/ &
      2._ark*y6**2*y7*y8)*y4
    dF(2,214) = (-y5*y6**3+y5**3*y6-2._ark/3._ark*sqrt(3._ark)*y5**2*y6**2)*y1+(2._ark/ &
      3._ark*sqrt(3._ark)*y5**2*y6**2+y5*y6**3-y5**3*y6)*y2+(2._ark/ &
      3._ark*sqrt(3._ark)*y5**2*y6**2+y5*y6**3-y5**3*y6)*y3+(-y5*y6**3+y5**3*y6-2._ark/ &
      3._ark*sqrt(3._ark)*y5**2*y6**2)*y4
    dF(2,215) = (-y5**2*y6**2-2._ark*sqrt(3._ark)*y5*y6**3+y5**4)*y1+(-y5**4+ &
      2._ark*sqrt(3._ark)*y5*y6**3+y5**2*y6**2)*y2+(-y5**4+2._ark*sqrt(3._ark)*y5*y6**3+ &
      y5**2*y6**2)*y3+(-y5**2*y6**2-2._ark*sqrt(3._ark)*y5*y6**3+y5**4)*y4
    dF(2,216) = (y5**2*y8+y6**2*y8)*y4*y1+(y5**2*y8+y6**2*y8)*y3*y2
    dF(2,217) = (y5*y6*y8-sqrt(3._ark)*y5**2*y8/3._ark)*y4*y1+(y5*y6*y8- &
      sqrt(3._ark)*y5**2*y8/3._ark)*y3*y2
    dF(2,218) = ((y5**2*y8/4._ark+3._ark/4._ark*y6**2*y8-sqrt(3._ark)*y5*y6*y8/2._ark)*y2+ &
      y3*y5**2*y8)*y1+y2*y4*y5**2*y8+(y5**2*y8/4._ark+3._ark/4._ark*y6**2*y8- &
      sqrt(3._ark)*y5*y6*y8/2._ark)*y4*y3
    dF(2,219) = (y5**3-3._ark*y5*y6**2)*y4*y1+(-y5**3+3._ark*y5*y6**2)*y3*y2
    dF(2,220) = y1*y4**2*y8**2+y1**2*y4*y8**2-y2*y3**2*y8**2-y2**2*y3*y8**2
    dF(2,221) = (-y2*y9**2-y3*y7**2)*y1**2+(y3**2*y7**2+y2**2*y9**2)*y1- &
      y3*y4**2*y9**2+y2**2*y4*y7**2+y3**2*y4*y9**2-y2*y4**2*y7**2
    dF(2,222) = ((y6**2*y8/4._ark+3._ark/4._ark*y5**2*y8+sqrt(3._ark)*y5*y6*y8/2._ark)*y2+ &
      y3*y6**2*y8)*y1+y2*y4*y6**2*y8+(y6**2*y8/4._ark+3._ark/4._ark*y5**2*y8+ &
      sqrt(3._ark)*y5*y6*y8/2._ark)*y4*y3
    dF(2,223) = ((y6**2*y7+y5**2*y7)*y2+(y6**2*y9+y5**2*y9)*y3)*y1+(-y6**2*y9- &
      y5**2*y9)*y4*y2+(-y6**2*y7-y5**2*y7)*y4*y3
    dF(2,224) = (-4._ark/3._ark*sqrt(3._ark)*y5*y6**2-y5**2*y6-y6**3)*y4*y1+(y5**2*y6+ &
      y6**3+4._ark/3._ark*sqrt(3._ark)*y5*y6**2)*y3*y2
    dF(2,225) = (-2._ark*y2*y5*y7*y9+(-sqrt(3._ark)*y6*y7*y9+y5*y7*y9)*y3)*y1+(- &
      sqrt(3._ark)*y6*y7*y9+y5*y7*y9)*y4*y2-2._ark*y3*y4*y5*y7*y9
    dF(2,226) = (sqrt(3._ark)*y6*y8**2+y5*y8**2)*y4*y1+(-sqrt(3._ark)*y6*y8**2- &
      y5*y8**2)*y3*y2
    dF(2,227) = ((y5*y6*y7+sqrt(3._ark)*y6**2*y7/2._ark+sqrt(3._ark)*y5**2*y7/6._ark)*y2+ &
      (y5*y6*y9+sqrt(3._ark)*y6**2*y9/2._ark+sqrt(3._ark)*y5**2*y9/6._ark)*y3)*y1+(- &
      sqrt(3._ark)*y5**2*y9/6._ark-y5*y6*y9-sqrt(3._ark)*y6**2*y9/2._ark)*y4*y2+(-y5*y6*y7- &
      sqrt(3._ark)*y6**2*y7/2._ark-sqrt(3._ark)*y5**2*y7/6._ark)*y4*y3
    dF(2,228) = (((y8*y9+y7*y8)*y3+(y8*y9-y7*y8)*y4)*y2+(y7*y8-y8*y9)*y4*y3)*y1+(- &
      y7*y8-y8*y9)*y4*y3*y2
    dF(2,229) = ((y4*y8**2-y3*y8**2)*y2+y3*y4*y8**2)*y1-y2*y3*y4*y8**2
    dF(2,230) = (((y6**2+y5**2)*y3+(-y6**2-y5**2)*y4)*y2+(-y6**2-y5**2)*y4*y3)*y1+ &
      (y6**2+y5**2)*y4*y3*y2
    dF(2,231) = (y9+y7)*y3*y2*y1**2+((-y9+y7)*y4*y2**2+(-y7+y9)*y4*y3**2)*y1+(-y9- &
      y7)*y4**2*y3*y2
    dF(2,232) = ((-y6**2-y5**2)*y2+(-y6**2-y5**2)*y3)*y1**2+((y6**2+y5**2)*y2**2+ &
      (y6**2+y5**2)*y3**2)*y1+(y6**2+y5**2)*y4*y2**2+(-y6**2-y5**2)*y4**2*y2+(y6**2+ &
      y5**2)*y4*y3**2+(-y6**2-y5**2)*y4**2*y3
    dF(2,233) = (y2*y8+y3*y8)*y1**3+(y3**3*y8+y2**3*y8)*y1+y2*y4**3*y8+y2**3*y4*y8+ &
      y3*y4**3*y8+y3**3*y4*y8
    dF(2,234) = ((sqrt(3._ark)*y5*y7*y9-y6*y7*y9)*y2+(y6*y7*y9- &
      sqrt(3._ark)*y5*y7*y9)*y3)*y1+(y6*y7*y9-sqrt(3._ark)*y5*y7*y9)*y4*y2+ &
      (sqrt(3._ark)*y5*y7*y9-y6*y7*y9)*y4*y3
    dF(2,235) = (-sqrt(3._ark)*y6*y7*y9/2._ark-y5*y7*y9/2._ark)*y4*y1+(- &
      sqrt(3._ark)*y6*y7*y9/2._ark-y5*y7*y9/2._ark)*y3*y2
    dF(2,236) = (((-y9**2-y7**2)*y3+(y9**2+y7**2)*y4)*y2+(y9**2+y7**2)*y4*y3)*y1+(- &
      y9**2-y7**2)*y4*y3*y2
    dF(2,237) = ((y3*y7*y9+y4*y7*y9)*y2+y3*y4*y7*y9)*y1+y2*y3*y4*y7*y9
    dF(2,238) = ((((sqrt(3._ark)*y9/2._ark-sqrt(3._ark)*y7/2._ark)*y5+(y7/2._ark-y9/ &
      2._ark)*y6)*y3+((-sqrt(3._ark)*y9/2._ark-sqrt(3._ark)*y7/2._ark)*y5+(y9/2._ark+y7/ &
      2._ark)*y6)*y4)*y2+((sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y9/2._ark)*y5+(-y9/2._ark-y7/ &
      2._ark)*y6)*y4*y3)*y1+((sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y9/2._ark)*y5+(y9/2._ark-y7/ &
      2._ark)*y6)*y4*y3*y2
    dF(2,239) = (((y5*y8+sqrt(3._ark)*y6*y8)*y3+(y5*y8+sqrt(3._ark)*y6*y8)*y4)*y2+ &
      (y5*y8+sqrt(3._ark)*y6*y8)*y4*y3)*y1+(y5*y8+sqrt(3._ark)*y6*y8)*y4*y3*y2
    dF(2,240) = ((((y9/2._ark+y7/2._ark)*y5+(sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y9/ &
      2._ark)*y6)*y3+((y7/2._ark-y9/2._ark)*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y9/ &
      2._ark)*y6)*y4)*y2+((y9/2._ark-y7/2._ark)*y5+(sqrt(3._ark)*y9/2._ark-sqrt(3._ark)*y7/ &
      2._ark)*y6)*y4*y3)*y1+((-y9/2._ark-y7/2._ark)*y5+(-sqrt(3._ark)*y9/2._ark-sqrt(3._ark)*y7/ &
      2._ark)*y6)*y4*y3*y2
    dF(2,241) = (((-sqrt(3._ark)*y5**2/2._ark+y5*y6-sqrt(3._ark)*y6**2/6._ark)*y3+ &
      (sqrt(3._ark)*y5**2/2._ark-y5*y6+sqrt(3._ark)*y6**2/6._ark)*y4)*y2+(sqrt(3._ark)*y5**2/ &
      2._ark-y5*y6+sqrt(3._ark)*y6**2/6._ark)*y4*y3)*y1+(-sqrt(3._ark)*y5**2/2._ark+y5*y6- &
      sqrt(3._ark)*y6**2/6._ark)*y4*y3*y2
    dF(2,242) = y1*y2*y3*y4*y8
    dF(2,243) = (-y3*y4*y9-y2*y4*y7)*y1**2+(-y2**2*y3*y7+(-y3**2*y9+y4**2*y9)*y2+ &
      y3*y4**2*y7)*y1+y2**2*y3*y4*y9+y2*y3**2*y4*y7
    dF(2,244) = y2**2*y8**3+y1**2*y8**3+y3**2*y8**3+y4**2*y8**3
    dF(2,245) = (y7**3+y9**3)*y1**2+(y7**3-y9**3)*y2**2+(y9**3-y7**3)*y3**2+(-y7**3- &
      y9**3)*y4**2
    dF(2,246) = y4**2*y7*y8*y9+y1**2*y7*y8*y9-y3**2*y7*y8*y9-y2**2*y7*y8*y9
    dF(2,247) = ((sqrt(3._ark)*y9**2-sqrt(3._ark)*y7**2)*y5+(y7**2-y9**2)*y6)*y1**2+ &
      ((sqrt(3._ark)*y7**2-sqrt(3._ark)*y9**2)*y5+(-y7**2+y9**2)*y6)*y2**2+ &
      ((sqrt(3._ark)*y7**2-sqrt(3._ark)*y9**2)*y5+(-y7**2+y9**2)*y6)*y3**2+ &
      ((sqrt(3._ark)*y9**2-sqrt(3._ark)*y7**2)*y5+(y7**2-y9**2)*y6)*y4**2
    dF(2,248) = (-sqrt(3._ark)*y5*y8**2/3._ark-y6*y8**2)*y1**2+(sqrt(3._ark)*y5*y8**2/ &
      3._ark+y6*y8**2)*y2**2+(sqrt(3._ark)*y5*y8**2/3._ark+y6*y8**2)*y3**2+(- &
      sqrt(3._ark)*y5*y8**2/3._ark-y6*y8**2)*y4**2
    dF(2,249) = (sqrt(3._ark)*y5**3-5._ark/3._ark*sqrt(3._ark)*y5*y6**2+y5**2*y6+ &
      y6**3)*y1**2+(-sqrt(3._ark)*y5**3-y6**3+5._ark/3._ark*sqrt(3._ark)*y5*y6**2- &
      y5**2*y6)*y2**2+(-sqrt(3._ark)*y5**3-y6**3+5._ark/3._ark*sqrt(3._ark)*y5*y6**2- &
      y5**2*y6)*y3**2+(sqrt(3._ark)*y5**3-5._ark/3._ark*sqrt(3._ark)*y5*y6**2+y5**2*y6+ &
      y6**3)*y4**2
    dF(2,250) = (-sqrt(3._ark)*y6*y7*y9/2._ark-y5*y7*y9/2._ark)*y1**2+(- &
      sqrt(3._ark)*y6*y7*y9/2._ark-y5*y7*y9/2._ark)*y2**2+(-sqrt(3._ark)*y6*y7*y9/2._ark- &
      y5*y7*y9/2._ark)*y3**2+(-sqrt(3._ark)*y6*y7*y9/2._ark-y5*y7*y9/2._ark)*y4**2
    dF(2,251) = ((sqrt(3._ark)*y7/6._ark+sqrt(3._ark)*y9/6._ark)*y5**2+(y9+y7)*y6*y5+ &
      (sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y9/2._ark)*y6**2)*y1**2+((-sqrt(3._ark)*y9/6._ark+ &
      sqrt(3._ark)*y7/6._ark)*y5**2+(-y9+y7)*y6*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y9/ &
      2._ark)*y6**2)*y2**2+((sqrt(3._ark)*y9/6._ark-sqrt(3._ark)*y7/6._ark)*y5**2+(-y7+ &
      y9)*y6*y5+(sqrt(3._ark)*y9/2._ark-sqrt(3._ark)*y7/2._ark)*y6**2)*y3**2+((-sqrt(3._ark)*y7/ &
      6._ark-sqrt(3._ark)*y9/6._ark)*y5**2+(-y9-y7)*y6*y5+(-sqrt(3._ark)*y9/2._ark- &
      sqrt(3._ark)*y7/2._ark)*y6**2)*y4**2
    dF(2,252) = (-sqrt(3._ark)*y6**2*y8/3._ark-y5*y6*y8)*y1**2+(-sqrt(3._ark)*y6**2*y8/ &
      3._ark-y5*y6*y8)*y2**2+(-sqrt(3._ark)*y6**2*y8/3._ark-y5*y6*y8)*y3**2+(- &
      sqrt(3._ark)*y6**2*y8/3._ark-y5*y6*y8)*y4**2
    dF(2,253) = ((-sqrt(3._ark)*y9/6._ark+sqrt(3._ark)*y7/3._ark)*y5**2+y5*y6*y7+ &
      sqrt(3._ark)*y6**2*y9/2._ark)*y1**2+((sqrt(3._ark)*y7/3._ark+sqrt(3._ark)*y9/6._ark)*y5**2+ &
      y5*y6*y7-sqrt(3._ark)*y6**2*y9/2._ark)*y2**2+((-sqrt(3._ark)*y7/3._ark-sqrt(3._ark)*y9/ &
      6._ark)*y5**2-y5*y6*y7+sqrt(3._ark)*y6**2*y9/2._ark)*y3**2+((sqrt(3._ark)*y9/6._ark- &
      sqrt(3._ark)*y7/3._ark)*y5**2-y5*y6*y7-sqrt(3._ark)*y6**2*y9/2._ark)*y4**2
    dF(2,254) = (y5**2*y8+y6**2*y8)*y1**2+(y5**2*y8+y6**2*y8)*y2**2+(y5**2*y8+ &
      y6**2*y8)*y3**2+(y5**2*y8+y6**2*y8)*y4**2
    dF(2,255) = (y5**3-3._ark*y5*y6**2)*y1**2+(-y5**3+3._ark*y5*y6**2)*y2**2+(-y5**3+ &
      3._ark*y5*y6**2)*y3**2+(y5**3-3._ark*y5*y6**2)*y4**2
    dF(2,256) = y2**2*y3**2*y8+y1**2*y4**2*y8
    dF(2,257) = (y3**2*y8+y2**2*y8)*y1**2+y3**2*y4**2*y8+y2**2*y4**2*y8
    dF(2,258) = (y3**2*y9+y2**2*y7)*y1**2-y3**2*y4**2*y7-y2**2*y4**2*y9
    dF(2,259) = (y6+sqrt(3._ark)*y5/3._ark)*y4**2*y1**2+(-sqrt(3._ark)*y5/3._ark- &
      y6)*y3**2*y2**2
    dF(2,260) = ((sqrt(3._ark)*y7-sqrt(3._ark)*y9)*y5+(-y7+y9)*y6)*y4*y1**2+((- &
      sqrt(3._ark)*y7+sqrt(3._ark)*y9)*y5+(-y9+y7)*y6)*y4**2*y1+((sqrt(3._ark)*y7+ &
      sqrt(3._ark)*y9)*y5+(-y9-y7)*y6)*y3*y2**2+((-sqrt(3._ark)*y9-sqrt(3._ark)*y7)*y5+(y9+ &
      y7)*y6)*y3**2*y2
    dF(2,261) = ((sqrt(3._ark)*y5*y8/2._ark+y6*y8/2._ark)*y2+y3*y6*y8)*y1**2+ &
      ((sqrt(3._ark)*y5*y8/2._ark+y6*y8/2._ark)*y2**2+y3**2*y6*y8)*y1+y2**2*y4*y6*y8+ &
      y2*y4**2*y6*y8+(sqrt(3._ark)*y5*y8/2._ark+y6*y8/2._ark)*y4*y3**2+(sqrt(3._ark)*y5*y8/ &
      2._ark+y6*y8/2._ark)*y4**2*y3
    dF(2,262) = (-sqrt(3._ark)*y6**2/3._ark-y5*y6)*y4*y1**2+(-sqrt(3._ark)*y6**2/3._ark- &
      y5*y6)*y4**2*y1+(y5*y6+sqrt(3._ark)*y6**2/3._ark)*y3*y2**2+(y5*y6+sqrt(3._ark)*y6**2/ &
      3._ark)*y3**2*y2
    dF(2,263) = y1**2*y2*y3*y8+(y3**2*y4*y8+y2**2*y4*y8)*y1+y2*y3*y4**2*y8
    dF(2,264) = (y3*y4**2+y2*y4**2)*y1**2-y1*y2**2*y3**2-y2**2*y3**2*y4
    dF(2,265) = ((-y4+y3)*y2**2+y2*y3**2-y3**2*y4)*y1**2+(-y2**2*y4**2- &
      y3**2*y4**2)*y1+y2**2*y3*y4**2+y2*y3**2*y4**2
    dF(2,266) = ((y5*y9-sqrt(3._ark)*y6*y9)*y2-2._ark*y3*y5*y7)*y1**2+ &
      ((sqrt(3._ark)*y6*y9-y5*y9)*y2**2+2._ark*y3**2*y5*y7)*y1-2._ark*y2**2*y4*y5*y7+ &
      2._ark*y2*y4**2*y5*y7+(y5*y9-sqrt(3._ark)*y6*y9)*y4*y3**2+(sqrt(3._ark)*y6*y9- &
      y5*y9)*y4**2*y3
    dF(2,267) = (-y5*y8/2._ark-sqrt(3._ark)*y6*y8/2._ark)*y4*y1**2+(-y5*y8/2._ark- &
      sqrt(3._ark)*y6*y8/2._ark)*y4**2*y1+(-y5*y8/2._ark-sqrt(3._ark)*y6*y8/2._ark)*y3*y2**2+(- &
      y5*y8/2._ark-sqrt(3._ark)*y6*y8/2._ark)*y3**2*y2
    dF(2,268) = ((-sqrt(3._ark)*y6*y7+y5*y7)*y2-2._ark*y3*y5*y9)*y1**2+((- &
      sqrt(3._ark)*y6*y7+y5*y7)*y2**2-2._ark*y3**2*y5*y9)*y1+2._ark*y2**2*y4*y5*y9+ &
      2._ark*y2*y4**2*y5*y9+(-y5*y7+sqrt(3._ark)*y6*y7)*y4*y3**2+(-y5*y7+ &
      sqrt(3._ark)*y6*y7)*y4**2*y3
    dF(2,269) = ((y9-2._ark*y7)*y5-sqrt(3._ark)*y6*y9)*y4*y1**2+((-y9+2._ark*y7)*y5+ &
      sqrt(3._ark)*y6*y9)*y4**2*y1+((-y9-2._ark*y7)*y5+sqrt(3._ark)*y6*y9)*y3*y2**2+((y9+ &
      2._ark*y7)*y5-sqrt(3._ark)*y6*y9)*y3**2*y2
    dF(2,270) = (-2._ark*y2*y4*y5+(y5-sqrt(3._ark)*y6)*y4*y3)*y1**2+(2._ark*y2**2*y3*y5+ &
      ((sqrt(3._ark)*y6-y5)*y3**2+(y5-sqrt(3._ark)*y6)*y4**2)*y2-2._ark*y3*y4**2*y5)*y1+ &
      (sqrt(3._ark)*y6-y5)*y4*y3*y2**2+2._ark*y2*y3**2*y4*y5
    dF(2,271) = y1**2*y2*y3*y4+(-y2**2*y3*y4+(y3*y4**2-y3**2*y4)*y2)*y1
    dF(2,272) = (-y6*y8-sqrt(3._ark)*y5*y8/3._ark)*y1**3+(-y6*y8-sqrt(3._ark)*y5*y8/ &
      3._ark)*y2**3+(-y6*y8-sqrt(3._ark)*y5*y8/3._ark)*y3**3+(-y6*y8-sqrt(3._ark)*y5*y8/ &
      3._ark)*y4**3
    dF(2,273) = (y5*y6-sqrt(3._ark)*y5**2/3._ark)*y1**3+(sqrt(3._ark)*y5**2/3._ark- &
      y5*y6)*y2**3+(sqrt(3._ark)*y5**2/3._ark-y5*y6)*y3**3+(y5*y6-sqrt(3._ark)*y5**2/ &
      3._ark)*y4**3
    dF(2,274) = y1**3*y4**2+y1**2*y4**3-y2**2*y3**3-y2**3*y3**2
    dF(2,275) = (y3**2+y2**2)*y1**3+(-y3**3-y2**3)*y1**2-y2**3*y4**2+y2**2*y4**3+ &
      y3**2*y4**3-y3**3*y4**2
    dF(2,276) = (-y2*y6+(-y6/2._ark-sqrt(3._ark)*y5/2._ark)*y3)*y1**3+(y2**3*y6+ &
      (sqrt(3._ark)*y5/2._ark+y6/2._ark)*y3**3)*y1+(sqrt(3._ark)*y5/2._ark+y6/2._ark)*y4*y2**3+(- &
      y6/2._ark-sqrt(3._ark)*y5/2._ark)*y4**3*y2-y3*y4**3*y6+y3**3*y4*y6
    dF(2,277) = y1**4*y8+y4**4*y8+y3**4*y8+y2**4*y8
    dF(2,278) = (y9+y7)*y1**4+(-y9+y7)*y2**4+(-y7+y9)*y3**4+(-y9-y7)*y4**4
    dF(2,279) = (-sqrt(3._ark)*y6/2._ark-y5/2._ark)*y1**4+(sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y2**4+(sqrt(3._ark)*y6/2._ark+y5/2._ark)*y3**4+(-sqrt(3._ark)*y6/2._ark-y5/ &
      2._ark)*y4**4
    dF(2,280) = y1*y4**4+y1**4*y4-y2*y3**4-y2**4*y3
    dF(2,281) = y1**5-y3**5-y2**5+y4**5
    dF(2,282) = y7*y9**5+y7**5*y9
    dF(2,283) = y7**3*y9**3
    dF(2,284) = y7**3*y8**2*y9+y7*y8**2*y9**3
    dF(2,285) = y7*y8**4*y9
    dF(2,286) = -sqrt(3._ark)*y5*y7**2*y8**3+(-y7**2*y8**3-2._ark*y8**3*y9**2)*y6
    dF(2,287) = (sqrt(3._ark)*y8*y9**4-sqrt(3._ark)*y7**4*y8)*y5+(y7**4*y8-y8*y9**4)*y6
    dF(2,288) = -4._ark/3._ark*sqrt(3._ark)*y5*y6**2*y8**3-y5**2*y6*y8**3-y6**3*y8**3
    dF(2,289) = -sqrt(3._ark)*y6*y7**2*y8*y9**2/2._ark-y5*y7**2*y8*y9**2/2._ark
    dF(2,290) = (-2._ark*y7**2*y8**3+y8**3*y9**2)*y5-sqrt(3._ark)*y6*y8**3*y9**2
    dF(2,291) = (y8*y9**4-2._ark*y7**4*y8)*y5-sqrt(3._ark)*y6*y8*y9**4
    dF(2,292) = y5*y8**5+sqrt(3._ark)*y6*y8**5
    dF(2,293) = (-sqrt(3._ark)*y7*y9**3/2._ark-sqrt(3._ark)*y7**3*y9/2._ark)*y5**2+ &
      (y7*y9**3+y7**3*y9)*y6*y5+(-sqrt(3._ark)*y7*y9**3/6._ark-sqrt(3._ark)*y7**3*y9/ &
      6._ark)*y6**2
    dF(2,294) = -sqrt(3._ark)*y6**2*y7*y8**2*y9/3._ark-y5*y6*y7*y8**2*y9
    dF(2,295) = -sqrt(3._ark)*y5**2*y7*y9**3/2._ark+y5*y6*y7**3*y9+(- &
      sqrt(3._ark)*y7**3*y9/3._ark+sqrt(3._ark)*y7*y9**3/6._ark)*y6**2
    dF(2,296) = (-y8*y9**2/3._ark-4._ark/3._ark*y7**2*y8)*y5**3- &
      sqrt(3._ark)*y5**2*y6*y8*y9**2+y5*y6**2*y8*y9**2+(-5._ark/9._ark*sqrt(3._ark)*y8*y9**2- &
      4._ark/9._ark*sqrt(3._ark)*y7**2*y8)*y6**3
    dF(2,297) = y5**3*y7**2*y8+sqrt(3._ark)*y5**2*y6*y8*y9**2+y5*y6**2*y7**2*y8+(5._ark/ &
      9._ark*sqrt(3._ark)*y8*y9**2+4._ark/9._ark*sqrt(3._ark)*y7**2*y8)*y6**3
    dF(2,298) = y6**2*y7*y8**2*y9+y5**2*y7*y8**2*y9
    dF(2,299) = (y7*y9**3+y7**3*y9)*y5**2+(y7*y9**3+y7**3*y9)*y6**2
    dF(2,300) = (4._ark/9._ark*sqrt(3._ark)*y8*y9**2+4._ark/9._ark*sqrt(3._ark)*y7**2*y8)*y5**3+ &
      2._ark*y5**2*y6*y8*y9**2+(2._ark/3._ark*y8*y9**2+4._ark/3._ark*y7**2*y8)*y6**3
    dF(2,301) = (4._ark/9._ark*sqrt(3._ark)*y8*y9**2+4._ark/9._ark*sqrt(3._ark)*y7**2*y8)*y5**3+ &
      (y8*y9**2+y7**2*y8)*y6*y5**2+(y8*y9**2+y7**2*y8)*y6**3
    dF(2,302) = -sqrt(3._ark)*y5*y6**3*y7*y9/3._ark+y5**2*y6**2*y7*y9/4._ark+5._ark/ &
      24._ark*y6**4*y7*y9+3._ark/8._ark*y5**4*y7*y9
    dF(2,303) = 17._ark/6._ark*sqrt(3._ark)*y5*y6**4*y8-2._ark*y6**5*y8-3._ark*y5**4*y6*y8- &
      3._ark/2._ark*sqrt(3._ark)*y5**5*y8-y5**2*y6**3*y8
    dF(2,304) = -3._ark*y5*y6**2*y8**3+y5**3*y8**3
    dF(2,305) = y5**3*y6*y7*y9-sqrt(3._ark)*y6**4*y7*y9/36._ark+y5*y6**3*y7*y9/3._ark- &
      sqrt(3._ark)*y5**2*y6**2*y7*y9/2._ark-sqrt(3._ark)*y5**4*y7*y9/4._ark
    dF(2,306) = -6._ark*y5*y6**4*y8+5._ark/2._ark*sqrt(3._ark)*y5**4*y6*y8+y5**3*y6**2*y8+ &
      3._ark/2._ark*sqrt(3._ark)*y6**5*y8+3._ark*y5**5*y8
    dF(2,307) = 7._ark/12._ark*y6**4*y7*y9+y5**4*y7*y9/4._ark+2._ark/ &
      3._ark*sqrt(3._ark)*y5*y6**3*y7*y9+3._ark/2._ark*y5**2*y6**2*y7*y9
    dF(2,308) = -15._ark*y5*y6**4*y8+5._ark*sqrt(3._ark)*y5**4*y6*y8+ &
      3._ark*sqrt(3._ark)*y6**5*y8+7._ark*y5**5*y8
    dF(2,309) = (-y8**4*y9-y7*y8**4)*y1+(y8**4*y9-y7*y8**4)*y2+(y7*y8**4- &
      y8**4*y9)*y3+(y7*y8**4+y8**4*y9)*y4
    dF(2,310) = ((-sqrt(3._ark)*y8**2*y9**2+sqrt(3._ark)*y7**2*y8**2)*y5+(y8**2*y9**2- &
      y7**2*y8**2)*y6)*y1+((sqrt(3._ark)*y8**2*y9**2-sqrt(3._ark)*y7**2*y8**2)*y5+ &
      (y7**2*y8**2-y8**2*y9**2)*y6)*y2+((sqrt(3._ark)*y8**2*y9**2- &
      sqrt(3._ark)*y7**2*y8**2)*y5+(y7**2*y8**2-y8**2*y9**2)*y6)*y3+((- &
      sqrt(3._ark)*y8**2*y9**2+sqrt(3._ark)*y7**2*y8**2)*y5+(y8**2*y9**2- &
      y7**2*y8**2)*y6)*y4
    dF(2,311) = (sqrt(3._ark)*y5*y7*y9**3/2._ark+(y7**3*y9+y7*y9**3/2._ark)*y6)*y1+ &
      (sqrt(3._ark)*y5*y7*y9**3/2._ark+(y7**3*y9+y7*y9**3/2._ark)*y6)*y2+ &
      (sqrt(3._ark)*y5*y7*y9**3/2._ark+(y7**3*y9+y7*y9**3/2._ark)*y6)*y3+ &
      (sqrt(3._ark)*y5*y7*y9**3/2._ark+(y7**3*y9+y7*y9**3/2._ark)*y6)*y4
    dF(2,312) = ((-y8**3*y9/2._ark-y7*y8**3/2._ark)*y5+(-sqrt(3._ark)*y8**3*y9/2._ark- &
      sqrt(3._ark)*y7*y8**3/2._ark)*y6)*y1+((-y8**3*y9/2._ark+y7*y8**3/2._ark)*y5+(- &
      sqrt(3._ark)*y8**3*y9/2._ark+sqrt(3._ark)*y7*y8**3/2._ark)*y6)*y2+((y8**3*y9/2._ark- &
      y7*y8**3/2._ark)*y5+(sqrt(3._ark)*y8**3*y9/2._ark-sqrt(3._ark)*y7*y8**3/2._ark)*y6)*y3+ &
      ((y8**3*y9/2._ark+y7*y8**3/2._ark)*y5+(sqrt(3._ark)*y8**3*y9/2._ark+sqrt(3._ark)*y7*y8**3/ &
      2._ark)*y6)*y4
    dF(2,313) = ((y7**3*y9-y7*y9**3/2._ark)*y5+sqrt(3._ark)*y6*y7*y9**3/2._ark)*y1+ &
      ((y7**3*y9-y7*y9**3/2._ark)*y5+sqrt(3._ark)*y6*y7*y9**3/2._ark)*y2+((y7**3*y9- &
      y7*y9**3/2._ark)*y5+sqrt(3._ark)*y6*y7*y9**3/2._ark)*y3+((y7**3*y9-y7*y9**3/2._ark)*y5+ &
      sqrt(3._ark)*y6*y7*y9**3/2._ark)*y4
    dF(2,314) = ((-sqrt(3._ark)*y7*y9**2/6._ark-2._ark/3._ark*sqrt(3._ark)*y7**2*y9)*y5**2+ &
      y5*y6*y7*y9**2-sqrt(3._ark)*y6**2*y7*y9**2/2._ark)*y1+((2._ark/ &
      3._ark*sqrt(3._ark)*y7**2*y9-sqrt(3._ark)*y7*y9**2/6._ark)*y5**2+y5*y6*y7*y9**2- &
      sqrt(3._ark)*y6**2*y7*y9**2/2._ark)*y2+((sqrt(3._ark)*y7*y9**2/6._ark-2._ark/ &
      3._ark*sqrt(3._ark)*y7**2*y9)*y5**2-y5*y6*y7*y9**2+sqrt(3._ark)*y6**2*y7*y9**2/ &
      2._ark)*y3+((sqrt(3._ark)*y7*y9**2/6._ark+2._ark/3._ark*sqrt(3._ark)*y7**2*y9)*y5**2- &
      y5*y6*y7*y9**2+sqrt(3._ark)*y6**2*y7*y9**2/2._ark)*y4
    dF(2,315) = ((-y7*y8**2-y8**2*y9)*y5**2+(-y7*y8**2-y8**2*y9)*y6**2)*y1+((- &
      y7*y8**2+y8**2*y9)*y5**2+(-y7*y8**2+y8**2*y9)*y6**2)*y2+((-y8**2*y9+ &
      y7*y8**2)*y5**2+(-y8**2*y9+y7*y8**2)*y6**2)*y3+((y8**2*y9+y7*y8**2)*y5**2+ &
      (y8**2*y9+y7*y8**2)*y6**2)*y4
    dF(2,316) = ((y8*y9**2+y7**2*y8)*y5**2+(y8*y9**2+y7**2*y8)*y6**2)*y1+((y8*y9**2+ &
      y7**2*y8)*y5**2+(y8*y9**2+y7**2*y8)*y6**2)*y2+((y8*y9**2+y7**2*y8)*y5**2+ &
      (y8*y9**2+y7**2*y8)*y6**2)*y3+((y8*y9**2+y7**2*y8)*y5**2+(y8*y9**2+ &
      y7**2*y8)*y6**2)*y4
    dF(2,317) = (6._ark*y5**4*y7-4._ark*sqrt(3._ark)*y5**3*y6*y9+(10._ark*y7+ &
      10._ark*y9)*y6**2*y5**2+4._ark*sqrt(3._ark)*y5*y6**3*y7+6._ark*y6**4*y9)*y1+ &
      (6._ark*y5**4*y7+4._ark*sqrt(3._ark)*y5**3*y6*y9+(-10._ark*y9+10._ark*y7)*y6**2*y5**2+ &
      4._ark*sqrt(3._ark)*y5*y6**3*y7-6._ark*y6**4*y9)*y2+(-6._ark*y5**4*y7- &
      4._ark*sqrt(3._ark)*y5**3*y6*y9+(-10._ark*y7+10._ark*y9)*y6**2*y5**2- &
      4._ark*sqrt(3._ark)*y5*y6**3*y7+6._ark*y6**4*y9)*y3+(-6._ark*y5**4*y7+ &
      4._ark*sqrt(3._ark)*y5**3*y6*y9+(-10._ark*y9-10._ark*y7)*y6**2*y5**2- &
      4._ark*sqrt(3._ark)*y5*y6**3*y7-6._ark*y6**4*y9)*y4
    dF(2,318) = y4**2*y7*y8**2*y9+y3**2*y7*y8**2*y9+y1**2*y7*y8**2*y9+ &
      y2**2*y7*y8**2*y9
    dF(2,319) = ((sqrt(3._ark)*y7**2*y9-sqrt(3._ark)*y7*y9**2)*y5+(-y7**2*y9+ &
      y7*y9**2)*y6)*y1**2+((-sqrt(3._ark)*y7**2*y9-sqrt(3._ark)*y7*y9**2)*y5+(y7*y9**2+ &
      y7**2*y9)*y6)*y2**2+((sqrt(3._ark)*y7*y9**2+sqrt(3._ark)*y7**2*y9)*y5+(-y7**2*y9- &
      y7*y9**2)*y6)*y3**2+((sqrt(3._ark)*y7*y9**2-sqrt(3._ark)*y7**2*y9)*y5+(y7**2*y9- &
      y7*y9**2)*y6)*y4**2
    dF(2,320) = (y5**2*y6*y8+y6**3*y8+4._ark/3._ark*sqrt(3._ark)*y5*y6**2*y8)*y1**2+ &
      (y5**2*y6*y8+y6**3*y8+4._ark/3._ark*sqrt(3._ark)*y5*y6**2*y8)*y2**2+(y5**2*y6*y8+ &
      y6**3*y8+4._ark/3._ark*sqrt(3._ark)*y5*y6**2*y8)*y3**2+(y5**2*y6*y8+y6**3*y8+4._ark/ &
      3._ark*sqrt(3._ark)*y5*y6**2*y8)*y4**2
    dF(2,321) = (y5**3*y8-3._ark*y5*y6**2*y8)*y1**2+(y5**3*y8-3._ark*y5*y6**2*y8)*y2**2+ &
      (y5**3*y8-3._ark*y5*y6**2*y8)*y3**2+(y5**3*y8-3._ark*y5*y6**2*y8)*y4**2
    dF(2,322) = (y5**2*y8+y6**2*y8)*y1**3+(y5**2*y8+y6**2*y8)*y2**3+(y5**2*y8+ &
      y6**2*y8)*y3**3+(y5**2*y8+y6**2*y8)*y4**3
    dF(2,323) = ((2._ark*y7+2._ark*y9)*y5**2+(-2._ark*sqrt(3._ark)*y9- &
      2._ark*sqrt(3._ark)*y7)*y6*y5)*y1**3+((2._ark*y7-2._ark*y9)*y5**2+(2._ark*sqrt(3._ark)*y9- &
      2._ark*sqrt(3._ark)*y7)*y6*y5)*y2**3+((2._ark*y9-2._ark*y7)*y5**2+(2._ark*sqrt(3._ark)*y7- &
      2._ark*sqrt(3._ark)*y9)*y6*y5)*y3**3+((-2._ark*y9-2._ark*y7)*y5**2+(2._ark*sqrt(3._ark)*y7+ &
      2._ark*sqrt(3._ark)*y9)*y6*y5)*y4**3
    dF(2,324) = ((-y7+y9/2._ark)*y5-sqrt(3._ark)*y6*y9/2._ark)*y1**4+((-y7-y9/2._ark)*y5+ &
      sqrt(3._ark)*y6*y9/2._ark)*y2**4+((y9/2._ark+y7)*y5-sqrt(3._ark)*y6*y9/2._ark)*y3**4+((y7- &
      y9/2._ark)*y5+sqrt(3._ark)*y6*y9/2._ark)*y4**4
    dF(2,325) = (-sqrt(3._ark)*y5**2/2._ark+y5*y6-sqrt(3._ark)*y6**2/6._ark)*y1**4+ &
      (sqrt(3._ark)*y5**2/2._ark-y5*y6+sqrt(3._ark)*y6**2/6._ark)*y2**4+(sqrt(3._ark)*y5**2/ &
      2._ark-y5*y6+sqrt(3._ark)*y6**2/6._ark)*y3**4+(-sqrt(3._ark)*y5**2/2._ark+y5*y6- &
      sqrt(3._ark)*y6**2/6._ark)*y4**4
    dF(2,326) = (-y6**2-y5**2)*y1**4+(y6**2+y5**2)*y2**4+(y6**2+y5**2)*y3**4+(- &
      y6**2-y5**2)*y4**4
    dF(2,327) = y2**5*y8+y4**5*y8+y3**5*y8+y1**5*y8
    dF(2,328) = (-y9-y7)*y1**5+(-y7+y9)*y2**5+(-y9+y7)*y3**5+(y9+y7)*y4**5
    dF(2,329) = (-y7**5-y9**5)*y1+(y9**5-y7**5)*y2+(y7**5-y9**5)*y3+(y7**5+y9**5)*y4
    dF(2,330) = (-y7**3*y8**2-y8**2*y9**3)*y1+(-y7**3*y8**2+y8**2*y9**3)*y2+ &
      (y7**3*y8**2-y8**2*y9**3)*y3+(y8**2*y9**3+y7**3*y8**2)*y4
    dF(2,331) = (y7*y8**2*y9**2+y7**2*y8**2*y9)*y1+(y7*y8**2*y9**2- &
      y7**2*y8**2*y9)*y2+(-y7*y8**2*y9**2+y7**2*y8**2*y9)*y3+(-y7**2*y8**2*y9- &
      y7*y8**2*y9**2)*y4
    dF(2,332) = (y6*y7*y8**2*y9+sqrt(3._ark)*y5*y7*y8**2*y9/3._ark)*y1+(y6*y7*y8**2*y9+ &
      sqrt(3._ark)*y5*y7*y8**2*y9/3._ark)*y2+(y6*y7*y8**2*y9+sqrt(3._ark)*y5*y7*y8**2*y9/ &
      3._ark)*y3+(y6*y7*y8**2*y9+sqrt(3._ark)*y5*y7*y8**2*y9/3._ark)*y4
    dF(2,333) = (-3._ark*sqrt(3._ark)*y5**3*y7*y8-9._ark*y5**2*y6*y7*y8- &
      3._ark*sqrt(3._ark)*y5*y6**2*y7*y8+(-y7*y8-8._ark*y8*y9)*y6**3)*y1+ &
      (3._ark*sqrt(3._ark)*y5**3*y7*y8+9._ark*y5**2*y6*y7*y8+3._ark*sqrt(3._ark)*y5*y6**2*y7*y8+ &
      (y7*y8-8._ark*y8*y9)*y6**3)*y2+(-3._ark*sqrt(3._ark)*y5**3*y7*y8-9._ark*y5**2*y6*y7*y8- &
      3._ark*sqrt(3._ark)*y5*y6**2*y7*y8+(-y7*y8+8._ark*y8*y9)*y6**3)*y3+ &
      (3._ark*sqrt(3._ark)*y5**3*y7*y8+9._ark*y5**2*y6*y7*y8+3._ark*sqrt(3._ark)*y5*y6**2*y7*y8+ &
      (8._ark*y8*y9+y7*y8)*y6**3)*y4
    dF(2,334) = ((y9**4/2._ark+y7**4/2._ark)*y5+(sqrt(3._ark)*y9**4/2._ark+sqrt(3._ark)*y7**4/ &
      2._ark)*y6)*y1+((-y9**4/2._ark-y7**4/2._ark)*y5+(-sqrt(3._ark)*y7**4/2._ark- &
      sqrt(3._ark)*y9**4/2._ark)*y6)*y2+((-y9**4/2._ark-y7**4/2._ark)*y5+(-sqrt(3._ark)*y7**4/ &
      2._ark-sqrt(3._ark)*y9**4/2._ark)*y6)*y3+((y9**4/2._ark+y7**4/2._ark)*y5+ &
      (sqrt(3._ark)*y9**4/2._ark+sqrt(3._ark)*y7**4/2._ark)*y6)*y4
    dF(2,335) = (-sqrt(3._ark)*y5**2*y7**3/2._ark+y5*y6*y9**3+(sqrt(3._ark)*y7**3/6._ark- &
      sqrt(3._ark)*y9**3/3._ark)*y6**2)*y1+(-sqrt(3._ark)*y5**2*y7**3/2._ark-y5*y6*y9**3+ &
      (sqrt(3._ark)*y7**3/6._ark+sqrt(3._ark)*y9**3/3._ark)*y6**2)*y2+(sqrt(3._ark)*y5**2*y7**3/ &
      2._ark+y5*y6*y9**3+(-sqrt(3._ark)*y9**3/3._ark-sqrt(3._ark)*y7**3/6._ark)*y6**2)*y3+ &
      (sqrt(3._ark)*y5**2*y7**3/2._ark-y5*y6*y9**3+(-sqrt(3._ark)*y7**3/6._ark+ &
      sqrt(3._ark)*y9**3/3._ark)*y6**2)*y4
    dF(2,336) = (sqrt(3._ark)*y5**2*y7*y8**2/2._ark+y5*y6*y7*y8**2+(2._ark/ &
      3._ark*sqrt(3._ark)*y8**2*y9+sqrt(3._ark)*y7*y8**2/6._ark)*y6**2)*y1+ &
      (sqrt(3._ark)*y5**2*y7*y8**2/2._ark+y5*y6*y7*y8**2+(-2._ark/3._ark*sqrt(3._ark)*y8**2*y9+ &
      sqrt(3._ark)*y7*y8**2/6._ark)*y6**2)*y2+(-sqrt(3._ark)*y5**2*y7*y8**2/2._ark- &
      y5*y6*y7*y8**2+(2._ark/3._ark*sqrt(3._ark)*y8**2*y9-sqrt(3._ark)*y7*y8**2/ &
      6._ark)*y6**2)*y3+(-sqrt(3._ark)*y5**2*y7*y8**2/2._ark-y5*y6*y7*y8**2+(- &
      sqrt(3._ark)*y7*y8**2/6._ark-2._ark/3._ark*sqrt(3._ark)*y8**2*y9)*y6**2)*y4
    dF(2,337) = (-y5*y6**2*y8**2+y5**3*y8**2/3._ark)*y1+(y5*y6**2*y8**2-y5**3*y8**2/ &
      3._ark)*y2+(y5*y6**2*y8**2-y5**3*y8**2/3._ark)*y3+(-y5*y6**2*y8**2+y5**3*y8**2/ &
      3._ark)*y4
    dF(2,338) = (y5*y6**2*y7*y9-y5**3*y7*y9/3._ark)*y1+(y5*y6**2*y7*y9-y5**3*y7*y9/ &
      3._ark)*y2+(y5*y6**2*y7*y9-y5**3*y7*y9/3._ark)*y3+(y5*y6**2*y7*y9-y5**3*y7*y9/ &
      3._ark)*y4
    dF(2,339) = (12._ark/5._ark*sqrt(3._ark)*y5**4*y7-27._ark/5._ark*y5**3*y6*y9+(24._ark/ &
      5._ark*sqrt(3._ark)*y7+21._ark/5._ark*sqrt(3._ark)*y9)*y6**2*y5**2+(28._ark/5._ark*y7+ &
      y9)*y6**3*y5+13._ark/5._ark*sqrt(3._ark)*y6**4*y9)*y1+(12._ark/5._ark*sqrt(3._ark)*y5**4*y7+ &
      27._ark/5._ark*y5**3*y6*y9+(-21._ark/5._ark*sqrt(3._ark)*y9+24._ark/ &
      5._ark*sqrt(3._ark)*y7)*y6**2*y5**2+(28._ark/5._ark*y7-y9)*y6**3*y5-13._ark/ &
      5._ark*sqrt(3._ark)*y6**4*y9)*y2+(-12._ark/5._ark*sqrt(3._ark)*y5**4*y7-27._ark/ &
      5._ark*y5**3*y6*y9+(-24._ark/5._ark*sqrt(3._ark)*y7+21._ark/ &
      5._ark*sqrt(3._ark)*y9)*y6**2*y5**2+(y9-28._ark/5._ark*y7)*y6**3*y5+13._ark/ &
      5._ark*sqrt(3._ark)*y6**4*y9)*y3+(-12._ark/5._ark*sqrt(3._ark)*y5**4*y7+27._ark/ &
      5._ark*y5**3*y6*y9+(-21._ark/5._ark*sqrt(3._ark)*y9-24._ark/ &
      5._ark*sqrt(3._ark)*y7)*y6**2*y5**2+(-28._ark/5._ark*y7-y9)*y6**3*y5-13._ark/ &
      5._ark*sqrt(3._ark)*y6**4*y9)*y4
    dF(2,340) = (y5*y6**3*y8-sqrt(3._ark)*y5**4*y8/6._ark+sqrt(3._ark)*y5**2*y6**2*y8/ &
      6._ark)*y1+(y5*y6**3*y8-sqrt(3._ark)*y5**4*y8/6._ark+sqrt(3._ark)*y5**2*y6**2*y8/ &
      6._ark)*y2+(y5*y6**3*y8-sqrt(3._ark)*y5**4*y8/6._ark+sqrt(3._ark)*y5**2*y6**2*y8/ &
      6._ark)*y3+(y5*y6**3*y8-sqrt(3._ark)*y5**4*y8/6._ark+sqrt(3._ark)*y5**2*y6**2*y8/ &
      6._ark)*y4
    dF(2,341) = (6._ark*y5*y6**4-3._ark/2._ark*sqrt(3._ark)*y6**5-y5**3*y6**2-3._ark*y5**5- &
      5._ark/2._ark*sqrt(3._ark)*y5**4*y6)*y1+(-6._ark*y5*y6**4+3._ark*y5**5+y5**3*y6**2+3._ark/ &
      2._ark*sqrt(3._ark)*y6**5+5._ark/2._ark*sqrt(3._ark)*y5**4*y6)*y2+(-6._ark*y5*y6**4+ &
      3._ark*y5**5+y5**3*y6**2+3._ark/2._ark*sqrt(3._ark)*y6**5+5._ark/ &
      2._ark*sqrt(3._ark)*y5**4*y6)*y3+(6._ark*y5*y6**4-3._ark/2._ark*sqrt(3._ark)*y6**5- &
      y5**3*y6**2-3._ark*y5**5-5._ark/2._ark*sqrt(3._ark)*y5**4*y6)*y4
    dF(2,342) = y2*y3*y8**4-y1*y4*y8**4
    dF(2,343) = ((-y6*y8*y9**2-sqrt(3._ark)*y5*y8*y9**2)*y2-2._ark*y3*y6*y7**2*y8)*y1- &
      2._ark*y2*y4*y6*y7**2*y8+(-y6*y8*y9**2-sqrt(3._ark)*y5*y8*y9**2)*y4*y3
    dF(2,344) = (-y2*y6*y7*y9**2+(-y6*y7**2*y9/2._ark-sqrt(3._ark)*y5*y7**2*y9/ &
      2._ark)*y3)*y1+(y6*y7**2*y9/2._ark+sqrt(3._ark)*y5*y7**2*y9/2._ark)*y4*y2+ &
      y3*y4*y6*y7*y9**2
    dF(2,345) = (-y2*y6*y8**3+(-y6*y8**3/2._ark-sqrt(3._ark)*y5*y8**3/2._ark)*y3)*y1+(- &
      y6*y8**3/2._ark-sqrt(3._ark)*y5*y8**3/2._ark)*y4*y2-y3*y4*y6*y8**3
    dF(2,346) = (-y2*y6*y7**3+(-sqrt(3._ark)*y5*y9**3/2._ark-y6*y9**3/2._ark)*y3)*y1+ &
      (sqrt(3._ark)*y5*y9**3/2._ark+y6*y9**3/2._ark)*y4*y2+y3*y4*y6*y7**3
    dF(2,347) = ((-y5**2*y8*y9-y6**2*y8*y9)*y2+(-y6**2*y7*y8-y5**2*y7*y8)*y3)*y1+ &
      (y5**2*y7*y8+y6**2*y7*y8)*y4*y2+(y6**2*y8*y9+y5**2*y8*y9)*y4*y3
    dF(2,348) = (y6**2*y7*y9+y5**2*y7*y9)*y4*y1+(y6**2*y7*y9+y5**2*y7*y9)*y3*y2
    dF(2,349) = (-y2*y5*y7*y8**2+(-sqrt(3._ark)*y6*y8**2*y9/2._ark+y5*y8**2*y9/ &
      2._ark)*y3)*y1+(-y5*y8**2*y9/2._ark+sqrt(3._ark)*y6*y8**2*y9/2._ark)*y4*y2+ &
      y3*y4*y5*y7*y8**2
    dF(2,350) = ((sqrt(3._ark)*y6*y7**2*y8+y5*y7**2*y8)*y2+(y5*y8*y9**2+ &
      sqrt(3._ark)*y6*y8*y9**2)*y3)*y1+(y5*y8*y9**2+sqrt(3._ark)*y6*y8*y9**2)*y4*y2+ &
      (sqrt(3._ark)*y6*y7**2*y8+y5*y7**2*y8)*y4*y3
    dF(2,351) = (y2*y5*y8**3+(-y5*y8**3/2._ark+sqrt(3._ark)*y6*y8**3/2._ark)*y3)*y1+(- &
      y5*y8**3/2._ark+sqrt(3._ark)*y6*y8**3/2._ark)*y4*y2+y3*y4*y5*y8**3
    dF(2,352) = (-y2*y5*y7**3+(-sqrt(3._ark)*y6*y9**3/2._ark+y5*y9**3/2._ark)*y3)*y1+ &
      (sqrt(3._ark)*y6*y9**3/2._ark-y5*y9**3/2._ark)*y4*y2+y3*y4*y5*y7**3
    dF(2,353) = (sqrt(3._ark)*y5**2*y7**2/2._ark-y5*y6*y9**2+(-sqrt(3._ark)*y7**2/6._ark+ &
      sqrt(3._ark)*y9**2/3._ark)*y6**2)*y4*y1+(-sqrt(3._ark)*y5**2*y7**2/2._ark+y5*y6*y9**2+(- &
      sqrt(3._ark)*y9**2/3._ark+sqrt(3._ark)*y7**2/6._ark)*y6**2)*y3*y2
    dF(2,354) = ((sqrt(3._ark)*y5**2*y8*y9/3._ark+y5*y6*y8*y9)*y2+ &
      (sqrt(3._ark)*y6**2*y7*y8/2._ark-sqrt(3._ark)*y5**2*y7*y8/6._ark)*y3)*y1+(- &
      sqrt(3._ark)*y6**2*y7*y8/2._ark+sqrt(3._ark)*y5**2*y7*y8/6._ark)*y4*y2+(- &
      sqrt(3._ark)*y5**2*y8*y9/3._ark-y5*y6*y8*y9)*y4*y3
    dF(2,355) = (sqrt(3._ark)*y5*y6**3-y5**4/8._ark+5._ark/4._ark*y5**2*y6**2+3._ark/ &
      8._ark*y6**4)*y4*y1+(y5**4/8._ark-sqrt(3._ark)*y5*y6**3-5._ark/4._ark*y5**2*y6**2-3._ark/ &
      8._ark*y6**4)*y3*y2
    dF(2,356) = ((y6*y7*y9-sqrt(3._ark)*y5*y7*y9)*y2+(sqrt(3._ark)*y5*y7*y9- &
      y6*y7*y9)*y3)*y1**2+((y6*y7*y9-sqrt(3._ark)*y5*y7*y9)*y2**2+(sqrt(3._ark)*y5*y7*y9- &
      y6*y7*y9)*y3**2)*y1+(sqrt(3._ark)*y5*y7*y9-y6*y7*y9)*y4*y2**2+ &
      (sqrt(3._ark)*y5*y7*y9-y6*y7*y9)*y4**2*y2+(y6*y7*y9-sqrt(3._ark)*y5*y7*y9)*y4*y3**2+ &
      (y6*y7*y9-sqrt(3._ark)*y5*y7*y9)*y4**2*y3
    dF(2,357) = ((-sqrt(3._ark)*y5**3/9._ark-11._ark/9._ark*y5**2*y6+y6**3-5._ark/ &
      9._ark*sqrt(3._ark)*y5*y6**2)*y2+(-2._ark/3._ark*y6**3+14._ark/9._ark*y5**2*y6+4._ark/ &
      9._ark*sqrt(3._ark)*y5**3)*y3)*y1**2+((5._ark/9._ark*sqrt(3._ark)*y5*y6**2+ &
      sqrt(3._ark)*y5**3/9._ark+11._ark/9._ark*y5**2*y6-y6**3)*y2**2+(-14._ark/9._ark*y5**2*y6- &
      4._ark/9._ark*sqrt(3._ark)*y5**3+2._ark/3._ark*y6**3)*y3**2)*y1+(-14._ark/9._ark*y5**2*y6- &
      4._ark/9._ark*sqrt(3._ark)*y5**3+2._ark/3._ark*y6**3)*y4*y2**2+(-2._ark/3._ark*y6**3+14._ark/ &
      9._ark*y5**2*y6+4._ark/9._ark*sqrt(3._ark)*y5**3)*y4**2*y2+(5._ark/ &
      9._ark*sqrt(3._ark)*y5*y6**2+sqrt(3._ark)*y5**3/9._ark+11._ark/9._ark*y5**2*y6- &
      y6**3)*y4*y3**2+(-sqrt(3._ark)*y5**3/9._ark-11._ark/9._ark*y5**2*y6+y6**3-5._ark/ &
      9._ark*sqrt(3._ark)*y5*y6**2)*y4**2*y3
    dF(2,358) = -y3**2*y7**2*y9**2-y2**2*y7**2*y9**2+y1**2*y7**2*y9**2+ &
      y4**2*y7**2*y9**2
    dF(2,359) = (-y7**2*y8**2-y8**2*y9**2)*y1**2+(y7**2*y8**2+y8**2*y9**2)*y2**2+ &
      (y7**2*y8**2+y8**2*y9**2)*y3**2+(-y7**2*y8**2-y8**2*y9**2)*y4**2
    dF(2,360) = (-y7*y8*y9**2-y7**2*y8*y9)*y1**2+(-y7**2*y8*y9+y7*y8*y9**2)*y2**2+ &
      (y7**2*y8*y9-y7*y8*y9**2)*y3**2+(y7*y8*y9**2+y7**2*y8*y9)*y4**2
    dF(2,361) = (y7*y9**3+y7**3*y9)*y1**2+(y7*y9**3+y7**3*y9)*y2**2+(y7*y9**3+ &
      y7**3*y9)*y3**2+(y7*y9**3+y7**3*y9)*y4**2
    dF(2,362) = (y8*y9**3+y7**3*y8)*y1**2+(y8*y9**3-y7**3*y8)*y2**2+(y7**3*y8- &
      y8*y9**3)*y3**2+(-y7**3*y8-y8*y9**3)*y4**2
    dF(2,363) = (-y5**2*y6**2-y6**4/2._ark-y5**4/2._ark)*y1**2+(y5**4/2._ark+y5**2*y6**2+ &
      y6**4/2._ark)*y2**2+(y5**4/2._ark+y5**2*y6**2+y6**4/2._ark)*y3**2+(-y5**2*y6**2-y6**4/ &
      2._ark-y5**4/2._ark)*y4**2
    dF(2,364) = ((-y9-y7)*y5**3+(-3._ark*sqrt(3._ark)*y7-3._ark*sqrt(3._ark)*y9)*y6*y5**2+(- &
      9._ark*y9-9._ark*y7)*y6**2*y5+(-3._ark*sqrt(3._ark)*y7-3._ark*sqrt(3._ark)*y9)*y6**3)*y1**2+ &
      ((-y7+y9)*y5**3+(-3._ark*sqrt(3._ark)*y7+3._ark*sqrt(3._ark)*y9)*y6*y5**2+(-9._ark*y7+ &
      9._ark*y9)*y6**2*y5+(-3._ark*sqrt(3._ark)*y7+3._ark*sqrt(3._ark)*y9)*y6**3)*y2**2+((-y9+ &
      y7)*y5**3+(3._ark*sqrt(3._ark)*y7-3._ark*sqrt(3._ark)*y9)*y6*y5**2+(-9._ark*y9+ &
      9._ark*y7)*y6**2*y5+(3._ark*sqrt(3._ark)*y7-3._ark*sqrt(3._ark)*y9)*y6**3)*y3**2+((y9+ &
      y7)*y5**3+(3._ark*sqrt(3._ark)*y9+3._ark*sqrt(3._ark)*y7)*y6*y5**2+(9._ark*y9+ &
      9._ark*y7)*y6**2*y5+(3._ark*sqrt(3._ark)*y9+3._ark*sqrt(3._ark)*y7)*y6**3)*y4**2
    dF(2,365) = (-y3*y7*y8**2-y2*y8**2*y9)*y1**2+(y2**2*y8**2*y9+y3**2*y7*y8**2)*y1- &
      y2**2*y4*y7*y8**2+y2*y4**2*y7*y8**2-y3**2*y4*y8**2*y9+y3*y4**2*y8**2*y9
    dF(2,366) = (y3*y7*y9**2+y2*y7**2*y9)*y1**2+(-y2**2*y7**2*y9-y3**2*y7*y9**2)*y1+ &
      y2**2*y4*y7*y9**2+y3**2*y4*y7**2*y9-y2*y4**2*y7*y9**2-y3*y4**2*y7**2*y9
    dF(2,367) = (-y2**2*y8*y9-y3**2*y7*y8)*y1**2+y3**2*y4**2*y8*y9+y2**2*y4**2*y7*y8
    dF(2,368) = y1**2*y4**2*y7*y9+y2**2*y3**2*y7*y9
    dF(2,369) = (-y9**2-y7**2)*y4**2*y1**2+(y9**2+y7**2)*y3**2*y2**2
    dF(2,370) = (-y2**2*y6*y8+(-y6*y8/2._ark-sqrt(3._ark)*y5*y8/2._ark)*y3**2)*y1**2+(- &
      y6*y8/2._ark-sqrt(3._ark)*y5*y8/2._ark)*y4**2*y2**2-y3**2*y4**2*y6*y8
    dF(2,371) = ((-y5*y7+sqrt(3._ark)*y6*y7)*y2**2+2._ark*y3**2*y5*y9)*y1**2- &
      2._ark*y2**2*y4**2*y5*y9+(-sqrt(3._ark)*y6*y7+y5*y7)*y4**2*y3**2
    dF(2,372) = (y2**2*y5*y8+(sqrt(3._ark)*y6*y8/2._ark-y5*y8/2._ark)*y3**2)*y1**2+ &
      (sqrt(3._ark)*y6*y8/2._ark-y5*y8/2._ark)*y4**2*y2**2+y3**2*y4**2*y5*y8
    dF(2,373) = (sqrt(3._ark)*y6/2._ark+y5/2._ark)*y4**2*y1**3+(sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y4**3*y1**2+(-sqrt(3._ark)*y6/2._ark-y5/2._ark)*y3**2*y2**3+(-sqrt(3._ark)*y6/ &
      2._ark-y5/2._ark)*y3**3*y2**2
    dF(2,374) = (-y3**2-y2**2)*y1**4+(y3**4+y2**4)*y1**2+y3**4*y4**2-y3**2*y4**4- &
      y2**2*y4**4+y2**4*y4**2
    dF(2,375) = ((-3._ark*y7-3._ark*y9)*y5**2+(2._ark*sqrt(3._ark)*y7+ &
      2._ark*sqrt(3._ark)*y9)*y6*y5+(-y9-y7)*y6**2)*y1**3+((-3._ark*y7+3._ark*y9)*y5**2+ &
      (2._ark*sqrt(3._ark)*y7-2._ark*sqrt(3._ark)*y9)*y6*y5+(-y7+y9)*y6**2)*y2**3+((-3._ark*y9+ &
      3._ark*y7)*y5**2+(2._ark*sqrt(3._ark)*y9-2._ark*sqrt(3._ark)*y7)*y6*y5+(-y9+ &
      y7)*y6**2)*y3**3+((3._ark*y9+3._ark*y7)*y5**2+(-2._ark*sqrt(3._ark)*y9- &
      2._ark*sqrt(3._ark)*y7)*y6*y5+(y9+y7)*y6**2)*y4**3
    dF(2,376) = (-y5**3+3._ark*y5*y6**2)*y1**3+(y5**3-3._ark*y5*y6**2)*y2**3+(y5**3- &
      3._ark*y5*y6**2)*y3**3+(-y5**3+3._ark*y5*y6**2)*y4**3
    dF(2,377) = (y3*y7*y9+y2*y7*y9)*y1**3+(y3**3*y7*y9+y2**3*y7*y9)*y1+ &
      y3*y4**3*y7*y9+y3**3*y4*y7*y9+y2*y4**3*y7*y9+y2**3*y4*y7*y9
    dF(2,378) = (-y6**2-y5**2)*y4*y1**3+(-y6**2-y5**2)*y4**3*y1+(y6**2+ &
      y5**2)*y3*y2**3+(y6**2+y5**2)*y3**3*y2
    dF(2,379) = ((sqrt(3._ark)*y6*y9+y5*y9)*y2+(sqrt(3._ark)*y6*y7+y5*y7)*y3)*y1**3+((- &
      sqrt(3._ark)*y6*y9-y5*y9)*y2**3+(-y5*y7-sqrt(3._ark)*y6*y7)*y3**3)*y1+ &
      (sqrt(3._ark)*y6*y7+y5*y7)*y4*y2**3+(-y5*y7-sqrt(3._ark)*y6*y7)*y4**3*y2+ &
      (sqrt(3._ark)*y6*y9+y5*y9)*y4*y3**3+(-sqrt(3._ark)*y6*y9-y5*y9)*y4**3*y3
    dF(2,380) = (-y5*y8/2._ark-sqrt(3._ark)*y6*y8/2._ark)*y4*y1**3+(-y5*y8/2._ark- &
      sqrt(3._ark)*y6*y8/2._ark)*y4**3*y1+(-y5*y8/2._ark-sqrt(3._ark)*y6*y8/2._ark)*y3*y2**3+(- &
      y5*y8/2._ark-sqrt(3._ark)*y6*y8/2._ark)*y3**3*y2
    dF(2,381) = ((-sqrt(3._ark)*y6**2/3._ark-y5*y6)*y2+(-sqrt(3._ark)*y6**2/3._ark- &
      y5*y6)*y3)*y1**3+((y5*y6+sqrt(3._ark)*y6**2/3._ark)*y2**3+(y5*y6+sqrt(3._ark)*y6**2/ &
      3._ark)*y3**3)*y1+(y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4*y2**3+(-sqrt(3._ark)*y6**2/3._ark- &
      y5*y6)*y4**3*y2+(y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4*y3**3+(-sqrt(3._ark)*y6**2/3._ark- &
      y5*y6)*y4**3*y3
    dF(2,382) = ((-y6**2-y5**2)*y2+(-y6**2-y5**2)*y3)*y1**3+((y6**2+y5**2)*y2**3+ &
      (y6**2+y5**2)*y3**3)*y1+(y6**2+y5**2)*y4*y2**3+(-y6**2-y5**2)*y4**3*y2+(y6**2+ &
      y5**2)*y4*y3**3+(-y6**2-y5**2)*y4**3*y3
    dF(2,383) = y1**3*y4**2*y8+y2**3*y3**2*y8+y2**2*y3**3*y8+y1**2*y4**3*y8
    dF(2,384) = (-y2**2*y7-y3**2*y9)*y1**3+(-y3**3*y9-y2**3*y7)*y1**2+ &
      y2**2*y4**3*y9+y2**3*y4**2*y9+y3**2*y4**3*y7+y3**3*y4**2*y7
    dF(2,385) = (y9+y7)*y4**2*y1**3+(-y9-y7)*y4**3*y1**2+(-y9+y7)*y3**2*y2**3+(-y7+ &
      y9)*y3**3*y2**2
    dF(2,386) = ((-y6-sqrt(3._ark)*y5)*y2**2-2._ark*y3**2*y6)*y1**3+((y6+ &
      sqrt(3._ark)*y5)*y2**3+2._ark*y3**3*y6)*y1**2+2._ark*y2**3*y4**2*y6- &
      2._ark*y2**2*y4**3*y6+(y6+sqrt(3._ark)*y5)*y4**2*y3**3+(-y6- &
      sqrt(3._ark)*y5)*y4**3*y3**2
    dF(2,387) = -y1**3*y4**3+y2**3*y3**3
    dF(2,388) = (y9**2+y7**2)*y1**4+(-y9**2-y7**2)*y2**4+(-y9**2-y7**2)*y3**4+ &
      (y9**2+y7**2)*y4**4
    dF(2,389) = (-sqrt(3._ark)*y5/3._ark-y6)*y1**5+(y6+sqrt(3._ark)*y5/3._ark)*y2**5+(y6+ &
      sqrt(3._ark)*y5/3._ark)*y3**5+(-sqrt(3._ark)*y5/3._ark-y6)*y4**5
    dF(2,390) = (-y7**3*y9**2-y7**2*y9**3)*y1+(y7**2*y9**3-y7**3*y9**2)*y2+ &
      (y7**3*y9**2-y7**2*y9**3)*y3+(y7**3*y9**2+y7**2*y9**3)*y4
    dF(2,391) = (y7**2*y8**3+y8**3*y9**2)*y1+(y7**2*y8**3+y8**3*y9**2)*y2+ &
      (y7**2*y8**3+y8**3*y9**2)*y3+(y7**2*y8**3+y8**3*y9**2)*y4
    dF(2,392) = (y7**4*y8+y8*y9**4)*y1+(y7**4*y8+y8*y9**4)*y2+(y7**4*y8+ &
      y8*y9**4)*y3+(y7**4*y8+y8*y9**4)*y4
    dF(2,393) = (y7*y9**4+y7**4*y9)*y1+(-y7**4*y9+y7*y9**4)*y2+(y7**4*y9- &
      y7*y9**4)*y3+(-y7**4*y9-y7*y9**4)*y4
    dF(2,394) = y4*y7*y8**3*y9-y2*y7*y8**3*y9-y3*y7*y8**3*y9+y1*y7*y8**3*y9
    dF(2,395) = ((-sqrt(3._ark)*y9**4/2._ark+sqrt(3._ark)*y7**4/2._ark)*y5+(y9**4/2._ark- &
      y7**4/2._ark)*y6)*y1+((sqrt(3._ark)*y9**4/2._ark-sqrt(3._ark)*y7**4/2._ark)*y5+(y7**4/ &
      2._ark-y9**4/2._ark)*y6)*y2+((sqrt(3._ark)*y9**4/2._ark-sqrt(3._ark)*y7**4/2._ark)*y5+ &
      (y7**4/2._ark-y9**4/2._ark)*y6)*y3+((-sqrt(3._ark)*y9**4/2._ark+sqrt(3._ark)*y7**4/ &
      2._ark)*y5+(y9**4/2._ark-y7**4/2._ark)*y6)*y4
    dF(2,396) = ((-sqrt(3._ark)*y8**3*y9/2._ark+sqrt(3._ark)*y7*y8**3/2._ark)*y5+(y8**3*y9/ &
      2._ark-y7*y8**3/2._ark)*y6)*y1+((-sqrt(3._ark)*y8**3*y9/2._ark-sqrt(3._ark)*y7*y8**3/ &
      2._ark)*y5+(y8**3*y9/2._ark+y7*y8**3/2._ark)*y6)*y2+((sqrt(3._ark)*y8**3*y9/2._ark+ &
      sqrt(3._ark)*y7*y8**3/2._ark)*y5+(-y8**3*y9/2._ark-y7*y8**3/2._ark)*y6)*y3+ &
      ((sqrt(3._ark)*y8**3*y9/2._ark-sqrt(3._ark)*y7*y8**3/2._ark)*y5+(-y8**3*y9/2._ark+ &
      y7*y8**3/2._ark)*y6)*y4
    dF(2,397) = (sqrt(3._ark)*y5*y7*y8*y9**2+(y7*y8*y9**2+2._ark*y7**2*y8*y9)*y6)*y1+(- &
      sqrt(3._ark)*y5*y7*y8*y9**2+(2._ark*y7**2*y8*y9-y7*y8*y9**2)*y6)*y2+ &
      (sqrt(3._ark)*y5*y7*y8*y9**2+(y7*y8*y9**2-2._ark*y7**2*y8*y9)*y6)*y3+(- &
      sqrt(3._ark)*y5*y7*y8*y9**2+(-2._ark*y7**2*y8*y9-y7*y8*y9**2)*y6)*y4
    dF(2,398) = (y6**4*y8+2._ark*y5**2*y6**2*y8+y5**4*y8)*y1+(y6**4*y8+ &
      2._ark*y5**2*y6**2*y8+y5**4*y8)*y2+(y6**4*y8+2._ark*y5**2*y6**2*y8+y5**4*y8)*y3+ &
      (y6**4*y8+2._ark*y5**2*y6**2*y8+y5**4*y8)*y4
    dF(2,399) = (-27._ark/5._ark*y5**4*y7+24._ark/5._ark*sqrt(3._ark)*y5**3*y6*y9+(-36._ark/ &
      5._ark*y9-54._ark/5._ark*y7)*y6**2*y5**2-16._ark/5._ark*sqrt(3._ark)*y5*y6**3*y7+(y7-28._ark/ &
      5._ark*y9)*y6**4)*y1+(-27._ark/5._ark*y5**4*y7-24._ark/5._ark*sqrt(3._ark)*y5**3*y6*y9+(- &
      54._ark/5._ark*y7+36._ark/5._ark*y9)*y6**2*y5**2-16._ark/5._ark*sqrt(3._ark)*y5*y6**3*y7+ &
      (28._ark/5._ark*y9+y7)*y6**4)*y2+(27._ark/5._ark*y5**4*y7+24._ark/ &
      5._ark*sqrt(3._ark)*y5**3*y6*y9+(-36._ark/5._ark*y9+54._ark/5._ark*y7)*y6**2*y5**2+16._ark/ &
      5._ark*sqrt(3._ark)*y5*y6**3*y7+(-28._ark/5._ark*y9-y7)*y6**4)*y3+(27._ark/5._ark*y5**4*y7- &
      24._ark/5._ark*sqrt(3._ark)*y5**3*y6*y9+(36._ark/5._ark*y9+54._ark/5._ark*y7)*y6**2*y5**2+ &
      16._ark/5._ark*sqrt(3._ark)*y5*y6**3*y7+(28._ark/5._ark*y9-y7)*y6**4)*y4
    dF(2,400) = (-sqrt(3._ark)*y6*y8**4/2._ark-y5*y8**4/2._ark)*y1+(y5*y8**4/2._ark+ &
      sqrt(3._ark)*y6*y8**4/2._ark)*y2+(y5*y8**4/2._ark+sqrt(3._ark)*y6*y8**4/2._ark)*y3+(- &
      sqrt(3._ark)*y6*y8**4/2._ark-y5*y8**4/2._ark)*y4
    dF(2,401) = ((-sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7**2*y8/2._ark)*y5**2+ &
      (y8*y9**2+y7**2*y8)*y6*y5+(-sqrt(3._ark)*y8*y9**2/6._ark-sqrt(3._ark)*y7**2*y8/ &
      6._ark)*y6**2)*y1+((-sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7**2*y8/2._ark)*y5**2+ &
      (y8*y9**2+y7**2*y8)*y6*y5+(-sqrt(3._ark)*y8*y9**2/6._ark-sqrt(3._ark)*y7**2*y8/ &
      6._ark)*y6**2)*y2+((-sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7**2*y8/2._ark)*y5**2+ &
      (y8*y9**2+y7**2*y8)*y6*y5+(-sqrt(3._ark)*y8*y9**2/6._ark-sqrt(3._ark)*y7**2*y8/ &
      6._ark)*y6**2)*y3+((-sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7**2*y8/2._ark)*y5**2+ &
      (y8*y9**2+y7**2*y8)*y6*y5+(-sqrt(3._ark)*y8*y9**2/6._ark-sqrt(3._ark)*y7**2*y8/ &
      6._ark)*y6**2)*y4
    dF(2,402) = ((-sqrt(3._ark)*y7*y9**2/3._ark-sqrt(3._ark)*y7**2*y9/3._ark)*y5**2+ &
      (y7*y9**2+y7**2*y9)*y6*y5)*y1+((sqrt(3._ark)*y7**2*y9/3._ark-sqrt(3._ark)*y7*y9**2/ &
      3._ark)*y5**2+(-y7**2*y9+y7*y9**2)*y6*y5)*y2+((-sqrt(3._ark)*y7**2*y9/3._ark+ &
      sqrt(3._ark)*y7*y9**2/3._ark)*y5**2+(y7**2*y9-y7*y9**2)*y6*y5)*y3+ &
      ((sqrt(3._ark)*y7*y9**2/3._ark+sqrt(3._ark)*y7**2*y9/3._ark)*y5**2+(-y7**2*y9- &
      y7*y9**2)*y6*y5)*y4
    dF(2,403) = (sqrt(3._ark)*y5**2*y7*y8**2/2._ark-y5*y6*y8**2*y9+(- &
      sqrt(3._ark)*y7*y8**2/6._ark+sqrt(3._ark)*y8**2*y9/3._ark)*y6**2)*y1+ &
      (sqrt(3._ark)*y5**2*y7*y8**2/2._ark+y5*y6*y8**2*y9+(-sqrt(3._ark)*y8**2*y9/3._ark- &
      sqrt(3._ark)*y7*y8**2/6._ark)*y6**2)*y2+(-sqrt(3._ark)*y5**2*y7*y8**2/2._ark- &
      y5*y6*y8**2*y9+(sqrt(3._ark)*y8**2*y9/3._ark+sqrt(3._ark)*y7*y8**2/6._ark)*y6**2)*y3+(- &
      sqrt(3._ark)*y5**2*y7*y8**2/2._ark+y5*y6*y8**2*y9+(-sqrt(3._ark)*y8**2*y9/3._ark+ &
      sqrt(3._ark)*y7*y8**2/6._ark)*y6**2)*y4
    dF(2,404) = (sqrt(3._ark)*y5**2*y8*y9**2/2._ark-y5*y6*y7**2*y8+(- &
      sqrt(3._ark)*y8*y9**2/6._ark+sqrt(3._ark)*y7**2*y8/3._ark)*y6**2)*y1+ &
      (sqrt(3._ark)*y5**2*y8*y9**2/2._ark-y5*y6*y7**2*y8+(-sqrt(3._ark)*y8*y9**2/6._ark+ &
      sqrt(3._ark)*y7**2*y8/3._ark)*y6**2)*y2+(sqrt(3._ark)*y5**2*y8*y9**2/2._ark- &
      y5*y6*y7**2*y8+(-sqrt(3._ark)*y8*y9**2/6._ark+sqrt(3._ark)*y7**2*y8/3._ark)*y6**2)*y3+ &
      (sqrt(3._ark)*y5**2*y8*y9**2/2._ark-y5*y6*y7**2*y8+(-sqrt(3._ark)*y8*y9**2/6._ark+ &
      sqrt(3._ark)*y7**2*y8/3._ark)*y6**2)*y4
    dF(2,405) = (sqrt(3._ark)*y5**2*y9**3/2._ark-y5*y6*y7**3+(sqrt(3._ark)*y7**3/3._ark- &
      sqrt(3._ark)*y9**3/6._ark)*y6**2)*y1+(-sqrt(3._ark)*y5**2*y9**3/2._ark-y5*y6*y7**3+ &
      (sqrt(3._ark)*y7**3/3._ark+sqrt(3._ark)*y9**3/6._ark)*y6**2)*y2+(sqrt(3._ark)*y5**2*y9**3/ &
      2._ark+y5*y6*y7**3+(-sqrt(3._ark)*y9**3/6._ark-sqrt(3._ark)*y7**3/3._ark)*y6**2)*y3+(- &
      sqrt(3._ark)*y5**2*y9**3/2._ark+y5*y6*y7**3+(-sqrt(3._ark)*y7**3/3._ark+ &
      sqrt(3._ark)*y9**3/6._ark)*y6**2)*y4
    dF(2,406) = (sqrt(3._ark)*y6**2*y8**3/2._ark+sqrt(3._ark)*y5**2*y8**3/6._ark+ &
      y5*y6*y8**3)*y1+(sqrt(3._ark)*y6**2*y8**3/2._ark+sqrt(3._ark)*y5**2*y8**3/6._ark+ &
      y5*y6*y8**3)*y2+(sqrt(3._ark)*y6**2*y8**3/2._ark+sqrt(3._ark)*y5**2*y8**3/6._ark+ &
      y5*y6*y8**3)*y3+(sqrt(3._ark)*y6**2*y8**3/2._ark+sqrt(3._ark)*y5**2*y8**3/6._ark+ &
      y5*y6*y8**3)*y4
    dF(2,407) = (-y5*y6*y7*y8*y9-sqrt(3._ark)*y5**2*y7*y8*y9/6._ark- &
      sqrt(3._ark)*y6**2*y7*y8*y9/2._ark)*y1+(y5*y6*y7*y8*y9+sqrt(3._ark)*y5**2*y7*y8*y9/ &
      6._ark+sqrt(3._ark)*y6**2*y7*y8*y9/2._ark)*y2+(y5*y6*y7*y8*y9+ &
      sqrt(3._ark)*y5**2*y7*y8*y9/6._ark+sqrt(3._ark)*y6**2*y7*y8*y9/2._ark)*y3+(- &
      y5*y6*y7*y8*y9-sqrt(3._ark)*y5**2*y7*y8*y9/6._ark-sqrt(3._ark)*y6**2*y7*y8*y9/2._ark)*y4
    dF(2,408) = (3._ark/16._ark*y5**3*y7**2-5._ark/16._ark*sqrt(3._ark)*y5**2*y6*y9**2+(- &
      y9**2-13._ark/16._ark*y7**2)*y6**2*y5+(-sqrt(3._ark)*y7**2/4._ark-sqrt(3._ark)*y9**2/ &
      16._ark)*y6**3)*y1+(-3._ark/16._ark*y5**3*y7**2+5._ark/16._ark*sqrt(3._ark)*y5**2*y6*y9**2+ &
      (13._ark/16._ark*y7**2+y9**2)*y6**2*y5+(sqrt(3._ark)*y7**2/4._ark+sqrt(3._ark)*y9**2/ &
      16._ark)*y6**3)*y2+(-3._ark/16._ark*y5**3*y7**2+5._ark/16._ark*sqrt(3._ark)*y5**2*y6*y9**2+ &
      (13._ark/16._ark*y7**2+y9**2)*y6**2*y5+(sqrt(3._ark)*y7**2/4._ark+sqrt(3._ark)*y9**2/ &
      16._ark)*y6**3)*y3+(3._ark/16._ark*y5**3*y7**2-5._ark/16._ark*sqrt(3._ark)*y5**2*y6*y9**2+(- &
      y9**2-13._ark/16._ark*y7**2)*y6**2*y5+(-sqrt(3._ark)*y7**2/4._ark-sqrt(3._ark)*y9**2/ &
      16._ark)*y6**3)*y4
    dF(2,409) = ((-y7**3-y9**3)*y5**2+(-y7**3-y9**3)*y6**2)*y1+((y9**3-y7**3)*y5**2+ &
      (y9**3-y7**3)*y6**2)*y2+((y7**3-y9**3)*y5**2+(y7**3-y9**3)*y6**2)*y3+((y7**3+ &
      y9**3)*y5**2+(y7**3+y9**3)*y6**2)*y4
    dF(2,410) = (sqrt(3._ark)*y5**3*y7**2/4._ark+(-y9**2/4._ark+y7**2)*y6*y5**2+ &
      sqrt(3._ark)*y5*y6**2*y7**2/4._ark+3._ark/4._ark*y6**3*y9**2)*y1+(- &
      sqrt(3._ark)*y5**3*y7**2/4._ark+(-y7**2+y9**2/4._ark)*y6*y5**2- &
      sqrt(3._ark)*y5*y6**2*y7**2/4._ark-3._ark/4._ark*y6**3*y9**2)*y2+(- &
      sqrt(3._ark)*y5**3*y7**2/4._ark+(-y7**2+y9**2/4._ark)*y6*y5**2- &
      sqrt(3._ark)*y5*y6**2*y7**2/4._ark-3._ark/4._ark*y6**3*y9**2)*y3+ &
      (sqrt(3._ark)*y5**3*y7**2/4._ark+(-y9**2/4._ark+y7**2)*y6*y5**2+ &
      sqrt(3._ark)*y5*y6**2*y7**2/4._ark+3._ark/4._ark*y6**3*y9**2)*y4
    dF(2,411) = (-sqrt(3._ark)*y5*y6**2*y8**2-y6**3*y8**2-y5**2*y6*y8**2- &
      sqrt(3._ark)*y5**3*y8**2/9._ark)*y1+(sqrt(3._ark)*y5*y6**2*y8**2+y5**2*y6*y8**2+ &
      y6**3*y8**2+sqrt(3._ark)*y5**3*y8**2/9._ark)*y2+(sqrt(3._ark)*y5*y6**2*y8**2+ &
      y5**2*y6*y8**2+y6**3*y8**2+sqrt(3._ark)*y5**3*y8**2/9._ark)*y3+(- &
      sqrt(3._ark)*y5*y6**2*y8**2-y6**3*y8**2-y5**2*y6*y8**2-sqrt(3._ark)*y5**3*y8**2/ &
      9._ark)*y4
    dF(2,412) = (sqrt(3._ark)*y5**3*y7*y8+(4._ark*y7*y8-y8*y9)*y6*y5**2+ &
      sqrt(3._ark)*y5*y6**2*y7*y8+3._ark*y6**3*y8*y9)*y1+(-sqrt(3._ark)*y5**3*y7*y8+(- &
      4._ark*y7*y8-y8*y9)*y6*y5**2-sqrt(3._ark)*y5*y6**2*y7*y8+3._ark*y6**3*y8*y9)*y2+ &
      (sqrt(3._ark)*y5**3*y7*y8+(4._ark*y7*y8+y8*y9)*y6*y5**2+sqrt(3._ark)*y5*y6**2*y7*y8- &
      3._ark*y6**3*y8*y9)*y3+(-sqrt(3._ark)*y5**3*y7*y8+(-4._ark*y7*y8+y8*y9)*y6*y5**2- &
      sqrt(3._ark)*y5*y6**2*y7*y8-3._ark*y6**3*y8*y9)*y4
    dF(2,413) = ((-y9**2-7._ark/16._ark*y7**2)*y5**3-15._ark/ &
      16._ark*sqrt(3._ark)*y5**2*y6*y9**2+9._ark/16._ark*y5*y6**2*y7**2+(-3._ark/ &
      16._ark*sqrt(3._ark)*y9**2-3._ark/4._ark*sqrt(3._ark)*y7**2)*y6**3)*y1+((y9**2+7._ark/ &
      16._ark*y7**2)*y5**3+15._ark/16._ark*sqrt(3._ark)*y5**2*y6*y9**2-9._ark/ &
      16._ark*y5*y6**2*y7**2+(3._ark/4._ark*sqrt(3._ark)*y7**2+3._ark/ &
      16._ark*sqrt(3._ark)*y9**2)*y6**3)*y2+((y9**2+7._ark/16._ark*y7**2)*y5**3+15._ark/ &
      16._ark*sqrt(3._ark)*y5**2*y6*y9**2-9._ark/16._ark*y5*y6**2*y7**2+(3._ark/ &
      4._ark*sqrt(3._ark)*y7**2+3._ark/16._ark*sqrt(3._ark)*y9**2)*y6**3)*y3+((-y9**2-7._ark/ &
      16._ark*y7**2)*y5**3-15._ark/16._ark*sqrt(3._ark)*y5**2*y6*y9**2+9._ark/ &
      16._ark*y5*y6**2*y7**2+(-3._ark/16._ark*sqrt(3._ark)*y9**2-3._ark/ &
      4._ark*sqrt(3._ark)*y7**2)*y6**3)*y4
    dF(2,414) = (13._ark/5._ark*sqrt(3._ark)*y5**4*y7+(-28._ark/5._ark*y9-y7)*y6*y5**3+(21._ark/ &
      5._ark*sqrt(3._ark)*y7+24._ark/5._ark*sqrt(3._ark)*y9)*y6**2*y5**2+27._ark/5._ark*y5*y6**3*y7+ &
      12._ark/5._ark*sqrt(3._ark)*y6**4*y9)*y1+(13._ark/5._ark*sqrt(3._ark)*y5**4*y7+(28._ark/ &
      5._ark*y9-y7)*y6*y5**3+(-24._ark/5._ark*sqrt(3._ark)*y9+21._ark/ &
      5._ark*sqrt(3._ark)*y7)*y6**2*y5**2+27._ark/5._ark*y5*y6**3*y7-12._ark/ &
      5._ark*sqrt(3._ark)*y6**4*y9)*y2+(-13._ark/5._ark*sqrt(3._ark)*y5**4*y7+(y7-28._ark/ &
      5._ark*y9)*y6*y5**3+(24._ark/5._ark*sqrt(3._ark)*y9-21._ark/ &
      5._ark*sqrt(3._ark)*y7)*y6**2*y5**2-27._ark/5._ark*y5*y6**3*y7+12._ark/ &
      5._ark*sqrt(3._ark)*y6**4*y9)*y3+(-13._ark/5._ark*sqrt(3._ark)*y5**4*y7+(28._ark/5._ark*y9+ &
      y7)*y6*y5**3+(-24._ark/5._ark*sqrt(3._ark)*y9-21._ark/5._ark*sqrt(3._ark)*y7)*y6**2*y5**2- &
      27._ark/5._ark*y5*y6**3*y7-12._ark/5._ark*sqrt(3._ark)*y6**4*y9)*y4
    dF(2,415) = (-sqrt(3._ark)*y5**2*y6**2*y8/2._ark+y5**3*y6*y8-sqrt(3._ark)*y5**4*y8/ &
      6._ark)*y1+(-sqrt(3._ark)*y5**2*y6**2*y8/2._ark+y5**3*y6*y8-sqrt(3._ark)*y5**4*y8/ &
      6._ark)*y2+(-sqrt(3._ark)*y5**2*y6**2*y8/2._ark+y5**3*y6*y8-sqrt(3._ark)*y5**4*y8/ &
      6._ark)*y3+(-sqrt(3._ark)*y5**2*y6**2*y8/2._ark+y5**3*y6*y8-sqrt(3._ark)*y5**4*y8/ &
      6._ark)*y4
    dF(2,416) = ((28._ark/5._ark*y7-y9)*y5**4-16._ark/5._ark*sqrt(3._ark)*y5**3*y6*y9+(54._ark/ &
      5._ark*y9+36._ark/5._ark*y7)*y6**2*y5**2+24._ark/5._ark*sqrt(3._ark)*y5*y6**3*y7+27._ark/ &
      5._ark*y6**4*y9)*y1+((28._ark/5._ark*y7+y9)*y5**4+16._ark/5._ark*sqrt(3._ark)*y5**3*y6*y9+ &
      (36._ark/5._ark*y7-54._ark/5._ark*y9)*y6**2*y5**2+24._ark/5._ark*sqrt(3._ark)*y5*y6**3*y7- &
      27._ark/5._ark*y6**4*y9)*y2+((-28._ark/5._ark*y7-y9)*y5**4-16._ark/ &
      5._ark*sqrt(3._ark)*y5**3*y6*y9+(54._ark/5._ark*y9-36._ark/5._ark*y7)*y6**2*y5**2-24._ark/ &
      5._ark*sqrt(3._ark)*y5*y6**3*y7+27._ark/5._ark*y6**4*y9)*y3+((y9-28._ark/5._ark*y7)*y5**4+ &
      16._ark/5._ark*sqrt(3._ark)*y5**3*y6*y9+(-36._ark/5._ark*y7-54._ark/5._ark*y9)*y6**2*y5**2- &
      24._ark/5._ark*sqrt(3._ark)*y5*y6**3*y7-27._ark/5._ark*y6**4*y9)*y4
    dF(2,417) = (5._ark*sqrt(3._ark)*y5**4*y6+3._ark*sqrt(3._ark)*y6**5-15._ark*y5*y6**4+ &
      7._ark*y5**5)*y1+(-5._ark*sqrt(3._ark)*y5**4*y6-7._ark*y5**5-3._ark*sqrt(3._ark)*y6**5+ &
      15._ark*y5*y6**4)*y2+(-5._ark*sqrt(3._ark)*y5**4*y6-7._ark*y5**5-3._ark*sqrt(3._ark)*y6**5+ &
      15._ark*y5*y6**4)*y3+(5._ark*sqrt(3._ark)*y5**4*y6+3._ark*sqrt(3._ark)*y6**5- &
      15._ark*y5*y6**4+7._ark*y5**5)*y4
    dF(2,418) = (-2._ark*y2*y6*y7**2*y8+(-y6*y8*y9**2-sqrt(3._ark)*y5*y8*y9**2)*y3)*y1+ &
      (-y6*y8*y9**2-sqrt(3._ark)*y5*y8*y9**2)*y4*y2-2._ark*y3*y4*y6*y7**2*y8
    dF(2,419) = ((-3._ark/5._ark*sqrt(3._ark)*y6**3*y8-3._ark/5._ark*sqrt(3._ark)*y5**2*y6*y8- &
      4._ark/5._ark*y5**3*y8)*y2+(-3._ark/5._ark*sqrt(3._ark)*y6**3*y8-3._ark/ &
      5._ark*sqrt(3._ark)*y5**2*y6*y8-4._ark/5._ark*y5**3*y8)*y3)*y1+(-3._ark/ &
      5._ark*sqrt(3._ark)*y6**3*y8-3._ark/5._ark*sqrt(3._ark)*y5**2*y6*y8-4._ark/ &
      5._ark*y5**3*y8)*y4*y2+(-3._ark/5._ark*sqrt(3._ark)*y6**3*y8-3._ark/ &
      5._ark*sqrt(3._ark)*y5**2*y6*y8-4._ark/5._ark*y5**3*y8)*y4*y3
    dF(2,420) = ((9._ark/4._ark*y5*y6**2*y7+5._ark/4._ark*y5**3*y7)*y2+(3._ark/ &
      4._ark*sqrt(3._ark)*y5**2*y6*y9+3._ark/4._ark*sqrt(3._ark)*y6**3*y9-y5**3*y9)*y3)*y1+(- &
      3._ark/4._ark*sqrt(3._ark)*y6**3*y9-3._ark/4._ark*sqrt(3._ark)*y5**2*y6*y9+y5**3*y9)*y4*y2+ &
      (-5._ark/4._ark*y5**3*y7-9._ark/4._ark*y5*y6**2*y7)*y4*y3
    dF(2,421) = (y7**3+y9**3)*y4*y1**2+(-y7**3-y9**3)*y4**2*y1+(y7**3- &
      y9**3)*y3*y2**2+(y9**3-y7**3)*y3**2*y2
    dF(2,422) = (-y2*y7**3-y3*y9**3)*y1**2+(-y2**2*y7**3-y3**2*y9**3)*y1+ &
      y2*y4**2*y9**3+y3*y4**2*y7**3+y3**2*y4*y7**3+y2**2*y4*y9**3
    dF(2,423) = (y2*y8**3+y3*y8**3)*y1**2+(y2**2*y8**3+y3**2*y8**3)*y1+ &
      y3**2*y4*y8**3+y2*y4**2*y8**3+y3*y4**2*y8**3+y2**2*y4*y8**3
    dF(2,424) = (-y7**2*y9-y7*y9**2)*y4*y1**2+(y7*y9**2+y7**2*y9)*y4**2*y1+ &
      (y7**2*y9-y7*y9**2)*y3*y2**2+(-y7**2*y9+y7*y9**2)*y3**2*y2
    dF(2,425) = ((sqrt(3._ark)*y6*y8**2/2._ark-y5*y8**2/2._ark)*y2+y3*y5*y8**2)*y1**2+ &
      ((y5*y8**2/2._ark-sqrt(3._ark)*y6*y8**2/2._ark)*y2**2-y3**2*y5*y8**2)*y1- &
      y2**2*y4*y5*y8**2+y2*y4**2*y5*y8**2+(y5*y8**2/2._ark-sqrt(3._ark)*y6*y8**2/ &
      2._ark)*y4*y3**2+(sqrt(3._ark)*y6*y8**2/2._ark-y5*y8**2/2._ark)*y4**2*y3
    dF(2,426) = ((-y5*y7**2+sqrt(3._ark)*y6*y7**2)*y2+2._ark*y3*y5*y9**2)*y1**2+((- &
      sqrt(3._ark)*y6*y7**2+y5*y7**2)*y2**2-2._ark*y3**2*y5*y9**2)*y1- &
      2._ark*y2**2*y4*y5*y9**2+2._ark*y2*y4**2*y5*y9**2+(-sqrt(3._ark)*y6*y7**2+ &
      y5*y7**2)*y4*y3**2+(-y5*y7**2+sqrt(3._ark)*y6*y7**2)*y4**2*y3
    dF(2,427) = ((-y6-sqrt(3._ark)*y5)*y2-2._ark*y3*y6)*y1**4+((y6+sqrt(3._ark)*y5)*y2**4+ &
      2._ark*y3**4*y6)*y1+2._ark*y2**4*y4*y6-2._ark*y2*y4**4*y6+(y6+sqrt(3._ark)*y5)*y4*y3**4+ &
      (-y6-sqrt(3._ark)*y5)*y4**4*y3
    dF(2,428) = ((sqrt(3._ark)*y6+y5)*y2+(sqrt(3._ark)*y6+y5)*y3)*y1**4+((-y5- &
      sqrt(3._ark)*y6)*y2**4+(-y5-sqrt(3._ark)*y6)*y3**4)*y1+(-y5-sqrt(3._ark)*y6)*y4*y2**4+ &
      (sqrt(3._ark)*y6+y5)*y4**4*y2+(-y5-sqrt(3._ark)*y6)*y4*y3**4+(sqrt(3._ark)*y6+ &
      y5)*y4**4*y3
    dF(2,429) = (y2*y7*y9**3+y3*y7**3*y9)*y1+y3*y4*y7*y9**3+y2*y4*y7**3*y9
    dF(2,430) = (-y3*y7**3*y8-y2*y8*y9**3)*y1+y3*y4*y8*y9**3+y2*y4*y7**3*y8
    dF(2,431) = y1*y4*y7*y8**2*y9+y2*y3*y7*y8**2*y9
    dF(2,432) = ((y8*y9**2+y7**2*y8)*y5+(sqrt(3._ark)*y8*y9**2+ &
      sqrt(3._ark)*y7**2*y8)*y6)*y4*y1+((y8*y9**2+y7**2*y8)*y5+(sqrt(3._ark)*y8*y9**2+ &
      sqrt(3._ark)*y7**2*y8)*y6)*y3*y2
    dF(2,433) = (-2._ark*y2*y5*y8*y9**2+(-sqrt(3._ark)*y6*y7**2*y8+y5*y7**2*y8)*y3)*y1+ &
      (-sqrt(3._ark)*y6*y7**2*y8+y5*y7**2*y8)*y4*y2-2._ark*y3*y4*y5*y8*y9**2
    dF(2,434) = (sqrt(3._ark)*y5**2*y8**2/2._ark+sqrt(3._ark)*y6**2*y8**2/6._ark- &
      y5*y6*y8**2)*y4*y1+(-sqrt(3._ark)*y5**2*y8**2/2._ark-sqrt(3._ark)*y6**2*y8**2/6._ark+ &
      y5*y6*y8**2)*y3*y2
    dF(2,435) = (((-y7*y8**2-y8**2*y9)*y3+(-y7*y8**2+y8**2*y9)*y4)*y2+(-y8**2*y9+ &
      y7*y8**2)*y4*y3)*y1+(y8**2*y9+y7*y8**2)*y4*y3*y2
    dF(2,436) = (((-y7**2*y9-y7*y9**2)*y3+(y7**2*y9-y7*y9**2)*y4)*y2+(-y7**2*y9+ &
      y7*y9**2)*y4*y3)*y1+(y7*y9**2+y7**2*y9)*y4*y3*y2
    dF(2,437) = (((y8*y9**2+y7**2*y8)*y3+(y8*y9**2+y7**2*y8)*y4)*y2+(y8*y9**2+ &
      y7**2*y8)*y4*y3)*y1+(y8*y9**2+y7**2*y8)*y4*y3*y2
    dF(2,438) = ((((sqrt(3._ark)*y7**2-sqrt(3._ark)*y9**2)*y5+(-y7**2+y9**2)*y6)*y3+ &
      ((sqrt(3._ark)*y9**2-sqrt(3._ark)*y7**2)*y5+(y7**2-y9**2)*y6)*y4)*y2+ &
      ((sqrt(3._ark)*y9**2-sqrt(3._ark)*y7**2)*y5+(y7**2-y9**2)*y6)*y4*y3)*y1+ &
      ((sqrt(3._ark)*y7**2-sqrt(3._ark)*y9**2)*y5+(-y7**2+y9**2)*y6)*y4*y3*y2
    dF(2,439) = y1**3*y2*y3*y8+(y3**3*y4*y8+y2**3*y4*y8)*y1+y2*y3*y4**3*y8
    dF(2,440) = (y3*y4*y7+y2*y4*y9)*y1**3+(-y2**3*y3*y9+(-y3**3*y7-y4**3*y7)*y2- &
      y3*y4**3*y9)*y1+y2*y3**3*y4*y9+y2**3*y3*y4*y7
    dF(2,441) = ((-sqrt(3._ark)*y5+y6)*y4*y2+(-y6+sqrt(3._ark)*y5)*y4*y3)*y1**3+((-y6+ &
      sqrt(3._ark)*y5)*y3*y2**3+((-sqrt(3._ark)*y5+y6)*y3**3+(-y6+ &
      sqrt(3._ark)*y5)*y4**3)*y2+(-sqrt(3._ark)*y5+y6)*y4**3*y3)*y1+(-sqrt(3._ark)*y5+ &
      y6)*y4*y3*y2**3+(-y6+sqrt(3._ark)*y5)*y4*y3**3*y2
    dF(2,442) = (2._ark*y2*y4*y5+(sqrt(3._ark)*y6-y5)*y4*y3)*y1**3+(-2._ark*y2**3*y3*y5+ &
      ((y5-sqrt(3._ark)*y6)*y3**3+(sqrt(3._ark)*y6-y5)*y4**3)*y2+2._ark*y3*y4**3*y5)*y1+(y5- &
      sqrt(3._ark)*y6)*y4*y3*y2**3-2._ark*y2*y3**3*y4*y5
    dF(2,443) = ((y6*y7**2-sqrt(3._ark)*y5*y7**2)*y2+(sqrt(3._ark)*y5*y9**2- &
      y6*y9**2)*y3)*y1**2+((-y6*y7**2+sqrt(3._ark)*y5*y7**2)*y2**2+(y6*y9**2- &
      sqrt(3._ark)*y5*y9**2)*y3**2)*y1+(y6*y9**2-sqrt(3._ark)*y5*y9**2)*y4*y2**2+ &
      (sqrt(3._ark)*y5*y9**2-y6*y9**2)*y4**2*y2+(-y6*y7**2+ &
      sqrt(3._ark)*y5*y7**2)*y4*y3**2+(y6*y7**2-sqrt(3._ark)*y5*y7**2)*y4**2*y3
    dF(2,444) = (-sqrt(3._ark)*y5*y9**2+(-y9**2-2._ark*y7**2)*y6)*y4*y1**2+(- &
      sqrt(3._ark)*y5*y9**2+(-y9**2-2._ark*y7**2)*y6)*y4**2*y1+(sqrt(3._ark)*y5*y9**2+ &
      (2._ark*y7**2+y9**2)*y6)*y3*y2**2+(sqrt(3._ark)*y5*y9**2+(2._ark*y7**2+ &
      y9**2)*y6)*y3**2*y2
    dF(2,445) = ((-sqrt(3._ark)*y5*y8**2/2._ark-y6*y8**2/2._ark)*y2-y3*y6*y8**2)*y1**2+ &
      ((sqrt(3._ark)*y5*y8**2/2._ark+y6*y8**2/2._ark)*y2**2+y3**2*y6*y8**2)*y1+ &
      y2**2*y4*y6*y8**2-y2*y4**2*y6*y8**2+(sqrt(3._ark)*y5*y8**2/2._ark+y6*y8**2/ &
      2._ark)*y4*y3**2+(-sqrt(3._ark)*y5*y8**2/2._ark-y6*y8**2/2._ark)*y4**2*y3
    dF(2,446) = ((-y6**2*y7-y5**2*y7)*y2+(-y6**2*y9-y5**2*y9)*y3)*y1**2+((-y6**2*y7- &
      y5**2*y7)*y2**2+(-y6**2*y9-y5**2*y9)*y3**2)*y1+(y6**2*y9+y5**2*y9)*y4*y2**2+ &
      (y6**2*y9+y5**2*y9)*y4**2*y2+(y6**2*y7+y5**2*y7)*y4*y3**2+(y6**2*y7+ &
      y5**2*y7)*y4**2*y3
    dF(2,447) = (-4._ark/9._ark*sqrt(3._ark)*y5**3-y5**2*y6-y6**3)*y4*y1**2+(-4._ark/ &
      9._ark*sqrt(3._ark)*y5**3-y5**2*y6-y6**3)*y4**2*y1+(4._ark/9._ark*sqrt(3._ark)*y5**3+ &
      y5**2*y6+y6**3)*y3*y2**2+(4._ark/9._ark*sqrt(3._ark)*y5**3+y5**2*y6+y6**3)*y3**2*y2
    dF(2,448) = ((4._ark/3._ark*sqrt(3._ark)*y5**2*y6+y5*y6**2+y5**3)*y2+(- &
      sqrt(3._ark)*y5**2*y6/3._ark+sqrt(3._ark)*y6**3)*y3)*y1**2+((-y5*y6**2-y5**3-4._ark/ &
      3._ark*sqrt(3._ark)*y5**2*y6)*y2**2+(sqrt(3._ark)*y5**2*y6/3._ark- &
      sqrt(3._ark)*y6**3)*y3**2)*y1+(sqrt(3._ark)*y5**2*y6/3._ark- &
      sqrt(3._ark)*y6**3)*y4*y2**2+(-sqrt(3._ark)*y5**2*y6/3._ark+ &
      sqrt(3._ark)*y6**3)*y4**2*y2+(-y5*y6**2-y5**3-4._ark/ &
      3._ark*sqrt(3._ark)*y5**2*y6)*y4*y3**2+(4._ark/3._ark*sqrt(3._ark)*y5**2*y6+y5*y6**2+ &
      y5**3)*y4**2*y3
    dF(2,449) = (-y2*y4*y8**2-y3*y4*y8**2)*y1**2+(y2**2*y3*y8**2+(y3**2*y8**2- &
      y4**2*y8**2)*y2-y3*y4**2*y8**2)*y1+y2*y3**2*y4*y8**2+y2**2*y3*y4*y8**2
    dF(2,450) = (y2*y4*y7*y8+y3*y4*y8*y9)*y1**2+(-y2**2*y3*y7*y8+(-y4**2*y8*y9- &
      y3**2*y8*y9)*y2-y3*y4**2*y7*y8)*y1+y2*y3**2*y4*y7*y8+y2**2*y3*y4*y8*y9
    dF(2,451) = ((-y3*y7-y4*y7)*y2**2-y2*y3**2*y9-y3**2*y4*y9)*y1**2+ &
      (y2**2*y4**2*y9+y3**2*y4**2*y7)*y1+y2*y3**2*y4**2*y7+y2**2*y3*y4**2*y9
    dF(2,452) = ((-y7+y9)*y4**2*y2+(-y9+y7)*y4**2*y3)*y1**2+(-y9-y7)*y3**2*y2**2*y1+ &
      (y9+y7)*y4*y3**2*y2**2
    dF(2,453) = ((sqrt(3._ark)*y6/2._ark+y5/2._ark)*y4**2*y2+(sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y4**2*y3)*y1**2+(-sqrt(3._ark)*y6/2._ark-y5/2._ark)*y3**2*y2**2*y1+(- &
      sqrt(3._ark)*y6/2._ark-y5/2._ark)*y4*y3**2*y2**2
    dF(2,454) = (y3*y8**2+y2*y8**2)*y1**3+(-y3**3*y8**2-y2**3*y8**2)*y1+ &
      y3*y4**3*y8**2+y2*y4**3*y8**2-y3**3*y4*y8**2-y2**3*y4*y8**2
    dF(2,455) = (y3+y2)*y1**5+(-y2**5-y3**5)*y1+y2*y4**5-y3**5*y4-y2**5*y4+y3*y4**5
    dF(2,456) = (-sqrt(3._ark)*y5*y7*y8**2/2._ark+(-y7*y8**2/2._ark-y8**2*y9)*y6)*y1**2+(- &
      sqrt(3._ark)*y5*y7*y8**2/2._ark+(-y7*y8**2/2._ark+y8**2*y9)*y6)*y2**2+ &
      (sqrt(3._ark)*y5*y7*y8**2/2._ark+(y7*y8**2/2._ark-y8**2*y9)*y6)*y3**2+ &
      (sqrt(3._ark)*y5*y7*y8**2/2._ark+(y7*y8**2/2._ark+y8**2*y9)*y6)*y4**2
    dF(2,457) = (sqrt(3._ark)*y5*y7**2*y8+(y7**2*y8+2._ark*y8*y9**2)*y6)*y1**2+ &
      (sqrt(3._ark)*y5*y7**2*y8+(y7**2*y8+2._ark*y8*y9**2)*y6)*y2**2+ &
      (sqrt(3._ark)*y5*y7**2*y8+(y7**2*y8+2._ark*y8*y9**2)*y6)*y3**2+ &
      (sqrt(3._ark)*y5*y7**2*y8+(y7**2*y8+2._ark*y8*y9**2)*y6)*y4**2
    dF(2,458) = ((y9**2+y7**2)*y5**2+(y9**2+y7**2)*y6**2)*y1**2+((-y9**2- &
      y7**2)*y5**2+(-y9**2-y7**2)*y6**2)*y2**2+((-y9**2-y7**2)*y5**2+(-y9**2- &
      y7**2)*y6**2)*y3**2+((y9**2+y7**2)*y5**2+(y9**2+y7**2)*y6**2)*y4**2
    dF(2,459) = (y6**2*y7*y9+y5**2*y7*y9)*y1**2+(y6**2*y7*y9+y5**2*y7*y9)*y2**2+ &
      (y6**2*y7*y9+y5**2*y7*y9)*y3**2+(y6**2*y7*y9+y5**2*y7*y9)*y4**2
    dF(2,460) = ((y8*y9**2+y7**2*y8)*y5+(sqrt(3._ark)*y8*y9**2+ &
      sqrt(3._ark)*y7**2*y8)*y6)*y1**2+((y8*y9**2+y7**2*y8)*y5+(sqrt(3._ark)*y8*y9**2+ &
      sqrt(3._ark)*y7**2*y8)*y6)*y2**2+((y8*y9**2+y7**2*y8)*y5+(sqrt(3._ark)*y8*y9**2+ &
      sqrt(3._ark)*y7**2*y8)*y6)*y3**2+((y8*y9**2+y7**2*y8)*y5+(sqrt(3._ark)*y8*y9**2+ &
      sqrt(3._ark)*y7**2*y8)*y6)*y4**2
    dF(2,461) = ((2._ark*y7**2*y9-y7*y9**2)*y5+sqrt(3._ark)*y6*y7*y9**2)*y1**2+((- &
      2._ark*y7**2*y9-y7*y9**2)*y5+sqrt(3._ark)*y6*y7*y9**2)*y2**2+((2._ark*y7**2*y9+ &
      y7*y9**2)*y5-sqrt(3._ark)*y6*y7*y9**2)*y3**2+((-2._ark*y7**2*y9+y7*y9**2)*y5- &
      sqrt(3._ark)*y6*y7*y9**2)*y4**2
    dF(2,462) = ((y7*y8**2/2._ark-y8**2*y9)*y5-sqrt(3._ark)*y6*y7*y8**2/2._ark)*y1**2+ &
      ((y7*y8**2/2._ark+y8**2*y9)*y5-sqrt(3._ark)*y6*y7*y8**2/2._ark)*y2**2+((-y7*y8**2/ &
      2._ark-y8**2*y9)*y5+sqrt(3._ark)*y6*y7*y8**2/2._ark)*y3**2+((-y7*y8**2/2._ark+ &
      y8**2*y9)*y5+sqrt(3._ark)*y6*y7*y8**2/2._ark)*y4**2
    dF(2,463) = (sqrt(3._ark)*y5**2*y8*y9/2._ark-y5*y6*y7*y8+(sqrt(3._ark)*y7*y8/3._ark- &
      sqrt(3._ark)*y8*y9/6._ark)*y6**2)*y1**2+(sqrt(3._ark)*y5**2*y8*y9/2._ark+y5*y6*y7*y8+(- &
      sqrt(3._ark)*y8*y9/6._ark-sqrt(3._ark)*y7*y8/3._ark)*y6**2)*y2**2+(- &
      sqrt(3._ark)*y5**2*y8*y9/2._ark-y5*y6*y7*y8+(sqrt(3._ark)*y8*y9/6._ark+sqrt(3._ark)*y7*y8/ &
      3._ark)*y6**2)*y3**2+(-sqrt(3._ark)*y5**2*y8*y9/2._ark+y5*y6*y7*y8+(-sqrt(3._ark)*y7*y8/ &
      3._ark+sqrt(3._ark)*y8*y9/6._ark)*y6**2)*y4**2
    dF(2,464) = ((-sqrt(3._ark)*y9**2/6._ark-sqrt(3._ark)*y7**2/6._ark)*y5**2+(-y9**2- &
      y7**2)*y6*y5+(-sqrt(3._ark)*y7**2/2._ark-sqrt(3._ark)*y9**2/2._ark)*y6**2)*y1**2+ &
      ((sqrt(3._ark)*y7**2/6._ark+sqrt(3._ark)*y9**2/6._ark)*y5**2+(y9**2+y7**2)*y6*y5+ &
      (sqrt(3._ark)*y9**2/2._ark+sqrt(3._ark)*y7**2/2._ark)*y6**2)*y2**2+((sqrt(3._ark)*y7**2/ &
      6._ark+sqrt(3._ark)*y9**2/6._ark)*y5**2+(y9**2+y7**2)*y6*y5+(sqrt(3._ark)*y9**2/2._ark+ &
      sqrt(3._ark)*y7**2/2._ark)*y6**2)*y3**2+((-sqrt(3._ark)*y9**2/6._ark-sqrt(3._ark)*y7**2/ &
      6._ark)*y5**2+(-y9**2-y7**2)*y6*y5+(-sqrt(3._ark)*y7**2/2._ark-sqrt(3._ark)*y9**2/ &
      2._ark)*y6**2)*y4**2
    dF(2,465) = ((-sqrt(3._ark)*y9**2/3._ark+sqrt(3._ark)*y7**2/6._ark)*y5**2-y5*y6*y9**2- &
      sqrt(3._ark)*y6**2*y7**2/2._ark)*y1**2+((-sqrt(3._ark)*y7**2/6._ark+sqrt(3._ark)*y9**2/ &
      3._ark)*y5**2+y5*y6*y9**2+sqrt(3._ark)*y6**2*y7**2/2._ark)*y2**2+((-sqrt(3._ark)*y7**2/ &
      6._ark+sqrt(3._ark)*y9**2/3._ark)*y5**2+y5*y6*y9**2+sqrt(3._ark)*y6**2*y7**2/ &
      2._ark)*y3**2+((-sqrt(3._ark)*y9**2/3._ark+sqrt(3._ark)*y7**2/6._ark)*y5**2-y5*y6*y9**2- &
      sqrt(3._ark)*y6**2*y7**2/2._ark)*y4**2
    dF(2,466) = (y8**2*y9+y7*y8**2)*y4*y1**2+(-y7*y8**2-y8**2*y9)*y4**2*y1+(- &
      y8**2*y9+y7*y8**2)*y3*y2**2+(-y7*y8**2+y8**2*y9)*y3**2*y2
    dF(2,467) = (-y2*y7*y8**2-y3*y8**2*y9)*y1**2+(-y3**2*y8**2*y9- &
      y2**2*y7*y8**2)*y1+y3*y4**2*y7*y8**2+y2*y4**2*y8**2*y9+y3**2*y4*y7*y8**2+ &
      y2**2*y4*y8**2*y9
    dF(2,468) = (y2*y7**2*y8+y3*y8*y9**2)*y1**2+(y3**2*y8*y9**2+y2**2*y7**2*y8)*y1+ &
      y2*y4**2*y8*y9**2+y3**2*y4*y7**2*y8+y2**2*y4*y8*y9**2+y3*y4**2*y7**2*y8
    dF(2,469) = (-sqrt(3._ark)*y5*y8**2/3._ark-y6*y8**2)*y4*y1**2+(-sqrt(3._ark)*y5*y8**2/ &
      3._ark-y6*y8**2)*y4**2*y1+(sqrt(3._ark)*y5*y8**2/3._ark+y6*y8**2)*y3*y2**2+ &
      (sqrt(3._ark)*y5*y8**2/3._ark+y6*y8**2)*y3**2*y2
    dF(2,470) = ((-sqrt(3._ark)*y5*y9**2-y6*y9**2)*y2-2._ark*y3*y6*y7**2)*y1**2+ &
      ((sqrt(3._ark)*y5*y9**2+y6*y9**2)*y2**2+2._ark*y3**2*y6*y7**2)*y1+ &
      2._ark*y2**2*y4*y6*y7**2-2._ark*y2*y4**2*y6*y7**2+(sqrt(3._ark)*y5*y9**2+ &
      y6*y9**2)*y4*y3**2+(-sqrt(3._ark)*y5*y9**2-y6*y9**2)*y4**2*y3
    dF(2,471) = ((-y6*y7*y8/2._ark+sqrt(3._ark)*y5*y7*y8/2._ark)*y2+(y6*y8*y9/2._ark- &
      sqrt(3._ark)*y5*y8*y9/2._ark)*y3)*y1**2+((y6*y7*y8/2._ark-sqrt(3._ark)*y5*y7*y8/ &
      2._ark)*y2**2+(-y6*y8*y9/2._ark+sqrt(3._ark)*y5*y8*y9/2._ark)*y3**2)*y1+(y6*y8*y9/2._ark- &
      sqrt(3._ark)*y5*y8*y9/2._ark)*y4*y2**2+(-y6*y8*y9/2._ark+sqrt(3._ark)*y5*y8*y9/ &
      2._ark)*y4**2*y2+(-y6*y7*y8/2._ark+sqrt(3._ark)*y5*y7*y8/2._ark)*y4*y3**2+(y6*y7*y8/ &
      2._ark-sqrt(3._ark)*y5*y7*y8/2._ark)*y4**2*y3
    dF(2,472) = ((sqrt(3._ark)*y5*y6*y9-y6**2*y9)*y2+(-3._ark/2._ark*y5**2*y7+y6**2*y7/ &
      2._ark)*y3)*y1**2+((-sqrt(3._ark)*y5*y6*y9+y6**2*y9)*y2**2+(3._ark/2._ark*y5**2*y7- &
      y6**2*y7/2._ark)*y3**2)*y1+(-3._ark/2._ark*y5**2*y7+y6**2*y7/2._ark)*y4*y2**2+(3._ark/ &
      2._ark*y5**2*y7-y6**2*y7/2._ark)*y4**2*y2+(sqrt(3._ark)*y5*y6*y9-y6**2*y9)*y4*y3**2+(- &
      sqrt(3._ark)*y5*y6*y9+y6**2*y9)*y4**2*y3
    dF(2,473) = ((-2._ark*sqrt(3._ark)*y5*y6*y8+3._ark*y5**2*y8+y6**2*y8)*y2+(- &
      2._ark*sqrt(3._ark)*y5*y6*y8+3._ark*y5**2*y8+y6**2*y8)*y3)*y1**2+((- &
      2._ark*sqrt(3._ark)*y5*y6*y8+3._ark*y5**2*y8+y6**2*y8)*y2**2+(- &
      2._ark*sqrt(3._ark)*y5*y6*y8+3._ark*y5**2*y8+y6**2*y8)*y3**2)*y1+(- &
      2._ark*sqrt(3._ark)*y5*y6*y8+3._ark*y5**2*y8+y6**2*y8)*y4*y2**2+(- &
      2._ark*sqrt(3._ark)*y5*y6*y8+3._ark*y5**2*y8+y6**2*y8)*y4**2*y2+(- &
      2._ark*sqrt(3._ark)*y5*y6*y8+3._ark*y5**2*y8+y6**2*y8)*y4*y3**2+(- &
      2._ark*sqrt(3._ark)*y5*y6*y8+3._ark*y5**2*y8+y6**2*y8)*y4**2*y3
    dF(2,474) = ((-y5*y9**2-sqrt(3._ark)*y6*y9**2)*y2+(-y5*y7**2- &
      sqrt(3._ark)*y6*y7**2)*y3)*y1**2+((y5*y9**2+sqrt(3._ark)*y6*y9**2)*y2**2+(y5*y7**2+ &
      sqrt(3._ark)*y6*y7**2)*y3**2)*y1+(y5*y7**2+sqrt(3._ark)*y6*y7**2)*y4*y2**2+(- &
      y5*y7**2-sqrt(3._ark)*y6*y7**2)*y4**2*y2+(y5*y9**2+sqrt(3._ark)*y6*y9**2)*y4*y3**2+ &
      (-y5*y9**2-sqrt(3._ark)*y6*y9**2)*y4**2*y3
    dF(2,475) = (-2._ark*y2*y5*y7*y9+(-sqrt(3._ark)*y6*y7*y9+y5*y7*y9)*y3)*y1**2+(- &
      2._ark*y2**2*y5*y7*y9+(-sqrt(3._ark)*y6*y7*y9+y5*y7*y9)*y3**2)*y1+(- &
      sqrt(3._ark)*y6*y7*y9+y5*y7*y9)*y4*y2**2+(-sqrt(3._ark)*y6*y7*y9+y5*y7*y9)*y4**2*y2- &
      2._ark*y3*y4**2*y5*y7*y9-2._ark*y3**2*y4*y5*y7*y9
    dF(2,476) = (-y2*y4*y9**2-y3*y4*y7**2)*y1**2+(y2**2*y3*y9**2+(y3**2*y7**2- &
      y4**2*y7**2)*y2-y3*y4**2*y9**2)*y1+y2**2*y3*y4*y7**2+y2*y3**2*y4*y9**2
    dF(2,477) = (y9**2+y7**2)*y3*y2*y1**2+((-y9**2-y7**2)*y4*y2**2+(-y9**2- &
      y7**2)*y4*y3**2)*y1+(y9**2+y7**2)*y4**2*y3*y2
    dF(2,478) = (-2._ark*y2*y4*y6*y9+(-y6*y7-sqrt(3._ark)*y5*y7)*y4*y3)*y1**2+ &
      (2._ark*y2**2*y3*y6*y9+((sqrt(3._ark)*y5*y7+y6*y7)*y3**2+(sqrt(3._ark)*y5*y7+ &
      y6*y7)*y4**2)*y2+2._ark*y3*y4**2*y6*y9)*y1+(-y6*y7-sqrt(3._ark)*y5*y7)*y4*y3*y2**2- &
      2._ark*y2*y3**2*y4*y6*y9
    dF(2,479) = (2._ark*y2**2*y6*y7+(y6*y9+sqrt(3._ark)*y5*y9)*y3**2)*y1**2+(- &
      sqrt(3._ark)*y5*y9-y6*y9)*y4**2*y2**2-2._ark*y3**2*y4**2*y6*y7
    dF(2,480) = (-y6**2-y5**2)*y4**2*y1**2+(y6**2+y5**2)*y3**2*y2**2
    dF(2,481) = (y5*y8+sqrt(3._ark)*y6*y8)*y4**2*y1**2+(y5*y8+ &
      sqrt(3._ark)*y6*y8)*y3**2*y2**2
    dF(2,482) = ((y4*y8+y3*y8)*y2**2+y2*y3**2*y8+y3**2*y4*y8)*y1**2+(y3**2*y4**2*y8+ &
      y2**2*y4**2*y8)*y1+y2*y3**2*y4**2*y8+y2**2*y3*y4**2*y8
    dF(2,483) = (-2._ark*y2**2*y5+(y5-sqrt(3._ark)*y6)*y3**2)*y1**3+(2._ark*y2**3*y5+ &
      (sqrt(3._ark)*y6-y5)*y3**3)*y1**2+(sqrt(3._ark)*y6-y5)*y4**2*y2**3+(y5- &
      sqrt(3._ark)*y6)*y4**3*y2**2-2._ark*y3**2*y4**3*y5+2._ark*y3**3*y4**2*y5
    dF(2,484) = (-y7**3-y9**3)*y1**3+(y9**3-y7**3)*y2**3+(y7**3-y9**3)*y3**3+(y7**3+ &
      y9**3)*y4**3
    dF(2,485) = y1**3*y8**3+y4**3*y8**3+y3**3*y8**3+y2**3*y8**3
    dF(2,486) = (y8*y9**2+y7**2*y8)*y1**3+(y8*y9**2+y7**2*y8)*y2**3+(y8*y9**2+ &
      y7**2*y8)*y3**3+(y8*y9**2+y7**2*y8)*y4**3
    dF(2,487) = (y7*y9**2+y7**2*y9)*y1**3+(-y7**2*y9+y7*y9**2)*y2**3+(y7**2*y9- &
      y7*y9**2)*y3**3+(-y7**2*y9-y7*y9**2)*y4**3
    dF(2,488) = (-y7*y8**2-y8**2*y9)*y1**3+(-y7*y8**2+y8**2*y9)*y2**3+(-y8**2*y9+ &
      y7*y8**2)*y3**3+(y8**2*y9+y7*y8**2)*y4**3
    dF(2,489) = -y4**3*y7*y8*y9+y3**3*y7*y8*y9-y1**3*y7*y8*y9+y2**3*y7*y8*y9
    dF(2,490) = (-sqrt(3._ark)*y5*y9**2+(-y9**2-2._ark*y7**2)*y6)*y1**3+ &
      (sqrt(3._ark)*y5*y9**2+(2._ark*y7**2+y9**2)*y6)*y2**3+(sqrt(3._ark)*y5*y9**2+ &
      (2._ark*y7**2+y9**2)*y6)*y3**3+(-sqrt(3._ark)*y5*y9**2+(-y9**2-2._ark*y7**2)*y6)*y4**3
    dF(2,491) = ((sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y8/2._ark)*y5+(y7*y8/2._ark-y8*y9/ &
      2._ark)*y6)*y1**3+((sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y8/2._ark)*y5+(-y7*y8/2._ark- &
      y8*y9/2._ark)*y6)*y2**3+((-sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y8/2._ark)*y5+(y7*y8/ &
      2._ark+y8*y9/2._ark)*y6)*y3**3+((-sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y8/2._ark)*y5+ &
      (y8*y9/2._ark-y7*y8/2._ark)*y6)*y4**3
    dF(2,492) = ((y7*y8/2._ark+y8*y9/2._ark)*y5+(sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y8/ &
      2._ark)*y6)*y1**3+((y8*y9/2._ark-y7*y8/2._ark)*y5+(sqrt(3._ark)*y8*y9/2._ark- &
      sqrt(3._ark)*y7*y8/2._ark)*y6)*y2**3+((y7*y8/2._ark-y8*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/ &
      2._ark+sqrt(3._ark)*y7*y8/2._ark)*y6)*y3**3+((-y7*y8/2._ark-y8*y9/2._ark)*y5+(- &
      sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y8/2._ark)*y6)*y4**3
    dF(2,493) = ((-y7**2+2._ark*y9**2)*y5+sqrt(3._ark)*y6*y7**2)*y1**3+((y7**2- &
      2._ark*y9**2)*y5-sqrt(3._ark)*y6*y7**2)*y2**3+((y7**2-2._ark*y9**2)*y5- &
      sqrt(3._ark)*y6*y7**2)*y3**3+((-y7**2+2._ark*y9**2)*y5+sqrt(3._ark)*y6*y7**2)*y4**3
    dF(2,494) = (-y3*y7*y8-y2*y8*y9)*y1**3+(-y3**3*y7*y8-y2**3*y8*y9)*y1+ &
      y3**3*y4*y8*y9+y2*y4**3*y7*y8+y2**3*y4*y7*y8+y3*y4**3*y8*y9
    dF(2,495) = (y3*y7**2+y2*y9**2)*y1**3+(-y3**3*y7**2-y2**3*y9**2)*y1+ &
      y2*y4**3*y7**2-y3**3*y4*y9**2+y3*y4**3*y9**2-y2**3*y4*y7**2
    dF(2,496) = (y3**2*y8+y2**2*y8)*y1**3+(y3**3*y8+y2**3*y8)*y1**2+y3**2*y4**3*y8+ &
      y2**2*y4**3*y8+y2**3*y4**2*y8+y3**3*y4**2*y8
    dF(2,497) = (y2**2*y9+y3**2*y7)*y1**3+(-y2**3*y9-y3**3*y7)*y1**2-y2**2*y4**3*y7+ &
      y2**3*y4**2*y7-y3**2*y4**3*y9+y3**3*y4**2*y9
    dF(2,498) = y1*y8**5+y4*y8**5+y3*y8**5+y2*y8**5
    dF(2,499) = (y7*y8*y9**3+y7**3*y8*y9)*y1+(-y7**3*y8*y9-y7*y8*y9**3)*y2+(- &
      y7**3*y8*y9-y7*y8*y9**3)*y3+(y7*y8*y9**3+y7**3*y8*y9)*y4
    dF(2,500) = y4*y7**2*y8*y9**2+y3*y7**2*y8*y9**2+y2*y7**2*y8*y9**2+ &
      y1*y7**2*y8*y9**2
    dF(2,501) = (sqrt(3._ark)*y5*y7**3*y8+(y7**3*y8+2._ark*y8*y9**3)*y6)*y1+(- &
      sqrt(3._ark)*y5*y7**3*y8+(-y7**3*y8+2._ark*y8*y9**3)*y6)*y2+(sqrt(3._ark)*y5*y7**3*y8+ &
      (-2._ark*y8*y9**3+y7**3*y8)*y6)*y3+(-sqrt(3._ark)*y5*y7**3*y8+(-y7**3*y8- &
      2._ark*y8*y9**3)*y6)*y4
    dF(2,502) = (y6**2*y8**3+y5**2*y8**3)*y1+(y6**2*y8**3+y5**2*y8**3)*y2+ &
      (y6**2*y8**3+y5**2*y8**3)*y3+(y6**2*y8**3+y5**2*y8**3)*y4
    dF(2,503) = ((y7*y9**2+y7**2*y9)*y5**2+(y7*y9**2+y7**2*y9)*y6**2)*y1+((- &
      y7**2*y9+y7*y9**2)*y5**2+(-y7**2*y9+y7*y9**2)*y6**2)*y2+((y7**2*y9- &
      y7*y9**2)*y5**2+(y7**2*y9-y7*y9**2)*y6**2)*y3+((-y7**2*y9-y7*y9**2)*y5**2+(- &
      y7**2*y9-y7*y9**2)*y6**2)*y4
    dF(2,504) = (y5**2*y7*y8*y9+y6**2*y7*y8*y9)*y1+(-y5**2*y7*y8*y9- &
      y6**2*y7*y8*y9)*y2+(-y5**2*y7*y8*y9-y6**2*y7*y8*y9)*y3+(y5**2*y7*y8*y9+ &
      y6**2*y7*y8*y9)*y4
    dF(2,505) = (y5**2*y6*y7*y9+4._ark/9._ark*sqrt(3._ark)*y5**3*y7*y9+y6**3*y7*y9)*y1+ &
      (y5**2*y6*y7*y9+4._ark/9._ark*sqrt(3._ark)*y5**3*y7*y9+y6**3*y7*y9)*y2+ &
      (y5**2*y6*y7*y9+4._ark/9._ark*sqrt(3._ark)*y5**3*y7*y9+y6**3*y7*y9)*y3+ &
      (y5**2*y6*y7*y9+4._ark/9._ark*sqrt(3._ark)*y5**3*y7*y9+y6**3*y7*y9)*y4
    dF(2,506) = (-3._ark/4._ark*sqrt(3._ark)*y5**3*y7**2-9._ark/4._ark*y5**2*y6*y9**2-3._ark/ &
      4._ark*sqrt(3._ark)*y5*y6**2*y7**2+(-5._ark/4._ark*y9**2-y7**2)*y6**3)*y1+(3._ark/ &
      4._ark*sqrt(3._ark)*y5**3*y7**2+9._ark/4._ark*y5**2*y6*y9**2+3._ark/ &
      4._ark*sqrt(3._ark)*y5*y6**2*y7**2+(y7**2+5._ark/4._ark*y9**2)*y6**3)*y2+(3._ark/ &
      4._ark*sqrt(3._ark)*y5**3*y7**2+9._ark/4._ark*y5**2*y6*y9**2+3._ark/ &
      4._ark*sqrt(3._ark)*y5*y6**2*y7**2+(y7**2+5._ark/4._ark*y9**2)*y6**3)*y3+(-3._ark/ &
      4._ark*sqrt(3._ark)*y5**3*y7**2-9._ark/4._ark*y5**2*y6*y9**2-3._ark/ &
      4._ark*sqrt(3._ark)*y5*y6**2*y7**2+(-5._ark/4._ark*y9**2-y7**2)*y6**3)*y4
    dF(2,507) = ((y7**2*y8**2-2._ark*y8**2*y9**2)*y5-sqrt(3._ark)*y6*y7**2*y8**2)*y1+ &
      ((2._ark*y8**2*y9**2-y7**2*y8**2)*y5+sqrt(3._ark)*y6*y7**2*y8**2)*y2+ &
      ((2._ark*y8**2*y9**2-y7**2*y8**2)*y5+sqrt(3._ark)*y6*y7**2*y8**2)*y3+((y7**2*y8**2- &
      2._ark*y8**2*y9**2)*y5-sqrt(3._ark)*y6*y7**2*y8**2)*y4
    dF(2,508) = ((y8*y9**3-2._ark*y7**3*y8)*y5-sqrt(3._ark)*y6*y8*y9**3)*y1+ &
      ((2._ark*y7**3*y8+y8*y9**3)*y5-sqrt(3._ark)*y6*y8*y9**3)*y2+((-2._ark*y7**3*y8- &
      y8*y9**3)*y5+sqrt(3._ark)*y6*y8*y9**3)*y3+((-y8*y9**3+2._ark*y7**3*y8)*y5+ &
      sqrt(3._ark)*y6*y8*y9**3)*y4
    dF(2,509) = (sqrt(3._ark)*y6*y7**2*y9**2+y5*y7**2*y9**2)*y1+(-y5*y7**2*y9**2- &
      sqrt(3._ark)*y6*y7**2*y9**2)*y2+(-y5*y7**2*y9**2-sqrt(3._ark)*y6*y7**2*y9**2)*y3+ &
      (sqrt(3._ark)*y6*y7**2*y9**2+y5*y7**2*y9**2)*y4
    dF(2,510) = ((y7**2*y8*y9-2._ark*y7*y8*y9**2)*y5-sqrt(3._ark)*y6*y7**2*y8*y9)*y1+ &
      ((y7**2*y8*y9+2._ark*y7*y8*y9**2)*y5-sqrt(3._ark)*y6*y7**2*y8*y9)*y2+((- &
      2._ark*y7*y8*y9**2-y7**2*y8*y9)*y5+sqrt(3._ark)*y6*y7**2*y8*y9)*y3+ &
      ((2._ark*y7*y8*y9**2-y7**2*y8*y9)*y5+sqrt(3._ark)*y6*y7**2*y8*y9)*y4
    dF(2,511) = (3._ark/4._ark*y5**3*y7*y8+5._ark/4._ark*sqrt(3._ark)*y5**2*y6*y7*y8+(7._ark/ &
      4._ark*y7*y8+y8*y9)*y6**2*y5+(sqrt(3._ark)*y7*y8/4._ark+sqrt(3._ark)*y8*y9)*y6**3)*y1+(- &
      3._ark/4._ark*y5**3*y7*y8-5._ark/4._ark*sqrt(3._ark)*y5**2*y6*y7*y8+(y8*y9-7._ark/ &
      4._ark*y7*y8)*y6**2*y5+(-sqrt(3._ark)*y7*y8/4._ark+sqrt(3._ark)*y8*y9)*y6**3)*y2+(3._ark/ &
      4._ark*y5**3*y7*y8+5._ark/4._ark*sqrt(3._ark)*y5**2*y6*y7*y8+(-y8*y9+7._ark/ &
      4._ark*y7*y8)*y6**2*y5+(sqrt(3._ark)*y7*y8/4._ark-sqrt(3._ark)*y8*y9)*y6**3)*y3+(-3._ark/ &
      4._ark*y5**3*y7*y8-5._ark/4._ark*sqrt(3._ark)*y5**2*y6*y7*y8+(-7._ark/4._ark*y7*y8- &
      y8*y9)*y6**2*y5+(-sqrt(3._ark)*y7*y8/4._ark-sqrt(3._ark)*y8*y9)*y6**3)*y4
    dF(2,512) = (17._ark/6._ark*sqrt(3._ark)*y5*y6**4-3._ark/2._ark*sqrt(3._ark)*y5**5- &
      3._ark*y5**4*y6-y5**2*y6**3-2._ark*y6**5)*y1+(3._ark/2._ark*sqrt(3._ark)*y5**5+2._ark*y6**5- &
      17._ark/6._ark*sqrt(3._ark)*y5*y6**4+3._ark*y5**4*y6+y5**2*y6**3)*y2+(3._ark/ &
      2._ark*sqrt(3._ark)*y5**5+2._ark*y6**5-17._ark/6._ark*sqrt(3._ark)*y5*y6**4+3._ark*y5**4*y6+ &
      y5**2*y6**3)*y3+(17._ark/6._ark*sqrt(3._ark)*y5*y6**4-3._ark/2._ark*sqrt(3._ark)*y5**5- &
      3._ark*y5**4*y6-y5**2*y6**3-2._ark*y6**5)*y4
    dF(2,513) = ((y8*y9+13._ark/4._ark*y7*y8)*y5**3+15._ark/ &
      4._ark*sqrt(3._ark)*y5**2*y6*y7*y8+9._ark/4._ark*y5*y6**2*y7*y8+(3._ark/ &
      4._ark*sqrt(3._ark)*y7*y8+3._ark*sqrt(3._ark)*y8*y9)*y6**3)*y1+((-13._ark/4._ark*y7*y8+ &
      y8*y9)*y5**3-15._ark/4._ark*sqrt(3._ark)*y5**2*y6*y7*y8-9._ark/4._ark*y5*y6**2*y7*y8+(- &
      3._ark/4._ark*sqrt(3._ark)*y7*y8+3._ark*sqrt(3._ark)*y8*y9)*y6**3)*y2+((-y8*y9+13._ark/ &
      4._ark*y7*y8)*y5**3+15._ark/4._ark*sqrt(3._ark)*y5**2*y6*y7*y8+9._ark/4._ark*y5*y6**2*y7*y8+ &
      (-3._ark*sqrt(3._ark)*y8*y9+3._ark/4._ark*sqrt(3._ark)*y7*y8)*y6**3)*y3+((-y8*y9-13._ark/ &
      4._ark*y7*y8)*y5**3-15._ark/4._ark*sqrt(3._ark)*y5**2*y6*y7*y8-9._ark/4._ark*y5*y6**2*y7*y8+ &
      (-3._ark*sqrt(3._ark)*y8*y9-3._ark/4._ark*sqrt(3._ark)*y7*y8)*y6**3)*y4
    dF(2,514) = (y2*y7*y8**2*y9+y3*y7*y8**2*y9)*y1+y3*y4*y7*y8**2*y9+ &
      y2*y4*y7*y8**2*y9
    dF(2,515) = ((-3._ark/20._ark*sqrt(3._ark)*y5**3*y8-3._ark/4._ark*sqrt(3._ark)*y5*y6**2*y8+ &
      27._ark/20._ark*y5**2*y6*y8+7._ark/20._ark*y6**3*y8)*y2+(y6**3*y8/10._ark+3._ark/ &
      5._ark*sqrt(3._ark)*y5**3*y8-9._ark/10._ark*y5**2*y6*y8)*y3)*y1+(y6**3*y8/10._ark+3._ark/ &
      5._ark*sqrt(3._ark)*y5**3*y8-9._ark/10._ark*y5**2*y6*y8)*y4*y2+(-3._ark/ &
      20._ark*sqrt(3._ark)*y5**3*y8-3._ark/4._ark*sqrt(3._ark)*y5*y6**2*y8+27._ark/ &
      20._ark*y5**2*y6*y8+7._ark/20._ark*y6**3*y8)*y4*y3
    dF(2,516) = ((y5*y6**2*y8-sqrt(3._ark)*y6**3*y8/5._ark-3._ark/5._ark*y5**3*y8- &
      sqrt(3._ark)*y5**2*y6*y8/5._ark)*y2+(y5*y6**2*y8-sqrt(3._ark)*y6**3*y8/5._ark-3._ark/ &
      5._ark*y5**3*y8-sqrt(3._ark)*y5**2*y6*y8/5._ark)*y3)*y1+(y5*y6**2*y8- &
      sqrt(3._ark)*y6**3*y8/5._ark-3._ark/5._ark*y5**3*y8-sqrt(3._ark)*y5**2*y6*y8/5._ark)*y4*y2+ &
      (y5*y6**2*y8-sqrt(3._ark)*y6**3*y8/5._ark-3._ark/5._ark*y5**3*y8-sqrt(3._ark)*y5**2*y6*y8/ &
      5._ark)*y4*y3
    dF(2,517) = ((y5*y6**2*y7/4._ark-3._ark/4._ark*y5**3*y7)*y2+(-sqrt(3._ark)*y5**2*y6*y9/ &
      4._ark+y9*y6**2*y5-sqrt(3._ark)*y6**3*y9/4._ark)*y3)*y1+(sqrt(3._ark)*y5**2*y6*y9/4._ark- &
      y9*y6**2*y5+sqrt(3._ark)*y6**3*y9/4._ark)*y4*y2+(-y5*y6**2*y7/4._ark+3._ark/ &
      4._ark*y5**3*y7)*y4*y3
    dF(2,518) = ((y6**2*y7*y9+y5**2*y7*y9)*y2+(y6**2*y7*y9+y5**2*y7*y9)*y3)*y1+ &
      (y6**2*y7*y9+y5**2*y7*y9)*y4*y2+(y6**2*y7*y9+y5**2*y7*y9)*y4*y3
    dF(2,519) = -y2*y3**5+y1**5*y4+y1*y4**5-y2**5*y3
    dF(2,520) = -y2*y3*y7**2*y9**2+y1*y4*y7**2*y9**2
    dF(2,521) = (y2*y7**3*y9+y3*y7*y9**3)*y1+y2*y4*y7*y9**3+y3*y4*y7**3*y9
    dF(2,522) = (y2*y8**3*y9+y3*y7*y8**3)*y1-y2*y4*y7*y8**3-y3*y4*y8**3*y9
    dF(2,523) = (y7**2*y8**2+y8**2*y9**2)*y4*y1+(-y7**2*y8**2-y8**2*y9**2)*y3*y2
    dF(2,524) = (-sqrt(3._ark)*y5*y6**2*y8-y5**2*y6*y8-y6**3*y8-sqrt(3._ark)*y5**3*y8/ &
      9._ark)*y4*y1+(-sqrt(3._ark)*y5*y6**2*y8-y5**2*y6*y8-y6**3*y8-sqrt(3._ark)*y5**3*y8/ &
      9._ark)*y3*y2
    dF(2,525) = (-y5*y6*y7*y9-sqrt(3._ark)*y5**2*y7*y9/6._ark-sqrt(3._ark)*y6**2*y7*y9/ &
      2._ark)*y4*y1+(-y5*y6*y7*y9-sqrt(3._ark)*y5**2*y7*y9/6._ark-sqrt(3._ark)*y6**2*y7*y9/ &
      2._ark)*y3*y2
    dF(2,526) = ((sqrt(3._ark)*y9**2/2._ark+sqrt(3._ark)*y7**2/2._ark)*y5**2+(-y9**2- &
      y7**2)*y6*y5+(sqrt(3._ark)*y7**2/6._ark+sqrt(3._ark)*y9**2/6._ark)*y6**2)*y4*y1+((- &
      sqrt(3._ark)*y7**2/2._ark-sqrt(3._ark)*y9**2/2._ark)*y5**2+(y9**2+y7**2)*y6*y5+(- &
      sqrt(3._ark)*y9**2/6._ark-sqrt(3._ark)*y7**2/6._ark)*y6**2)*y3*y2
    dF(2,527) = (-2._ark/3._ark*sqrt(3._ark)*y2*y5**2*y8*y9+(-sqrt(3._ark)*y6**2*y7*y8/2._ark+ &
      y5*y6*y7*y8-sqrt(3._ark)*y5**2*y7*y8/6._ark)*y3)*y1+(sqrt(3._ark)*y5**2*y7*y8/6._ark+ &
      sqrt(3._ark)*y6**2*y7*y8/2._ark-y5*y6*y7*y8)*y4*y2+2._ark/ &
      3._ark*sqrt(3._ark)*y3*y4*y5**2*y8*y9
    dF(2,528) = ((-sqrt(3._ark)*y5**2*y7*y9/2._ark+sqrt(3._ark)*y6**2*y7*y9/6._ark)*y2+ &
      (y5*y6*y7*y9-sqrt(3._ark)*y6**2*y7*y9/3._ark)*y3)*y1+(y5*y6*y7*y9- &
      sqrt(3._ark)*y6**2*y7*y9/3._ark)*y4*y2+(-sqrt(3._ark)*y5**2*y7*y9/2._ark+ &
      sqrt(3._ark)*y6**2*y7*y9/6._ark)*y4*y3
    dF(2,529) = ((-9._ark/20._ark*y6**3*y8+sqrt(3._ark)*y5**3*y8/20._ark+ &
      sqrt(3._ark)*y5*y6**2*y8/4._ark+11._ark/20._ark*y5**2*y6*y8)*y2+(3._ark/10._ark*y6**3*y8- &
      sqrt(3._ark)*y5**3*y8/5._ark-7._ark/10._ark*y5**2*y6*y8)*y3)*y1+(3._ark/10._ark*y6**3*y8- &
      sqrt(3._ark)*y5**3*y8/5._ark-7._ark/10._ark*y5**2*y6*y8)*y4*y2+(-9._ark/20._ark*y6**3*y8+ &
      sqrt(3._ark)*y5**3*y8/20._ark+sqrt(3._ark)*y5*y6**2*y8/4._ark+11._ark/ &
      20._ark*y5**2*y6*y8)*y4*y3
    dF(2,530) = (-3._ark/2._ark*sqrt(3._ark)*y5**2*y6**2-y5**3*y6-3._ark*y5*y6**3-3._ark/ &
      4._ark*sqrt(3._ark)*y6**4-sqrt(3._ark)*y5**4/12._ark)*y4*y1+(3._ark/4._ark*sqrt(3._ark)*y6**4+ &
      y5**3*y6+3._ark*y5*y6**3+3._ark/2._ark*sqrt(3._ark)*y5**2*y6**2+sqrt(3._ark)*y5**4/ &
      12._ark)*y3*y2
    dF(2,531) = (sqrt(3._ark)*y5*y7*y8+(y7*y8+2._ark*y8*y9)*y6)*y4*y1**2+(- &
      sqrt(3._ark)*y5*y7*y8+(-2._ark*y8*y9-y7*y8)*y6)*y4**2*y1+(-sqrt(3._ark)*y5*y7*y8+(- &
      y7*y8+2._ark*y8*y9)*y6)*y3*y2**2+(sqrt(3._ark)*y5*y7*y8+(y7*y8- &
      2._ark*y8*y9)*y6)*y3**2*y2
    dF(2,532) = ((-y6*y8*y9+sqrt(3._ark)*y5*y8*y9)*y2+(-sqrt(3._ark)*y5*y7*y8+ &
      y6*y7*y8)*y3)*y1**2+((-y6*y8*y9+sqrt(3._ark)*y5*y8*y9)*y2**2+(- &
      sqrt(3._ark)*y5*y7*y8+y6*y7*y8)*y3**2)*y1+(-y6*y7*y8+ &
      sqrt(3._ark)*y5*y7*y8)*y4*y2**2+(-y6*y7*y8+sqrt(3._ark)*y5*y7*y8)*y4**2*y2+ &
      (y6*y8*y9-sqrt(3._ark)*y5*y8*y9)*y4*y3**2+(y6*y8*y9-sqrt(3._ark)*y5*y8*y9)*y4**2*y3
    dF(2,533) = ((16._ark/9._ark*y5**2*y6-sqrt(3._ark)*y5**3/9._ark-5._ark/ &
      9._ark*sqrt(3._ark)*y5*y6**2)*y2+(-13._ark/9._ark*y5**2*y6+y6**3/3._ark+4._ark/ &
      9._ark*sqrt(3._ark)*y5**3)*y3)*y1**2+((-16._ark/9._ark*y5**2*y6+5._ark/ &
      9._ark*sqrt(3._ark)*y5*y6**2+sqrt(3._ark)*y5**3/9._ark)*y2**2+(-4._ark/ &
      9._ark*sqrt(3._ark)*y5**3-y6**3/3._ark+13._ark/9._ark*y5**2*y6)*y3**2)*y1+(-4._ark/ &
      9._ark*sqrt(3._ark)*y5**3-y6**3/3._ark+13._ark/9._ark*y5**2*y6)*y4*y2**2+(-13._ark/ &
      9._ark*y5**2*y6+y6**3/3._ark+4._ark/9._ark*sqrt(3._ark)*y5**3)*y4**2*y2+(-16._ark/ &
      9._ark*y5**2*y6+5._ark/9._ark*sqrt(3._ark)*y5*y6**2+sqrt(3._ark)*y5**3/9._ark)*y4*y3**2+ &
      (16._ark/9._ark*y5**2*y6-sqrt(3._ark)*y5**3/9._ark-5._ark/ &
      9._ark*sqrt(3._ark)*y5*y6**2)*y4**2*y3
    dF(2,534) = ((2._ark*y7*y8-y8*y9)*y5+sqrt(3._ark)*y6*y8*y9)*y4*y1**2+((y8*y9- &
      2._ark*y7*y8)*y5-sqrt(3._ark)*y6*y8*y9)*y4**2*y1+((-2._ark*y7*y8-y8*y9)*y5+ &
      sqrt(3._ark)*y6*y8*y9)*y3*y2**2+((2._ark*y7*y8+y8*y9)*y5- &
      sqrt(3._ark)*y6*y8*y9)*y3**2*y2
    dF(2,535) = (sqrt(3._ark)*y5**2*y8/2._ark+sqrt(3._ark)*y6**2*y8/6._ark- &
      y5*y6*y8)*y4*y1**2+(sqrt(3._ark)*y5**2*y8/2._ark+sqrt(3._ark)*y6**2*y8/6._ark- &
      y5*y6*y8)*y4**2*y1+(sqrt(3._ark)*y5**2*y8/2._ark+sqrt(3._ark)*y6**2*y8/6._ark- &
      y5*y6*y8)*y3*y2**2+(sqrt(3._ark)*y5**2*y8/2._ark+sqrt(3._ark)*y6**2*y8/6._ark- &
      y5*y6*y8)*y3**2*y2
    dF(2,536) = ((-sqrt(3._ark)*y6**2*y9/2._ark-y5*y6*y9/2._ark)*y2+(-sqrt(3._ark)*y6**2*y7/ &
      4._ark-y5*y6*y7-sqrt(3._ark)*y5**2*y7/4._ark)*y3)*y1**2+((sqrt(3._ark)*y6**2*y9/2._ark+ &
      y5*y6*y9/2._ark)*y2**2+(sqrt(3._ark)*y6**2*y7/4._ark+sqrt(3._ark)*y5**2*y7/4._ark+ &
      y5*y6*y7)*y3**2)*y1+(-sqrt(3._ark)*y6**2*y7/4._ark-y5*y6*y7-sqrt(3._ark)*y5**2*y7/ &
      4._ark)*y4*y2**2+(sqrt(3._ark)*y6**2*y7/4._ark+sqrt(3._ark)*y5**2*y7/4._ark+ &
      y5*y6*y7)*y4**2*y2+(-sqrt(3._ark)*y6**2*y9/2._ark-y5*y6*y9/2._ark)*y4*y3**2+ &
      (sqrt(3._ark)*y6**2*y9/2._ark+y5*y6*y9/2._ark)*y4**2*y3
    dF(2,537) = (y3*y4*y7*y9+y2*y4*y7*y9)*y1**2+(y2**2*y3*y7*y9+(y3**2*y7*y9+ &
      y4**2*y7*y9)*y2+y3*y4**2*y7*y9)*y1+y2*y3**2*y4*y7*y9+y2**2*y3*y4*y7*y9
    dF(2,538) = (-y7*y8-y8*y9)*y3*y2*y1**2+((y7*y8-y8*y9)*y4*y2**2+(y8*y9- &
      y7*y8)*y4*y3**2)*y1+(y8*y9+y7*y8)*y4**2*y3*y2
    dF(2,539) = -y1**2*y2*y3*y8**2+(y2**2*y4*y8**2+y3**2*y4*y8**2)*y1- &
      y2*y3*y4**2*y8**2
    dF(2,540) = (-y3*y4*y9**2-y2*y4*y7**2)*y1**2+(y2**2*y3*y7**2+(y3**2*y9**2- &
      y4**2*y9**2)*y2-y3*y4**2*y7**2)*y1+y2**2*y3*y4*y9**2+y2*y3**2*y4*y7**2
    dF(2,541) = (y2*y7**2*y8*y9+y3*y7*y8*y9**2)*y1-y2*y4*y7*y8*y9**2- &
      y3*y4*y7**2*y8*y9
    dF(2,542) = (y7*y9**3+y7**3*y9)*y4*y1+(y7*y9**3+y7**3*y9)*y3*y2
    dF(2,543) = (y9**4+y7**4)*y4*y1+(-y9**4-y7**4)*y3*y2
    dF(2,544) = (-sqrt(3._ark)*y5*y7**2*y8+(-y7**2*y8-2._ark*y8*y9**2)*y6)*y4*y1+(- &
      sqrt(3._ark)*y5*y7**2*y8+(-y7**2*y8-2._ark*y8*y9**2)*y6)*y3*y2
    dF(2,545) = (-y2*y6*y7*y8**2+(-y6*y8**2*y9/2._ark-sqrt(3._ark)*y5*y8**2*y9/ &
      2._ark)*y3)*y1+(sqrt(3._ark)*y5*y8**2*y9/2._ark+y6*y8**2*y9/2._ark)*y4*y2+ &
      y3*y4*y6*y7*y8**2
    dF(2,546) = ((-3._ark/4._ark*sqrt(3._ark)*y5*y6**2*y7-3._ark/4._ark*sqrt(3._ark)*y5**3*y7+ &
      y6**3*y7)*y2+(3._ark/4._ark*sqrt(3._ark)*y5*y6**2*y9+3._ark/4._ark*sqrt(3._ark)*y5**3*y9- &
      y6**3*y9)*y3)*y1+(y6**3*y9-3._ark/4._ark*sqrt(3._ark)*y5**3*y9-3._ark/ &
      4._ark*sqrt(3._ark)*y5*y6**2*y9)*y4*y2+(3._ark/4._ark*sqrt(3._ark)*y5**3*y7-y6**3*y7+3._ark/ &
      4._ark*sqrt(3._ark)*y5*y6**2*y7)*y4*y3
    dF(2,547) = (9._ark/2._ark*y5**2*y6**2+2._ark*sqrt(3._ark)*y5*y6**3+7._ark/4._ark*y6**4+ &
      3._ark/4._ark*y5**4)*y4*y1+(-3._ark/4._ark*y5**4-7._ark/4._ark*y6**4-9._ark/2._ark*y5**2*y6**2- &
      2._ark*sqrt(3._ark)*y5*y6**3)*y3*y2
    dF(2,548) = (-y5*y8**3/2._ark-sqrt(3._ark)*y6*y8**3/2._ark)*y4*y1+(-y5*y8**3/2._ark- &
      sqrt(3._ark)*y6*y8**3/2._ark)*y3*y2
    dF(2,549) = (y2*y5*y7*y9**2+(-y5*y7**2*y9/2._ark+sqrt(3._ark)*y6*y7**2*y9/ &
      2._ark)*y3)*y1+(y5*y7**2*y9/2._ark-sqrt(3._ark)*y6*y7**2*y9/2._ark)*y4*y2- &
      y3*y4*y5*y7*y9**2
    dF(2,550) = (-sqrt(3._ark)*y6*y7*y8*y9/2._ark-y5*y7*y8*y9/2._ark)*y4*y1+(y5*y7*y8*y9/ &
      2._ark+sqrt(3._ark)*y6*y7*y8*y9/2._ark)*y3*y2
    dF(2,551) = ((y5*y6*y7*y9-sqrt(3._ark)*y5**2*y7*y9/2._ark-sqrt(3._ark)*y6**2*y7*y9/ &
      6._ark)*y2+(y5*y6*y7*y9-sqrt(3._ark)*y5**2*y7*y9/2._ark-sqrt(3._ark)*y6**2*y7*y9/ &
      6._ark)*y3)*y1+(y5*y6*y7*y9-sqrt(3._ark)*y5**2*y7*y9/2._ark-sqrt(3._ark)*y6**2*y7*y9/ &
      6._ark)*y4*y2+(y5*y6*y7*y9-sqrt(3._ark)*y5**2*y7*y9/2._ark-sqrt(3._ark)*y6**2*y7*y9/ &
      6._ark)*y4*y3
    dF(2,552) = (-y5**3*y8/3._ark+y5*y6**2*y8)*y4*y1+(-y5**3*y8/3._ark+ &
      y5*y6**2*y8)*y3*y2
    dF(2,553) = (y6**2*y8**2+y5**2*y8**2)*y4*y1+(-y6**2*y8**2-y5**2*y8**2)*y3*y2
    dF(2,554) = ((y9**2+y7**2)*y5**2+(y9**2+y7**2)*y6**2)*y4*y1+((-y9**2- &
      y7**2)*y5**2+(-y9**2-y7**2)*y6**2)*y3*y2
    dF(2,555) = ((-sqrt(3._ark)*y5*y6**2*y7/4._ark-sqrt(3._ark)*y5**3*y7/4._ark+ &
      y5**2*y6*y7)*y2+(sqrt(3._ark)*y5**3*y9/4._ark-y5**2*y6*y9+sqrt(3._ark)*y5*y6**2*y9/ &
      4._ark)*y3)*y1+(y5**2*y6*y9-sqrt(3._ark)*y5**3*y9/4._ark-sqrt(3._ark)*y5*y6**2*y9/ &
      4._ark)*y4*y2+(-y5**2*y6*y7+sqrt(3._ark)*y5**3*y7/4._ark+sqrt(3._ark)*y5*y6**2*y7/ &
      4._ark)*y4*y3
    dF(2,556) = ((y3*y8**3+y4*y8**3)*y2+y3*y4*y8**3)*y1+y2*y3*y4*y8**3
    dF(2,557) = (((-sqrt(3._ark)*y5*y8*y9/2._ark+(-y8*y9/2._ark-y7*y8)*y6)*y3+(- &
      sqrt(3._ark)*y5*y8*y9/2._ark+(-y8*y9/2._ark+y7*y8)*y6)*y4)*y2+(sqrt(3._ark)*y5*y8*y9/ &
      2._ark+(y8*y9/2._ark-y7*y8)*y6)*y4*y3)*y1+(sqrt(3._ark)*y5*y8*y9/2._ark+(y7*y8+y8*y9/ &
      2._ark)*y6)*y4*y3*y2
    dF(2,558) = ((((-y9**2+2._ark*y7**2)*y5+sqrt(3._ark)*y6*y9**2)*y3+((-2._ark*y7**2+ &
      y9**2)*y5-sqrt(3._ark)*y6*y9**2)*y4)*y2+((-2._ark*y7**2+y9**2)*y5- &
      sqrt(3._ark)*y6*y9**2)*y4*y3)*y1+((-y9**2+2._ark*y7**2)*y5+ &
      sqrt(3._ark)*y6*y9**2)*y4*y3*y2
    dF(2,559) = ((((-y8*y9/2._ark+y7*y8)*y5+sqrt(3._ark)*y6*y8*y9/2._ark)*y3+((-y8*y9/ &
      2._ark-y7*y8)*y5+sqrt(3._ark)*y6*y8*y9/2._ark)*y4)*y2+((y7*y8+y8*y9/2._ark)*y5- &
      sqrt(3._ark)*y6*y8*y9/2._ark)*y4*y3)*y1+((y8*y9/2._ark-y7*y8)*y5-sqrt(3._ark)*y6*y8*y9/ &
      2._ark)*y4*y3*y2
    dF(2,560) = (((-sqrt(3._ark)*y5**2*y9/2._ark+(-y7-2._ark*y9)*y6*y5+(-sqrt(3._ark)*y7- &
      sqrt(3._ark)*y9/2._ark)*y6**2)*y3+(sqrt(3._ark)*y5**2*y9/2._ark+(2._ark*y9-y7)*y6*y5+ &
      (sqrt(3._ark)*y9/2._ark-sqrt(3._ark)*y7)*y6**2)*y4)*y2+(-sqrt(3._ark)*y5**2*y9/2._ark+(- &
      2._ark*y9+y7)*y6*y5+(-sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y7)*y6**2)*y4*y3)*y1+ &
      (sqrt(3._ark)*y5**2*y9/2._ark+(y7+2._ark*y9)*y6*y5+(sqrt(3._ark)*y9/2._ark+ &
      sqrt(3._ark)*y7)*y6**2)*y4*y3*y2
    dF(2,561) = (((y5**3/3._ark-y5*y6**2)*y3+(y5*y6**2-y5**3/3._ark)*y4)*y2+(y5*y6**2- &
      y5**3/3._ark)*y4*y3)*y1+(y5**3/3._ark-y5*y6**2)*y4*y3*y2
    dF(2,562) = (-y3*y4*y9-y2*y4*y7)*y1**3+(-y2**3*y3*y7+(y4**3*y9-y3**3*y9)*y2+ &
      y3*y4**3*y7)*y1+y2**3*y3*y4*y9+y2*y3**3*y4*y7
    dF(2,563) = (sqrt(3._ark)*y6/2._ark+y5/2._ark)*y3*y2*y1**3+((-sqrt(3._ark)*y6/2._ark-y5/ &
      2._ark)*y4*y2**3+(-sqrt(3._ark)*y6/2._ark-y5/2._ark)*y4*y3**3)*y1+(sqrt(3._ark)*y6/2._ark+ &
      y5/2._ark)*y4**3*y3*y2
    dF(2,564) = (((y7**3+y9**3)*y3+(y7**3-y9**3)*y4)*y2+(y9**3-y7**3)*y4*y3)*y1+(- &
      y7**3-y9**3)*y4*y3*y2
    dF(2,565) = ((y4*y7*y8*y9-y3*y7*y8*y9)*y2+y3*y4*y7*y8*y9)*y1-y2*y3*y4*y7*y8*y9
    dF(2,566) = (((y5**2*y8+y6**2*y8)*y3+(y5**2*y8+y6**2*y8)*y4)*y2+(y5**2*y8+ &
      y6**2*y8)*y4*y3)*y1+(y5**2*y8+y6**2*y8)*y4*y3*y2
    dF(2,567) = ((((2._ark*sqrt(3._ark)*y7+2._ark*sqrt(3._ark)*y9)*y6*y5+(2._ark*y7+ &
      2._ark*y9)*y6**2)*y3+((2._ark*sqrt(3._ark)*y7-2._ark*sqrt(3._ark)*y9)*y6*y5+(2._ark*y7- &
      2._ark*y9)*y6**2)*y4)*y2+((2._ark*sqrt(3._ark)*y9-2._ark*sqrt(3._ark)*y7)*y6*y5+(2._ark*y9- &
      2._ark*y7)*y6**2)*y4*y3)*y1+((-2._ark*sqrt(3._ark)*y9-2._ark*sqrt(3._ark)*y7)*y6*y5+(- &
      2._ark*y9-2._ark*y7)*y6**2)*y4*y3*y2
    dF(2,568) = (((sqrt(3._ark)*y5*y6**2+sqrt(3._ark)*y5**3/9._ark+y6**3+y5**2*y6)*y3+(- &
      y5**2*y6-sqrt(3._ark)*y5*y6**2-y6**3-sqrt(3._ark)*y5**3/9._ark)*y4)*y2+(-y5**2*y6- &
      sqrt(3._ark)*y5*y6**2-y6**3-sqrt(3._ark)*y5**3/9._ark)*y4*y3)*y1+(sqrt(3._ark)*y5*y6**2+ &
      sqrt(3._ark)*y5**3/9._ark+y6**3+y5**2*y6)*y4*y3*y2
    dF(2,569) = (((-sqrt(3._ark)*y6*y7*y9/2._ark-y5*y7*y9/2._ark)*y3+(- &
      sqrt(3._ark)*y6*y7*y9/2._ark-y5*y7*y9/2._ark)*y4)*y2+(-sqrt(3._ark)*y6*y7*y9/2._ark- &
      y5*y7*y9/2._ark)*y4*y3)*y1+(-sqrt(3._ark)*y6*y7*y9/2._ark-y5*y7*y9/2._ark)*y4*y3*y2
    dF(2,570) = (((-sqrt(3._ark)*y6*y8**2-y5*y8**2)*y3+(sqrt(3._ark)*y6*y8**2+ &
      y5*y8**2)*y4)*y2+(sqrt(3._ark)*y6*y8**2+y5*y8**2)*y4*y3)*y1+(-sqrt(3._ark)*y6*y8**2- &
      y5*y8**2)*y4*y3*y2
    dF(2,571) = (((-sqrt(3._ark)*y5**2*y8/6._ark-sqrt(3._ark)*y6**2*y8/2._ark-y5*y6*y8)*y3+ &
      (-sqrt(3._ark)*y5**2*y8/6._ark-sqrt(3._ark)*y6**2*y8/2._ark-y5*y6*y8)*y4)*y2+(- &
      sqrt(3._ark)*y5**2*y8/6._ark-sqrt(3._ark)*y6**2*y8/2._ark-y5*y6*y8)*y4*y3)*y1+(- &
      sqrt(3._ark)*y5**2*y8/6._ark-sqrt(3._ark)*y6**2*y8/2._ark-y5*y6*y8)*y4*y3*y2
    dF(2,572) = ((((y9+y7)*y5**2+(2._ark*sqrt(3._ark)*y7+2._ark*sqrt(3._ark)*y9)*y6*y5+ &
      (3._ark*y9+3._ark*y7)*y6**2)*y3+((-y9+y7)*y5**2+(2._ark*sqrt(3._ark)*y7- &
      2._ark*sqrt(3._ark)*y9)*y6*y5+(-3._ark*y9+3._ark*y7)*y6**2)*y4)*y2+((-y7+y9)*y5**2+ &
      (2._ark*sqrt(3._ark)*y9-2._ark*sqrt(3._ark)*y7)*y6*y5+(-3._ark*y7+ &
      3._ark*y9)*y6**2)*y4*y3)*y1+((-y9-y7)*y5**2+(-2._ark*sqrt(3._ark)*y9- &
      2._ark*sqrt(3._ark)*y7)*y6*y5+(-3._ark*y7-3._ark*y9)*y6**2)*y4*y3*y2
    dF(2,573) = y1*y2*y3*y4*y7*y9
    dF(2,574) = (-y6*y8-sqrt(3._ark)*y5*y8/3._ark)*y4*y3*y2*y1
    dF(2,575) = (y3*y4*y8+y2*y4*y8)*y1**3+(y2**3*y3*y8+(y3**3*y8+y4**3*y8)*y2+ &
      y3*y4**3*y8)*y1+y2*y3**3*y4*y8+y2**3*y3*y4*y8
    dF(2,576) = (-y3*y4-y2*y4)*y1**4+(y2**4*y3+(-y4**4+y3**4)*y2-y3*y4**4)*y1+ &
      y2**4*y3*y4+y2*y3**4*y4
    dF(2,577) = -y1**4*y2*y3+(y3**4*y4+y2**4*y4)*y1-y2*y3*y4**4
    dF(2,578) = (y3*y7**2*y8+y2*y8*y9**2)*y1**2+(y3**2*y7**2*y8+y2**2*y8*y9**2)*y1+ &
      y2**2*y4*y7**2*y8+y2*y4**2*y7**2*y8+y3*y4**2*y8*y9**2+y3**2*y4*y8*y9**2
    dF(2,579) = ((y5*y6+sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y4*y2+(y5*y6+ &
      sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y4*y3)*y1**2+((-sqrt(3._ark)*y6**2/ &
      2._ark-sqrt(3._ark)*y5**2/6._ark-y5*y6)*y3*y2**2+((-sqrt(3._ark)*y6**2/2._ark- &
      sqrt(3._ark)*y5**2/6._ark-y5*y6)*y3**2+(y5*y6+sqrt(3._ark)*y6**2/2._ark+ &
      sqrt(3._ark)*y5**2/6._ark)*y4**2)*y2+(y5*y6+sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/ &
      6._ark)*y4**2*y3)*y1+(-sqrt(3._ark)*y6**2/2._ark-sqrt(3._ark)*y5**2/6._ark- &
      y5*y6)*y4*y3*y2**2+(-sqrt(3._ark)*y6**2/2._ark-sqrt(3._ark)*y5**2/6._ark- &
      y5*y6)*y4*y3**2*y2
    dF(2,580) = (-sqrt(3._ark)*y6**2/3._ark-y5*y6)*y3*y2*y1**2+((y5*y6+sqrt(3._ark)*y6**2/ &
      3._ark)*y4*y2**2+(y5*y6+sqrt(3._ark)*y6**2/3._ark)*y4*y3**2)*y1+(-sqrt(3._ark)*y6**2/ &
      3._ark-y5*y6)*y4**2*y3*y2
    dF(2,581) = (y2*y4**2*y8+y3*y4**2*y8)*y1**2+y1*y2**2*y3**2*y8+y2**2*y3**2*y4*y8
    dF(2,582) = ((-y4*y9+y3*y9)*y2**2+y2*y3**2*y7-y3**2*y4*y7)*y1**2+ &
      (y2**2*y4**2*y7+y3**2*y4**2*y9)*y1-y2**2*y3*y4**2*y7-y2*y3**2*y4**2*y9
    dF(2,583) = (((sqrt(3._ark)*y6-y5)*y3+(y5-sqrt(3._ark)*y6)*y4)*y2**2+ &
      2._ark*y2*y3**2*y5-2._ark*y3**2*y4*y5)*y1**2+(-2._ark*y2**2*y4**2*y5+(y5- &
      sqrt(3._ark)*y6)*y4**2*y3**2)*y1+2._ark*y2**2*y3*y4**2*y5+(sqrt(3._ark)*y6- &
      y5)*y4**2*y3**2*y2
    dF(2,584) = (y3**2*y4+y2**2*y4)*y1**3+(-y2**3*y3-y2*y3**3)*y1**2+(y2**2*y4**3+ &
      y3**2*y4**3)*y1-y2**3*y3*y4**2-y2*y3**3*y4**2
    dF(2,585) = y1**2*y2*y3*y7*y9+(y3**2*y4*y7*y9+y2**2*y4*y7*y9)*y1+ &
      y2*y3*y4**2*y7*y9
    dF(2,586) = (sqrt(3._ark)*y5*y9+(y9+2._ark*y7)*y6)*y3*y2*y1**2+((-sqrt(3._ark)*y5*y9+ &
      (-y9+2._ark*y7)*y6)*y4*y2**2+(sqrt(3._ark)*y5*y9+(y9-2._ark*y7)*y6)*y4*y3**2)*y1+(- &
      sqrt(3._ark)*y5*y9+(-y9-2._ark*y7)*y6)*y4**2*y3*y2
    dF(2,587) = (-y9-y7)*y4*y3*y2*y1**2+((-y7+y9)*y4*y3*y2**2+((-y9+y7)*y4*y3**2+ &
      (y9+y7)*y4**2*y3)*y2)*y1
    dF(2,588) = y1**2*y2*y3*y4*y8+(y2**2*y3*y4*y8+(y3*y4**2*y8+y3**2*y4*y8)*y2)*y1
    dF(2,589) = (-y5-sqrt(3._ark)*y6)*y4*y3*y2*y1**2+((sqrt(3._ark)*y6+y5)*y4*y3*y2**2+ &
      ((sqrt(3._ark)*y6+y5)*y4*y3**2+(-y5-sqrt(3._ark)*y6)*y4**2*y3)*y2)*y1
    dF(2,590) = y1**2*y8**4+y4**2*y8**4-y3**2*y8**4-y2**2*y8**4
    dF(2,591) = (y8**3*y9+y7*y8**3)*y1**2+(y8**3*y9-y7*y8**3)*y2**2+(y7*y8**3- &
      y8**3*y9)*y3**2+(-y7*y8**3-y8**3*y9)*y4**2
    dF(2,592) = (y9**4+y7**4)*y1**2+(-y9**4-y7**4)*y2**2+(-y9**4-y7**4)*y3**2+ &
      (y9**4+y7**4)*y4**2
    dF(2,593) = ((sqrt(3._ark)*y9**3-sqrt(3._ark)*y7**3)*y5+(y7**3-y9**3)*y6)*y1**2+((- &
      sqrt(3._ark)*y9**3-sqrt(3._ark)*y7**3)*y5+(y7**3+y9**3)*y6)*y2**2+ &
      ((sqrt(3._ark)*y7**3+sqrt(3._ark)*y9**3)*y5+(-y7**3-y9**3)*y6)*y3**2+ &
      ((sqrt(3._ark)*y7**3-sqrt(3._ark)*y9**3)*y5+(y9**3-y7**3)*y6)*y4**2
    dF(2,594) = (-3._ark/4._ark*sqrt(3._ark)*y5**3*y7+(-9._ark/2._ark*y9-9._ark/ &
      4._ark*y7)*y6*y5**2+(-3._ark*sqrt(3._ark)*y9-15._ark/4._ark*sqrt(3._ark)*y7)*y6**2*y5+(- &
      7._ark/2._ark*y9-13._ark/4._ark*y7)*y6**3)*y1**2+(-3._ark/4._ark*sqrt(3._ark)*y5**3*y7+(-9._ark/ &
      4._ark*y7+9._ark/2._ark*y9)*y6*y5**2+(3._ark*sqrt(3._ark)*y9-15._ark/ &
      4._ark*sqrt(3._ark)*y7)*y6**2*y5+(7._ark/2._ark*y9-13._ark/4._ark*y7)*y6**3)*y2**2+(3._ark/ &
      4._ark*sqrt(3._ark)*y5**3*y7+(9._ark/4._ark*y7-9._ark/2._ark*y9)*y6*y5**2+(15._ark/ &
      4._ark*sqrt(3._ark)*y7-3._ark*sqrt(3._ark)*y9)*y6**2*y5+(13._ark/4._ark*y7-7._ark/ &
      2._ark*y9)*y6**3)*y3**2+(3._ark/4._ark*sqrt(3._ark)*y5**3*y7+(9._ark/4._ark*y7+9._ark/ &
      2._ark*y9)*y6*y5**2+(15._ark/4._ark*sqrt(3._ark)*y7+3._ark*sqrt(3._ark)*y9)*y6**2*y5+(13._ark/ &
      4._ark*y7+7._ark/2._ark*y9)*y6**3)*y4**2
    dF(2,595) = (-sqrt(3._ark)*y6*y7*y8*y9/2._ark-y5*y7*y8*y9/2._ark)*y1**2+(y5*y7*y8*y9/ &
      2._ark+sqrt(3._ark)*y6*y7*y8*y9/2._ark)*y2**2+(y5*y7*y8*y9/2._ark+ &
      sqrt(3._ark)*y6*y7*y8*y9/2._ark)*y3**2+(-sqrt(3._ark)*y6*y7*y8*y9/2._ark-y5*y7*y8*y9/ &
      2._ark)*y4**2
    dF(2,596) = (y5*y8**3+sqrt(3._ark)*y6*y8**3)*y1**2+(y5*y8**3+ &
      sqrt(3._ark)*y6*y8**3)*y2**2+(y5*y8**3+sqrt(3._ark)*y6*y8**3)*y3**2+(y5*y8**3+ &
      sqrt(3._ark)*y6*y8**3)*y4**2
    dF(2,597) = ((y7**3-2._ark*y9**3)*y5-sqrt(3._ark)*y6*y7**3)*y1**2+((y7**3+ &
      2._ark*y9**3)*y5-sqrt(3._ark)*y6*y7**3)*y2**2+((-y7**3-2._ark*y9**3)*y5+ &
      sqrt(3._ark)*y6*y7**3)*y3**2+((-y7**3+2._ark*y9**3)*y5+sqrt(3._ark)*y6*y7**3)*y4**2
    dF(2,598) = (-sqrt(3._ark)*y5**2*y8**2/2._ark-sqrt(3._ark)*y6**2*y8**2/6._ark+ &
      y5*y6*y8**2)*y1**2+(sqrt(3._ark)*y5**2*y8**2/2._ark+sqrt(3._ark)*y6**2*y8**2/6._ark- &
      y5*y6*y8**2)*y2**2+(sqrt(3._ark)*y5**2*y8**2/2._ark+sqrt(3._ark)*y6**2*y8**2/6._ark- &
      y5*y6*y8**2)*y3**2+(-sqrt(3._ark)*y5**2*y8**2/2._ark-sqrt(3._ark)*y6**2*y8**2/6._ark+ &
      y5*y6*y8**2)*y4**2
    dF(2,599) = ((sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y8/2._ark)*y5**2+(-y7*y8- &
      y8*y9)*y6*y5+(sqrt(3._ark)*y8*y9/6._ark+sqrt(3._ark)*y7*y8/6._ark)*y6**2)*y1**2+ &
      ((sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y8/2._ark)*y5**2+(y7*y8-y8*y9)*y6*y5+ &
      (sqrt(3._ark)*y8*y9/6._ark-sqrt(3._ark)*y7*y8/6._ark)*y6**2)*y2**2+((-sqrt(3._ark)*y8*y9/ &
      2._ark+sqrt(3._ark)*y7*y8/2._ark)*y5**2+(y8*y9-y7*y8)*y6*y5+(sqrt(3._ark)*y7*y8/6._ark- &
      sqrt(3._ark)*y8*y9/6._ark)*y6**2)*y3**2+((-sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y8/ &
      2._ark)*y5**2+(y8*y9+y7*y8)*y6*y5+(-sqrt(3._ark)*y8*y9/6._ark-sqrt(3._ark)*y7*y8/ &
      6._ark)*y6**2)*y4**2
    dF(2,600) = (sqrt(3._ark)*y5**2*y7*y9/6._ark+sqrt(3._ark)*y6**2*y7*y9/2._ark+ &
      y5*y6*y7*y9)*y1**2+(sqrt(3._ark)*y5**2*y7*y9/6._ark+sqrt(3._ark)*y6**2*y7*y9/2._ark+ &
      y5*y6*y7*y9)*y2**2+(sqrt(3._ark)*y5**2*y7*y9/6._ark+sqrt(3._ark)*y6**2*y7*y9/2._ark+ &
      y5*y6*y7*y9)*y3**2+(sqrt(3._ark)*y5**2*y7*y9/6._ark+sqrt(3._ark)*y6**2*y7*y9/2._ark+ &
      y5*y6*y7*y9)*y4**2
    dF(2,601) = ((sqrt(3._ark)*y7+sqrt(3._ark)*y9)*y6*y5**2+(4._ark*y7+4._ark*y9)*y6**2*y5+ &
      (sqrt(3._ark)*y7+sqrt(3._ark)*y9)*y6**3)*y1**2+((sqrt(3._ark)*y7- &
      sqrt(3._ark)*y9)*y6*y5**2+(4._ark*y7-4._ark*y9)*y6**2*y5+(sqrt(3._ark)*y7- &
      sqrt(3._ark)*y9)*y6**3)*y2**2+((-sqrt(3._ark)*y7+sqrt(3._ark)*y9)*y6*y5**2+(-4._ark*y7+ &
      4._ark*y9)*y6**2*y5+(-sqrt(3._ark)*y7+sqrt(3._ark)*y9)*y6**3)*y3**2+((-sqrt(3._ark)*y9- &
      sqrt(3._ark)*y7)*y6*y5**2+(-4._ark*y7-4._ark*y9)*y6**2*y5+(-sqrt(3._ark)*y9- &
      sqrt(3._ark)*y7)*y6**3)*y4**2
    dF(2,602) = (-y5*y6**3-7._ark/24._ark*sqrt(3._ark)*y6**4-sqrt(3._ark)*y5**4/8._ark-3._ark/ &
      4._ark*sqrt(3._ark)*y5**2*y6**2)*y1**2+(y5*y6**3+7._ark/24._ark*sqrt(3._ark)*y6**4+3._ark/ &
      4._ark*sqrt(3._ark)*y5**2*y6**2+sqrt(3._ark)*y5**4/8._ark)*y2**2+(y5*y6**3+7._ark/ &
      24._ark*sqrt(3._ark)*y6**4+3._ark/4._ark*sqrt(3._ark)*y5**2*y6**2+sqrt(3._ark)*y5**4/ &
      8._ark)*y3**2+(-y5*y6**3-7._ark/24._ark*sqrt(3._ark)*y6**4-sqrt(3._ark)*y5**4/8._ark-3._ark/ &
      4._ark*sqrt(3._ark)*y5**2*y6**2)*y4**2
    dF(2,603) = (y6**2*y8**2+y5**2*y8**2)*y1**2+(-y6**2*y8**2-y5**2*y8**2)*y2**2+(- &
      y6**2*y8**2-y5**2*y8**2)*y3**2+(y6**2*y8**2+y5**2*y8**2)*y4**2
    dF(2,604) = ((y8*y9+y7*y8)*y5**2+(y8*y9+y7*y8)*y6**2)*y1**2+((y8*y9- &
      y7*y8)*y5**2+(y8*y9-y7*y8)*y6**2)*y2**2+((y7*y8-y8*y9)*y5**2+(y7*y8- &
      y8*y9)*y6**2)*y3**2+((-y7*y8-y8*y9)*y5**2+(-y7*y8-y8*y9)*y6**2)*y4**2
    dF(2,605) = (-sqrt(3._ark)*y5**3*y7/4._ark+(-7._ark/4._ark*y7-y9/2._ark)*y6*y5**2+(-5._ark/ &
      4._ark*sqrt(3._ark)*y7-sqrt(3._ark)*y9)*y6**2*y5+(-3._ark/4._ark*y7-3._ark/ &
      2._ark*y9)*y6**3)*y1**2+(-sqrt(3._ark)*y5**3*y7/4._ark+(-7._ark/4._ark*y7+y9/ &
      2._ark)*y6*y5**2+(-5._ark/4._ark*sqrt(3._ark)*y7+sqrt(3._ark)*y9)*y6**2*y5+(3._ark/2._ark*y9- &
      3._ark/4._ark*y7)*y6**3)*y2**2+(sqrt(3._ark)*y5**3*y7/4._ark+(-y9/2._ark+7._ark/ &
      4._ark*y7)*y6*y5**2+(5._ark/4._ark*sqrt(3._ark)*y7-sqrt(3._ark)*y9)*y6**2*y5+(3._ark/ &
      4._ark*y7-3._ark/2._ark*y9)*y6**3)*y3**2+(sqrt(3._ark)*y5**3*y7/4._ark+(y9/2._ark+7._ark/ &
      4._ark*y7)*y6*y5**2+(sqrt(3._ark)*y9+5._ark/4._ark*sqrt(3._ark)*y7)*y6**2*y5+(3._ark/ &
      2._ark*y9+3._ark/4._ark*y7)*y6**3)*y4**2
    dF(2,606) = (-sqrt(3._ark)*y6**4/8._ark-3._ark/4._ark*sqrt(3._ark)*y5**2*y6**2-7._ark/ &
      24._ark*sqrt(3._ark)*y5**4+y5**3*y6)*y1**2+(-y5**3*y6+sqrt(3._ark)*y6**4/8._ark+3._ark/ &
      4._ark*sqrt(3._ark)*y5**2*y6**2+7._ark/24._ark*sqrt(3._ark)*y5**4)*y2**2+(-y5**3*y6+ &
      sqrt(3._ark)*y6**4/8._ark+3._ark/4._ark*sqrt(3._ark)*y5**2*y6**2+7._ark/ &
      24._ark*sqrt(3._ark)*y5**4)*y3**2+(-sqrt(3._ark)*y6**4/8._ark-3._ark/ &
      4._ark*sqrt(3._ark)*y5**2*y6**2-7._ark/24._ark*sqrt(3._ark)*y5**4+y5**3*y6)*y4**2
    dF(2,607) = y1**2*y4*y8**3+y1*y4**2*y8**3+y2**2*y3*y8**3+y2*y3**2*y8**3
    dF(2,608) = (y2*y9**3+y3*y7**3)*y1**2+(-y2**2*y9**3-y3**2*y7**3)*y1- &
      y3*y4**2*y9**3-y2*y4**2*y7**3+y2**2*y4*y7**3+y3**2*y4*y9**3
    dF(2,609) = y1**2*y4*y7*y8*y9-y2*y3**2*y7*y8*y9-y2**2*y3*y7*y8*y9+ &
      y1*y4**2*y7*y8*y9
    dF(2,610) = (y2*y7*y8*y9+y3*y7*y8*y9)*y1**2+(-y2**2*y7*y8*y9-y3**2*y7*y8*y9)*y1- &
      y2**2*y4*y7*y8*y9+y3*y4**2*y7*y8*y9-y3**2*y4*y7*y8*y9+y2*y4**2*y7*y8*y9
    dF(2,611) = (y2*y7*y9**2+y3*y7**2*y9)*y1**2+(y3**2*y7**2*y9+y2**2*y7*y9**2)*y1- &
      y2*y4**2*y7**2*y9-y3**2*y4*y7*y9**2-y2**2*y4*y7**2*y9-y3*y4**2*y7*y9**2
    dF(2,612) = (sqrt(3._ark)*y5*y7*y9/3._ark+y6*y7*y9)*y4*y1**2+(sqrt(3._ark)*y5*y7*y9/ &
      3._ark+y6*y7*y9)*y4**2*y1+(sqrt(3._ark)*y5*y7*y9/3._ark+y6*y7*y9)*y3*y2**2+ &
      (sqrt(3._ark)*y5*y7*y9/3._ark+y6*y7*y9)*y3**2*y2
    dF(2,613) = (3._ark/4._ark*y5**2*y9+sqrt(3._ark)*y5*y6*y9/2._ark+(y9/4._ark+ &
      y7)*y6**2)*y4*y1**2+(-3._ark/4._ark*y5**2*y9-sqrt(3._ark)*y5*y6*y9/2._ark+(-y9/4._ark- &
      y7)*y6**2)*y4**2*y1+(-3._ark/4._ark*y5**2*y9-sqrt(3._ark)*y5*y6*y9/2._ark+(-y9/4._ark+ &
      y7)*y6**2)*y3*y2**2+(3._ark/4._ark*y5**2*y9+sqrt(3._ark)*y5*y6*y9/2._ark+(y9/4._ark- &
      y7)*y6**2)*y3**2*y2
    dF(2,614) = ((-sqrt(3._ark)*y6*y7*y8/2._ark-y5*y7*y8/2._ark)*y2+(-y5*y8*y9/2._ark- &
      sqrt(3._ark)*y6*y8*y9/2._ark)*y3)*y1**2+((y5*y7*y8/2._ark+sqrt(3._ark)*y6*y7*y8/ &
      2._ark)*y2**2+(sqrt(3._ark)*y6*y8*y9/2._ark+y5*y8*y9/2._ark)*y3**2)*y1+(-y5*y8*y9/2._ark- &
      sqrt(3._ark)*y6*y8*y9/2._ark)*y4*y2**2+(sqrt(3._ark)*y6*y8*y9/2._ark+y5*y8*y9/ &
      2._ark)*y4**2*y2+(-sqrt(3._ark)*y6*y7*y8/2._ark-y5*y7*y8/2._ark)*y4*y3**2+(y5*y7*y8/ &
      2._ark+sqrt(3._ark)*y6*y7*y8/2._ark)*y4**2*y3
    dF(2,615) = (y5*y6**2-y5**3/3._ark)*y4*y1**2+(y5*y6**2-y5**3/3._ark)*y4**2*y1+ &
      (y5**3/3._ark-y5*y6**2)*y3*y2**2+(y5**3/3._ark-y5*y6**2)*y3**2*y2
    dF(2,616) = ((y7**2-2._ark*y9**2)*y5-sqrt(3._ark)*y6*y7**2)*y4*y1**2+((y7**2- &
      2._ark*y9**2)*y5-sqrt(3._ark)*y6*y7**2)*y4**2*y1+((-y7**2+2._ark*y9**2)*y5+ &
      sqrt(3._ark)*y6*y7**2)*y3*y2**2+((-y7**2+2._ark*y9**2)*y5+ &
      sqrt(3._ark)*y6*y7**2)*y3**2*y2
    dF(2,617) = ((y5*y8*y9-sqrt(3._ark)*y6*y8*y9)*y2-2._ark*y3*y5*y7*y8)*y1**2+ &
      ((y5*y8*y9-sqrt(3._ark)*y6*y8*y9)*y2**2-2._ark*y3**2*y5*y7*y8)*y1+ &
      2._ark*y2**2*y4*y5*y7*y8+2._ark*y2*y4**2*y5*y7*y8+(-y5*y8*y9+ &
      sqrt(3._ark)*y6*y8*y9)*y4*y3**2+(-y5*y8*y9+sqrt(3._ark)*y6*y8*y9)*y4**2*y3
    dF(2,618) = ((y5*y6*y7-sqrt(3._ark)*y6**2*y7/6._ark-sqrt(3._ark)*y5**2*y7/2._ark)*y2+ &
      (y5*y6*y9-sqrt(3._ark)*y5**2*y9/2._ark-sqrt(3._ark)*y6**2*y9/6._ark)*y3)*y1**2+ &
      ((y5*y6*y7-sqrt(3._ark)*y6**2*y7/6._ark-sqrt(3._ark)*y5**2*y7/2._ark)*y2**2+(y5*y6*y9- &
      sqrt(3._ark)*y5**2*y9/2._ark-sqrt(3._ark)*y6**2*y9/6._ark)*y3**2)*y1+(-y5*y6*y9+ &
      sqrt(3._ark)*y6**2*y9/6._ark+sqrt(3._ark)*y5**2*y9/2._ark)*y4*y2**2+(-y5*y6*y9+ &
      sqrt(3._ark)*y6**2*y9/6._ark+sqrt(3._ark)*y5**2*y9/2._ark)*y4**2*y2+ &
      (sqrt(3._ark)*y5**2*y7/2._ark+sqrt(3._ark)*y6**2*y7/6._ark-y5*y6*y7)*y4*y3**2+ &
      (sqrt(3._ark)*y5**2*y7/2._ark+sqrt(3._ark)*y6**2*y7/6._ark-y5*y6*y7)*y4**2*y3
    dF(2,619) = ((sqrt(3._ark)*y5**2*y8/2._ark-2._ark*y5*y6*y8+sqrt(3._ark)*y6**2*y8/ &
      2._ark)*y2+(sqrt(3._ark)*y5**2*y8-y5*y6*y8)*y3)*y1**2+((sqrt(3._ark)*y5**2*y8/2._ark- &
      2._ark*y5*y6*y8+sqrt(3._ark)*y6**2*y8/2._ark)*y2**2+(sqrt(3._ark)*y5**2*y8- &
      y5*y6*y8)*y3**2)*y1+(sqrt(3._ark)*y5**2*y8-y5*y6*y8)*y4*y2**2+ &
      (sqrt(3._ark)*y5**2*y8-y5*y6*y8)*y4**2*y2+(sqrt(3._ark)*y5**2*y8/2._ark-2._ark*y5*y6*y8+ &
      sqrt(3._ark)*y6**2*y8/2._ark)*y4*y3**2+(sqrt(3._ark)*y5**2*y8/2._ark-2._ark*y5*y6*y8+ &
      sqrt(3._ark)*y6**2*y8/2._ark)*y4**2*y3
    dF(2,620) = (sqrt(3._ark)*y5**2*y9/4._ark+(-y7-y9/2._ark)*y6*y5-sqrt(3._ark)*y6**2*y9/ &
      4._ark)*y4*y1**2+(-sqrt(3._ark)*y5**2*y9/4._ark+(y9/2._ark+y7)*y6*y5+ &
      sqrt(3._ark)*y6**2*y9/4._ark)*y4**2*y1+(-sqrt(3._ark)*y5**2*y9/4._ark+(-y7+y9/ &
      2._ark)*y6*y5+sqrt(3._ark)*y6**2*y9/4._ark)*y3*y2**2+(sqrt(3._ark)*y5**2*y9/4._ark+(y7-y9/ &
      2._ark)*y6*y5-sqrt(3._ark)*y6**2*y9/4._ark)*y3**2*y2
    dF(2,621) = ((sqrt(3._ark)*y6**2*y7/3._ark-y5*y6*y7)*y2+(-sqrt(3._ark)*y6**2*y9/6._ark+ &
      sqrt(3._ark)*y5**2*y9/2._ark)*y3)*y1**2+((sqrt(3._ark)*y6**2*y7/3._ark-y5*y6*y7)*y2**2+ &
      (-sqrt(3._ark)*y6**2*y9/6._ark+sqrt(3._ark)*y5**2*y9/2._ark)*y3**2)*y1+ &
      (sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y4*y2**2+ &
      (sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y4**2*y2+(- &
      sqrt(3._ark)*y6**2*y7/3._ark+y5*y6*y7)*y4*y3**2+(-sqrt(3._ark)*y6**2*y7/3._ark+ &
      y5*y6*y7)*y4**2*y3
    dF(2,622) = ((2._ark/3._ark*y5*y6**2-4._ark/9._ark*sqrt(3._ark)*y5**2*y6-2._ark/ &
      3._ark*y5**3)*y2+(sqrt(3._ark)*y5**2*y6/9._ark+y5*y6**2-y5**3/3._ark-sqrt(3._ark)*y6**3/ &
      3._ark)*y3)*y1**2+((4._ark/9._ark*sqrt(3._ark)*y5**2*y6+2._ark/3._ark*y5**3-2._ark/ &
      3._ark*y5*y6**2)*y2**2+(y5**3/3._ark-sqrt(3._ark)*y5**2*y6/9._ark-y5*y6**2+ &
      sqrt(3._ark)*y6**3/3._ark)*y3**2)*y1+(y5**3/3._ark-sqrt(3._ark)*y5**2*y6/9._ark-y5*y6**2+ &
      sqrt(3._ark)*y6**3/3._ark)*y4*y2**2+(sqrt(3._ark)*y5**2*y6/9._ark+y5*y6**2-y5**3/3._ark- &
      sqrt(3._ark)*y6**3/3._ark)*y4**2*y2+(4._ark/9._ark*sqrt(3._ark)*y5**2*y6+2._ark/3._ark*y5**3- &
      2._ark/3._ark*y5*y6**2)*y4*y3**2+(2._ark/3._ark*y5*y6**2-4._ark/9._ark*sqrt(3._ark)*y5**2*y6- &
      2._ark/3._ark*y5**3)*y4**2*y3
    dF(2,623) = ((-2._ark*y5**2*y8+2._ark*sqrt(3._ark)*y5*y6*y8)*y2+(-2._ark*y5**2*y8+ &
      2._ark*sqrt(3._ark)*y5*y6*y8)*y3)*y1**2+((-2._ark*y5**2*y8+ &
      2._ark*sqrt(3._ark)*y5*y6*y8)*y2**2+(-2._ark*y5**2*y8+ &
      2._ark*sqrt(3._ark)*y5*y6*y8)*y3**2)*y1+(-2._ark*y5**2*y8+ &
      2._ark*sqrt(3._ark)*y5*y6*y8)*y4*y2**2+(-2._ark*y5**2*y8+ &
      2._ark*sqrt(3._ark)*y5*y6*y8)*y4**2*y2+(-2._ark*y5**2*y8+ &
      2._ark*sqrt(3._ark)*y5*y6*y8)*y4*y3**2+(-2._ark*y5**2*y8+ &
      2._ark*sqrt(3._ark)*y5*y6*y8)*y4**2*y3
    dF(2,624) = ((y9/4._ark+y7)*y5**2-sqrt(3._ark)*y5*y6*y9/2._ark+3._ark/ &
      4._ark*y6**2*y9)*y4*y1**2+((-y9/4._ark-y7)*y5**2+sqrt(3._ark)*y5*y6*y9/2._ark-3._ark/ &
      4._ark*y6**2*y9)*y4**2*y1+((-y9/4._ark+y7)*y5**2+sqrt(3._ark)*y5*y6*y9/2._ark-3._ark/ &
      4._ark*y6**2*y9)*y3*y2**2+((y9/4._ark-y7)*y5**2-sqrt(3._ark)*y5*y6*y9/2._ark+3._ark/ &
      4._ark*y6**2*y9)*y3**2*y2
    dF(2,625) = (y2**2*y7*y9+y3**2*y7*y9)*y1**2+y3**2*y4**2*y7*y9+y2**2*y4**2*y7*y9
    dF(2,626) = (sqrt(3._ark)*y5**2/3._ark-y5*y6)*y4**2*y1**2+(y5*y6-sqrt(3._ark)*y5**2/ &
      3._ark)*y3**2*y2**2
    dF(2,627) = (y8*y9**2+y7**2*y8)*y4*y1**2+(y8*y9**2+y7**2*y8)*y4**2*y1+(y8*y9**2+ &
      y7**2*y8)*y3*y2**2+(y8*y9**2+y7**2*y8)*y3**2*y2
    dF(2,628) = ((sqrt(3._ark)*y5*y6*y9+y5**2*y9)*y2+(-y5**2*y7/2._ark+3._ark/ &
      2._ark*y6**2*y7)*y3)*y1**2+((-y5**2*y9-sqrt(3._ark)*y5*y6*y9)*y2**2+(-3._ark/ &
      2._ark*y6**2*y7+y5**2*y7/2._ark)*y3**2)*y1+(-y5**2*y7/2._ark+3._ark/ &
      2._ark*y6**2*y7)*y4*y2**2+(-3._ark/2._ark*y6**2*y7+y5**2*y7/2._ark)*y4**2*y2+ &
      (sqrt(3._ark)*y5*y6*y9+y5**2*y9)*y4*y3**2+(-y5**2*y9-sqrt(3._ark)*y5*y6*y9)*y4**2*y3
    dF(2,629) = (y5**2*y8+y6**2*y8)*y4*y1**2+(y5**2*y8+y6**2*y8)*y4**2*y1+(y5**2*y8+ &
      y6**2*y8)*y3*y2**2+(y5**2*y8+y6**2*y8)*y3**2*y2
    dF(2,630) = (y3*y4*y7*y8+y2*y4*y8*y9)*y1**2+(y2**2*y3*y8*y9+(y3**2*y7*y8- &
      y4**2*y7*y8)*y2-y3*y4**2*y8*y9)*y1-y2**2*y3*y4*y7*y8-y2*y3**2*y4*y8*y9
    dF(2,631) = ((y6*y7/2._ark-sqrt(3._ark)*y5*y7/2._ark)*y4*y2+(-y6*y9/2._ark+ &
      sqrt(3._ark)*y5*y9/2._ark)*y4*y3)*y1**2+((y6*y7/2._ark-sqrt(3._ark)*y5*y7/ &
      2._ark)*y3*y2**2+((-y6*y9/2._ark+sqrt(3._ark)*y5*y9/2._ark)*y3**2+(y6*y9/2._ark- &
      sqrt(3._ark)*y5*y9/2._ark)*y4**2)*y2+(-y6*y7/2._ark+sqrt(3._ark)*y5*y7/ &
      2._ark)*y4**2*y3)*y1+(y6*y9/2._ark-sqrt(3._ark)*y5*y9/2._ark)*y4*y3*y2**2+(-y6*y7/2._ark+ &
      sqrt(3._ark)*y5*y7/2._ark)*y4*y3**2*y2
    dF(2,632) = ((-2._ark*y9+y7)*y5-sqrt(3._ark)*y6*y7)*y3*y2*y1**2+(((y7+2._ark*y9)*y5- &
      sqrt(3._ark)*y6*y7)*y4*y2**2+((-y7-2._ark*y9)*y5+sqrt(3._ark)*y6*y7)*y4*y3**2)*y1+ &
      ((2._ark*y9-y7)*y5+sqrt(3._ark)*y6*y7)*y4**2*y3*y2
    dF(2,633) = ((y5*y8+sqrt(3._ark)*y6*y8)*y4*y2+(y5*y8+ &
      sqrt(3._ark)*y6*y8)*y4*y3)*y1**2+((y5*y8+sqrt(3._ark)*y6*y8)*y3*y2**2+((y5*y8+ &
      sqrt(3._ark)*y6*y8)*y3**2+(y5*y8+sqrt(3._ark)*y6*y8)*y4**2)*y2+(y5*y8+ &
      sqrt(3._ark)*y6*y8)*y4**2*y3)*y1+(y5*y8+sqrt(3._ark)*y6*y8)*y4*y3*y2**2+(y5*y8+ &
      sqrt(3._ark)*y6*y8)*y4*y3**2*y2
    dF(2,634) = (y5*y8+sqrt(3._ark)*y6*y8)*y3*y2*y1**2+((y5*y8+ &
      sqrt(3._ark)*y6*y8)*y4*y2**2+(y5*y8+sqrt(3._ark)*y6*y8)*y4*y3**2)*y1+(y5*y8+ &
      sqrt(3._ark)*y6*y8)*y4**2*y3*y2
    dF(2,635) = ((-sqrt(3._ark)*y6*y7/2._ark-y5*y7/2._ark)*y4*y2+(-y5*y9/2._ark- &
      sqrt(3._ark)*y6*y9/2._ark)*y4*y3)*y1**2+((-sqrt(3._ark)*y6*y7/2._ark-y5*y7/ &
      2._ark)*y3*y2**2+((-y5*y9/2._ark-sqrt(3._ark)*y6*y9/2._ark)*y3**2+(sqrt(3._ark)*y6*y9/ &
      2._ark+y5*y9/2._ark)*y4**2)*y2+(sqrt(3._ark)*y6*y7/2._ark+y5*y7/2._ark)*y4**2*y3)*y1+ &
      (sqrt(3._ark)*y6*y9/2._ark+y5*y9/2._ark)*y4*y3*y2**2+(sqrt(3._ark)*y6*y7/2._ark+y5*y7/ &
      2._ark)*y4*y3**2*y2
    dF(2,636) = ((-sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y4*y2+(-y5*y6- &
      sqrt(3._ark)*y5**2/3._ark)*y4*y3)*y1**2+((-sqrt(3._ark)*y5**2/6._ark+sqrt(3._ark)*y6**2/ &
      2._ark)*y3*y2**2+((y5*y6+sqrt(3._ark)*y5**2/3._ark)*y3**2+(-y5*y6-sqrt(3._ark)*y5**2/ &
      3._ark)*y4**2)*y2+(-sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y4**2*y3)*y1+ &
      (y5*y6+sqrt(3._ark)*y5**2/3._ark)*y4*y3*y2**2+(-sqrt(3._ark)*y5**2/6._ark+ &
      sqrt(3._ark)*y6**2/2._ark)*y4*y3**2*y2
    dF(2,637) = (y6**2+y5**2)*y3*y2*y1**2+((-y6**2-y5**2)*y4*y2**2+(-y6**2- &
      y5**2)*y4*y3**2)*y1+(y6**2+y5**2)*y4**2*y3*y2
    dF(2,638) = (((-sqrt(3._ark)*y5+y6)*y3+(-y6+sqrt(3._ark)*y5)*y4)*y2**2+(-y6+ &
      sqrt(3._ark)*y5)*y3**2*y2+(-sqrt(3._ark)*y5+y6)*y4*y3**2)*y1**2+((-sqrt(3._ark)*y5+ &
      y6)*y4**2*y2**2+(-y6+sqrt(3._ark)*y5)*y4**2*y3**2)*y1+(-y6+ &
      sqrt(3._ark)*y5)*y4**2*y3*y2**2+(-sqrt(3._ark)*y5+y6)*y4**2*y3**2*y2
    dF(2,639) = (y3*y4**2+y2*y4**2)*y1**3+(y2*y4**3+y3*y4**3)*y1**2+(-y2**2*y3**3- &
      y2**3*y3**2)*y1-y2**3*y3**2*y4-y2**2*y3**3*y4
    dF(2,640) = (-y2*y3**2-y2**2*y3)*y1**3+(y3**3*y4+y2**3*y4)*y1**2+(y2**3*y4**2+ &
      y3**3*y4**2)*y1-y2*y3**2*y4**3-y2**2*y3*y4**3
    dF(2,641) = (-2._ark*y2*y4*y6*y8+(-y6*y8-sqrt(3._ark)*y5*y8)*y4*y3)*y1**2+(- &
      2._ark*y2**2*y3*y6*y8+((-y6*y8-sqrt(3._ark)*y5*y8)*y3**2+(-y6*y8- &
      sqrt(3._ark)*y5*y8)*y4**2)*y2-2._ark*y3*y4**2*y6*y8)*y1+(-y6*y8- &
      sqrt(3._ark)*y5*y8)*y4*y3*y2**2-2._ark*y2*y3**2*y4*y6*y8
    dF(2,642) = ((y6**2+y5**2)*y4*y2+(y6**2+y5**2)*y4*y3)*y1**2+((-y6**2- &
      y5**2)*y3*y2**2+((-y6**2-y5**2)*y3**2+(y6**2+y5**2)*y4**2)*y2+(y6**2+ &
      y5**2)*y4**2*y3)*y1+(-y6**2-y5**2)*y4*y3*y2**2+(-y6**2-y5**2)*y4*y3**2*y2
    dF(2,643) = ((y5*y9-sqrt(3._ark)*y6*y9)*y4*y2-2._ark*y3*y4*y5*y7)*y1**2+ &
      ((sqrt(3._ark)*y6*y9-y5*y9)*y3*y2**2+(2._ark*y4**2*y5*y7+2._ark*y3**2*y5*y7)*y2+ &
      (sqrt(3._ark)*y6*y9-y5*y9)*y4**2*y3)*y1-2._ark*y2**2*y3*y4*y5*y7+(y5*y9- &
      sqrt(3._ark)*y6*y9)*y4*y3**2*y2
    dF(2,644) = -y1*y2**2*y3**2*y4+y1**2*y2*y3*y4**2
    dF(2,645) = y1**2*y4**2*y8**2-y2**2*y3**2*y8**2
    dF(2,646) = ((-y3**2+y4**2)*y2**2+y3**2*y4**2)*y1**2-y2**2*y3**2*y4**2
    dF(2,647) = (-sqrt(3._ark)*y5*y8**2/3._ark-y6*y8**2)*y1**3+(sqrt(3._ark)*y5*y8**2/ &
      3._ark+y6*y8**2)*y2**3+(sqrt(3._ark)*y5*y8**2/3._ark+y6*y8**2)*y3**3+(- &
      sqrt(3._ark)*y5*y8**2/3._ark-y6*y8**2)*y4**3
    dF(2,648) = (-sqrt(3._ark)*y6*y7*y9/2._ark-y5*y7*y9/2._ark)*y1**3+(- &
      sqrt(3._ark)*y6*y7*y9/2._ark-y5*y7*y9/2._ark)*y2**3+(-sqrt(3._ark)*y6*y7*y9/2._ark- &
      y5*y7*y9/2._ark)*y3**3+(-sqrt(3._ark)*y6*y7*y9/2._ark-y5*y7*y9/2._ark)*y4**3
    dF(2,649) = (-sqrt(3._ark)*y6**2*y8/3._ark-y5*y6*y8)*y1**3+(-sqrt(3._ark)*y6**2*y8/ &
      3._ark-y5*y6*y8)*y2**3+(-sqrt(3._ark)*y6**2*y8/3._ark-y5*y6*y8)*y3**3+(- &
      sqrt(3._ark)*y6**2*y8/3._ark-y5*y6*y8)*y4**3
    dF(2,650) = ((sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y7)*y5**2+(-y7-2._ark*y9)*y6*y5+ &
      sqrt(3._ark)*y6**2*y9/2._ark)*y1**3+((-sqrt(3._ark)*y9/2._ark+sqrt(3._ark)*y7)*y5**2+ &
      (2._ark*y9-y7)*y6*y5-sqrt(3._ark)*y6**2*y9/2._ark)*y2**3+((sqrt(3._ark)*y9/2._ark- &
      sqrt(3._ark)*y7)*y5**2+(-2._ark*y9+y7)*y6*y5+sqrt(3._ark)*y6**2*y9/2._ark)*y3**3+((- &
      sqrt(3._ark)*y7-sqrt(3._ark)*y9/2._ark)*y5**2+(y7+2._ark*y9)*y6*y5-sqrt(3._ark)*y6**2*y9/ &
      2._ark)*y4**3
    dF(2,651) = (sqrt(3._ark)*y5**3-5._ark/3._ark*sqrt(3._ark)*y5*y6**2+y5**2*y6+ &
      y6**3)*y1**3+(-sqrt(3._ark)*y5**3-y6**3+5._ark/3._ark*sqrt(3._ark)*y5*y6**2- &
      y5**2*y6)*y2**3+(-sqrt(3._ark)*y5**3-y6**3+5._ark/3._ark*sqrt(3._ark)*y5*y6**2- &
      y5**2*y6)*y3**3+(sqrt(3._ark)*y5**3-5._ark/3._ark*sqrt(3._ark)*y5*y6**2+y5**2*y6+ &
      y6**3)*y4**3
    dF(2,652) = y1**3*y4*y7*y9+y2*y3**3*y7*y9+y2**3*y3*y7*y9+y1*y4**3*y7*y9
    dF(2,653) = (y3*y8*y9+y2*y7*y8)*y1**3+(-y3**3*y8*y9-y2**3*y7*y8)*y1+ &
      y3**3*y4*y7*y8+y2**3*y4*y8*y9-y3*y4**3*y7*y8-y2*y4**3*y8*y9
    dF(2,654) = (y8*y9+y7*y8)*y4*y1**3+(-y7*y8-y8*y9)*y4**3*y1+(y8*y9- &
      y7*y8)*y3*y2**3+(y7*y8-y8*y9)*y3**3*y2
    dF(2,655) = (y2*y7**2+y3*y9**2)*y1**3+(-y3**3*y9**2-y2**3*y7**2)*y1+ &
      y3*y4**3*y7**2+y2*y4**3*y9**2-y2**3*y4*y9**2-y3**3*y4*y7**2
    dF(2,656) = (y9**2+y7**2)*y4*y1**3+(y9**2+y7**2)*y4**3*y1+(-y9**2- &
      y7**2)*y3*y2**3+(-y9**2-y7**2)*y3**3*y2
    dF(2,657) = (-sqrt(3._ark)*y5*y9/2._ark+(-y7-y9/2._ark)*y6)*y4*y1**3+ &
      (sqrt(3._ark)*y5*y9/2._ark+(y9/2._ark+y7)*y6)*y4**3*y1+(sqrt(3._ark)*y5*y9/2._ark+(-y7+y9/ &
      2._ark)*y6)*y3*y2**3+(-sqrt(3._ark)*y5*y9/2._ark+(y7-y9/2._ark)*y6)*y3**3*y2
    dF(2,658) = ((sqrt(3._ark)*y5*y8/2._ark+y6*y8/2._ark)*y2+y3*y6*y8)*y1**3+ &
      ((sqrt(3._ark)*y5*y8/2._ark+y6*y8/2._ark)*y2**3+y3**3*y6*y8)*y1+y2**3*y4*y6*y8+ &
      y2*y4**3*y6*y8+(sqrt(3._ark)*y5*y8/2._ark+y6*y8/2._ark)*y4*y3**3+(sqrt(3._ark)*y5*y8/ &
      2._ark+y6*y8/2._ark)*y4**3*y3
    dF(2,659) = ((-y6*y7/2._ark+sqrt(3._ark)*y5*y7/2._ark)*y2+(y6*y9/2._ark- &
      sqrt(3._ark)*y5*y9/2._ark)*y3)*y1**3+((-y6*y7/2._ark+sqrt(3._ark)*y5*y7/2._ark)*y2**3+ &
      (y6*y9/2._ark-sqrt(3._ark)*y5*y9/2._ark)*y3**3)*y1+(-y6*y9/2._ark+sqrt(3._ark)*y5*y9/ &
      2._ark)*y4*y2**3+(-y6*y9/2._ark+sqrt(3._ark)*y5*y9/2._ark)*y4**3*y2+(y6*y7/2._ark- &
      sqrt(3._ark)*y5*y7/2._ark)*y4*y3**3+(y6*y7/2._ark-sqrt(3._ark)*y5*y7/2._ark)*y4**3*y3
    dF(2,660) = (-2._ark*y2*y6*y9+(-y6*y7-sqrt(3._ark)*y5*y7)*y3)*y1**3+ &
      (2._ark*y2**3*y6*y9+(sqrt(3._ark)*y5*y7+y6*y7)*y3**3)*y1+(-y6*y7- &
      sqrt(3._ark)*y5*y7)*y4*y2**3+(sqrt(3._ark)*y5*y7+y6*y7)*y4**3*y2- &
      2._ark*y3**3*y4*y6*y9+2._ark*y3*y4**3*y6*y9
    dF(2,661) = ((y7-y9/2._ark)*y5+sqrt(3._ark)*y6*y9/2._ark)*y4*y1**3+((-y7+y9/2._ark)*y5- &
      sqrt(3._ark)*y6*y9/2._ark)*y4**3*y1+((y9/2._ark+y7)*y5-sqrt(3._ark)*y6*y9/ &
      2._ark)*y3*y2**3+((-y7-y9/2._ark)*y5+sqrt(3._ark)*y6*y9/2._ark)*y3**3*y2
    dF(2,662) = ((sqrt(3._ark)*y6*y8/2._ark-y5*y8/2._ark)*y2+y3*y5*y8)*y1**3+ &
      ((sqrt(3._ark)*y6*y8/2._ark-y5*y8/2._ark)*y2**3+y3**3*y5*y8)*y1+y2**3*y4*y5*y8+ &
      y2*y4**3*y5*y8+(sqrt(3._ark)*y6*y8/2._ark-y5*y8/2._ark)*y4*y3**3+(sqrt(3._ark)*y6*y8/ &
      2._ark-y5*y8/2._ark)*y4**3*y3
    dF(2,663) = ((-sqrt(3._ark)*y6*y7/2._ark-y5*y7/2._ark)*y2+(-y5*y9/2._ark- &
      sqrt(3._ark)*y6*y9/2._ark)*y3)*y1**3+((-sqrt(3._ark)*y6*y7/2._ark-y5*y7/2._ark)*y2**3+(- &
      y5*y9/2._ark-sqrt(3._ark)*y6*y9/2._ark)*y3**3)*y1+(sqrt(3._ark)*y6*y9/2._ark+y5*y9/ &
      2._ark)*y4*y2**3+(sqrt(3._ark)*y6*y9/2._ark+y5*y9/2._ark)*y4**3*y2+(sqrt(3._ark)*y6*y7/ &
      2._ark+y5*y7/2._ark)*y4*y3**3+(sqrt(3._ark)*y6*y7/2._ark+y5*y7/2._ark)*y4**3*y3
    dF(2,664) = (-sqrt(3._ark)*y6**2/2._ark-sqrt(3._ark)*y5**2/6._ark-y5*y6)*y4*y1**3+(- &
      sqrt(3._ark)*y6**2/2._ark-sqrt(3._ark)*y5**2/6._ark-y5*y6)*y4**3*y1+(y5*y6+ &
      sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y3*y2**3+(y5*y6+sqrt(3._ark)*y6**2/ &
      2._ark+sqrt(3._ark)*y5**2/6._ark)*y3**3*y2
    dF(2,665) = (2._ark/3._ark*sqrt(3._ark)*y2*y6**2+(y5*y6+sqrt(3._ark)*y5**2/2._ark+ &
      sqrt(3._ark)*y6**2/6._ark)*y3)*y1**3+(-2._ark/3._ark*sqrt(3._ark)*y2**3*y6**2+(- &
      sqrt(3._ark)*y6**2/6._ark-y5*y6-sqrt(3._ark)*y5**2/2._ark)*y3**3)*y1+(-sqrt(3._ark)*y6**2/ &
      6._ark-y5*y6-sqrt(3._ark)*y5**2/2._ark)*y4*y2**3+(y5*y6+sqrt(3._ark)*y5**2/2._ark+ &
      sqrt(3._ark)*y6**2/6._ark)*y4**3*y2-2._ark/3._ark*sqrt(3._ark)*y3**3*y4*y6**2+2._ark/ &
      3._ark*sqrt(3._ark)*y3*y4**3*y6**2
    dF(2,666) = (y9+y7)*y3*y2*y1**3+((-y9+y7)*y4*y2**3+(-y7+y9)*y4*y3**3)*y1+(-y9- &
      y7)*y4**3*y3*y2
    dF(2,667) = y1*y4**3*y8**2+y1**3*y4*y8**2-y2*y3**3*y8**2-y2**3*y3*y8**2
    dF(2,668) = y1**3*y2*y3*y4+(-y2**3*y3*y4+(y3*y4**3-y3**3*y4)*y2)*y1
    dF(2,669) = (y8*y9+y7*y8)*y1**4+(y8*y9-y7*y8)*y2**4+(y7*y8-y8*y9)*y3**4+(-y7*y8- &
      y8*y9)*y4**4
    dF(2,670) = y4**4*y7*y9+y3**4*y7*y9+y2**4*y7*y9+y1**4*y7*y9
    dF(2,671) = y4**4*y8**2+y1**4*y8**2-y2**4*y8**2-y3**4*y8**2
    dF(2,672) = (-y6*y8-sqrt(3._ark)*y5*y8/3._ark)*y1**4+(-y6*y8-sqrt(3._ark)*y5*y8/ &
      3._ark)*y2**4+(-y6*y8-sqrt(3._ark)*y5*y8/3._ark)*y3**4+(-y6*y8-sqrt(3._ark)*y5*y8/ &
      3._ark)*y4**4
    dF(2,673) = (sqrt(3._ark)*y5*y9/2._ark+(y9/2._ark+y7)*y6)*y1**4+(-sqrt(3._ark)*y5*y9/ &
      2._ark+(y7-y9/2._ark)*y6)*y2**4+(sqrt(3._ark)*y5*y9/2._ark+(-y7+y9/2._ark)*y6)*y3**4+(- &
      sqrt(3._ark)*y5*y9/2._ark+(-y7-y9/2._ark)*y6)*y4**4
    dF(2,674) = (y3*y7+y2*y9)*y1**4+(-y3**4*y7-y2**4*y9)*y1+y3**4*y4*y9-y2*y4**4*y7+ &
      y2**4*y4*y7-y3*y4**4*y9
    dF(2,675) = y1**4*y4*y8+y1*y4**4*y8+y2*y3**4*y8+y2**4*y3*y8
    dF(2,676) = (y2*y8+y3*y8)*y1**4+(y2**4*y8+y3**4*y8)*y1+y3*y4**4*y8+y2*y4**4*y8+ &
      y2**4*y4*y8+y3**4*y4*y8
    dF(2,677) = (y3*y9+y2*y7)*y1**4+(y3**4*y9+y2**4*y7)*y1-y3**4*y4*y7-y2**4*y4*y9- &
      y3*y4**4*y7-y2*y4**4*y9
    dF(2,678) = (y9+y7)*y4*y1**4+(-y9-y7)*y4**4*y1+(-y9+y7)*y3*y2**4+(-y7+ &
      y9)*y3**4*y2
    dF(2,679) = (-sqrt(3._ark)*y6/2._ark-y5/2._ark)*y4*y1**4+(-sqrt(3._ark)*y6/2._ark-y5/ &
      2._ark)*y4**4*y1+(sqrt(3._ark)*y6/2._ark+y5/2._ark)*y3*y2**4+(sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y3**4*y2
    dF(2,680) = y1**4*y4**2+y1**2*y4**4-y2**4*y3**2-y2**2*y3**4
    dF(2,681) = y4**6+y1**6-y3**6-y2**6
    dF(3,1) = 0._ark
    dF(3,2) = y9
    dF(3,3) = -y4+y2+y1-y3
    dF(3,4) = y7*y8
    dF(3,5) = -2._ark*y5*y9
    dF(3,6) = y1**2+y2**2-y4**2-y3**2
    dF(3,7) = y3*y9+y1*y9+y2*y9+y4*y9
    dF(3,8) = (y8+y7)*y1+(-y7-y8)*y2+(y8-y7)*y3+(-y8+y7)*y4
    dF(3,9) = -y4*y5+y2*y5-y3*y5+y1*y5
    dF(3,10) = y1*y2-y3*y4
    dF(3,11) = y8**2*y9+y7**2*y9
    dF(3,12) = y9**3
    dF(3,13) = y5*y7*y8
    dF(3,14) = sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark
    dF(3,15) = y6**2*y9+y5**2*y9
    dF(3,16) = (-y7*y9-y8*y9)*y1+(y7*y9+y8*y9)*y2+(-y7*y9+y8*y9)*y3+(-y8*y9+ &
      y7*y9)*y4
    dF(3,17) = y1*y7*y8+y4*y7*y8+y3*y7*y8+y2*y7*y8
    dF(3,18) = y3*y9**2+y4*y9**2-y2*y9**2-y1*y9**2
    dF(3,19) = (-y8**2-y7**2)*y1+(-y8**2-y7**2)*y2+(y8**2+y7**2)*y3+(y8**2+y7**2)*y4
    dF(3,20) = 2._ark/3._ark*sqrt(3._ark)*y4*y5*y9+2._ark/3._ark*sqrt(3._ark)*y1*y5*y9+2._ark/ &
      3._ark*sqrt(3._ark)*y2*y5*y9+2._ark/3._ark*sqrt(3._ark)*y3*y5*y9
    dF(3,21) = (-y8+y7)*y6*y1+(y8-y7)*y6*y2+(-y7-y8)*y6*y3+(y8+y7)*y6*y4
    dF(3,22) = (y8+y7)*y5*y1+(-y7-y8)*y5*y2+(y8-y7)*y5*y3+(-y8+y7)*y5*y4
    dF(3,23) = (-sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y1+(-sqrt(3._ark)*y6**2/ &
      2._ark+sqrt(3._ark)*y5**2/6._ark)*y2+(-sqrt(3._ark)*y5**2/6._ark+sqrt(3._ark)*y6**2/ &
      2._ark)*y3+(-sqrt(3._ark)*y5**2/6._ark+sqrt(3._ark)*y6**2/2._ark)*y4
    dF(3,24) = -y1**2*y5-y2**2*y5+y3**2*y5+y4**2*y5
    dF(3,25) = -y2**3+y3**3-y1**3+y4**3
    dF(3,26) = y3*y4*y9+y1*y2*y9
    dF(3,27) = (y4*y9+y3*y9)*y1+(y4*y9+y3*y9)*y2
    dF(3,28) = (-y4*y7-y3*y8)*y1+(y4*y8+y3*y7)*y2
    dF(3,29) = y3*y4*y5-y1*y2*y5
    dF(3,30) = (y4+y3)*y1**2+(-y4**2-y3**2)*y1+(y4+y3)*y2**2+(-y4**2-y3**2)*y2
    dF(3,31) = y3*y4**2-y1**2*y2+y3**2*y4-y1*y2**2
    dF(3,32) = (-y7-y8)*y1**2+(y8+y7)*y2**2+(-y8+y7)*y3**2+(y8-y7)*y4**2
    dF(3,33) = y3**2*y9+y4**2*y9+y2**2*y9+y1**2*y9
    dF(3,34) = (-y6**2-y5**2)*y1+(-y6**2-y5**2)*y2+(y6**2+y5**2)*y3+(y6**2+y5**2)*y4
    dF(3,35) = ((-y4-y3)*y2+y3*y4)*y1+y2*y3*y4
    dF(3,36) = y7*y8*y9**2
    dF(3,37) = y7*y8**3+y7**3*y8
    dF(3,38) = (-2._ark*y7**2*y9+2._ark*y8**2*y9)*y6
    dF(3,39) = (y8**2*y9+y7**2*y9)*y5+(sqrt(3._ark)*y7**2*y9-sqrt(3._ark)*y8**2*y9)*y6
    dF(3,40) = -2._ark*y5*y9**3
    dF(3,41) = -2._ark/3._ark*sqrt(3._ark)*y6**2*y7*y8
    dF(3,42) = y5**2*y7*y8+y6**2*y7*y8
    dF(3,43) = 8._ark/3._ark*sqrt(3._ark)*y5*y6**2*y9
    dF(3,44) = y5**3*y9-3._ark*y9*y6**2*y5
    dF(3,45) = y4*y9**3+y3*y9**3+y2*y9**3+y1*y9**3
    dF(3,46) = ((-sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y7/3._ark)*y5**2+(y8-y7)*y6*y5)*y1+ &
      ((sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/3._ark)*y5**2+(-y8+y7)*y6*y5)*y2+((- &
      sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/3._ark)*y5**2+(y8+y7)*y6*y5)*y3+((sqrt(3._ark)*y8/ &
      3._ark-sqrt(3._ark)*y7/3._ark)*y5**2+(-y7-y8)*y6*y5)*y4
    dF(3,47) = (y8-y7)*y6*y1**2+(-y8+y7)*y6*y2**2+(y8+y7)*y6*y3**2+(-y7-y8)*y6*y4**2
    dF(3,48) = (-y8**2-y7**2)*y1**2+(-y8**2-y7**2)*y2**2+(y8**2+y7**2)*y3**2+(y8**2+ &
      y7**2)*y4**2
    dF(3,49) = y2**3*y5-y4**3*y5-y3**3*y5+y1**3*y5
    dF(3,50) = -y2**4+y3**4+y4**4-y1**4
    dF(3,51) = (y8**2*y9+y7**2*y9)*y1+(y8**2*y9+y7**2*y9)*y2+(y8**2*y9+y7**2*y9)*y3+ &
      (y8**2*y9+y7**2*y9)*y4
    dF(3,52) = (y7*y8**2+y7**2*y8)*y1+(-y7**2*y8-y7*y8**2)*y2+(y7**2*y8- &
      y7*y8**2)*y3+(y7*y8**2-y7**2*y8)*y4
    dF(3,53) = (y8*y9**2+y7*y9**2)*y1+(-y7*y9**2-y8*y9**2)*y2+(y8*y9**2- &
      y7*y9**2)*y3+(y7*y9**2-y8*y9**2)*y4
    dF(3,54) = (y7**3+y8**3)*y1+(-y7**3-y8**3)*y2+(y8**3-y7**3)*y3+(y7**3-y8**3)*y4
    dF(3,55) = y1*y7*y8*y9+y2*y7*y8*y9-y4*y7*y8*y9-y3*y7*y8*y9
    dF(3,56) = ((-sqrt(3._ark)*y8**2/2._ark-sqrt(3._ark)*y7**2/2._ark)*y5+(y8**2/2._ark-y7**2/ &
      2._ark)*y6)*y1+((-sqrt(3._ark)*y8**2/2._ark-sqrt(3._ark)*y7**2/2._ark)*y5+(y8**2/2._ark- &
      y7**2/2._ark)*y6)*y2+((sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y7**2/2._ark)*y5+(-y8**2/ &
      2._ark+y7**2/2._ark)*y6)*y3+((sqrt(3._ark)*y8**2/2._ark+sqrt(3._ark)*y7**2/2._ark)*y5+(- &
      y8**2/2._ark+y7**2/2._ark)*y6)*y4
    dF(3,57) = 2._ark/3._ark*sqrt(3._ark)*y4*y5*y7*y8+2._ark/3._ark*sqrt(3._ark)*y3*y5*y7*y8+ &
      2._ark/3._ark*sqrt(3._ark)*y2*y5*y7*y8+2._ark/3._ark*sqrt(3._ark)*y1*y5*y7*y8
    dF(3,58) = ((-sqrt(3._ark)*y8*y9-sqrt(3._ark)*y7*y9)*y5+(-y7*y9+y8*y9)*y6)*y1+ &
      ((sqrt(3._ark)*y7*y9+sqrt(3._ark)*y8*y9)*y5+(-y8*y9+y7*y9)*y6)*y2+ &
      ((sqrt(3._ark)*y8*y9-sqrt(3._ark)*y7*y9)*y5+(-y7*y9-y8*y9)*y6)*y3+ &
      ((sqrt(3._ark)*y7*y9-sqrt(3._ark)*y8*y9)*y5+(y7*y9+y8*y9)*y6)*y4
    dF(3,59) = ((2._ark/3._ark*y8+2._ark/3._ark*y7)*y5**2+(2._ark/3._ark*sqrt(3._ark)*y8-2._ark/ &
      3._ark*sqrt(3._ark)*y7)*y6*y5)*y1+((-2._ark/3._ark*y8-2._ark/3._ark*y7)*y5**2+(-2._ark/ &
      3._ark*sqrt(3._ark)*y8+2._ark/3._ark*sqrt(3._ark)*y7)*y6*y5)*y2+((2._ark/3._ark*y8-2._ark/ &
      3._ark*y7)*y5**2+(2._ark/3._ark*sqrt(3._ark)*y8+2._ark/3._ark*sqrt(3._ark)*y7)*y6*y5)*y3+((- &
      2._ark/3._ark*y8+2._ark/3._ark*y7)*y5**2+(-2._ark/3._ark*sqrt(3._ark)*y8-2._ark/ &
      3._ark*sqrt(3._ark)*y7)*y6*y5)*y4
    dF(3,60) = ((y8/3._ark+y7/3._ark)*y5**2+(-2._ark/3._ark*sqrt(3._ark)*y8+2._ark/ &
      3._ark*sqrt(3._ark)*y7)*y6*y5+(y8+y7)*y6**2)*y1+((-y8/3._ark-y7/3._ark)*y5**2+(2._ark/ &
      3._ark*sqrt(3._ark)*y8-2._ark/3._ark*sqrt(3._ark)*y7)*y6*y5+(-y7-y8)*y6**2)*y2+((y8/3._ark- &
      y7/3._ark)*y5**2+(-2._ark/3._ark*sqrt(3._ark)*y8-2._ark/3._ark*sqrt(3._ark)*y7)*y6*y5+(y8- &
      y7)*y6**2)*y3+((-y8/3._ark+y7/3._ark)*y5**2+(2._ark/3._ark*sqrt(3._ark)*y8+2._ark/ &
      3._ark*sqrt(3._ark)*y7)*y6*y5+(-y8+y7)*y6**2)*y4
    dF(3,61) = ((-y8**2/2._ark-y7**2/2._ark)*y5+(sqrt(3._ark)*y7**2/2._ark-sqrt(3._ark)*y8**2/ &
      2._ark)*y6)*y1+((-y8**2/2._ark-y7**2/2._ark)*y5+(sqrt(3._ark)*y7**2/2._ark- &
      sqrt(3._ark)*y8**2/2._ark)*y6)*y2+((y8**2/2._ark+y7**2/2._ark)*y5+(-sqrt(3._ark)*y7**2/ &
      2._ark+sqrt(3._ark)*y8**2/2._ark)*y6)*y3+((y8**2/2._ark+y7**2/2._ark)*y5+(- &
      sqrt(3._ark)*y7**2/2._ark+sqrt(3._ark)*y8**2/2._ark)*y6)*y4
    dF(3,62) = -y4*y5*y9**2+y1*y5*y9**2+y2*y5*y9**2-y3*y5*y9**2
    dF(3,63) = ((y7*y9+y8*y9)*y5+(sqrt(3._ark)*y7*y9-sqrt(3._ark)*y8*y9)*y6)*y1+((- &
      y7*y9-y8*y9)*y5+(sqrt(3._ark)*y8*y9-sqrt(3._ark)*y7*y9)*y6)*y2+((-y8*y9+y7*y9)*y5+ &
      (sqrt(3._ark)*y7*y9+sqrt(3._ark)*y8*y9)*y6)*y3+((-y7*y9+y8*y9)*y5+(- &
      sqrt(3._ark)*y8*y9-sqrt(3._ark)*y7*y9)*y6)*y4
    dF(3,64) = (sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y1+ &
      (sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y2+(sqrt(3._ark)*y6**2*y9/ &
      6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y3+(sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/ &
      2._ark)*y4
    dF(3,65) = (y5*y6**2-y5**3/3._ark)*y1+(y5*y6**2-y5**3/3._ark)*y2+(y5**3/3._ark- &
      y5*y6**2)*y3+(y5**3/3._ark-y5*y6**2)*y4
    dF(3,66) = (y6**2*y9+y5**2*y9)*y1+(y6**2*y9+y5**2*y9)*y2+(y6**2*y9+y5**2*y9)*y3+ &
      (y6**2*y9+y5**2*y9)*y4
    dF(3,67) = 8._ark/9._ark*sqrt(3._ark)*y3*y5**3+8._ark/9._ark*sqrt(3._ark)*y4*y5**3-8._ark/ &
      9._ark*sqrt(3._ark)*y2*y5**3-8._ark/9._ark*sqrt(3._ark)*y1*y5**3
    dF(3,68) = -2._ark/3._ark*sqrt(3._ark)*y3*y4*y5**2+2._ark/3._ark*sqrt(3._ark)*y1*y2*y5**2
    dF(3,69) = 2._ark/3._ark*sqrt(3._ark)*y3*y4**2*y5+2._ark/3._ark*sqrt(3._ark)*y3**2*y4*y5- &
      2._ark/3._ark*sqrt(3._ark)*y1**2*y2*y5-2._ark/3._ark*sqrt(3._ark)*y1*y2**2*y5
    dF(3,70) = y1*y2**3+y1**3*y2-y3**3*y4-y3*y4**3
    dF(3,71) = (-y4-y3)*y1**3+(y4**3+y3**3)*y1+(-y4-y3)*y2**3+(y4**3+y3**3)*y2
    dF(3,72) = y1*y2*y9**2-y3*y4*y9**2
    dF(3,73) = (y3*y7*y9+y4*y8*y9)*y1+(-y3*y8*y9-y4*y7*y9)*y2
    dF(3,74) = (y8**2+y7**2)*y2*y1+(-y8**2-y7**2)*y4*y3
    dF(3,75) = (y3*y7*y8+y4*y7*y8)*y1+(y3*y7*y8+y4*y7*y8)*y2
    dF(3,76) = y1*y2*y7*y8+y3*y4*y7*y8
    dF(3,77) = ((y6*y9/2._ark+sqrt(3._ark)*y5*y9/2._ark)*y3+(-y6*y9/2._ark+sqrt(3._ark)*y5*y9/ &
      2._ark)*y4)*y1+((-y6*y9/2._ark+sqrt(3._ark)*y5*y9/2._ark)*y3+(y6*y9/2._ark+ &
      sqrt(3._ark)*y5*y9/2._ark)*y4)*y2
    dF(3,78) = ((-y6*y8/2._ark-sqrt(3._ark)*y5*y8/2._ark)*y3+(y6*y7/2._ark-sqrt(3._ark)*y5*y7/ &
      2._ark)*y4)*y1+((-y6*y7/2._ark+sqrt(3._ark)*y5*y7/2._ark)*y3+(sqrt(3._ark)*y5*y8/2._ark+ &
      y6*y8/2._ark)*y4)*y2
    dF(3,79) = (y6**2+y5**2)*y2*y1+(-y6**2-y5**2)*y4*y3
    dF(3,80) = y1*y2*y5*y9+y3*y4*y5*y9
    dF(3,81) = ((sqrt(3._ark)*y6*y9/2._ark-y5*y9/2._ark)*y3+(-y5*y9/2._ark-sqrt(3._ark)*y6*y9/ &
      2._ark)*y4)*y1+((-y5*y9/2._ark-sqrt(3._ark)*y6*y9/2._ark)*y3+(sqrt(3._ark)*y6*y9/2._ark- &
      y5*y9/2._ark)*y4)*y2
    dF(3,82) = ((sqrt(3._ark)*y6*y8/2._ark-y5*y8/2._ark)*y3+(-sqrt(3._ark)*y6*y7/2._ark-y5*y7/ &
      2._ark)*y4)*y1+((sqrt(3._ark)*y6*y7/2._ark+y5*y7/2._ark)*y3+(-sqrt(3._ark)*y6*y8/2._ark+ &
      y5*y8/2._ark)*y4)*y2
    dF(3,83) = ((y4*y9+y3*y9)*y2+y3*y4*y9)*y1+y2*y3*y4*y9
    dF(3,84) = ((y4*y5+y3*y5)*y2-y3*y4*y5)*y1-y2*y3*y4*y5
    dF(3,85) = (((y8-y7)*y3+(-y8+y7)*y4)*y2+(y8+y7)*y4*y3)*y1+(-y7-y8)*y4*y3*y2
    dF(3,86) = ((-sqrt(3._ark)*y6/2._ark+y5/2._ark)*y3+(sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y4)*y1**2+((sqrt(3._ark)*y6/2._ark-y5/2._ark)*y3**2+(-sqrt(3._ark)*y6/2._ark-y5/ &
      2._ark)*y4**2)*y1+((sqrt(3._ark)*y6/2._ark+y5/2._ark)*y3+(-sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y4)*y2**2+((-sqrt(3._ark)*y6/2._ark-y5/2._ark)*y3**2+(sqrt(3._ark)*y6/2._ark-y5/ &
      2._ark)*y4**2)*y2
    dF(3,87) = (y7*y9+y8*y9)*y1**2+(-y7*y9-y8*y9)*y2**2+(-y8*y9+y7*y9)*y3**2+(- &
      y7*y9+y8*y9)*y4**2
    dF(3,88) = y4**2*y7*y8+y3**2*y7*y8+y2**2*y7*y8+y1**2*y7*y8
    dF(3,89) = y2**2*y9**2+y1**2*y9**2-y3**2*y9**2-y4**2*y9**2
    dF(3,90) = 2._ark/3._ark*sqrt(3._ark)*y2**2*y5*y9+2._ark/3._ark*sqrt(3._ark)*y4**2*y5*y9+ &
      2._ark/3._ark*sqrt(3._ark)*y3**2*y5*y9+2._ark/3._ark*sqrt(3._ark)*y1**2*y5*y9
    dF(3,91) = (y6**2+y5**2)*y1**2+(y6**2+y5**2)*y2**2+(-y6**2-y5**2)*y3**2+(-y6**2- &
      y5**2)*y4**2
    dF(3,92) = (y8+y7)*y5*y1**2+(-y7-y8)*y5*y2**2+(y8-y7)*y5*y3**2+(-y8+y7)*y5*y4**2
    dF(3,93) = (-sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y1**2+(- &
      sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y2**2+(-sqrt(3._ark)*y5**2/6._ark+ &
      sqrt(3._ark)*y6**2/2._ark)*y3**2+(-sqrt(3._ark)*y5**2/6._ark+sqrt(3._ark)*y6**2/ &
      2._ark)*y4**2
    dF(3,94) = (y3*y8+y4*y7)*y1**2+(y3**2*y8+y4**2*y7)*y1+(-y3*y7-y4*y8)*y2**2+(- &
      y4**2*y8-y3**2*y7)*y2
    dF(3,95) = (y8+y7)*y2*y1**2+(-y7-y8)*y2**2*y1+(y8-y7)*y4*y3**2+(-y8+y7)*y4**2*y3
    dF(3,96) = (y4*y9+y3*y9)*y1**2+(y3**2*y9+y4**2*y9)*y1+(y4*y9+y3*y9)*y2**2+ &
      (y3**2*y9+y4**2*y9)*y2
    dF(3,97) = y3*y4**2*y9+y3**2*y4*y9+y1**2*y2*y9+y1*y2**2*y9
    dF(3,98) = (y4*y8+y3*y7)*y1**2+(-y4**2*y8-y3**2*y7)*y1+(-y4*y7-y3*y8)*y2**2+ &
      (y3**2*y8+y4**2*y7)*y2
    dF(3,99) = y1**2*y2**2-y3**2*y4**2
    dF(3,100) = ((sqrt(3._ark)*y5/2._ark+y6/2._ark)*y3+(-y6/2._ark+sqrt(3._ark)*y5/ &
      2._ark)*y4)*y1**2+((-y6/2._ark-sqrt(3._ark)*y5/2._ark)*y3**2+(-sqrt(3._ark)*y5/2._ark+y6/ &
      2._ark)*y4**2)*y1+((-y6/2._ark+sqrt(3._ark)*y5/2._ark)*y3+(sqrt(3._ark)*y5/2._ark+y6/ &
      2._ark)*y4)*y2**2+((-sqrt(3._ark)*y5/2._ark+y6/2._ark)*y3**2+(-y6/2._ark-sqrt(3._ark)*y5/ &
      2._ark)*y4**2)*y2
    dF(3,101) = (y4+y3)*y2*y1**2+((y4+y3)*y2**2-y3*y4**2-y3**2*y4)*y1+(-y3*y4**2- &
      y3**2*y4)*y2
    dF(3,102) = y1**2*y3*y4+(-y4**2-y3**2)*y2*y1+y2**2*y3*y4
    dF(3,103) = (y8+y7)*y1**3+(-y7-y8)*y2**3+(y8-y7)*y3**3+(-y8+y7)*y4**3
    dF(3,104) = y4**3*y9+y3**3*y9+y2**3*y9+y1**3*y9
    dF(3,105) = y8**2*y9**3+y7**2*y9**3
    dF(3,106) = y7**2*y8**2*y9
    dF(3,107) = y8**4*y9+y7**4*y9
    dF(3,108) = y9**5
    dF(3,109) = (sqrt(3._ark)*y7**3*y8+sqrt(3._ark)*y7*y8**3)*y5+(y7**3*y8-y7*y8**3)*y6
    dF(3,110) = (3._ark*y7**2*y9+3._ark*y8**2*y9)*y5**2+(-2._ark*sqrt(3._ark)*y7**2*y9+ &
      2._ark*sqrt(3._ark)*y8**2*y9)*y6*y5+(y8**2*y9+y7**2*y9)*y6**2
    dF(3,111) = -y6**4*y9/2._ark+5._ark/6._ark*y5**4*y9+3._ark*y5**2*y6**2*y9
    dF(3,112) = (-2._ark*y7**3*y8-2._ark*y7*y8**3)*y5
    dF(3,113) = -2._ark*y5*y7*y8*y9**2
    dF(3,114) = (-sqrt(3._ark)*y8**2*y9/2._ark-sqrt(3._ark)*y7**2*y9/2._ark)*y5**2+ &
      (2._ark*y7**2*y9-2._ark*y8**2*y9)*y6*y5+(-sqrt(3._ark)*y8**2*y9/2._ark- &
      sqrt(3._ark)*y7**2*y9/2._ark)*y6**2
    dF(3,115) = sqrt(3._ark)*y6**2*y9**3/6._ark-sqrt(3._ark)*y5**2*y9**3/2._ark
    dF(3,116) = sqrt(3._ark)*y5**2*y6**2*y9/2._ark+sqrt(3._ark)*y6**4*y9/4._ark-7._ark/ &
      36._ark*sqrt(3._ark)*y5**4*y9
    dF(3,117) = (-2._ark*y7**2*y9-2._ark*y8**2*y9)*y5**2+(2._ark*sqrt(3._ark)*y7**2*y9- &
      2._ark*sqrt(3._ark)*y8**2*y9)*y6*y5
    dF(3,118) = y6**2*y9**3+y5**2*y9**3
    dF(3,119) = -8._ark/3._ark*sqrt(3._ark)*y5*y6**2*y7*y8
    dF(3,120) = 3._ark/4._ark*y6**4*y9+y5**4*y9/12._ark-y5**2*y6**2*y9/2._ark
    dF(3,121) = y5**3*y7*y8-3._ark*y5*y6**2*y7*y8
    dF(3,122) = (y7**4+y8**4)*y1+(y7**4+y8**4)*y2+(-y8**4-y7**4)*y3+(-y8**4- &
      y7**4)*y4
    dF(3,123) = (-y7**3*y9-y8**3*y9)*y1+(y7**3*y9+y8**3*y9)*y2+(-y7**3*y9+ &
      y8**3*y9)*y3+(-y8**3*y9+y7**3*y9)*y4
    dF(3,124) = ((-sqrt(3._ark)*y8**3/2._ark-sqrt(3._ark)*y7**3/2._ark)*y5+(-y8**3/2._ark+ &
      y7**3/2._ark)*y6)*y1+((sqrt(3._ark)*y8**3/2._ark+sqrt(3._ark)*y7**3/2._ark)*y5+(-y7**3/ &
      2._ark+y8**3/2._ark)*y6)*y2+((sqrt(3._ark)*y7**3/2._ark-sqrt(3._ark)*y8**3/2._ark)*y5+(- &
      y8**3/2._ark-y7**3/2._ark)*y6)*y3+((-sqrt(3._ark)*y7**3/2._ark+sqrt(3._ark)*y8**3/ &
      2._ark)*y5+(y8**3/2._ark+y7**3/2._ark)*y6)*y4
    dF(3,125) = ((-2._ark*sqrt(3._ark)*y8-2._ark*sqrt(3._ark)*y7)*y5**3+(-2._ark*y8+ &
      2._ark*y7)*y6*y5**2)*y1+((2._ark*sqrt(3._ark)*y8+2._ark*sqrt(3._ark)*y7)*y5**3+(-2._ark*y7+ &
      2._ark*y8)*y6*y5**2)*y2+((-2._ark*sqrt(3._ark)*y8+2._ark*sqrt(3._ark)*y7)*y5**3+(-2._ark*y8- &
      2._ark*y7)*y6*y5**2)*y3+((-2._ark*sqrt(3._ark)*y7+2._ark*sqrt(3._ark)*y8)*y5**3+(2._ark*y8+ &
      2._ark*y7)*y6*y5**2)*y4
    dF(3,126) = y1**3*y7*y8+y4**3*y7*y8+y3**3*y7*y8+y2**3*y7*y8
    dF(3,127) = ((y8**3/2._ark+y7**3/2._ark)*y5+(sqrt(3._ark)*y7**3/2._ark-sqrt(3._ark)*y8**3/ &
      2._ark)*y6)*y1+((-y8**3/2._ark-y7**3/2._ark)*y5+(-sqrt(3._ark)*y7**3/2._ark+ &
      sqrt(3._ark)*y8**3/2._ark)*y6)*y2+((-y7**3/2._ark+y8**3/2._ark)*y5+(-sqrt(3._ark)*y8**3/ &
      2._ark-sqrt(3._ark)*y7**3/2._ark)*y6)*y3+((-y8**3/2._ark+y7**3/2._ark)*y5+ &
      (sqrt(3._ark)*y8**3/2._ark+sqrt(3._ark)*y7**3/2._ark)*y6)*y4
    dF(3,128) = ((-y8**2-y7**2)*y5**2+(-y8**2-y7**2)*y6**2)*y1+((-y8**2- &
      y7**2)*y5**2+(-y8**2-y7**2)*y6**2)*y2+((y8**2+y7**2)*y5**2+(y8**2+ &
      y7**2)*y6**2)*y3+((y8**2+y7**2)*y5**2+(y8**2+y7**2)*y6**2)*y4
    dF(3,129) = (y5**3*y9-3._ark*y9*y6**2*y5)*y1+(y5**3*y9-3._ark*y9*y6**2*y5)*y2+ &
      (y5**3*y9-3._ark*y9*y6**2*y5)*y3+(y5**3*y9-3._ark*y9*y6**2*y5)*y4
    dF(3,130) = ((-10._ark*y8-10._ark*y7)*y5**3+(-6._ark*sqrt(3._ark)*y8+ &
      6._ark*sqrt(3._ark)*y7)*y6*y5**2)*y1+((10._ark*y8+10._ark*y7)*y5**3+(6._ark*sqrt(3._ark)*y8- &
      6._ark*sqrt(3._ark)*y7)*y6*y5**2)*y2+((-10._ark*y8+10._ark*y7)*y5**3+(- &
      6._ark*sqrt(3._ark)*y7-6._ark*sqrt(3._ark)*y8)*y6*y5**2)*y3+((10._ark*y8-10._ark*y7)*y5**3+ &
      (6._ark*sqrt(3._ark)*y8+6._ark*sqrt(3._ark)*y7)*y6*y5**2)*y4
    dF(3,131) = ((-y8**2-y7**2)*y5+(-sqrt(3._ark)*y8**2+sqrt(3._ark)*y7**2)*y6)*y2*y1+ &
      ((y8**2+y7**2)*y5+(sqrt(3._ark)*y8**2-sqrt(3._ark)*y7**2)*y6)*y4*y3
    dF(3,132) = ((sqrt(3._ark)*y6/2._ark-y5/2._ark)*y3+(-sqrt(3._ark)*y6/2._ark-y5/ &
      2._ark)*y4)*y1**3+((-sqrt(3._ark)*y6/2._ark+y5/2._ark)*y3**3+(sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y4**3)*y1+((-sqrt(3._ark)*y6/2._ark-y5/2._ark)*y3+(sqrt(3._ark)*y6/2._ark-y5/ &
      2._ark)*y4)*y2**3+((sqrt(3._ark)*y6/2._ark+y5/2._ark)*y3**3+(-sqrt(3._ark)*y6/2._ark+y5/ &
      2._ark)*y4**3)*y2
    dF(3,133) = y3*y4**3*y5-y1**3*y2*y5-y1*y2**3*y5+y3**3*y4*y5
    dF(3,134) = ((y8+y7)*y5**2+(y8+y7)*y6**2)*y1**2+((-y7-y8)*y5**2+(-y7- &
      y8)*y6**2)*y2**2+((y8-y7)*y5**2+(y8-y7)*y6**2)*y3**2+((-y8+y7)*y5**2+(-y8+ &
      y7)*y6**2)*y4**2
    dF(3,135) = (y8-y7)*y6*y1**3+(-y8+y7)*y6*y2**3+(y8+y7)*y6*y3**3+(-y7- &
      y8)*y6*y4**3
    dF(3,136) = (y6**2+y5**2)*y1**3+(y6**2+y5**2)*y2**3+(-y6**2-y5**2)*y3**3+(- &
      y6**2-y5**2)*y4**3
    dF(3,137) = ((-sqrt(3._ark)*y7*y8**2/2._ark-sqrt(3._ark)*y7**2*y8/2._ark)*y5+(y7*y8**2/ &
      2._ark-y7**2*y8/2._ark)*y6)*y1+((sqrt(3._ark)*y7*y8**2/2._ark+sqrt(3._ark)*y7**2*y8/ &
      2._ark)*y5+(-y7*y8**2/2._ark+y7**2*y8/2._ark)*y6)*y2+((sqrt(3._ark)*y7*y8**2/2._ark- &
      sqrt(3._ark)*y7**2*y8/2._ark)*y5+(-y7**2*y8/2._ark-y7*y8**2/2._ark)*y6)*y3+ &
      ((sqrt(3._ark)*y7**2*y8/2._ark-sqrt(3._ark)*y7*y8**2/2._ark)*y5+(y7*y8**2/2._ark+y7**2*y8/ &
      2._ark)*y6)*y4
    dF(3,138) = -2._ark/3._ark*sqrt(3._ark)*y4*y5*y7*y8*y9-2._ark/ &
      3._ark*sqrt(3._ark)*y3*y5*y7*y8*y9+2._ark/3._ark*sqrt(3._ark)*y2*y5*y7*y8*y9+2._ark/ &
      3._ark*sqrt(3._ark)*y1*y5*y7*y8*y9
    dF(3,139) = ((-y8**2+y7**2)*y6*y5+(-sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/ &
      3._ark)*y6**2)*y1+((-y8**2+y7**2)*y6*y5+(-sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/ &
      3._ark)*y6**2)*y2+((-y7**2+y8**2)*y6*y5+(sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y7**2/ &
      3._ark)*y6**2)*y3+((-y7**2+y8**2)*y6*y5+(sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y7**2/ &
      3._ark)*y6**2)*y4
    dF(3,140) = (y4*y7*y9+y3*y8*y9)*y1**2+(-y4**2*y7*y9-y3**2*y8*y9)*y1+(-y4*y8*y9- &
      y3*y7*y9)*y2**2+(y4**2*y8*y9+y3**2*y7*y9)*y2
    dF(3,141) = (-y8**2-y7**2)*y2*y1**2+(-y8**2-y7**2)*y2**2*y1+(y8**2+ &
      y7**2)*y4*y3**2+(y8**2+y7**2)*y4**2*y3
    dF(3,142) = (y3*y9**2+y4*y9**2)*y1**2+(-y4**2*y9**2-y3**2*y9**2)*y1+(y3*y9**2+ &
      y4*y9**2)*y2**2+(-y4**2*y9**2-y3**2*y9**2)*y2
    dF(3,143) = (-y3*y8**2-y4*y7**2)*y1**2+(y3**2*y8**2+y4**2*y7**2)*y1+(-y3*y7**2- &
      y4*y8**2)*y2**2+(y4**2*y8**2+y3**2*y7**2)*y2
    dF(3,144) = (y7*y9+y8*y9)*y2*y1**2+(-y7*y9-y8*y9)*y2**2*y1+(-y8*y9+ &
      y7*y9)*y4*y3**2+(-y7*y9+y8*y9)*y4**2*y3
    dF(3,145) = (-y4*y8*y9-y3*y7*y9)*y1**2+(-y3**2*y7*y9-y4**2*y8*y9)*y1+(y4*y7*y9+ &
      y3*y8*y9)*y2**2+(y4**2*y7*y9+y3**2*y8*y9)*y2
    dF(3,146) = (y3*y7*y8+y4*y7*y8)*y1**2+(y3**2*y7*y8+y4**2*y7*y8)*y1+(y3*y7*y8+ &
      y4*y7*y8)*y2**2+(y3**2*y7*y8+y4**2*y7*y8)*y2
    dF(3,147) = ((sqrt(3._ark)*y5*y7+y6*y7)*y3+(sqrt(3._ark)*y5*y8-y6*y8)*y4)*y1**2+((- &
      y6*y7-sqrt(3._ark)*y5*y7)*y3**2+(-sqrt(3._ark)*y5*y8+y6*y8)*y4**2)*y1+((- &
      sqrt(3._ark)*y5*y8+y6*y8)*y3+(-y6*y7-sqrt(3._ark)*y5*y7)*y4)*y2**2+ &
      ((sqrt(3._ark)*y5*y8-y6*y8)*y3**2+(sqrt(3._ark)*y5*y7+y6*y7)*y4**2)*y2
    dF(3,148) = ((-y5*y9/2._ark-sqrt(3._ark)*y6*y9/2._ark)*y3+(sqrt(3._ark)*y6*y9/2._ark- &
      y5*y9/2._ark)*y4)*y1**2+((-y5*y9/2._ark-sqrt(3._ark)*y6*y9/2._ark)*y3**2+ &
      (sqrt(3._ark)*y6*y9/2._ark-y5*y9/2._ark)*y4**2)*y1+((sqrt(3._ark)*y6*y9/2._ark-y5*y9/ &
      2._ark)*y3+(-y5*y9/2._ark-sqrt(3._ark)*y6*y9/2._ark)*y4)*y2**2+((sqrt(3._ark)*y6*y9/2._ark- &
      y5*y9/2._ark)*y3**2+(-y5*y9/2._ark-sqrt(3._ark)*y6*y9/2._ark)*y4**2)*y2
    dF(3,149) = (-2._ark/3._ark*sqrt(3._ark)*y4*y6**2-2._ark/ &
      3._ark*sqrt(3._ark)*y3*y6**2)*y1**2+(2._ark/3._ark*sqrt(3._ark)*y3**2*y6**2+2._ark/ &
      3._ark*sqrt(3._ark)*y4**2*y6**2)*y1+(-2._ark/3._ark*sqrt(3._ark)*y4*y6**2-2._ark/ &
      3._ark*sqrt(3._ark)*y3*y6**2)*y2**2+(2._ark/3._ark*sqrt(3._ark)*y3**2*y6**2+2._ark/ &
      3._ark*sqrt(3._ark)*y4**2*y6**2)*y2
    dF(3,150) = ((-sqrt(3._ark)*y6**2/3._ark-y5*y6)*y3+(y5*y6-sqrt(3._ark)*y6**2/ &
      3._ark)*y4)*y1**2+((y5*y6+sqrt(3._ark)*y6**2/3._ark)*y3**2+(sqrt(3._ark)*y6**2/3._ark- &
      y5*y6)*y4**2)*y1+((y5*y6-sqrt(3._ark)*y6**2/3._ark)*y3+(-sqrt(3._ark)*y6**2/3._ark- &
      y5*y6)*y4)*y2**2+((sqrt(3._ark)*y6**2/3._ark-y5*y6)*y3**2+(y5*y6+sqrt(3._ark)*y6**2/ &
      3._ark)*y4**2)*y2
    dF(3,151) = (-y6**2-y5**2)*y2*y1**2+(-y6**2-y5**2)*y2**2*y1+(y6**2+ &
      y5**2)*y4*y3**2+(y6**2+y5**2)*y4**2*y3
    dF(3,152) = (y8+y7)*y2*y1**3+(-y7-y8)*y2**3*y1+(y8-y7)*y4*y3**3+(-y8+ &
      y7)*y4**3*y3
    dF(3,153) = (-y4*y7-y3*y8)*y1**3+(-y3**3*y8-y4**3*y7)*y1+(y4*y8+y3*y7)*y2**3+ &
      (y4**3*y8+y3**3*y7)*y2
    dF(3,154) = y1**3*y2*y9+y3*y4**3*y9+y3**3*y4*y9+y1*y2**3*y9
    dF(3,155) = (y4*y8+y3*y7)*y1**3+(-y4**3*y8-y3**3*y7)*y1+(-y4*y7-y3*y8)*y2**3+ &
      (y3**3*y8+y4**3*y7)*y2
    dF(3,156) = (-y4-y3)*y1**4+(y3**4+y4**4)*y1+(-y4-y3)*y2**4+(y3**4+y4**4)*y2
    dF(3,157) = (y3*y9**3+y4*y9**3)*y1+(y3*y9**3+y4*y9**3)*y2
    dF(3,158) = (y8**2*y9+y7**2*y9)*y2*y1+(y8**2*y9+y7**2*y9)*y4*y3
    dF(3,159) = (-y4*y7*y8**2-y3*y7**2*y8)*y1+(y4*y7**2*y8+y3*y7*y8**2)*y2
    dF(3,160) = y3*y4*y9**3+y1*y2*y9**3
    dF(3,161) = y3*y4*y7*y8*y9-y1*y2*y7*y8*y9
    dF(3,162) = (y3*y7**2*y9+y4*y8**2*y9)*y1+(y4*y7**2*y9+y3*y8**2*y9)*y2
    dF(3,163) = (-y4*y7*y9**2-y3*y8*y9**2)*y1+(y4*y8*y9**2+y3*y7*y9**2)*y2
    dF(3,164) = (y4*y7**2*y9+y3*y8**2*y9)*y1+(y3*y7**2*y9+y4*y8**2*y9)*y2
    dF(3,165) = (-y4*y7**3-y3*y8**3)*y1+(y3*y7**3+y4*y8**3)*y2
    dF(3,166) = (-2._ark*y8**2+2._ark*y7**2)*y6*y2*y1+(2._ark*y8**2-2._ark*y7**2)*y6*y4*y3
    dF(3,167) = ((y6*y7*y9-sqrt(3._ark)*y5*y7*y9)*y3+(-y6*y8*y9- &
      sqrt(3._ark)*y5*y8*y9)*y4)*y1+((sqrt(3._ark)*y5*y8*y9+y6*y8*y9)*y3+ &
      (sqrt(3._ark)*y5*y7*y9-y6*y7*y9)*y4)*y2
    dF(3,168) = (2._ark*y4*y5*y8*y9+2._ark*y3*y5*y7*y9)*y1+(-2._ark*y3*y5*y8*y9- &
      2._ark*y4*y5*y7*y9)*y2
    dF(3,169) = ((-sqrt(3._ark)*y5**2*y9/4._ark-y5*y6*y9/2._ark+sqrt(3._ark)*y6**2*y9/ &
      4._ark)*y3+(sqrt(3._ark)*y6**2*y9/4._ark-sqrt(3._ark)*y5**2*y9/4._ark+y5*y6*y9/ &
      2._ark)*y4)*y1+((sqrt(3._ark)*y6**2*y9/4._ark-sqrt(3._ark)*y5**2*y9/4._ark+y5*y6*y9/ &
      2._ark)*y3+(-sqrt(3._ark)*y5**2*y9/4._ark-y5*y6*y9/2._ark+sqrt(3._ark)*y6**2*y9/ &
      4._ark)*y4)*y2
    dF(3,170) = ((y5*y6*y8-sqrt(3._ark)*y5**2*y8/3._ark)*y3+(-y5*y6*y7- &
      sqrt(3._ark)*y5**2*y7/3._ark)*y4)*y1+((sqrt(3._ark)*y5**2*y7/3._ark+y5*y6*y7)*y3+ &
      (sqrt(3._ark)*y5**2*y8/3._ark-y5*y6*y8)*y4)*y2
    dF(3,171) = -2._ark/3._ark*sqrt(3._ark)*y1**2*y3*y4*y5+(2._ark/3._ark*sqrt(3._ark)*y3**2*y5+ &
      2._ark/3._ark*sqrt(3._ark)*y4**2*y5)*y2*y1-2._ark/3._ark*sqrt(3._ark)*y2**2*y3*y4*y5
    dF(3,172) = (y4*y9+y3*y9)*y2*y1**2+((y4*y9+y3*y9)*y2**2+y3*y4**2*y9+ &
      y3**2*y4*y9)*y1+(y3*y4**2*y9+y3**2*y4*y9)*y2
    dF(3,173) = (y4*y8+y3*y7)*y2*y1**2+((-y4*y7-y3*y8)*y2**2-y3**2*y4*y7- &
      y3*y4**2*y8)*y1+(y3**2*y4*y8+y3*y4**2*y7)*y2
    dF(3,174) = (2._ark*y3*y6-2._ark*y4*y6)*y2*y1**2+((2._ark*y4*y6-2._ark*y3*y6)*y2**2- &
      2._ark*y3**2*y4*y6+2._ark*y3*y4**2*y6)*y1+(-2._ark*y3*y4**2*y6+2._ark*y3**2*y4*y6)*y2
    dF(3,175) = y1**3*y3*y4+(-y3**3-y4**3)*y2*y1+y2**3*y3*y4
    dF(3,176) = (-y4-y3)*y2*y1**3+((-y4-y3)*y2**3+y3**3*y4+y3*y4**3)*y1+(y3**3*y4+ &
      y3*y4**3)*y2
    dF(3,177) = (y7*y8**2+y7**2*y8)*y1**2+(-y7**2*y8-y7*y8**2)*y2**2+(y7**2*y8- &
      y7*y8**2)*y3**2+(y7*y8**2-y7**2*y8)*y4**2
    dF(3,178) = (y8**2*y9+y7**2*y9)*y1**2+(y8**2*y9+y7**2*y9)*y2**2+(y8**2*y9+ &
      y7**2*y9)*y3**2+(y8**2*y9+y7**2*y9)*y4**2
    dF(3,179) = (-y7*y9**2-y8*y9**2)*y1**2+(y8*y9**2+y7*y9**2)*y2**2+(y7*y9**2- &
      y8*y9**2)*y3**2+(y8*y9**2-y7*y9**2)*y4**2
    dF(3,180) = ((-sqrt(3._ark)*y7*y9/2._ark-sqrt(3._ark)*y8*y9/2._ark)*y5+(-y7*y9/2._ark+ &
      y8*y9/2._ark)*y6)*y1**2+((sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5+(-y8*y9/ &
      2._ark+y7*y9/2._ark)*y6)*y2**2+((sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y5+(- &
      y8*y9/2._ark-y7*y9/2._ark)*y6)*y3**2+((-sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/ &
      2._ark)*y5+(y8*y9/2._ark+y7*y9/2._ark)*y6)*y4**2
    dF(3,181) = ((y8**2+y7**2)*y5+(sqrt(3._ark)*y8**2-sqrt(3._ark)*y7**2)*y6)*y1**2+ &
      ((y8**2+y7**2)*y5+(sqrt(3._ark)*y8**2-sqrt(3._ark)*y7**2)*y6)*y2**2+((-y8**2- &
      y7**2)*y5+(-sqrt(3._ark)*y8**2+sqrt(3._ark)*y7**2)*y6)*y3**2+((-y8**2-y7**2)*y5+(- &
      sqrt(3._ark)*y8**2+sqrt(3._ark)*y7**2)*y6)*y4**2
    dF(3,182) = ((-y8*y9/2._ark-y7*y9/2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark+ &
      sqrt(3._ark)*y7*y9/2._ark)*y6)*y1**2+((y8*y9/2._ark+y7*y9/2._ark)*y5+(sqrt(3._ark)*y8*y9/ &
      2._ark-sqrt(3._ark)*y7*y9/2._ark)*y6)*y2**2+((-y7*y9/2._ark+y8*y9/2._ark)*y5+ &
      (sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y3**2+((-y8*y9/2._ark+y7*y9/ &
      2._ark)*y5+(-sqrt(3._ark)*y7*y9/2._ark-sqrt(3._ark)*y8*y9/2._ark)*y6)*y4**2
    dF(3,183) = y3*y4**2*y7*y8+y3**2*y4*y7*y8+y1*y2**2*y7*y8+y1**2*y2*y7*y8
    dF(3,184) = ((y6*y8+sqrt(3._ark)*y5*y8)*y3+(sqrt(3._ark)*y5*y7-y6*y7)*y4)*y1**2+ &
      ((y6*y8+sqrt(3._ark)*y5*y8)*y3**2+(sqrt(3._ark)*y5*y7-y6*y7)*y4**2)*y1+((y6*y7- &
      sqrt(3._ark)*y5*y7)*y3+(-y6*y8-sqrt(3._ark)*y5*y8)*y4)*y2**2+((y6*y7- &
      sqrt(3._ark)*y5*y7)*y3**2+(-y6*y8-sqrt(3._ark)*y5*y8)*y4**2)*y2
    dF(3,185) = y4**3*y9**2+y3**3*y9**2-y1**3*y9**2-y2**3*y9**2
    dF(3,186) = (y8**2+y7**2)*y1**3+(y8**2+y7**2)*y2**3+(-y8**2-y7**2)*y3**3+(- &
      y8**2-y7**2)*y4**3
    dF(3,187) = (-y7*y9-y8*y9)*y1**3+(y7*y9+y8*y9)*y2**3+(-y7*y9+y8*y9)*y3**3+(- &
      y8*y9+y7*y9)*y4**3
    dF(3,188) = (-y7-y8)*y5*y1**3+(y8+y7)*y5*y2**3+(-y8+y7)*y5*y3**3+(y8- &
      y7)*y5*y4**3
    dF(3,189) = -y3*y7**2*y8**2-y4*y7**2*y8**2+y2*y7**2*y8**2+y1*y7**2*y8**2
    dF(3,190) = (y7**2*y9**2+y8**2*y9**2)*y1+(y7**2*y9**2+y8**2*y9**2)*y2+(- &
      y7**2*y9**2-y8**2*y9**2)*y3+(-y7**2*y9**2-y8**2*y9**2)*y4
    dF(3,191) = (y8*y9**3+y7*y9**3)*y1+(-y7*y9**3-y8*y9**3)*y2+(y7*y9**3- &
      y8*y9**3)*y3+(y8*y9**3-y7*y9**3)*y4
    dF(3,192) = y2*y9**4+y1*y9**4-y4*y9**4-y3*y9**4
    dF(3,193) = (y7*y8**3+y7**3*y8)*y1+(y7*y8**3+y7**3*y8)*y2+(y7*y8**3+ &
      y7**3*y8)*y3+(y7*y8**3+y7**3*y8)*y4
    dF(3,194) = y1*y7*y8*y9**2+y2*y7*y8*y9**2+y3*y7*y8*y9**2+y4*y7*y8*y9**2
    dF(3,195) = (y7*y8**2*y9+y7**2*y8*y9)*y1+(-y7**2*y8*y9-y7*y8**2*y9)*y2+ &
      (y7*y8**2*y9-y7**2*y8*y9)*y3+(-y7*y8**2*y9+y7**2*y8*y9)*y4
    dF(3,196) = ((-sqrt(3._ark)*y8**2*y9-sqrt(3._ark)*y7**2*y9)*y5+(y8**2*y9- &
      y7**2*y9)*y6)*y1+((-sqrt(3._ark)*y8**2*y9-sqrt(3._ark)*y7**2*y9)*y5+(y8**2*y9- &
      y7**2*y9)*y6)*y2+((-sqrt(3._ark)*y8**2*y9-sqrt(3._ark)*y7**2*y9)*y5+(y8**2*y9- &
      y7**2*y9)*y6)*y3+((-sqrt(3._ark)*y8**2*y9-sqrt(3._ark)*y7**2*y9)*y5+(y8**2*y9- &
      y7**2*y9)*y6)*y4
    dF(3,197) = ((sqrt(3._ark)*y7*y9**2+sqrt(3._ark)*y8*y9**2)*y5+(y7*y9**2- &
      y8*y9**2)*y6)*y1+((-sqrt(3._ark)*y8*y9**2-sqrt(3._ark)*y7*y9**2)*y5+(y8*y9**2- &
      y7*y9**2)*y6)*y2+((sqrt(3._ark)*y8*y9**2-sqrt(3._ark)*y7*y9**2)*y5+(-y7*y9**2- &
      y8*y9**2)*y6)*y3+((-sqrt(3._ark)*y8*y9**2+sqrt(3._ark)*y7*y9**2)*y5+(y8*y9**2+ &
      y7*y9**2)*y6)*y4
    dF(3,198) = (y6**2*y9**2+y5**2*y9**2)*y1+(y6**2*y9**2+y5**2*y9**2)*y2+(- &
      y6**2*y9**2-y5**2*y9**2)*y3+(-y6**2*y9**2-y5**2*y9**2)*y4
    dF(3,199) = (y5**2*y7*y8+y6**2*y7*y8)*y1+(y5**2*y7*y8+y6**2*y7*y8)*y2+ &
      (y5**2*y7*y8+y6**2*y7*y8)*y3+(y5**2*y7*y8+y6**2*y7*y8)*y4
    dF(3,200) = ((sqrt(3._ark)*y7*y9-sqrt(3._ark)*y8*y9)*y6*y5+(y7*y9+y8*y9)*y6**2)*y1+ &
      ((sqrt(3._ark)*y8*y9-sqrt(3._ark)*y7*y9)*y6*y5+(-y7*y9-y8*y9)*y6**2)*y2+ &
      ((sqrt(3._ark)*y7*y9+sqrt(3._ark)*y8*y9)*y6*y5+(-y8*y9+y7*y9)*y6**2)*y3+((- &
      sqrt(3._ark)*y8*y9-sqrt(3._ark)*y7*y9)*y6*y5+(-y7*y9+y8*y9)*y6**2)*y4
    dF(3,201) = (-sqrt(3._ark)*y5*y6**2*y9/3._ark+sqrt(3._ark)*y5**3*y9)*y1+(- &
      sqrt(3._ark)*y5*y6**2*y9/3._ark+sqrt(3._ark)*y5**3*y9)*y2+(-sqrt(3._ark)*y5*y6**2*y9/ &
      3._ark+sqrt(3._ark)*y5**3*y9)*y3+(-sqrt(3._ark)*y5*y6**2*y9/3._ark+ &
      sqrt(3._ark)*y5**3*y9)*y4
    dF(3,202) = ((6._ark*sqrt(3._ark)*y8+6._ark*sqrt(3._ark)*y7)*y5**3+(-9._ark*y7+ &
      9._ark*y8)*y6*y5**2+(-y8+y7)*y6**3)*y1+((-6._ark*sqrt(3._ark)*y7- &
      6._ark*sqrt(3._ark)*y8)*y5**3+(-9._ark*y8+9._ark*y7)*y6*y5**2+(y8-y7)*y6**3)*y2+ &
      ((6._ark*sqrt(3._ark)*y8-6._ark*sqrt(3._ark)*y7)*y5**3+(9._ark*y7+9._ark*y8)*y6*y5**2+(-y7- &
      y8)*y6**3)*y3+((-6._ark*sqrt(3._ark)*y8+6._ark*sqrt(3._ark)*y7)*y5**3+(-9._ark*y7- &
      9._ark*y8)*y6*y5**2+(y8+y7)*y6**3)*y4
    dF(3,203) = (y6**4/4._ark+9._ark/4._ark*y5**4-3._ark/2._ark*y5**2*y6**2)*y1+(y6**4/4._ark+ &
      9._ark/4._ark*y5**4-3._ark/2._ark*y5**2*y6**2)*y2+(3._ark/2._ark*y5**2*y6**2-9._ark/ &
      4._ark*y5**4-y6**4/4._ark)*y3+(3._ark/2._ark*y5**2*y6**2-9._ark/4._ark*y5**4-y6**4/4._ark)*y4
    dF(3,204) = y1*y5*y9**3+y4*y5*y9**3+y3*y5*y9**3+y2*y5*y9**3
    dF(3,205) = ((y8*y9**2+y7*y9**2)*y5+(-sqrt(3._ark)*y8*y9**2+ &
      sqrt(3._ark)*y7*y9**2)*y6)*y1+((-y7*y9**2-y8*y9**2)*y5+(sqrt(3._ark)*y8*y9**2- &
      sqrt(3._ark)*y7*y9**2)*y6)*y2+((y8*y9**2-y7*y9**2)*y5+(-sqrt(3._ark)*y8*y9**2- &
      sqrt(3._ark)*y7*y9**2)*y6)*y3+((y7*y9**2-y8*y9**2)*y5+(sqrt(3._ark)*y7*y9**2+ &
      sqrt(3._ark)*y8*y9**2)*y6)*y4
    dF(3,206) = (-2._ark*y7**2*y9-2._ark*y8**2*y9)*y5*y1+(-2._ark*y7**2*y9- &
      2._ark*y8**2*y9)*y5*y2+(-2._ark*y7**2*y9-2._ark*y8**2*y9)*y5*y3+(-2._ark*y7**2*y9- &
      2._ark*y8**2*y9)*y5*y4
    dF(3,207) = ((-y7**2*y8/2._ark-y7*y8**2/2._ark)*y5+(sqrt(3._ark)*y7**2*y8/2._ark- &
      sqrt(3._ark)*y7*y8**2/2._ark)*y6)*y1+((y7*y8**2/2._ark+y7**2*y8/2._ark)*y5+ &
      (sqrt(3._ark)*y7*y8**2/2._ark-sqrt(3._ark)*y7**2*y8/2._ark)*y6)*y2+((y7*y8**2/2._ark- &
      y7**2*y8/2._ark)*y5+(sqrt(3._ark)*y7*y8**2/2._ark+sqrt(3._ark)*y7**2*y8/2._ark)*y6)*y3+((- &
      y7*y8**2/2._ark+y7**2*y8/2._ark)*y5+(-sqrt(3._ark)*y7*y8**2/2._ark-sqrt(3._ark)*y7**2*y8/ &
      2._ark)*y6)*y4
    dF(3,208) = (-2._ark/3._ark*sqrt(3._ark)*y8**2-2._ark/3._ark*sqrt(3._ark)*y7**2)*y6**2*y1+(- &
      2._ark/3._ark*sqrt(3._ark)*y8**2-2._ark/3._ark*sqrt(3._ark)*y7**2)*y6**2*y2+(2._ark/ &
      3._ark*sqrt(3._ark)*y8**2+2._ark/3._ark*sqrt(3._ark)*y7**2)*y6**2*y3+(2._ark/ &
      3._ark*sqrt(3._ark)*y8**2+2._ark/3._ark*sqrt(3._ark)*y7**2)*y6**2*y4
    dF(3,209) = (-sqrt(3._ark)*y6**2*y7*y8/2._ark+sqrt(3._ark)*y5**2*y7*y8/6._ark)*y1+(- &
      sqrt(3._ark)*y6**2*y7*y8/2._ark+sqrt(3._ark)*y5**2*y7*y8/6._ark)*y2+(- &
      sqrt(3._ark)*y6**2*y7*y8/2._ark+sqrt(3._ark)*y5**2*y7*y8/6._ark)*y3+(- &
      sqrt(3._ark)*y6**2*y7*y8/2._ark+sqrt(3._ark)*y5**2*y7*y8/6._ark)*y4
    dF(3,210) = (-2._ark*y8*y9+2._ark*y7*y9)*y6*y5*y1+(2._ark*y8*y9-2._ark*y7*y9)*y6*y5*y2+ &
      (2._ark*y8*y9+2._ark*y7*y9)*y6*y5*y3+(-2._ark*y7*y9-2._ark*y8*y9)*y6*y5*y4
    dF(3,211) = (sqrt(3._ark)*y6**2*y9**2/2._ark-sqrt(3._ark)*y5**2*y9**2/6._ark)*y1+ &
      (sqrt(3._ark)*y6**2*y9**2/2._ark-sqrt(3._ark)*y5**2*y9**2/6._ark)*y2+(- &
      sqrt(3._ark)*y6**2*y9**2/2._ark+sqrt(3._ark)*y5**2*y9**2/6._ark)*y3+(- &
      sqrt(3._ark)*y6**2*y9**2/2._ark+sqrt(3._ark)*y5**2*y9**2/6._ark)*y4
    dF(3,212) = ((3._ark*y8+3._ark*y7)*y5**3+(-2._ark*sqrt(3._ark)*y7+ &
      2._ark*sqrt(3._ark)*y8)*y6*y5**2+(y8+y7)*y6**2*y5)*y1+((-3._ark*y7-3._ark*y8)*y5**3+(- &
      2._ark*sqrt(3._ark)*y8+2._ark*sqrt(3._ark)*y7)*y6*y5**2+(-y7-y8)*y6**2*y5)*y2+((3._ark*y8- &
      3._ark*y7)*y5**3+(2._ark*sqrt(3._ark)*y8+2._ark*sqrt(3._ark)*y7)*y6*y5**2+(y8- &
      y7)*y6**2*y5)*y3+((-3._ark*y8+3._ark*y7)*y5**3+(-2._ark*sqrt(3._ark)*y8- &
      2._ark*sqrt(3._ark)*y7)*y6*y5**2+(-y8+y7)*y6**2*y5)*y4
    dF(3,213) = ((y7*y9+y8*y9)*y5**2+(sqrt(3._ark)*y8*y9-sqrt(3._ark)*y7*y9)*y6*y5)*y1+ &
      ((-y7*y9-y8*y9)*y5**2+(sqrt(3._ark)*y7*y9-sqrt(3._ark)*y8*y9)*y6*y5)*y2+((-y8*y9+ &
      y7*y9)*y5**2+(-sqrt(3._ark)*y8*y9-sqrt(3._ark)*y7*y9)*y6*y5)*y3+((-y7*y9+ &
      y8*y9)*y5**2+(sqrt(3._ark)*y7*y9+sqrt(3._ark)*y8*y9)*y6*y5)*y4
    dF(3,214) = (5._ark/6._ark*sqrt(3._ark)*y5**2*y6**2-sqrt(3._ark)*y5**4/4._ark- &
      sqrt(3._ark)*y6**4/4._ark)*y1+(5._ark/6._ark*sqrt(3._ark)*y5**2*y6**2-sqrt(3._ark)*y5**4/ &
      4._ark-sqrt(3._ark)*y6**4/4._ark)*y2+(sqrt(3._ark)*y5**4/4._ark+sqrt(3._ark)*y6**4/4._ark- &
      5._ark/6._ark*sqrt(3._ark)*y5**2*y6**2)*y3+(sqrt(3._ark)*y5**4/4._ark+sqrt(3._ark)*y6**4/ &
      4._ark-5._ark/6._ark*sqrt(3._ark)*y5**2*y6**2)*y4
    dF(3,215) = (-5._ark/4._ark*y5**4+3._ark/4._ark*y6**4+7._ark/2._ark*y5**2*y6**2)*y1+(-5._ark/ &
      4._ark*y5**4+3._ark/4._ark*y6**4+7._ark/2._ark*y5**2*y6**2)*y2+(-3._ark/4._ark*y6**4+5._ark/ &
      4._ark*y5**4-7._ark/2._ark*y5**2*y6**2)*y3+(-3._ark/4._ark*y6**4+5._ark/4._ark*y5**4-7._ark/ &
      2._ark*y5**2*y6**2)*y4
    dF(3,216) = (y6**2*y9+y5**2*y9)*y2*y1+(y6**2*y9+y5**2*y9)*y4*y3
    dF(3,217) = (-sqrt(3._ark)*y6**2*y9/2._ark+sqrt(3._ark)*y5**2*y9/6._ark)*y2*y1+(- &
      sqrt(3._ark)*y6**2*y9/2._ark+sqrt(3._ark)*y5**2*y9/6._ark)*y4*y3
    dF(3,218) = ((sqrt(3._ark)*y5*y6*y9/2._ark+y5**2*y9/4._ark+3._ark/4._ark*y6**2*y9)*y3+(- &
      sqrt(3._ark)*y5*y6*y9/2._ark+y5**2*y9/4._ark+3._ark/4._ark*y6**2*y9)*y4)*y1+((- &
      sqrt(3._ark)*y5*y6*y9/2._ark+y5**2*y9/4._ark+3._ark/4._ark*y6**2*y9)*y3+ &
      (sqrt(3._ark)*y5*y6*y9/2._ark+y5**2*y9/4._ark+3._ark/4._ark*y6**2*y9)*y4)*y2
    dF(3,219) = (y5**3-3._ark*y5*y6**2)*y2*y1+(-y5**3+3._ark*y5*y6**2)*y4*y3
    dF(3,220) = y1**2*y2*y9**2+y1*y2**2*y9**2-y3*y4**2*y9**2-y3**2*y4*y9**2
    dF(3,221) = (-y3*y7**2-y4*y8**2)*y1**2+(y4**2*y8**2+y3**2*y7**2)*y1+(-y3*y8**2- &
      y4*y7**2)*y2**2+(y3**2*y8**2+y4**2*y7**2)*y2
    dF(3,222) = ((-sqrt(3._ark)*y5*y6*y9/2._ark+3._ark/4._ark*y5**2*y9+y6**2*y9/4._ark)*y3+ &
      (sqrt(3._ark)*y5*y6*y9/2._ark+3._ark/4._ark*y5**2*y9+y6**2*y9/4._ark)*y4)*y1+ &
      ((sqrt(3._ark)*y5*y6*y9/2._ark+3._ark/4._ark*y5**2*y9+y6**2*y9/4._ark)*y3+(- &
      sqrt(3._ark)*y5*y6*y9/2._ark+3._ark/4._ark*y5**2*y9+y6**2*y9/4._ark)*y4)*y2
    dF(3,223) = ((y5**2*y8+y6**2*y8)*y3+(y6**2*y7+y5**2*y7)*y4)*y1+((-y6**2*y7- &
      y5**2*y7)*y3+(-y6**2*y8-y5**2*y8)*y4)*y2
    dF(3,224) = (-sqrt(3._ark)*y5*y6**2/3._ark+sqrt(3._ark)*y5**3)*y2*y1+ &
      (sqrt(3._ark)*y5*y6**2/3._ark-sqrt(3._ark)*y5**3)*y4*y3
    dF(3,225) = ((y5*y7*y8-sqrt(3._ark)*y6*y7*y8)*y3+(sqrt(3._ark)*y6*y7*y8+ &
      y5*y7*y8)*y4)*y1+((sqrt(3._ark)*y6*y7*y8+y5*y7*y8)*y3+(y5*y7*y8- &
      sqrt(3._ark)*y6*y7*y8)*y4)*y2
    dF(3,226) = 2._ark*y3*y4*y5*y9**2-2._ark*y1*y2*y5*y9**2
    dF(3,227) = (2._ark/3._ark*sqrt(3._ark)*y3*y5**2*y8+2._ark/ &
      3._ark*sqrt(3._ark)*y4*y5**2*y7)*y1+(-2._ark/3._ark*sqrt(3._ark)*y3*y5**2*y7-2._ark/ &
      3._ark*sqrt(3._ark)*y4*y5**2*y8)*y2
    dF(3,228) = (((-y8*y9+y7*y9)*y3+(-y7*y9+y8*y9)*y4)*y2+(y7*y9+y8*y9)*y4*y3)*y1+(- &
      y7*y9-y8*y9)*y4*y3*y2
    dF(3,229) = ((y3*y9**2+y4*y9**2)*y2-y3*y4*y9**2)*y1-y2*y3*y4*y9**2
    dF(3,230) = (((-y6**2-y5**2)*y3+(-y6**2-y5**2)*y4)*y2+(y6**2+y5**2)*y4*y3)*y1+ &
      (y6**2+y5**2)*y4*y3*y2
    dF(3,231) = (y8+y7)*y4*y3*y1**2+((y8-y7)*y3**2+(-y8+y7)*y4**2)*y2*y1+(-y7- &
      y8)*y4*y3*y2**2
    dF(3,232) = ((-y6**2-y5**2)*y3+(-y6**2-y5**2)*y4)*y1**2+((y6**2+y5**2)*y3**2+ &
      (y6**2+y5**2)*y4**2)*y1+((-y6**2-y5**2)*y3+(-y6**2-y5**2)*y4)*y2**2+((y6**2+ &
      y5**2)*y3**2+(y6**2+y5**2)*y4**2)*y2
    dF(3,233) = (y4*y9+y3*y9)*y1**3+(y4**3*y9+y3**3*y9)*y1+(y4*y9+y3*y9)*y2**3+ &
      (y4**3*y9+y3**3*y9)*y2
    dF(3,234) = (2._ark*y3*y6*y7*y8-2._ark*y4*y6*y7*y8)*y1+(-2._ark*y3*y6*y7*y8+ &
      2._ark*y4*y6*y7*y8)*y2
    dF(3,235) = y1*y2*y5*y7*y8+y3*y4*y5*y7*y8
    dF(3,236) = (((y8**2+y7**2)*y3+(y8**2+y7**2)*y4)*y2+(-y8**2-y7**2)*y4*y3)*y1+(- &
      y8**2-y7**2)*y4*y3*y2
    dF(3,237) = ((y3*y7*y8+y4*y7*y8)*y2+y3*y4*y7*y8)*y1+y2*y3*y4*y7*y8
    dF(3,238) = (((-y7-y8)*y6*y3+(y8+y7)*y6*y4)*y2+(-y8+y7)*y6*y4*y3)*y1+(y8- &
      y7)*y6*y4*y3*y2
    dF(3,239) = ((-2._ark*y3*y5*y9-2._ark*y4*y5*y9)*y2-2._ark*y3*y4*y5*y9)*y1- &
      2._ark*y2*y3*y4*y5*y9
    dF(3,240) = (((-y8+y7)*y5*y3+(y8-y7)*y5*y4)*y2+(-y7-y8)*y5*y4*y3)*y1+(y8+ &
      y7)*y5*y4*y3*y2
    dF(3,241) = ((2._ark/3._ark*sqrt(3._ark)*y3*y6**2+2._ark/3._ark*sqrt(3._ark)*y4*y6**2)*y2- &
      2._ark/3._ark*sqrt(3._ark)*y3*y4*y6**2)*y1-2._ark/3._ark*sqrt(3._ark)*y2*y3*y4*y6**2
    dF(3,242) = y1*y2*y3*y4*y9
    dF(3,243) = (-y4*y7-y3*y8)*y2*y1**2+((y4*y8+y3*y7)*y2**2-y3**2*y4*y8- &
      y3*y4**2*y7)*y1+(y3*y4**2*y8+y3**2*y4*y7)*y2
    dF(3,244) = y1**2*y9**3+y4**2*y9**3+y3**2*y9**3+y2**2*y9**3
    dF(3,245) = (y7**3+y8**3)*y1**2+(-y7**3-y8**3)*y2**2+(y8**3-y7**3)*y3**2+(y7**3- &
      y8**3)*y4**2
    dF(3,246) = -y4**2*y7*y8*y9+y1**2*y7*y8*y9-y3**2*y7*y8*y9+y2**2*y7*y8*y9
    dF(3,247) = (-2._ark*y8**2+2._ark*y7**2)*y6*y1**2+(-2._ark*y8**2+2._ark*y7**2)*y6*y2**2+ &
      (2._ark*y8**2-2._ark*y7**2)*y6*y3**2+(2._ark*y8**2-2._ark*y7**2)*y6*y4**2
    dF(3,248) = -2._ark/3._ark*sqrt(3._ark)*y4**2*y5*y9**2-2._ark/ &
      3._ark*sqrt(3._ark)*y3**2*y5*y9**2+2._ark/3._ark*sqrt(3._ark)*y2**2*y5*y9**2+2._ark/ &
      3._ark*sqrt(3._ark)*y1**2*y5*y9**2
    dF(3,249) = -8._ark/3._ark*sqrt(3._ark)*y1**2*y5*y6**2+8._ark/ &
      3._ark*sqrt(3._ark)*y4**2*y5*y6**2+8._ark/3._ark*sqrt(3._ark)*y3**2*y5*y6**2-8._ark/ &
      3._ark*sqrt(3._ark)*y2**2*y5*y6**2
    dF(3,250) = y4**2*y5*y7*y8+y3**2*y5*y7*y8+y2**2*y5*y7*y8+y1**2*y5*y7*y8
    dF(3,251) = (2._ark/3._ark*sqrt(3._ark)*y8+2._ark/3._ark*sqrt(3._ark)*y7)*y5**2*y1**2+(- &
      2._ark/3._ark*sqrt(3._ark)*y8-2._ark/3._ark*sqrt(3._ark)*y7)*y5**2*y2**2+(2._ark/ &
      3._ark*sqrt(3._ark)*y8-2._ark/3._ark*sqrt(3._ark)*y7)*y5**2*y3**2+(-2._ark/ &
      3._ark*sqrt(3._ark)*y8+2._ark/3._ark*sqrt(3._ark)*y7)*y5**2*y4**2
    dF(3,252) = (sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y1**2+ &
      (sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y2**2+(sqrt(3._ark)*y6**2*y9/ &
      6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y3**2+(sqrt(3._ark)*y6**2*y9/6._ark- &
      sqrt(3._ark)*y5**2*y9/2._ark)*y4**2
    dF(3,253) = ((sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/3._ark)*y5**2+(-y8+ &
      y7)*y6*y5)*y1**2+((-sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y7/3._ark)*y5**2+(y8- &
      y7)*y6*y5)*y2**2+((sqrt(3._ark)*y8/3._ark-sqrt(3._ark)*y7/3._ark)*y5**2+(-y7- &
      y8)*y6*y5)*y3**2+((-sqrt(3._ark)*y8/3._ark+sqrt(3._ark)*y7/3._ark)*y5**2+(y8+ &
      y7)*y6*y5)*y4**2
    dF(3,254) = (y6**2*y9+y5**2*y9)*y1**2+(y6**2*y9+y5**2*y9)*y2**2+(y6**2*y9+ &
      y5**2*y9)*y3**2+(y6**2*y9+y5**2*y9)*y4**2
    dF(3,255) = (y5**3-3._ark*y5*y6**2)*y1**2+(y5**3-3._ark*y5*y6**2)*y2**2+(-y5**3+ &
      3._ark*y5*y6**2)*y3**2+(-y5**3+3._ark*y5*y6**2)*y4**2
    dF(3,256) = y1**2*y2**2*y9+y3**2*y4**2*y9
    dF(3,257) = (y3**2*y9+y4**2*y9)*y1**2+(y3**2*y9+y4**2*y9)*y2**2
    dF(3,258) = (y3**2*y8+y4**2*y7)*y1**2+(-y4**2*y8-y3**2*y7)*y2**2
    dF(3,259) = -2._ark/3._ark*sqrt(3._ark)*y1**2*y2**2*y5+2._ark/ &
      3._ark*sqrt(3._ark)*y3**2*y4**2*y5
    dF(3,260) = (-2._ark*y7+2._ark*y8)*y6*y2*y1**2+(-2._ark*y8+2._ark*y7)*y6*y2**2*y1+ &
      (2._ark*y8+2._ark*y7)*y6*y4*y3**2+(-2._ark*y8-2._ark*y7)*y6*y4**2*y3
    dF(3,261) = ((y6*y9/2._ark-sqrt(3._ark)*y5*y9/2._ark)*y3+(-sqrt(3._ark)*y5*y9/2._ark- &
      y6*y9/2._ark)*y4)*y1**2+((y6*y9/2._ark-sqrt(3._ark)*y5*y9/2._ark)*y3**2+(- &
      sqrt(3._ark)*y5*y9/2._ark-y6*y9/2._ark)*y4**2)*y1+((-sqrt(3._ark)*y5*y9/2._ark-y6*y9/ &
      2._ark)*y3+(y6*y9/2._ark-sqrt(3._ark)*y5*y9/2._ark)*y4)*y2**2+((-sqrt(3._ark)*y5*y9/2._ark- &
      y6*y9/2._ark)*y3**2+(y6*y9/2._ark-sqrt(3._ark)*y5*y9/2._ark)*y4**2)*y2
    dF(3,262) = (sqrt(3._ark)*y6**2/6._ark-sqrt(3._ark)*y5**2/2._ark)*y2*y1**2+ &
      (sqrt(3._ark)*y6**2/6._ark-sqrt(3._ark)*y5**2/2._ark)*y2**2*y1+(sqrt(3._ark)*y5**2/2._ark- &
      sqrt(3._ark)*y6**2/6._ark)*y4*y3**2+(sqrt(3._ark)*y5**2/2._ark-sqrt(3._ark)*y6**2/ &
      6._ark)*y4**2*y3
    dF(3,263) = y1**2*y3*y4*y9+(y3**2*y9+y4**2*y9)*y2*y1+y2**2*y3*y4*y9
    dF(3,264) = (y4+y3)*y2**2*y1**2-y1*y3**2*y4**2-y2*y3**2*y4**2
    dF(3,265) = ((-y4**2-y3**2)*y2+y3**2*y4+y3*y4**2)*y1**2+(-y4**2-y3**2)*y2**2*y1+ &
      (y3**2*y4+y3*y4**2)*y2**2
    dF(3,266) = ((sqrt(3._ark)*y6*y7+y5*y7)*y3+(-sqrt(3._ark)*y6*y8+y5*y8)*y4)*y1**2+((- &
      y5*y7-sqrt(3._ark)*y6*y7)*y3**2+(sqrt(3._ark)*y6*y8-y5*y8)*y4**2)*y1+ &
      ((sqrt(3._ark)*y6*y8-y5*y8)*y3+(-y5*y7-sqrt(3._ark)*y6*y7)*y4)*y2**2+((- &
      sqrt(3._ark)*y6*y8+y5*y8)*y3**2+(sqrt(3._ark)*y6*y7+y5*y7)*y4**2)*y2
    dF(3,267) = y1**2*y2*y5*y9+y3*y4**2*y5*y9+y3**2*y4*y5*y9+y1*y2**2*y5*y9
    dF(3,268) = ((y5*y8+sqrt(3._ark)*y6*y8)*y3+(-sqrt(3._ark)*y6*y7+y5*y7)*y4)*y1**2+ &
      ((y5*y8+sqrt(3._ark)*y6*y8)*y3**2+(-sqrt(3._ark)*y6*y7+y5*y7)*y4**2)*y1+((-y5*y7+ &
      sqrt(3._ark)*y6*y7)*y3+(-y5*y8-sqrt(3._ark)*y6*y8)*y4)*y2**2+((-y5*y7+ &
      sqrt(3._ark)*y6*y7)*y3**2+(-y5*y8-sqrt(3._ark)*y6*y8)*y4**2)*y2
    dF(3,269) = ((y8+y7)*y5+(sqrt(3._ark)*y7-sqrt(3._ark)*y8)*y6)*y2*y1**2+((-y7-y8)*y5+ &
      (sqrt(3._ark)*y8-sqrt(3._ark)*y7)*y6)*y2**2*y1+((y8-y7)*y5+(-sqrt(3._ark)*y7- &
      sqrt(3._ark)*y8)*y6)*y4*y3**2+((-y8+y7)*y5+(sqrt(3._ark)*y8+ &
      sqrt(3._ark)*y7)*y6)*y4**2*y3
    dF(3,270) = ((y5-sqrt(3._ark)*y6)*y3+(sqrt(3._ark)*y6+y5)*y4)*y2*y1**2+ &
      (((sqrt(3._ark)*y6+y5)*y3+(y5-sqrt(3._ark)*y6)*y4)*y2**2+(sqrt(3._ark)*y6- &
      y5)*y4*y3**2+(-y5-sqrt(3._ark)*y6)*y4**2*y3)*y1+((-y5-sqrt(3._ark)*y6)*y4*y3**2+ &
      (sqrt(3._ark)*y6-y5)*y4**2*y3)*y2
    dF(3,271) = y1**2*y2*y3*y4+(y2**2*y3*y4+(-y3*y4**2-y3**2*y4)*y2)*y1
    dF(3,272) = 2._ark/3._ark*sqrt(3._ark)*y1**3*y5*y9+2._ark/3._ark*sqrt(3._ark)*y4**3*y5*y9+ &
      2._ark/3._ark*sqrt(3._ark)*y3**3*y5*y9+2._ark/3._ark*sqrt(3._ark)*y2**3*y5*y9
    dF(3,273) = (-sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y1**3+(- &
      sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y2**3+(-sqrt(3._ark)*y5**2/6._ark+ &
      sqrt(3._ark)*y6**2/2._ark)*y3**3+(-sqrt(3._ark)*y5**2/6._ark+sqrt(3._ark)*y6**2/ &
      2._ark)*y4**3
    dF(3,274) = y1**2*y2**3+y1**3*y2**2-y3**2*y4**3-y3**3*y4**2
    dF(3,275) = (y4**2+y3**2)*y1**3+(-y3**3-y4**3)*y1**2+(y4**2+y3**2)*y2**3+(- &
      y3**3-y4**3)*y2**2
    dF(3,276) = ((sqrt(3._ark)*y5/2._ark+y6/2._ark)*y3+(-y6/2._ark+sqrt(3._ark)*y5/ &
      2._ark)*y4)*y1**3+((-y6/2._ark-sqrt(3._ark)*y5/2._ark)*y3**3+(-sqrt(3._ark)*y5/2._ark+y6/ &
      2._ark)*y4**3)*y1+((-y6/2._ark+sqrt(3._ark)*y5/2._ark)*y3+(sqrt(3._ark)*y5/2._ark+y6/ &
      2._ark)*y4)*y2**3+((-sqrt(3._ark)*y5/2._ark+y6/2._ark)*y3**3+(-y6/2._ark-sqrt(3._ark)*y5/ &
      2._ark)*y4**3)*y2
    dF(3,277) = y1**4*y9+y3**4*y9+y4**4*y9+y2**4*y9
    dF(3,278) = (y8+y7)*y1**4+(-y7-y8)*y2**4+(y8-y7)*y3**4+(-y8+y7)*y4**4
    dF(3,279) = y1**4*y5-y3**4*y5-y4**4*y5+y2**4*y5
    dF(3,280) = y1**4*y2+y1*y2**4-y3*y4**4-y3**4*y4
    dF(3,281) = y1**5-y3**5+y2**5-y4**5
    dF(3,282) = y7**5*y8+y7*y8**5
    dF(3,283) = y7**3*y8**3
    dF(3,284) = y7*y8**3*y9**2+y7**3*y8*y9**2
    dF(3,285) = y7*y8*y9**4
    dF(3,286) = (sqrt(3._ark)*y7**2*y9**3+sqrt(3._ark)*y8**2*y9**3)*y5+(-y8**2*y9**3+ &
      y7**2*y9**3)*y6
    dF(3,287) = (2._ark*y7**4*y9-2._ark*y8**4*y9)*y6
    dF(3,288) = -sqrt(3._ark)*y5*y6**2*y9**3/3._ark+sqrt(3._ark)*y5**3*y9**3
    dF(3,289) = y5*y7**2*y8**2*y9
    dF(3,290) = (y8**2*y9**3+y7**2*y9**3)*y5+(sqrt(3._ark)*y7**2*y9**3- &
      sqrt(3._ark)*y8**2*y9**3)*y6
    dF(3,291) = (y8**4*y9+y7**4*y9)*y5+(-sqrt(3._ark)*y8**4*y9+sqrt(3._ark)*y7**4*y9)*y6
    dF(3,292) = -2._ark*y5*y9**5
    dF(3,293) = (-2._ark/3._ark*sqrt(3._ark)*y7*y8**3-2._ark/3._ark*sqrt(3._ark)*y7**3*y8)*y6**2
    dF(3,294) = sqrt(3._ark)*y6**2*y7*y8*y9**2/6._ark-sqrt(3._ark)*y5**2*y7*y8*y9**2/2._ark
    dF(3,295) = (y7**3*y8-y7*y8**3)*y6*y5+(-sqrt(3._ark)*y7*y8**3/3._ark- &
      sqrt(3._ark)*y7**3*y8/3._ark)*y6**2
    dF(3,296) = (2._ark/3._ark*y8**2*y9+2._ark/3._ark*y7**2*y9)*y5**3+(2._ark*y8**2*y9+ &
      2._ark*y7**2*y9)*y6**2*y5+(4._ark/9._ark*sqrt(3._ark)*y7**2*y9-4._ark/ &
      9._ark*sqrt(3._ark)*y8**2*y9)*y6**3
    dF(3,297) = (-y8**2*y9-y7**2*y9)*y5**3+(-y8**2*y9-y7**2*y9)*y6**2*y5+(-4._ark/ &
      9._ark*sqrt(3._ark)*y7**2*y9+4._ark/9._ark*sqrt(3._ark)*y8**2*y9)*y6**3
    dF(3,298) = y6**2*y7*y8*y9**2+y5**2*y7*y8*y9**2
    dF(3,299) = (y7*y8**3+y7**3*y8)*y5**2+(y7*y8**3+y7**3*y8)*y6**2
    dF(3,300) = (-5._ark/9._ark*sqrt(3._ark)*y7**2*y9-5._ark/ &
      9._ark*sqrt(3._ark)*y8**2*y9)*y5**3+(y7**2*y9-y8**2*y9)*y6*y5**2+(- &
      sqrt(3._ark)*y8**2*y9-sqrt(3._ark)*y7**2*y9)*y6**2*y5+(y8**2*y9/3._ark-y7**2*y9/ &
      3._ark)*y6**3
    dF(3,301) = (-5._ark/9._ark*sqrt(3._ark)*y7**2*y9-5._ark/ &
      9._ark*sqrt(3._ark)*y8**2*y9)*y5**3+(-sqrt(3._ark)*y8**2*y9- &
      sqrt(3._ark)*y7**2*y9)*y6**2*y5
    dF(3,302) = y5**2*y6**2*y7*y8+y6**4*y7*y8/3._ark
    dF(3,303) = 10._ark/3._ark*sqrt(3._ark)*y5*y6**4*y9+6._ark*sqrt(3._ark)*y5**3*y6**2*y9
    dF(3,304) = y5**3*y9**3-3._ark*y5*y6**2*y9**3
    dF(3,305) = -4._ark/9._ark*sqrt(3._ark)*y6**4*y7*y8
    dF(3,306) = -6._ark*y5*y6**4*y9-14._ark*y5**3*y6**2*y9
    dF(3,307) = y5**4*y7*y8+y6**4*y7*y8/3._ark
    dF(3,308) = y5**5*y9-15._ark*y5*y6**4*y9-30._ark*y5**3*y6**2*y9
    dF(3,309) = (-y8*y9**4-y7*y9**4)*y1+(y7*y9**4+y8*y9**4)*y2+(y7*y9**4- &
      y8*y9**4)*y3+(y8*y9**4-y7*y9**4)*y4
    dF(3,310) = (2._ark*y8**2*y9**2-2._ark*y7**2*y9**2)*y6*y1+(2._ark*y8**2*y9**2- &
      2._ark*y7**2*y9**2)*y6*y2+(2._ark*y7**2*y9**2-2._ark*y8**2*y9**2)*y6*y3+ &
      (2._ark*y7**2*y9**2-2._ark*y8**2*y9**2)*y6*y4
    dF(3,311) = ((-sqrt(3._ark)*y7*y8**3/2._ark-sqrt(3._ark)*y7**3*y8/2._ark)*y5+(-y7*y8**3/ &
      2._ark+y7**3*y8/2._ark)*y6)*y1+((-sqrt(3._ark)*y7*y8**3/2._ark-sqrt(3._ark)*y7**3*y8/ &
      2._ark)*y5+(-y7*y8**3/2._ark+y7**3*y8/2._ark)*y6)*y2+((-sqrt(3._ark)*y7*y8**3/2._ark- &
      sqrt(3._ark)*y7**3*y8/2._ark)*y5+(-y7*y8**3/2._ark+y7**3*y8/2._ark)*y6)*y3+((- &
      sqrt(3._ark)*y7*y8**3/2._ark-sqrt(3._ark)*y7**3*y8/2._ark)*y5+(-y7*y8**3/2._ark+y7**3*y8/ &
      2._ark)*y6)*y4
    dF(3,312) = (y8*y9**3+y7*y9**3)*y5*y1+(-y7*y9**3-y8*y9**3)*y5*y2+(y7*y9**3- &
      y8*y9**3)*y5*y3+(y8*y9**3-y7*y9**3)*y5*y4
    dF(3,313) = ((-y7*y8**3/2._ark-y7**3*y8/2._ark)*y5+(sqrt(3._ark)*y7*y8**3/2._ark- &
      sqrt(3._ark)*y7**3*y8/2._ark)*y6)*y1+((-y7*y8**3/2._ark-y7**3*y8/2._ark)*y5+ &
      (sqrt(3._ark)*y7*y8**3/2._ark-sqrt(3._ark)*y7**3*y8/2._ark)*y6)*y2+((-y7*y8**3/2._ark- &
      y7**3*y8/2._ark)*y5+(sqrt(3._ark)*y7*y8**3/2._ark-sqrt(3._ark)*y7**3*y8/2._ark)*y6)*y3+((- &
      y7*y8**3/2._ark-y7**3*y8/2._ark)*y5+(sqrt(3._ark)*y7*y8**3/2._ark-sqrt(3._ark)*y7**3*y8/ &
      2._ark)*y6)*y4
    dF(3,314) = ((-sqrt(3._ark)*y7**2*y8/6._ark-sqrt(3._ark)*y7*y8**2/6._ark)*y5**2+ &
      (y7*y8**2-y7**2*y8)*y6*y5+(-sqrt(3._ark)*y7*y8**2/2._ark-sqrt(3._ark)*y7**2*y8/ &
      2._ark)*y6**2)*y1+((sqrt(3._ark)*y7**2*y8/6._ark+sqrt(3._ark)*y7*y8**2/6._ark)*y5**2+ &
      (y7**2*y8-y7*y8**2)*y6*y5+(sqrt(3._ark)*y7*y8**2/2._ark+sqrt(3._ark)*y7**2*y8/ &
      2._ark)*y6**2)*y2+((sqrt(3._ark)*y7*y8**2/6._ark-sqrt(3._ark)*y7**2*y8/6._ark)*y5**2+(- &
      y7**2*y8-y7*y8**2)*y6*y5+(sqrt(3._ark)*y7*y8**2/2._ark-sqrt(3._ark)*y7**2*y8/ &
      2._ark)*y6**2)*y3+((sqrt(3._ark)*y7**2*y8/6._ark-sqrt(3._ark)*y7*y8**2/6._ark)*y5**2+ &
      (y7*y8**2+y7**2*y8)*y6*y5+(sqrt(3._ark)*y7**2*y8/2._ark-sqrt(3._ark)*y7*y8**2/ &
      2._ark)*y6**2)*y4
    dF(3,315) = ((-y7*y9**2-y8*y9**2)*y5**2+(-y7*y9**2-y8*y9**2)*y6**2)*y1+ &
      ((y8*y9**2+y7*y9**2)*y5**2+(y8*y9**2+y7*y9**2)*y6**2)*y2+((y7*y9**2- &
      y8*y9**2)*y5**2+(y7*y9**2-y8*y9**2)*y6**2)*y3+((y8*y9**2-y7*y9**2)*y5**2+ &
      (y8*y9**2-y7*y9**2)*y6**2)*y4
    dF(3,316) = ((y8**2*y9+y7**2*y9)*y5**2+(y8**2*y9+y7**2*y9)*y6**2)*y1+((y8**2*y9+ &
      y7**2*y9)*y5**2+(y8**2*y9+y7**2*y9)*y6**2)*y2+((y8**2*y9+y7**2*y9)*y5**2+ &
      (y8**2*y9+y7**2*y9)*y6**2)*y3+((y8**2*y9+y7**2*y9)*y5**2+(y8**2*y9+ &
      y7**2*y9)*y6**2)*y4
    dF(3,317) = ((9._ark/2._ark*y7+9._ark/2._ark*y8)*y5**4+(-4._ark*sqrt(3._ark)*y8+ &
      4._ark*sqrt(3._ark)*y7)*y6*y5**3+(y8+y7)*y6**2*y5**2+(-4._ark*sqrt(3._ark)*y8+ &
      4._ark*sqrt(3._ark)*y7)*y6**3*y5+(9._ark/2._ark*y7+9._ark/2._ark*y8)*y6**4)*y1+((-9._ark/ &
      2._ark*y7-9._ark/2._ark*y8)*y5**4+(-4._ark*sqrt(3._ark)*y7+4._ark*sqrt(3._ark)*y8)*y6*y5**3+(- &
      y7-y8)*y6**2*y5**2+(-4._ark*sqrt(3._ark)*y7+4._ark*sqrt(3._ark)*y8)*y6**3*y5+(-9._ark/ &
      2._ark*y7-9._ark/2._ark*y8)*y6**4)*y2+((9._ark/2._ark*y8-9._ark/2._ark*y7)*y5**4+(- &
      4._ark*sqrt(3._ark)*y7-4._ark*sqrt(3._ark)*y8)*y6*y5**3+(y8-y7)*y6**2*y5**2+(- &
      4._ark*sqrt(3._ark)*y7-4._ark*sqrt(3._ark)*y8)*y6**3*y5+(9._ark/2._ark*y8-9._ark/ &
      2._ark*y7)*y6**4)*y3+((9._ark/2._ark*y7-9._ark/2._ark*y8)*y5**4+(4._ark*sqrt(3._ark)*y7+ &
      4._ark*sqrt(3._ark)*y8)*y6*y5**3+(-y8+y7)*y6**2*y5**2+(4._ark*sqrt(3._ark)*y7+ &
      4._ark*sqrt(3._ark)*y8)*y6**3*y5+(9._ark/2._ark*y7-9._ark/2._ark*y8)*y6**4)*y4
    dF(3,318) = y4**2*y7*y8*y9**2+y1**2*y7*y8*y9**2+y3**2*y7*y8*y9**2+ &
      y2**2*y7*y8*y9**2
    dF(3,319) = (2._ark*y7*y8**2-2._ark*y7**2*y8)*y6*y1**2+(2._ark*y7**2*y8- &
      2._ark*y7*y8**2)*y6*y2**2+(-2._ark*y7*y8**2-2._ark*y7**2*y8)*y6*y3**2+(2._ark*y7*y8**2+ &
      2._ark*y7**2*y8)*y6*y4**2
    dF(3,320) = (sqrt(3._ark)*y5*y6**2*y9/3._ark-sqrt(3._ark)*y5**3*y9)*y1**2+ &
      (sqrt(3._ark)*y5*y6**2*y9/3._ark-sqrt(3._ark)*y5**3*y9)*y2**2+(sqrt(3._ark)*y5*y6**2*y9/ &
      3._ark-sqrt(3._ark)*y5**3*y9)*y3**2+(sqrt(3._ark)*y5*y6**2*y9/3._ark- &
      sqrt(3._ark)*y5**3*y9)*y4**2
    dF(3,321) = (y5**3*y9-3._ark*y9*y6**2*y5)*y1**2+(y5**3*y9-3._ark*y9*y6**2*y5)*y2**2+ &
      (y5**3*y9-3._ark*y9*y6**2*y5)*y3**2+(y5**3*y9-3._ark*y9*y6**2*y5)*y4**2
    dF(3,322) = (y6**2*y9+y5**2*y9)*y1**3+(y6**2*y9+y5**2*y9)*y2**3+(y6**2*y9+ &
      y5**2*y9)*y3**3+(y6**2*y9+y5**2*y9)*y4**3
    dF(3,323) = ((-y7-y8)*y5**2+(3._ark*y8+3._ark*y7)*y6**2)*y1**3+((y8+y7)*y5**2+(- &
      3._ark*y7-3._ark*y8)*y6**2)*y2**3+((-y8+y7)*y5**2+(3._ark*y8-3._ark*y7)*y6**2)*y3**3+ &
      ((y8-y7)*y5**2+(-3._ark*y8+3._ark*y7)*y6**2)*y4**3
    dF(3,324) = ((y7/2._ark+y8/2._ark)*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
      2._ark)*y6)*y1**4+((-y8/2._ark-y7/2._ark)*y5+(sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y7/ &
      2._ark)*y6)*y2**4+((y8/2._ark-y7/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
      2._ark)*y6)*y3**4+((-y8/2._ark+y7/2._ark)*y5+(sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/ &
      2._ark)*y6)*y4**4
    dF(3,325) = -2._ark/3._ark*sqrt(3._ark)*y2**4*y6**2-2._ark/3._ark*sqrt(3._ark)*y1**4*y6**2+ &
      2._ark/3._ark*sqrt(3._ark)*y4**4*y6**2+2._ark/3._ark*sqrt(3._ark)*y3**4*y6**2
    dF(3,326) = (-y6**2-y5**2)*y1**4+(-y6**2-y5**2)*y2**4+(y6**2+y5**2)*y3**4+ &
      (y6**2+y5**2)*y4**4
    dF(3,327) = y4**5*y9+y1**5*y9+y2**5*y9+y3**5*y9
    dF(3,328) = (-y7-y8)*y1**5+(y8+y7)*y2**5+(-y8+y7)*y3**5+(y8-y7)*y4**5
    dF(3,329) = (-y7**5-y8**5)*y1+(y7**5+y8**5)*y2+(-y8**5+y7**5)*y3+(-y7**5+ &
      y8**5)*y4
    dF(3,330) = (-y8**3*y9**2-y7**3*y9**2)*y1+(y7**3*y9**2+y8**3*y9**2)*y2+(- &
      y8**3*y9**2+y7**3*y9**2)*y3+(y8**3*y9**2-y7**3*y9**2)*y4
    dF(3,331) = (y7*y8**2*y9**2+y7**2*y8*y9**2)*y1+(-y7**2*y8*y9**2- &
      y7*y8**2*y9**2)*y2+(-y7*y8**2*y9**2+y7**2*y8*y9**2)*y3+(-y7**2*y8*y9**2+ &
      y7*y8**2*y9**2)*y4
    dF(3,332) = -2._ark/3._ark*sqrt(3._ark)*y4*y5*y7*y8*y9**2-2._ark/ &
      3._ark*sqrt(3._ark)*y3*y5*y7*y8*y9**2-2._ark/3._ark*sqrt(3._ark)*y1*y5*y7*y8*y9**2-2._ark/ &
      3._ark*sqrt(3._ark)*y2*y5*y7*y8*y9**2
    dF(3,333) = ((3._ark*sqrt(3._ark)*y7*y9+3._ark*sqrt(3._ark)*y8*y9)*y5**3+(9._ark*y7*y9- &
      9._ark*y8*y9)*y6*y5**2+(3._ark*sqrt(3._ark)*y7*y9+3._ark*sqrt(3._ark)*y8*y9)*y6**2*y5+(- &
      y8*y9+y7*y9)*y6**3)*y1+((-3._ark*sqrt(3._ark)*y7*y9-3._ark*sqrt(3._ark)*y8*y9)*y5**3+(- &
      9._ark*y7*y9+9._ark*y8*y9)*y6*y5**2+(-3._ark*sqrt(3._ark)*y7*y9- &
      3._ark*sqrt(3._ark)*y8*y9)*y6**2*y5+(-y7*y9+y8*y9)*y6**3)*y2+((- &
      3._ark*sqrt(3._ark)*y8*y9+3._ark*sqrt(3._ark)*y7*y9)*y5**3+(9._ark*y7*y9+ &
      9._ark*y8*y9)*y6*y5**2+(-3._ark*sqrt(3._ark)*y8*y9+3._ark*sqrt(3._ark)*y7*y9)*y6**2*y5+ &
      (y7*y9+y8*y9)*y6**3)*y3+((3._ark*sqrt(3._ark)*y8*y9-3._ark*sqrt(3._ark)*y7*y9)*y5**3+(- &
      9._ark*y7*y9-9._ark*y8*y9)*y6*y5**2+(3._ark*sqrt(3._ark)*y8*y9- &
      3._ark*sqrt(3._ark)*y7*y9)*y6**2*y5+(-y7*y9-y8*y9)*y6**3)*y4
    dF(3,334) = (-y8**4-y7**4)*y5*y1+(-y8**4-y7**4)*y5*y2+(y7**4+y8**4)*y5*y3+ &
      (y7**4+y8**4)*y5*y4
    dF(3,335) = ((y8**3-y7**3)*y6*y5+(-sqrt(3._ark)*y7**3/3._ark-sqrt(3._ark)*y8**3/ &
      3._ark)*y6**2)*y1+((y7**3-y8**3)*y6*y5+(sqrt(3._ark)*y8**3/3._ark+sqrt(3._ark)*y7**3/ &
      3._ark)*y6**2)*y2+((y7**3+y8**3)*y6*y5+(-sqrt(3._ark)*y8**3/3._ark+sqrt(3._ark)*y7**3/ &
      3._ark)*y6**2)*y3+((-y7**3-y8**3)*y6*y5+(-sqrt(3._ark)*y7**3/3._ark+sqrt(3._ark)*y8**3/ &
      3._ark)*y6**2)*y4
    dF(3,336) = ((sqrt(3._ark)*y7*y9**2/2._ark+sqrt(3._ark)*y8*y9**2/2._ark)*y5**2+ &
      (y7*y9**2-y8*y9**2)*y6*y5+(sqrt(3._ark)*y7*y9**2/6._ark+sqrt(3._ark)*y8*y9**2/ &
      6._ark)*y6**2)*y1+((-sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7*y9**2/2._ark)*y5**2+ &
      (y8*y9**2-y7*y9**2)*y6*y5+(-sqrt(3._ark)*y8*y9**2/6._ark-sqrt(3._ark)*y7*y9**2/ &
      6._ark)*y6**2)*y2+((sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7*y9**2/2._ark)*y5**2+(- &
      y7*y9**2-y8*y9**2)*y6*y5+(sqrt(3._ark)*y8*y9**2/6._ark-sqrt(3._ark)*y7*y9**2/ &
      6._ark)*y6**2)*y3+((sqrt(3._ark)*y7*y9**2/2._ark-sqrt(3._ark)*y8*y9**2/2._ark)*y5**2+ &
      (y8*y9**2+y7*y9**2)*y6*y5+(sqrt(3._ark)*y7*y9**2/6._ark-sqrt(3._ark)*y8*y9**2/ &
      6._ark)*y6**2)*y4
    dF(3,337) = (-y5*y6**2*y9**2+y5**3*y9**2/3._ark)*y1+(-y5*y6**2*y9**2+y5**3*y9**2/ &
      3._ark)*y2+(y5*y6**2*y9**2-y5**3*y9**2/3._ark)*y3+(y5*y6**2*y9**2-y5**3*y9**2/ &
      3._ark)*y4
    dF(3,338) = (-y5**3*y7*y8/3._ark+y5*y6**2*y7*y8)*y1+(-y5**3*y7*y8/3._ark+ &
      y5*y6**2*y7*y8)*y2+(-y5**3*y7*y8/3._ark+y5*y6**2*y7*y8)*y3+(-y5**3*y7*y8/3._ark+ &
      y5*y6**2*y7*y8)*y4
    dF(3,339) = ((21._ark/10._ark*sqrt(3._ark)*y7+21._ark/10._ark*sqrt(3._ark)*y8)*y5**4+(27._ark/ &
      5._ark*y7-27._ark/5._ark*y8)*y6*y5**3+(23._ark/5._ark*y7-23._ark/5._ark*y8)*y6**3*y5+(19._ark/ &
      10._ark*sqrt(3._ark)*y7+19._ark/10._ark*sqrt(3._ark)*y8)*y6**4)*y1+((-21._ark/ &
      10._ark*sqrt(3._ark)*y7-21._ark/10._ark*sqrt(3._ark)*y8)*y5**4+(-27._ark/5._ark*y7+27._ark/ &
      5._ark*y8)*y6*y5**3+(-23._ark/5._ark*y7+23._ark/5._ark*y8)*y6**3*y5+(-19._ark/ &
      10._ark*sqrt(3._ark)*y7-19._ark/10._ark*sqrt(3._ark)*y8)*y6**4)*y2+((21._ark/ &
      10._ark*sqrt(3._ark)*y8-21._ark/10._ark*sqrt(3._ark)*y7)*y5**4+(-27._ark/5._ark*y7-27._ark/ &
      5._ark*y8)*y6*y5**3+(-23._ark/5._ark*y7-23._ark/5._ark*y8)*y6**3*y5+(-19._ark/ &
      10._ark*sqrt(3._ark)*y7+19._ark/10._ark*sqrt(3._ark)*y8)*y6**4)*y3+((21._ark/ &
      10._ark*sqrt(3._ark)*y7-21._ark/10._ark*sqrt(3._ark)*y8)*y5**4+(27._ark/5._ark*y7+27._ark/ &
      5._ark*y8)*y6*y5**3+(23._ark/5._ark*y7+23._ark/5._ark*y8)*y6**3*y5+(19._ark/ &
      10._ark*sqrt(3._ark)*y7-19._ark/10._ark*sqrt(3._ark)*y8)*y6**4)*y4
    dF(3,340) = (-sqrt(3._ark)*y6**4*y9/8._ark+5._ark/24._ark*sqrt(3._ark)*y5**4*y9-7._ark/ &
      12._ark*sqrt(3._ark)*y5**2*y6**2*y9)*y1+(-sqrt(3._ark)*y6**4*y9/8._ark+5._ark/ &
      24._ark*sqrt(3._ark)*y5**4*y9-7._ark/12._ark*sqrt(3._ark)*y5**2*y6**2*y9)*y2+(- &
      sqrt(3._ark)*y6**4*y9/8._ark+5._ark/24._ark*sqrt(3._ark)*y5**4*y9-7._ark/ &
      12._ark*sqrt(3._ark)*y5**2*y6**2*y9)*y3+(-sqrt(3._ark)*y6**4*y9/8._ark+5._ark/ &
      24._ark*sqrt(3._ark)*y5**4*y9-7._ark/12._ark*sqrt(3._ark)*y5**2*y6**2*y9)*y4
    dF(3,341) = (14._ark*y5**3*y6**2+6._ark*y5*y6**4)*y1+(14._ark*y5**3*y6**2+ &
      6._ark*y5*y6**4)*y2+(-14._ark*y5**3*y6**2-6._ark*y5*y6**4)*y3+(-14._ark*y5**3*y6**2- &
      6._ark*y5*y6**4)*y4
    dF(3,342) = y3*y4*y9**4-y1*y2*y9**4
    dF(3,343) = ((sqrt(3._ark)*y5*y7**2*y9-y6*y7**2*y9)*y3+(y6*y8**2*y9+ &
      sqrt(3._ark)*y5*y8**2*y9)*y4)*y1+((y6*y8**2*y9+sqrt(3._ark)*y5*y8**2*y9)*y3+ &
      (sqrt(3._ark)*y5*y7**2*y9-y6*y7**2*y9)*y4)*y2
    dF(3,344) = ((y6*y7**2*y8/2._ark+sqrt(3._ark)*y5*y7**2*y8/2._ark)*y3+ &
      (sqrt(3._ark)*y5*y7*y8**2/2._ark-y6*y7*y8**2/2._ark)*y4)*y1+((-sqrt(3._ark)*y5*y7*y8**2/ &
      2._ark+y6*y7*y8**2/2._ark)*y3+(-sqrt(3._ark)*y5*y7**2*y8/2._ark-y6*y7**2*y8/2._ark)*y4)*y2
    dF(3,345) = ((sqrt(3._ark)*y5*y9**3/2._ark+y6*y9**3/2._ark)*y3+(sqrt(3._ark)*y5*y9**3/ &
      2._ark-y6*y9**3/2._ark)*y4)*y1+((sqrt(3._ark)*y5*y9**3/2._ark-y6*y9**3/2._ark)*y3+ &
      (sqrt(3._ark)*y5*y9**3/2._ark+y6*y9**3/2._ark)*y4)*y2
    dF(3,346) = ((sqrt(3._ark)*y5*y8**3/2._ark+y6*y8**3/2._ark)*y3+(-y6*y7**3/2._ark+ &
      sqrt(3._ark)*y5*y7**3/2._ark)*y4)*y1+((-sqrt(3._ark)*y5*y7**3/2._ark+y6*y7**3/2._ark)*y3+ &
      (-y6*y8**3/2._ark-sqrt(3._ark)*y5*y8**3/2._ark)*y4)*y2
    dF(3,347) = ((-y6**2*y7*y9-y5**2*y7*y9)*y3+(-y5**2*y8*y9-y6**2*y8*y9)*y4)*y1+ &
      ((y6**2*y8*y9+y5**2*y8*y9)*y3+(y6**2*y7*y9+y5**2*y7*y9)*y4)*y2
    dF(3,348) = (y5**2*y7*y8+y6**2*y7*y8)*y2*y1+(y5**2*y7*y8+y6**2*y7*y8)*y4*y3
    dF(3,349) = ((-sqrt(3._ark)*y6*y8*y9**2/2._ark+y5*y8*y9**2/2._ark)*y3+(y5*y7*y9**2/ &
      2._ark+sqrt(3._ark)*y6*y7*y9**2/2._ark)*y4)*y1+((-y5*y7*y9**2/2._ark- &
      sqrt(3._ark)*y6*y7*y9**2/2._ark)*y3+(-y5*y8*y9**2/2._ark+sqrt(3._ark)*y6*y8*y9**2/ &
      2._ark)*y4)*y2
    dF(3,350) = (-2._ark*y4*y5*y7**2*y9-2._ark*y3*y5*y8**2*y9)*y1+(-2._ark*y3*y5*y7**2*y9- &
      2._ark*y4*y5*y8**2*y9)*y2
    dF(3,351) = ((sqrt(3._ark)*y6*y9**3/2._ark-y5*y9**3/2._ark)*y3+(-sqrt(3._ark)*y6*y9**3/ &
      2._ark-y5*y9**3/2._ark)*y4)*y1+((-sqrt(3._ark)*y6*y9**3/2._ark-y5*y9**3/2._ark)*y3+ &
      (sqrt(3._ark)*y6*y9**3/2._ark-y5*y9**3/2._ark)*y4)*y2
    dF(3,352) = ((y5*y8**3/2._ark-sqrt(3._ark)*y6*y8**3/2._ark)*y3+(sqrt(3._ark)*y6*y7**3/ &
      2._ark+y5*y7**3/2._ark)*y4)*y1+((-y5*y7**3/2._ark-sqrt(3._ark)*y6*y7**3/2._ark)*y3+(- &
      y5*y8**3/2._ark+sqrt(3._ark)*y6*y8**3/2._ark)*y4)*y2
    dF(3,353) = ((-y8**2+y7**2)*y6*y5+(sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y7**2/ &
      3._ark)*y6**2)*y2*y1+((-y7**2+y8**2)*y6*y5+(-sqrt(3._ark)*y8**2/3._ark- &
      sqrt(3._ark)*y7**2/3._ark)*y6**2)*y4*y3
    dF(3,354) = ((sqrt(3._ark)*y5**2*y7*y9/3._ark-y5*y6*y7*y9)*y3+ &
      (sqrt(3._ark)*y5**2*y8*y9/3._ark+y5*y6*y8*y9)*y4)*y1+((-sqrt(3._ark)*y5**2*y8*y9/3._ark- &
      y5*y6*y8*y9)*y3+(-sqrt(3._ark)*y5**2*y7*y9/3._ark+y5*y6*y7*y9)*y4)*y2
    dF(3,355) = (y5**4-y5**2*y6**2)*y2*y1+(-y5**4+y5**2*y6**2)*y4*y3
    dF(3,356) = (-2._ark*y3*y6*y7*y8+2._ark*y4*y6*y7*y8)*y1**2+(-2._ark*y3**2*y6*y7*y8+ &
      2._ark*y4**2*y6*y7*y8)*y1+(2._ark*y3*y6*y7*y8-2._ark*y4*y6*y7*y8)*y2**2+ &
      (2._ark*y3**2*y6*y7*y8-2._ark*y4**2*y6*y7*y8)*y2
    dF(3,357) = ((-20._ark/9._ark*y5**2*y6-4._ark/9._ark*sqrt(3._ark)*y5*y6**2)*y3+(-4._ark/ &
      9._ark*sqrt(3._ark)*y5*y6**2+20._ark/9._ark*y5**2*y6)*y4)*y1**2+((4._ark/ &
      9._ark*sqrt(3._ark)*y5*y6**2+20._ark/9._ark*y5**2*y6)*y3**2+(4._ark/ &
      9._ark*sqrt(3._ark)*y5*y6**2-20._ark/9._ark*y5**2*y6)*y4**2)*y1+((-4._ark/ &
      9._ark*sqrt(3._ark)*y5*y6**2+20._ark/9._ark*y5**2*y6)*y3+(-20._ark/9._ark*y5**2*y6-4._ark/ &
      9._ark*sqrt(3._ark)*y5*y6**2)*y4)*y2**2+((4._ark/9._ark*sqrt(3._ark)*y5*y6**2-20._ark/ &
      9._ark*y5**2*y6)*y3**2+(4._ark/9._ark*sqrt(3._ark)*y5*y6**2+20._ark/ &
      9._ark*y5**2*y6)*y4**2)*y2
    dF(3,358) = y2**2*y7**2*y8**2+y1**2*y7**2*y8**2-y4**2*y7**2*y8**2- &
      y3**2*y7**2*y8**2
    dF(3,359) = (-y7**2*y9**2-y8**2*y9**2)*y1**2+(-y7**2*y9**2-y8**2*y9**2)*y2**2+ &
      (y7**2*y9**2+y8**2*y9**2)*y3**2+(y7**2*y9**2+y8**2*y9**2)*y4**2
    dF(3,360) = (-y7**2*y8*y9-y7*y8**2*y9)*y1**2+(y7*y8**2*y9+y7**2*y8*y9)*y2**2+(- &
      y7*y8**2*y9+y7**2*y8*y9)*y3**2+(y7*y8**2*y9-y7**2*y8*y9)*y4**2
    dF(3,361) = (y7*y8**3+y7**3*y8)*y1**2+(y7*y8**3+y7**3*y8)*y2**2+(y7*y8**3+ &
      y7**3*y8)*y3**2+(y7*y8**3+y7**3*y8)*y4**2
    dF(3,362) = (y7**3*y9+y8**3*y9)*y1**2+(-y7**3*y9-y8**3*y9)*y2**2+(-y8**3*y9+ &
      y7**3*y9)*y3**2+(-y7**3*y9+y8**3*y9)*y4**2
    dF(3,363) = (-y5**2*y6**2-y6**4/2._ark-y5**4/2._ark)*y1**2+(-y5**2*y6**2-y6**4/2._ark- &
      y5**4/2._ark)*y2**2+(y5**4/2._ark+y5**2*y6**2+y6**4/2._ark)*y3**2+(y5**4/2._ark+ &
      y5**2*y6**2+y6**4/2._ark)*y4**2
    dF(3,364) = (8._ark*y7+8._ark*y8)*y5**3*y1**2+(-8._ark*y8-8._ark*y7)*y5**3*y2**2+(- &
      8._ark*y7+8._ark*y8)*y5**3*y3**2+(-8._ark*y8+8._ark*y7)*y5**3*y4**2
    dF(3,365) = (-y3*y7*y9**2-y4*y8*y9**2)*y1**2+(y3**2*y7*y9**2+y4**2*y8*y9**2)*y1+ &
      (y3*y8*y9**2+y4*y7*y9**2)*y2**2+(-y3**2*y8*y9**2-y4**2*y7*y9**2)*y2
    dF(3,366) = (y4*y7**2*y8+y3*y7*y8**2)*y1**2+(-y3**2*y7*y8**2-y4**2*y7**2*y8)*y1+ &
      (-y4*y7*y8**2-y3*y7**2*y8)*y2**2+(y4**2*y7*y8**2+y3**2*y7**2*y8)*y2
    dF(3,367) = (-y3**2*y7*y9-y4**2*y8*y9)*y1**2+(y4**2*y7*y9+y3**2*y8*y9)*y2**2
    dF(3,368) = y1**2*y2**2*y7*y8+y3**2*y4**2*y7*y8
    dF(3,369) = (-y8**2-y7**2)*y2**2*y1**2+(y8**2+y7**2)*y4**2*y3**2
    dF(3,370) = ((y6*y9/2._ark+sqrt(3._ark)*y5*y9/2._ark)*y3**2+(-y6*y9/2._ark+ &
      sqrt(3._ark)*y5*y9/2._ark)*y4**2)*y1**2+((-y6*y9/2._ark+sqrt(3._ark)*y5*y9/2._ark)*y3**2+ &
      (y6*y9/2._ark+sqrt(3._ark)*y5*y9/2._ark)*y4**2)*y2**2
    dF(3,371) = ((-y5*y8-sqrt(3._ark)*y6*y8)*y3**2+(-y5*y7+ &
      sqrt(3._ark)*y6*y7)*y4**2)*y1**2+((-sqrt(3._ark)*y6*y7+y5*y7)*y3**2+(y5*y8+ &
      sqrt(3._ark)*y6*y8)*y4**2)*y2**2
    dF(3,372) = ((sqrt(3._ark)*y6*y9/2._ark-y5*y9/2._ark)*y3**2+(-y5*y9/2._ark- &
      sqrt(3._ark)*y6*y9/2._ark)*y4**2)*y1**2+((-y5*y9/2._ark-sqrt(3._ark)*y6*y9/2._ark)*y3**2+ &
      (sqrt(3._ark)*y6*y9/2._ark-y5*y9/2._ark)*y4**2)*y2**2
    dF(3,373) = y3**2*y4**3*y5+y3**3*y4**2*y5-y1**3*y2**2*y5-y1**2*y2**3*y5
    dF(3,374) = (-y4**2-y3**2)*y1**4+(y3**4+y4**4)*y1**2+(-y4**2-y3**2)*y2**4+ &
      (y3**4+y4**4)*y2**2
    dF(3,375) = (-4._ark*y8-4._ark*y7)*y6**2*y1**3+(4._ark*y8+4._ark*y7)*y6**2*y2**3+(- &
      4._ark*y8+4._ark*y7)*y6**2*y3**3+(4._ark*y8-4._ark*y7)*y6**2*y4**3
    dF(3,376) = (-y5**3+3._ark*y5*y6**2)*y1**3+(-y5**3+3._ark*y5*y6**2)*y2**3+(y5**3- &
      3._ark*y5*y6**2)*y3**3+(y5**3-3._ark*y5*y6**2)*y4**3
    dF(3,377) = (y3*y7*y8+y4*y7*y8)*y1**3+(y3**3*y7*y8+y4**3*y7*y8)*y1+(y3*y7*y8+ &
      y4*y7*y8)*y2**3+(y3**3*y7*y8+y4**3*y7*y8)*y2
    dF(3,378) = (-y6**2-y5**2)*y2*y1**3+(-y6**2-y5**2)*y2**3*y1+(y6**2+ &
      y5**2)*y4*y3**3+(y6**2+y5**2)*y4**3*y3
    dF(3,379) = (-2._ark*y3*y5*y7-2._ark*y4*y5*y8)*y1**3+(2._ark*y4**3*y5*y8+ &
      2._ark*y3**3*y5*y7)*y1+(2._ark*y3*y5*y8+2._ark*y4*y5*y7)*y2**3+(-2._ark*y3**3*y5*y8- &
      2._ark*y4**3*y5*y7)*y2
    dF(3,380) = y3**3*y4*y5*y9+y3*y4**3*y5*y9+y1*y2**3*y5*y9+y1**3*y2*y5*y9
    dF(3,381) = ((sqrt(3._ark)*y6**2/6._ark-sqrt(3._ark)*y5**2/2._ark)*y3+(sqrt(3._ark)*y6**2/ &
      6._ark-sqrt(3._ark)*y5**2/2._ark)*y4)*y1**3+((sqrt(3._ark)*y5**2/2._ark-sqrt(3._ark)*y6**2/ &
      6._ark)*y3**3+(sqrt(3._ark)*y5**2/2._ark-sqrt(3._ark)*y6**2/6._ark)*y4**3)*y1+ &
      ((sqrt(3._ark)*y6**2/6._ark-sqrt(3._ark)*y5**2/2._ark)*y3+(sqrt(3._ark)*y6**2/6._ark- &
      sqrt(3._ark)*y5**2/2._ark)*y4)*y2**3+((sqrt(3._ark)*y5**2/2._ark-sqrt(3._ark)*y6**2/ &
      6._ark)*y3**3+(sqrt(3._ark)*y5**2/2._ark-sqrt(3._ark)*y6**2/6._ark)*y4**3)*y2
    dF(3,382) = ((-y6**2-y5**2)*y3+(-y6**2-y5**2)*y4)*y1**3+((y6**2+y5**2)*y3**3+ &
      (y6**2+y5**2)*y4**3)*y1+((-y6**2-y5**2)*y3+(-y6**2-y5**2)*y4)*y2**3+((y6**2+ &
      y5**2)*y3**3+(y6**2+y5**2)*y4**3)*y2
    dF(3,383) = y1**3*y2**2*y9+y3**2*y4**3*y9+y1**2*y2**3*y9+y3**3*y4**2*y9
    dF(3,384) = (-y4**2*y7-y3**2*y8)*y1**3+(-y3**3*y8-y4**3*y7)*y1**2+(y3**2*y7+ &
      y4**2*y8)*y2**3+(y4**3*y8+y3**3*y7)*y2**2
    dF(3,385) = (y8+y7)*y2**2*y1**3+(-y7-y8)*y2**3*y1**2+(y8-y7)*y4**2*y3**3+(-y8+ &
      y7)*y4**3*y3**2
    dF(3,386) = ((-y6+sqrt(3._ark)*y5)*y3**2+(y6+sqrt(3._ark)*y5)*y4**2)*y1**3+((- &
      sqrt(3._ark)*y5+y6)*y3**3+(-y6-sqrt(3._ark)*y5)*y4**3)*y1**2+((y6+ &
      sqrt(3._ark)*y5)*y3**2+(-y6+sqrt(3._ark)*y5)*y4**2)*y2**3+((-y6- &
      sqrt(3._ark)*y5)*y3**3+(-sqrt(3._ark)*y5+y6)*y4**3)*y2**2
    dF(3,387) = -y1**3*y2**3+y3**3*y4**3
    dF(3,388) = (y8**2+y7**2)*y1**4+(y8**2+y7**2)*y2**4+(-y8**2-y7**2)*y3**4+(- &
      y8**2-y7**2)*y4**4
    dF(3,389) = 2._ark/3._ark*sqrt(3._ark)*y2**5*y5-2._ark/3._ark*sqrt(3._ark)*y4**5*y5+2._ark/ &
      3._ark*sqrt(3._ark)*y1**5*y5-2._ark/3._ark*sqrt(3._ark)*y3**5*y5
    dF(3,390) = (-y7**2*y8**3-y7**3*y8**2)*y1+(y7**2*y8**3+y7**3*y8**2)*y2+(- &
      y7**2*y8**3+y7**3*y8**2)*y3+(-y7**3*y8**2+y7**2*y8**3)*y4
    dF(3,391) = (y8**2*y9**3+y7**2*y9**3)*y1+(y8**2*y9**3+y7**2*y9**3)*y2+ &
      (y8**2*y9**3+y7**2*y9**3)*y3+(y8**2*y9**3+y7**2*y9**3)*y4
    dF(3,392) = (y8**4*y9+y7**4*y9)*y1+(y8**4*y9+y7**4*y9)*y2+(y8**4*y9+ &
      y7**4*y9)*y3+(y8**4*y9+y7**4*y9)*y4
    dF(3,393) = (y7*y8**4+y7**4*y8)*y1+(-y7*y8**4-y7**4*y8)*y2+(y7**4*y8- &
      y7*y8**4)*y3+(y7*y8**4-y7**4*y8)*y4
    dF(3,394) = y1*y7*y8*y9**3-y4*y7*y8*y9**3-y3*y7*y8*y9**3+y2*y7*y8*y9**3
    dF(3,395) = (-y7**4+y8**4)*y6*y1+(-y7**4+y8**4)*y6*y2+(-y8**4+y7**4)*y6*y3+(- &
      y8**4+y7**4)*y6*y4
    dF(3,396) = (y8*y9**3-y7*y9**3)*y6*y1+(y7*y9**3-y8*y9**3)*y6*y2+(-y7*y9**3- &
      y8*y9**3)*y6*y3+(y8*y9**3+y7*y9**3)*y6*y4
    dF(3,397) = ((-sqrt(3._ark)*y7**2*y8*y9-sqrt(3._ark)*y7*y8**2*y9)*y5+(-y7*y8**2*y9+ &
      y7**2*y8*y9)*y6)*y1+((sqrt(3._ark)*y7**2*y8*y9+sqrt(3._ark)*y7*y8**2*y9)*y5+ &
      (y7*y8**2*y9-y7**2*y8*y9)*y6)*y2+((sqrt(3._ark)*y7**2*y8*y9- &
      sqrt(3._ark)*y7*y8**2*y9)*y5+(-y7**2*y8*y9-y7*y8**2*y9)*y6)*y3+ &
      ((sqrt(3._ark)*y7*y8**2*y9-sqrt(3._ark)*y7**2*y8*y9)*y5+(y7*y8**2*y9+ &
      y7**2*y8*y9)*y6)*y4
    dF(3,398) = (y6**4*y9+2._ark*y5**2*y6**2*y9+y5**4*y9)*y1+(y6**4*y9+ &
      2._ark*y5**2*y6**2*y9+y5**4*y9)*y2+(y6**4*y9+2._ark*y5**2*y6**2*y9+y5**4*y9)*y3+ &
      (y6**4*y9+2._ark*y5**2*y6**2*y9+y5**4*y9)*y4
    dF(3,399) = ((-18._ark/5._ark*y7-18._ark/5._ark*y8)*y5**4+(-24._ark/5._ark*sqrt(3._ark)*y7+ &
      24._ark/5._ark*sqrt(3._ark)*y8)*y6*y5**3+(-16._ark/5._ark*sqrt(3._ark)*y7+16._ark/ &
      5._ark*sqrt(3._ark)*y8)*y6**3*y5+(-22._ark/5._ark*y7-22._ark/5._ark*y8)*y6**4)*y1+((18._ark/ &
      5._ark*y7+18._ark/5._ark*y8)*y5**4+(-24._ark/5._ark*sqrt(3._ark)*y8+24._ark/ &
      5._ark*sqrt(3._ark)*y7)*y6*y5**3+(-16._ark/5._ark*sqrt(3._ark)*y8+16._ark/ &
      5._ark*sqrt(3._ark)*y7)*y6**3*y5+(22._ark/5._ark*y8+22._ark/5._ark*y7)*y6**4)*y2+((-18._ark/ &
      5._ark*y8+18._ark/5._ark*y7)*y5**4+(24._ark/5._ark*sqrt(3._ark)*y8+24._ark/ &
      5._ark*sqrt(3._ark)*y7)*y6*y5**3+(16._ark/5._ark*sqrt(3._ark)*y7+16._ark/ &
      5._ark*sqrt(3._ark)*y8)*y6**3*y5+(-22._ark/5._ark*y8+22._ark/5._ark*y7)*y6**4)*y3+((-18._ark/ &
      5._ark*y7+18._ark/5._ark*y8)*y5**4+(-24._ark/5._ark*sqrt(3._ark)*y7-24._ark/ &
      5._ark*sqrt(3._ark)*y8)*y6*y5**3+(-16._ark/5._ark*sqrt(3._ark)*y7-16._ark/ &
      5._ark*sqrt(3._ark)*y8)*y6**3*y5+(22._ark/5._ark*y8-22._ark/5._ark*y7)*y6**4)*y4
    dF(3,400) = y1*y5*y9**4+y2*y5*y9**4-y4*y5*y9**4-y3*y5*y9**4
    dF(3,401) = (-2._ark/3._ark*sqrt(3._ark)*y7**2*y9-2._ark/ &
      3._ark*sqrt(3._ark)*y8**2*y9)*y6**2*y1+(-2._ark/3._ark*sqrt(3._ark)*y7**2*y9-2._ark/ &
      3._ark*sqrt(3._ark)*y8**2*y9)*y6**2*y2+(-2._ark/3._ark*sqrt(3._ark)*y7**2*y9-2._ark/ &
      3._ark*sqrt(3._ark)*y8**2*y9)*y6**2*y3+(-2._ark/3._ark*sqrt(3._ark)*y7**2*y9-2._ark/ &
      3._ark*sqrt(3._ark)*y8**2*y9)*y6**2*y4
    dF(3,402) = ((sqrt(3._ark)*y7**2*y8/6._ark+sqrt(3._ark)*y7*y8**2/6._ark)*y5**2+(- &
      sqrt(3._ark)*y7*y8**2/2._ark-sqrt(3._ark)*y7**2*y8/2._ark)*y6**2)*y1+((- &
      sqrt(3._ark)*y7**2*y8/6._ark-sqrt(3._ark)*y7*y8**2/6._ark)*y5**2+(sqrt(3._ark)*y7*y8**2/ &
      2._ark+sqrt(3._ark)*y7**2*y8/2._ark)*y6**2)*y2+((sqrt(3._ark)*y7**2*y8/6._ark- &
      sqrt(3._ark)*y7*y8**2/6._ark)*y5**2+(sqrt(3._ark)*y7*y8**2/2._ark-sqrt(3._ark)*y7**2*y8/ &
      2._ark)*y6**2)*y3+((sqrt(3._ark)*y7*y8**2/6._ark-sqrt(3._ark)*y7**2*y8/6._ark)*y5**2+ &
      (sqrt(3._ark)*y7**2*y8/2._ark-sqrt(3._ark)*y7*y8**2/2._ark)*y6**2)*y4
    dF(3,403) = ((y7*y9**2-y8*y9**2)*y6*y5+(sqrt(3._ark)*y8*y9**2/3._ark+ &
      sqrt(3._ark)*y7*y9**2/3._ark)*y6**2)*y1+((y8*y9**2-y7*y9**2)*y6*y5+(- &
      sqrt(3._ark)*y8*y9**2/3._ark-sqrt(3._ark)*y7*y9**2/3._ark)*y6**2)*y2+((-y7*y9**2- &
      y8*y9**2)*y6*y5+(-sqrt(3._ark)*y7*y9**2/3._ark+sqrt(3._ark)*y8*y9**2/3._ark)*y6**2)*y3+ &
      ((y8*y9**2+y7*y9**2)*y6*y5+(sqrt(3._ark)*y7*y9**2/3._ark-sqrt(3._ark)*y8*y9**2/ &
      3._ark)*y6**2)*y4
    dF(3,404) = ((y8**2*y9-y7**2*y9)*y6*y5+(sqrt(3._ark)*y7**2*y9/3._ark+ &
      sqrt(3._ark)*y8**2*y9/3._ark)*y6**2)*y1+((y8**2*y9-y7**2*y9)*y6*y5+ &
      (sqrt(3._ark)*y7**2*y9/3._ark+sqrt(3._ark)*y8**2*y9/3._ark)*y6**2)*y2+((y8**2*y9- &
      y7**2*y9)*y6*y5+(sqrt(3._ark)*y7**2*y9/3._ark+sqrt(3._ark)*y8**2*y9/3._ark)*y6**2)*y3+ &
      ((y8**2*y9-y7**2*y9)*y6*y5+(sqrt(3._ark)*y7**2*y9/3._ark+sqrt(3._ark)*y8**2*y9/ &
      3._ark)*y6**2)*y4
    dF(3,405) = ((y8**3-y7**3)*y6*y5+(sqrt(3._ark)*y8**3/3._ark+sqrt(3._ark)*y7**3/ &
      3._ark)*y6**2)*y1+((y7**3-y8**3)*y6*y5+(-sqrt(3._ark)*y7**3/3._ark-sqrt(3._ark)*y8**3/ &
      3._ark)*y6**2)*y2+((y7**3+y8**3)*y6*y5+(-sqrt(3._ark)*y7**3/3._ark+sqrt(3._ark)*y8**3/ &
      3._ark)*y6**2)*y3+((-y7**3-y8**3)*y6*y5+(-sqrt(3._ark)*y8**3/3._ark+sqrt(3._ark)*y7**3/ &
      3._ark)*y6**2)*y4
    dF(3,406) = 2._ark/3._ark*sqrt(3._ark)*y4*y5**2*y9**3+2._ark/ &
      3._ark*sqrt(3._ark)*y1*y5**2*y9**3+2._ark/3._ark*sqrt(3._ark)*y3*y5**2*y9**3+2._ark/ &
      3._ark*sqrt(3._ark)*y2*y5**2*y9**3
    dF(3,407) = 2._ark/3._ark*sqrt(3._ark)*y4*y5**2*y7*y8*y9-2._ark/ &
      3._ark*sqrt(3._ark)*y2*y5**2*y7*y8*y9-2._ark/3._ark*sqrt(3._ark)*y1*y5**2*y7*y8*y9+2._ark/ &
      3._ark*sqrt(3._ark)*y3*y5**2*y7*y8*y9
    dF(3,408) = ((9._ark/16._ark*y7**2+9._ark/16._ark*y8**2)*y5**3+(sqrt(3._ark)*y8**2/4._ark- &
      sqrt(3._ark)*y7**2/4._ark)*y6*y5**2+(-7._ark/16._ark*y7**2-7._ark/ &
      16._ark*y8**2)*y6**2*y5)*y1+((9._ark/16._ark*y7**2+9._ark/16._ark*y8**2)*y5**3+ &
      (sqrt(3._ark)*y8**2/4._ark-sqrt(3._ark)*y7**2/4._ark)*y6*y5**2+(-7._ark/16._ark*y7**2-7._ark/ &
      16._ark*y8**2)*y6**2*y5)*y2+((-9._ark/16._ark*y8**2-9._ark/16._ark*y7**2)*y5**3+ &
      (sqrt(3._ark)*y7**2/4._ark-sqrt(3._ark)*y8**2/4._ark)*y6*y5**2+(7._ark/16._ark*y7**2+7._ark/ &
      16._ark*y8**2)*y6**2*y5)*y3+((-9._ark/16._ark*y8**2-9._ark/16._ark*y7**2)*y5**3+ &
      (sqrt(3._ark)*y7**2/4._ark-sqrt(3._ark)*y8**2/4._ark)*y6*y5**2+(7._ark/16._ark*y7**2+7._ark/ &
      16._ark*y8**2)*y6**2*y5)*y4
    dF(3,409) = ((-y7**3-y8**3)*y5**2+(-y7**3-y8**3)*y6**2)*y1+((y7**3+y8**3)*y5**2+ &
      (y7**3+y8**3)*y6**2)*y2+((y7**3-y8**3)*y5**2+(y7**3-y8**3)*y6**2)*y3+((y8**3- &
      y7**3)*y5**2+(y8**3-y7**3)*y6**2)*y4
    dF(3,410) = ((-sqrt(3._ark)*y8**2/4._ark-sqrt(3._ark)*y7**2/4._ark)*y5**3+(-y7**2+ &
      y8**2)*y6*y5**2+(-sqrt(3._ark)*y8**2/4._ark-sqrt(3._ark)*y7**2/4._ark)*y6**2*y5)*y1+((- &
      sqrt(3._ark)*y8**2/4._ark-sqrt(3._ark)*y7**2/4._ark)*y5**3+(-y7**2+y8**2)*y6*y5**2+(- &
      sqrt(3._ark)*y8**2/4._ark-sqrt(3._ark)*y7**2/4._ark)*y6**2*y5)*y2+((sqrt(3._ark)*y7**2/ &
      4._ark+sqrt(3._ark)*y8**2/4._ark)*y5**3+(-y8**2+y7**2)*y6*y5**2+(sqrt(3._ark)*y7**2/ &
      4._ark+sqrt(3._ark)*y8**2/4._ark)*y6**2*y5)*y3+((sqrt(3._ark)*y7**2/4._ark+ &
      sqrt(3._ark)*y8**2/4._ark)*y5**3+(-y8**2+y7**2)*y6*y5**2+(sqrt(3._ark)*y7**2/4._ark+ &
      sqrt(3._ark)*y8**2/4._ark)*y6**2*y5)*y4
    dF(3,411) = 8._ark/9._ark*sqrt(3._ark)*y2*y5**3*y9**2+8._ark/ &
      9._ark*sqrt(3._ark)*y1*y5**3*y9**2-8._ark/9._ark*sqrt(3._ark)*y4*y5**3*y9**2-8._ark/ &
      9._ark*sqrt(3._ark)*y3*y5**3*y9**2
    dF(3,412) = ((-sqrt(3._ark)*y8*y9-sqrt(3._ark)*y7*y9)*y5**3+(4._ark*y8*y9- &
      4._ark*y7*y9)*y6*y5**2+(-sqrt(3._ark)*y8*y9-sqrt(3._ark)*y7*y9)*y6**2*y5)*y1+ &
      ((sqrt(3._ark)*y7*y9+sqrt(3._ark)*y8*y9)*y5**3+(4._ark*y7*y9-4._ark*y8*y9)*y6*y5**2+ &
      (sqrt(3._ark)*y7*y9+sqrt(3._ark)*y8*y9)*y6**2*y5)*y2+((sqrt(3._ark)*y8*y9- &
      sqrt(3._ark)*y7*y9)*y5**3+(-4._ark*y7*y9-4._ark*y8*y9)*y6*y5**2+(sqrt(3._ark)*y8*y9- &
      sqrt(3._ark)*y7*y9)*y6**2*y5)*y3+((sqrt(3._ark)*y7*y9-sqrt(3._ark)*y8*y9)*y5**3+ &
      (4._ark*y7*y9+4._ark*y8*y9)*y6*y5**2+(sqrt(3._ark)*y7*y9- &
      sqrt(3._ark)*y8*y9)*y6**2*y5)*y4
    dF(3,413) = ((11._ark/16._ark*y8**2+11._ark/16._ark*y7**2)*y5**3+(3._ark/ &
      4._ark*sqrt(3._ark)*y8**2-3._ark/4._ark*sqrt(3._ark)*y7**2)*y6*y5**2+(27._ark/16._ark*y8**2+ &
      27._ark/16._ark*y7**2)*y6**2*y5)*y1+((11._ark/16._ark*y8**2+11._ark/16._ark*y7**2)*y5**3+ &
      (3._ark/4._ark*sqrt(3._ark)*y8**2-3._ark/4._ark*sqrt(3._ark)*y7**2)*y6*y5**2+(27._ark/ &
      16._ark*y8**2+27._ark/16._ark*y7**2)*y6**2*y5)*y2+((-11._ark/16._ark*y8**2-11._ark/ &
      16._ark*y7**2)*y5**3+(3._ark/4._ark*sqrt(3._ark)*y7**2-3._ark/ &
      4._ark*sqrt(3._ark)*y8**2)*y6*y5**2+(-27._ark/16._ark*y8**2-27._ark/ &
      16._ark*y7**2)*y6**2*y5)*y3+((-11._ark/16._ark*y8**2-11._ark/16._ark*y7**2)*y5**3+(3._ark/ &
      4._ark*sqrt(3._ark)*y7**2-3._ark/4._ark*sqrt(3._ark)*y8**2)*y6*y5**2+(-27._ark/16._ark*y8**2- &
      27._ark/16._ark*y7**2)*y6**2*y5)*y4
    dF(3,414) = ((19._ark/10._ark*sqrt(3._ark)*y7+19._ark/10._ark*sqrt(3._ark)*y8)*y5**4+(23._ark/ &
      5._ark*y7-23._ark/5._ark*y8)*y6*y5**3+(27._ark/5._ark*y7-27._ark/5._ark*y8)*y6**3*y5+(21._ark/ &
      10._ark*sqrt(3._ark)*y7+21._ark/10._ark*sqrt(3._ark)*y8)*y6**4)*y1+((-19._ark/ &
      10._ark*sqrt(3._ark)*y7-19._ark/10._ark*sqrt(3._ark)*y8)*y5**4+(-23._ark/5._ark*y7+23._ark/ &
      5._ark*y8)*y6*y5**3+(-27._ark/5._ark*y7+27._ark/5._ark*y8)*y6**3*y5+(-21._ark/ &
      10._ark*sqrt(3._ark)*y7-21._ark/10._ark*sqrt(3._ark)*y8)*y6**4)*y2+((-19._ark/ &
      10._ark*sqrt(3._ark)*y7+19._ark/10._ark*sqrt(3._ark)*y8)*y5**4+(-23._ark/5._ark*y7-23._ark/ &
      5._ark*y8)*y6*y5**3+(-27._ark/5._ark*y7-27._ark/5._ark*y8)*y6**3*y5+(21._ark/ &
      10._ark*sqrt(3._ark)*y8-21._ark/10._ark*sqrt(3._ark)*y7)*y6**4)*y3+((19._ark/ &
      10._ark*sqrt(3._ark)*y7-19._ark/10._ark*sqrt(3._ark)*y8)*y5**4+(23._ark/5._ark*y7+23._ark/ &
      5._ark*y8)*y6*y5**3+(27._ark/5._ark*y7+27._ark/5._ark*y8)*y6**3*y5+(21._ark/ &
      10._ark*sqrt(3._ark)*y7-21._ark/10._ark*sqrt(3._ark)*y8)*y6**4)*y4
    dF(3,415) = (-3._ark/8._ark*sqrt(3._ark)*y6**4*y9+sqrt(3._ark)*y5**2*y6**2*y9/4._ark- &
      sqrt(3._ark)*y5**4*y9/24._ark)*y1+(-3._ark/8._ark*sqrt(3._ark)*y6**4*y9+ &
      sqrt(3._ark)*y5**2*y6**2*y9/4._ark-sqrt(3._ark)*y5**4*y9/24._ark)*y2+(-3._ark/ &
      8._ark*sqrt(3._ark)*y6**4*y9+sqrt(3._ark)*y5**2*y6**2*y9/4._ark-sqrt(3._ark)*y5**4*y9/ &
      24._ark)*y3+(-3._ark/8._ark*sqrt(3._ark)*y6**4*y9+sqrt(3._ark)*y5**2*y6**2*y9/4._ark- &
      sqrt(3._ark)*y5**4*y9/24._ark)*y4
    dF(3,416) = ((22._ark/5._ark*y8+22._ark/5._ark*y7)*y5**4+(-16._ark/5._ark*sqrt(3._ark)*y8+ &
      16._ark/5._ark*sqrt(3._ark)*y7)*y6*y5**3+(-24._ark/5._ark*sqrt(3._ark)*y8+24._ark/ &
      5._ark*sqrt(3._ark)*y7)*y6**3*y5+(18._ark/5._ark*y7+18._ark/5._ark*y8)*y6**4)*y1+((-22._ark/ &
      5._ark*y7-22._ark/5._ark*y8)*y5**4+(-16._ark/5._ark*sqrt(3._ark)*y7+16._ark/ &
      5._ark*sqrt(3._ark)*y8)*y6*y5**3+(-24._ark/5._ark*sqrt(3._ark)*y7+24._ark/ &
      5._ark*sqrt(3._ark)*y8)*y6**3*y5+(-18._ark/5._ark*y7-18._ark/5._ark*y8)*y6**4)*y2+((22._ark/ &
      5._ark*y8-22._ark/5._ark*y7)*y5**4+(-16._ark/5._ark*sqrt(3._ark)*y7-16._ark/ &
      5._ark*sqrt(3._ark)*y8)*y6*y5**3+(-24._ark/5._ark*sqrt(3._ark)*y7-24._ark/ &
      5._ark*sqrt(3._ark)*y8)*y6**3*y5+(-18._ark/5._ark*y7+18._ark/5._ark*y8)*y6**4)*y3+((-22._ark/ &
      5._ark*y8+22._ark/5._ark*y7)*y5**4+(16._ark/5._ark*sqrt(3._ark)*y7+16._ark/ &
      5._ark*sqrt(3._ark)*y8)*y6*y5**3+(24._ark/5._ark*sqrt(3._ark)*y8+24._ark/ &
      5._ark*sqrt(3._ark)*y7)*y6**3*y5+(-18._ark/5._ark*y8+18._ark/5._ark*y7)*y6**4)*y4
    dF(3,417) = (y5**5-30._ark*y5**3*y6**2-15._ark*y5*y6**4)*y1+(y5**5- &
      30._ark*y5**3*y6**2-15._ark*y5*y6**4)*y2+(15._ark*y5*y6**4+30._ark*y5**3*y6**2- &
      y5**5)*y3+(15._ark*y5*y6**4+30._ark*y5**3*y6**2-y5**5)*y4
    dF(3,418) = ((y6*y8**2*y9+sqrt(3._ark)*y5*y8**2*y9)*y3+(sqrt(3._ark)*y5*y7**2*y9- &
      y6*y7**2*y9)*y4)*y1+((sqrt(3._ark)*y5*y7**2*y9-y6*y7**2*y9)*y3+(y6*y8**2*y9+ &
      sqrt(3._ark)*y5*y8**2*y9)*y4)*y2
    dF(3,419) = ((y5**3*y9+9._ark/5._ark*y9*y6**2*y5)*y3+(y5**3*y9+9._ark/ &
      5._ark*y9*y6**2*y5)*y4)*y1+((y5**3*y9+9._ark/5._ark*y9*y6**2*y5)*y3+(y5**3*y9+9._ark/ &
      5._ark*y9*y6**2*y5)*y4)*y2
    dF(3,420) = ((3._ark/4._ark*sqrt(3._ark)*y5**2*y6*y8+3._ark/4._ark*sqrt(3._ark)*y6**3*y8- &
      y5**3*y8)*y3+(-3._ark/4._ark*sqrt(3._ark)*y6**3*y7-y5**3*y7-3._ark/ &
      4._ark*sqrt(3._ark)*y5**2*y6*y7)*y4)*y1+((y5**3*y7+3._ark/4._ark*sqrt(3._ark)*y5**2*y6*y7+ &
      3._ark/4._ark*sqrt(3._ark)*y6**3*y7)*y3+(-3._ark/4._ark*sqrt(3._ark)*y6**3*y8-3._ark/ &
      4._ark*sqrt(3._ark)*y5**2*y6*y8+y5**3*y8)*y4)*y2
    dF(3,421) = (y7**3+y8**3)*y2*y1**2+(-y7**3-y8**3)*y2**2*y1+(y8**3- &
      y7**3)*y4*y3**2+(y7**3-y8**3)*y4**2*y3
    dF(3,422) = (-y4*y7**3-y3*y8**3)*y1**2+(-y3**2*y8**3-y4**2*y7**3)*y1+(y3*y7**3+ &
      y4*y8**3)*y2**2+(y4**2*y8**3+y3**2*y7**3)*y2
    dF(3,423) = (y3*y9**3+y4*y9**3)*y1**2+(y4**2*y9**3+y3**2*y9**3)*y1+(y3*y9**3+ &
      y4*y9**3)*y2**2+(y4**2*y9**3+y3**2*y9**3)*y2
    dF(3,424) = (-y7**2*y8-y7*y8**2)*y2*y1**2+(y7*y8**2+y7**2*y8)*y2**2*y1+ &
      (y7*y8**2-y7**2*y8)*y4*y3**2+(y7**2*y8-y7*y8**2)*y4**2*y3
    dF(3,425) = ((-y5*y9**2/2._ark-sqrt(3._ark)*y6*y9**2/2._ark)*y3+(sqrt(3._ark)*y6*y9**2/ &
      2._ark-y5*y9**2/2._ark)*y4)*y1**2+((y5*y9**2/2._ark+sqrt(3._ark)*y6*y9**2/2._ark)*y3**2+(- &
      sqrt(3._ark)*y6*y9**2/2._ark+y5*y9**2/2._ark)*y4**2)*y1+((sqrt(3._ark)*y6*y9**2/2._ark- &
      y5*y9**2/2._ark)*y3+(-y5*y9**2/2._ark-sqrt(3._ark)*y6*y9**2/2._ark)*y4)*y2**2+((- &
      sqrt(3._ark)*y6*y9**2/2._ark+y5*y9**2/2._ark)*y3**2+(y5*y9**2/2._ark+ &
      sqrt(3._ark)*y6*y9**2/2._ark)*y4**2)*y2
    dF(3,426) = ((-sqrt(3._ark)*y6*y8**2-y5*y8**2)*y3+(-y5*y7**2+ &
      sqrt(3._ark)*y6*y7**2)*y4)*y1**2+((sqrt(3._ark)*y6*y8**2+y5*y8**2)*y3**2+(- &
      sqrt(3._ark)*y6*y7**2+y5*y7**2)*y4**2)*y1+((-y5*y7**2+sqrt(3._ark)*y6*y7**2)*y3+(- &
      sqrt(3._ark)*y6*y8**2-y5*y8**2)*y4)*y2**2+((-sqrt(3._ark)*y6*y7**2+y5*y7**2)*y3**2+ &
      (sqrt(3._ark)*y6*y8**2+y5*y8**2)*y4**2)*y2
    dF(3,427) = ((-y6+sqrt(3._ark)*y5)*y3+(y6+sqrt(3._ark)*y5)*y4)*y1**4+((- &
      sqrt(3._ark)*y5+y6)*y3**4+(-y6-sqrt(3._ark)*y5)*y4**4)*y1+((y6+sqrt(3._ark)*y5)*y3+(- &
      y6+sqrt(3._ark)*y5)*y4)*y2**4+((-y6-sqrt(3._ark)*y5)*y3**4+(-sqrt(3._ark)*y5+ &
      y6)*y4**4)*y2
    dF(3,428) = (-2._ark*y3*y5-2._ark*y4*y5)*y1**4+(2._ark*y3**4*y5+2._ark*y4**4*y5)*y1+(- &
      2._ark*y3*y5-2._ark*y4*y5)*y2**4+(2._ark*y3**4*y5+2._ark*y4**4*y5)*y2
    dF(3,429) = (y4*y7*y8**3+y3*y7**3*y8)*y1+(y4*y7**3*y8+y3*y7*y8**3)*y2
    dF(3,430) = (-y4*y8**3*y9-y3*y7**3*y9)*y1+(y4*y7**3*y9+y3*y8**3*y9)*y2
    dF(3,431) = y1*y2*y7*y8*y9**2+y3*y4*y7*y8*y9**2
    dF(3,432) = (-2._ark*y7**2*y9-2._ark*y8**2*y9)*y5*y2*y1+(-2._ark*y7**2*y9- &
      2._ark*y8**2*y9)*y5*y4*y3
    dF(3,433) = ((y5*y7**2*y9-sqrt(3._ark)*y6*y7**2*y9)*y3+(y5*y8**2*y9+ &
      sqrt(3._ark)*y6*y8**2*y9)*y4)*y1+((y5*y8**2*y9+sqrt(3._ark)*y6*y8**2*y9)*y3+ &
      (y5*y7**2*y9-sqrt(3._ark)*y6*y7**2*y9)*y4)*y2
    dF(3,434) = 2._ark/3._ark*sqrt(3._ark)*y1*y2*y6**2*y9**2-2._ark/ &
      3._ark*sqrt(3._ark)*y3*y4*y6**2*y9**2
    dF(3,435) = (((y7*y9**2-y8*y9**2)*y3+(y8*y9**2-y7*y9**2)*y4)*y2+(-y7*y9**2- &
      y8*y9**2)*y4*y3)*y1+(y8*y9**2+y7*y9**2)*y4*y3*y2
    dF(3,436) = (((y7*y8**2-y7**2*y8)*y3+(y7**2*y8-y7*y8**2)*y4)*y2+(-y7**2*y8- &
      y7*y8**2)*y4*y3)*y1+(y7*y8**2+y7**2*y8)*y4*y3*y2
    dF(3,437) = (((y8**2*y9+y7**2*y9)*y3+(y8**2*y9+y7**2*y9)*y4)*y2+(y8**2*y9+ &
      y7**2*y9)*y4*y3)*y1+(y8**2*y9+y7**2*y9)*y4*y3*y2
    dF(3,438) = (((-2._ark*y8**2+2._ark*y7**2)*y6*y3+(-2._ark*y8**2+2._ark*y7**2)*y6*y4)*y2+ &
      (2._ark*y8**2-2._ark*y7**2)*y6*y4*y3)*y1+(2._ark*y8**2-2._ark*y7**2)*y6*y4*y3*y2
    dF(3,439) = y1**3*y3*y4*y9+(y4**3*y9+y3**3*y9)*y2*y1+y2**3*y3*y4*y9
    dF(3,440) = (y4*y8+y3*y7)*y2*y1**3+((-y4*y7-y3*y8)*y2**3-y3**3*y4*y7- &
      y3*y4**3*y8)*y1+(y3*y4**3*y7+y3**3*y4*y8)*y2
    dF(3,441) = (2._ark*y4*y6-2._ark*y3*y6)*y2*y1**3+((2._ark*y3*y6-2._ark*y4*y6)*y2**3- &
      2._ark*y3*y4**3*y6+2._ark*y3**3*y4*y6)*y1+(-2._ark*y3**3*y4*y6+2._ark*y3*y4**3*y6)*y2
    dF(3,442) = ((sqrt(3._ark)*y6-y5)*y3+(-y5-sqrt(3._ark)*y6)*y4)*y2*y1**3+(((-y5- &
      sqrt(3._ark)*y6)*y3+(sqrt(3._ark)*y6-y5)*y4)*y2**3+(y5-sqrt(3._ark)*y6)*y4*y3**3+ &
      (sqrt(3._ark)*y6+y5)*y4**3*y3)*y1+((sqrt(3._ark)*y6+y5)*y4*y3**3+(y5- &
      sqrt(3._ark)*y6)*y4**3*y3)*y2
    dF(3,443) = (2._ark*y4*y6*y7**2-2._ark*y3*y6*y8**2)*y1**2+(-2._ark*y4**2*y6*y7**2+ &
      2._ark*y3**2*y6*y8**2)*y1+(-2._ark*y4*y6*y8**2+2._ark*y3*y6*y7**2)*y2**2+ &
      (2._ark*y4**2*y6*y8**2-2._ark*y3**2*y6*y7**2)*y2
    dF(3,444) = ((sqrt(3._ark)*y7**2+sqrt(3._ark)*y8**2)*y5+(-y7**2+y8**2)*y6)*y2*y1**2+ &
      ((sqrt(3._ark)*y7**2+sqrt(3._ark)*y8**2)*y5+(-y7**2+y8**2)*y6)*y2**2*y1+((- &
      sqrt(3._ark)*y7**2-sqrt(3._ark)*y8**2)*y5+(-y8**2+y7**2)*y6)*y4*y3**2+((- &
      sqrt(3._ark)*y7**2-sqrt(3._ark)*y8**2)*y5+(-y8**2+y7**2)*y6)*y4**2*y3
    dF(3,445) = ((-y6*y9**2/2._ark+sqrt(3._ark)*y5*y9**2/2._ark)*y3+(y6*y9**2/2._ark+ &
      sqrt(3._ark)*y5*y9**2/2._ark)*y4)*y1**2+((y6*y9**2/2._ark-sqrt(3._ark)*y5*y9**2/ &
      2._ark)*y3**2+(-y6*y9**2/2._ark-sqrt(3._ark)*y5*y9**2/2._ark)*y4**2)*y1+((y6*y9**2/2._ark+ &
      sqrt(3._ark)*y5*y9**2/2._ark)*y3+(-y6*y9**2/2._ark+sqrt(3._ark)*y5*y9**2/ &
      2._ark)*y4)*y2**2+((-y6*y9**2/2._ark-sqrt(3._ark)*y5*y9**2/2._ark)*y3**2+(y6*y9**2/2._ark- &
      sqrt(3._ark)*y5*y9**2/2._ark)*y4**2)*y2
    dF(3,446) = ((-y6**2*y8-y5**2*y8)*y3+(-y6**2*y7-y5**2*y7)*y4)*y1**2+((-y6**2*y8- &
      y5**2*y8)*y3**2+(-y6**2*y7-y5**2*y7)*y4**2)*y1+((y6**2*y7+y5**2*y7)*y3+ &
      (y5**2*y8+y6**2*y8)*y4)*y2**2+((y6**2*y7+y5**2*y7)*y3**2+(y5**2*y8+ &
      y6**2*y8)*y4**2)*y2
    dF(3,447) = (5._ark/9._ark*sqrt(3._ark)*y5**3+sqrt(3._ark)*y5*y6**2)*y2*y1**2+(5._ark/ &
      9._ark*sqrt(3._ark)*y5**3+sqrt(3._ark)*y5*y6**2)*y2**2*y1+(-sqrt(3._ark)*y5*y6**2-5._ark/ &
      9._ark*sqrt(3._ark)*y5**3)*y4*y3**2+(-sqrt(3._ark)*y5*y6**2-5._ark/ &
      9._ark*sqrt(3._ark)*y5**3)*y4**2*y3
    dF(3,448) = ((4._ark/3._ark*sqrt(3._ark)*y5**2*y6-y5*y6**2-y5**3)*y3+(-y5*y6**2-y5**3- &
      4._ark/3._ark*sqrt(3._ark)*y5**2*y6)*y4)*y1**2+((y5**3+y5*y6**2-4._ark/ &
      3._ark*sqrt(3._ark)*y5**2*y6)*y3**2+(4._ark/3._ark*sqrt(3._ark)*y5**2*y6+y5*y6**2+ &
      y5**3)*y4**2)*y1+((-y5*y6**2-y5**3-4._ark/3._ark*sqrt(3._ark)*y5**2*y6)*y3+(4._ark/ &
      3._ark*sqrt(3._ark)*y5**2*y6-y5*y6**2-y5**3)*y4)*y2**2+((4._ark/ &
      3._ark*sqrt(3._ark)*y5**2*y6+y5*y6**2+y5**3)*y3**2+(y5**3+y5*y6**2-4._ark/ &
      3._ark*sqrt(3._ark)*y5**2*y6)*y4**2)*y2
    dF(3,449) = (-y3*y9**2-y4*y9**2)*y2*y1**2+((-y3*y9**2-y4*y9**2)*y2**2+ &
      y3**2*y4*y9**2+y3*y4**2*y9**2)*y1+(y3**2*y4*y9**2+y3*y4**2*y9**2)*y2
    dF(3,450) = (y4*y7*y9+y3*y8*y9)*y2*y1**2+((-y4*y8*y9-y3*y7*y9)*y2**2- &
      y3**2*y4*y8*y9-y3*y4**2*y7*y9)*y1+(y3**2*y4*y7*y9+y3*y4**2*y8*y9)*y2
    dF(3,451) = ((-y4**2*y7-y3**2*y8)*y2-y3*y4**2*y7-y3**2*y4*y8)*y1**2+(y3**2*y7+ &
      y4**2*y8)*y2**2*y1+(y3*y4**2*y8+y3**2*y4*y7)*y2**2
    dF(3,452) = ((-y8+y7)*y3+(y8-y7)*y4)*y2**2*y1**2+(-y7-y8)*y4**2*y3**2*y1+(y8+ &
      y7)*y4**2*y3**2*y2
    dF(3,453) = (-y4*y5-y3*y5)*y2**2*y1**2+y1*y3**2*y4**2*y5+y2*y3**2*y4**2*y5
    dF(3,454) = (y3*y9**2+y4*y9**2)*y1**3+(-y4**3*y9**2-y3**3*y9**2)*y1+(y3*y9**2+ &
      y4*y9**2)*y2**3+(-y4**3*y9**2-y3**3*y9**2)*y2
    dF(3,455) = (y4+y3)*y1**5+(-y3**5-y4**5)*y1+(y4+y3)*y2**5+(-y3**5-y4**5)*y2
    dF(3,456) = ((sqrt(3._ark)*y7*y9**2/2._ark+sqrt(3._ark)*y8*y9**2/2._ark)*y5+(-y8*y9**2/ &
      2._ark+y7*y9**2/2._ark)*y6)*y1**2+((-sqrt(3._ark)*y8*y9**2/2._ark-sqrt(3._ark)*y7*y9**2/ &
      2._ark)*y5+(y8*y9**2/2._ark-y7*y9**2/2._ark)*y6)*y2**2+((sqrt(3._ark)*y8*y9**2/2._ark- &
      sqrt(3._ark)*y7*y9**2/2._ark)*y5+(-y8*y9**2/2._ark-y7*y9**2/2._ark)*y6)*y3**2+ &
      ((sqrt(3._ark)*y7*y9**2/2._ark-sqrt(3._ark)*y8*y9**2/2._ark)*y5+(y8*y9**2/2._ark+y7*y9**2/ &
      2._ark)*y6)*y4**2
    dF(3,457) = ((-sqrt(3._ark)*y8**2*y9-sqrt(3._ark)*y7**2*y9)*y5+(y8**2*y9- &
      y7**2*y9)*y6)*y1**2+((-sqrt(3._ark)*y8**2*y9-sqrt(3._ark)*y7**2*y9)*y5+(y8**2*y9- &
      y7**2*y9)*y6)*y2**2+((-sqrt(3._ark)*y8**2*y9-sqrt(3._ark)*y7**2*y9)*y5+(y8**2*y9- &
      y7**2*y9)*y6)*y3**2+((-sqrt(3._ark)*y8**2*y9-sqrt(3._ark)*y7**2*y9)*y5+(y8**2*y9- &
      y7**2*y9)*y6)*y4**2
    dF(3,458) = ((y8**2+y7**2)*y5**2+(y8**2+y7**2)*y6**2)*y1**2+((y8**2+ &
      y7**2)*y5**2+(y8**2+y7**2)*y6**2)*y2**2+((-y8**2-y7**2)*y5**2+(-y8**2- &
      y7**2)*y6**2)*y3**2+((-y8**2-y7**2)*y5**2+(-y8**2-y7**2)*y6**2)*y4**2
    dF(3,459) = (y5**2*y7*y8+y6**2*y7*y8)*y1**2+(y5**2*y7*y8+y6**2*y7*y8)*y2**2+ &
      (y5**2*y7*y8+y6**2*y7*y8)*y3**2+(y5**2*y7*y8+y6**2*y7*y8)*y4**2
    dF(3,460) = (-2._ark*y7**2*y9-2._ark*y8**2*y9)*y5*y1**2+(-2._ark*y7**2*y9- &
      2._ark*y8**2*y9)*y5*y2**2+(-2._ark*y7**2*y9-2._ark*y8**2*y9)*y5*y3**2+(-2._ark*y7**2*y9- &
      2._ark*y8**2*y9)*y5*y4**2
    dF(3,461) = ((-y7**2*y8-y7*y8**2)*y5+(sqrt(3._ark)*y7*y8**2- &
      sqrt(3._ark)*y7**2*y8)*y6)*y1**2+((y7*y8**2+y7**2*y8)*y5+(sqrt(3._ark)*y7**2*y8- &
      sqrt(3._ark)*y7*y8**2)*y6)*y2**2+((y7*y8**2-y7**2*y8)*y5+(-sqrt(3._ark)*y7**2*y8- &
      sqrt(3._ark)*y7*y8**2)*y6)*y3**2+((y7**2*y8-y7*y8**2)*y5+(sqrt(3._ark)*y7*y8**2+ &
      sqrt(3._ark)*y7**2*y8)*y6)*y4**2
    dF(3,462) = ((y8*y9**2/2._ark+y7*y9**2/2._ark)*y5+(sqrt(3._ark)*y8*y9**2/2._ark- &
      sqrt(3._ark)*y7*y9**2/2._ark)*y6)*y1**2+((-y8*y9**2/2._ark-y7*y9**2/2._ark)*y5+ &
      (sqrt(3._ark)*y7*y9**2/2._ark-sqrt(3._ark)*y8*y9**2/2._ark)*y6)*y2**2+((y8*y9**2/2._ark- &
      y7*y9**2/2._ark)*y5+(sqrt(3._ark)*y7*y9**2/2._ark+sqrt(3._ark)*y8*y9**2/2._ark)*y6)*y3**2+ &
      ((-y8*y9**2/2._ark+y7*y9**2/2._ark)*y5+(-sqrt(3._ark)*y8*y9**2/2._ark- &
      sqrt(3._ark)*y7*y9**2/2._ark)*y6)*y4**2
    dF(3,463) = ((-y7*y9+y8*y9)*y6*y5+(sqrt(3._ark)*y8*y9/3._ark+sqrt(3._ark)*y7*y9/ &
      3._ark)*y6**2)*y1**2+((-y8*y9+y7*y9)*y6*y5+(-sqrt(3._ark)*y7*y9/3._ark- &
      sqrt(3._ark)*y8*y9/3._ark)*y6**2)*y2**2+((-y7*y9-y8*y9)*y6*y5+(-sqrt(3._ark)*y8*y9/ &
      3._ark+sqrt(3._ark)*y7*y9/3._ark)*y6**2)*y3**2+((y7*y9+y8*y9)*y6*y5+(- &
      sqrt(3._ark)*y7*y9/3._ark+sqrt(3._ark)*y8*y9/3._ark)*y6**2)*y4**2
    dF(3,464) = (-2._ark/3._ark*sqrt(3._ark)*y8**2-2._ark/ &
      3._ark*sqrt(3._ark)*y7**2)*y5**2*y1**2+(-2._ark/3._ark*sqrt(3._ark)*y8**2-2._ark/ &
      3._ark*sqrt(3._ark)*y7**2)*y5**2*y2**2+(2._ark/3._ark*sqrt(3._ark)*y8**2+2._ark/ &
      3._ark*sqrt(3._ark)*y7**2)*y5**2*y3**2+(2._ark/3._ark*sqrt(3._ark)*y8**2+2._ark/ &
      3._ark*sqrt(3._ark)*y7**2)*y5**2*y4**2
    dF(3,465) = ((-sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y5**2+(-y8**2+ &
      y7**2)*y6*y5)*y1**2+((-sqrt(3._ark)*y8**2/3._ark-sqrt(3._ark)*y7**2/3._ark)*y5**2+(- &
      y8**2+y7**2)*y6*y5)*y2**2+((sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y7**2/3._ark)*y5**2+ &
      (-y7**2+y8**2)*y6*y5)*y3**2+((sqrt(3._ark)*y8**2/3._ark+sqrt(3._ark)*y7**2/ &
      3._ark)*y5**2+(-y7**2+y8**2)*y6*y5)*y4**2
    dF(3,466) = (y8*y9**2+y7*y9**2)*y2*y1**2+(-y7*y9**2-y8*y9**2)*y2**2*y1+ &
      (y8*y9**2-y7*y9**2)*y4*y3**2+(y7*y9**2-y8*y9**2)*y4**2*y3
    dF(3,467) = (-y4*y7*y9**2-y3*y8*y9**2)*y1**2+(-y3**2*y8*y9**2- &
      y4**2*y7*y9**2)*y1+(y4*y8*y9**2+y3*y7*y9**2)*y2**2+(y3**2*y7*y9**2+ &
      y4**2*y8*y9**2)*y2
    dF(3,468) = (y4*y7**2*y9+y3*y8**2*y9)*y1**2+(y3**2*y8**2*y9+y4**2*y7**2*y9)*y1+ &
      (y3*y7**2*y9+y4*y8**2*y9)*y2**2+(y4**2*y8**2*y9+y3**2*y7**2*y9)*y2
    dF(3,469) = -2._ark/3._ark*sqrt(3._ark)*y3**2*y4*y5*y9**2+2._ark/ &
      3._ark*sqrt(3._ark)*y1*y2**2*y5*y9**2+2._ark/3._ark*sqrt(3._ark)*y1**2*y2*y5*y9**2-2._ark/ &
      3._ark*sqrt(3._ark)*y3*y4**2*y5*y9**2
    dF(3,470) = ((-y6*y7**2+sqrt(3._ark)*y5*y7**2)*y3+(sqrt(3._ark)*y5*y8**2+ &
      y6*y8**2)*y4)*y1**2+((y6*y7**2-sqrt(3._ark)*y5*y7**2)*y3**2+(-y6*y8**2- &
      sqrt(3._ark)*y5*y8**2)*y4**2)*y1+((sqrt(3._ark)*y5*y8**2+y6*y8**2)*y3+(-y6*y7**2+ &
      sqrt(3._ark)*y5*y7**2)*y4)*y2**2+((-y6*y8**2-sqrt(3._ark)*y5*y8**2)*y3**2+(y6*y7**2- &
      sqrt(3._ark)*y5*y7**2)*y4**2)*y2
    dF(3,471) = (y3*y6*y8*y9-y4*y6*y7*y9)*y1**2+(y4**2*y6*y7*y9-y3**2*y6*y8*y9)*y1+ &
      (-y4*y6*y8*y9+y3*y6*y7*y9)*y2**2+(y4**2*y6*y8*y9-y3**2*y6*y7*y9)*y2
    dF(3,472) = ((-sqrt(3._ark)*y5*y6*y7-y6**2*y7)*y3+(sqrt(3._ark)*y5*y6*y8- &
      y6**2*y8)*y4)*y1**2+((sqrt(3._ark)*y5*y6*y7+y6**2*y7)*y3**2+(-sqrt(3._ark)*y5*y6*y8+ &
      y6**2*y8)*y4**2)*y1+((-sqrt(3._ark)*y5*y6*y8+y6**2*y8)*y3+(sqrt(3._ark)*y5*y6*y7+ &
      y6**2*y7)*y4)*y2**2+((sqrt(3._ark)*y5*y6*y8-y6**2*y8)*y3**2+(-sqrt(3._ark)*y5*y6*y7- &
      y6**2*y7)*y4**2)*y2
    dF(3,473) = (4._ark*y4*y6**2*y9+4._ark*y3*y6**2*y9)*y1**2+(4._ark*y3**2*y6**2*y9+ &
      4._ark*y4**2*y6**2*y9)*y1+(4._ark*y4*y6**2*y9+4._ark*y3*y6**2*y9)*y2**2+ &
      (4._ark*y3**2*y6**2*y9+4._ark*y4**2*y6**2*y9)*y2
    dF(3,474) = (2._ark*y3*y5*y7**2+2._ark*y4*y5*y8**2)*y1**2+(-2._ark*y3**2*y5*y7**2- &
      2._ark*y4**2*y5*y8**2)*y1+(2._ark*y3*y5*y8**2+2._ark*y4*y5*y7**2)*y2**2+(- &
      2._ark*y3**2*y5*y8**2-2._ark*y4**2*y5*y7**2)*y2
    dF(3,475) = ((y5*y7*y8-sqrt(3._ark)*y6*y7*y8)*y3+(sqrt(3._ark)*y6*y7*y8+ &
      y5*y7*y8)*y4)*y1**2+((y5*y7*y8-sqrt(3._ark)*y6*y7*y8)*y3**2+(sqrt(3._ark)*y6*y7*y8+ &
      y5*y7*y8)*y4**2)*y1+((sqrt(3._ark)*y6*y7*y8+y5*y7*y8)*y3+(y5*y7*y8- &
      sqrt(3._ark)*y6*y7*y8)*y4)*y2**2+((sqrt(3._ark)*y6*y7*y8+y5*y7*y8)*y3**2+(y5*y7*y8- &
      sqrt(3._ark)*y6*y7*y8)*y4**2)*y2
    dF(3,476) = (-y3*y7**2-y4*y8**2)*y2*y1**2+((-y3*y8**2-y4*y7**2)*y2**2+ &
      y3**2*y4*y7**2+y3*y4**2*y8**2)*y1+(y3*y4**2*y7**2+y3**2*y4*y8**2)*y2
    dF(3,477) = (y8**2+y7**2)*y4*y3*y1**2+((-y8**2-y7**2)*y3**2+(-y8**2- &
      y7**2)*y4**2)*y2*y1+(y8**2+y7**2)*y4*y3*y2**2
    dF(3,478) = ((sqrt(3._ark)*y5*y7+y6*y7)*y3+(sqrt(3._ark)*y5*y8-y6*y8)*y4)*y2*y1**2+ &
      (((-sqrt(3._ark)*y5*y8+y6*y8)*y3+(-y6*y7-sqrt(3._ark)*y5*y7)*y4)*y2**2+(-y6*y7- &
      sqrt(3._ark)*y5*y7)*y4*y3**2+(-sqrt(3._ark)*y5*y8+y6*y8)*y4**2*y3)*y1+ &
      ((sqrt(3._ark)*y5*y8-y6*y8)*y4*y3**2+(sqrt(3._ark)*y5*y7+y6*y7)*y4**2*y3)*y2
    dF(3,479) = ((-y6*y8-sqrt(3._ark)*y5*y8)*y3**2+(y6*y7- &
      sqrt(3._ark)*y5*y7)*y4**2)*y1**2+((sqrt(3._ark)*y5*y7-y6*y7)*y3**2+(y6*y8+ &
      sqrt(3._ark)*y5*y8)*y4**2)*y2**2
    dF(3,480) = (-y6**2-y5**2)*y2**2*y1**2+(y6**2+y5**2)*y4**2*y3**2
    dF(3,481) = -2._ark*y3**2*y4**2*y5*y9-2._ark*y1**2*y2**2*y5*y9
    dF(3,482) = ((y3**2*y9+y4**2*y9)*y2+y3*y4**2*y9+y3**2*y4*y9)*y1**2+(y3**2*y9+ &
      y4**2*y9)*y2**2*y1+(y3*y4**2*y9+y3**2*y4*y9)*y2**2
    dF(3,483) = ((y5-sqrt(3._ark)*y6)*y3**2+(sqrt(3._ark)*y6+y5)*y4**2)*y1**3+ &
      ((sqrt(3._ark)*y6-y5)*y3**3+(-y5-sqrt(3._ark)*y6)*y4**3)*y1**2+((sqrt(3._ark)*y6+ &
      y5)*y3**2+(y5-sqrt(3._ark)*y6)*y4**2)*y2**3+((-y5-sqrt(3._ark)*y6)*y3**3+ &
      (sqrt(3._ark)*y6-y5)*y4**3)*y2**2
    dF(3,484) = (-y7**3-y8**3)*y1**3+(y7**3+y8**3)*y2**3+(y7**3-y8**3)*y3**3+(y8**3- &
      y7**3)*y4**3
    dF(3,485) = y4**3*y9**3+y1**3*y9**3+y3**3*y9**3+y2**3*y9**3
    dF(3,486) = (y8**2*y9+y7**2*y9)*y1**3+(y8**2*y9+y7**2*y9)*y2**3+(y8**2*y9+ &
      y7**2*y9)*y3**3+(y8**2*y9+y7**2*y9)*y4**3
    dF(3,487) = (y7*y8**2+y7**2*y8)*y1**3+(-y7**2*y8-y7*y8**2)*y2**3+(y7**2*y8- &
      y7*y8**2)*y3**3+(y7*y8**2-y7**2*y8)*y4**3
    dF(3,488) = (-y7*y9**2-y8*y9**2)*y1**3+(y8*y9**2+y7*y9**2)*y2**3+(y7*y9**2- &
      y8*y9**2)*y3**3+(y8*y9**2-y7*y9**2)*y4**3
    dF(3,489) = y4**3*y7*y8*y9+y3**3*y7*y8*y9-y1**3*y7*y8*y9-y2**3*y7*y8*y9
    dF(3,490) = ((sqrt(3._ark)*y7**2+sqrt(3._ark)*y8**2)*y5+(-y7**2+y8**2)*y6)*y1**3+ &
      ((sqrt(3._ark)*y7**2+sqrt(3._ark)*y8**2)*y5+(-y7**2+y8**2)*y6)*y2**3+((- &
      sqrt(3._ark)*y7**2-sqrt(3._ark)*y8**2)*y5+(-y8**2+y7**2)*y6)*y3**3+((- &
      sqrt(3._ark)*y7**2-sqrt(3._ark)*y8**2)*y5+(-y8**2+y7**2)*y6)*y4**3
    dF(3,491) = (-y8*y9+y7*y9)*y6*y1**3+(-y7*y9+y8*y9)*y6*y2**3+(y7*y9+ &
      y8*y9)*y6*y3**3+(-y7*y9-y8*y9)*y6*y4**3
    dF(3,492) = (-y7*y9-y8*y9)*y5*y1**3+(y7*y9+y8*y9)*y5*y2**3+(-y7*y9+ &
      y8*y9)*y5*y3**3+(-y8*y9+y7*y9)*y5*y4**3
    dF(3,493) = ((-y8**2-y7**2)*y5+(-sqrt(3._ark)*y8**2+sqrt(3._ark)*y7**2)*y6)*y1**3+ &
      ((-y8**2-y7**2)*y5+(-sqrt(3._ark)*y8**2+sqrt(3._ark)*y7**2)*y6)*y2**3+((y8**2+ &
      y7**2)*y5+(sqrt(3._ark)*y8**2-sqrt(3._ark)*y7**2)*y6)*y3**3+((y8**2+y7**2)*y5+ &
      (sqrt(3._ark)*y8**2-sqrt(3._ark)*y7**2)*y6)*y4**3
    dF(3,494) = (-y4*y8*y9-y3*y7*y9)*y1**3+(-y3**3*y7*y9-y4**3*y8*y9)*y1+(y4*y7*y9+ &
      y3*y8*y9)*y2**3+(y3**3*y8*y9+y4**3*y7*y9)*y2
    dF(3,495) = (y3*y7**2+y4*y8**2)*y1**3+(-y4**3*y8**2-y3**3*y7**2)*y1+(y3*y8**2+ &
      y4*y7**2)*y2**3+(-y4**3*y7**2-y3**3*y8**2)*y2
    dF(3,496) = (y3**2*y9+y4**2*y9)*y1**3+(y4**3*y9+y3**3*y9)*y1**2+(y3**2*y9+ &
      y4**2*y9)*y2**3+(y4**3*y9+y3**3*y9)*y2**2
    dF(3,497) = (y3**2*y7+y4**2*y8)*y1**3+(-y4**3*y8-y3**3*y7)*y1**2+(-y4**2*y7- &
      y3**2*y8)*y2**3+(y3**3*y8+y4**3*y7)*y2**2
    dF(3,498) = y4*y9**5+y3*y9**5+y2*y9**5+y1*y9**5
    dF(3,499) = (y7**3*y8*y9+y7*y8**3*y9)*y1+(y7**3*y8*y9+y7*y8**3*y9)*y2+(- &
      y7**3*y8*y9-y7*y8**3*y9)*y3+(-y7**3*y8*y9-y7*y8**3*y9)*y4
    dF(3,500) = y4*y7**2*y8**2*y9+y3*y7**2*y8**2*y9+y2*y7**2*y8**2*y9+ &
      y1*y7**2*y8**2*y9
    dF(3,501) = ((-sqrt(3._ark)*y7**3*y9-sqrt(3._ark)*y8**3*y9)*y5+(-y7**3*y9+ &
      y8**3*y9)*y6)*y1+((sqrt(3._ark)*y8**3*y9+sqrt(3._ark)*y7**3*y9)*y5+(-y8**3*y9+ &
      y7**3*y9)*y6)*y2+((sqrt(3._ark)*y8**3*y9-sqrt(3._ark)*y7**3*y9)*y5+(-y7**3*y9- &
      y8**3*y9)*y6)*y3+((sqrt(3._ark)*y7**3*y9-sqrt(3._ark)*y8**3*y9)*y5+(y7**3*y9+ &
      y8**3*y9)*y6)*y4
    dF(3,502) = (y6**2*y9**3+y5**2*y9**3)*y1+(y6**2*y9**3+y5**2*y9**3)*y2+ &
      (y6**2*y9**3+y5**2*y9**3)*y3+(y6**2*y9**3+y5**2*y9**3)*y4
    dF(3,503) = ((y7*y8**2+y7**2*y8)*y5**2+(y7*y8**2+y7**2*y8)*y6**2)*y1+((- &
      y7**2*y8-y7*y8**2)*y5**2+(-y7**2*y8-y7*y8**2)*y6**2)*y2+((y7**2*y8- &
      y7*y8**2)*y5**2+(y7**2*y8-y7*y8**2)*y6**2)*y3+((y7*y8**2-y7**2*y8)*y5**2+ &
      (y7*y8**2-y7**2*y8)*y6**2)*y4
    dF(3,504) = (y5**2*y7*y8*y9+y6**2*y7*y8*y9)*y1+(y5**2*y7*y8*y9+ &
      y6**2*y7*y8*y9)*y2+(-y5**2*y7*y8*y9-y6**2*y7*y8*y9)*y3+(-y5**2*y7*y8*y9- &
      y6**2*y7*y8*y9)*y4
    dF(3,505) = (-sqrt(3._ark)*y5*y6**2*y7*y8-5._ark/9._ark*sqrt(3._ark)*y5**3*y7*y8)*y1+(- &
      sqrt(3._ark)*y5*y6**2*y7*y8-5._ark/9._ark*sqrt(3._ark)*y5**3*y7*y8)*y2+(- &
      sqrt(3._ark)*y5*y6**2*y7*y8-5._ark/9._ark*sqrt(3._ark)*y5**3*y7*y8)*y3+(- &
      sqrt(3._ark)*y5*y6**2*y7*y8-5._ark/9._ark*sqrt(3._ark)*y5**3*y7*y8)*y4
    dF(3,506) = ((3._ark/4._ark*sqrt(3._ark)*y8**2+3._ark/4._ark*sqrt(3._ark)*y7**2)*y5**3+ &
      (3._ark/4._ark*sqrt(3._ark)*y8**2+3._ark/4._ark*sqrt(3._ark)*y7**2)*y6**2*y5+(-y8**2+ &
      y7**2)*y6**3)*y1+((3._ark/4._ark*sqrt(3._ark)*y8**2+3._ark/4._ark*sqrt(3._ark)*y7**2)*y5**3+ &
      (3._ark/4._ark*sqrt(3._ark)*y8**2+3._ark/4._ark*sqrt(3._ark)*y7**2)*y6**2*y5+(-y8**2+ &
      y7**2)*y6**3)*y2+((-3._ark/4._ark*sqrt(3._ark)*y8**2-3._ark/ &
      4._ark*sqrt(3._ark)*y7**2)*y5**3+(-3._ark/4._ark*sqrt(3._ark)*y8**2-3._ark/ &
      4._ark*sqrt(3._ark)*y7**2)*y6**2*y5+(-y7**2+y8**2)*y6**3)*y3+((-3._ark/ &
      4._ark*sqrt(3._ark)*y8**2-3._ark/4._ark*sqrt(3._ark)*y7**2)*y5**3+(-3._ark/ &
      4._ark*sqrt(3._ark)*y8**2-3._ark/4._ark*sqrt(3._ark)*y7**2)*y6**2*y5+(-y7**2+ &
      y8**2)*y6**3)*y4
    dF(3,507) = ((y7**2*y9**2+y8**2*y9**2)*y5+(sqrt(3._ark)*y8**2*y9**2- &
      sqrt(3._ark)*y7**2*y9**2)*y6)*y1+((y7**2*y9**2+y8**2*y9**2)*y5+ &
      (sqrt(3._ark)*y8**2*y9**2-sqrt(3._ark)*y7**2*y9**2)*y6)*y2+((-y7**2*y9**2- &
      y8**2*y9**2)*y5+(-sqrt(3._ark)*y8**2*y9**2+sqrt(3._ark)*y7**2*y9**2)*y6)*y3+((- &
      y7**2*y9**2-y8**2*y9**2)*y5+(-sqrt(3._ark)*y8**2*y9**2+ &
      sqrt(3._ark)*y7**2*y9**2)*y6)*y4
    dF(3,508) = ((y7**3*y9+y8**3*y9)*y5+(sqrt(3._ark)*y7**3*y9- &
      sqrt(3._ark)*y8**3*y9)*y6)*y1+((-y7**3*y9-y8**3*y9)*y5+(sqrt(3._ark)*y8**3*y9- &
      sqrt(3._ark)*y7**3*y9)*y6)*y2+((-y8**3*y9+y7**3*y9)*y5+(sqrt(3._ark)*y8**3*y9+ &
      sqrt(3._ark)*y7**3*y9)*y6)*y3+((-y7**3*y9+y8**3*y9)*y5+(-sqrt(3._ark)*y7**3*y9- &
      sqrt(3._ark)*y8**3*y9)*y6)*y4
    dF(3,509) = 2._ark*y4*y5*y7**2*y8**2-2._ark*y1*y5*y7**2*y8**2+ &
      2._ark*y3*y5*y7**2*y8**2-2._ark*y2*y5*y7**2*y8**2
    dF(3,510) = ((y7*y8**2*y9+y7**2*y8*y9)*y5+(sqrt(3._ark)*y7*y8**2*y9- &
      sqrt(3._ark)*y7**2*y8*y9)*y6)*y1+((-y7**2*y8*y9-y7*y8**2*y9)*y5+ &
      (sqrt(3._ark)*y7**2*y8*y9-sqrt(3._ark)*y7*y8**2*y9)*y6)*y2+((y7*y8**2*y9- &
      y7**2*y8*y9)*y5+(sqrt(3._ark)*y7**2*y8*y9+sqrt(3._ark)*y7*y8**2*y9)*y6)*y3+((- &
      y7*y8**2*y9+y7**2*y8*y9)*y5+(-sqrt(3._ark)*y7**2*y8*y9- &
      sqrt(3._ark)*y7*y8**2*y9)*y6)*y4
    dF(3,511) = ((-3._ark/2._ark*y8*y9-3._ark/2._ark*y7*y9)*y5**3+(sqrt(3._ark)*y8*y9- &
      sqrt(3._ark)*y7*y9)*y6*y5**2+(-y8*y9/2._ark-y7*y9/2._ark)*y6**2*y5)*y1+((3._ark/ &
      2._ark*y8*y9+3._ark/2._ark*y7*y9)*y5**3+(sqrt(3._ark)*y7*y9-sqrt(3._ark)*y8*y9)*y6*y5**2+ &
      (y8*y9/2._ark+y7*y9/2._ark)*y6**2*y5)*y2+((-3._ark/2._ark*y7*y9+3._ark/2._ark*y8*y9)*y5**3+ &
      (-sqrt(3._ark)*y8*y9-sqrt(3._ark)*y7*y9)*y6*y5**2+(-y7*y9/2._ark+y8*y9/ &
      2._ark)*y6**2*y5)*y3+((-3._ark/2._ark*y8*y9+3._ark/2._ark*y7*y9)*y5**3+(sqrt(3._ark)*y7*y9+ &
      sqrt(3._ark)*y8*y9)*y6*y5**2+(-y8*y9/2._ark+y7*y9/2._ark)*y6**2*y5)*y4
    dF(3,512) = (6._ark*sqrt(3._ark)*y5**3*y6**2+10._ark/3._ark*sqrt(3._ark)*y5*y6**4)*y1+ &
      (6._ark*sqrt(3._ark)*y5**3*y6**2+10._ark/3._ark*sqrt(3._ark)*y5*y6**4)*y2+(- &
      6._ark*sqrt(3._ark)*y5**3*y6**2-10._ark/3._ark*sqrt(3._ark)*y5*y6**4)*y3+(- &
      6._ark*sqrt(3._ark)*y5**3*y6**2-10._ark/3._ark*sqrt(3._ark)*y5*y6**4)*y4
    dF(3,513) = ((-7._ark/2._ark*y8*y9-7._ark/2._ark*y7*y9)*y5**3+(3._ark*sqrt(3._ark)*y8*y9- &
      3._ark*sqrt(3._ark)*y7*y9)*y6*y5**2+(-9._ark/2._ark*y8*y9-9._ark/2._ark*y7*y9)*y6**2*y5)*y1+ &
      ((7._ark/2._ark*y8*y9+7._ark/2._ark*y7*y9)*y5**3+(-3._ark*sqrt(3._ark)*y8*y9+ &
      3._ark*sqrt(3._ark)*y7*y9)*y6*y5**2+(9._ark/2._ark*y8*y9+9._ark/2._ark*y7*y9)*y6**2*y5)*y2+ &
      ((-7._ark/2._ark*y7*y9+7._ark/2._ark*y8*y9)*y5**3+(-3._ark*sqrt(3._ark)*y7*y9- &
      3._ark*sqrt(3._ark)*y8*y9)*y6*y5**2+(-9._ark/2._ark*y7*y9+9._ark/2._ark*y8*y9)*y6**2*y5)*y3+ &
      ((-7._ark/2._ark*y8*y9+7._ark/2._ark*y7*y9)*y5**3+(3._ark*sqrt(3._ark)*y7*y9+ &
      3._ark*sqrt(3._ark)*y8*y9)*y6*y5**2+(-9._ark/2._ark*y8*y9+9._ark/2._ark*y7*y9)*y6**2*y5)*y4
    dF(3,514) = (y3*y7*y8*y9**2+y4*y7*y8*y9**2)*y1+(y3*y7*y8*y9**2+ &
      y4*y7*y8*y9**2)*y2
    dF(3,515) = ((-y6**3*y9-3._ark/5._ark*sqrt(3._ark)*y5*y6**2*y9)*y3+(y6**3*y9-3._ark/ &
      5._ark*sqrt(3._ark)*y5*y6**2*y9)*y4)*y1+((y6**3*y9-3._ark/ &
      5._ark*sqrt(3._ark)*y5*y6**2*y9)*y3+(-y6**3*y9-3._ark/ &
      5._ark*sqrt(3._ark)*y5*y6**2*y9)*y4)*y2
    dF(3,516) = (8._ark/5._ark*y4*y5*y6**2*y9+8._ark/5._ark*y3*y5*y6**2*y9)*y1+(8._ark/ &
      5._ark*y4*y5*y6**2*y9+8._ark/5._ark*y3*y5*y6**2*y9)*y2
    dF(3,517) = ((-sqrt(3._ark)*y6**3*y8/4._ark-sqrt(3._ark)*y5**2*y6*y8/4._ark+ &
      y5*y6**2*y8)*y3+(sqrt(3._ark)*y6**3*y7/4._ark+sqrt(3._ark)*y5**2*y6*y7/4._ark+ &
      y5*y6**2*y7)*y4)*y1+((-sqrt(3._ark)*y5**2*y6*y7/4._ark-y5*y6**2*y7- &
      sqrt(3._ark)*y6**3*y7/4._ark)*y3+(-y5*y6**2*y8+sqrt(3._ark)*y6**3*y8/4._ark+ &
      sqrt(3._ark)*y5**2*y6*y8/4._ark)*y4)*y2
    dF(3,518) = ((y5**2*y7*y8+y6**2*y7*y8)*y3+(y5**2*y7*y8+y6**2*y7*y8)*y4)*y1+ &
      ((y5**2*y7*y8+y6**2*y7*y8)*y3+(y5**2*y7*y8+y6**2*y7*y8)*y4)*y2
    dF(3,519) = y1**5*y2+y1*y2**5-y3*y4**5-y3**5*y4
    dF(3,520) = y1*y2*y7**2*y8**2-y3*y4*y7**2*y8**2
    dF(3,521) = (y4*y7**3*y8+y3*y7*y8**3)*y1+(y4*y7*y8**3+y3*y7**3*y8)*y2
    dF(3,522) = (y4*y8*y9**3+y3*y7*y9**3)*y1+(-y4*y7*y9**3-y3*y8*y9**3)*y2
    dF(3,523) = (y7**2*y9**2+y8**2*y9**2)*y2*y1+(-y7**2*y9**2-y8**2*y9**2)*y4*y3
    dF(3,524) = 8._ark/9._ark*sqrt(3._ark)*y1*y2*y5**3*y9+8._ark/ &
      9._ark*sqrt(3._ark)*y3*y4*y5**3*y9
    dF(3,525) = -2._ark/3._ark*sqrt(3._ark)*y3*y4*y5**2*y7*y8-2._ark/ &
      3._ark*sqrt(3._ark)*y1*y2*y5**2*y7*y8
    dF(3,526) = (2._ark/3._ark*sqrt(3._ark)*y8**2+2._ark/3._ark*sqrt(3._ark)*y7**2)*y6**2*y2*y1+ &
      (-2._ark/3._ark*sqrt(3._ark)*y8**2-2._ark/3._ark*sqrt(3._ark)*y7**2)*y6**2*y4*y3
    dF(3,527) = ((y5*y6*y7*y9-sqrt(3._ark)*y6**2*y7*y9/2._ark-sqrt(3._ark)*y5**2*y7*y9/ &
      6._ark)*y3+(-sqrt(3._ark)*y5**2*y8*y9/6._ark-y5*y6*y8*y9-sqrt(3._ark)*y6**2*y8*y9/ &
      2._ark)*y4)*y1+((sqrt(3._ark)*y6**2*y8*y9/2._ark+sqrt(3._ark)*y5**2*y8*y9/6._ark+ &
      y5*y6*y8*y9)*y3+(sqrt(3._ark)*y5**2*y7*y9/6._ark-y5*y6*y7*y9+sqrt(3._ark)*y6**2*y7*y9/ &
      2._ark)*y4)*y2
    dF(3,528) = ((-sqrt(3._ark)*y6**2*y7*y8/3._ark+y5*y6*y7*y8)*y3+(-y5*y6*y7*y8- &
      sqrt(3._ark)*y6**2*y7*y8/3._ark)*y4)*y1+((-y5*y6*y7*y8-sqrt(3._ark)*y6**2*y7*y8/ &
      3._ark)*y3+(-sqrt(3._ark)*y6**2*y7*y8/3._ark+y5*y6*y7*y8)*y4)*y2
    dF(3,529) = ((y5**2*y6*y9+sqrt(3._ark)*y5*y6**2*y9/5._ark)*y3+ &
      (sqrt(3._ark)*y5*y6**2*y9/5._ark-y5**2*y6*y9)*y4)*y1+((sqrt(3._ark)*y5*y6**2*y9/5._ark- &
      y5**2*y6*y9)*y3+(y5**2*y6*y9+sqrt(3._ark)*y5*y6**2*y9/5._ark)*y4)*y2
    dF(3,530) = -4._ark/3._ark*sqrt(3._ark)*y1*y2*y5**4+4._ark/3._ark*sqrt(3._ark)*y3*y4*y5**4
    dF(3,531) = ((-sqrt(3._ark)*y8*y9-sqrt(3._ark)*y7*y9)*y5+(-y7*y9+ &
      y8*y9)*y6)*y2*y1**2+((sqrt(3._ark)*y7*y9+sqrt(3._ark)*y8*y9)*y5+(-y8*y9+ &
      y7*y9)*y6)*y2**2*y1+((sqrt(3._ark)*y8*y9-sqrt(3._ark)*y7*y9)*y5+(-y7*y9- &
      y8*y9)*y6)*y4*y3**2+((sqrt(3._ark)*y7*y9-sqrt(3._ark)*y8*y9)*y5+(y7*y9+ &
      y8*y9)*y6)*y4**2*y3
    dF(3,532) = (-2._ark*y4*y6*y8*y9+2._ark*y3*y6*y7*y9)*y1**2+(2._ark*y3**2*y6*y7*y9- &
      2._ark*y4**2*y6*y8*y9)*y1+(2._ark*y3*y6*y8*y9-2._ark*y4*y6*y7*y9)*y2**2+(- &
      2._ark*y4**2*y6*y7*y9+2._ark*y3**2*y6*y8*y9)*y2
    dF(3,533) = ((-y6**3+7._ark/9._ark*y5**2*y6-4._ark/9._ark*sqrt(3._ark)*y5*y6**2)*y3+(- &
      4._ark/9._ark*sqrt(3._ark)*y5*y6**2-7._ark/9._ark*y5**2*y6+y6**3)*y4)*y1**2+((-7._ark/ &
      9._ark*y5**2*y6+y6**3+4._ark/9._ark*sqrt(3._ark)*y5*y6**2)*y3**2+(7._ark/9._ark*y5**2*y6- &
      y6**3+4._ark/9._ark*sqrt(3._ark)*y5*y6**2)*y4**2)*y1+((-4._ark/9._ark*sqrt(3._ark)*y5*y6**2- &
      7._ark/9._ark*y5**2*y6+y6**3)*y3+(-y6**3+7._ark/9._ark*y5**2*y6-4._ark/ &
      9._ark*sqrt(3._ark)*y5*y6**2)*y4)*y2**2+((7._ark/9._ark*y5**2*y6-y6**3+4._ark/ &
      9._ark*sqrt(3._ark)*y5*y6**2)*y3**2+(-7._ark/9._ark*y5**2*y6+y6**3+4._ark/ &
      9._ark*sqrt(3._ark)*y5*y6**2)*y4**2)*y2
    dF(3,534) = ((-y7*y9-y8*y9)*y5+(sqrt(3._ark)*y8*y9-sqrt(3._ark)*y7*y9)*y6)*y2*y1**2+ &
      ((y7*y9+y8*y9)*y5+(sqrt(3._ark)*y7*y9-sqrt(3._ark)*y8*y9)*y6)*y2**2*y1+((-y7*y9+ &
      y8*y9)*y5+(-sqrt(3._ark)*y8*y9-sqrt(3._ark)*y7*y9)*y6)*y4*y3**2+((-y8*y9+y7*y9)*y5+ &
      (sqrt(3._ark)*y7*y9+sqrt(3._ark)*y8*y9)*y6)*y4**2*y3
    dF(3,535) = 2._ark/3._ark*sqrt(3._ark)*y1*y2**2*y6**2*y9+2._ark/ &
      3._ark*sqrt(3._ark)*y3*y4**2*y6**2*y9+2._ark/3._ark*sqrt(3._ark)*y3**2*y4*y6**2*y9+2._ark/ &
      3._ark*sqrt(3._ark)*y1**2*y2*y6**2*y9
    dF(3,536) = ((-y5*y6*y7/2._ark-sqrt(3._ark)*y5**2*y7/2._ark)*y3+(-sqrt(3._ark)*y5**2*y8/ &
      2._ark+y5*y6*y8/2._ark)*y4)*y1**2+((y5*y6*y7/2._ark+sqrt(3._ark)*y5**2*y7/2._ark)*y3**2+(- &
      y5*y6*y8/2._ark+sqrt(3._ark)*y5**2*y8/2._ark)*y4**2)*y1+((-y5*y6*y8/2._ark+ &
      sqrt(3._ark)*y5**2*y8/2._ark)*y3+(y5*y6*y7/2._ark+sqrt(3._ark)*y5**2*y7/2._ark)*y4)*y2**2+ &
      ((-sqrt(3._ark)*y5**2*y8/2._ark+y5*y6*y8/2._ark)*y3**2+(-y5*y6*y7/2._ark- &
      sqrt(3._ark)*y5**2*y7/2._ark)*y4**2)*y2
    dF(3,537) = (y3*y7*y8+y4*y7*y8)*y2*y1**2+((y3*y7*y8+y4*y7*y8)*y2**2+ &
      y3*y4**2*y7*y8+y3**2*y4*y7*y8)*y1+(y3*y4**2*y7*y8+y3**2*y4*y7*y8)*y2
    dF(3,538) = (-y7*y9-y8*y9)*y4*y3*y1**2+((-y7*y9+y8*y9)*y3**2+(-y8*y9+ &
      y7*y9)*y4**2)*y2*y1+(y7*y9+y8*y9)*y4*y3*y2**2
    dF(3,539) = -y1**2*y3*y4*y9**2+(y4**2*y9**2+y3**2*y9**2)*y2*y1-y2**2*y3*y4*y9**2
    dF(3,540) = (-y3*y8**2-y4*y7**2)*y2*y1**2+((-y3*y7**2-y4*y8**2)*y2**2+ &
      y3*y4**2*y7**2+y3**2*y4*y8**2)*y1+(y3*y4**2*y8**2+y3**2*y4*y7**2)*y2
    dF(3,541) = (y3*y7*y8**2*y9+y4*y7**2*y8*y9)*y1+(-y3*y7**2*y8*y9- &
      y4*y7*y8**2*y9)*y2
    dF(3,542) = (y7*y8**3+y7**3*y8)*y2*y1+(y7*y8**3+y7**3*y8)*y4*y3
    dF(3,543) = (y7**4+y8**4)*y2*y1+(-y8**4-y7**4)*y4*y3
    dF(3,544) = ((sqrt(3._ark)*y8**2*y9+sqrt(3._ark)*y7**2*y9)*y5+(y7**2*y9- &
      y8**2*y9)*y6)*y2*y1+((sqrt(3._ark)*y8**2*y9+sqrt(3._ark)*y7**2*y9)*y5+(y7**2*y9- &
      y8**2*y9)*y6)*y4*y3
    dF(3,545) = ((sqrt(3._ark)*y5*y8*y9**2/2._ark+y6*y8*y9**2/2._ark)*y3+ &
      (sqrt(3._ark)*y5*y7*y9**2/2._ark-y6*y7*y9**2/2._ark)*y4)*y1+((y6*y7*y9**2/2._ark- &
      sqrt(3._ark)*y5*y7*y9**2/2._ark)*y3+(-y6*y8*y9**2/2._ark-sqrt(3._ark)*y5*y8*y9**2/ &
      2._ark)*y4)*y2
    dF(3,546) = ((-5._ark/4._ark*y6**3*y8-9._ark/4._ark*y5**2*y6*y8)*y3+(5._ark/4._ark*y6**3*y7+ &
      9._ark/4._ark*y5**2*y6*y7)*y4)*y1+((-9._ark/4._ark*y5**2*y6*y7-5._ark/4._ark*y6**3*y7)*y3+ &
      (5._ark/4._ark*y6**3*y8+9._ark/4._ark*y5**2*y6*y8)*y4)*y2
    dF(3,547) = (3._ark*y5**4+y6**4)*y2*y1+(-3._ark*y5**4-y6**4)*y4*y3
    dF(3,548) = y3*y4*y5*y9**3+y1*y2*y5*y9**3
    dF(3,549) = ((sqrt(3._ark)*y6*y7**2*y8/2._ark-y5*y7**2*y8/2._ark)*y3+(-y5*y7*y8**2/ &
      2._ark-sqrt(3._ark)*y6*y7*y8**2/2._ark)*y4)*y1+((sqrt(3._ark)*y6*y7*y8**2/2._ark+ &
      y5*y7*y8**2/2._ark)*y3+(y5*y7**2*y8/2._ark-sqrt(3._ark)*y6*y7**2*y8/2._ark)*y4)*y2
    dF(3,550) = -y3*y4*y5*y7*y8*y9+y1*y2*y5*y7*y8*y9
    dF(3,551) = (-2._ark/3._ark*sqrt(3._ark)*y3*y6**2*y7*y8-2._ark/ &
      3._ark*sqrt(3._ark)*y4*y6**2*y7*y8)*y1+(-2._ark/3._ark*sqrt(3._ark)*y3*y6**2*y7*y8-2._ark/ &
      3._ark*sqrt(3._ark)*y4*y6**2*y7*y8)*y2
    dF(3,552) = (-y5**3*y9/3._ark+y9*y6**2*y5)*y2*y1+(-y5**3*y9/3._ark+ &
      y9*y6**2*y5)*y4*y3
    dF(3,553) = (y6**2*y9**2+y5**2*y9**2)*y2*y1+(-y6**2*y9**2-y5**2*y9**2)*y4*y3
    dF(3,554) = ((y8**2+y7**2)*y5**2+(y8**2+y7**2)*y6**2)*y2*y1+((-y8**2- &
      y7**2)*y5**2+(-y8**2-y7**2)*y6**2)*y4*y3
    dF(3,555) = ((y5**2*y6*y8/4._ark-3._ark/4._ark*y6**3*y8)*y3+(-y5**2*y6*y7/4._ark+3._ark/ &
      4._ark*y6**3*y7)*y4)*y1+((-3._ark/4._ark*y6**3*y7+y5**2*y6*y7/4._ark)*y3+(-y5**2*y6*y8/ &
      4._ark+3._ark/4._ark*y6**3*y8)*y4)*y2
    dF(3,556) = ((y3*y9**3+y4*y9**3)*y2+y3*y4*y9**3)*y1+y2*y3*y4*y9**3
    dF(3,557) = ((((-sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5+(-y8*y9/2._ark- &
      y7*y9/2._ark)*y6)*y3+((sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y5+(y8*y9/ &
      2._ark+y7*y9/2._ark)*y6)*y4)*y2+((sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y5+(- &
      y7*y9/2._ark+y8*y9/2._ark)*y6)*y4*y3)*y1+((-sqrt(3._ark)*y7*y9/2._ark-sqrt(3._ark)*y8*y9/ &
      2._ark)*y5+(-y8*y9/2._ark+y7*y9/2._ark)*y6)*y4*y3*y2
    dF(3,558) = ((((y8**2+y7**2)*y5+(-sqrt(3._ark)*y8**2+sqrt(3._ark)*y7**2)*y6)*y3+ &
      ((y8**2+y7**2)*y5+(-sqrt(3._ark)*y8**2+sqrt(3._ark)*y7**2)*y6)*y4)*y2+((-y8**2- &
      y7**2)*y5+(sqrt(3._ark)*y8**2-sqrt(3._ark)*y7**2)*y6)*y4*y3)*y1+((-y8**2-y7**2)*y5+ &
      (sqrt(3._ark)*y8**2-sqrt(3._ark)*y7**2)*y6)*y4*y3*y2
    dF(3,559) = ((((-y7*y9/2._ark+y8*y9/2._ark)*y5+(-sqrt(3._ark)*y7*y9/2._ark- &
      sqrt(3._ark)*y8*y9/2._ark)*y6)*y3+((-y8*y9/2._ark+y7*y9/2._ark)*y5+(sqrt(3._ark)*y8*y9/ &
      2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y4)*y2+((-y8*y9/2._ark-y7*y9/2._ark)*y5+ &
      (sqrt(3._ark)*y8*y9/2._ark-sqrt(3._ark)*y7*y9/2._ark)*y6)*y4*y3)*y1+((y8*y9/2._ark+y7*y9/ &
      2._ark)*y5+(-sqrt(3._ark)*y8*y9/2._ark+sqrt(3._ark)*y7*y9/2._ark)*y6)*y4*y3*y2
    dF(3,560) = ((((sqrt(3._ark)*y7-sqrt(3._ark)*y8)*y5**2+(-y7-y8)*y6*y5)*y3+ &
      ((sqrt(3._ark)*y8-sqrt(3._ark)*y7)*y5**2+(y8+y7)*y6*y5)*y4)*y2+((-sqrt(3._ark)*y7- &
      sqrt(3._ark)*y8)*y5**2+(-y8+y7)*y6*y5)*y4*y3)*y1+((sqrt(3._ark)*y8+ &
      sqrt(3._ark)*y7)*y5**2+(y8-y7)*y6*y5)*y4*y3*y2
    dF(3,561) = (((y5*y6**2-y5**3/3._ark)*y3+(y5*y6**2-y5**3/3._ark)*y4)*y2+(y5**3/3._ark- &
      y5*y6**2)*y4*y3)*y1+(y5**3/3._ark-y5*y6**2)*y4*y3*y2
    dF(3,562) = (-y4*y7-y3*y8)*y2*y1**3+((y4*y8+y3*y7)*y2**3-y3*y4**3*y7- &
      y3**3*y4*y8)*y1+(y3**3*y4*y7+y3*y4**3*y8)*y2
    dF(3,563) = -y1**3*y3*y4*y5+(y3**3*y5+y4**3*y5)*y2*y1-y2**3*y3*y4*y5
    dF(3,564) = (((y8**3-y7**3)*y3+(y7**3-y8**3)*y4)*y2+(y7**3+y8**3)*y4*y3)*y1+(- &
      y7**3-y8**3)*y4*y3*y2
    dF(3,565) = ((y4*y7*y8*y9+y3*y7*y8*y9)*y2-y3*y4*y7*y8*y9)*y1-y2*y3*y4*y7*y8*y9
    dF(3,566) = (((y6**2*y9+y5**2*y9)*y3+(y6**2*y9+y5**2*y9)*y4)*y2+(y6**2*y9+ &
      y5**2*y9)*y4*y3)*y1+(y6**2*y9+y5**2*y9)*y4*y3*y2
    dF(3,567) = ((((3._ark*y8-3._ark*y7)*y5**2+(-y8+y7)*y6**2)*y3+((-3._ark*y8+ &
      3._ark*y7)*y5**2+(y8-y7)*y6**2)*y4)*y2+((3._ark*y8+3._ark*y7)*y5**2+(-y7- &
      y8)*y6**2)*y4*y3)*y1+((-3._ark*y7-3._ark*y8)*y5**2+(y8+y7)*y6**2)*y4*y3*y2
    dF(3,568) = ((8._ark/9._ark*sqrt(3._ark)*y3*y5**3+8._ark/9._ark*sqrt(3._ark)*y4*y5**3)*y2- &
      8._ark/9._ark*sqrt(3._ark)*y3*y4*y5**3)*y1-8._ark/9._ark*sqrt(3._ark)*y2*y3*y4*y5**3
    dF(3,569) = ((y4*y5*y7*y8+y3*y5*y7*y8)*y2+y3*y4*y5*y7*y8)*y1+y2*y3*y4*y5*y7*y8
    dF(3,570) = ((-2._ark*y3*y5*y9**2-2._ark*y4*y5*y9**2)*y2+2._ark*y3*y4*y5*y9**2)*y1+ &
      2._ark*y2*y3*y4*y5*y9**2
    dF(3,571) = ((-2._ark/3._ark*sqrt(3._ark)*y3*y5**2*y9-2._ark/ &
      3._ark*sqrt(3._ark)*y4*y5**2*y9)*y2-2._ark/3._ark*sqrt(3._ark)*y3*y4*y5**2*y9)*y1-2._ark/ &
      3._ark*sqrt(3._ark)*y2*y3*y4*y5**2*y9
    dF(3,572) = (((4._ark*y8-4._ark*y7)*y5**2*y3+(-4._ark*y8+4._ark*y7)*y5**2*y4)*y2+ &
      (4._ark*y8+4._ark*y7)*y5**2*y4*y3)*y1+(-4._ark*y8-4._ark*y7)*y5**2*y4*y3*y2
    dF(3,573) = y1*y2*y3*y4*y7*y8
    dF(3,574) = 2._ark/3._ark*sqrt(3._ark)*y1*y2*y3*y4*y5*y9
    dF(3,575) = (y4*y9+y3*y9)*y2*y1**3+((y4*y9+y3*y9)*y2**3+y3**3*y4*y9+ &
      y3*y4**3*y9)*y1+(y3**3*y4*y9+y3*y4**3*y9)*y2
    dF(3,576) = (-y4-y3)*y2*y1**4+((-y4-y3)*y2**4+y3*y4**4+y3**4*y4)*y1+(y3*y4**4+ &
      y3**4*y4)*y2
    dF(3,577) = -y1**4*y3*y4+(y3**4+y4**4)*y2*y1-y2**4*y3*y4
    dF(3,578) = (y3*y7**2*y9+y4*y8**2*y9)*y1**2+(y4**2*y8**2*y9+y3**2*y7**2*y9)*y1+ &
      (y4*y7**2*y9+y3*y8**2*y9)*y2**2+(y3**2*y8**2*y9+y4**2*y7**2*y9)*y2
    dF(3,579) = (2._ark/3._ark*sqrt(3._ark)*y4*y5**2+2._ark/ &
      3._ark*sqrt(3._ark)*y3*y5**2)*y2*y1**2+((2._ark/3._ark*sqrt(3._ark)*y4*y5**2+2._ark/ &
      3._ark*sqrt(3._ark)*y3*y5**2)*y2**2-2._ark/3._ark*sqrt(3._ark)*y3**2*y4*y5**2-2._ark/ &
      3._ark*sqrt(3._ark)*y3*y4**2*y5**2)*y1+(-2._ark/3._ark*sqrt(3._ark)*y3**2*y4*y5**2-2._ark/ &
      3._ark*sqrt(3._ark)*y3*y4**2*y5**2)*y2
    dF(3,580) = (sqrt(3._ark)*y6**2/6._ark-sqrt(3._ark)*y5**2/2._ark)*y4*y3*y1**2+ &
      ((sqrt(3._ark)*y5**2/2._ark-sqrt(3._ark)*y6**2/6._ark)*y3**2+(sqrt(3._ark)*y5**2/2._ark- &
      sqrt(3._ark)*y6**2/6._ark)*y4**2)*y2*y1+(sqrt(3._ark)*y6**2/6._ark-sqrt(3._ark)*y5**2/ &
      2._ark)*y4*y3*y2**2
    dF(3,581) = (y4*y9+y3*y9)*y2**2*y1**2+y1*y3**2*y4**2*y9+y2*y3**2*y4**2*y9
    dF(3,582) = ((-y4**2*y8-y3**2*y7)*y2+y3*y4**2*y8+y3**2*y4*y7)*y1**2+(y3**2*y8+ &
      y4**2*y7)*y2**2*y1+(-y3*y4**2*y7-y3**2*y4*y8)*y2**2
    dF(3,583) = (((sqrt(3._ark)*y6+y5)*y3**2+(y5-sqrt(3._ark)*y6)*y4**2)*y2+(-y5- &
      sqrt(3._ark)*y6)*y4*y3**2+(sqrt(3._ark)*y6-y5)*y4**2*y3)*y1**2+((y5- &
      sqrt(3._ark)*y6)*y3**2+(sqrt(3._ark)*y6+y5)*y4**2)*y2**2*y1+((sqrt(3._ark)*y6- &
      y5)*y4*y3**2+(-y5-sqrt(3._ark)*y6)*y4**2*y3)*y2**2
    dF(3,584) = (y4**2+y3**2)*y2*y1**3+(-y3**3*y4-y3*y4**3)*y1**2+(y4**2+ &
      y3**2)*y2**3*y1+(-y3**3*y4-y3*y4**3)*y2**2
    dF(3,585) = y1**2*y3*y4*y7*y8+(y3**2*y7*y8+y4**2*y7*y8)*y2*y1+y2**2*y3*y4*y7*y8
    dF(3,586) = ((-sqrt(3._ark)*y7-sqrt(3._ark)*y8)*y5+(-y8+y7)*y6)*y4*y3*y1**2+ &
      (((sqrt(3._ark)*y7-sqrt(3._ark)*y8)*y5+(-y7-y8)*y6)*y3**2+((sqrt(3._ark)*y8- &
      sqrt(3._ark)*y7)*y5+(y8+y7)*y6)*y4**2)*y2*y1+((sqrt(3._ark)*y8+sqrt(3._ark)*y7)*y5+ &
      (y8-y7)*y6)*y4*y3*y2**2
    dF(3,587) = (-y7-y8)*y4*y3*y2*y1**2+((y8+y7)*y4*y3*y2**2+((-y8+y7)*y4*y3**2+(y8- &
      y7)*y4**2*y3)*y2)*y1
    dF(3,588) = y1**2*y2*y3*y4*y9+(y2**2*y3*y4*y9+(y3*y4**2*y9+y3**2*y4*y9)*y2)*y1
    dF(3,589) = 2._ark*y1**2*y2*y3*y4*y5+(2._ark*y2**2*y3*y4*y5+(-2._ark*y3*y4**2*y5- &
      2._ark*y3**2*y4*y5)*y2)*y1
    dF(3,590) = y1**2*y9**4+y2**2*y9**4-y4**2*y9**4-y3**2*y9**4
    dF(3,591) = (y8*y9**3+y7*y9**3)*y1**2+(-y7*y9**3-y8*y9**3)*y2**2+(y7*y9**3- &
      y8*y9**3)*y3**2+(y8*y9**3-y7*y9**3)*y4**2
    dF(3,592) = (y7**4+y8**4)*y1**2+(y7**4+y8**4)*y2**2+(-y8**4-y7**4)*y3**2+(- &
      y8**4-y7**4)*y4**2
    dF(3,593) = (-2._ark*y8**3+2._ark*y7**3)*y6*y1**2+(2._ark*y8**3-2._ark*y7**3)*y6*y2**2+ &
      (-2._ark*y7**3-2._ark*y8**3)*y6*y3**2+(2._ark*y7**3+2._ark*y8**3)*y6*y4**2
    dF(3,594) = ((3._ark*sqrt(3._ark)*y7+3._ark*sqrt(3._ark)*y8)*y5**3+(-y8+ &
      y7)*y6**3)*y1**2+((-3._ark*sqrt(3._ark)*y8-3._ark*sqrt(3._ark)*y7)*y5**3+(y8- &
      y7)*y6**3)*y2**2+((-3._ark*sqrt(3._ark)*y7+3._ark*sqrt(3._ark)*y8)*y5**3+(-y7- &
      y8)*y6**3)*y3**2+((-3._ark*sqrt(3._ark)*y8+3._ark*sqrt(3._ark)*y7)*y5**3+(y8+ &
      y7)*y6**3)*y4**2
    dF(3,595) = -y4**2*y5*y7*y8*y9-y3**2*y5*y7*y8*y9+y2**2*y5*y7*y8*y9+ &
      y1**2*y5*y7*y8*y9
    dF(3,596) = -2._ark*y3**2*y5*y9**3-2._ark*y2**2*y5*y9**3-2._ark*y1**2*y5*y9**3- &
      2._ark*y4**2*y5*y9**3
    dF(3,597) = ((y7**3+y8**3)*y5+(sqrt(3._ark)*y8**3-sqrt(3._ark)*y7**3)*y6)*y1**2+((- &
      y7**3-y8**3)*y5+(-sqrt(3._ark)*y8**3+sqrt(3._ark)*y7**3)*y6)*y2**2+((y8**3- &
      y7**3)*y5+(sqrt(3._ark)*y7**3+sqrt(3._ark)*y8**3)*y6)*y3**2+((y7**3-y8**3)*y5+(- &
      sqrt(3._ark)*y7**3-sqrt(3._ark)*y8**3)*y6)*y4**2
    dF(3,598) = -2._ark/3._ark*sqrt(3._ark)*y1**2*y6**2*y9**2+2._ark/ &
      3._ark*sqrt(3._ark)*y4**2*y6**2*y9**2+2._ark/3._ark*sqrt(3._ark)*y3**2*y6**2*y9**2-2._ark/ &
      3._ark*sqrt(3._ark)*y2**2*y6**2*y9**2
    dF(3,599) = (2._ark/3._ark*sqrt(3._ark)*y8*y9+2._ark/3._ark*sqrt(3._ark)*y7*y9)*y6**2*y1**2+ &
      (-2._ark/3._ark*sqrt(3._ark)*y7*y9-2._ark/3._ark*sqrt(3._ark)*y8*y9)*y6**2*y2**2+(2._ark/ &
      3._ark*sqrt(3._ark)*y7*y9-2._ark/3._ark*sqrt(3._ark)*y8*y9)*y6**2*y3**2+(-2._ark/ &
      3._ark*sqrt(3._ark)*y7*y9+2._ark/3._ark*sqrt(3._ark)*y8*y9)*y6**2*y4**2
    dF(3,600) = 2._ark/3._ark*sqrt(3._ark)*y4**2*y5**2*y7*y8+2._ark/ &
      3._ark*sqrt(3._ark)*y3**2*y5**2*y7*y8+2._ark/3._ark*sqrt(3._ark)*y2**2*y5**2*y7*y8+2._ark/ &
      3._ark*sqrt(3._ark)*y1**2*y5**2*y7*y8
    dF(3,601) = ((-3._ark*y7-3._ark*y8)*y5**3+(y8+y7)*y6**2*y5)*y1**2+((3._ark*y8+ &
      3._ark*y7)*y5**3+(-y7-y8)*y6**2*y5)*y2**2+((-3._ark*y8+3._ark*y7)*y5**3+(y8- &
      y7)*y6**2*y5)*y3**2+((3._ark*y8-3._ark*y7)*y5**3+(-y8+y7)*y6**2*y5)*y4**2
    dF(3,602) = (-sqrt(3._ark)*y6**4/6._ark-sqrt(3._ark)*y5**4/2._ark)*y1**2+(- &
      sqrt(3._ark)*y6**4/6._ark-sqrt(3._ark)*y5**4/2._ark)*y2**2+(sqrt(3._ark)*y6**4/6._ark+ &
      sqrt(3._ark)*y5**4/2._ark)*y3**2+(sqrt(3._ark)*y6**4/6._ark+sqrt(3._ark)*y5**4/2._ark)*y4**2
    dF(3,603) = (y6**2*y9**2+y5**2*y9**2)*y1**2+(y6**2*y9**2+y5**2*y9**2)*y2**2+(- &
      y6**2*y9**2-y5**2*y9**2)*y3**2+(-y6**2*y9**2-y5**2*y9**2)*y4**2
    dF(3,604) = ((y7*y9+y8*y9)*y5**2+(y7*y9+y8*y9)*y6**2)*y1**2+((-y7*y9- &
      y8*y9)*y5**2+(-y7*y9-y8*y9)*y6**2)*y2**2+((-y8*y9+y7*y9)*y5**2+(-y8*y9+ &
      y7*y9)*y6**2)*y3**2+((-y7*y9+y8*y9)*y5**2+(-y7*y9+y8*y9)*y6**2)*y4**2
    dF(3,605) = ((sqrt(3._ark)*y8+sqrt(3._ark)*y7)*y5**3+(-y8+y7)*y6*y5**2)*y1**2+((- &
      sqrt(3._ark)*y7-sqrt(3._ark)*y8)*y5**3+(y8-y7)*y6*y5**2)*y2**2+((sqrt(3._ark)*y8- &
      sqrt(3._ark)*y7)*y5**3+(-y7-y8)*y6*y5**2)*y3**2+((sqrt(3._ark)*y7- &
      sqrt(3._ark)*y8)*y5**3+(y8+y7)*y6*y5**2)*y4**2
    dF(3,606) = (-sqrt(3._ark)*y5**4/6._ark-sqrt(3._ark)*y6**4/2._ark)*y1**2+(- &
      sqrt(3._ark)*y5**4/6._ark-sqrt(3._ark)*y6**4/2._ark)*y2**2+(sqrt(3._ark)*y6**4/2._ark+ &
      sqrt(3._ark)*y5**4/6._ark)*y3**2+(sqrt(3._ark)*y6**4/2._ark+sqrt(3._ark)*y5**4/6._ark)*y4**2
    dF(3,607) = y3**2*y4*y9**3+y1*y2**2*y9**3+y1**2*y2*y9**3+y3*y4**2*y9**3
    dF(3,608) = (y3*y7**3+y4*y8**3)*y1**2+(-y4**2*y8**3-y3**2*y7**3)*y1+(-y4*y7**3- &
      y3*y8**3)*y2**2+(y4**2*y7**3+y3**2*y8**3)*y2
    dF(3,609) = -y3*y4**2*y7*y8*y9+y1**2*y2*y7*y8*y9+y1*y2**2*y7*y8*y9- &
      y3**2*y4*y7*y8*y9
    dF(3,610) = (y4*y7*y8*y9+y3*y7*y8*y9)*y1**2+(-y4**2*y7*y8*y9-y3**2*y7*y8*y9)*y1+ &
      (y4*y7*y8*y9+y3*y7*y8*y9)*y2**2+(-y4**2*y7*y8*y9-y3**2*y7*y8*y9)*y2
    dF(3,611) = (y3*y7**2*y8+y4*y7*y8**2)*y1**2+(y4**2*y7*y8**2+y3**2*y7**2*y8)*y1+ &
      (-y3*y7*y8**2-y4*y7**2*y8)*y2**2+(-y3**2*y7*y8**2-y4**2*y7**2*y8)*y2
    dF(3,612) = -2._ark/3._ark*sqrt(3._ark)*y1*y2**2*y5*y7*y8-2._ark/ &
      3._ark*sqrt(3._ark)*y1**2*y2*y5*y7*y8-2._ark/3._ark*sqrt(3._ark)*y3**2*y4*y5*y7*y8-2._ark/ &
      3._ark*sqrt(3._ark)*y3*y4**2*y5*y7*y8
    dF(3,613) = ((3._ark/4._ark*y8+3._ark/4._ark*y7)*y5**2+(sqrt(3._ark)*y8/2._ark- &
      sqrt(3._ark)*y7/2._ark)*y6*y5+(y8/4._ark+y7/4._ark)*y6**2)*y2*y1**2+((-3._ark/4._ark*y7- &
      3._ark/4._ark*y8)*y5**2+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y6*y5+(-y7/4._ark-y8/ &
      4._ark)*y6**2)*y2**2*y1+((-3._ark/4._ark*y7+3._ark/4._ark*y8)*y5**2+(sqrt(3._ark)*y7/2._ark+ &
      sqrt(3._ark)*y8/2._ark)*y6*y5+(-y7/4._ark+y8/4._ark)*y6**2)*y4*y3**2+((3._ark/4._ark*y7- &
      3._ark/4._ark*y8)*y5**2+(-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y6*y5+(y7/4._ark-y8/ &
      4._ark)*y6**2)*y4**2*y3
    dF(3,614) = (y3*y5*y8*y9+y4*y5*y7*y9)*y1**2+(-y4**2*y5*y7*y9-y3**2*y5*y8*y9)*y1+ &
      (-y4*y5*y8*y9-y3*y5*y7*y9)*y2**2+(y4**2*y5*y8*y9+y3**2*y5*y7*y9)*y2
    dF(3,615) = (y5*y6**2-y5**3/3._ark)*y2*y1**2+(y5*y6**2-y5**3/3._ark)*y2**2*y1+ &
      (y5**3/3._ark-y5*y6**2)*y4*y3**2+(y5**3/3._ark-y5*y6**2)*y4**2*y3
    dF(3,616) = ((y8**2+y7**2)*y5+(sqrt(3._ark)*y8**2-sqrt(3._ark)*y7**2)*y6)*y2*y1**2+ &
      ((y8**2+y7**2)*y5+(sqrt(3._ark)*y8**2-sqrt(3._ark)*y7**2)*y6)*y2**2*y1+((-y8**2- &
      y7**2)*y5+(-sqrt(3._ark)*y8**2+sqrt(3._ark)*y7**2)*y6)*y4*y3**2+((-y8**2-y7**2)*y5+ &
      (-sqrt(3._ark)*y8**2+sqrt(3._ark)*y7**2)*y6)*y4**2*y3
    dF(3,617) = ((y5*y7*y9+sqrt(3._ark)*y6*y7*y9)*y3+(y5*y8*y9- &
      sqrt(3._ark)*y6*y8*y9)*y4)*y1**2+((y5*y7*y9+sqrt(3._ark)*y6*y7*y9)*y3**2+(y5*y8*y9- &
      sqrt(3._ark)*y6*y8*y9)*y4**2)*y1+((-y5*y8*y9+sqrt(3._ark)*y6*y8*y9)*y3+(- &
      sqrt(3._ark)*y6*y7*y9-y5*y7*y9)*y4)*y2**2+((-y5*y8*y9+sqrt(3._ark)*y6*y8*y9)*y3**2+ &
      (-sqrt(3._ark)*y6*y7*y9-y5*y7*y9)*y4**2)*y2
    dF(3,618) = (-2._ark/3._ark*sqrt(3._ark)*y3*y6**2*y8-2._ark/ &
      3._ark*sqrt(3._ark)*y4*y6**2*y7)*y1**2+(-2._ark/3._ark*sqrt(3._ark)*y3**2*y6**2*y8-2._ark/ &
      3._ark*sqrt(3._ark)*y4**2*y6**2*y7)*y1+(2._ark/3._ark*sqrt(3._ark)*y4*y6**2*y8+2._ark/ &
      3._ark*sqrt(3._ark)*y3*y6**2*y7)*y2**2+(2._ark/3._ark*sqrt(3._ark)*y4**2*y6**2*y8+2._ark/ &
      3._ark*sqrt(3._ark)*y3**2*y6**2*y7)*y2
    dF(3,619) = ((y5*y6*y9+sqrt(3._ark)*y6**2*y9)*y3+(sqrt(3._ark)*y6**2*y9- &
      y5*y6*y9)*y4)*y1**2+((y5*y6*y9+sqrt(3._ark)*y6**2*y9)*y3**2+(sqrt(3._ark)*y6**2*y9- &
      y5*y6*y9)*y4**2)*y1+((sqrt(3._ark)*y6**2*y9-y5*y6*y9)*y3+(y5*y6*y9+ &
      sqrt(3._ark)*y6**2*y9)*y4)*y2**2+((sqrt(3._ark)*y6**2*y9-y5*y6*y9)*y3**2+(y5*y6*y9+ &
      sqrt(3._ark)*y6**2*y9)*y4**2)*y2
    dF(3,620) = ((-sqrt(3._ark)*y8/4._ark-sqrt(3._ark)*y7/4._ark)*y5**2+(y8/2._ark-y7/ &
      2._ark)*y6*y5+(sqrt(3._ark)*y7/4._ark+sqrt(3._ark)*y8/4._ark)*y6**2)*y2*y1**2+ &
      ((sqrt(3._ark)*y7/4._ark+sqrt(3._ark)*y8/4._ark)*y5**2+(-y8/2._ark+y7/2._ark)*y6*y5+(- &
      sqrt(3._ark)*y8/4._ark-sqrt(3._ark)*y7/4._ark)*y6**2)*y2**2*y1+((-sqrt(3._ark)*y8/4._ark+ &
      sqrt(3._ark)*y7/4._ark)*y5**2+(y7/2._ark+y8/2._ark)*y6*y5+(-sqrt(3._ark)*y7/4._ark+ &
      sqrt(3._ark)*y8/4._ark)*y6**2)*y4*y3**2+((-sqrt(3._ark)*y7/4._ark+sqrt(3._ark)*y8/ &
      4._ark)*y5**2+(-y8/2._ark-y7/2._ark)*y6*y5+(-sqrt(3._ark)*y8/4._ark+sqrt(3._ark)*y7/ &
      4._ark)*y6**2)*y4**2*y3
    dF(3,621) = ((y5*y6*y8+sqrt(3._ark)*y6**2*y8/3._ark)*y3+(sqrt(3._ark)*y6**2*y7/3._ark- &
      y5*y6*y7)*y4)*y1**2+((y5*y6*y8+sqrt(3._ark)*y6**2*y8/3._ark)*y3**2+ &
      (sqrt(3._ark)*y6**2*y7/3._ark-y5*y6*y7)*y4**2)*y1+((-sqrt(3._ark)*y6**2*y7/3._ark+ &
      y5*y6*y7)*y3+(-sqrt(3._ark)*y6**2*y8/3._ark-y5*y6*y8)*y4)*y2**2+((- &
      sqrt(3._ark)*y6**2*y7/3._ark+y5*y6*y7)*y3**2+(-sqrt(3._ark)*y6**2*y8/3._ark- &
      y5*y6*y8)*y4**2)*y2
    dF(3,622) = ((4._ark/3._ark*y5*y6**2-4._ark/9._ark*sqrt(3._ark)*y5**2*y6)*y3+(4._ark/ &
      3._ark*y5*y6**2+4._ark/9._ark*sqrt(3._ark)*y5**2*y6)*y4)*y1**2+((-4._ark/3._ark*y5*y6**2+ &
      4._ark/9._ark*sqrt(3._ark)*y5**2*y6)*y3**2+(-4._ark/3._ark*y5*y6**2-4._ark/ &
      9._ark*sqrt(3._ark)*y5**2*y6)*y4**2)*y1+((4._ark/3._ark*y5*y6**2+4._ark/ &
      9._ark*sqrt(3._ark)*y5**2*y6)*y3+(4._ark/3._ark*y5*y6**2-4._ark/ &
      9._ark*sqrt(3._ark)*y5**2*y6)*y4)*y2**2+((-4._ark/3._ark*y5*y6**2-4._ark/ &
      9._ark*sqrt(3._ark)*y5**2*y6)*y3**2+(-4._ark/3._ark*y5*y6**2+4._ark/ &
      9._ark*sqrt(3._ark)*y5**2*y6)*y4**2)*y2
    dF(3,623) = ((y5**2*y9-3._ark*y6**2*y9)*y3+(y5**2*y9-3._ark*y6**2*y9)*y4)*y1**2+ &
      ((y5**2*y9-3._ark*y6**2*y9)*y3**2+(y5**2*y9-3._ark*y6**2*y9)*y4**2)*y1+((y5**2*y9- &
      3._ark*y6**2*y9)*y3+(y5**2*y9-3._ark*y6**2*y9)*y4)*y2**2+((y5**2*y9- &
      3._ark*y6**2*y9)*y3**2+(y5**2*y9-3._ark*y6**2*y9)*y4**2)*y2
    dF(3,624) = ((y8/4._ark+y7/4._ark)*y5**2+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
      2._ark)*y6*y5+(3._ark/4._ark*y8+3._ark/4._ark*y7)*y6**2)*y2*y1**2+((-y7/4._ark-y8/ &
      4._ark)*y5**2+(sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y7/2._ark)*y6*y5+(-3._ark/4._ark*y7-3._ark/ &
      4._ark*y8)*y6**2)*y2**2*y1+((-y7/4._ark+y8/4._ark)*y5**2+(-sqrt(3._ark)*y7/2._ark- &
      sqrt(3._ark)*y8/2._ark)*y6*y5+(-3._ark/4._ark*y7+3._ark/4._ark*y8)*y6**2)*y4*y3**2+((y7/ &
      4._ark-y8/4._ark)*y5**2+(sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y6*y5+(3._ark/4._ark*y7- &
      3._ark/4._ark*y8)*y6**2)*y4**2*y3
    dF(3,625) = (y3**2*y7*y8+y4**2*y7*y8)*y1**2+(y3**2*y7*y8+y4**2*y7*y8)*y2**2
    dF(3,626) = (-sqrt(3._ark)*y5**2/6._ark+sqrt(3._ark)*y6**2/2._ark)*y2**2*y1**2+(- &
      sqrt(3._ark)*y6**2/2._ark+sqrt(3._ark)*y5**2/6._ark)*y4**2*y3**2
    dF(3,627) = (y8**2*y9+y7**2*y9)*y2*y1**2+(y8**2*y9+y7**2*y9)*y2**2*y1+(y8**2*y9+ &
      y7**2*y9)*y4*y3**2+(y8**2*y9+y7**2*y9)*y4**2*y3
    dF(3,628) = ((-sqrt(3._ark)*y5*y6*y7+y5**2*y7)*y3+(y5**2*y8+ &
      sqrt(3._ark)*y5*y6*y8)*y4)*y1**2+((sqrt(3._ark)*y5*y6*y7-y5**2*y7)*y3**2+(- &
      sqrt(3._ark)*y5*y6*y8-y5**2*y8)*y4**2)*y1+((-sqrt(3._ark)*y5*y6*y8-y5**2*y8)*y3+ &
      (sqrt(3._ark)*y5*y6*y7-y5**2*y7)*y4)*y2**2+((y5**2*y8+sqrt(3._ark)*y5*y6*y8)*y3**2+ &
      (-sqrt(3._ark)*y5*y6*y7+y5**2*y7)*y4**2)*y2
    dF(3,629) = (y6**2*y9+y5**2*y9)*y2*y1**2+(y6**2*y9+y5**2*y9)*y2**2*y1+(y6**2*y9+ &
      y5**2*y9)*y4*y3**2+(y6**2*y9+y5**2*y9)*y4**2*y3
    dF(3,630) = (y3*y7*y9+y4*y8*y9)*y2*y1**2+((-y3*y8*y9-y4*y7*y9)*y2**2+ &
      y3**2*y4*y7*y9+y3*y4**2*y8*y9)*y1+(-y3*y4**2*y7*y9-y3**2*y4*y8*y9)*y2
    dF(3,631) = (-y3*y6*y8+y4*y6*y7)*y2*y1**2+((y4*y6*y8-y3*y6*y7)*y2**2+ &
      y3*y4**2*y6*y7-y3**2*y4*y6*y8)*y1+(y3*y4**2*y6*y8-y3**2*y4*y6*y7)*y2
    dF(3,632) = ((y8+y7)*y5+(sqrt(3._ark)*y8-sqrt(3._ark)*y7)*y6)*y4*y3*y1**2+(((y8- &
      y7)*y5+(sqrt(3._ark)*y8+sqrt(3._ark)*y7)*y6)*y3**2+((-y8+y7)*y5+(-sqrt(3._ark)*y7- &
      sqrt(3._ark)*y8)*y6)*y4**2)*y2*y1+((-y7-y8)*y5+(sqrt(3._ark)*y7- &
      sqrt(3._ark)*y8)*y6)*y4*y3*y2**2
    dF(3,633) = (-2._ark*y3*y5*y9-2._ark*y4*y5*y9)*y2*y1**2+((-2._ark*y3*y5*y9- &
      2._ark*y4*y5*y9)*y2**2-2._ark*y3**2*y4*y5*y9-2._ark*y3*y4**2*y5*y9)*y1+(- &
      2._ark*y3**2*y4*y5*y9-2._ark*y3*y4**2*y5*y9)*y2
    dF(3,634) = -2._ark*y1**2*y3*y4*y5*y9+(-2._ark*y3**2*y5*y9-2._ark*y4**2*y5*y9)*y2*y1- &
      2._ark*y2**2*y3*y4*y5*y9
    dF(3,635) = (y3*y5*y8+y4*y5*y7)*y2*y1**2+((-y3*y5*y7-y4*y5*y8)*y2**2+ &
      y3*y4**2*y5*y7+y3**2*y4*y5*y8)*y1+(-y3**2*y4*y5*y7-y3*y4**2*y5*y8)*y2
    dF(3,636) = ((-y5*y6-sqrt(3._ark)*y5**2/3._ark)*y3+(y5*y6-sqrt(3._ark)*y5**2/ &
      3._ark)*y4)*y2*y1**2+(((y5*y6-sqrt(3._ark)*y5**2/3._ark)*y3+(-y5*y6-sqrt(3._ark)*y5**2/ &
      3._ark)*y4)*y2**2+(y5*y6+sqrt(3._ark)*y5**2/3._ark)*y4*y3**2+(sqrt(3._ark)*y5**2/3._ark- &
      y5*y6)*y4**2*y3)*y1+((sqrt(3._ark)*y5**2/3._ark-y5*y6)*y4*y3**2+(y5*y6+ &
      sqrt(3._ark)*y5**2/3._ark)*y4**2*y3)*y2
    dF(3,637) = (y6**2+y5**2)*y4*y3*y1**2+((-y6**2-y5**2)*y3**2+(-y6**2- &
      y5**2)*y4**2)*y2*y1+(y6**2+y5**2)*y4*y3*y2**2
    dF(3,638) = ((2._ark*y3**2*y6-2._ark*y4**2*y6)*y2+2._ark*y3*y4**2*y6- &
      2._ark*y3**2*y4*y6)*y1**2+(2._ark*y4**2*y6-2._ark*y3**2*y6)*y2**2*y1+(- &
      2._ark*y3*y4**2*y6+2._ark*y3**2*y4*y6)*y2**2
    dF(3,639) = (y4+y3)*y2**2*y1**3+(y4+y3)*y2**3*y1**2+(-y3**3*y4**2- &
      y3**2*y4**3)*y1+(-y3**3*y4**2-y3**2*y4**3)*y2
    dF(3,640) = (-y3*y4**2-y3**2*y4)*y1**3+(y4**3+y3**3)*y2*y1**2+(y4**3+ &
      y3**3)*y2**2*y1+(-y3*y4**2-y3**2*y4)*y2**3
    dF(3,641) = ((y6*y9+sqrt(3._ark)*y5*y9)*y3+(-y6*y9+sqrt(3._ark)*y5*y9)*y4)*y2*y1**2+ &
      (((-y6*y9+sqrt(3._ark)*y5*y9)*y3+(y6*y9+sqrt(3._ark)*y5*y9)*y4)*y2**2+(y6*y9+ &
      sqrt(3._ark)*y5*y9)*y4*y3**2+(-y6*y9+sqrt(3._ark)*y5*y9)*y4**2*y3)*y1+((-y6*y9+ &
      sqrt(3._ark)*y5*y9)*y4*y3**2+(y6*y9+sqrt(3._ark)*y5*y9)*y4**2*y3)*y2
    dF(3,642) = ((y6**2+y5**2)*y3+(y6**2+y5**2)*y4)*y2*y1**2+(((y6**2+y5**2)*y3+ &
      (y6**2+y5**2)*y4)*y2**2+(-y6**2-y5**2)*y4*y3**2+(-y6**2-y5**2)*y4**2*y3)*y1+((- &
      y6**2-y5**2)*y4*y3**2+(-y6**2-y5**2)*y4**2*y3)*y2
    dF(3,643) = ((sqrt(3._ark)*y6*y7+y5*y7)*y3+(-sqrt(3._ark)*y6*y8+y5*y8)*y4)*y2*y1**2+ &
      (((sqrt(3._ark)*y6*y8-y5*y8)*y3+(-y5*y7-sqrt(3._ark)*y6*y7)*y4)*y2**2+(-y5*y7- &
      sqrt(3._ark)*y6*y7)*y4*y3**2+(sqrt(3._ark)*y6*y8-y5*y8)*y4**2*y3)*y1+((- &
      sqrt(3._ark)*y6*y8+y5*y8)*y4*y3**2+(sqrt(3._ark)*y6*y7+y5*y7)*y4**2*y3)*y2
    dF(3,644) = y1**2*y2**2*y3*y4-y1*y2*y3**2*y4**2
    dF(3,645) = -y3**2*y4**2*y9**2+y1**2*y2**2*y9**2
    dF(3,646) = ((y4**2+y3**2)*y2**2-y3**2*y4**2)*y1**2-y2**2*y3**2*y4**2
    dF(3,647) = 2._ark/3._ark*sqrt(3._ark)*y2**3*y5*y9**2+2._ark/ &
      3._ark*sqrt(3._ark)*y1**3*y5*y9**2-2._ark/3._ark*sqrt(3._ark)*y4**3*y5*y9**2-2._ark/ &
      3._ark*sqrt(3._ark)*y3**3*y5*y9**2
    dF(3,648) = y4**3*y5*y7*y8+y3**3*y5*y7*y8+y2**3*y5*y7*y8+y1**3*y5*y7*y8
    dF(3,649) = (sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y1**3+ &
      (sqrt(3._ark)*y6**2*y9/6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y2**3+(sqrt(3._ark)*y6**2*y9/ &
      6._ark-sqrt(3._ark)*y5**2*y9/2._ark)*y3**3+(sqrt(3._ark)*y6**2*y9/6._ark- &
      sqrt(3._ark)*y5**2*y9/2._ark)*y4**3
    dF(3,650) = ((-y8+y7)*y6*y5+(sqrt(3._ark)*y8+sqrt(3._ark)*y7)*y6**2)*y1**3+((y8- &
      y7)*y6*y5+(-sqrt(3._ark)*y7-sqrt(3._ark)*y8)*y6**2)*y2**3+((-y7-y8)*y6*y5+ &
      (sqrt(3._ark)*y8-sqrt(3._ark)*y7)*y6**2)*y3**3+((y8+y7)*y6*y5+(sqrt(3._ark)*y7- &
      sqrt(3._ark)*y8)*y6**2)*y4**3
    dF(3,651) = -8._ark/3._ark*sqrt(3._ark)*y1**3*y5*y6**2+8._ark/ &
      3._ark*sqrt(3._ark)*y3**3*y5*y6**2+8._ark/3._ark*sqrt(3._ark)*y4**3*y5*y6**2-8._ark/ &
      3._ark*sqrt(3._ark)*y2**3*y5*y6**2
    dF(3,652) = y3*y4**3*y7*y8+y3**3*y4*y7*y8+y1*y2**3*y7*y8+y1**3*y2*y7*y8
    dF(3,653) = (y4*y7*y9+y3*y8*y9)*y1**3+(-y3**3*y8*y9-y4**3*y7*y9)*y1+(-y4*y8*y9- &
      y3*y7*y9)*y2**3+(y3**3*y7*y9+y4**3*y8*y9)*y2
    dF(3,654) = (y7*y9+y8*y9)*y2*y1**3+(-y7*y9-y8*y9)*y2**3*y1+(-y8*y9+ &
      y7*y9)*y4*y3**3+(-y7*y9+y8*y9)*y4**3*y3
    dF(3,655) = (y3*y8**2+y4*y7**2)*y1**3+(-y4**3*y7**2-y3**3*y8**2)*y1+(y3*y7**2+ &
      y4*y8**2)*y2**3+(-y4**3*y8**2-y3**3*y7**2)*y2
    dF(3,656) = (y8**2+y7**2)*y2*y1**3+(y8**2+y7**2)*y2**3*y1+(-y8**2- &
      y7**2)*y4*y3**3+(-y8**2-y7**2)*y4**3*y3
    dF(3,657) = ((sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y8/2._ark-y7/ &
      2._ark)*y6)*y2*y1**3+((-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y5+(-y8/2._ark+y7/ &
      2._ark)*y6)*y2**3*y1+((sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y7/2._ark)*y5+(y7/2._ark+y8/ &
      2._ark)*y6)*y4*y3**3+((sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y5+(-y8/2._ark-y7/ &
      2._ark)*y6)*y4**3*y3
    dF(3,658) = ((y6*y9/2._ark-sqrt(3._ark)*y5*y9/2._ark)*y3+(-sqrt(3._ark)*y5*y9/2._ark- &
      y6*y9/2._ark)*y4)*y1**3+((y6*y9/2._ark-sqrt(3._ark)*y5*y9/2._ark)*y3**3+(- &
      sqrt(3._ark)*y5*y9/2._ark-y6*y9/2._ark)*y4**3)*y1+((-sqrt(3._ark)*y5*y9/2._ark-y6*y9/ &
      2._ark)*y3+(y6*y9/2._ark-sqrt(3._ark)*y5*y9/2._ark)*y4)*y2**3+((-sqrt(3._ark)*y5*y9/2._ark- &
      y6*y9/2._ark)*y3**3+(y6*y9/2._ark-sqrt(3._ark)*y5*y9/2._ark)*y4**3)*y2
    dF(3,659) = (-y4*y6*y7+y3*y6*y8)*y1**3+(-y4**3*y6*y7+y3**3*y6*y8)*y1+(-y4*y6*y8+ &
      y3*y6*y7)*y2**3+(-y4**3*y6*y8+y3**3*y6*y7)*y2
    dF(3,660) = ((sqrt(3._ark)*y5*y7+y6*y7)*y3+(sqrt(3._ark)*y5*y8-y6*y8)*y4)*y1**3+((- &
      y6*y7-sqrt(3._ark)*y5*y7)*y3**3+(-sqrt(3._ark)*y5*y8+y6*y8)*y4**3)*y1+((- &
      sqrt(3._ark)*y5*y8+y6*y8)*y3+(-y6*y7-sqrt(3._ark)*y5*y7)*y4)*y2**3+ &
      ((sqrt(3._ark)*y5*y8-y6*y8)*y3**3+(sqrt(3._ark)*y5*y7+y6*y7)*y4**3)*y2
    dF(3,661) = ((-y8/2._ark-y7/2._ark)*y5+(sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y7/ &
      2._ark)*y6)*y2*y1**3+((y7/2._ark+y8/2._ark)*y5+(sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
      2._ark)*y6)*y2**3*y1+((-y8/2._ark+y7/2._ark)*y5+(sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/ &
      2._ark)*y6)*y4*y3**3+((y8/2._ark-y7/2._ark)*y5+(-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/ &
      2._ark)*y6)*y4**3*y3
    dF(3,662) = ((-y5*y9/2._ark-sqrt(3._ark)*y6*y9/2._ark)*y3+(sqrt(3._ark)*y6*y9/2._ark- &
      y5*y9/2._ark)*y4)*y1**3+((-y5*y9/2._ark-sqrt(3._ark)*y6*y9/2._ark)*y3**3+ &
      (sqrt(3._ark)*y6*y9/2._ark-y5*y9/2._ark)*y4**3)*y1+((sqrt(3._ark)*y6*y9/2._ark-y5*y9/ &
      2._ark)*y3+(-y5*y9/2._ark-sqrt(3._ark)*y6*y9/2._ark)*y4)*y2**3+((sqrt(3._ark)*y6*y9/2._ark- &
      y5*y9/2._ark)*y3**3+(-y5*y9/2._ark-sqrt(3._ark)*y6*y9/2._ark)*y4**3)*y2
    dF(3,663) = (y3*y5*y8+y4*y5*y7)*y1**3+(y4**3*y5*y7+y3**3*y5*y8)*y1+(-y3*y5*y7- &
      y4*y5*y8)*y2**3+(-y4**3*y5*y8-y3**3*y5*y7)*y2
    dF(3,664) = 2._ark/3._ark*sqrt(3._ark)*y3*y4**3*y5**2+2._ark/ &
      3._ark*sqrt(3._ark)*y3**3*y4*y5**2-2._ark/3._ark*sqrt(3._ark)*y1*y2**3*y5**2-2._ark/ &
      3._ark*sqrt(3._ark)*y1**3*y2*y5**2
    dF(3,665) = ((y5*y6+sqrt(3._ark)*y5**2/2._ark+sqrt(3._ark)*y6**2/6._ark)*y3+ &
      (sqrt(3._ark)*y5**2/2._ark-y5*y6+sqrt(3._ark)*y6**2/6._ark)*y4)*y1**3+((- &
      sqrt(3._ark)*y6**2/6._ark-y5*y6-sqrt(3._ark)*y5**2/2._ark)*y3**3+(-sqrt(3._ark)*y5**2/ &
      2._ark+y5*y6-sqrt(3._ark)*y6**2/6._ark)*y4**3)*y1+((sqrt(3._ark)*y5**2/2._ark-y5*y6+ &
      sqrt(3._ark)*y6**2/6._ark)*y3+(y5*y6+sqrt(3._ark)*y5**2/2._ark+sqrt(3._ark)*y6**2/ &
      6._ark)*y4)*y2**3+((-sqrt(3._ark)*y5**2/2._ark+y5*y6-sqrt(3._ark)*y6**2/6._ark)*y3**3+(- &
      sqrt(3._ark)*y6**2/6._ark-y5*y6-sqrt(3._ark)*y5**2/2._ark)*y4**3)*y2
    dF(3,666) = (y8+y7)*y4*y3*y1**3+((y8-y7)*y3**3+(-y8+y7)*y4**3)*y2*y1+(-y7- &
      y8)*y4*y3*y2**3
    dF(3,667) = -y3**3*y4*y9**2-y3*y4**3*y9**2+y1**3*y2*y9**2+y1*y2**3*y9**2
    dF(3,668) = y1**3*y2*y3*y4+(y2**3*y3*y4+(-y3**3*y4-y3*y4**3)*y2)*y1
    dF(3,669) = (y7*y9+y8*y9)*y1**4+(-y7*y9-y8*y9)*y2**4+(-y8*y9+y7*y9)*y3**4+(- &
      y7*y9+y8*y9)*y4**4
    dF(3,670) = y3**4*y7*y8+y2**4*y7*y8+y1**4*y7*y8+y4**4*y7*y8
    dF(3,671) = -y3**4*y9**2+y2**4*y9**2+y1**4*y9**2-y4**4*y9**2
    dF(3,672) = 2._ark/3._ark*sqrt(3._ark)*y4**4*y5*y9+2._ark/3._ark*sqrt(3._ark)*y3**4*y5*y9+ &
      2._ark/3._ark*sqrt(3._ark)*y2**4*y5*y9+2._ark/3._ark*sqrt(3._ark)*y1**4*y5*y9
    dF(3,673) = ((-sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y5+(-y8/2._ark+y7/ &
      2._ark)*y6)*y1**4+((sqrt(3._ark)*y7/2._ark+sqrt(3._ark)*y8/2._ark)*y5+(y8/2._ark-y7/ &
      2._ark)*y6)*y2**4+((sqrt(3._ark)*y7/2._ark-sqrt(3._ark)*y8/2._ark)*y5+(-y8/2._ark-y7/ &
      2._ark)*y6)*y3**4+((sqrt(3._ark)*y8/2._ark-sqrt(3._ark)*y7/2._ark)*y5+(y7/2._ark+y8/ &
      2._ark)*y6)*y4**4
    dF(3,674) = (y4*y8+y3*y7)*y1**4+(-y4**4*y8-y3**4*y7)*y1+(-y4*y7-y3*y8)*y2**4+ &
      (y3**4*y8+y4**4*y7)*y2
    dF(3,675) = y3**4*y4*y9+y1**4*y2*y9+y1*y2**4*y9+y3*y4**4*y9
    dF(3,676) = (y4*y9+y3*y9)*y1**4+(y3**4*y9+y4**4*y9)*y1+(y4*y9+y3*y9)*y2**4+ &
      (y3**4*y9+y4**4*y9)*y2
    dF(3,677) = (y3*y8+y4*y7)*y1**4+(y3**4*y8+y4**4*y7)*y1+(-y3*y7-y4*y8)*y2**4+(- &
      y4**4*y8-y3**4*y7)*y2
    dF(3,678) = (y8+y7)*y2*y1**4+(-y7-y8)*y2**4*y1+(y8-y7)*y4*y3**4+(-y8+ &
      y7)*y4**4*y3
    dF(3,679) = -y3*y4**4*y5+y1**4*y2*y5-y3**4*y4*y5+y1*y2**4*y5
    dF(3,680) = y1**2*y2**4+y1**4*y2**2-y3**2*y4**4-y3**4*y4**2
    dF(3,681) = -y4**6+y1**6-y3**6+y2**6
      !
 end subroutine dipch4_diff_mu


end module pot_user
