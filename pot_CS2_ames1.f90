!
!  This unit is for a user defined potential of CH3OH from Bowman (2007)
!
module pot_user
  use accuracy
  use moltype

  implicit none

  public MLdipole,MLpoten,ML_MEP,MLpoten_name

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
   !f = MLdipole_xy4_dF(xyz)
   !
 end subroutine MLdipole
 !
 ! Check the potential name 
 subroutine MLpoten_name(name)
   !
   character(len=cl),intent(in) ::  name
   character(len=cl),parameter ::  poten_name = 'CS2_AMES1'
   ! 
   if (poten_name/=trim(name)) then
     write(out,"(a,a,a,a)") 'Wrong Potential ',trim(name),'; should be ',trim(poten_name)
   endif
   !
   write(out,"(a,a)") '  Using USER-type PES ',trim(poten_name)
   !
 end subroutine MLpoten_name
 !
 ! Defining potential energy function)
 function MLpoten(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   f = MLpoten_cs2_ames1(ncoords,natoms,local,xyz,force)
   !
 end function MLpoten

 
 !
  function MLpoten_cs2_ames1(ncoords,natoms,local,xyz,force) result(f)
   !
   implicit none
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   integer(ik)          :: i,k,Ncoe
   real(ark)            :: alpha1,alpha2,de1,de2,rcs1,rcs2,rcsref,rrefa,rrefb,rcsd1,rcsd2,ang2
   !
   real(ark)            :: emin,etmp2,sumx0,ax,enetmp2,etmp1,enetmp3
   !
   real(ark)            :: str1,str2,enetmp1,atp1,edamp11,edamp12,edamp1,V
   real(ark)            :: v0,r1,r2,a3,xx1,xx2,ang,sumstr2,sumstr4,angref,angx,bdamp2,bdamp4
   real(ark)            :: edp1,edp2,edp3,edp4,edp5,edp6,De_1,De_2,De3,De_3
   !
   real(ark)            ::   v1(3),v2(3),x1,x2,cosalpha1
        !
        rcsref  = force(1)
        alpha1  = force(2)
        alpha2  = alpha1
        De1     = force(3)
        De2     = De1
        De_1    = force(4)
        De_2    = De_1
        De3     = force(5)
        De_3    = force(6)
        edp1    = force(7)
        edp2    = force(8)
        edp3    = force(9)
        edp4    = force(10)
        edp5    = force(11)
        edp6    = force(12)
        Emin    = force(13)*219474.63067_ark
        !
        rcs1 = local(1)
        rcs2 = local(2)
        !
        ang = local(3)
        !
        r1=rcs1-rcsref
        r2=rcs2-rcsref
        !
        a3=1.0_ark+cos(ang)
        !
        Ncoe = molec%parmax
        !
        v0=0
        do i=14,Ncoe
           !
           v0=v0+force(i)*r1**molec%pot_ind(1,i)*r2**molec%pot_ind(2,i)*a3**(molec%pot_ind(3,i))
           if(molec%pot_ind(1,i).ne.molec%pot_ind(2,i))then
             v0=v0+force(i)*r2**molec%pot_ind(1,i)*r1**molec%pot_ind(2,i)*a3**(molec%pot_ind(3,i))
           end if
           !
        end do
        !
        rrefa=rcsref; rrefb=rcsref
        rcsd1=rcs1-rrefa; rcsd2=rcs2-rrefb
        rcsd1=rcsd1*rcsd1; rcsd2=rcsd2*rcsd2
        sumstr2=rcsd1+rcsd2
        sumstr4=rcsd1*rcsd1+rcsd2*rcsd2
        !
        r1=rcs1-rcsref
        r2=rcs2-rcsref
        !
        ang2=pi-ang
        angref=150.d0*pi/180.d0-pi
        angx=min(-abs(ang2)-angref,0.0_ark)
        bdamp2=angx**2;  bdamp4=angx**4
        etmp2=exp(edp1*sumstr2+edp2*sumstr4+edp3*bdamp2+edp4*bdamp4)
        !
        V0=V0*etmp2
        !
        enetmp1=De1*(1.0_ark-exp(-alpha1*(rcs1-rcsref)))**2 + De_1*(1.0_ark-exp(-alpha1*(rcs1-rcsref)))**4
        enetmp2=De2*(1.0_ark-exp(-alpha2*(rcs2-rcsref)))**2 + De_2*(1.0_ark-exp(-alpha2*(rcs2-rcsref)))**4
        !
        ax=(pi-ang)*0.5_ark
        enetmp3=De3*sin(ax)**2 + De_3*sin(ax)**4  !bending simulation
        edamp1=exp(edp5*sumstr2+1.d0*edp6*sumstr4)
        enetmp3=edamp1*enetmp3
        !
        etmp1=enetmp1+enetmp2+enetmp3 
        !
        V=V0-emin+etmp1
        !
        f = V
        !
  end function MLpoten_cs2_ames1

end module pot_user
