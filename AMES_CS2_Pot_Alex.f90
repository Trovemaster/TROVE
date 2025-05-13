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
   real(ark)            :: rssref,rcsref,alpha1,alpha2,De1,De2,De3,edamp1,edamp2,edamp3,edamp4,edamp5,fdamp
   !
   real(ark)            :: emin,rmin,rminbohr,alpha,rref,a2b
   !
   real(ark)            :: str1,str2,enetmp1,atp1,enetmp2,edamp11,edamp12,V
   real(ark)            :: v0,rrcs1,rrcs2,r1,r2,a3,xx1,xx2,rcs1,rcs2,ang,ang2,dstr1,dstr2,sumstr2,sumstr4,angref,angx,bdamp2,bdamp4
   !
   real(ark)            ::   v1(3),v2(3),x1,x2,cosalpha1
   character(len=cl)    :: txt = 'MLpoten_c3_R_theta'
        !
        rssref  = force(1)
        rcsref  = force(2)
        alpha1  = force(3)
        alpha2  = force(4)
        De1     = force(5)
        De2     = force(6)
        De3     = force(7)
        edamp1  = force(8)
        edamp2  = force(9)
        edamp3  = force(10)
        edamp4  = force(11)
        edamp5  = force(12)
        Emin    = force(13)
        !
        rcs1 = local(1)
        rcs2 = local(2)
        !
        rrcs1=rcs1/0.529177249_ark
        rrcs2=rcs2/0.529177249_ark
        !
        !   ang = pi-local(3)
        ang = local(3)
        !
        if  (molec%Ndihedrals>1) then
          !
          ang = asin( sqrt( sin(local(3))**2+sin(local(4))**2 ))
          !
          v1(:) = xyz(2,:)-xyz(1,:)
          v2(:) = xyz(3,:)-xyz(1,:)
          !
          x1 = sqrt(sum(v1(:)**2))
          x2 = sqrt(sum(v2(:)**2))
          !
          cosalpha1 = sum( v1(:)*v2(:))/(x1*x2)
          !
          !alpha = aacos(cosalpha1,txt)
          ang  = aacos(cosalpha1,txt)-pi
          !
        endif
        !
      !   r1=1._ark-exp(-alpha1*(rcs1-rcsref))
      !   r2=1._ark-exp(-alpha2*(rcs2-rcsref))
        r1=rcs1-rcsref
        r2=rcs2-rcsref
        a3=1._ark + cos(ang)
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
        xx1=rcs1; xx2=rcs2
        dstr1=(xx1-rcsref)
        dstr2=(xx2-rcsref)
        sumstr2=dstr1**2+dstr2**2
        sumstr4=dstr1**4+dstr2**4
        !
        angref=150.0_ark*pi/180.0_ark-pi
        !
        ang2=pi-ang
        angx=min(-abs(ang2)-angref,0.0_ark)
        bdamp2=angx**2;bdamp4=angx**4
      !   bdamp2=ang**2; bdamp4=ang**4
        !
        !angx=min(-abs(ang)+angref,0.d0)
        !bdamp2=angx**2;bdamp4=angx**4
        !bdamp2=ang**2; bdamp4=ang**4
        !
        V0=V0*exp(edamp1*sumstr2+edamp2*sumstr4+edamp3*bdamp2+edamp4*bdamp4)
      !   write(*, '(A, F20.10, F20.10, F20.10, F20.10)') "Vshort ", local(1), local(2), local(3), V0
        !
        a2b=0.529177249_ark
        str1=1._ark-exp(-alpha1*(xx1-rcsref))
        str2=1._ark-exp(-alpha2*(xx2-rcsref))
        enetmp1=De1*(str1**2+str2**2)+De2*(str1**4+str2**4)
        enetmp1=enetmp1
        !
        atp1=(acos(-1._ark)-abs(ang))/2
        !   enetmp2=Ae1*(1._ark+cos(atp1))**2+Ae2*(1._ark+cos(atp1))**4
        enetmp2=De3*sin(atp1)**2!+Ae2*(1._ark+cos(atp1))**4
        enetmp2=enetmp2
        !
        sumstr2=(xx1-rcsref)**2+(xx2-rcsref)**2
      !   sumstr4=(xx1-rcsref)**4+(xx2-rcsref)**4
        !   edamp11=-0.5_ark; edamp12=-0.5_ark
        fdamp=exp(edamp5*sumstr2)
        enetmp2=fdamp*enetmp2
        !
        V0=V0+enetmp1+enetmp2
        !
        V=V0/219474.63067_ark-emin
        !
        f = V*219474.63067_ark
      !   write(*, '(A, F20.10, F20.10, F20.10, F20.10)') "enetmp1 ", local(1), local(2), local(3), enetmp1
      !   write(*, '(A, F20.10, F20.10, F20.10, F20.10)') "enetmp2 ", local(1), local(2), local(3), enetmp2
      !   write(*, '(A, F20.10, F20.10, F20.10, F20.10)') "grep ", local(1), local(2), local(3), f
        !
  end function MLpoten_cs2_ames1


end module pot_user
