!
!  This unit is for a user defined potential 
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
 ! Check the potential name 
 subroutine MLpoten_name(name)
   !
   character(len=cl),intent(in) ::  name
   character(len=cl),parameter ::  poten_name = 'GENERAL'
   ! 
   if (poten_name/=trim(name)) then
     write(out,"(a,a,a,a)") 'Wrong Potential ',trim(name),'; should be ',trim(poten_name)
   endif
   !
   write(out,"(a,a)") '  Using USER-type PES ',trim(poten_name)
   !
 end subroutine MLpoten_name
 !


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
   f = MLpoten_so2_damp(ncoords,natoms,local,xyz,force)
   !
 end function MLpoten
 !

 !
 function MLpoten_so2_damp(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
      !
      integer(ik)          :: iterm,k_ind(3),ifst
      real(ark)            :: y(3),xi(3),vshort,vlong,vdamp
      real(ark)            :: De,ae,re,beta,damp2,damp4,th1,vshort2,temp
      !
      re     = molec%req(1)
      ae     = molec%alphaeq(1)
      !
      De     = force(1)
      beta   = force(2)
      damp2  = force(3)
      damp4  = force(4)
      !
      ifst   = 5
      !
      y(1) = local(1)-re
      y(2) = local(2)-re
      y(3) = local(3)-ae
      !
      vlong = De*(1.0_ark-exp(-beta*y(1)))**2 + De*(1.0_ark-exp(-beta*y(2)))**2
      vdamp = exp( -damp2*( y(1)**2+y(2)**2)-damp4*( y(1)**4+y(2)**4 ) )
      !
      vshort2 = 0.0_ark
      !
      vshort = 0.0_ark
      !
      do iterm=ifst,size(force)
       !
       xi(1:3) = y(1:3)**molec%pot_ind(1:3,iterm)
       !
       temp = force(iterm)*product(xi(1:3))
       !
       if (molec%pot_ind(1,iterm)/=molec%pot_ind(2,iterm)) then
         !
         k_ind(1) = molec%pot_ind(2,iterm)
         k_ind(2) = molec%pot_ind(1,iterm)
         k_ind(3) = molec%pot_ind(3,iterm)
         !
         xi(1:3) = y(1:3)**k_ind(1:3)
         !
         temp = temp + force(iterm)*product(xi(1:3))
         !
       endif
       !
       if (sum(molec%pot_ind(1:3,iterm))>2) then 
         !
         vshort = vshort + temp
         !
       else 
         !
         vshort2 = vshort2 + temp
         !
       endif
       !
      enddo
      !
      temp = vshort2*vdamp+vlong
      !
      th1 = 0.5d0*( 1.0d0-tanh( 0.0001_ark*( temp-40000.0_ark ) ) )
      !
      f = temp + vshort*th1*vdamp
      !
      !f = vlong + vshort*vdamp
      !
  end function MLpoten_so2_damp

end module pot_user
