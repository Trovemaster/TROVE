!
!  This unit is for a user defined potential: N2O PES 
!
module pot_user
  use accuracy
  use moltype
  use pot_xy2, only : MLpoten_xy2_morse_powers,MLpoten_co2_ames1,MLloc2pqr_xyz

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



 !
 ! Defining potential energy function (built for SO2)

 function MLpoten(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force( :)
   real(ark)              ::  f,f1,f2
   integer(ik)            ::  nparam,ntot
   !
   !
   nparam = int(force(1))
   ntot = size(force)
   !
   f1 = MLpoten_co2_ames1(ncoords,natoms,local,xyz,force(2:nparam+1))
   !
   f2 = MLpoten_xy2_morse_powers(ncoords,natoms,local,xyz,force(nparam+2:ntot))
   !
   f = f1+f2
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
    stop 'MLdipole N2O is undefined'
    !
  end subroutine MLdipole






end module pot_user
