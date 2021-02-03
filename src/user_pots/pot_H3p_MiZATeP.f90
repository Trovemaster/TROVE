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
   f = 0
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
   f = MLpoten_xy2_dmbe(ncoords,natoms,local,xyz,force)
   !
 end function MLpoten
 !
 !
  !
  ! Defining potential energy function 
  !
  ! This is a Varandas type (DMBE) PES for XY2 molecules. It simply uses the external routine,
  ! defined in J. Chem. Phys. 118, 2637 (2003).
  !
  function MLpoten_xy2_dmbe(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f
   !
   real(ark)               :: r12,r32,alpha,xcos,v
   real(ark),parameter     :: tocm = 219474.63067_ark
   !
   if (verbose>=6) write(out,"('MLpoten_xy2_tennys/start')") 


     !write (out,"('MLpoten_xy2_dmbe: is turned off')")
     !stop 'MLpoten_xy2_dmbe: PES is turned off'


     r12 = local(1)/bohr ; r32 = local(2)/bohr ;  alpha = local(3)
     !
     if (molec%AtomMasses(1)>15.8_rk.and.molec%AtomMasses(1)<18.0_rk) then
       !
       ! whater
       ! 
       xcos = cos(alpha)
       !
       call potv(v,r12,r32,xcos)
       !v = 0
       v = v*tocm
       !
       !v = v + MLpoten_xy2_bubukina(ncoords,natoms,local,xyz,force)
       !
     elseif(molec%AtomMasses(1)>31.9_rk.and.molec%AtomMasses(1)<36.0_rk) then 
       !
       ! h2s
       !
       !call potv_h2s(v,r12,r32,alpha)
       v = 0 
       !
     else
        !
        write (out,"('MLpoten_xy2_dmbe: for these atoms ',3f12.6,' PES is not provided.')") molec%AtomMasses(1:3)
        stop 'MLpoten_xy2_dmbe: PES is not provided'
        !
     endif
     ! 
     f = v
     !
     if (verbose>=6) write(out,"('MLpoten_xy2_dmbe/end')") 
 
 end function MLpoten_xy2_dmbe


end module pot_user
