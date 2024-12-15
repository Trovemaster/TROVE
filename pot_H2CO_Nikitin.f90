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
   character(len=cl),parameter ::  poten_name = 'H2CO_NIKITIN_PES'
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
   f = MLpoten_H2CO_Nikitin_PES(ncoords,natoms,local,xyz,force)
   !
 end function MLpoten
 !
 !!!!start change
  function ML_MEP_OH3P(x)  result(dst)

   real(ark),intent(in)   ::  x
   real(ark)              ::  dst
   real(ark)              ::  rf(0:6),y
   integer(ik)          ::  k(0:6) = (/0,1,2,3,4,5,6/)
     !
     rf(0:6)  = molec%mep_params(1:7) 
     !
     y=sin(x)-1.0_ark
     !
     dst =  sum(rf(:)*y**k(:))
     ! 
  end function ML_MEP_OH3P
 !!!!end change
 !
  ! Defining potential energy function 
  !
function MLpoten_H2CO_Nikitin_PES(ncoords,natoms,local,xyz,force) result(f) 
  !
  integer(ik),intent(in) ::  ncoords,natoms
  real(ark),intent(in)   ::  local(ncoords)
  real(ark),intent(in)   ::  xyz(natoms,3)
  real(ark),intent(in)   ::  force(:)
  real(ark)              ::  f
  !
  integer(ik)          :: iterm,k_ind(6)
  real(ark)            :: a,y(6),S(6),xi(6),req(3),alphaeq(2),cosrho,rho,xieq(6),eq(6), point(6), Spoint(6), ypoint(6), xipoint(6), fpoint
  !
  rho = local(6)
  !
  cosrho = cos(rho) + 1.0_ark

  eq(1) = 1.287759196_ark
  eq(2) = 1.078383659_ark
  eq(3) = 1.078383659_ark
  eq(4) = 119.723652929_ark*pi/180.0_ark
  eq(5) = 119.723652929_ark*pi/180.0_ark

  point(1) = 1.29006214_ark
  point(2) = 1.08000288_ark
  point(3) = 1.08000288_ark
  point(4) = 2.09125828_ark
  point(5) = 2.09125828_ark
  point(6) = 3.14159265_ark
  !
  ! reference
  !
  xieq(1:6) = eq(1:6)
  !
  ! expansion functions
  !
  a = 1.9_ark
  S(1:3) = 1.0_ark-exp(-a*(local(1:3) - eq(1:3)))
  Spoint(1:3) = 1.0_ark-exp(-a*(point(1:3) - eq(1:3)))
  !  y(1:3) = 1.0_ark-exp(-(local(1:3)-xieq(1:3)))
  y(1) = S(2)
  y(2) = S(1)
  y(3) = S(3)
  ypoint(1) = Spoint(2)
  ypoint(2) = Spoint(1)
  ypoint(3) = Spoint(3)
  S(4:5) = sin(local(4:5)-eq(4:5))
  Spoint(4:5) = sin(point(4:5)-eq(4:5))
  !  y(4:5) = local(4:5)-xieq(4:5)
  y(5) = S(4)
  y(6) = S(5)
  y(4) = cos(local(6)/2.0_ark)
  ypoint(5) = Spoint(4)
  ypoint(6) = Spoint(5)
  ypoint(4) = cos(point(6)/2.0_ark)
  !
  !  y(6)   = cosrho
  !
  ! def potential energy
  !
  f = 0.0_ark
  fpoint = 0.0_ark
  !
  do iterm = 1,size(force)
    !
    xi(1:6) = y(1:6)**molec%pot_ind(1:6,iterm)
    xipoint(1:6) = ypoint(1:6)**molec%pot_ind(1:6,iterm)
    !
    f = f + force(iterm)*product(xi(1:6))
    fpoint = fpoint + force(iterm)*product(xipoint(1:6))
    !
  enddo

  f = (f - (-114.5624285449_ark)) * 219474.63_ark
  fpoint = (fpoint - (-114.5624285449_ark)) * 219474.63_ark
  !
end function MLpoten_H2CO_Nikitin_PES


end module pot_user
