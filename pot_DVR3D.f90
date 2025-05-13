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
 !
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
 !
 recursive subroutine MLdipole(rank,ncoords,natoms,local,xyz,f)
   !
   integer(ik),intent(in) ::  rank,ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
   real(ark),intent(out)  ::  f(rank)
   !
   call MLdms2pqr_dvr3d(rank,ncoords,natoms,local,xyz,f)
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
   f = MLpoten_xy2_DVR3D(ncoords,natoms,local,xyz,force)
   !
 end function MLpoten
 !
 !
  !
  ! Defining potential energy function 
  !
  function MLpoten_xy2_DVR3D(ncoords,natoms,local,xyz,force) result(f) 
   !
   integer(ik),intent(in) ::  ncoords,natoms
   real(ark),intent(in)   ::  local(ncoords)
   real(ark),intent(in)   ::  xyz(natoms,3)
   real(ark),intent(in)   ::  force(:)
   real(ark)              ::  f,r1,r2,r3,cosalpha
   !
   real(ark)               :: r12,r32,alpha,xcos,v
   real(ark),parameter     :: tocm = 219474.63067_ark
   !
   if (verbose>=6) write(out,"('MLpoten_DVR3D/start')") 
     !
     select case(trim(molec%coords_transform))
     case default
       !write (out,"('MLpoten_xy2_DVR3D: coord. type ',a,' unknown')") trim(molec%coords_transform)
       !stop 'MLpoten_xy2_DVR3D - bad coord. type'
       !
     !case('R-ALPHA')
       !
       r12 = local(1)/bohr ; r32 = local(2)/bohr ;  alpha = local(3)
       ! 
       xcos = cos(alpha)
       !
     case('R-R-R')
       !
       r1 = local(1) ; r2 = local(2) ;   alpha = local(3)
       !
       !cosalpha = (r1**2+r2**2-r3**2)/(2.0_ark*r1*r2)
       !
       r12 = r1/bohr ; r32 = r2/bohr ;  xcos = cos(alpha)
       !
     end select 
     !
     call potv(v,r12,r32,xcos)
     f = v*tocm
     !
     if (verbose>=6) write(out,"('MLpoten_DVR3D/end')") 
 
 end function MLpoten_xy2_DVR3D


 recursive subroutine MLdms2pqr_dvr3d(rank,ncoords,natoms,local,xyz,f)

    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
    real(ark),intent(out)  ::  f(rank)
    !
    integer(ik)           :: i,ik(1:3)
    real(ark)             :: b(molec%ncoords), mu(3),u1(3),u2(3),u3(3),tmat(3,3),n1(3),n2(3),x(2,3),alpha
    real(rk)              :: r1,r2,xcos,dipp,dipq
    !
    x(1,:) = xyz(2,:) - xyz(1,:)
    x(2,:) = xyz(3,:) - xyz(1,:)
    !
    r1 = sqrt(sum(x(1,:)**2))
    r2 = sqrt(sum(x(2,:)**2))
    !
    n1 = x(1,:) / r1
    n2 = x(2,:) / r2
    !
    alpha = acos(sum(n1*n2))
    !
    u1 = n1 + n2
    u2 = n2 - n1
    !
    u1 = u1 / sqrt(sum(u1(:)**2))
    u2 = u2 / sqrt(sum(u2(:)**2))
    !
    u3 = MLvector_product(u1,u2)
    !
    tmat(1, :) = u1
    tmat(2, :) = u2
    tmat(3, :) = u3
    !
    r1 = local(1)/bohr ; r2 = local(2)/bohr ;  alpha = local(3)
    ! 
    xcos = cos(alpha)
    !
    !call dipd(dipp,r1,r2,xcos,1)
    !
    !call dipd(dipq,r1,r2,xcos,0)
    !
    mu(1) = dipq
    mu(2) = dipp
    mu(3) = 0 
    !
    f(1:3) = matmul(mu,tmat)
    !
 end subroutine MLdms2pqr_dvr3d
 !

end module pot_user
