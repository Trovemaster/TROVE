module prop_xy2_spinrot
  use accuracy
  use moltype
  use pot_xy2, only : MLloc2pqr_xy2

  implicit none

  public prop_xy2_spin_rotation_bisector


contains


subroutine prop_xy2_spin_rotation_bisector(rank, ncoords, natoms, local, xyz, f)

  implicit none 
  !
  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)
  integer(ik) :: iatom
  !
  real(ark) :: xyz0(3), xyz_(natoms,3),r1,r2,alpha,rho,rho_over_sinrho,rho2_over_sinrho2
  real(ark) :: c(3,3), mat(3,3), c_out(5,-2:0), e1(3), e2(3), e3(3),x(natoms,3)
  real(ark),parameter  :: rho_threshold = 0.01_rk
  integer(ik) :: icentre
  !
  if (rank/=9) then
    write(out, '(/a,1x,i3,1x,a)') 'prop_xy2_spin_rotation_bisector: rank of the dipole moment vector =', rank, ', expected 9'
    stop
  endif
  !
  ! xyz are undefined for the local case
  if (all(abs(xyz)<small_)) then
    !
    select case(trim(molec%coords_transform))
    case default
       write (out,"('prop_xy2_qmom_bisect_frame: coord. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'prop_xy2_qmom_bisect_frame - bad coord. type'
    case('R-RHO-Z')
       !
       x = MLloc2pqr_xy2(local)
       !
    end select
    !
  else
    !
    x = xyz
    !
  endif

  xyz0 = x(1,:)
  do iatom=1, natoms
    xyz_(iatom,:) = x(iatom,:) - xyz0(:)
  enddo
  !
  r1 = sqrt(sum(xyz_(2,:)**2))
  r2 = sqrt(sum(xyz_(3,:)**2))
  !
  e1 = xyz_(2,:)/r1
  e2 = xyz_(3,:)/r2
  !
  alpha = aacos(sum(e1*e2))
  !
  rho = pi-alpha
  !
  icentre = extF%coef(2,1) ! = 1 or 2, determines on which Y-atom spin-rotation tensor is centered
  !
  ! fitted tensor elements
  c = 0
  !
  if (icentre==1) then
      c(1,1) = fit_xy2_sr(extF%nterms(1)-2, extF%coef(3:extF%nterms(1),1), (/r1, r2, alpha/))
      c(2,2) = fit_xy2_sr(extF%nterms(2)-2, extF%coef(3:extF%nterms(2),2), (/r1, r2, alpha/))
      c(3,3) = fit_xy2_sr(extF%nterms(3)-2, extF%coef(3:extF%nterms(3),3), (/r1, r2, alpha/))
      c(1,3) = fit_xy2_sr(extF%nterms(4)-2, extF%coef(3:extF%nterms(4),4), (/r1, r2, alpha/))
      c(3,1) = fit_xy2_sr(extF%nterms(5)-2, extF%coef(3:extF%nterms(5),5), (/r1, r2, alpha/))
  elseif (icentre==2) then
      c(1,1) = fit_xy2_sr(extF%nterms(1)-2, extF%coef(3:extF%nterms(1),1), (/r2, r1, alpha/))
      c(2,2) = fit_xy2_sr(extF%nterms(2)-2, extF%coef(3:extF%nterms(2),2), (/r2, r1, alpha/))
      c(3,3) = fit_xy2_sr(extF%nterms(3)-2, extF%coef(3:extF%nterms(3),3), (/r2, r1, alpha/))
      c(1,3) = -fit_xy2_sr(extF%nterms(4)-2, extF%coef(3:extF%nterms(4),4), (/r2, r1, alpha/))
      c(3,1) = -fit_xy2_sr(extF%nterms(5)-2, extF%coef(3:extF%nterms(5),5), (/r2, r1, alpha/))
  else
      write(out, '(a,1x,i3)') 'prop_xy2_qmom_bisect_frame error: icentre /= (1 or 2)'
      stop 'prop_xy2_qmom_bisect_frame - bad icentre parameter' 
  endif
  !
  ! inverse transform
  !
  mat = 0
  !
  !
  if (rho>rho_threshold) then
    !
    ! rho/sin(rho)
    rho_over_sinrho = rho/sin(rho)
    rho2_over_sinrho2 = rho**2/sin(0.5_ark*rho)**2
    !
  else
    ! up to rho^7
    rho_over_sinrho = 1.0_ark+1.0_ark/6.0_ark*rho**2+7.0_ark/360.0_ark*rho**4+31.0_ark/15120.0_ark*rho**6
    rho2_over_sinrho2 = 4.0_ark+1.0_ark/3.0_ark*rho**2+1.0_ark/60.0_ark*rho**4+1.0_ark/1512.0_ark*rho**6
    !
  endif
  !
  mat(1,1) = 1.0_ark/cos(0.5_ark*rho)**2 * (-1.0_ark)/(3.0_ark*(r1-r2)**2-4.0_ark*r1*r2)
  mat(2,2) = 1.0_ark
  mat(3,3) = rho2_over_sinrho2 * (-1.0_ark)/(3.0_ark*(r1-r2)**2-4.0_ark*r1*r2)
  mat(1,3) = rho_over_sinrho * ( 4.0_ark*(r1-r2) )/( (r1+r2)*(3.0_ark*(r1-r2)**2-4.0_ark*r1*r2) )
  mat(3,1) = mat(1,3) 
  !
  !1,1
  c_out(1,0)  = mat(1,1)*c(1,1)
  c_out(1,-1) = mat(1,3)*c(3,1)
  c_out(1,-2) = 0
  !
  !1,3
  c_out(2,0)  = mat(1,1)*c(1,3)
  c_out(2,-1) = mat(1,3)*c(3,3)
  c_out(2,-2) = 0
  !
  !2,2
  c_out(3,0)  = mat(2,2)*c(2,2)
  c_out(3,-1) = 0
  c_out(3,-2) = 0
  !
  !3,1
  c_out(4,0)  = 0
  c_out(4,-1) = mat(1,3)*c(1,1)
  c_out(4,-2) = mat(3,3)*c(3,1)
  !
  !3,3
  c_out(5,0)  = 0
  c_out(5,-1) = mat(1,3)*c(1,3)
  c_out(5,-2) = mat(3,3)*c(3,3)
  !
  f = (/c_out(1,0),c_out(1,-1),c_out(2,0),c_out(2,-1),c_out(3,0),c_out(4,-1),c_out(4,-2),c_out(5,-1),c_out(5,-2)/)
  !
end subroutine prop_xy2_spin_rotation_bisector


!###################################################################################################


function fit_xy2_sr(nparams, params, coords) result(f)

  integer(ik), intent(in) :: nparams 
  real(ark), intent(in) :: params(nparams), coords(3)
  real(ark) :: f

  real(ark) :: y1, y2, y3, alpha, rho, rad, f0, f1, f2, f3, req, alphaeq, beta

  rad = pi/180.0_ark

  req     = params(1)
  alphaeq = params(2)*rad ! obsolete
  beta    = params(3)

  y1     = (coords(1)-req) *exp(-beta*(coords(1)-req)**2)
  y2     = (coords(2)-req) *exp(-beta*(coords(2)-req)**2)
  alpha  = coords(3)
  rho    = pi - alpha
  y3     = 1.0_ark - cos(rho)

  f0 = params(4)
  f1 = xy2_func_n1_d6( (/y1,y2,y3/), params(5:22)  )  ! nparams = 18
  f2 = xy2_func_n2_d6( (/y1,y2,y3/), params(23:67) )  ! nparams = 45
  f3 = xy2_func_n3_d6( (/y1,y2,y3/), params(68:87) )  ! nparams = 20

  f = f0 + f1 + f2 + f3

end function fit_xy2_sr


!###################################################################################################


function xy2_func_n1_d6(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(18)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(18)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = r1
f(2) = r1**2
f(3) = r1**3
f(4) = r1**4
f(5) = r1**5
f(6) = r1**6
f(7) = r2
f(8) = r2**2
f(9) = r2**3
f(10) = r2**4
f(11) = r2**5
f(12) = r2**6
f(13) = a1
f(14) = a1**2
f(15) = a1**3
f(16) = a1**4
f(17) = a1**5
f(18) = a1**6
v = sum(f*params)
end function xy2_func_n1_d6


!###################################################################################################


function xy2_func_n2_d6(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(45)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(45)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = r1*r2
f(2) = r1*r2**2
f(3) = r1*r2**3
f(4) = r1*r2**4
f(5) = r1*r2**5
f(6) = r1**2*r2
f(7) = r1**2*r2**2
f(8) = r1**2*r2**3
f(9) = r1**2*r2**4
f(10) = r1**3*r2
f(11) = r1**3*r2**2
f(12) = r1**3*r2**3
f(13) = r1**4*r2
f(14) = r1**4*r2**2
f(15) = r1**5*r2
f(16) = a1*r1
f(17) = a1**2*r1
f(18) = a1**3*r1
f(19) = a1**4*r1
f(20) = a1**5*r1
f(21) = a1*r1**2
f(22) = a1**2*r1**2
f(23) = a1**3*r1**2
f(24) = a1**4*r1**2
f(25) = a1*r1**3
f(26) = a1**2*r1**3
f(27) = a1**3*r1**3
f(28) = a1*r1**4
f(29) = a1**2*r1**4
f(30) = a1*r1**5
f(31) = a1*r2
f(32) = a1**2*r2
f(33) = a1**3*r2
f(34) = a1**4*r2
f(35) = a1**5*r2
f(36) = a1*r2**2
f(37) = a1**2*r2**2
f(38) = a1**3*r2**2
f(39) = a1**4*r2**2
f(40) = a1*r2**3
f(41) = a1**2*r2**3
f(42) = a1**3*r2**3
f(43) = a1*r2**4
f(44) = a1**2*r2**4
f(45) = a1*r2**5
v = sum(f*params)
end function xy2_func_n2_d6


!###################################################################################################


function xy2_func_n3_d6(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(20)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(20)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = a1*r1*r2
f(2) = a1**2*r1*r2
f(3) = a1**3*r1*r2
f(4) = a1**4*r1*r2
f(5) = a1*r1*r2**2
f(6) = a1**2*r1*r2**2
f(7) = a1**3*r1*r2**2
f(8) = a1*r1*r2**3
f(9) = a1**2*r1*r2**3
f(10) = a1*r1*r2**4
f(11) = a1*r1**2*r2
f(12) = a1**2*r1**2*r2
f(13) = a1**3*r1**2*r2
f(14) = a1*r1**2*r2**2
f(15) = a1**2*r1**2*r2**2
f(16) = a1*r1**2*r2**3
f(17) = a1*r1**3*r2
f(18) = a1**2*r1**3*r2
f(19) = a1*r1**3*r2**2
f(20) = a1*r1**4*r2
v = sum(f*params)
end function xy2_func_n3_d6


end module prop_xy2_spinrot
