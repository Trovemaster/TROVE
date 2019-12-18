module prop_xy2_quad
  use accuracy
  use moltype
  use timer
  use pot_xy2, only : MLloc2pqr_xy2

  implicit none

  public prop_xy2_qmom_bisect_frame, TEST_xy2_qmom_bisect_frame

  private

  integer(ik), parameter :: verbose     = 4                          ! Verbosity level


  contains


!###################################################################################################################################


recursive subroutine prop_xy2_qmom_bisect_frame(rank, ncoords, natoms, local, xyz, f)

  implicit none
  !
  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)

  integer(ik) :: iatom,ierr
  real(ark) :: xyz0(3), xyz_(natoms,3), r1, r2, alpha1, e1(3), e2(3), coords(3), x(3,3), q11_A1, q22_A1, q33_A1, q11_min_q22_A1, q13_B2

  if (rank/=6) then
    write(out, '(/a,1x,i3,1x,a)') 'xy2_qmom_bisect_frame error: rank of the dipole moment vector =', rank, ', expected 6'
    stop
  endif


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

  r1 = sqrt(sum(xyz_(2,:)**2))
  r2 = sqrt(sum(xyz_(3,:)**2))

  e1 = xyz_(2,:)/r1
  e2 = xyz_(3,:)/r2

  alpha1 = acos(sum(e1*e2))

  coords = (/r1,r2,alpha1/)

  q33_A1         = fit_xy2_quad_A1(extF%nterms(1), extF%coef(1:extF%nterms(1),1), coords)
  q11_min_q22_A1 = fit_xy2_quad_A1(extF%nterms(2), extF%coef(1:extF%nterms(2),2), coords)
  q13_B2         = fit_xy2_quad_B2(extF%nterms(3), extF%coef(1:extF%nterms(3),3), coords)

  q11_A1 =  (q11_min_q22_A1-q33_A1)*0.5_ark
  q22_A1 = -(q11_min_q22_A1+q33_A1)*0.5_ark

  f(1:6) = (/ q11_A1, 0.0_ark, q13_B2, q22_A1, 0.0_ark, q33_A1 /)

end subroutine prop_xy2_qmom_bisect_frame


!###################################################################################################################################


! Subroutine to test xy2_qmom_bisect_frame if it is able to reproduce the original ab initio Cartesian components

subroutine TEST_xy2_qmom_bisect_frame(rank, ncoords, natoms, local, xyz, f)

  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)

  integer(ik) :: iounit, iatom, npoints, ipoint, info, i, ithread
  real(ark) :: xyz_(natoms,3), quad(6), quad0(6), diff(6), rms(6), angs
  character(cl) :: fname
  integer(ik), external :: omp_get_thread_num

  ithread = omp_get_thread_num()

  if (ithread==0) then

  write(out, '(/a)') 'TEST_xy2_qmom_bisect_frame/start: test polarizability tensor transformation'

  angs = 0.529177209_ark
  fname = 'QMOM.dat'

  call IOstart(fname, iounit)

  open(iounit,form='formatted',file=fname,iostat=info)
  if (info/=0) then
    write(out, '(/a,a,a)') 'TEST_xy2_qmom_bisect_frame error: file "', trim(fname), '" not found'
    stop
  endif

  ipoint = 0
  rms = 0
  do
    ipoint = ipoint + 1
    read(iounit,*,iostat=info) (xyz_(iatom,1:3), iatom=1, natoms), quad0(1:6)
    if (info/=0) exit
    !xyz_ = xyz_ * angs
    call prop_xy2_qmom_bisect_frame(rank, ncoords, natoms, local, xyz_, quad(1:6))
    diff = quad-quad0
    write(out, '(1x,i6,3x,6(4x,f12.4,1x,f12.4))') ipoint, (quad0(i), diff(i), i=1, 6)
    rms = rms + diff**2
  enddo
  npoints = ipoint-1
  close(iounit)

  rms = sqrt(rms/real(npoints,ark))

  write(out, '(/1x,a,6(1x,f12.4))') 'rms =', rms

  write(out, '(/a)') 'TEST_xy2_qmom_bisect_frame/done'

  endif

  !$omp barrier
  stop

end subroutine TEST_xy2_qmom_bisect_frame


!###################################################################################################################################


function fit_xy2_quad_A1(nparams, params, coords) result(f)

  integer(ik), intent(in) :: nparams
  real(ark), intent(in) :: params(nparams), coords(:)
  real(ark) :: f

  real(ark) :: y1, y2, y3, rho, alpha1, rad, f0, f1, f2, f3, req, alphaeq, beta

  rad = pi/180.0_ark

  req     = params(1)
  alphaeq = params(2)*rad
  beta    = params(3)

  y1     = (coords(1)-req) *exp(-beta*(coords(1)-req)**2)
  y2     = (coords(2)-req) *exp(-beta*(coords(2)-req)**2)
  rho    = alphaeq - coords(3)
  y3      = rho

  f0 = params(4)
  f1 = xy2_func_a1_n1_d6( (/y1,y2,y3/), params(5:16)  )  ! nparams = 12
  f2 = xy2_func_a1_n2_d6( (/y1,y2,y3/), params(17:40) )  ! nparams = 24
  f3 = xy2_func_a1_n3_d6( (/y1,y2,y3/), params(41:53) )  ! nparams = 13

  f = f0 + f1 + f2 + f3

end function fit_xy2_quad_A1


!###################################################################################################################################


function fit_xy2_quad_B2(nparams, params, coords) result(f)

  integer(ik), intent(in) :: nparams
  real(ark), intent(in) :: params(nparams), coords(:)
  real(ark) :: f

  real(ark) :: y1, y2, y3, rho, r1, r2, alpha1, rad, f0, f1, f2, f3, req, alphaeq, beta

  rad = pi/180.0_ark

  req     = params(1)
  alphaeq = params(2)*rad
  beta    = params(3)

  y1     = (coords(1)-req) *exp(-beta*(coords(1)-req)**2)
  y2     = (coords(2)-req) *exp(-beta*(coords(2)-req)**2)
  rho    = alphaeq - coords(3)
  y3      = rho

  f0 = params(4)
  f1 = xy2_func_b2_n1_d6((/y1,y2,y3/), params(5:10))    ! nparams = 6
  f2 = xy2_func_b2_n2_d6((/y1,y2,y3/), params(11:32))   ! nparams = 22
  f3 = xy2_func_b2_n3_d6((/y1,y2,y3/), params(33:44))   ! nparams = 12

  f = f0 + f1 + f2 + f3

end function fit_xy2_quad_B2


!###################################################################################################################################


function xy2_func_a1_n1_d6(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(12)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(12)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = r1 + r2
f(2) = r1**2 + r2**2
f(3) = r1**3 + r2**3
f(4) = r1**4 + r2**4
f(5) = r1**5 + r2**5
f(6) = r1**6 + r2**6
f(7) = a1
f(8) = a1**2
f(9) = a1**3
f(10) = a1**4
f(11) = a1**5
f(12) = a1**6
v = sum(f*params)
end function xy2_func_a1_n1_d6


!###################################################################################################################################


function xy2_func_a1_n2_d6(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(24)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(24)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = r1*r2
f(2) = r1*r2*(r1 + r2)
f(3) = r1*r2*(r1**2 + r2**2)
f(4) = r1*r2*(r1**3 + r2**3)
f(5) = r1*r2*(r1**4 + r2**4)
f(6) = r1**2*r2**2
f(7) = r1**2*r2**2*(r1 + r2)
f(8) = r1**2*r2**2*(r1**2 + r2**2)
f(9) = r1**3*r2**3
f(10) = a1*(r1 + r2)
f(11) = a1**2*(r1 + r2)
f(12) = a1**3*(r1 + r2)
f(13) = a1**4*(r1 + r2)
f(14) = a1**5*(r1 + r2)
f(15) = a1*(r1**2 + r2**2)
f(16) = a1**2*(r1**2 + r2**2)
f(17) = a1**3*(r1**2 + r2**2)
f(18) = a1**4*(r1**2 + r2**2)
f(19) = a1*(r1**3 + r2**3)
f(20) = a1**2*(r1**3 + r2**3)
f(21) = a1**3*(r1**3 + r2**3)
f(22) = a1*(r1**4 + r2**4)
f(23) = a1**2*(r1**4 + r2**4)
f(24) = a1*(r1**5 + r2**5)
v = sum(f*params)
end function xy2_func_a1_n2_d6


!###################################################################################################################################


function xy2_func_a1_n3_d6(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(13)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(13)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = a1*r1*r2
f(2) = a1**2*r1*r2
f(3) = a1**3*r1*r2
f(4) = a1**4*r1*r2
f(5) = a1*r1*r2*(r1 + r2)
f(6) = a1**2*r1*r2*(r1 + r2)
f(7) = a1**3*r1*r2*(r1 + r2)
f(8) = a1*r1*r2*(r1**2 + r2**2)
f(9) = a1**2*r1*r2*(r1**2 + r2**2)
f(10) = a1*r1*r2*(r1**3 + r2**3)
f(11) = a1*r1**2*r2**2
f(12) = a1**2*r1**2*r2**2
f(13) = a1*r1**2*r2**2*(r1 + r2)
v = sum(f*params)
end function xy2_func_a1_n3_d6


!###################################################################################################################################


function xy2_func_b2_n1_d6(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(6)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(6)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = -r1 + r2
f(2) = r1**2 - r2**2
f(3) = r1**3 - r2**3
f(4) = -r1**4 + r2**4
f(5) = r1**5 - r2**5
f(6) = -r1**6 + r2**6
v = sum(f*params)
end function xy2_func_b2_n1_d6


!###################################################################################################################################


function xy2_func_b2_n2_d6(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(22)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(22)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = r1*r2*(r1 - r2)
f(2) = r1*r2*(-r1**2 + r2**2)
f(3) = r1*r2*(r1**3 - r2**3)
f(4) = r1*r2*(r1**4 - r2**4)
f(5) = r1**2*r2**2*(r1 - r2)
f(6) = r1**2*r2**2*(r1**2 - r2**2)
f(7) = r1**2*r2**2*(-r1**2 + r2**2)
f(8) = a1*(-r1 + r2)
f(9) = a1**2*(-r1 + r2)
f(10) = a1**3*(r1 - r2)
f(11) = a1**4*(r1 - r2)
f(12) = a1**5*(r1 - r2)
f(13) = a1*(-r1**2 + r2**2)
f(14) = a1**2*(-r1**2 + r2**2)
f(15) = a1**3*(-r1**2 + r2**2)
f(16) = a1**4*(r1**2 - r2**2)
f(17) = a1*(-r1**3 + r2**3)
f(18) = a1**2*(r1**3 - r2**3)
f(19) = a1**3*(-r1**3 + r2**3)
f(20) = a1*(r1**4 - r2**4)
f(21) = a1**2*(-r1**4 + r2**4)
f(22) = a1*(r1**5 - r2**5)
v = sum(f*params)
end function xy2_func_b2_n2_d6


!###################################################################################################################################


function xy2_func_b2_n3_d6(coords, params) result(v)
real(ark), intent(in) :: coords(3)
real(ark), intent(in) :: params(12)
real(ark) :: v
real(ark) :: r1,r2,a1
real(ark) :: f(12)
r1 = coords(1)
r2 = coords(2)
a1 = coords(3)
f(1) = a1*r1*r2*(r1 - r2)
f(2) = a1**2*r1*r2*(r1 - r2)
f(3) = a1**3*r1*r2*(r1 - r2)
f(4) = a1*r1*r2*(r1**2 - r2**2)
f(5) = a1**2*r1*r2*(r1**2 - r2**2)
f(6) = a1*r1*r2*(-r1**3 + r2**3)
f(7) = a1*r1*r2*(-r1 + r2)
f(8) = a1**2*r1*r2*(-r1 + r2)
f(9) = a1**3*r1*r2*(-r1 + r2)
f(10) = a1*r1**2*r2**2*(-r1 + r2)
f(11) = a1*r1*r2*(-r1**2 + r2**2)
f(12) = a1**2*r1*r2*(-r1**2 + r2**2)
v = sum(f*params)
end function xy2_func_b2_n3_d6


end module prop_xy2_quad