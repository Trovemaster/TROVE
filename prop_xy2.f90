!
!  This unit defines all specific routines for a triatomic molecule of XY2 type
!
module prop_xy2
  use accuracy
  use moltype
  use timer
  use pot_xy2, only : MLloc2pqr_xy2

  implicit none

  public prop_xy2_qmom_sym, MLdipole_h2o_lpt2011, prop_xy2_sr

  private
 
  integer(ik), parameter :: verbose     = 4                          ! Verbosity level


  contains


recursive subroutine prop_xy2_sr(rank, ncoords, natoms, local, xyz, f)

  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)

  integer(ik) :: iatom, icentre
  real(ark) :: xyz0(3), xyz_(natoms,3), r1, r2, alpha1, e1(3), e2(3), e3(3), n1(3), n2(3), n3(3), tmat(3,3), &!
               coords(3), c_mb(3,3), c_xyz(3,3), c_xyz_(3,3), tmat_inv(3,3)

  icentre = nint(extf%coef(1,1)) ! icentre=1 means centre on Y1 and icentre=2 - on Y2

  if (rank/=9) then
    write(out, '(/a,1x,i3,1x,a)') 'prop_xy2_sr error: rank of the dipole moment vector =', rank, ', expected 9'
    stop
  endif

  xyz0 = xyz(1,:)
  do iatom=1, natoms
    xyz_(iatom,:) = xyz(iatom,:) - xyz0(:)
  enddo

  r1 = sqrt(sum(xyz_(2,:)**2))
  r2 = sqrt(sum(xyz_(3,:)**2))

  e1 = xyz_(2,:)/r1
  e2 = xyz_(3,:)/r2

  alpha1 = acos(sum(e1*e2))

  if (abs(alpha1-pi)<0.0001) stop 'prop_xy2_sr error: valence bond angle is 180 degrees (does not work for linear molecule)'

  coords = (/r1,r2,alpha1/)

  e3 = vector_product(e1,e2)

  tmat(1,:) = e1
  tmat(2,:) = e2
  tmat(3,:) = e3

  ! compute spin-rotation tensor in molecular-bond frame

  c_mb = 0
  c_mb(1,1) = fit_xy2_nosym(extF%nterms(1)-1, extF%coef(2:extF%nterms(1),1), coords) ! first parameter defines centre on Y1 or Y2, see above
  c_mb(2,2) = fit_xy2_nosym(extF%nterms(2), extF%coef(1:extF%nterms(2),2), coords)
  c_mb(3,3) = fit_xy2_nosym(extF%nterms(3), extF%coef(1:extF%nterms(3),3), coords)
  c_mb(1,2) = fit_xy2_nosym(extF%nterms(4), extF%coef(1:extF%nterms(4),4), coords)
  c_mb(2,1) = fit_xy2_nosym(extF%nterms(5), extF%coef(1:extF%nterms(5),5), coords)

  ! transform tensor into coordinate frame defined by the input xyz coordinates

  call matrix_pseudoinverse_ark( 3, 3, tmat(1:3,1:3), tmat_inv(1:3,1:3) )

  c_xyz_(1:3,1:3) = matmul(tmat_inv(1:3,1:3), c_mb(1:3,1:3))
  c_xyz(1:3,1:3) = matmul(c_xyz_(1:3,1:3), transpose(tmat_inv(1:3,1:3)))

  f(1:9) = (/ c_xyz(1,1), c_xyz(1,2), c_xyz(1,3), c_xyz(2,1), c_xyz(2,2), c_xyz(2,3), c_xyz(3,1), c_xyz(3,2), c_xyz(3,3) /)

  contains

  function vector_product(v1,v2) result (v)
    real(ark), intent(in) :: v1(3), v2(3)
    real(ark) :: v(3)
    v(1) = v1(2)*v2(3)-v1(3)*v2(2)
    v(2) = v1(3)*v2(1)-v1(1)*v2(3)
    v(3) = v1(1)*v2(2)-v1(2)*v2(1)
  end function vector_product

end subroutine prop_xy2_sr


function fit_xy2_nosym(nparams, params, coords) result(f)

  integer(ik), intent(in) :: nparams
  real(ark), intent(in) :: params(nparams), coords(:)
  real(ark) :: f

  real(ark) :: r1, r2, alpha1, rad, f0, f1, f2, f3, req, alphaeq, beta

  rad = real(pi,ark)/180.0_ark

  req     = params(1)
  alphaeq = params(2)*rad
  beta    = params(3)

  r1     = (coords(1)-req) * exp(-beta*(coords(1)-req)**2)
  r2     = (coords(2)-req) * exp(-beta*(coords(2)-req)**2)
  alpha1 = coords(3)-alphaeq

  f0 = params(4)
  f1 = xy2_func_n1_d6( (/r1,r2,alpha1/), params(5:22)  )  ! nparams = 18
  f2 = xy2_func_n2_d6( (/r1,r2,alpha1/), params(23:67) )  ! nparams = 45
  f3 = xy2_func_n3_d6( (/r1,r2,alpha1/), params(68:87) )  ! nparams = 20

  f = f0 + f1 + f2 + f3

end function fit_xy2_nosym


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


!###################################################################################################################################


recursive subroutine xy2_dipole_sym(rank, ncoords, natoms, local, xyz, f)

  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)

  integer(ik) :: iatom,ierr
  real(ark) :: xyz0(3), xyz_(natoms,3), r1, r2, alpha1, e1(3), e2(3), n1(3), n2(3), n3(3), tmat(3,3), coords(3), mu_mb(3), tmat_inv(3,3)

  if (rank/=3) then
    write(out, '(/a,1x,i3,1x,a)') 'xy2_dipole_sym error: rank of the dipole moment vector =', rank, ', expected 3'
    stop
  endif

  xyz0 = xyz(1,:)
  do iatom=1, natoms
    xyz_(iatom,:) = xyz(iatom,:) - xyz0(:)
  enddo

  r1 = sqrt(sum(xyz_(2,:)**2))
  r2 = sqrt(sum(xyz_(3,:)**2))

  e1 = xyz_(2,:)/r1
  e2 = xyz_(3,:)/r2

  alpha1 = acos(sum(e1*e2))

  n1 = e1 + e2
  n2 = e1 - e2

  n1 = n1/sqrt(sum(n1(:)**2))
  n2 = n2/sqrt(sum(n2(:)**2))

  n3 = vector_product(e1,e2)

  n3 = n3/sqrt(sum(n3(:)**2))

  tmat(1,:) = n1
  tmat(2,:) = n2
  tmat(3,:) = n3

  coords = (/r1,r2,alpha1/)

  mu_mb(1) = fit_xy2_dipole_A1(extF%nterms(1), extF%coef(1:extF%nterms(1),1), coords)
  mu_mb(2) = fit_xy2_dipole_B2(extF%nterms(2), extF%coef(1:extF%nterms(2),2), coords)
  mu_mb(3) = 0.0

  call matrix_pseudoinverse_ark( 3, 3, tmat(1:3,1:3), tmat_inv(1:3,1:3) )
  !
  f(1:3) = matmul(tmat_inv(1:3,1:3), mu_mb(1:3))
  !
  contains

  function vector_product(v1,v2) result (v)
    real(ark), intent(in) :: v1(3), v2(3)
    real(ark) :: v(3)
    v(1) = v1(2)*v2(3)-v1(3)*v2(2)
    v(2) = v1(3)*v2(1)-v1(1)*v2(3)
    v(3) = v1(1)*v2(2)-v1(2)*v2(1)
  end function vector_product

end subroutine xy2_dipole_sym


  subroutine MLinvmatark(al,ai,dimen,ierr)
  integer,intent(in)   :: dimen
  real(ark),intent(in)  :: al(:,:)
  real(ark),intent(out) :: ai(:,:)
  integer(ik),intent(out) :: ierr
  real(ark)             :: h(dimen),p,q
  integer(ik)           :: i1,i2,k,i,j,k8,k9
      

    ierr = 0
    ai=al
 
    do i1=1,dimen
      k=dimen-i1+1
      p=ai(1,1)
      do i=2,dimen
        q=ai(i,1)
        !
        if (abs(p)<small_a) then 
          !
          ierr = i
          !
          return
          !
        endif 
        !
        h(i)=q/p
        if(i.le.k) h(i)=-q/p
        do j=2,i
          k8=i-1
          k9=j-1
          ai(k8,k9)=ai(i,j)+q*h(j)
        enddo 
      enddo 
      ai(dimen,dimen)=1.0_ark/p
      do i=2,dimen
        k8=i-1
        ai(dimen,k8)=h(i)
      enddo 
   end do 
   do i=1,dimen
     k8=i-1
     do j=1,k8
       ai(j,i)=ai(i,j)
     enddo 
   enddo 
   return
 end subroutine MLinvmatark


! Subroutine to test xy2_dipole_sym if it is able to reproduce the original ab initio Cartesian components

subroutine TEST_xy2_dipole_sym(rank, ncoords, natoms, local, xyz, f)

  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)

  integer(ik) :: iounit, iatom, npoints, ipoint, info, ithread
  real(ark) :: xyz_(natoms,3), mu_xyz(3), mu_xyz0(3), rms(3), angs
  character(cl) :: fname
  integer(ik), external :: omp_get_thread_num

  f = 0
  ithread = omp_get_thread_num()
  print*, ithread

  if (ithread==0) then

  write(out, '(/a)') 'TEST_xy2_dipole_sym/start: test dipole moment transformation'

  angs = 0.529177209_ark
  fname = 'MU.dat'

  call IOstart(fname, iounit)

  open(iounit,form='formatted',file=fname,iostat=info)
  if (info/=0) then
    write(out, '(/a,a,a)') 'TEST_xy2_dipole_sym error: file "', trim(fname), '" not found'
    stop
  endif

  ipoint = 0
  rms = 0
  do
    ipoint = ipoint + 1
    read(iounit,*,iostat=info) (xyz_(iatom,1:3), iatom=1, natoms), mu_xyz0(1:3)
    if (info/=0) exit
    xyz_ = xyz_ * angs
    call xy2_dipole_sym(rank, ncoords, natoms, local, xyz_, mu_xyz(1:3))
    write(out, '(1x,i6,9(1x,f12.4))') ipoint, mu_xyz0(1:3), mu_xyz(1:3),  mu_xyz0(1:3)-mu_xyz(1:3)
    rms(1:3) = rms(1:3) + (mu_xyz0(1:3)-mu_xyz(1:3))**2
  enddo
  npoints = ipoint-1
  close(iounit)

  rms = sqrt(rms/real(npoints,ark))

  write(out, '(/1x,a,3(1x,f12.4))') 'rms =', rms

  write(out, '(/a)') 'TEST_xy2_dipole_sym/done'

  endif

  !$omp barrier
  stop

end subroutine TEST_xy2_dipole_sym


!###################################################################################################################################

! Electric field gradient tensor for XY2-type molecule, centered on Y atom

recursive subroutine xy2_efg_y(rank, ncoords, natoms, local, xyz, f)

  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)

  integer(ik) :: iatom, icentre
  real(ark) :: xyz0(3), xyz_(natoms,3), r1, r2, alpha1, e1(3), e2(3), e3(3), n1(3), n2(3), tmat(3,3), &!
               coords(3), efg_mb1(3,3), efg_mb2(3,3), efg_A1, efg_B2, efg_xyz(3,3,2), efg_xyz_(3,3), tmat_inv(3,3)

  icentre = nint(extf%coef(1,1)) ! icentre=1 means centre on Y1 and icentre=2 - on Y2

  if (rank/=6) then
    write(out, '(/a,1x,i3,1x,a)') 'xy2_efg_y error: rank of the dipole moment vector =', rank, ', expected 6'
    stop
  endif

  xyz0 = xyz(1,:)
  do iatom=1, natoms
    xyz_(iatom,:) = xyz(iatom,:) - xyz0(:)
  enddo

  r1 = sqrt(sum(xyz_(2,:)**2))
  r2 = sqrt(sum(xyz_(3,:)**2))

  e1 = xyz_(2,:)/r1
  e2 = xyz_(3,:)/r2

  alpha1 = acos(sum(e1*e2))

  if (abs(alpha1-pi)<0.0001) stop 'xy2_efg_y error: valence bond angle is 180 degrees (does not work for linear molecule)'

  coords = (/r1,r2,alpha1/)

  ! construct projection vectors

  if (icentre==1) then

    xyz0 = xyz(2,:)
    do iatom=1, natoms
      xyz_(iatom,:) = xyz(iatom,:) - xyz0(:)
    enddo

    n1 = xyz_(1,:)/sqrt(sum(xyz_(1,:)**2))
    n2 = xyz_(3,:)/sqrt(sum(xyz_(3,:)**2))

  elseif (icentre==2) then

    xyz0 = xyz(3,:)
    do iatom=1, natoms
      xyz_(iatom,:) = xyz(iatom,:) - xyz0(:)
    enddo

    n1 = xyz_(1,:)/sqrt(sum(xyz_(1,:)**2))
    n2 = xyz_(2,:)/sqrt(sum(xyz_(2,:)**2))

  endif

  e1 = n1 + n2
  e2 = n1 - n2

  e1 = e1/sqrt(sum(e1(:)**2))
  e2 = e2/sqrt(sum(e2(:)**2))

  e3 = vector_product(e1,e2)
  e3 = e3/sqrt(sum(e3(:)**2))

  tmat(1,:) = e1
  tmat(2,:) = e2
  tmat(3,:) = e3

  ! compute EFG tensor

  efg_mb1 = 0
  efg_mb2 = 0

  efg_A1 = fit_xy2_dipole_A1(extF%nterms(1)-1, extF%coef(2:extF%nterms(1),1), coords) ! first parameter defines centre on Y1 or Y2, see above
  efg_B2 = fit_xy2_dipole_B2(extF%nterms(2), extF%coef(1:extF%nterms(2),2), coords)

  efg_mb1(1,1) = (efg_A1 + efg_B2)*0.5_ark
  efg_mb2(1,1) = (efg_A1 - efg_B2)*0.5_ark

  efg_A1 = fit_xy2_dipole_A1(extF%nterms(3), extF%coef(1:extF%nterms(3),3), coords)
  efg_B2 = fit_xy2_dipole_B2(extF%nterms(4), extF%coef(1:extF%nterms(4),4), coords)

  efg_mb1(2,2) = (efg_A1 + efg_B2)*0.5_ark
  efg_mb2(2,2) = (efg_A1 - efg_B2)*0.5_ark

  efg_A1 = fit_xy2_dipole_A1(extF%nterms(5), extF%coef(1:extF%nterms(5),5), coords)
  efg_B2 = fit_xy2_dipole_B2(extF%nterms(6), extF%coef(1:extF%nterms(6),6), coords)

  efg_mb1(1,2) = (efg_A1 + efg_B2)*0.5_ark
  efg_mb2(1,2) = (efg_A1 - efg_B2)*0.5_ark

  efg_mb1(3,3) = -( efg_mb1(1,1) + efg_mb1(2,2) )
  efg_mb1(2,1) = efg_mb1(1,2)

  efg_mb2(3,3) = -( efg_mb2(1,1) + efg_mb2(2,2) )
  efg_mb2(2,1) = efg_mb2(1,2)

  ! transform EFG tensor into coordinate frame defined by the input xyz coordinates

  call matrix_pseudoinverse_ark( 3, 3, tmat(1:3,1:3), tmat_inv(1:3,1:3) )

  efg_xyz_(1:3,1:3) = matmul(tmat_inv(1:3,1:3), efg_mb1(1:3,1:3))
  efg_xyz(1:3,1:3,1) = matmul(efg_xyz_(1:3,1:3), transpose(tmat_inv(1:3,1:3)))

  efg_xyz_(1:3,1:3) = matmul(tmat_inv(1:3,1:3), efg_mb2(1:3,1:3))
  efg_xyz(1:3,1:3,2) = matmul(efg_xyz_(1:3,1:3), transpose(tmat_inv(1:3,1:3)))

  f(1:6) = (/ efg_xyz(1,1,icentre), efg_xyz(1,2,icentre), efg_xyz(1,3,icentre), &!
              efg_xyz(2,2,icentre), efg_xyz(2,3,icentre), efg_xyz(3,3,icentre) /)

  contains

  function vector_product(v1,v2) result (v)
    real(ark), intent(in) :: v1(3), v2(3)
    real(ark) :: v(3)
    v(1) = v1(2)*v2(3)-v1(3)*v2(2)
    v(2) = v1(3)*v2(1)-v1(1)*v2(3)
    v(3) = v1(1)*v2(2)-v1(2)*v2(1)
  end function vector_product

end subroutine xy2_efg_y


! Subroutine to test xy2_efg_y1 or xy2_efg_y2  if it is able to reproduce the original ab initio Cartesian components

subroutine TEST_xy2_efg_y(rank, ncoords, natoms, local, xyz, f)

  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)

  integer(ik) :: iounit, iatom, npoints, ipoint, info, i, j, k, l, ii, ithread
  real(ark) :: xyz_(natoms,3), mu_xyz(6), mu_xyz0(3,3), rms(6), angs, dmu(6)
  character(cl) :: fname
  integer(ik), external :: omp_get_thread_num

  f = 0
  ithread = omp_get_thread_num()

  if (ithread==0) then

  write(out, '(/a)') 'TEST_xy2_efg_y/start: test polarizability tensor transformation'

  angs = 0.529177209_ark
  fname = 'EFG.dat'

  call IOstart(fname, iounit)

  open(iounit,form='formatted',file=fname,iostat=info)
  if (info/=0) then
    write(out, '(/a,a,a)') 'xy2_efg_y error: file "', trim(fname), '" not found'
    stop
  endif

  ipoint = 0
  rms = 0
  do
    ipoint = ipoint + 1
    read(iounit,*,iostat=info) (xyz_(iatom,1:3), iatom=1, natoms), (mu_xyz0(i,1:3), i=1, 3)
    if (info/=0) exit
    xyz_ = xyz_ * angs
    call xy2_efg_y(rank, ncoords, natoms, local, xyz_, mu_xyz(1:6))
    ii = 0
    dmu = 0
    do i=1, 3
      do j=i, 3
        ii = ii + 1
        dmu(ii) = mu_xyz0(i,j) - mu_xyz(ii)
      enddo
    enddo
    write(out, '(1x,i6,3x,6(4x,f12.4,1x,f12.4))') ipoint, (dmu(i)+mu_xyz(i), dmu(i), i=1, 6)
    rms = rms + dmu**2
  enddo
  npoints = ipoint-1
  close(iounit)

  rms = sqrt(rms/real(npoints,ark))

  write(out, '(/1x,a,6(1x,f12.4))') 'rms =', rms

  write(out, '(/a)') 'TEST_xy2_efg_y/done'

  endif

  !$omp barrier
  stop

end subroutine TEST_xy2_efg_y


!###################################################################################################################################


recursive subroutine prop_xy2_qmom_sym(rank, ncoords, natoms, local, xyz, f)

  implicit none
  !
  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)

  integer(ik) :: iatom,ierr
  real(ark) :: xyz0(3), xyz_(natoms,3), r1, r2, alpha1, e1(3), e2(3), n1(3), n2(3), n3(3), tmat(3,3), &
               coords(3), qmom_mb(3,3), qmom_xyz(3,3), qmom_xyz_(3,3), tmat_inv(3,3), x(3,3)

  if (rank/=6) then
    write(out, '(/a,1x,i3,1x,a)') 'xy2_qmom_sym error: rank of the dipole moment vector =', rank, ', expected 6'
    stop
  endif


  ! xyz are undefined for the local case
  if (all(abs(xyz)<small_)) then
    !
    select case(trim(molec%coords_transform))
    case default
       write (out,"('prop_xy2_qmom_sym: coord. type ',a,' unknown')") trim(molec%coords_transform)
       stop 'prop_xy2_qmom_sym - bad coord. type'
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

  n1 = e1 + e2
  n2 = e1 - e2

  n1 = n1/sqrt(sum(n1(:)**2))
  n2 = n2/sqrt(sum(n2(:)**2))

  n3 = vector_product(e1,e2)

  n3 = n3/sqrt(sum(n3(:)**2))

  tmat(1,:) = n1
  tmat(2,:) = n2
  tmat(3,:) = n3

  coords = (/r1,r2,alpha1/)

  qmom_mb = 0
  qmom_mb(1,1) = fit_xy2_dipole_A1(extF%nterms(1), extF%coef(1:extF%nterms(1),1), coords)
  qmom_mb(2,2) = fit_xy2_dipole_A1(extF%nterms(2), extF%coef(1:extF%nterms(2),2), coords)
  qmom_mb(1,2) = fit_xy2_dipole_B2(extF%nterms(3), extF%coef(1:extF%nterms(3),3), coords)

  qmom_mb(3,3) = -( qmom_mb(1,1) + qmom_mb(2,2) )
  qmom_mb(2,1) = qmom_mb(1,2)

  call matrix_pseudoinverse_ark( 3, 3, tmat(1:3,1:3), tmat_inv(1:3,1:3) )
  !
  !call MLinvmatark(tmat,tmat_inv,3,ierr)
  !
  !if (ierr/=0) then
  !  write(out,"('xy2_dipole_sym error: failed inverse',i0)") ierr
  !  stop 'xy2_dipole_sym error: failed inverse'
  !endif
  !
  qmom_xyz_(1:3,1:3) = matmul(tmat_inv(1:3,1:3), qmom_mb(1:3,1:3))
  qmom_xyz(1:3,1:3) = matmul(qmom_xyz_(1:3,1:3), transpose(tmat_inv(1:3,1:3)))

  f(1:6) = (/ qmom_xyz(1,1), qmom_xyz(1,2), qmom_xyz(1,3), qmom_xyz(2,2), qmom_xyz(2,3), qmom_xyz(3,3) /)

  contains

  function vector_product(v1,v2) result (v)
    real(ark), intent(in) :: v1(3), v2(3)
    real(ark) :: v(3)
    v(1) = v1(2)*v2(3)-v1(3)*v2(2)
    v(2) = v1(3)*v2(1)-v1(1)*v2(3)
    v(3) = v1(1)*v2(2)-v1(2)*v2(1)
  end function vector_product

end subroutine prop_xy2_qmom_sym


! Subroutine to test xy2_qmom_sym if it is able to reproduce the original ab initio Cartesian components

subroutine TEST_xy2_qmom_sym(rank, ncoords, natoms, local, xyz, f)

  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)

  integer(ik) :: iounit, iatom, npoints, ipoint, info, i, j, k, l, ii, ithread
  real(ark) :: xyz_(natoms,3), mu_xyz(6), mu_xyz0(3,3), rms(6), angs, dmu(6)
  character(cl) :: fname
  integer(ik), external :: omp_get_thread_num

  ithread = omp_get_thread_num()

  if (ithread==0) then

  write(out, '(/a)') 'TEST_xy2_qmom_sym/start: test polarizability tensor transformation'

  angs = 0.529177209_ark
  fname = 'QMOM.dat'

  call IOstart(fname, iounit)

  open(iounit,form='formatted',file=fname,iostat=info)
  if (info/=0) then
    write(out, '(/a,a,a)') 'TEST_xy2_qmom_sym error: file "', trim(fname), '" not found'
    stop
  endif

  ipoint = 0
  rms = 0
  do
    ipoint = ipoint + 1
    read(iounit,*,iostat=info) (xyz_(iatom,1:3), iatom=1, natoms), (mu_xyz0(i,1:3), i=1, 3)
    if (info/=0) exit
    xyz_ = xyz_ * angs
    call prop_xy2_qmom_sym(rank, ncoords, natoms, local, xyz_, mu_xyz(1:6))
    ii = 0
    dmu = 0
    do i=1, 3
      do j=i, 3
        ii = ii + 1
        dmu(ii) = mu_xyz0(i,j) - mu_xyz(ii)
      enddo
    enddo
    write(out, '(1x,i6,3x,6(4x,f12.4,1x,f12.4))') ipoint, (dmu(i)+mu_xyz(i), dmu(i), i=1, 6)
    rms = rms + dmu**2
  enddo
  npoints = ipoint-1
  close(iounit)

  rms = sqrt(rms/real(npoints,ark))

  write(out, '(/1x,a,6(1x,f12.4))') 'rms =', rms

  write(out, '(/a)') 'TEST_xy2_qmom_sym/done'

  endif

  !$omp barrier
  stop

end subroutine TEST_xy2_qmom_sym


!###################################################################################################################################


recursive subroutine xy2_alpha_sym(rank, ncoords, natoms, local, xyz, f)

  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)

  integer(ik) :: iatom
  real(ark) :: xyz0(3), xyz_(natoms,3), r1, r2, alpha1, e1(3), e2(3), n1(3), n2(3), n3(3), tmat(3,3), &!
               coords(3), alpha_mb(3,3), alpha_xyz(3,3), alpha_xyz_(3,3), tmat_inv(3,3)

  if (rank/=6) then
    write(out, '(/a,1x,i3,1x,a)') 'xy2_alpha_sym error: rank of the dipole moment vector =', rank, ', expected 6'
    stop
  endif

  xyz0 = xyz(1,:)
  do iatom=1, natoms
    xyz_(iatom,:) = xyz(iatom,:) - xyz0(:)
  enddo

  r1 = sqrt(sum(xyz_(2,:)**2))
  r2 = sqrt(sum(xyz_(3,:)**2))

  e1 = xyz_(2,:)/r1
  e2 = xyz_(3,:)/r2

  alpha1 = acos(sum(e1*e2))

  n1 = e1 + e2
  n2 = e1 - e2

  n1 = n1/sqrt(sum(n1(:)**2))
  n2 = n2/sqrt(sum(n2(:)**2))

  n3 = vector_product(e1,e2)

  n3 = n3/sqrt(sum(n3(:)**2))

  tmat(1,:) = n1
  tmat(2,:) = n2
  tmat(3,:) = n3

  coords = (/r1,r2,alpha1/)

  alpha_mb = 0
  alpha_mb(1,1) = fit_xy2_dipole_A1(extF%nterms(1), extF%coef(1:extF%nterms(1),1), coords)
  alpha_mb(2,2) = fit_xy2_dipole_A1(extF%nterms(2), extF%coef(1:extF%nterms(2),2), coords)
  alpha_mb(3,3) = fit_xy2_dipole_A1(extF%nterms(3), extF%coef(1:extF%nterms(3),3), coords)
  alpha_mb(1,2) = fit_xy2_dipole_B2(extF%nterms(4), extF%coef(1:extF%nterms(4),4), coords)
  alpha_mb(2,1) = alpha_mb(1,2)

  call matrix_pseudoinverse_ark( 3, 3, tmat(1:3,1:3), tmat_inv(1:3,1:3) )

  alpha_xyz_(1:3,1:3) = matmul(tmat_inv(1:3,1:3), alpha_mb(1:3,1:3))
  alpha_xyz(1:3,1:3) = matmul(alpha_xyz_(1:3,1:3), transpose(tmat_inv(1:3,1:3)))

  f(1:6) = (/ alpha_xyz(1,1), alpha_xyz(1,2), alpha_xyz(1,3), alpha_xyz(2,2), alpha_xyz(2,3), alpha_xyz(3,3) /)

  contains

  function vector_product(v1,v2) result (v)
    real(ark), intent(in) :: v1(3), v2(3)
    real(ark) :: v(3)
    v(1) = v1(2)*v2(3)-v1(3)*v2(2)
    v(2) = v1(3)*v2(1)-v1(1)*v2(3)
    v(3) = v1(1)*v2(2)-v1(2)*v2(1)
  end function vector_product

end subroutine xy2_alpha_sym


! Subroutine to test xy2_alpha_sym if it is able to reproduce the original ab initio Cartesian components

subroutine TEST_xy2_alpha_sym(rank, ncoords, natoms, local, xyz, f)

  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)

  integer(ik) :: iounit, iatom, npoints, ipoint, info, i, j, k, l, ii, ithread
  real(ark) :: xyz_(natoms,3), mu_xyz(6), mu_xyz0(3,3), rms(6), angs, dmu(6)
  character(cl) :: fname
  integer(ik), external :: omp_get_thread_num

  ithread = omp_get_thread_num()

  if (ithread==0) then

  write(out, '(/a)') 'TEST_xy2_alpha_sym/start: test polarizability tensor transformation'

  angs = 0.529177209_ark
  fname = 'ALPHA.dat'

  call IOstart(fname, iounit)

  open(iounit,form='formatted',file=fname,iostat=info)
  if (info/=0) then
    write(out, '(/a,a,a)') 'TEST_xy2_alpha_sym error: file "', trim(fname), '" not found'
    stop
  endif

  ipoint = 0
  rms = 0
  do
    ipoint = ipoint + 1
    read(iounit,*,iostat=info) (xyz_(iatom,1:3), iatom=1, natoms), (mu_xyz0(i,1:3), i=1, 3)
    if (info/=0) exit
    xyz_ = xyz_ * angs
    call xy2_alpha_sym(rank, ncoords, natoms, local, xyz_, mu_xyz(1:6))
    ii = 0
    dmu = 0
    do i=1, 3
      do j=i, 3
        ii = ii + 1
        dmu(ii) = mu_xyz0(i,j) - mu_xyz(ii)
      enddo
    enddo
    write(out, '(1x,i6,3x,6(4x,f12.4,1x,f12.4))') ipoint, (dmu(i)+mu_xyz(i), dmu(i), i=1, 6)
    rms = rms + dmu**2
  enddo
  npoints = ipoint-1
  close(iounit)

  rms = sqrt(rms/real(npoints,ark))

  write(out, '(/1x,a,6(1x,f12.4))') 'rms =', rms

  write(out, '(/a)') 'TEST_xy2_alpha_sym/done'

  endif

  !$omp barrier
  stop

end subroutine TEST_xy2_alpha_sym


!###################################################################################################################################


function fit_xy2_dipole_A1(nparams, params, coords) result(f)

  integer(ik), intent(in) :: nparams
  real(ark), intent(in) :: params(nparams), coords(:)
  real(ark) :: f

  real(ark) :: r1, r2, alpha1, rad, f0, f1, f2, f3, req, alphaeq, beta

  rad = real(pi,ark)/180.0_ark

  req     = params(1)
  alphaeq = params(2)*rad
  beta    = params(3)

  r1     = (coords(1)-req) *exp(-beta*(coords(1)-req)**2)
  r2     = (coords(2)-req) *exp(-beta*(coords(2)-req)**2)
  alpha1 = coords(3)-alphaeq

  f0 = params(4)
  f1 = xy2_func_a1_n1_d6( (/r1,r2,alpha1/), params(5:16)  )  ! nparams = 12
  f2 = xy2_func_a1_n2_d6( (/r1,r2,alpha1/), params(17:40) )  ! nparams = 24
  f3 = xy2_func_a1_n3_d6( (/r1,r2,alpha1/), params(41:53) )  ! nparams = 13

  f = f0 + f1 + f2 + f3

end function fit_xy2_dipole_A1


!###################################################################################################################################


function fit_xy2_dipole_B2(nparams, params, coords) result(f)

  integer(ik), intent(in) :: nparams
  real(ark), intent(in) :: params(nparams), coords(:)
  real(ark) :: f

  real(ark) :: r1, r2, alpha1, rad, f0, f1, f2, f3, req, alphaeq, beta

  rad = real(pi,ark)/180.0_ark

  req     = params(1)
  alphaeq = params(2)*rad
  beta    = params(3)

  r1     = (coords(1)-req) *exp(-beta*(coords(1)-req)**2)
  r2     = (coords(2)-req) *exp(-beta*(coords(2)-req)**2)
  alpha1 = coords(3)-alphaeq

  f0 = params(4)
  f1 = xy2_func_b2_n1_d6((/r1,r2,alpha1/), params(5:10))    ! nparams = 6
  f2 = xy2_func_b2_n2_d6((/r1,r2,alpha1/), params(11:32))   ! nparams = 22
  f3 = xy2_func_b2_n3_d6((/r1,r2,alpha1/), params(33:44))   ! nparams = 12

  f = f0 + f1 + f2 + f3

end function fit_xy2_dipole_B2


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


subroutine matrix_pseudoinverse_ark(m, n, mat, invmat, info)

  integer(ik), intent(in) :: m, n
  real(ark), intent(in) :: mat(m,n)
  real(ark), intent(out) :: invmat(n,m)
  integer(ik), intent(out), optional :: info

  integer(ik) lwork, info_, i, j
  double precision work1(1), matd(m,n), matu(m,m), matvt(n,n), invmatd(n,m), mat_d(n,m), sv(n), tmat(n,m), tol
  double precision, allocatable :: work(:)

  tol = 1.0d-08 !epsilon(1.0d0)

  matd = dble(mat)

  lwork = -1
  call dgesvd('A', 'A', m, n, matd, m, sv, matu, m, matvt, n, work1, lwork, info_)
  lwork = int(work1(1), kind=ik)

  allocate(work(lwork), stat=info_)
  if (info_/=0) then
    write(out, '(/a,1x,i6)') 'matrix_pseudoinverse_ark error: failed to allocate workspace array required for SVD, size =', lwork
    stop
  endif

  call dgesvd('A', 'A', m, n, matd, m, sv, matu, m, matvt, n, work, lwork, info_)

  if (info_/=0) then
    write(out, '(/a,1x,i6)') 'matrix_pseudoinverse_ark error: SVD failed, info =', info_
    stop
  endif

  if (present(info)) info = 0

  mat_d = 0.0
  do i=1, n
    if (sv(i)>=tol) then
      mat_d(i,i) = 1.0d0/sv(i)
    else
      if (present(info)) then
        info = i
        mat_d(i,i) = 0.0
      else
        write(out, '(/a)') 'matrix_pseudoinverse_ark error: matrix is singular'
        write(out, '(a)') 'matrix:'
        do j=1, m
          write(out, '(<n>(1x,f10.6))') mat(j,1:n)
        enddo
        write(out, '(a)') 'singular elements:'
        do j=1, n
          write(out, '(1x,f)') sv(j)
        enddo
        stop 'STOP'
      endif
    endif
  enddo

  call dgemm('N', 'T', n, m, m, 1.0d0, mat_d(1:n,1:m), n, matu(1:m,1:m), m, 0.0d0, tmat(1:n,1:m), n)
  call dgemm('T', 'N', n, m, n, 1.0d0, matvt(1:n,1:n), n, tmat(1:n,1:m), n, 0.0d0, invmatd(1:n,1:m), n)

  invmat = real(invmatd,kind=ark)

  deallocate(work)

end subroutine matrix_pseudoinverse_ark


subroutine DIPS_LPT2011(muY, muX, r1, r2, cos_theta)
implicit none
!==============================================================================================================
! Lorenzo Lodi, 12 May 2011
! DIPS is Fortran 90 subroutine containing four
! dipole-moment surfaces (DMS) for the ground-state water molecule: LTP2011, LTP2011P, LTP2011NR, LTP2011S.
! Modify the character variable "dipole_surface" to select the DMS of choice (see below).
! REF: L. Lodi, J. Tennyson and O.L Polyansky, "A global, high accuracy ab initio dipole moment
! surface for the electronic ground state of the water molecule", J. Chem. Phys. (submitted) (2011).
!
! INPUT:  r1/bohrs, r2/bohrs, cos( theta )
! OUTPUT: mu_Y(a.u.), mu_X(a.u.)
! mu_Y is the "small" component, perpendicular to the bond-angle bisector  (=0 a.u. at equilibrium)
! mu_X is the "large" component, parallel      to the bond-angle bisector  (=~0.730 a.u. at equilibrium)
! See the ref. for a more precise description of the orientation of the axes.
!    M O R E     D E T A I L S
! All surfaces are based on IC-MRCI+Q[8,10] energy-derivative dipoles in the aug-cc-pCV6Z basis set,
! and incorporates a damping function to ensure a reasonable behaviour upon dissociation r->+inf.
! A part from LTP2001NR, they include relativistic corrections.
! The Y component of the dipole has an irregular behaviour for linear H-O-H geometries when one bond is
! stretched beyond about 3 bohrs (due to surface intersection).
!
! SPEED: on a Pentium4 3 GHz compiled with Intel Fortran compiler 11.1
! each call to DIPS takes ~2e-6 seconds. Could be made faster by implementing Horner polynomial evaluation.
! All the coefficients are hard-coded to make it self-contained. A DATA block would require
! too many continuation lines. As it is the routine takes about 40 kbytes of memory.
!
! NOTE: Numerically this fit is rather ill-conditioned, about 9 digits of accuracy are lost
!       in the summation because of floating-point truncation error.
!       Ugly as it is, it shouldn't be a problem in practice, though.
!==============================================================================================================
!integer , parameter :: dp=kind(1.d0)
real(ark), parameter :: PI=3.14159265358979323846264338328_ark
real(ark), intent(in) :: r1, r2, cos_theta
real(ark), intent(out) :: muX, muY
real(ark) :: theta
!******************************************************
logical :: zLoaded = .false.
integer :: n_r_X, n_theta_X
integer :: n_r_Y, n_theta_Y
integer, parameter :: n_max_coeffs=225       ! Statically allocate arrays for the fit coeffs.
integer, parameter :: max_exponent=9
real(ark) :: d_X(n_max_coeffs)                ! fit coefficients X
real(ark) :: d_Y(n_max_coeffs)                ! fit coefficients Y
integer  :: index_X(3,n_max_coeffs)          ! index table containing the exponents of the X fit
integer  :: index_Y(3,n_max_coeffs)          ! index table containing the exponents of the Y fit
real(ark) :: powers_x1(0:max_exponent), powers_x2(0:max_exponent), powers_x3(0:max_exponent)
integer :: ii, ncoeffs_X, ncoeffs_Y
real(ark) :: x1, x2, x3
real(ark) :: damp1, damp2, value1, value2, cos_theta_half, sin_theta_half
!real(ark) :: damp, muOH
!**********************************************************************************************
! Possible choices for dipole_surface are: LTP2011 (default), LTP2011P, LTP2011NR, LTP2001S
! All are based on aug-cc-pCV6Z IC-MRCI[8, 10] calculations
! The differ in the energies used to compute the dipoles (Davidson or Pople)
! in whether they include the MVD1 relativistic correction
! and in the number of fitting parameters included.
! LTP2011 should be used for studies involving very high energy (>25000 cm-1) states.
!
!  NAME              DIPOLES    RELATIVISTIC    FIT ORDER muX/ muY
!  LTP2011 (default)  Davidson   Yes             [8,8]/[9,8]
!  LTP2011P           Pople      Yes             [8,8]/[9,8]
!  LTP2011NR          Davidson   No              [8,8]/[9,8]
!  LTP2011S           Davidson   Yes             [6,7]/[7,7]
character(len=8), parameter :: dipole_surface = "LTP2011"
!**********************************************************************************************
!save n_r_X, n_theta_X, ncoeffs_X, index_X, d_X
!save n_r_Y, n_theta_Y, ncoeffs_Y, index_Y, d_Y
!save zLoaded


! Let's load the opportune coefficients for muX and muY
if(zLoaded .eqv. .false.) then
   d_X = 0._ark
   d_Y = 0._ark

 if(dipole_surface == "LTP2011") then
   n_r_X = 8; n_theta_X = 8; ncoeffs_X = 200
   index_X(1,001)= 0; index_X(2,001)= 0; index_X(3,001)= 1; d_X(001)= -4.3606625464226454E+01_ark
   index_X(1,002)= 0; index_X(2,002)= 0; index_X(3,002)= 2; d_X(002)=  3.3027196259894445E+03_ark
   index_X(1,003)= 0; index_X(2,003)= 0; index_X(3,003)= 3; d_X(003)= -1.2727332561105239E+04_ark
   index_X(1,004)= 0; index_X(2,004)= 0; index_X(3,004)= 4; d_X(004)=  2.0478327360204032E+04_ark
   index_X(1,005)= 0; index_X(2,005)= 0; index_X(3,005)= 5; d_X(005)= -1.7320391032929223E+04_ark
   index_X(1,006)= 0; index_X(2,006)= 0; index_X(3,006)= 6; d_X(006)=  8.0857296986251222E+03_ark
   index_X(1,007)= 0; index_X(2,007)= 0; index_X(3,007)= 7; d_X(007)= -1.9717556345513731E+03_ark
   index_X(1,008)= 0; index_X(2,008)= 0; index_X(3,008)= 8; d_X(008)=  1.9607508538273396E+02_ark
   index_X(1,009)= 0; index_X(2,009)= 2; index_X(3,009)= 1; d_X(009)=  3.0335075152822174E+02_ark
   index_X(1,010)= 0; index_X(2,010)= 2; index_X(3,010)= 2; d_X(010)= -1.1709356490854843E+03_ark
   index_X(1,011)= 0; index_X(2,011)= 2; index_X(3,011)= 3; d_X(011)=  1.9956778959599906E+03_ark
   index_X(1,012)= 0; index_X(2,012)= 2; index_X(3,012)= 4; d_X(012)= -1.9712161529225559E+03_ark
   index_X(1,013)= 0; index_X(2,013)= 2; index_X(3,013)= 5; d_X(013)=  1.2340220565548516E+03_ark
   index_X(1,014)= 0; index_X(2,014)= 2; index_X(3,014)= 6; d_X(014)= -4.8970182591694174E+02_ark
   index_X(1,015)= 0; index_X(2,015)= 2; index_X(3,015)= 7; d_X(015)=  1.1206177641027170E+02_ark
   index_X(1,016)= 0; index_X(2,016)= 2; index_X(3,016)= 8; d_X(016)= -1.1105175444354245E+01_ark
   index_X(1,017)= 0; index_X(2,017)= 4; index_X(3,017)= 1; d_X(017)=  1.6925709138547063E+01_ark
   index_X(1,018)= 0; index_X(2,018)= 4; index_X(3,018)= 2; d_X(018)= -4.3912662965647542E+01_ark
   index_X(1,019)= 0; index_X(2,019)= 4; index_X(3,019)= 3; d_X(019)=  2.4288671317884109E+01_ark
   index_X(1,020)= 0; index_X(2,020)= 4; index_X(3,020)= 4; d_X(020)=  3.4365336338240013E+01_ark
   index_X(1,021)= 0; index_X(2,021)= 4; index_X(3,021)= 5; d_X(021)= -5.5705027501924633E+01_ark
   index_X(1,022)= 0; index_X(2,022)= 4; index_X(3,022)= 6; d_X(022)=  3.1458201995923446E+01_ark
   index_X(1,023)= 0; index_X(2,023)= 4; index_X(3,023)= 7; d_X(023)= -8.0773946006338520E+00_ark
   index_X(1,024)= 0; index_X(2,024)= 4; index_X(3,024)= 8; d_X(024)=  7.8402541108334844E-01_ark
   index_X(1,025)= 0; index_X(2,025)= 6; index_X(3,025)= 1; d_X(025)=  1.3694677677176514E-01_ark
   index_X(1,026)= 0; index_X(2,026)= 6; index_X(3,026)= 2; d_X(026)= -1.2184309628182746E+00_ark
   index_X(1,027)= 0; index_X(2,027)= 6; index_X(3,027)= 3; d_X(027)=  3.9898597473516020E+00_ark
   index_X(1,028)= 0; index_X(2,028)= 6; index_X(3,028)= 4; d_X(028)= -6.4466157734350418E+00_ark
   index_X(1,029)= 0; index_X(2,029)= 6; index_X(3,029)= 5; d_X(029)=  5.6885223634108115E+00_ark
   index_X(1,030)= 0; index_X(2,030)= 6; index_X(3,030)= 6; d_X(030)= -2.7918382966599893E+00_ark
   index_X(1,031)= 0; index_X(2,031)= 6; index_X(3,031)= 7; d_X(031)=  7.1398484320161515E-01_ark
   index_X(1,032)= 0; index_X(2,032)= 6; index_X(3,032)= 8; d_X(032)= -7.3916051669584704E-02_ark
   index_X(1,033)= 0; index_X(2,033)= 8; index_X(3,033)= 1; d_X(033)=  2.6790328253412099E-04_ark
   index_X(1,034)= 0; index_X(2,034)= 8; index_X(3,034)= 2; d_X(034)= -4.0414235845958046E-03_ark
   index_X(1,035)= 0; index_X(2,035)= 8; index_X(3,035)= 3; d_X(035)=  1.5002012491777350E-02_ark
   index_X(1,036)= 0; index_X(2,036)= 8; index_X(3,036)= 4; d_X(036)= -2.5395644888476454E-02_ark
   index_X(1,037)= 0; index_X(2,037)= 8; index_X(3,037)= 5; d_X(037)=  2.3034273439407116E-02_ark
   index_X(1,038)= 0; index_X(2,038)= 8; index_X(3,038)= 6; d_X(038)= -1.1575981466648955E-02_ark
   index_X(1,039)= 0; index_X(2,039)= 8; index_X(3,039)= 7; d_X(039)=  3.0378956698768889E-03_ark
   index_X(1,040)= 0; index_X(2,040)= 8; index_X(3,040)= 8; d_X(040)= -3.2474704221385764E-04_ark
   index_X(1,041)= 1; index_X(2,041)= 0; index_X(3,041)= 1; d_X(041)=  1.7490126454435995E+02_ark
   index_X(1,042)= 1; index_X(2,042)= 0; index_X(3,042)= 2; d_X(042)= -6.6062814319734498E+03_ark
   index_X(1,043)= 1; index_X(2,043)= 0; index_X(3,043)= 3; d_X(043)=  2.4632070084248153E+04_ark
   index_X(1,044)= 1; index_X(2,044)= 0; index_X(3,044)= 4; d_X(044)= -3.9170928839035245E+04_ark
   index_X(1,045)= 1; index_X(2,045)= 0; index_X(3,045)= 5; d_X(045)=  3.2936121345909283E+04_ark
   index_X(1,046)= 1; index_X(2,046)= 0; index_X(3,046)= 6; d_X(046)= -1.5323433860934820E+04_ark
   index_X(1,047)= 1; index_X(2,047)= 0; index_X(3,047)= 7; d_X(047)=  3.7287418976416811E+03_ark
   index_X(1,048)= 1; index_X(2,048)= 0; index_X(3,048)= 8; d_X(048)= -3.7025900128080684E+02_ark
   index_X(1,049)= 1; index_X(2,049)= 2; index_X(3,049)= 1; d_X(049)= -3.8204879799564151E+02_ark
   index_X(1,050)= 1; index_X(2,050)= 2; index_X(3,050)= 2; d_X(050)=  1.4816809813710424E+03_ark
   index_X(1,051)= 1; index_X(2,051)= 2; index_X(3,051)= 3; d_X(051)= -2.5410134350081353E+03_ark
   index_X(1,052)= 1; index_X(2,052)= 2; index_X(3,052)= 4; d_X(052)=  2.5251484275723342E+03_ark
   index_X(1,053)= 1; index_X(2,053)= 2; index_X(3,053)= 5; d_X(053)= -1.5866169642140449E+03_ark
   index_X(1,054)= 1; index_X(2,054)= 2; index_X(3,054)= 6; d_X(054)=  6.2936423488600121E+02_ark
   index_X(1,055)= 1; index_X(2,055)= 2; index_X(3,055)= 7; d_X(055)= -1.4335551639070036E+02_ark
   index_X(1,056)= 1; index_X(2,056)= 2; index_X(3,056)= 8; d_X(056)=  1.4095056124002440E+01_ark
   index_X(1,057)= 1; index_X(2,057)= 4; index_X(3,057)= 1; d_X(057)= -1.3412212213906969E+01_ark
   index_X(1,058)= 1; index_X(2,058)= 4; index_X(3,058)= 2; d_X(058)=  3.4861250143381767E+01_ark
   index_X(1,059)= 1; index_X(2,059)= 4; index_X(3,059)= 3; d_X(059)= -1.9057659324804263E+01_ark
   index_X(1,060)= 1; index_X(2,060)= 4; index_X(3,060)= 4; d_X(060)= -2.8186887285377452E+01_ark
   index_X(1,061)= 1; index_X(2,061)= 4; index_X(3,061)= 5; d_X(061)=  4.5464201636816142E+01_ark
   index_X(1,062)= 1; index_X(2,062)= 4; index_X(3,062)= 6; d_X(062)= -2.5773491567353631E+01_ark
   index_X(1,063)= 1; index_X(2,063)= 4; index_X(3,063)= 7; d_X(063)=  6.6598208994000743E+00_ark
   index_X(1,064)= 1; index_X(2,064)= 4; index_X(3,064)= 8; d_X(064)= -6.5198890475585358E-01_ark
   index_X(1,065)= 1; index_X(2,065)= 6; index_X(3,065)= 1; d_X(065)= -3.8146332383576009E-02_ark
   index_X(1,066)= 1; index_X(2,066)= 6; index_X(3,066)= 2; d_X(066)=  3.4491055845546725E-01_ark
   index_X(1,067)= 1; index_X(2,067)= 6; index_X(3,067)= 3; d_X(067)= -1.1504040370477924E+00_ark
   index_X(1,068)= 1; index_X(2,068)= 6; index_X(3,068)= 4; d_X(068)=  1.8871106806586795E+00_ark
   index_X(1,069)= 1; index_X(2,069)= 6; index_X(3,069)= 5; d_X(069)= -1.6847004003047914E+00_ark
   index_X(1,070)= 1; index_X(2,070)= 6; index_X(3,070)= 6; d_X(070)=  8.3441442563798773E-01_ark
   index_X(1,071)= 1; index_X(2,071)= 6; index_X(3,071)= 7; d_X(071)= -2.1497937103777076E-01_ark
   index_X(1,072)= 1; index_X(2,072)= 6; index_X(3,072)= 8; d_X(072)=  2.2393755010853056E-02_ark
   index_X(1,073)= 2; index_X(2,073)= 0; index_X(3,073)= 1; d_X(073)= -2.1892849125475004E+02_ark
   index_X(1,074)= 2; index_X(2,074)= 0; index_X(3,074)= 2; d_X(074)=  5.6962800355550971E+03_ark
   index_X(1,075)= 2; index_X(2,075)= 0; index_X(3,075)= 3; d_X(075)= -2.0588561699539132E+04_ark
   index_X(1,076)= 2; index_X(2,076)= 0; index_X(3,076)= 4; d_X(076)=  3.2364759842469120E+04_ark
   index_X(1,077)= 2; index_X(2,077)= 0; index_X(3,077)= 5; d_X(077)= -2.7053799841186432E+04_ark
   index_X(1,078)= 2; index_X(2,078)= 0; index_X(3,078)= 6; d_X(078)=  1.2543623220686441E+04_ark
   index_X(1,079)= 2; index_X(2,079)= 0; index_X(3,079)= 7; d_X(079)= -3.0456876672294966E+03_ark
   index_X(1,080)= 2; index_X(2,080)= 0; index_X(3,080)= 8; d_X(080)=  3.0198138552824821E+02_ark
   index_X(1,081)= 2; index_X(2,081)= 2; index_X(3,081)= 1; d_X(081)=  1.9542895348778075E+02_ark
   index_X(1,082)= 2; index_X(2,082)= 2; index_X(3,082)= 2; d_X(082)= -7.6245485172340705E+02_ark
   index_X(1,083)= 2; index_X(2,083)= 2; index_X(3,083)= 3; d_X(083)=  1.3183009633221081E+03_ark
   index_X(1,084)= 2; index_X(2,084)= 2; index_X(3,084)= 4; d_X(084)= -1.3213562538592669E+03_ark
   index_X(1,085)= 2; index_X(2,085)= 2; index_X(3,085)= 5; d_X(085)=  8.3564896815556131E+02_ark
   index_X(1,086)= 2; index_X(2,086)= 2; index_X(3,086)= 6; d_X(086)= -3.3223581216292450E+02_ark
   index_X(1,087)= 2; index_X(2,087)= 2; index_X(3,087)= 7; d_X(087)=  7.5487293678030255E+01_ark
   index_X(1,088)= 2; index_X(2,088)= 2; index_X(3,088)= 8; d_X(088)= -7.3736087346769636E+00_ark
   index_X(1,089)= 2; index_X(2,089)= 4; index_X(3,089)= 1; d_X(089)=  3.9449696341125673E+00_ark
   index_X(1,090)= 2; index_X(2,090)= 4; index_X(3,090)= 2; d_X(090)= -1.0362686658220809E+01_ark
   index_X(1,091)= 2; index_X(2,091)= 4; index_X(3,091)= 3; d_X(091)=  5.8361834900792928E+00_ark
   index_X(1,092)= 2; index_X(2,092)= 4; index_X(3,092)= 4; d_X(092)=  8.1898774691727567E+00_ark
   index_X(1,093)= 2; index_X(2,093)= 4; index_X(3,093)= 5; d_X(093)= -1.3488763912672312E+01_ark
   index_X(1,094)= 2; index_X(2,094)= 4; index_X(3,094)= 6; d_X(094)=  7.7195530796107050E+00_ark
   index_X(1,095)= 2; index_X(2,095)= 4; index_X(3,095)= 7; d_X(095)= -2.0111792595416773E+00_ark
   index_X(1,096)= 2; index_X(2,096)= 4; index_X(3,096)= 8; d_X(096)=  1.9856014939796296E-01_ark
   index_X(1,097)= 2; index_X(2,097)= 6; index_X(3,097)= 1; d_X(097)=  2.6992892313586481E-03_ark
   index_X(1,098)= 2; index_X(2,098)= 6; index_X(3,098)= 2; d_X(098)= -2.3169406981168095E-02_ark
   index_X(1,099)= 2; index_X(2,099)= 6; index_X(3,099)= 3; d_X(099)=  7.6920653129892269E-02_ark
   index_X(1,100)= 2; index_X(2,100)= 6; index_X(3,100)= 4; d_X(100)= -1.2699177163437980E-01_ark
   index_X(1,101)= 2; index_X(2,101)= 6; index_X(3,101)= 5; d_X(101)=  1.1419975744036037E-01_ark
   index_X(1,102)= 2; index_X(2,102)= 6; index_X(3,102)= 6; d_X(102)= -5.6911577327582563E-02_ark
   index_X(1,103)= 2; index_X(2,103)= 6; index_X(3,103)= 7; d_X(103)=  1.4731428730556217E-02_ark
   index_X(1,104)= 2; index_X(2,104)= 6; index_X(3,104)= 8; d_X(104)= -1.5391529138923943E-03_ark
   index_X(1,105)= 3; index_X(2,105)= 0; index_X(3,105)= 1; d_X(105)=  1.3622336107641650E+02_ark
   index_X(1,106)= 3; index_X(2,106)= 0; index_X(3,106)= 2; d_X(106)= -2.7620215237450120E+03_ark
   index_X(1,107)= 3; index_X(2,107)= 0; index_X(3,107)= 3; d_X(107)=  9.6961566179726942E+03_ark
   index_X(1,108)= 3; index_X(2,108)= 0; index_X(3,108)= 4; d_X(108)= -1.5071634521499458E+04_ark
   index_X(1,109)= 3; index_X(2,109)= 0; index_X(3,109)= 5; d_X(109)=  1.2525535485496084E+04_ark
   index_X(1,110)= 3; index_X(2,110)= 0; index_X(3,110)= 6; d_X(110)= -5.7877475885727608E+03_ark
   index_X(1,111)= 3; index_X(2,111)= 0; index_X(3,111)= 7; d_X(111)=  1.4022418853302088E+03_ark
   index_X(1,112)= 3; index_X(2,112)= 0; index_X(3,112)= 8; d_X(112)= -1.3882010495019813E+02_ark
   index_X(1,113)= 3; index_X(2,113)= 2; index_X(3,113)= 1; d_X(113)= -5.1784446024105819E+01_ark
   index_X(1,114)= 3; index_X(2,114)= 2; index_X(3,114)= 2; d_X(114)=  2.0347186203791898E+02_ark
   index_X(1,115)= 3; index_X(2,115)= 2; index_X(3,115)= 3; d_X(115)= -3.5517945792744649E+02_ark
   index_X(1,116)= 3; index_X(2,116)= 2; index_X(3,116)= 4; d_X(116)=  3.5968370839488307E+02_ark
   index_X(1,117)= 3; index_X(2,117)= 2; index_X(3,117)= 5; d_X(117)= -2.2941598322942446E+02_ark
   index_X(1,118)= 3; index_X(2,118)= 2; index_X(3,118)= 6; d_X(118)=  9.1615227223264810E+01_ark
   index_X(1,119)= 3; index_X(2,119)= 2; index_X(3,119)= 7; d_X(119)= -2.0803746205627249E+01_ark
   index_X(1,120)= 3; index_X(2,120)= 2; index_X(3,120)= 8; d_X(120)=  2.0215951358823077E+00_ark
   index_X(1,121)= 3; index_X(2,121)= 4; index_X(3,121)= 1; d_X(121)= -5.1223059493923984E-01_ark
   index_X(1,122)= 3; index_X(2,122)= 4; index_X(3,122)= 2; d_X(122)=  1.3793274425197524E+00_ark
   index_X(1,123)= 3; index_X(2,123)= 4; index_X(3,123)= 3; d_X(123)= -8.6111511230274118E-01_ark
   index_X(1,124)= 3; index_X(2,124)= 4; index_X(3,124)= 4; d_X(124)= -9.3618549735765555E-01_ark
   index_X(1,125)= 3; index_X(2,125)= 4; index_X(3,125)= 5; d_X(125)=  1.6723633643267704E+00_ark
   index_X(1,126)= 3; index_X(2,126)= 4; index_X(3,126)= 6; d_X(126)= -9.7717597698039071E-01_ark
   index_X(1,127)= 3; index_X(2,127)= 4; index_X(3,127)= 7; d_X(127)=  2.5754294537479439E-01_ark
   index_X(1,128)= 3; index_X(2,128)= 4; index_X(3,128)= 8; d_X(128)= -2.5639662240507732E-02_ark
   index_X(1,129)= 4; index_X(2,129)= 0; index_X(3,129)= 1; d_X(129)= -4.8664524863777615E+01_ark
   index_X(1,130)= 4; index_X(2,130)= 0; index_X(3,130)= 2; d_X(130)=  8.2253513436822686E+02_ark
   index_X(1,131)= 4; index_X(2,131)= 0; index_X(3,131)= 3; d_X(131)= -2.8104242401845022E+03_ark
   index_X(1,132)= 4; index_X(2,132)= 0; index_X(3,132)= 4; d_X(132)=  4.3215344250088574E+03_ark
   index_X(1,133)= 4; index_X(2,133)= 0; index_X(3,133)= 5; d_X(133)= -3.5712787276469517E+03_ark
   index_X(1,134)= 4; index_X(2,134)= 0; index_X(3,134)= 6; d_X(134)=  1.6446990040199869E+03_ark
   index_X(1,135)= 4; index_X(2,135)= 0; index_X(3,135)= 7; d_X(135)= -3.9761250447597183E+02_ark
   index_X(1,136)= 4; index_X(2,136)= 0; index_X(3,136)= 8; d_X(136)=  3.9302009433029525E+01_ark
   index_X(1,137)= 4; index_X(2,137)= 2; index_X(3,137)= 1; d_X(137)=  7.4642192204291860E+00_ark
   index_X(1,138)= 4; index_X(2,138)= 2; index_X(3,138)= 2; d_X(138)= -2.9549045720626509E+01_ark
   index_X(1,139)= 4; index_X(2,139)= 2; index_X(3,139)= 3; d_X(139)=  5.2096102541524488E+01_ark
   index_X(1,140)= 4; index_X(2,140)= 2; index_X(3,140)= 4; d_X(140)= -5.3340409897036352E+01_ark
   index_X(1,141)= 4; index_X(2,141)= 2; index_X(3,141)= 5; d_X(141)=  3.4353054530474992E+01_ark
   index_X(1,142)= 4; index_X(2,142)= 2; index_X(3,142)= 6; d_X(142)= -1.3800578378807131E+01_ark
   index_X(1,143)= 4; index_X(2,143)= 2; index_X(3,143)= 7; d_X(143)=  3.1368956403484844E+00_ark
   index_X(1,144)= 4; index_X(2,144)= 2; index_X(3,144)= 8; d_X(144)= -3.0360358130678833E-01_ark
   index_X(1,145)= 4; index_X(2,145)= 4; index_X(3,145)= 1; d_X(145)=  2.4705941986514546E-02_ark
   index_X(1,146)= 4; index_X(2,146)= 4; index_X(3,146)= 2; d_X(146)= -6.9117806702358564E-02_ark
   index_X(1,147)= 4; index_X(2,147)= 4; index_X(3,147)= 3; d_X(147)=  5.0433248065949954E-02_ark
   index_X(1,148)= 4; index_X(2,148)= 4; index_X(3,148)= 4; d_X(148)=  3.2425384490856857E-02_ark
   index_X(1,149)= 4; index_X(2,149)= 4; index_X(3,149)= 5; d_X(149)= -7.1011085503201343E-02_ark
   index_X(1,150)= 4; index_X(2,150)= 4; index_X(3,150)= 6; d_X(150)=  4.3066331371839794E-02_ark
   index_X(1,151)= 4; index_X(2,151)= 4; index_X(3,151)= 7; d_X(151)= -1.1522000939407917E-02_ark
   index_X(1,152)= 4; index_X(2,152)= 4; index_X(3,152)= 8; d_X(152)=  1.1545216566446470E-03_ark
   index_X(1,153)= 5; index_X(2,153)= 0; index_X(3,153)= 1; d_X(153)=  1.0455117395260515E+01_ark
   index_X(1,154)= 5; index_X(2,154)= 0; index_X(3,154)= 2; d_X(154)= -1.5379526273858755E+02_ark
   index_X(1,155)= 5; index_X(2,155)= 0; index_X(3,155)= 3; d_X(155)=  5.1262413047375753E+02_ark
   index_X(1,156)= 5; index_X(2,156)= 0; index_X(3,156)= 4; d_X(156)= -7.8025562897946952E+02_ark
   index_X(1,157)= 5; index_X(2,157)= 0; index_X(3,157)= 5; d_X(157)=  6.4134095430272794E+02_ark
   index_X(1,158)= 5; index_X(2,158)= 0; index_X(3,158)= 6; d_X(158)= -2.9441601383245484E+02_ark
   index_X(1,159)= 5; index_X(2,159)= 0; index_X(3,159)= 7; d_X(159)=  7.1026892747306960E+01_ark
   index_X(1,160)= 5; index_X(2,160)= 0; index_X(3,160)= 8; d_X(160)= -7.0097630764862515E+00_ark
   index_X(1,161)= 5; index_X(2,161)= 2; index_X(3,161)= 1; d_X(161)= -5.5181072854595925E-01_ark
   index_X(1,162)= 5; index_X(2,162)= 2; index_X(3,162)= 2; d_X(162)=  2.1988689794610252E+00_ark
   index_X(1,163)= 5; index_X(2,163)= 2; index_X(3,163)= 3; d_X(163)= -3.9118084258312251E+00_ark
   index_X(1,164)= 5; index_X(2,164)= 2; index_X(3,164)= 4; d_X(164)=  4.0477572843179388E+00_ark
   index_X(1,165)= 5; index_X(2,165)= 2; index_X(3,165)= 5; d_X(165)= -2.6330324078872316E+00_ark
   index_X(1,166)= 5; index_X(2,166)= 2; index_X(3,166)= 6; d_X(166)=  1.0650938147472431E+00_ark
   index_X(1,167)= 5; index_X(2,167)= 2; index_X(3,167)= 7; d_X(167)= -2.4261849843153982E-01_ark
   index_X(1,168)= 5; index_X(2,168)= 2; index_X(3,168)= 8; d_X(168)=  2.3404996829526681E-02_ark
   index_X(1,169)= 6; index_X(2,169)= 0; index_X(3,169)= 1; d_X(169)= -1.3330355566063190E+00_ark
   index_X(1,170)= 6; index_X(2,170)= 0; index_X(3,170)= 2; d_X(170)=  1.7597899896022227E+01_ark
   index_X(1,171)= 6; index_X(2,171)= 0; index_X(3,171)= 3; d_X(171)= -5.7368093491015827E+01_ark
   index_X(1,172)= 6; index_X(2,172)= 0; index_X(3,172)= 4; d_X(172)=  8.6504254282326798E+01_ark
   index_X(1,173)= 6; index_X(2,173)= 0; index_X(3,173)= 5; d_X(173)= -7.0750176940476308E+01_ark
   index_X(1,174)= 6; index_X(2,174)= 0; index_X(3,174)= 6; d_X(174)=  3.2382040446805647E+01_ark
   index_X(1,175)= 6; index_X(2,175)= 0; index_X(3,175)= 7; d_X(175)= -7.7965304444907169E+00_ark
   index_X(1,176)= 6; index_X(2,176)= 0; index_X(3,176)= 8; d_X(176)=  7.6827845837534881E-01_ark
   index_X(1,177)= 6; index_X(2,177)= 2; index_X(3,177)= 1; d_X(177)=  1.6229968998813238E-02_ark
   index_X(1,178)= 6; index_X(2,178)= 2; index_X(3,178)= 2; d_X(178)= -6.4908083799454452E-02_ark
   index_X(1,179)= 6; index_X(2,179)= 2; index_X(3,179)= 3; d_X(179)=  1.1618507169545289E-01_ark
   index_X(1,180)= 6; index_X(2,180)= 2; index_X(3,180)= 4; d_X(180)= -1.2124387765841604E-01_ark
   index_X(1,181)= 6; index_X(2,181)= 2; index_X(3,181)= 5; d_X(181)=  7.9581487575241194E-02_ark
   index_X(1,182)= 6; index_X(2,182)= 2; index_X(3,182)= 6; d_X(182)= -3.2415275442599034E-02_ark
   index_X(1,183)= 6; index_X(2,183)= 2; index_X(3,183)= 7; d_X(183)=  7.4020781937269575E-03_ark
   index_X(1,184)= 6; index_X(2,184)= 2; index_X(3,184)= 8; d_X(184)= -7.1146892846662979E-04_ark
   index_X(1,185)= 7; index_X(2,185)= 0; index_X(3,185)= 1; d_X(185)=  9.2649525643705655E-02_ark
   index_X(1,186)= 7; index_X(2,186)= 0; index_X(3,186)= 2; d_X(186)= -1.1242874381538730E+00_ark
   index_X(1,187)= 7; index_X(2,187)= 0; index_X(3,187)= 3; d_X(187)=  3.5950660241542902E+00_ark
   index_X(1,188)= 7; index_X(2,188)= 0; index_X(3,188)= 4; d_X(188)= -5.3760392178762846E+00_ark
   index_X(1,189)= 7; index_X(2,189)= 0; index_X(3,189)= 5; d_X(189)=  4.3774577026751516E+00_ark
   index_X(1,190)= 7; index_X(2,190)= 0; index_X(3,190)= 6; d_X(190)= -1.9981707526528956E+00_ark
   index_X(1,191)= 7; index_X(2,191)= 0; index_X(3,191)= 7; d_X(191)=  4.8021225340536289E-01_ark
   index_X(1,192)= 7; index_X(2,192)= 0; index_X(3,192)= 8; d_X(192)= -4.7250502079100087E-02_ark
   index_X(1,193)= 8; index_X(2,193)= 0; index_X(3,193)= 1; d_X(193)= -2.6924375660369038E-03_ark
   index_X(1,194)= 8; index_X(2,194)= 0; index_X(3,194)= 2; d_X(194)=  3.0639958403565398E-02_ark
   index_X(1,195)= 8; index_X(2,195)= 0; index_X(3,195)= 3; d_X(195)= -9.6417377890139655E-02_ark
   index_X(1,196)= 8; index_X(2,196)= 0; index_X(3,196)= 4; d_X(196)=  1.4317113229849029E-01_ark
   index_X(1,197)= 8; index_X(2,197)= 0; index_X(3,197)= 5; d_X(197)= -1.1613665480567770E-01_ark
   index_X(1,198)= 8; index_X(2,198)= 0; index_X(3,198)= 6; d_X(198)=  5.2889764456063852E-02_ark
   index_X(1,199)= 8; index_X(2,199)= 0; index_X(3,199)= 7; d_X(199)= -1.2689809454080206E-02_ark
   index_X(1,200)= 8; index_X(2,200)= 0; index_X(3,200)= 8; d_X(200)=  1.2468169536632245E-03_ark
!********************************************************************************************
   n_r_Y = 9; n_theta_Y = 8; ncoeffs_Y = 225
   index_Y(1,001)= 0; index_Y(2,001)= 1; index_Y(3,001)= 0; d_Y(001)= -6.0025296063971437E+03_ark
   index_Y(1,002)= 0; index_Y(2,002)= 1; index_Y(3,002)= 1; d_Y(002)=  3.9291648647678376E+04_ark
   index_Y(1,003)= 0; index_Y(2,003)= 1; index_Y(3,003)= 2; d_Y(003)= -1.0569378589871156E+05_ark
   index_Y(1,004)= 0; index_Y(2,004)= 1; index_Y(3,004)= 3; d_Y(004)=  1.5020926617361890E+05_ark
   index_Y(1,005)= 0; index_Y(2,005)= 1; index_Y(3,005)= 4; d_Y(005)= -1.2289367399421075E+05_ark
   index_Y(1,006)= 0; index_Y(2,006)= 1; index_Y(3,006)= 5; d_Y(006)=  5.9210063162055972E+04_ark
   index_Y(1,007)= 0; index_Y(2,007)= 1; index_Y(3,007)= 6; d_Y(007)= -1.6309896166413906E+04_ark
   index_Y(1,008)= 0; index_Y(2,008)= 1; index_Y(3,008)= 7; d_Y(008)=  2.3058900193641894E+03_ark
   index_Y(1,009)= 0; index_Y(2,009)= 1; index_Y(3,009)= 8; d_Y(009)= -1.2187808235920966E+02_ark
   index_Y(1,010)= 0; index_Y(2,010)= 3; index_Y(3,010)= 0; d_Y(010)= -1.0438762294771004E+02_ark
   index_Y(1,011)= 0; index_Y(2,011)= 3; index_Y(3,011)= 1; d_Y(011)= -2.6954003638798604E+03_ark
   index_Y(1,012)= 0; index_Y(2,012)= 3; index_Y(3,012)= 2; d_Y(012)=  1.6112550047979807E+04_ark
   index_Y(1,013)= 0; index_Y(2,013)= 3; index_Y(3,013)= 3; d_Y(013)= -3.5922088374701911E+04_ark
   index_Y(1,014)= 0; index_Y(2,014)= 3; index_Y(3,014)= 4; d_Y(014)=  4.1234354044479434E+04_ark
   index_Y(1,015)= 0; index_Y(2,015)= 3; index_Y(3,015)= 5; d_Y(015)= -2.6668830950832620E+04_ark
   index_Y(1,016)= 0; index_Y(2,016)= 3; index_Y(3,016)= 6; d_Y(016)=  9.7672387842139869E+03_ark
   index_Y(1,017)= 0; index_Y(2,017)= 3; index_Y(3,017)= 7; d_Y(017)= -1.8779090769474860E+03_ark
   index_Y(1,018)= 0; index_Y(2,018)= 3; index_Y(3,018)= 8; d_Y(018)=  1.4585975086590042E+02_ark
   index_Y(1,019)= 0; index_Y(2,019)= 5; index_Y(3,019)= 0; d_Y(019)= -2.1822372815887320E+01_ark
   index_Y(1,020)= 0; index_Y(2,020)= 5; index_Y(3,020)= 1; d_Y(020)= -2.2110734126872558E+01_ark
   index_Y(1,021)= 0; index_Y(2,021)= 5; index_Y(3,021)= 2; d_Y(021)=  3.7571119498786572E+02_ark
   index_Y(1,022)= 0; index_Y(2,022)= 5; index_Y(3,022)= 3; d_Y(022)= -8.9280859359521492E+02_ark
   index_Y(1,023)= 0; index_Y(2,023)= 5; index_Y(3,023)= 4; d_Y(023)=  9.9854772253553529E+02_ark
   index_Y(1,024)= 0; index_Y(2,024)= 5; index_Y(3,024)= 5; d_Y(024)= -6.1509045913653972E+02_ark
   index_Y(1,025)= 0; index_Y(2,025)= 5; index_Y(3,025)= 6; d_Y(025)=  2.1314747796376469E+02_ark
   index_Y(1,026)= 0; index_Y(2,026)= 5; index_Y(3,026)= 7; d_Y(026)= -3.8767443468488636E+01_ark
   index_Y(1,027)= 0; index_Y(2,027)= 5; index_Y(3,027)= 8; d_Y(027)=  2.8598573219642276E+00_ark
   index_Y(1,028)= 0; index_Y(2,028)= 7; index_Y(3,028)= 0; d_Y(028)= -6.9625815262816104E-01_ark
   index_Y(1,029)= 0; index_Y(2,029)= 7; index_Y(3,029)= 1; d_Y(029)=  6.6073020690482736E+00_ark
   index_Y(1,030)= 0; index_Y(2,030)= 7; index_Y(3,030)= 2; d_Y(030)= -2.7882436566596880E+01_ark
   index_Y(1,031)= 0; index_Y(2,031)= 7; index_Y(3,031)= 3; d_Y(031)=  6.1709270655274850E+01_ark
   index_Y(1,032)= 0; index_Y(2,032)= 7; index_Y(3,032)= 4; d_Y(032)= -7.7645710806545594E+01_ark
   index_Y(1,033)= 0; index_Y(2,033)= 7; index_Y(3,033)= 5; d_Y(033)=  5.7561024705354612E+01_ark
   index_Y(1,034)= 0; index_Y(2,034)= 7; index_Y(3,034)= 6; d_Y(034)= -2.4882164185464717E+01_ark
   index_Y(1,035)= 0; index_Y(2,035)= 7; index_Y(3,035)= 7; d_Y(035)=  5.7990635062305955E+00_ark
   index_Y(1,036)= 0; index_Y(2,036)= 7; index_Y(3,036)= 8; d_Y(036)= -5.6292930951167364E-01_ark
   index_Y(1,037)= 0; index_Y(2,037)= 9; index_Y(3,037)= 0; d_Y(037)=  1.0287016772296820E-03_ark
   index_Y(1,038)= 0; index_Y(2,038)= 9; index_Y(3,038)= 1; d_Y(038)= -1.0860989693355805E-02_ark
   index_Y(1,039)= 0; index_Y(2,039)= 9; index_Y(3,039)= 2; d_Y(039)=  4.7514016035023587E-02_ark
   index_Y(1,040)= 0; index_Y(2,040)= 9; index_Y(3,040)= 3; d_Y(040)= -1.0624429525728374E-01_ark
   index_Y(1,041)= 0; index_Y(2,041)= 9; index_Y(3,041)= 4; d_Y(041)=  1.3119592284215287E-01_ark
   index_Y(1,042)= 0; index_Y(2,042)= 9; index_Y(3,042)= 5; d_Y(042)= -9.2928998177512767E-02_ark
   index_Y(1,043)= 0; index_Y(2,043)= 9; index_Y(3,043)= 6; d_Y(043)=  3.7497752367926296E-02_ark
   index_Y(1,044)= 0; index_Y(2,044)= 9; index_Y(3,044)= 7; d_Y(044)= -7.9895325261531980E-03_ark
   index_Y(1,045)= 0; index_Y(2,045)= 9; index_Y(3,045)= 8; d_Y(045)=  6.9510952312157315E-04_ark
   index_Y(1,046)= 1; index_Y(2,046)= 1; index_Y(3,046)= 0; d_Y(046)=  1.1107348110055093E+04_ark
   index_Y(1,047)= 1; index_Y(2,047)= 1; index_Y(3,047)= 1; d_Y(047)= -7.3812023665533910E+04_ark
   index_Y(1,048)= 1; index_Y(2,048)= 1; index_Y(3,048)= 2; d_Y(048)=  2.0230408994767399E+05_ark
   index_Y(1,049)= 1; index_Y(2,049)= 1; index_Y(3,049)= 3; d_Y(049)= -2.9410739202800079E+05_ark
   index_Y(1,050)= 1; index_Y(2,050)= 1; index_Y(3,050)= 4; d_Y(050)=  2.4753687850269134E+05_ark
   index_Y(1,051)= 1; index_Y(2,051)= 1; index_Y(3,051)= 5; d_Y(051)= -1.2375368727427488E+05_ark
   index_Y(1,052)= 1; index_Y(2,052)= 1; index_Y(3,052)= 6; d_Y(052)=  3.5882111890754779E+04_ark
   index_Y(1,053)= 1; index_Y(2,053)= 1; index_Y(3,053)= 7; d_Y(053)= -5.4833204488015035E+03_ark
   index_Y(1,054)= 1; index_Y(2,054)= 1; index_Y(3,054)= 8; d_Y(054)=  3.3225098766572773E+02_ark
   index_Y(1,055)= 1; index_Y(2,055)= 3; index_Y(3,055)= 0; d_Y(055)=  1.1210052264027763E+02_ark
   index_Y(1,056)= 1; index_Y(2,056)= 3; index_Y(3,056)= 1; d_Y(056)=  3.4399633045124938E+03_ark
   index_Y(1,057)= 1; index_Y(2,057)= 3; index_Y(3,057)= 2; d_Y(057)= -2.0218419913111138E+04_ark
   index_Y(1,058)= 1; index_Y(2,058)= 3; index_Y(3,058)= 3; d_Y(058)=  4.4955127263088827E+04_ark
   index_Y(1,059)= 1; index_Y(2,059)= 3; index_Y(3,059)= 4; d_Y(059)= -5.1673000975081872E+04_ark
   index_Y(1,060)= 1; index_Y(2,060)= 3; index_Y(3,060)= 5; d_Y(060)=  3.3546952728580436E+04_ark
   index_Y(1,061)= 1; index_Y(2,061)= 3; index_Y(3,061)= 6; d_Y(061)= -1.2360422032248898E+04_ark
   index_Y(1,062)= 1; index_Y(2,062)= 3; index_Y(3,062)= 7; d_Y(062)=  2.3971254975670017E+03_ark
   index_Y(1,063)= 1; index_Y(2,063)= 3; index_Y(3,063)= 8; d_Y(063)= -1.8849571401080175E+02_ark
   index_Y(1,064)= 1; index_Y(2,064)= 5; index_Y(3,064)= 0; d_Y(064)=  1.5801327026620129E+01_ark
   index_Y(1,065)= 1; index_Y(2,065)= 5; index_Y(3,065)= 1; d_Y(065)=  3.0778086309463106E+01_ark
   index_Y(1,066)= 1; index_Y(2,066)= 5; index_Y(3,066)= 2; d_Y(066)= -3.5038283671477348E+02_ark
   index_Y(1,067)= 1; index_Y(2,067)= 5; index_Y(3,067)= 3; d_Y(067)=  8.1883162195943078E+02_ark
   index_Y(1,068)= 1; index_Y(2,068)= 5; index_Y(3,068)= 4; d_Y(068)= -9.2646558530730545E+02_ark
   index_Y(1,069)= 1; index_Y(2,069)= 5; index_Y(3,069)= 5; d_Y(069)=  5.8468471543774649E+02_ark
   index_Y(1,070)= 1; index_Y(2,070)= 5; index_Y(3,070)= 6; d_Y(070)= -2.0984500327254500E+02_ark
   index_Y(1,071)= 1; index_Y(2,071)= 5; index_Y(3,071)= 7; d_Y(071)=  4.0027073412467871E+01_ark
   index_Y(1,072)= 1; index_Y(2,072)= 5; index_Y(3,072)= 8; d_Y(072)= -3.1483292751508998E+00_ark
   index_Y(1,073)= 1; index_Y(2,073)= 7; index_Y(3,073)= 0; d_Y(073)=  2.8459849010778271E-01_ark
   index_Y(1,074)= 1; index_Y(2,074)= 7; index_Y(3,074)= 1; d_Y(074)= -2.6711004988846980E+00_ark
   index_Y(1,075)= 1; index_Y(2,075)= 7; index_Y(3,075)= 2; d_Y(075)=  1.1089900367103837E+01_ark
   index_Y(1,076)= 1; index_Y(2,076)= 7; index_Y(3,076)= 3; d_Y(076)= -2.4258537179810446E+01_ark
   index_Y(1,077)= 1; index_Y(2,077)= 7; index_Y(3,077)= 4; d_Y(077)=  3.0265838116033592E+01_ark
   index_Y(1,078)= 1; index_Y(2,078)= 7; index_Y(3,078)= 5; d_Y(078)= -2.2283342391293445E+01_ark
   index_Y(1,079)= 1; index_Y(2,079)= 7; index_Y(3,079)= 6; d_Y(079)=  9.5735821221169317E+00_ark
   index_Y(1,080)= 1; index_Y(2,080)= 7; index_Y(3,080)= 7; d_Y(080)= -2.2183009214204503E+00_ark
   index_Y(1,081)= 1; index_Y(2,081)= 7; index_Y(3,081)= 8; d_Y(081)=  2.1411987383908127E-01_ark
   index_Y(1,082)= 2; index_Y(2,082)= 1; index_Y(3,082)= 0; d_Y(082)= -8.8954806908923638E+03_ark
   index_Y(1,083)= 2; index_Y(2,083)= 1; index_Y(3,083)= 1; d_Y(083)=  6.0029031672360856E+04_ark
   index_Y(1,084)= 2; index_Y(2,084)= 1; index_Y(3,084)= 2; d_Y(084)= -1.6752066997622303E+05_ark
   index_Y(1,085)= 2; index_Y(2,085)= 1; index_Y(3,085)= 3; d_Y(085)=  2.4866939134875033E+05_ark
   index_Y(1,086)= 2; index_Y(2,086)= 1; index_Y(3,086)= 4; d_Y(086)= -2.1454483436766942E+05_ark
   index_Y(1,087)= 2; index_Y(2,087)= 1; index_Y(3,087)= 5; d_Y(087)=  1.1058141888443311E+05_ark
   index_Y(1,088)= 2; index_Y(2,088)= 1; index_Y(3,088)= 6; d_Y(088)= -3.3344258355232014E+04_ark
   index_Y(1,089)= 2; index_Y(2,089)= 1; index_Y(3,089)= 7; d_Y(089)=  5.3747705450634530E+03_ark
   index_Y(1,090)= 2; index_Y(2,090)= 1; index_Y(3,090)= 8; d_Y(090)= -3.5245815839042189E+02_ark
   index_Y(1,091)= 2; index_Y(2,091)= 3; index_Y(3,091)= 0; d_Y(091)= -4.7782641193291056E+01_ark
   index_Y(1,092)= 2; index_Y(2,092)= 3; index_Y(3,092)= 1; d_Y(092)= -1.7770035056206107E+03_ark
   index_Y(1,093)= 2; index_Y(2,093)= 3; index_Y(3,093)= 2; d_Y(093)=  1.0278888960119511E+04_ark
   index_Y(1,094)= 2; index_Y(2,094)= 3; index_Y(3,094)= 3; d_Y(094)= -2.2787458652314963E+04_ark
   index_Y(1,095)= 2; index_Y(2,095)= 3; index_Y(3,095)= 4; d_Y(095)=  2.6208096604709513E+04_ark
   index_Y(1,096)= 2; index_Y(2,096)= 3; index_Y(3,096)= 5; d_Y(096)= -1.7057952242937361E+04_ark
   index_Y(1,097)= 2; index_Y(2,097)= 3; index_Y(3,097)= 6; d_Y(097)=  6.3111728824057036E+03_ark
   index_Y(1,098)= 2; index_Y(2,098)= 3; index_Y(3,098)= 7; d_Y(098)= -1.2312067912722996E+03_ark
   index_Y(1,099)= 2; index_Y(2,099)= 3; index_Y(3,099)= 8; d_Y(099)=  9.7611219354381319E+01_ark
   index_Y(1,100)= 2; index_Y(2,100)= 5; index_Y(3,100)= 0; d_Y(100)= -3.9142006385295645E+00_ark
   index_Y(1,101)= 2; index_Y(2,101)= 5; index_Y(3,101)= 1; d_Y(101)= -1.5359912640847142E+01_ark
   index_Y(1,102)= 2; index_Y(2,102)= 5; index_Y(3,102)= 2; d_Y(102)=  1.2776397874860595E+02_ark
   index_Y(1,103)= 2; index_Y(2,103)= 5; index_Y(3,103)= 3; d_Y(103)= -2.9297032442760428E+02_ark
   index_Y(1,104)= 2; index_Y(2,104)= 5; index_Y(3,104)= 4; d_Y(104)=  3.3596087026596683E+02_ark
   index_Y(1,105)= 2; index_Y(2,105)= 5; index_Y(3,105)= 5; d_Y(105)= -2.1782781495290965E+02_ark
   index_Y(1,106)= 2; index_Y(2,106)= 5; index_Y(3,106)= 6; d_Y(106)=  8.1111138066315789E+01_ark
   index_Y(1,107)= 2; index_Y(2,107)= 5; index_Y(3,107)= 7; d_Y(107)= -1.6200769375291998E+01_ark
   index_Y(1,108)= 2; index_Y(2,108)= 5; index_Y(3,108)= 8; d_Y(108)=  1.3474396226224599E+00_ark
   index_Y(1,109)= 2; index_Y(2,109)= 7; index_Y(3,109)= 0; d_Y(109)= -2.8478794804414065E-02_ark
   index_Y(1,110)= 2; index_Y(2,110)= 7; index_Y(3,110)= 1; d_Y(110)=  2.6593452089218772E-01_ark
   index_Y(1,111)= 2; index_Y(2,111)= 7; index_Y(3,111)= 2; d_Y(111)= -1.0932653638288912E+00_ark
   index_Y(1,112)= 2; index_Y(2,112)= 7; index_Y(3,112)= 3; d_Y(112)=  2.3736962936600321E+00_ark
   index_Y(1,113)= 2; index_Y(2,113)= 7; index_Y(3,113)= 4; d_Y(113)= -2.9437369227216550E+00_ark
   index_Y(1,114)= 2; index_Y(2,114)= 7; index_Y(3,114)= 5; d_Y(114)=  2.1552315808708045E+00_ark
   index_Y(1,115)= 2; index_Y(2,115)= 7; index_Y(3,115)= 6; d_Y(115)= -9.2075081022630911E-01_ark
   index_Y(1,116)= 2; index_Y(2,116)= 7; index_Y(3,116)= 7; d_Y(116)=  2.1211273777407769E-01_ark
   index_Y(1,117)= 2; index_Y(2,117)= 7; index_Y(3,117)= 8; d_Y(117)= -2.0351663807559817E-02_ark
   index_Y(1,118)= 3; index_Y(2,118)= 1; index_Y(3,118)= 0; d_Y(118)=  4.0244729236675139E+03_ark
   index_Y(1,119)= 3; index_Y(2,119)= 1; index_Y(3,119)= 1; d_Y(119)= -2.7579194260415024E+04_ark
   index_Y(1,120)= 3; index_Y(2,120)= 1; index_Y(3,120)= 2; d_Y(120)=  7.8286179604332443E+04_ark
   index_Y(1,121)= 3; index_Y(2,121)= 1; index_Y(3,121)= 3; d_Y(121)= -1.1841790641238638E+05_ark
   index_Y(1,122)= 3; index_Y(2,122)= 1; index_Y(3,122)= 4; d_Y(122)=  1.0438173383473718E+05_ark
   index_Y(1,123)= 3; index_Y(2,123)= 1; index_Y(3,123)= 5; d_Y(123)= -5.5170356933581454E+04_ark
   index_Y(1,124)= 3; index_Y(2,124)= 1; index_Y(3,124)= 6; d_Y(124)=  1.7149894029817311E+04_ark
   index_Y(1,125)= 3; index_Y(2,125)= 1; index_Y(3,125)= 7; d_Y(125)= -2.8724756908097334E+03_ark
   index_Y(1,126)= 3; index_Y(2,126)= 1; index_Y(3,126)= 8; d_Y(126)=  1.9822578233739478E+02_ark
   index_Y(1,127)= 3; index_Y(2,127)= 3; index_Y(3,127)= 0; d_Y(127)=  1.0548770133245853E+01_ark
   index_Y(1,128)= 3; index_Y(2,128)= 3; index_Y(3,128)= 1; d_Y(128)=  4.7160384631398483E+02_ark
   index_Y(1,129)= 3; index_Y(2,129)= 3; index_Y(3,129)= 2; d_Y(129)= -2.6927523238212161E+03_ark
   index_Y(1,130)= 3; index_Y(2,130)= 3; index_Y(3,130)= 3; d_Y(130)=  5.9529358189060440E+03_ark
   index_Y(1,131)= 3; index_Y(2,131)= 3; index_Y(3,131)= 4; d_Y(131)= -6.8453123481617804E+03_ark
   index_Y(1,132)= 3; index_Y(2,132)= 3; index_Y(3,132)= 5; d_Y(132)=  4.4601981106036874E+03_ark
   index_Y(1,133)= 3; index_Y(2,133)= 3; index_Y(3,133)= 6; d_Y(133)= -1.6534222457115247E+03_ark
   index_Y(1,134)= 3; index_Y(2,134)= 3; index_Y(3,134)= 7; d_Y(134)=  3.2343774421433045E+02_ark
   index_Y(1,135)= 3; index_Y(2,135)= 3; index_Y(3,135)= 8; d_Y(135)= -2.5735116318981454E+01_ark
   index_Y(1,136)= 3; index_Y(2,136)= 5; index_Y(3,136)= 0; d_Y(136)=  3.6760240764022001E-01_ark
   index_Y(1,137)= 3; index_Y(2,137)= 5; index_Y(3,137)= 1; d_Y(137)=  3.1648925950707962E+00_ark
   index_Y(1,138)= 3; index_Y(2,138)= 5; index_Y(3,138)= 2; d_Y(138)= -2.1074204175226441E+01_ark
   index_Y(1,139)= 3; index_Y(2,139)= 5; index_Y(3,139)= 3; d_Y(139)=  4.7397684800605475E+01_ark
   index_Y(1,140)= 3; index_Y(2,140)= 5; index_Y(3,140)= 4; d_Y(140)= -5.4950297427654732E+01_ark
   index_Y(1,141)= 3; index_Y(2,141)= 5; index_Y(3,141)= 5; d_Y(141)=  3.6442449369556925E+01_ark
   index_Y(1,142)= 3; index_Y(2,142)= 5; index_Y(3,142)= 6; d_Y(142)= -1.3973705057789175E+01_ark
   index_Y(1,143)= 3; index_Y(2,143)= 5; index_Y(3,143)= 7; d_Y(143)=  2.8879156570819760E+00_ark
   index_Y(1,144)= 3; index_Y(2,144)= 5; index_Y(3,144)= 8; d_Y(144)= -2.4945287573154928E-01_ark
   index_Y(1,145)= 4; index_Y(2,145)= 1; index_Y(3,145)= 0; d_Y(145)= -1.1241905260007698E+03_ark
   index_Y(1,146)= 4; index_Y(2,146)= 1; index_Y(3,146)= 1; d_Y(146)=  7.8203937269640610E+03_ark
   index_Y(1,147)= 4; index_Y(2,147)= 1; index_Y(3,147)= 2; d_Y(147)= -2.2550301375257903E+04_ark
   index_Y(1,148)= 4; index_Y(2,148)= 1; index_Y(3,148)= 3; d_Y(148)=  3.4682736691500708E+04_ark
   index_Y(1,149)= 4; index_Y(2,149)= 1; index_Y(3,149)= 4; d_Y(149)= -3.1134366985278044E+04_ark
   index_Y(1,150)= 4; index_Y(2,150)= 1; index_Y(3,150)= 5; d_Y(150)=  1.6797128474809051E+04_ark
   index_Y(1,151)= 4; index_Y(2,151)= 1; index_Y(3,151)= 6; d_Y(151)= -5.3466702394040858E+03_ark
   index_Y(1,152)= 4; index_Y(2,152)= 1; index_Y(3,152)= 7; d_Y(152)=  9.2113318815936827E+02_ark
   index_Y(1,153)= 4; index_Y(2,153)= 1; index_Y(3,153)= 8; d_Y(153)= -6.5819873444136647E+01_ark
   index_Y(1,154)= 4; index_Y(2,154)= 3; index_Y(3,154)= 0; d_Y(154)= -1.3629199079396130E+00_ark
   index_Y(1,155)= 4; index_Y(2,155)= 3; index_Y(3,155)= 1; d_Y(155)= -6.6848897861893420E+01_ark
   index_Y(1,156)= 4; index_Y(2,156)= 3; index_Y(3,156)= 2; d_Y(156)=  3.7931898119598372E+02_ark
   index_Y(1,157)= 4; index_Y(2,157)= 3; index_Y(3,157)= 3; d_Y(157)= -8.3698743594581538E+02_ark
   index_Y(1,158)= 4; index_Y(2,158)= 3; index_Y(3,158)= 4; d_Y(158)=  9.6148168185669510E+02_ark
   index_Y(1,159)= 4; index_Y(2,159)= 3; index_Y(3,159)= 5; d_Y(159)= -6.2592613445075858E+02_ark
   index_Y(1,160)= 4; index_Y(2,160)= 3; index_Y(3,160)= 6; d_Y(160)=  2.3178244914996230E+02_ark
   index_Y(1,161)= 4; index_Y(2,161)= 3; index_Y(3,161)= 7; d_Y(161)= -4.5266862800091985E+01_ark
   index_Y(1,162)= 4; index_Y(2,162)= 3; index_Y(3,162)= 8; d_Y(162)=  3.5923807281924383E+00_ark
   index_Y(1,163)= 4; index_Y(2,163)= 5; index_Y(3,163)= 0; d_Y(163)= -8.8188292516697686E-03_ark
   index_Y(1,164)= 4; index_Y(2,164)= 5; index_Y(3,164)= 1; d_Y(164)= -2.2668415176542567E-01_ark
   index_Y(1,165)= 4; index_Y(2,165)= 5; index_Y(3,165)= 2; d_Y(165)=  1.2977989413177795E+00_ark
   index_Y(1,166)= 4; index_Y(2,166)= 5; index_Y(3,166)= 3; d_Y(166)= -2.8690230287728440E+00_ark
   index_Y(1,167)= 4; index_Y(2,167)= 5; index_Y(3,167)= 4; d_Y(167)=  3.3521400971742992E+00_ark
   index_Y(1,168)= 4; index_Y(2,168)= 5; index_Y(3,168)= 5; d_Y(168)= -2.2606551635301457E+00_ark
   index_Y(1,169)= 4; index_Y(2,169)= 5; index_Y(3,169)= 6; d_Y(169)=  8.8515537020313673E-01_ark
   index_Y(1,170)= 4; index_Y(2,170)= 5; index_Y(3,170)= 7; d_Y(170)= -1.8718353979136282E-01_ark
   index_Y(1,171)= 4; index_Y(2,171)= 5; index_Y(3,171)= 8; d_Y(171)=  1.6558293980750705E-02_ark
   index_Y(1,172)= 5; index_Y(2,172)= 1; index_Y(3,172)= 0; d_Y(172)=  1.9837545634243168E+02_ark
   index_Y(1,173)= 5; index_Y(2,173)= 1; index_Y(3,173)= 1; d_Y(173)= -1.3997062144072215E+03_ark
   index_Y(1,174)= 5; index_Y(2,174)= 1; index_Y(3,174)= 2; d_Y(174)=  4.0931169349134307E+03_ark
   index_Y(1,175)= 5; index_Y(2,175)= 1; index_Y(3,175)= 3; d_Y(175)= -6.3858066957452929E+03_ark
   index_Y(1,176)= 5; index_Y(2,176)= 1; index_Y(3,176)= 4; d_Y(176)=  5.8197297979717314E+03_ark
   index_Y(1,177)= 5; index_Y(2,177)= 1; index_Y(3,177)= 5; d_Y(177)= -3.1917314434027830E+03_ark
   index_Y(1,178)= 5; index_Y(2,178)= 1; index_Y(3,178)= 6; d_Y(178)=  1.0346443125482040E+03_ark
   index_Y(1,179)= 5; index_Y(2,179)= 1; index_Y(3,179)= 7; d_Y(179)= -1.8198254073487243E+02_ark
   index_Y(1,180)= 5; index_Y(2,180)= 1; index_Y(3,180)= 8; d_Y(180)=  1.3322791574812072E+01_ark
   index_Y(1,181)= 5; index_Y(2,181)= 3; index_Y(3,181)= 0; d_Y(181)=  1.1165936308050561E-01_ark
   index_Y(1,182)= 5; index_Y(2,182)= 3; index_Y(3,182)= 1; d_Y(182)=  4.6716401640872505E+00_ark
   index_Y(1,183)= 5; index_Y(2,183)= 3; index_Y(3,183)= 2; d_Y(183)= -2.6730580795888727E+01_ark
   index_Y(1,184)= 5; index_Y(2,184)= 3; index_Y(3,184)= 3; d_Y(184)=  5.9012280279981155E+01_ark
   index_Y(1,185)= 5; index_Y(2,185)= 3; index_Y(3,185)= 4; d_Y(185)= -6.7650478121726849E+01_ark
   index_Y(1,186)= 5; index_Y(2,186)= 3; index_Y(3,186)= 5; d_Y(186)=  4.3866386197446985E+01_ark
   index_Y(1,187)= 5; index_Y(2,187)= 3; index_Y(3,187)= 6; d_Y(187)= -1.6145447026892441E+01_ark
   index_Y(1,188)= 5; index_Y(2,188)= 3; index_Y(3,188)= 7; d_Y(188)=  3.1251545236183915E+00_ark
   index_Y(1,189)= 5; index_Y(2,189)= 3; index_Y(3,189)= 8; d_Y(189)= -2.4475078288219265E-01_ark
   index_Y(1,190)= 6; index_Y(2,190)= 1; index_Y(3,190)= 0; d_Y(190)= -2.1572413034670078E+01_ark
   index_Y(1,191)= 6; index_Y(2,191)= 1; index_Y(3,191)= 1; d_Y(191)=  1.5418085563717500E+02_ark
   index_Y(1,192)= 6; index_Y(2,192)= 1; index_Y(3,192)= 2; d_Y(192)= -4.5630573530583376E+02_ark
   index_Y(1,193)= 6; index_Y(2,193)= 1; index_Y(3,193)= 3; d_Y(193)=  7.2028687205061055E+02_ark
   index_Y(1,194)= 6; index_Y(2,194)= 1; index_Y(3,194)= 4; d_Y(194)= -6.6435370200832926E+02_ark
   index_Y(1,195)= 6; index_Y(2,195)= 1; index_Y(3,195)= 5; d_Y(195)=  3.6897713056122154E+02_ark
   index_Y(1,196)= 6; index_Y(2,196)= 1; index_Y(3,196)= 6; d_Y(196)= -1.2123709316119232E+02_ark
   index_Y(1,197)= 6; index_Y(2,197)= 1; index_Y(3,197)= 7; d_Y(197)=  2.1641503394243188E+01_ark
   index_Y(1,198)= 6; index_Y(2,198)= 1; index_Y(3,198)= 8; d_Y(198)= -1.6107065097376747E+00_ark
   index_Y(1,199)= 6; index_Y(2,199)= 3; index_Y(3,199)= 0; d_Y(199)= -4.8973249873083269E-03_ark
   index_Y(1,200)= 6; index_Y(2,200)= 3; index_Y(3,200)= 1; d_Y(200)= -1.1852809969847300E-01_ark
   index_Y(1,201)= 6; index_Y(2,201)= 3; index_Y(3,201)= 2; d_Y(201)=  7.0847529039930990E-01_ark
   index_Y(1,202)= 6; index_Y(2,202)= 3; index_Y(3,202)= 3; d_Y(202)= -1.5740582223622930E+00_ark
   index_Y(1,203)= 6; index_Y(2,203)= 3; index_Y(3,203)= 4; d_Y(203)=  1.7975475049534229E+00_ark
   index_Y(1,204)= 6; index_Y(2,204)= 3; index_Y(3,204)= 5; d_Y(204)= -1.1537548422919350E+00_ark
   index_Y(1,205)= 6; index_Y(2,205)= 3; index_Y(3,205)= 6; d_Y(205)=  4.1766798107112368E-01_ark
   index_Y(1,206)= 6; index_Y(2,206)= 3; index_Y(3,206)= 7; d_Y(206)= -7.8846735007957136E-02_ark
   index_Y(1,207)= 6; index_Y(2,207)= 3; index_Y(3,207)= 8; d_Y(207)=  5.9433544921461134E-03_ark
   index_Y(1,208)= 7; index_Y(2,208)= 1; index_Y(3,208)= 0; d_Y(208)=  1.3200627540448489E+00_ark
   index_Y(1,209)= 7; index_Y(2,209)= 1; index_Y(3,209)= 1; d_Y(209)= -9.5383998224375830E+00_ark
   index_Y(1,210)= 7; index_Y(2,210)= 1; index_Y(3,210)= 2; d_Y(210)=  2.8500648410562967E+01_ark
   index_Y(1,211)= 7; index_Y(2,211)= 1; index_Y(3,211)= 3; d_Y(211)= -4.5392119414633939E+01_ark
   index_Y(1,212)= 7; index_Y(2,212)= 1; index_Y(3,212)= 4; d_Y(212)=  4.2237271808823294E+01_ark
   index_Y(1,213)= 7; index_Y(2,213)= 1; index_Y(3,213)= 5; d_Y(213)= -2.3668495437070739E+01_ark
   index_Y(1,214)= 7; index_Y(2,214)= 1; index_Y(3,214)= 6; d_Y(214)=  7.8486672857396815E+00_ark
   index_Y(1,215)= 7; index_Y(2,215)= 1; index_Y(3,215)= 7; d_Y(215)= -1.4145069394761105E+00_ark
   index_Y(1,216)= 7; index_Y(2,216)= 1; index_Y(3,216)= 8; d_Y(216)=  1.0634681558898884E-01_ark
   index_Y(1,217)= 8; index_Y(2,217)= 1; index_Y(3,217)= 0; d_Y(217)= -3.4744476473319408E-02_ark
   index_Y(1,218)= 8; index_Y(2,218)= 1; index_Y(3,218)= 1; d_Y(218)=  2.5316013291309014E-01_ark
   index_Y(1,219)= 8; index_Y(2,219)= 1; index_Y(3,219)= 2; d_Y(219)= -7.6152525136575100E-01_ark
   index_Y(1,220)= 8; index_Y(2,220)= 1; index_Y(3,220)= 3; d_Y(220)=  1.2199452241563369E+00_ark
   index_Y(1,221)= 8; index_Y(2,221)= 1; index_Y(3,221)= 4; d_Y(221)= -1.1412953127550940E+00_ark
   index_Y(1,222)= 8; index_Y(2,222)= 1; index_Y(3,222)= 5; d_Y(222)=  6.4282295895550057E-01_ark
   index_Y(1,223)= 8; index_Y(2,223)= 1; index_Y(3,223)= 6; d_Y(223)= -2.1419872523397476E-01_ark
   index_Y(1,224)= 8; index_Y(2,224)= 1; index_Y(3,224)= 7; d_Y(224)=  3.8777767747524089E-02_ark
   index_Y(1,225)= 8; index_Y(2,225)= 1; index_Y(3,225)= 8; d_Y(225)= -2.9272691210419748E-03_ark

 elseif(dipole_surface == "LTP2011P") then
   n_r_X = 8; n_theta_X = 8; ncoeffs_X = 200
   index_X(1,001)= 0; index_X(2,001)= 0; index_X(3,001)= 1; d_X(001)= -4.5369603348199234E+01_ark
   index_X(1,002)= 0; index_X(2,002)= 0; index_X(3,002)= 2; d_X(002)=  3.3121642614652737E+03_ark
   index_X(1,003)= 0; index_X(2,003)= 0; index_X(3,003)= 3; d_X(003)= -1.2713239499539919E+04_ark
   index_X(1,004)= 0; index_X(2,004)= 0; index_X(3,004)= 4; d_X(004)=  2.0401024583471099E+04_ark
   index_X(1,005)= 0; index_X(2,005)= 0; index_X(3,005)= 5; d_X(005)= -1.7219554492854211E+04_ark
   index_X(1,006)= 0; index_X(2,006)= 0; index_X(3,006)= 6; d_X(006)=  8.0253517338231704E+03_ark
   index_X(1,007)= 0; index_X(2,007)= 0; index_X(3,007)= 7; d_X(007)= -1.9543461683175301E+03_ark
   index_X(1,008)= 0; index_X(2,008)= 0; index_X(3,008)= 8; d_X(008)=  1.9411738062995573E+02_ark
   index_X(1,009)= 0; index_X(2,009)= 2; index_X(3,009)= 1; d_X(009)=  3.0210459000648279E+02_ark
   index_X(1,010)= 0; index_X(2,010)= 2; index_X(3,010)= 2; d_X(010)= -1.1654314922083213E+03_ark
   index_X(1,011)= 0; index_X(2,011)= 2; index_X(3,011)= 3; d_X(011)=  1.9876355844482823E+03_ark
   index_X(1,012)= 0; index_X(2,012)= 2; index_X(3,012)= 4; d_X(012)= -1.9630023267338838E+03_ark
   index_X(1,013)= 0; index_X(2,013)= 2; index_X(3,013)= 5; d_X(013)=  1.2261452365222358E+03_ark
   index_X(1,014)= 0; index_X(2,014)= 2; index_X(3,014)= 6; d_X(014)= -4.8433223359058320E+02_ark
   index_X(1,015)= 0; index_X(2,015)= 2; index_X(3,015)= 7; d_X(015)=  1.1014236423987313E+02_ark
   index_X(1,016)= 0; index_X(2,016)= 2; index_X(3,016)= 8; d_X(016)= -1.0840530090525135E+01_ark
   index_X(1,017)= 0; index_X(2,017)= 4; index_X(3,017)= 1; d_X(017)=  1.6872145325226484E+01_ark
   index_X(1,018)= 0; index_X(2,018)= 4; index_X(3,018)= 2; d_X(018)= -4.3245483564356618E+01_ark
   index_X(1,019)= 0; index_X(2,019)= 4; index_X(3,019)= 3; d_X(019)=  2.2442896296735853E+01_ark
   index_X(1,020)= 0; index_X(2,020)= 4; index_X(3,020)= 4; d_X(020)=  3.6711273235951921E+01_ark
   index_X(1,021)= 0; index_X(2,021)= 4; index_X(3,021)= 5; d_X(021)= -5.7283495960546134E+01_ark
   index_X(1,022)= 0; index_X(2,022)= 4; index_X(3,022)= 6; d_X(022)=  3.2033186406068126E+01_ark
   index_X(1,023)= 0; index_X(2,023)= 4; index_X(3,023)= 7; d_X(023)= -8.1826520846589119E+00_ark
   index_X(1,024)= 0; index_X(2,024)= 4; index_X(3,024)= 8; d_X(024)=  7.9138613227405585E-01_ark
   index_X(1,025)= 0; index_X(2,025)= 6; index_X(3,025)= 1; d_X(025)=  1.3343250226432701E-01_ark
   index_X(1,026)= 0; index_X(2,026)= 6; index_X(3,026)= 2; d_X(026)= -1.1984880199825056E+00_ark
   index_X(1,027)= 0; index_X(2,027)= 6; index_X(3,027)= 3; d_X(027)=  3.9411367604175211E+00_ark
   index_X(1,028)= 0; index_X(2,028)= 6; index_X(3,028)= 4; d_X(028)= -6.3736382917932133E+00_ark
   index_X(1,029)= 0; index_X(2,029)= 6; index_X(3,029)= 5; d_X(029)=  5.6213359918610877E+00_ark
   index_X(1,030)= 0; index_X(2,030)= 6; index_X(3,030)= 6; d_X(030)= -2.7560199326962902E+00_ark
   index_X(1,031)= 0; index_X(2,031)= 6; index_X(3,031)= 7; d_X(031)=  7.0397132626931125E-01_ark
   index_X(1,032)= 0; index_X(2,032)= 6; index_X(3,032)= 8; d_X(032)= -7.2786127157087321E-02_ark
   index_X(1,033)= 0; index_X(2,033)= 8; index_X(3,033)= 1; d_X(033)=  2.7092784705473605E-04_ark
   index_X(1,034)= 0; index_X(2,034)= 8; index_X(3,034)= 2; d_X(034)= -4.0608925539658003E-03_ark
   index_X(1,035)= 0; index_X(2,035)= 8; index_X(3,035)= 3; d_X(035)=  1.5108618930639750E-02_ark
   index_X(1,036)= 0; index_X(2,036)= 8; index_X(3,036)= 4; d_X(036)= -2.5656002762275421E-02_ark
   index_X(1,037)= 0; index_X(2,037)= 8; index_X(3,037)= 5; d_X(037)=  2.3340877331293086E-02_ark
   index_X(1,038)= 0; index_X(2,038)= 8; index_X(3,038)= 6; d_X(038)= -1.1759937940041709E-02_ark
   index_X(1,039)= 0; index_X(2,039)= 8; index_X(3,039)= 7; d_X(039)=  3.0921277498237032E-03_ark
   index_X(1,040)= 0; index_X(2,040)= 8; index_X(3,040)= 8; d_X(040)= -3.3098012590926373E-04_ark
   index_X(1,041)= 1; index_X(2,041)= 0; index_X(3,041)= 1; d_X(041)=  1.7738874086750820E+02_ark
   index_X(1,042)= 1; index_X(2,042)= 0; index_X(3,042)= 2; d_X(042)= -6.6194392040028943E+03_ark
   index_X(1,043)= 1; index_X(2,043)= 0; index_X(3,043)= 3; d_X(043)=  2.4595307914140554E+04_ark
   index_X(1,044)= 1; index_X(2,044)= 0; index_X(3,044)= 4; d_X(044)= -3.9013219897221083E+04_ark
   index_X(1,045)= 1; index_X(2,045)= 0; index_X(3,045)= 5; d_X(045)=  3.2737811094170491E+04_ark
   index_X(1,046)= 1; index_X(2,046)= 0; index_X(3,046)= 6; d_X(046)= -1.5206307423432383E+04_ark
   index_X(1,047)= 1; index_X(2,047)= 0; index_X(3,047)= 7; d_X(047)=  3.6951974893005245E+03_ark
   index_X(1,048)= 1; index_X(2,048)= 0; index_X(3,048)= 8; d_X(048)= -3.6650121700513409E+02_ark
   index_X(1,049)= 1; index_X(2,049)= 2; index_X(3,049)= 1; d_X(049)= -3.8042524644143487E+02_ark
   index_X(1,050)= 1; index_X(2,050)= 2; index_X(3,050)= 2; d_X(050)=  1.4748129876645835E+03_ark
   index_X(1,051)= 1; index_X(2,051)= 2; index_X(3,051)= 3; d_X(051)= -2.5311905380668904E+03_ark
   index_X(1,052)= 1; index_X(2,052)= 2; index_X(3,052)= 4; d_X(052)=  2.5148712676982395E+03_ark
   index_X(1,053)= 1; index_X(2,053)= 2; index_X(3,053)= 5; d_X(053)= -1.5763630152998667E+03_ark
   index_X(1,054)= 1; index_X(2,054)= 2; index_X(3,054)= 6; d_X(054)=  6.2228967720275978E+02_ark
   index_X(1,055)= 1; index_X(2,055)= 2; index_X(3,055)= 7; d_X(055)= -1.4083146205414960E+02_ark
   index_X(1,056)= 1; index_X(2,056)= 2; index_X(3,056)= 8; d_X(056)=  1.3748680534394225E+01_ark
   index_X(1,057)= 1; index_X(2,057)= 4; index_X(3,057)= 1; d_X(057)= -1.3373702520790175E+01_ark
   index_X(1,058)= 1; index_X(2,058)= 4; index_X(3,058)= 2; d_X(058)=  3.4368349304128060E+01_ark
   index_X(1,059)= 1; index_X(2,059)= 4; index_X(3,059)= 3; d_X(059)= -1.7705662862770623E+01_ark
   index_X(1,060)= 1; index_X(2,060)= 4; index_X(3,060)= 4; d_X(060)= -2.9869929862723438E+01_ark
   index_X(1,061)= 1; index_X(2,061)= 4; index_X(3,061)= 5; d_X(061)=  4.6560137615433632E+01_ark
   index_X(1,062)= 1; index_X(2,062)= 4; index_X(3,062)= 6; d_X(062)= -2.6153549821350680E+01_ark
   index_X(1,063)= 1; index_X(2,063)= 4; index_X(3,063)= 7; d_X(063)=  6.7241769890988508E+00_ark
   index_X(1,064)= 1; index_X(2,064)= 4; index_X(3,064)= 8; d_X(064)= -6.5589005081710638E-01_ark
   index_X(1,065)= 1; index_X(2,065)= 6; index_X(3,065)= 1; d_X(065)= -3.6751387100821375E-02_ark
   index_X(1,066)= 1; index_X(2,066)= 6; index_X(3,066)= 2; d_X(066)=  3.3626411488512531E-01_ark
   index_X(1,067)= 1; index_X(2,067)= 6; index_X(3,067)= 3; d_X(067)= -1.1275661887038950E+00_ark
   index_X(1,068)= 1; index_X(2,068)= 6; index_X(3,068)= 4; d_X(068)=  1.8521532810364079E+00_ark
   index_X(1,069)= 1; index_X(2,069)= 6; index_X(3,069)= 5; d_X(069)= -1.6530192149984941E+00_ark
   index_X(1,070)= 1; index_X(2,070)= 6; index_X(3,070)= 6; d_X(070)=  8.1799051694179070E-01_ark
   index_X(1,071)= 1; index_X(2,071)= 6; index_X(3,071)= 7; d_X(071)= -2.1051645222178195E-01_ark
   index_X(1,072)= 1; index_X(2,072)= 6; index_X(3,072)= 8; d_X(072)=  2.1902739881625166E-02_ark
   index_X(1,073)= 2; index_X(2,073)= 0; index_X(3,073)= 1; d_X(073)= -2.2026341186145237E+02_ark
   index_X(1,074)= 2; index_X(2,074)= 0; index_X(3,074)= 2; d_X(074)=  5.7031785205919332E+03_ark
   index_X(1,075)= 2; index_X(2,075)= 0; index_X(3,075)= 3; d_X(075)= -2.0550148655123492E+04_ark
   index_X(1,076)= 2; index_X(2,076)= 0; index_X(3,076)= 4; d_X(076)=  3.2226106761960680E+04_ark
   index_X(1,077)= 2; index_X(2,077)= 0; index_X(3,077)= 5; d_X(077)= -2.6885061494203863E+04_ark
   index_X(1,078)= 2; index_X(2,078)= 0; index_X(3,078)= 6; d_X(078)=  1.2445227945902654E+04_ark
   index_X(1,079)= 2; index_X(2,079)= 0; index_X(3,079)= 7; d_X(079)= -3.0176887122543376E+03_ark
   index_X(1,080)= 2; index_X(2,080)= 0; index_X(3,080)= 8; d_X(080)=  2.9885626785513887E+02_ark
   index_X(1,081)= 2; index_X(2,081)= 2; index_X(3,081)= 1; d_X(081)=  1.9456185151329555E+02_ark
   index_X(1,082)= 2; index_X(2,082)= 2; index_X(3,082)= 2; d_X(082)= -7.5894117952407396E+02_ark
   index_X(1,083)= 2; index_X(2,083)= 2; index_X(3,083)= 3; d_X(083)=  1.3133167779591895E+03_ark
   index_X(1,084)= 2; index_X(2,084)= 2; index_X(3,084)= 4; d_X(084)= -1.3158697661721817E+03_ark
   index_X(1,085)= 2; index_X(2,085)= 2; index_X(3,085)= 5; d_X(085)=  8.2992177491243638E+02_ark
   index_X(1,086)= 2; index_X(2,086)= 2; index_X(3,086)= 6; d_X(086)= -3.2826729721757829E+02_ark
   index_X(1,087)= 2; index_X(2,087)= 2; index_X(3,087)= 7; d_X(087)=  7.4084063585236436E+01_ark
   index_X(1,088)= 2; index_X(2,088)= 2; index_X(3,088)= 8; d_X(088)= -7.1829130178812193E+00_ark
   index_X(1,089)= 2; index_X(2,089)= 4; index_X(3,089)= 1; d_X(089)=  3.9355717177472798E+00_ark
   index_X(1,090)= 2; index_X(2,090)= 4; index_X(3,090)= 2; d_X(090)= -1.0235047437898174E+01_ark
   index_X(1,091)= 2; index_X(2,091)= 4; index_X(3,091)= 3; d_X(091)=  5.4922692793288661E+00_ark
   index_X(1,092)= 2; index_X(2,092)= 4; index_X(3,092)= 4; d_X(092)=  8.5995847397580292E+00_ark
   index_X(1,093)= 2; index_X(2,093)= 4; index_X(3,093)= 5; d_X(093)= -1.3736088628860784E+01_ark
   index_X(1,094)= 2; index_X(2,094)= 4; index_X(3,094)= 6; d_X(094)=  7.7946534546117618E+00_ark
   index_X(1,095)= 2; index_X(2,095)= 4; index_X(3,095)= 7; d_X(095)= -2.0207988752208621E+00_ark
   index_X(1,096)= 2; index_X(2,096)= 4; index_X(3,096)= 8; d_X(096)=  1.9875343013700331E-01_ark
   index_X(1,097)= 2; index_X(2,097)= 6; index_X(3,097)= 1; d_X(097)=  2.5639698780892672E-03_ark
   index_X(1,098)= 2; index_X(2,098)= 6; index_X(3,098)= 2; d_X(098)= -2.2285863535500994E-02_ark
   index_X(1,099)= 2; index_X(2,099)= 6; index_X(3,099)= 3; d_X(099)=  7.4458721011075646E-02_ark
   index_X(1,100)= 2; index_X(2,100)= 6; index_X(3,100)= 4; d_X(100)= -1.2312898757835455E-01_ark
   index_X(1,101)= 2; index_X(2,101)= 6; index_X(3,101)= 5; d_X(101)=  1.1068434039498243E-01_ark
   index_X(1,102)= 2; index_X(2,102)= 6; index_X(3,102)= 6; d_X(102)= -5.5098921513717869E-02_ark
   index_X(1,103)= 2; index_X(2,103)= 6; index_X(3,103)= 7; d_X(103)=  1.4243265089135093E-02_ark
   index_X(1,104)= 2; index_X(2,104)= 6; index_X(3,104)= 8; d_X(104)= -1.4859788393550843E-03_ark
   index_X(1,105)= 3; index_X(2,105)= 0; index_X(3,105)= 1; d_X(105)=  1.3651440213032902E+02_ark
   index_X(1,106)= 3; index_X(2,106)= 0; index_X(3,106)= 2; d_X(106)= -2.7634175766984063E+03_ark
   index_X(1,107)= 3; index_X(2,107)= 0; index_X(3,107)= 3; d_X(107)=  9.6745152774145899E+03_ark
   index_X(1,108)= 3; index_X(2,108)= 0; index_X(3,108)= 4; d_X(108)= -1.5002983886899117E+04_ark
   index_X(1,109)= 3; index_X(2,109)= 0; index_X(3,109)= 5; d_X(109)=  1.2444388982658649E+04_ark
   index_X(1,110)= 3; index_X(2,110)= 0; index_X(3,110)= 6; d_X(110)= -5.7409860982858218E+03_ark
   index_X(1,111)= 3; index_X(2,111)= 0; index_X(3,111)= 7; d_X(111)=  1.3890161652223824E+03_ark
   index_X(1,112)= 3; index_X(2,112)= 0; index_X(3,112)= 8; d_X(112)= -1.3734903421730905E+02_ark
   index_X(1,113)= 3; index_X(2,113)= 2; index_X(3,113)= 1; d_X(113)= -5.1541226375674341E+01_ark
   index_X(1,114)= 3; index_X(2,114)= 2; index_X(3,114)= 2; d_X(114)=  2.0252430116602773E+02_ark
   index_X(1,115)= 3; index_X(2,115)= 2; index_X(3,115)= 3; d_X(115)= -3.5381621035763601E+02_ark
   index_X(1,116)= 3; index_X(2,116)= 2; index_X(3,116)= 4; d_X(116)=  3.5805847628593801E+02_ark
   index_X(1,117)= 3; index_X(2,117)= 2; index_X(3,117)= 5; d_X(117)= -2.2764688095956990E+02_ark
   index_X(1,118)= 3; index_X(2,118)= 2; index_X(3,118)= 6; d_X(118)=  9.0397709662287525E+01_ark
   index_X(1,119)= 3; index_X(2,119)= 2; index_X(3,119)= 7; d_X(119)= -2.0380481928821609E+01_ark
   index_X(1,120)= 3; index_X(2,120)= 2; index_X(3,120)= 8; d_X(120)=  1.9649408918976405E+00_ark
   index_X(1,121)= 3; index_X(2,121)= 4; index_X(3,121)= 1; d_X(121)= -5.1138966153507681E-01_ark
   index_X(1,122)= 3; index_X(2,122)= 4; index_X(3,122)= 2; d_X(122)=  1.3661441966523853E+00_ark
   index_X(1,123)= 3; index_X(2,123)= 4; index_X(3,123)= 3; d_X(123)= -8.2702172930817142E-01_ark
   index_X(1,124)= 3; index_X(2,124)= 4; index_X(3,124)= 4; d_X(124)= -9.7265421016186338E-01_ark
   index_X(1,125)= 3; index_X(2,125)= 4; index_X(3,125)= 5; d_X(125)=  1.6898253479760115E+00_ark
   index_X(1,126)= 3; index_X(2,126)= 4; index_X(3,126)= 6; d_X(126)= -9.7978204775955646E-01_ark
   index_X(1,127)= 3; index_X(2,127)= 4; index_X(3,127)= 7; d_X(127)=  2.5697994352849207E-01_ark
   index_X(1,128)= 3; index_X(2,128)= 4; index_X(3,128)= 8; d_X(128)= -2.5478771171322023E-02_ark
   index_X(1,129)= 4; index_X(2,129)= 0; index_X(3,129)= 1; d_X(129)= -4.8654387322359185E+01_ark
   index_X(1,130)= 4; index_X(2,130)= 0; index_X(3,130)= 2; d_X(130)=  8.2242190879531154E+02_ark
   index_X(1,131)= 4; index_X(2,131)= 0; index_X(3,131)= 3; d_X(131)= -2.8031241274627664E+03_ark
   index_X(1,132)= 4; index_X(2,132)= 0; index_X(3,132)= 4; d_X(132)=  4.3005863272196921E+03_ark
   index_X(1,133)= 4; index_X(2,133)= 0; index_X(3,133)= 5; d_X(133)= -3.5471505975257060E+03_ark
   index_X(1,134)= 4; index_X(2,134)= 0; index_X(3,134)= 6; d_X(134)=  1.6309452788358606E+03_ark
   index_X(1,135)= 4; index_X(2,135)= 0; index_X(3,135)= 7; d_X(135)= -3.9374442698865153E+02_ark
   index_X(1,136)= 4; index_X(2,136)= 0; index_X(3,136)= 8; d_X(136)=  3.8873173765189676E+01_ark
   index_X(1,137)= 4; index_X(2,137)= 2; index_X(3,137)= 1; d_X(137)=  7.4263482937953995E+00_ark
   index_X(1,138)= 4; index_X(2,138)= 2; index_X(3,138)= 2; d_X(138)= -2.9405856334703230E+01_ark
   index_X(1,139)= 4; index_X(2,139)= 2; index_X(3,139)= 3; d_X(139)=  5.1880127762536176E+01_ark
   index_X(1,140)= 4; index_X(2,140)= 2; index_X(3,140)= 4; d_X(140)= -5.3055228290167406E+01_ark
   index_X(1,141)= 4; index_X(2,141)= 2; index_X(3,141)= 5; d_X(141)=  3.4033537692254413E+01_ark
   index_X(1,142)= 4; index_X(2,142)= 2; index_X(3,142)= 6; d_X(142)= -1.3584827501465725E+01_ark
   index_X(1,143)= 4; index_X(2,143)= 2; index_X(3,143)= 7; d_X(143)=  3.0637068668156644E+00_ark
   index_X(1,144)= 4; index_X(2,144)= 2; index_X(3,144)= 8; d_X(144)= -2.9400400090571566E-01_ark
   index_X(1,145)= 4; index_X(2,145)= 4; index_X(3,145)= 1; d_X(145)=  2.4690080975599216E-02_ark
   index_X(1,146)= 4; index_X(2,146)= 4; index_X(3,146)= 2; d_X(146)= -6.8697766478130973E-02_ark
   index_X(1,147)= 4; index_X(2,147)= 4; index_X(3,147)= 3; d_X(147)=  4.9467067030398937E-02_ark
   index_X(1,148)= 4; index_X(2,148)= 4; index_X(3,148)= 4; d_X(148)=  3.3100127227321252E-02_ark
   index_X(1,149)= 4; index_X(2,149)= 4; index_X(3,149)= 5; d_X(149)= -7.0895292664445719E-02_ark
   index_X(1,150)= 4; index_X(2,150)= 4; index_X(3,150)= 6; d_X(150)=  4.2720818901763380E-02_ark
   index_X(1,151)= 4; index_X(2,151)= 4; index_X(3,151)= 7; d_X(151)= -1.1374359110760679E-02_ark
   index_X(1,152)= 4; index_X(2,152)= 4; index_X(3,152)= 8; d_X(152)=  1.1342970628049898E-03_ark
   index_X(1,153)= 5; index_X(2,153)= 0; index_X(3,153)= 1; d_X(153)=  1.0436540988380131E+01_ark
   index_X(1,154)= 5; index_X(2,154)= 0; index_X(3,154)= 2; d_X(154)= -1.5368270153630783E+02_ark
   index_X(1,155)= 5; index_X(2,155)= 0; index_X(3,155)= 3; d_X(155)=  5.1110097429431448E+02_ark
   index_X(1,156)= 5; index_X(2,156)= 0; index_X(3,156)= 4; d_X(156)= -7.7621947769184590E+02_ark
   index_X(1,157)= 5; index_X(2,157)= 0; index_X(3,157)= 5; d_X(157)=  6.3679719376412140E+02_ark
   index_X(1,158)= 5; index_X(2,158)= 0; index_X(3,158)= 6; d_X(158)= -2.9185137246518241E+02_ark
   index_X(1,159)= 5; index_X(2,159)= 0; index_X(3,159)= 7; d_X(159)=  7.0309352600092097E+01_ark
   index_X(1,160)= 5; index_X(2,160)= 0; index_X(3,160)= 8; d_X(160)= -6.9304526297042912E+00_ark
   index_X(1,161)= 5; index_X(2,161)= 2; index_X(3,161)= 1; d_X(161)= -5.4869428913236362E-01_ark
   index_X(1,162)= 5; index_X(2,162)= 2; index_X(3,162)= 2; d_X(162)=  2.1872350770547087E+00_ark
   index_X(1,163)= 5; index_X(2,163)= 2; index_X(3,163)= 3; d_X(163)= -3.8926626315604551E+00_ark
   index_X(1,164)= 5; index_X(2,164)= 2; index_X(3,164)= 4; d_X(164)=  4.0196063233264567E+00_ark
   index_X(1,165)= 5; index_X(2,165)= 2; index_X(3,165)= 5; d_X(165)= -2.6011652493055912E+00_ark
   index_X(1,166)= 5; index_X(2,166)= 2; index_X(3,166)= 6; d_X(166)=  1.0441917028073675E+00_ark
   index_X(1,167)= 5; index_X(2,167)= 2; index_X(3,167)= 7; d_X(167)= -2.3573818048228645E-01_ark
   index_X(1,168)= 5; index_X(2,168)= 2; index_X(3,168)= 8; d_X(168)=  2.2524187110434468E-02_ark
   index_X(1,169)= 6; index_X(2,169)= 0; index_X(3,169)= 1; d_X(169)= -1.3291220057371689E+00_ark
   index_X(1,170)= 6; index_X(2,170)= 0; index_X(3,170)= 2; d_X(170)=  1.7575143292323617E+01_ark
   index_X(1,171)= 6; index_X(2,171)= 0; index_X(3,171)= 3; d_X(171)= -5.7175015072088598E+01_ark
   index_X(1,172)= 6; index_X(2,172)= 0; index_X(3,172)= 4; d_X(172)=  8.6024404207820368E+01_ark
   index_X(1,173)= 6; index_X(2,173)= 0; index_X(3,173)= 5; d_X(173)= -7.0220684042135588E+01_ark
   index_X(1,174)= 6; index_X(2,174)= 0; index_X(3,174)= 6; d_X(174)=  3.2085804235700714E+01_ark
   index_X(1,175)= 6; index_X(2,175)= 0; index_X(3,175)= 7; d_X(175)= -7.7140363956595692E+00_ark
   index_X(1,176)= 6; index_X(2,176)= 0; index_X(3,176)= 8; d_X(176)=  7.5918527350713561E-01_ark
   index_X(1,177)= 6; index_X(2,177)= 2; index_X(3,177)= 1; d_X(177)=  1.6123411124608911E-02_ark
   index_X(1,178)= 6; index_X(2,178)= 2; index_X(3,178)= 2; d_X(178)= -6.4505379064557111E-02_ark
   index_X(1,179)= 6; index_X(2,179)= 2; index_X(3,179)= 3; d_X(179)=  1.1543644995943403E-01_ark
   index_X(1,180)= 6; index_X(2,180)= 2; index_X(3,180)= 4; d_X(180)= -1.2003620164744744E-01_ark
   index_X(1,181)= 6; index_X(2,181)= 2; index_X(3,181)= 5; d_X(181)=  7.8223395979861721E-02_ark
   index_X(1,182)= 6; index_X(2,182)= 2; index_X(3,182)= 6; d_X(182)= -3.1554559699329410E-02_ark
   index_X(1,183)= 6; index_X(2,183)= 2; index_X(3,183)= 7; d_X(183)=  7.1278330034481030E-03_ark
   index_X(1,184)= 6; index_X(2,184)= 2; index_X(3,184)= 8; d_X(184)= -6.7727103487680296E-04_ark
   index_X(1,185)= 7; index_X(2,185)= 0; index_X(3,185)= 1; d_X(185)=  9.2289664396406845E-02_ark
   index_X(1,186)= 7; index_X(2,186)= 0; index_X(3,186)= 2; d_X(186)= -1.1222172572779332E+00_ark
   index_X(1,187)= 7; index_X(2,187)= 0; index_X(3,187)= 3; d_X(187)=  3.5814121358859019E+00_ark
   index_X(1,188)= 7; index_X(2,188)= 0; index_X(3,188)= 4; d_X(188)= -5.3438253503695741E+00_ark
   index_X(1,189)= 7; index_X(2,189)= 0; index_X(3,189)= 5; d_X(189)=  4.3425221592129066E+00_ark
   index_X(1,190)= 7; index_X(2,190)= 0; index_X(3,190)= 6; d_X(190)= -1.9787764705962503E+00_ark
   index_X(1,191)= 7; index_X(2,191)= 0; index_X(3,191)= 7; d_X(191)=  4.7483379302160922E-01_ark
   index_X(1,192)= 7; index_X(2,192)= 0; index_X(3,192)= 8; d_X(192)= -4.6659087756172335E-02_ark
   index_X(1,193)= 8; index_X(2,193)= 0; index_X(3,193)= 1; d_X(193)= -2.6796358497062911E-03_ark
   index_X(1,194)= 8; index_X(2,194)= 0; index_X(3,194)= 2; d_X(194)=  3.0566325542951939E-02_ark
   index_X(1,195)= 8; index_X(2,195)= 0; index_X(3,195)= 3; d_X(195)= -9.6003394596242131E-02_ark
   index_X(1,196)= 8; index_X(2,196)= 0; index_X(3,196)= 4; d_X(196)=  1.4223492896628337E-01_ark
   index_X(1,197)= 8; index_X(2,197)= 0; index_X(3,197)= 5; d_X(197)= -1.1513633676733945E-01_ark
   index_X(1,198)= 8; index_X(2,198)= 0; index_X(3,198)= 6; d_X(198)=  5.2338159228265617E-02_ark
   index_X(1,199)= 8; index_X(2,199)= 0; index_X(3,199)= 7; d_X(199)= -1.2537385453518853E-02_ark
   index_X(1,200)= 8; index_X(2,200)= 0; index_X(3,200)= 8; d_X(200)=  1.2300920389397851E-03_ark
!********************************************************************************************
   n_r_Y = 9; n_theta_Y = 8; ncoeffs_Y = 225
   index_Y(1,001)= 0; index_Y(2,001)= 1; index_Y(3,001)= 0; d_Y(001)= -5.8428391777016677E+03_ark
   index_Y(1,002)= 0; index_Y(2,002)= 1; index_Y(3,002)= 1; d_Y(002)=  3.2573988953241562E+04_ark
   index_Y(1,003)= 0; index_Y(2,003)= 1; index_Y(3,003)= 2; d_Y(003)= -6.8553974380785585E+04_ark
   index_Y(1,004)= 0; index_Y(2,004)= 1; index_Y(3,004)= 3; d_Y(004)=  6.4678362290164412E+04_ark
   index_Y(1,005)= 0; index_Y(2,005)= 1; index_Y(3,005)= 4; d_Y(005)= -1.8451749483852953E+04_ark
   index_Y(1,006)= 0; index_Y(2,006)= 1; index_Y(3,006)= 5; d_Y(006)= -1.4308801916067372E+04_ark
   index_Y(1,007)= 0; index_Y(2,007)= 1; index_Y(3,007)= 6; d_Y(007)=  1.3696210464942968E+04_ark
   index_Y(1,008)= 0; index_Y(2,008)= 1; index_Y(3,008)= 7; d_Y(008)= -4.2940724990385352E+03_ark
   index_Y(1,009)= 0; index_Y(2,009)= 1; index_Y(3,009)= 8; d_Y(009)=  4.8326208650320768E+02_ark
   index_Y(1,010)= 0; index_Y(2,010)= 3; index_Y(3,010)= 0; d_Y(010)= -1.1156119678427785E+02_ark
   index_Y(1,011)= 0; index_Y(2,011)= 3; index_Y(3,011)= 1; d_Y(011)= -2.3423585488410899E+03_ark
   index_Y(1,012)= 0; index_Y(2,012)= 3; index_Y(3,012)= 2; d_Y(012)=  1.5284863314126444E+04_ark
   index_Y(1,013)= 0; index_Y(2,013)= 3; index_Y(3,013)= 3; d_Y(013)= -3.6735022429284116E+04_ark
   index_Y(1,014)= 0; index_Y(2,014)= 3; index_Y(3,014)= 4; d_Y(014)=  4.5666097600721812E+04_ark
   index_Y(1,015)= 0; index_Y(2,015)= 3; index_Y(3,015)= 5; d_Y(015)= -3.2257718204151926E+04_ark
   index_Y(1,016)= 0; index_Y(2,016)= 3; index_Y(3,016)= 6; d_Y(016)=  1.3053587657162425E+04_ark
   index_Y(1,017)= 0; index_Y(2,017)= 3; index_Y(3,017)= 7; d_Y(017)= -2.8172483116416261E+03_ark
   index_Y(1,018)= 0; index_Y(2,018)= 3; index_Y(3,018)= 8; d_Y(018)=  2.5118648520440911E+02_ark
   index_Y(1,019)= 0; index_Y(2,019)= 5; index_Y(3,019)= 0; d_Y(019)= -2.3228852260314397E+01_ark
   index_Y(1,020)= 0; index_Y(2,020)= 5; index_Y(3,020)= 1; d_Y(020)= -7.8394696388659213E-01_ark
   index_Y(1,021)= 0; index_Y(2,021)= 5; index_Y(3,021)= 2; d_Y(021)=  2.6306077703650351E+02_ark
   index_Y(1,022)= 0; index_Y(2,022)= 5; index_Y(3,022)= 3; d_Y(022)= -6.1891025835930122E+02_ark
   index_Y(1,023)= 0; index_Y(2,023)= 5; index_Y(3,023)= 4; d_Y(023)=  6.3630381759865122E+02_ark
   index_Y(1,024)= 0; index_Y(2,024)= 5; index_Y(3,024)= 5; d_Y(024)= -3.3458737650602416E+02_ark
   index_Y(1,025)= 0; index_Y(2,025)= 5; index_Y(3,025)= 6; d_Y(025)=  8.5441636123447097E+01_ark
   index_Y(1,026)= 0; index_Y(2,026)= 5; index_Y(3,026)= 7; d_Y(026)= -7.0972841656621313E+00_ark
   index_Y(1,027)= 0; index_Y(2,027)= 5; index_Y(3,027)= 8; d_Y(027)= -4.3291673791463836E-01_ark
   index_Y(1,028)= 0; index_Y(2,028)= 7; index_Y(3,028)= 0; d_Y(028)= -6.3949404915501873E-01_ark
   index_Y(1,029)= 0; index_Y(2,029)= 7; index_Y(3,029)= 1; d_Y(029)=  5.3697620428538357E+00_ark
   index_Y(1,030)= 0; index_Y(2,030)= 7; index_Y(3,030)= 2; d_Y(030)= -2.1246715925590252E+01_ark
   index_Y(1,031)= 0; index_Y(2,031)= 7; index_Y(3,031)= 3; d_Y(031)=  4.6150850531621927E+01_ark
   index_Y(1,032)= 0; index_Y(2,032)= 7; index_Y(3,032)= 4; d_Y(032)= -5.8164070783002217E+01_ark
   index_Y(1,033)= 0; index_Y(2,033)= 7; index_Y(3,033)= 5; d_Y(033)=  4.3540851665255104E+01_ark
   index_Y(1,034)= 0; index_Y(2,034)= 7; index_Y(3,034)= 6; d_Y(034)= -1.9065031598205678E+01_ark
   index_Y(1,035)= 0; index_Y(2,035)= 7; index_Y(3,035)= 7; d_Y(035)=  4.5054792348310002E+00_ark
   index_Y(1,036)= 0; index_Y(2,036)= 7; index_Y(3,036)= 8; d_Y(036)= -4.4354261085391045E-01_ark
   index_Y(1,037)= 0; index_Y(2,037)= 9; index_Y(3,037)= 0; d_Y(037)=  9.1961188862033083E-04_ark
   index_Y(1,038)= 0; index_Y(2,038)= 9; index_Y(3,038)= 1; d_Y(038)= -1.0495907744115129E-02_ark
   index_Y(1,039)= 0; index_Y(2,039)= 9; index_Y(3,039)= 2; d_Y(039)=  4.8841840636214329E-02_ark
   index_Y(1,040)= 0; index_Y(2,040)= 9; index_Y(3,040)= 3; d_Y(040)= -1.1364750491850373E-01_ark
   index_Y(1,041)= 0; index_Y(2,041)= 9; index_Y(3,041)= 4; d_Y(041)=  1.4418584905831722E-01_ark
   index_Y(1,042)= 0; index_Y(2,042)= 9; index_Y(3,042)= 5; d_Y(042)= -1.0421563094701014E-01_ark
   index_Y(1,043)= 0; index_Y(2,043)= 9; index_Y(3,043)= 6; d_Y(043)=  4.2765590164208334E-02_ark
   index_Y(1,044)= 0; index_Y(2,044)= 9; index_Y(3,044)= 7; d_Y(044)= -9.2557807083721855E-03_ark
   index_Y(1,045)= 0; index_Y(2,045)= 9; index_Y(3,045)= 8; d_Y(045)=  8.1841915448421787E-04_ark
   index_Y(1,046)= 1; index_Y(2,046)= 1; index_Y(3,046)= 0; d_Y(046)=  1.0828513080615849E+04_ark
   index_Y(1,047)= 1; index_Y(2,047)= 1; index_Y(3,047)= 1; d_Y(047)= -6.1958344528485381E+04_ark
   index_Y(1,048)= 1; index_Y(2,048)= 1; index_Y(3,048)= 2; d_Y(048)=  1.3732693343310739E+05_ark
   index_Y(1,049)= 1; index_Y(2,049)= 1; index_Y(3,049)= 3; d_Y(049)= -1.4602616110645336E+05_ark
   index_Y(1,050)= 1; index_Y(2,050)= 1; index_Y(3,050)= 4; d_Y(050)=  6.8811950508502021E+04_ark
   index_Y(1,051)= 1; index_Y(2,051)= 1; index_Y(3,051)= 5; d_Y(051)=  4.9475135291751940E+02_ark
   index_Y(1,052)= 1; index_Y(2,052)= 1; index_Y(3,052)= 6; d_Y(052)= -1.4177604144980898E+04_ark
   index_Y(1,053)= 1; index_Y(2,053)= 1; index_Y(3,053)= 7; d_Y(053)=  5.3843224207515595E+03_ark
   index_Y(1,054)= 1; index_Y(2,054)= 1; index_Y(3,054)= 8; d_Y(054)= -6.5125864102935884E+02_ark
   index_Y(1,055)= 1; index_Y(2,055)= 3; index_Y(3,055)= 0; d_Y(055)=  1.1608320410891611E+02_ark
   index_Y(1,056)= 1; index_Y(2,056)= 3; index_Y(3,056)= 1; d_Y(056)=  3.0841664668855374E+03_ark
   index_Y(1,057)= 1; index_Y(2,057)= 3; index_Y(3,057)= 2; d_Y(057)= -1.9599848175697436E+04_ark
   index_Y(1,058)= 1; index_Y(2,058)= 3; index_Y(3,058)= 3; d_Y(058)=  4.6880598293958697E+04_ark
   index_Y(1,059)= 1; index_Y(2,059)= 3; index_Y(3,059)= 4; d_Y(059)= -5.8280899986154342E+04_ark
   index_Y(1,060)= 1; index_Y(2,060)= 3; index_Y(3,060)= 5; d_Y(060)=  4.1256074369251844E+04_ark
   index_Y(1,061)= 1; index_Y(2,061)= 3; index_Y(3,061)= 6; d_Y(061)= -1.6751348811883654E+04_ark
   index_Y(1,062)= 1; index_Y(2,062)= 3; index_Y(3,062)= 7; d_Y(062)=  3.6309072374957323E+03_ark
   index_Y(1,063)= 1; index_Y(2,063)= 3; index_Y(3,063)= 8; d_Y(063)= -3.2538908818212803E+02_ark
   index_Y(1,064)= 1; index_Y(2,064)= 5; index_Y(3,064)= 0; d_Y(064)=  1.7127138269743227E+01_ark
   index_Y(1,065)= 1; index_Y(2,065)= 5; index_Y(3,065)= 1; d_Y(065)=  9.7780005264248757E+00_ark
   index_Y(1,066)= 1; index_Y(2,066)= 5; index_Y(3,066)= 2; d_Y(066)= -2.4288699679490583E+02_ark
   index_Y(1,067)= 1; index_Y(2,067)= 5; index_Y(3,067)= 3; d_Y(067)=  5.6674945667351130E+02_ark
   index_Y(1,068)= 1; index_Y(2,068)= 5; index_Y(3,068)= 4; d_Y(068)= -6.0452587426050741E+02_ark
   index_Y(1,069)= 1; index_Y(2,069)= 5; index_Y(3,069)= 5; d_Y(069)=  3.4361181457917701E+02_ark
   index_Y(1,070)= 1; index_Y(2,070)= 5; index_Y(3,070)= 6; d_Y(070)= -1.0348652397692786E+02_ark
   index_Y(1,071)= 1; index_Y(2,071)= 5; index_Y(3,071)= 7; d_Y(071)=  1.4383600738097812E+01_ark
   index_Y(1,072)= 1; index_Y(2,072)= 5; index_Y(3,072)= 8; d_Y(072)= -5.4561409435700625E-01_ark
   index_Y(1,073)= 1; index_Y(2,073)= 7; index_Y(3,073)= 0; d_Y(073)=  2.6381298131173025E-01_ark
   index_Y(1,074)= 1; index_Y(2,074)= 7; index_Y(3,074)= 1; d_Y(074)= -2.2513328844352145E+00_ark
   index_Y(1,075)= 1; index_Y(2,075)= 7; index_Y(3,075)= 2; d_Y(075)=  8.8800878948668469E+00_ark
   index_Y(1,076)= 1; index_Y(2,076)= 7; index_Y(3,076)= 3; d_Y(076)= -1.9109779860387789E+01_ark
   index_Y(1,077)= 1; index_Y(2,077)= 7; index_Y(3,077)= 4; d_Y(077)=  2.3835633796726142E+01_ark
   index_Y(1,078)= 1; index_Y(2,078)= 7; index_Y(3,078)= 5; d_Y(078)= -1.7662596569180096E+01_ark
   index_Y(1,079)= 1; index_Y(2,079)= 7; index_Y(3,079)= 6; d_Y(079)=  7.6588459639606299E+00_ark
   index_Y(1,080)= 1; index_Y(2,080)= 7; index_Y(3,080)= 7; d_Y(080)= -1.7931398108121357E+00_ark
   index_Y(1,081)= 1; index_Y(2,081)= 7; index_Y(3,081)= 8; d_Y(081)=  1.7495153230265714E-01_ark
   index_Y(1,082)= 2; index_Y(2,082)= 1; index_Y(3,082)= 0; d_Y(082)= -8.6874732719412368E+03_ark
   index_Y(1,083)= 2; index_Y(2,083)= 1; index_Y(3,083)= 1; d_Y(083)=  5.1021782358335586E+04_ark
   index_Y(1,084)= 2; index_Y(2,084)= 1; index_Y(3,084)= 2; d_Y(084)= -1.1856528606455067E+05_ark
   index_Y(1,085)= 2; index_Y(2,085)= 1; index_Y(3,085)= 3; d_Y(085)=  1.3832302249277182E+05_ark
   index_Y(1,086)= 2; index_Y(2,086)= 1; index_Y(3,086)= 4; d_Y(086)= -8.3027166901945369E+04_ark
   index_Y(1,087)= 2; index_Y(2,087)= 1; index_Y(3,087)= 5; d_Y(087)=  2.0398308682465227E+04_ark
   index_Y(1,088)= 2; index_Y(2,088)= 1; index_Y(3,088)= 6; d_Y(088)=  2.4654199999068514E+03_ark
   index_Y(1,089)= 2; index_Y(2,089)= 1; index_Y(3,089)= 7; d_Y(089)= -2.2827310582726786E+03_ark
   index_Y(1,090)= 2; index_Y(2,090)= 1; index_Y(3,090)= 8; d_Y(090)=  3.2990680808486650E+02_ark
   index_Y(1,091)= 2; index_Y(2,091)= 3; index_Y(3,091)= 0; d_Y(091)= -4.6754929801885737E+01_ark
   index_Y(1,092)= 2; index_Y(2,092)= 3; index_Y(3,092)= 1; d_Y(092)= -1.6468072940279089E+03_ark
   index_Y(1,093)= 2; index_Y(2,093)= 3; index_Y(3,093)= 2; d_Y(093)=  1.0216542603031499E+04_ark
   index_Y(1,094)= 2; index_Y(2,094)= 3; index_Y(3,094)= 3; d_Y(094)= -2.4336093649657036E+04_ark
   index_Y(1,095)= 2; index_Y(2,095)= 3; index_Y(3,095)= 4; d_Y(095)=  3.0262911938653473E+04_ark
   index_Y(1,096)= 2; index_Y(2,096)= 3; index_Y(3,096)= 5; d_Y(096)= -2.1468007124369775E+04_ark
   index_Y(1,097)= 2; index_Y(2,097)= 3; index_Y(3,097)= 6; d_Y(097)=  8.7443606750326289E+03_ark
   index_Y(1,098)= 2; index_Y(2,098)= 3; index_Y(3,098)= 7; d_Y(098)= -1.9027857210679285E+03_ark
   index_Y(1,099)= 2; index_Y(2,099)= 3; index_Y(3,099)= 8; d_Y(099)=  1.7129168135279906E+02_ark
   index_Y(1,100)= 2; index_Y(2,100)= 5; index_Y(3,100)= 0; d_Y(100)= -4.3764170026406646E+00_ark
   index_Y(1,101)= 2; index_Y(2,101)= 5; index_Y(3,101)= 1; d_Y(101)= -7.8373275434832976E+00_ark
   index_Y(1,102)= 2; index_Y(2,102)= 5; index_Y(3,102)= 2; d_Y(102)=  9.0024360359064303E+01_ark
   index_Y(1,103)= 2; index_Y(2,103)= 5; index_Y(3,103)= 3; d_Y(103)= -2.0651224054530849E+02_ark
   index_Y(1,104)= 2; index_Y(2,104)= 5; index_Y(3,104)= 4; d_Y(104)=  2.2812987148446541E+02_ark
   index_Y(1,105)= 2; index_Y(2,105)= 5; index_Y(3,105)= 5; d_Y(105)= -1.3903806777980640E+02_ark
   index_Y(1,106)= 2; index_Y(2,106)= 5; index_Y(3,106)= 6; d_Y(106)=  4.7211169573288316E+01_ark
   index_Y(1,107)= 2; index_Y(2,107)= 5; index_Y(3,107)= 7; d_Y(107)= -8.2255573207794441E+00_ark
   index_Y(1,108)= 2; index_Y(2,108)= 5; index_Y(3,108)= 8; d_Y(108)=  5.5621079584034305E-01_ark
   index_Y(1,109)= 2; index_Y(2,109)= 7; index_Y(3,109)= 0; d_Y(109)= -2.6543678818597982E-02_ark
   index_Y(1,110)= 2; index_Y(2,110)= 7; index_Y(3,110)= 1; d_Y(110)=  2.3036099053183534E-01_ark
   index_Y(1,111)= 2; index_Y(2,111)= 7; index_Y(3,111)= 2; d_Y(111)= -9.1111843931139447E-01_ark
   index_Y(1,112)= 2; index_Y(2,112)= 7; index_Y(3,112)= 3; d_Y(112)=  1.9544922946695351E+00_ark
   index_Y(1,113)= 2; index_Y(2,113)= 7; index_Y(3,113)= 4; d_Y(113)= -2.4238456777600277E+00_ark
   index_Y(1,114)= 2; index_Y(2,114)= 7; index_Y(3,114)= 5; d_Y(114)=  1.7833819855381989E+00_ark
   index_Y(1,115)= 2; index_Y(2,115)= 7; index_Y(3,115)= 6; d_Y(115)= -7.6722022467765782E-01_ark
   index_Y(1,116)= 2; index_Y(2,116)= 7; index_Y(3,116)= 7; d_Y(116)=  1.7812924495774496E-01_ark
   index_Y(1,117)= 2; index_Y(2,117)= 7; index_Y(3,117)= 8; d_Y(117)= -1.7230348931661865E-02_ark
   index_Y(1,118)= 3; index_Y(2,118)= 1; index_Y(3,118)= 0; d_Y(118)=  3.9383501583444358E+03_ark
   index_Y(1,119)= 3; index_Y(2,119)= 1; index_Y(3,119)= 1; d_Y(119)= -2.3736734872909554E+04_ark
   index_Y(1,120)= 3; index_Y(2,120)= 1; index_Y(3,120)= 2; d_Y(120)=  5.7577324569642020E+04_ark
   index_Y(1,121)= 3; index_Y(2,121)= 1; index_Y(3,121)= 3; d_Y(121)= -7.2284295018403602E+04_ark
   index_Y(1,122)= 3; index_Y(2,122)= 1; index_Y(3,122)= 4; d_Y(122)=  5.0149441167938770E+04_ark
   index_Y(1,123)= 3; index_Y(2,123)= 1; index_Y(3,123)= 5; d_Y(123)= -1.8552433495144898E+04_ark
   index_Y(1,124)= 3; index_Y(2,124)= 1; index_Y(3,124)= 6; d_Y(124)=  2.8520833363890270E+03_ark
   index_Y(1,125)= 3; index_Y(2,125)= 1; index_Y(3,125)= 7; d_Y(125)=  1.3054776132589905E+02_ark
   index_Y(1,126)= 3; index_Y(2,126)= 1; index_Y(3,126)= 8; d_Y(126)= -6.4350459709648931E+01_ark
   index_Y(1,127)= 3; index_Y(2,127)= 3; index_Y(3,127)= 0; d_Y(127)=  9.2717552110279939E+00_ark
   index_Y(1,128)= 3; index_Y(2,128)= 3; index_Y(3,128)= 1; d_Y(128)=  4.5338721075168723E+02_ark
   index_Y(1,129)= 3; index_Y(2,129)= 3; index_Y(3,129)= 2; d_Y(129)= -2.7579176608353373E+03_ark
   index_Y(1,130)= 3; index_Y(2,130)= 3; index_Y(3,130)= 3; d_Y(130)=  6.5508926793801147E+03_ark
   index_Y(1,131)= 3; index_Y(2,131)= 3; index_Y(3,131)= 4; d_Y(131)= -8.1530881467083163E+03_ark
   index_Y(1,132)= 3; index_Y(2,132)= 3; index_Y(3,132)= 5; d_Y(132)=  5.7967755981790360E+03_ark
   index_Y(1,133)= 3; index_Y(2,133)= 3; index_Y(3,133)= 6; d_Y(133)= -2.3683363245503715E+03_ark
   index_Y(1,134)= 3; index_Y(2,134)= 3; index_Y(3,134)= 7; d_Y(134)=  5.1719063463986822E+02_ark
   index_Y(1,135)= 3; index_Y(2,135)= 3; index_Y(3,135)= 8; d_Y(135)= -4.6742688110094605E+01_ark
   index_Y(1,136)= 3; index_Y(2,136)= 5; index_Y(3,136)= 0; d_Y(136)=  4.3784929804974126E-01_ark
   index_Y(1,137)= 3; index_Y(2,137)= 5; index_Y(3,137)= 1; d_Y(137)=  2.0124058596093164E+00_ark
   index_Y(1,138)= 3; index_Y(2,138)= 5; index_Y(3,138)= 2; d_Y(138)= -1.5382753893549307E+01_ark
   index_Y(1,139)= 3; index_Y(2,139)= 5; index_Y(3,139)= 3; d_Y(139)=  3.4566424204806481E+01_ark
   index_Y(1,140)= 3; index_Y(2,140)= 5; index_Y(3,140)= 4; d_Y(140)= -3.9202473142326880E+01_ark
   index_Y(1,141)= 3; index_Y(2,141)= 5; index_Y(3,141)= 5; d_Y(141)=  2.5132562380268382E+01_ark
   index_Y(1,142)= 3; index_Y(2,142)= 5; index_Y(3,142)= 6; d_Y(142)= -9.1976071586091166E+00_ark
   index_Y(1,143)= 3; index_Y(2,143)= 5; index_Y(3,143)= 7; d_Y(143)=  1.7859862244968099E+00_ark
   index_Y(1,144)= 3; index_Y(2,144)= 5; index_Y(3,144)= 8; d_Y(144)= -1.4221126718348387E-01_ark
   index_Y(1,145)= 4; index_Y(2,145)= 1; index_Y(3,145)= 0; d_Y(145)= -1.1027156149868561E+03_ark
   index_Y(1,146)= 4; index_Y(2,146)= 1; index_Y(3,146)= 1; d_Y(146)=  6.8161331466427855E+03_ark
   index_Y(1,147)= 4; index_Y(2,147)= 1; index_Y(3,147)= 2; d_Y(147)= -1.7182440121492898E+04_ark
   index_Y(1,148)= 4; index_Y(2,148)= 1; index_Y(3,148)= 3; d_Y(148)=  2.2875789474071840E+04_ark
   index_Y(1,149)= 4; index_Y(2,149)= 1; index_Y(3,149)= 4; d_Y(149)= -1.7467745131210548E+04_ark
   index_Y(1,150)= 4; index_Y(2,150)= 1; index_Y(3,150)= 5; d_Y(150)=  7.7326116020025074E+03_ark
   index_Y(1,151)= 4; index_Y(2,151)= 1; index_Y(3,151)= 6; d_Y(151)= -1.8776612942374941E+03_ark
   index_Y(1,152)= 4; index_Y(2,152)= 1; index_Y(3,152)= 7; d_Y(152)=  2.0853550685340088E+02_ark
   index_Y(1,153)= 4; index_Y(2,153)= 1; index_Y(3,153)= 8; d_Y(153)= -5.0149630309797431E+00_ark
   index_Y(1,154)= 4; index_Y(2,154)= 3; index_Y(3,154)= 0; d_Y(154)= -9.9114520260036443E-01_ark
   index_Y(1,155)= 4; index_Y(2,155)= 3; index_Y(3,155)= 1; d_Y(155)= -6.7136728118731298E+01_ark
   index_Y(1,156)= 4; index_Y(2,156)= 3; index_Y(3,156)= 2; d_Y(156)=  4.0357173145310298E+02_ark
   index_Y(1,157)= 4; index_Y(2,157)= 3; index_Y(3,157)= 3; d_Y(157)= -9.5823170142405161E+02_ark
   index_Y(1,158)= 4; index_Y(2,158)= 3; index_Y(3,158)= 4; d_Y(158)=  1.1948539306631246E+03_ark
   index_Y(1,159)= 4; index_Y(2,159)= 3; index_Y(3,159)= 5; d_Y(159)= -8.5181525608713855E+02_ark
   index_Y(1,160)= 4; index_Y(2,160)= 3; index_Y(3,160)= 6; d_Y(160)=  3.4908880621907701E+02_ark
   index_Y(1,161)= 4; index_Y(2,161)= 3; index_Y(3,161)= 7; d_Y(161)= -7.6484976427967922E+01_ark
   index_Y(1,162)= 4; index_Y(2,162)= 3; index_Y(3,162)= 8; d_Y(162)=  6.9364222380518186E+00_ark
   index_Y(1,163)= 4; index_Y(2,163)= 5; index_Y(3,163)= 0; d_Y(163)= -1.2746893289577343E-02_ark
   index_Y(1,164)= 4; index_Y(2,164)= 5; index_Y(3,164)= 1; d_Y(164)= -1.6283774396723061E-01_ark
   index_Y(1,165)= 4; index_Y(2,165)= 5; index_Y(3,165)= 2; d_Y(165)=  9.8732533594710503E-01_ark
   index_Y(1,166)= 4; index_Y(2,166)= 5; index_Y(3,166)= 3; d_Y(166)= -2.1780451753992907E+00_ark
   index_Y(1,167)= 4; index_Y(2,167)= 5; index_Y(3,167)= 4; d_Y(167)=  2.5143243279706127E+00_ark
   index_Y(1,168)= 4; index_Y(2,168)= 5; index_Y(3,168)= 5; d_Y(168)= -1.6666167524719526E+00_ark
   index_Y(1,169)= 4; index_Y(2,169)= 5; index_Y(3,169)= 6; d_Y(169)=  6.3782156683242874E-01_ark
   index_Y(1,170)= 4; index_Y(2,170)= 5; index_Y(3,170)= 7; d_Y(170)= -1.3098288697327121E-01_ark
   index_Y(1,171)= 4; index_Y(2,171)= 5; index_Y(3,171)= 8; d_Y(171)=  1.1173600748094259E-02_ark
   index_Y(1,172)= 5; index_Y(2,172)= 1; index_Y(3,172)= 0; d_Y(172)=  1.9511508123662236E+02_ark
   index_Y(1,173)= 5; index_Y(2,173)= 1; index_Y(3,173)= 1; d_Y(173)= -1.2354985857407928E+03_ark
   index_Y(1,174)= 5; index_Y(2,174)= 1; index_Y(3,174)= 2; d_Y(174)=  3.2224016656602016E+03_ark
   index_Y(1,175)= 5; index_Y(2,175)= 1; index_Y(3,175)= 3; d_Y(175)= -4.4974293239414765E+03_ark
   index_Y(1,176)= 5; index_Y(2,176)= 1; index_Y(3,176)= 4; d_Y(176)=  3.6726601751274707E+03_ark
   index_Y(1,177)= 5; index_Y(2,177)= 1; index_Y(3,177)= 5; d_Y(177)= -1.7978261556599946E+03_ark
   index_Y(1,178)= 5; index_Y(2,178)= 1; index_Y(3,178)= 6; d_Y(178)=  5.1439224865405663E+02_ark
   index_Y(1,179)= 5; index_Y(2,179)= 1; index_Y(3,179)= 7; d_Y(179)= -7.8173752132894776E+01_ark
   index_Y(1,180)= 5; index_Y(2,180)= 1; index_Y(3,180)= 8; d_Y(180)=  4.7582169005224273E+00_ark
   index_Y(1,181)= 5; index_Y(2,181)= 3; index_Y(3,181)= 0; d_Y(181)=  6.5238887102054832E-02_ark
   index_Y(1,182)= 5; index_Y(2,182)= 3; index_Y(3,182)= 1; d_Y(182)=  4.9742931408716800E+00_ark
   index_Y(1,183)= 5; index_Y(2,183)= 3; index_Y(3,183)= 2; d_Y(183)= -2.9972307095736369E+01_ark
   index_Y(1,184)= 5; index_Y(2,184)= 3; index_Y(3,184)= 3; d_Y(184)=  7.1459537936572701E+01_ark
   index_Y(1,185)= 5; index_Y(2,185)= 3; index_Y(3,185)= 4; d_Y(185)= -8.9454713279142197E+01_ark
   index_Y(1,186)= 5; index_Y(2,186)= 3; index_Y(3,186)= 5; d_Y(186)=  6.4003785143104835E+01_ark
   index_Y(1,187)= 5; index_Y(2,187)= 3; index_Y(3,187)= 6; d_Y(187)= -2.6319108118808060E+01_ark
   index_Y(1,188)= 5; index_Y(2,188)= 3; index_Y(3,188)= 7; d_Y(188)=  5.7850389033271199E+00_ark
   index_Y(1,189)= 5; index_Y(2,189)= 3; index_Y(3,189)= 8; d_Y(189)= -5.2623643022423039E-01_ark
   index_Y(1,190)= 6; index_Y(2,190)= 1; index_Y(3,190)= 0; d_Y(190)= -2.1284827840736583E+01_ark
   index_Y(1,191)= 6; index_Y(2,191)= 1; index_Y(3,191)= 1; d_Y(191)=  1.3783558054829109E+02_ark
   index_Y(1,192)= 6; index_Y(2,192)= 1; index_Y(3,192)= 2; d_Y(192)= -3.7028385993336815E+02_ark
   index_Y(1,193)= 6; index_Y(2,193)= 1; index_Y(3,193)= 3; d_Y(193)=  5.3670058217117173E+02_ark
   index_Y(1,194)= 6; index_Y(2,194)= 1; index_Y(3,194)= 4; d_Y(194)= -4.6005254619420379E+02_ark
   index_Y(1,195)= 6; index_Y(2,195)= 1; index_Y(3,195)= 5; d_Y(195)=  2.3986456471729525E+02_ark
   index_Y(1,196)= 6; index_Y(2,196)= 1; index_Y(3,196)= 6; d_Y(196)= -7.4621583549456147E+01_ark
   index_Y(1,197)= 6; index_Y(2,197)= 1; index_Y(3,197)= 7; d_Y(197)=  1.2713815864941807E+01_ark
   index_Y(1,198)= 6; index_Y(2,198)= 1; index_Y(3,198)= 8; d_Y(198)= -9.1098804145342172E-01_ark
   index_Y(1,199)= 6; index_Y(2,199)= 3; index_Y(3,199)= 0; d_Y(199)= -2.7266201597160311E-03_ark
   index_Y(1,200)= 6; index_Y(2,200)= 3; index_Y(3,200)= 1; d_Y(200)= -1.3890418236048951E-01_ark
   index_Y(1,201)= 6; index_Y(2,201)= 3; index_Y(3,201)= 2; d_Y(201)=  8.6302004193719739E-01_ark
   index_Y(1,202)= 6; index_Y(2,202)= 3; index_Y(3,202)= 3; d_Y(202)= -2.0841598203457945E+00_ark
   index_Y(1,203)= 6; index_Y(2,203)= 3; index_Y(3,203)= 4; d_Y(203)=  2.6294525324186946E+00_ark
   index_Y(1,204)= 6; index_Y(2,204)= 3; index_Y(3,204)= 5; d_Y(204)= -1.8918128255971371E+00_ark
   index_Y(1,205)= 6; index_Y(2,205)= 3; index_Y(3,205)= 6; d_Y(205)=  7.8128164115712462E-01_ark
   index_Y(1,206)= 6; index_Y(2,206)= 3; index_Y(3,206)= 7; d_Y(206)= -1.7232173401906614E-01_ark
   index_Y(1,207)= 6; index_Y(2,207)= 3; index_Y(3,207)= 8; d_Y(207)=  1.5718670163984427E-02_ark
   index_Y(1,208)= 7; index_Y(2,208)= 1; index_Y(3,208)= 0; d_Y(208)=  1.3072513891947799E+00_ark
   index_Y(1,209)= 7; index_Y(2,209)= 1; index_Y(3,209)= 1; d_Y(209)= -8.6372154169610553E+00_ark
   index_Y(1,210)= 7; index_Y(2,210)= 1; index_Y(3,210)= 2; d_Y(210)=  2.3789523722660253E+01_ark
   index_Y(1,211)= 7; index_Y(2,211)= 1; index_Y(3,211)= 3; d_Y(211)= -3.5526970635376472E+01_ark
   index_Y(1,212)= 7; index_Y(2,212)= 1; index_Y(3,212)= 4; d_Y(212)=  3.1552252207458725E+01_ark
   index_Y(1,213)= 7; index_Y(2,213)= 1; index_Y(3,213)= 5; d_Y(213)= -1.7154825873124022E+01_ark
   index_Y(1,214)= 7; index_Y(2,214)= 1; index_Y(3,214)= 6; d_Y(214)=  5.6066546500266927E+00_ark
   index_Y(1,215)= 7; index_Y(2,215)= 1; index_Y(3,215)= 7; d_Y(215)= -1.0120573737961336E+00_ark
   index_Y(1,216)= 7; index_Y(2,216)= 1; index_Y(3,216)= 8; d_Y(216)=  7.7570472862957707E-02_ark
   index_Y(1,217)= 8; index_Y(2,217)= 1; index_Y(3,217)= 0; d_Y(217)= -3.4555418066994373E-02_ark
   index_Y(1,218)= 8; index_Y(2,218)= 1; index_Y(3,218)= 1; d_Y(218)=  2.3223833862609206E-01_ark
   index_Y(1,219)= 8; index_Y(2,219)= 1; index_Y(3,219)= 2; d_Y(219)= -6.5271843232022964E-01_ark
   index_Y(1,220)= 8; index_Y(2,220)= 1; index_Y(3,220)= 3; d_Y(220)=  9.9740864964466303E-01_ark
   index_Y(1,221)= 8; index_Y(2,221)= 1; index_Y(3,221)= 4; d_Y(221)= -9.0888159851258443E-01_ark
   index_Y(1,222)= 8; index_Y(2,222)= 1; index_Y(3,222)= 5; d_Y(222)=  5.0839987505135653E-01_ark
   index_Y(1,223)= 8; index_Y(2,223)= 1; index_Y(3,223)= 6; d_Y(223)= -1.7139035613840276E-01_ark
   index_Y(1,224)= 8; index_Y(2,224)= 1; index_Y(3,224)= 7; d_Y(224)=  3.1985031944868778E-02_ark
   index_Y(1,225)= 8; index_Y(2,225)= 1; index_Y(3,225)= 8; d_Y(225)= -2.5392866543105371E-03_ark

 elseif(dipole_surface == "LTP2011NR") then
   n_r_X = 8; n_theta_X = 8; ncoeffs_X = 200
   index_X(1,001)= 0; index_X(2,001)= 0; index_X(3,001)= 1; d_X(001)= -3.1856223546828005E+01_ark
   index_X(1,002)= 0; index_X(2,002)= 0; index_X(3,002)= 2; d_X(002)=  3.1925671906666976E+03_ark
   index_X(1,003)= 0; index_X(2,003)= 0; index_X(3,003)= 3; d_X(003)= -1.2420131648763749E+04_ark
   index_X(1,004)= 0; index_X(2,004)= 0; index_X(3,004)= 4; d_X(004)=  2.0075642943638664E+04_ark
   index_X(1,005)= 0; index_X(2,005)= 0; index_X(3,005)= 5; d_X(005)= -1.7035726816875936E+04_ark
   index_X(1,006)= 0; index_X(2,006)= 0; index_X(3,006)= 6; d_X(006)=  7.9746213922284333E+03_ark
   index_X(1,007)= 0; index_X(2,007)= 0; index_X(3,007)= 7; d_X(007)= -1.9493867522180699E+03_ark
   index_X(1,008)= 0; index_X(2,008)= 0; index_X(3,008)= 8; d_X(008)=  1.9428164249905967E+02_ark
   index_X(1,009)= 0; index_X(2,009)= 2; index_X(3,009)= 1; d_X(009)=  2.9924522206937309E+02_ark
   index_X(1,010)= 0; index_X(2,010)= 2; index_X(3,010)= 2; d_X(010)= -1.1569756680752980E+03_ark
   index_X(1,011)= 0; index_X(2,011)= 2; index_X(3,011)= 3; d_X(011)=  1.9832648185496291E+03_ark
   index_X(1,012)= 0; index_X(2,012)= 2; index_X(3,012)= 4; d_X(012)= -1.9797615328614920E+03_ark
   index_X(1,013)= 0; index_X(2,013)= 2; index_X(3,013)= 5; d_X(013)=  1.2565441096005234E+03_ark
   index_X(1,014)= 0; index_X(2,014)= 2; index_X(3,014)= 6; d_X(014)= -5.0546976311040635E+02_ark
   index_X(1,015)= 0; index_X(2,015)= 2; index_X(3,015)= 7; d_X(015)=  1.1692673396212194E+02_ark
   index_X(1,016)= 0; index_X(2,016)= 2; index_X(3,016)= 8; d_X(016)= -1.1674410487112254E+01_ark
   index_X(1,017)= 0; index_X(2,017)= 4; index_X(3,017)= 1; d_X(017)=  1.6996178709804553E+01_ark
   index_X(1,018)= 0; index_X(2,018)= 4; index_X(3,018)= 2; d_X(018)= -4.4541118073244434E+01_ark
   index_X(1,019)= 0; index_X(2,019)= 4; index_X(3,019)= 3; d_X(019)=  2.6566083807165342E+01_ark
   index_X(1,020)= 0; index_X(2,020)= 4; index_X(3,020)= 4; d_X(020)=  3.0480209957349871E+01_ark
   index_X(1,021)= 0; index_X(2,021)= 4; index_X(3,021)= 5; d_X(021)= -5.2196416202423279E+01_ark
   index_X(1,022)= 0; index_X(2,022)= 4; index_X(3,022)= 6; d_X(022)=  2.9720320733140397E+01_ark
   index_X(1,023)= 0; index_X(2,023)= 4; index_X(3,023)= 7; d_X(023)= -7.6306119003202184E+00_ark
   index_X(1,024)= 0; index_X(2,024)= 4; index_X(3,024)= 8; d_X(024)=  7.3740685389020655E-01_ark
   index_X(1,025)= 0; index_X(2,025)= 6; index_X(3,025)= 1; d_X(025)=  1.6293169456628220E-01_ark
   index_X(1,026)= 0; index_X(2,026)= 6; index_X(3,026)= 2; d_X(026)= -1.4070879553987652E+00_ark
   index_X(1,027)= 0; index_X(2,027)= 6; index_X(3,027)= 3; d_X(027)=  4.4959166166463547E+00_ark
   index_X(1,028)= 0; index_X(2,028)= 6; index_X(3,028)= 4; d_X(028)= -7.1377949545012598E+00_ark
   index_X(1,029)= 0; index_X(2,029)= 6; index_X(3,029)= 5; d_X(029)=  6.2195188419691476E+00_ark
   index_X(1,030)= 0; index_X(2,030)= 6; index_X(3,030)= 6; d_X(030)= -3.0239509686398378E+00_ark
   index_X(1,031)= 0; index_X(2,031)= 6; index_X(3,031)= 7; d_X(031)=  7.6786069603986107E-01_ark
   index_X(1,032)= 0; index_X(2,032)= 6; index_X(3,032)= 8; d_X(032)= -7.9069733024880406E-02_ark
   index_X(1,033)= 0; index_X(2,033)= 8; index_X(3,033)= 1; d_X(033)=  2.8070402693458618E-04_ark
   index_X(1,034)= 0; index_X(2,034)= 8; index_X(3,034)= 2; d_X(034)= -4.1218232621815787E-03_ark
   index_X(1,035)= 0; index_X(2,035)= 8; index_X(3,035)= 3; d_X(035)=  1.5121787172574841E-02_ark
   index_X(1,036)= 0; index_X(2,036)= 8; index_X(3,036)= 4; d_X(036)= -2.5391562452227845E-02_ark
   index_X(1,037)= 0; index_X(2,037)= 8; index_X(3,037)= 5; d_X(037)=  2.2891065717885795E-02_ark
   index_X(1,038)= 0; index_X(2,038)= 8; index_X(3,038)= 6; d_X(038)= -1.1452798402387998E-02_ark
   index_X(1,039)= 0; index_X(2,039)= 8; index_X(3,039)= 7; d_X(039)=  2.9962127675275951E-03_ark
   index_X(1,040)= 0; index_X(2,040)= 8; index_X(3,040)= 8; d_X(040)= -3.1964150639396394E-04_ark
   index_X(1,041)= 1; index_X(2,041)= 0; index_X(3,041)= 1; d_X(041)=  1.5373003291362602E+02_ark
   index_X(1,042)= 1; index_X(2,042)= 0; index_X(3,042)= 2; d_X(042)= -6.4057280427812584E+03_ark
   index_X(1,043)= 1; index_X(2,043)= 0; index_X(3,043)= 3; d_X(043)=  2.4073340020083182E+04_ark
   index_X(1,044)= 1; index_X(2,044)= 0; index_X(3,044)= 4; d_X(044)= -3.8441814811672084E+04_ark
   index_X(1,045)= 1; index_X(2,045)= 0; index_X(3,045)= 5; d_X(045)=  3.2424005520887815E+04_ark
   index_X(1,046)= 1; index_X(2,046)= 0; index_X(3,046)= 6; d_X(046)= -1.5125172580457649E+04_ark
   index_X(1,047)= 1; index_X(2,047)= 0; index_X(3,047)= 7; d_X(047)=  3.6892392361755774E+03_ark
   index_X(1,048)= 1; index_X(2,048)= 0; index_X(3,048)= 8; d_X(048)= -3.6713571347937977E+02_ark
   index_X(1,049)= 1; index_X(2,049)= 2; index_X(3,049)= 1; d_X(049)= -3.7675889069706500E+02_ark
   index_X(1,050)= 1; index_X(2,050)= 2; index_X(3,050)= 2; d_X(050)=  1.4639580399207916E+03_ark
   index_X(1,051)= 1; index_X(2,051)= 2; index_X(3,051)= 3; d_X(051)= -2.5256132448586141E+03_ark
   index_X(1,052)= 1; index_X(2,052)= 2; index_X(3,052)= 4; d_X(052)=  2.5365317180957645E+03_ark
   index_X(1,053)= 1; index_X(2,053)= 2; index_X(3,053)= 5; d_X(053)= -1.6155044149936002E+03_ark
   index_X(1,054)= 1; index_X(2,054)= 2; index_X(3,054)= 6; d_X(054)=  6.4942157546253293E+02_ark
   index_X(1,055)= 1; index_X(2,055)= 2; index_X(3,055)= 7; d_X(055)= -1.4951414477602520E+02_ark
   index_X(1,056)= 1; index_X(2,056)= 2; index_X(3,056)= 8; d_X(056)=  1.4813093447301071E+01_ark
   index_X(1,057)= 1; index_X(2,057)= 4; index_X(3,057)= 1; d_X(057)= -1.3381456732675360E+01_ark
   index_X(1,058)= 1; index_X(2,058)= 4; index_X(3,058)= 2; d_X(058)=  3.4809161489940379E+01_ark
   index_X(1,059)= 1; index_X(2,059)= 4; index_X(3,059)= 3; d_X(059)= -1.9451065659275628E+01_ark
   index_X(1,060)= 1; index_X(2,060)= 4; index_X(3,060)= 4; d_X(060)= -2.7010078130702823E+01_ark
   index_X(1,061)= 1; index_X(2,061)= 4; index_X(3,061)= 5; d_X(061)=  4.4144473738557281E+01_ark
   index_X(1,062)= 1; index_X(2,062)= 4; index_X(3,062)= 6; d_X(062)= -2.5038527058686668E+01_ark
   index_X(1,063)= 1; index_X(2,063)= 4; index_X(3,063)= 7; d_X(063)=  6.4561519604467321E+00_ark
   index_X(1,064)= 1; index_X(2,064)= 4; index_X(3,064)= 8; d_X(064)= -6.2957802882010583E-01_ark
   index_X(1,065)= 1; index_X(2,065)= 6; index_X(3,065)= 1; d_X(065)= -4.7513303288951647E-02_ark
   index_X(1,066)= 1; index_X(2,066)= 6; index_X(3,066)= 2; d_X(066)=  4.1317954123815070E-01_ark
   index_X(1,067)= 1; index_X(2,067)= 6; index_X(3,067)= 3; d_X(067)= -1.3350380923966441E+00_ark
   index_X(1,068)= 1; index_X(2,068)= 6; index_X(3,068)= 4; d_X(068)=  2.1414445422810786E+00_ark
   index_X(1,069)= 1; index_X(2,069)= 6; index_X(3,069)= 5; d_X(069)= -1.8816184095858262E+00_ark
   index_X(1,070)= 1; index_X(2,070)= 6; index_X(3,070)= 6; d_X(070)=  9.2107846989347308E-01_ark
   index_X(1,071)= 1; index_X(2,071)= 6; index_X(3,071)= 7; d_X(071)= -2.3521121966768987E-01_ark
   index_X(1,072)= 1; index_X(2,072)= 6; index_X(3,072)= 8; d_X(072)=  2.4338397091923980E-02_ark
   index_X(1,073)= 2; index_X(2,073)= 0; index_X(3,073)= 1; d_X(073)= -2.0253432908627656E+02_ark
   index_X(1,074)= 2; index_X(2,074)= 0; index_X(3,074)= 2; d_X(074)=  5.5391320275401958E+03_ark
   index_X(1,075)= 2; index_X(2,075)= 0; index_X(3,075)= 3; d_X(075)= -2.0151367965662172E+04_ark
   index_X(1,076)= 2; index_X(2,076)= 0; index_X(3,076)= 4; d_X(076)=  3.1797329736150969E+04_ark
   index_X(1,077)= 2; index_X(2,077)= 0; index_X(3,077)= 5; d_X(077)= -2.6658369492740610E+04_ark
   index_X(1,078)= 2; index_X(2,078)= 0; index_X(3,078)= 6; d_X(078)=  1.2392088578403942E+04_ark
   index_X(1,079)= 2; index_X(2,079)= 0; index_X(3,079)= 7; d_X(079)= -3.0158955090298696E+03_ark
   index_X(1,080)= 2; index_X(2,080)= 0; index_X(3,080)= 8; d_X(080)=  2.9966922229060583E+02_ark
   index_X(1,081)= 2; index_X(2,081)= 2; index_X(3,081)= 1; d_X(081)=  1.9270779047346332E+02_ark
   index_X(1,082)= 2; index_X(2,082)= 2; index_X(3,082)= 2; d_X(082)= -7.5374240366585582E+02_ark
   index_X(1,083)= 2; index_X(2,083)= 2; index_X(3,083)= 3; d_X(083)=  1.3118567256053502E+03_ark
   index_X(1,084)= 2; index_X(2,084)= 2; index_X(3,084)= 4; d_X(084)= -1.3294912831195543E+03_ark
   index_X(1,085)= 2; index_X(2,085)= 2; index_X(3,085)= 5; d_X(085)=  8.5235981285723392E+02_ark
   index_X(1,086)= 2; index_X(2,086)= 2; index_X(3,086)= 6; d_X(086)= -3.4337758213098573E+02_ark
   index_X(1,087)= 2; index_X(2,087)= 2; index_X(3,087)= 7; d_X(087)=  7.8845389692316530E+01_ark
   index_X(1,088)= 2; index_X(2,088)= 2; index_X(3,088)= 8; d_X(088)= -7.7607080209563719E+00_ark
   index_X(1,089)= 2; index_X(2,089)= 4; index_X(3,089)= 1; d_X(089)=  3.9049971675840425E+00_ark
   index_X(1,090)= 2; index_X(2,090)= 4; index_X(3,090)= 2; d_X(090)= -1.0145874291866676E+01_ark
   index_X(1,091)= 2; index_X(2,091)= 4; index_X(3,091)= 3; d_X(091)=  5.4277729576847378E+00_ark
   index_X(1,092)= 2; index_X(2,092)= 4; index_X(3,092)= 4; d_X(092)=  8.5558531323423495E+00_ark
   index_X(1,093)= 2; index_X(2,093)= 4; index_X(3,093)= 5; d_X(093)= -1.3650188667025304E+01_ark
   index_X(1,094)= 2; index_X(2,094)= 4; index_X(3,094)= 6; d_X(094)=  7.7457206622493686E+00_ark
   index_X(1,095)= 2; index_X(2,095)= 4; index_X(3,095)= 7; d_X(095)= -2.0081043747613876E+00_ark
   index_X(1,096)= 2; index_X(2,096)= 4; index_X(3,096)= 8; d_X(096)=  1.9746219270382426E-01_ark
   index_X(1,097)= 2; index_X(2,097)= 6; index_X(3,097)= 1; d_X(097)=  3.5218514701327308E-03_ark
   index_X(1,098)= 2; index_X(2,098)= 6; index_X(3,098)= 2; d_X(098)= -2.9197748432068238E-02_ark
   index_X(1,099)= 2; index_X(2,099)= 6; index_X(3,099)= 3; d_X(099)=  9.3399784229820426E-02_ark
   index_X(1,100)= 2; index_X(2,100)= 6; index_X(3,100)= 4; d_X(100)= -1.4995339815152420E-01_ark
   index_X(1,101)= 2; index_X(2,101)= 6; index_X(3,101)= 5; d_X(101)=  1.3217254926189526E-01_ark
   index_X(1,102)= 2; index_X(2,102)= 6; index_X(3,102)= 6; d_X(102)= -6.4899411247552052E-02_ark
   index_X(1,103)= 2; index_X(2,103)= 6; index_X(3,103)= 7; d_X(103)=  1.6612182991138980E-02_ark
   index_X(1,104)= 2; index_X(2,104)= 6; index_X(3,104)= 8; d_X(104)= -1.7212415089034039E-03_ark
   index_X(1,105)= 3; index_X(2,105)= 0; index_X(3,105)= 1; d_X(105)=  1.2913183641677597E+02_ark
   index_X(1,106)= 3; index_X(2,106)= 0; index_X(3,106)= 2; d_X(106)= -2.6930578694100168E+03_ark
   index_X(1,107)= 3; index_X(2,107)= 0; index_X(3,107)= 3; d_X(107)=  9.5046349081418484E+03_ark
   index_X(1,108)= 3; index_X(2,108)= 0; index_X(3,108)= 4; d_X(108)= -1.4824732491712377E+04_ark
   index_X(1,109)= 3; index_X(2,109)= 0; index_X(3,109)= 5; d_X(109)=  1.2355183698162784E+04_ark
   index_X(1,110)= 3; index_X(2,110)= 0; index_X(3,110)= 6; d_X(110)= -5.7233299176658975E+03_ark
   index_X(1,111)= 3; index_X(2,111)= 0; index_X(3,111)= 7; d_X(111)=  1.3898040542378367E+03_ark
   index_X(1,112)= 3; index_X(2,112)= 0; index_X(3,112)= 8; d_X(112)= -1.3787988718284487E+02_ark
   index_X(1,113)= 3; index_X(2,113)= 2; index_X(3,113)= 1; d_X(113)= -5.1081129858614986E+01_ark
   index_X(1,114)= 3; index_X(2,114)= 2; index_X(3,114)= 2; d_X(114)=  2.0144326403617742E+02_ark
   index_X(1,115)= 3; index_X(2,115)= 2; index_X(3,115)= 3; d_X(115)= -3.5440994024697648E+02_ark
   index_X(1,116)= 3; index_X(2,116)= 2; index_X(3,116)= 4; d_X(116)=  3.6328484820175800E+02_ark
   index_X(1,117)= 3; index_X(2,117)= 2; index_X(3,117)= 5; d_X(117)= -2.3504449634318189E+02_ark
   index_X(1,118)= 3; index_X(2,118)= 2; index_X(3,118)= 6; d_X(118)=  9.5120552492981005E+01_ark
   index_X(1,119)= 3; index_X(2,119)= 2; index_X(3,119)= 7; d_X(119)= -2.1826541775715668E+01_ark
   index_X(1,120)= 3; index_X(2,120)= 2; index_X(3,120)= 8; d_X(120)=  2.1371898429588327E+00_ark
   index_X(1,121)= 3; index_X(2,121)= 4; index_X(3,121)= 1; d_X(121)= -5.0246929995154233E-01_ark
   index_X(1,122)= 3; index_X(2,122)= 4; index_X(3,122)= 2; d_X(122)=  1.3207642869515439E+00_ark
   index_X(1,123)= 3; index_X(2,123)= 4; index_X(3,123)= 3; d_X(123)= -7.2781535606350189E-01_ark
   index_X(1,124)= 3; index_X(2,124)= 4; index_X(3,124)= 4; d_X(124)= -1.0939424243233589E+00_ark
   index_X(1,125)= 3; index_X(2,125)= 4; index_X(3,125)= 5; d_X(125)=  1.7790921554027364E+00_ark
   index_X(1,126)= 3; index_X(2,126)= 4; index_X(3,126)= 6; d_X(126)= -1.0186242369028946E+00_ark
   index_X(1,127)= 3; index_X(2,127)= 4; index_X(3,127)= 7; d_X(127)=  2.6610726284081920E-01_ark
   index_X(1,128)= 3; index_X(2,128)= 4; index_X(3,128)= 8; d_X(128)= -2.6366449318516061E-02_ark
   index_X(1,129)= 4; index_X(2,129)= 0; index_X(3,129)= 1; d_X(129)= -4.6804510637386215E+01_ark
   index_X(1,130)= 4; index_X(2,130)= 0; index_X(3,130)= 2; d_X(130)=  8.0409703386131741E+02_ark
   index_X(1,131)= 4; index_X(2,131)= 0; index_X(3,131)= 3; d_X(131)= -2.7593303703850152E+03_ark
   index_X(1,132)= 4; index_X(2,132)= 0; index_X(3,132)= 4; d_X(132)=  4.2562421909951690E+03_ark
   index_X(1,133)= 4; index_X(2,133)= 0; index_X(3,133)= 5; d_X(133)= -3.5268283327869794E+03_ark
   index_X(1,134)= 4; index_X(2,134)= 0; index_X(3,134)= 6; d_X(134)=  1.6281985563157816E+03_ark
   index_X(1,135)= 4; index_X(2,135)= 0; index_X(3,135)= 7; d_X(135)= -3.9450935406747647E+02_ark
   index_X(1,136)= 4; index_X(2,136)= 0; index_X(3,136)= 8; d_X(136)=  3.9076814468183329E+01_ark
   index_X(1,137)= 4; index_X(2,137)= 2; index_X(3,137)= 1; d_X(137)=  7.3707943891934633E+00_ark
   index_X(1,138)= 4; index_X(2,138)= 2; index_X(3,138)= 2; d_X(138)= -2.9341221609692639E+01_ark
   index_X(1,139)= 4; index_X(2,139)= 2; index_X(3,139)= 3; d_X(139)=  5.2250417267000103E+01_ark
   index_X(1,140)= 4; index_X(2,140)= 2; index_X(3,140)= 4; d_X(140)= -5.4259653388165361E+01_ark
   index_X(1,141)= 4; index_X(2,141)= 2; index_X(3,141)= 5; d_X(141)=  3.5494594466831586E+01_ark
   index_X(1,142)= 4; index_X(2,142)= 2; index_X(3,142)= 6; d_X(142)= -1.4458071608977434E+01_ark
   index_X(1,143)= 4; index_X(2,143)= 2; index_X(3,143)= 7; d_X(143)=  3.3210993124366723E+00_ark
   index_X(1,144)= 4; index_X(2,144)= 2; index_X(3,144)= 8; d_X(144)= -3.2389237839811358E-01_ark
   index_X(1,145)= 4; index_X(2,145)= 4; index_X(3,145)= 1; d_X(145)=  2.3998610833167788E-02_ark
   index_X(1,146)= 4; index_X(2,146)= 4; index_X(3,146)= 2; d_X(146)= -6.4691426544641217E-02_ark
   index_X(1,147)= 4; index_X(2,147)= 4; index_X(3,147)= 3; d_X(147)=  3.9711989450999852E-02_ark
   index_X(1,148)= 4; index_X(2,148)= 4; index_X(3,148)= 4; d_X(148)=  4.6004022307073456E-02_ark
   index_X(1,149)= 4; index_X(2,149)= 4; index_X(3,149)= 5; d_X(149)= -8.0853422638398342E-02_ark
   index_X(1,150)= 4; index_X(2,150)= 4; index_X(3,150)= 6; d_X(150)=  4.7166682057593334E-02_ark
   index_X(1,151)= 4; index_X(2,151)= 4; index_X(3,151)= 7; d_X(151)= -1.2433618903614274E-02_ark
   index_X(1,152)= 4; index_X(2,152)= 4; index_X(3,152)= 8; d_X(152)=  1.2382037677980406E-03_ark
   index_X(1,153)= 5; index_X(2,153)= 0; index_X(3,153)= 1; d_X(153)=  1.0155984421619394E+01_ark
   index_X(1,154)= 5; index_X(2,154)= 0; index_X(3,154)= 2; d_X(154)= -1.5074553779922746E+02_ark
   index_X(1,155)= 5; index_X(2,155)= 0; index_X(3,155)= 3; d_X(155)=  5.0419538027665544E+02_ark
   index_X(1,156)= 5; index_X(2,156)= 0; index_X(3,156)= 4; d_X(156)= -7.6961446078857739E+02_ark
   index_X(1,157)= 5; index_X(2,157)= 0; index_X(3,157)= 5; d_X(157)=  6.3423459295398879E+02_ark
   index_X(1,158)= 5; index_X(2,158)= 0; index_X(3,158)= 6; d_X(158)= -2.9185084020960710E+02_ark
   index_X(1,159)= 5; index_X(2,159)= 0; index_X(3,159)= 7; d_X(159)=  7.0564545751516391E+01_ark
   index_X(1,160)= 5; index_X(2,160)= 0; index_X(3,160)= 8; d_X(160)= -6.9785608764290146E+00_ark
   index_X(1,161)= 5; index_X(2,161)= 2; index_X(3,161)= 1; d_X(161)= -5.4616773751840952E-01_ark
   index_X(1,162)= 5; index_X(2,162)= 2; index_X(3,162)= 2; d_X(162)=  2.1951780536307979E+00_ark
   index_X(1,163)= 5; index_X(2,163)= 2; index_X(3,163)= 3; d_X(163)= -3.9586131606004358E+00_ark
   index_X(1,164)= 5; index_X(2,164)= 2; index_X(3,164)= 4; d_X(164)=  4.1684562839570845E+00_ark
   index_X(1,165)= 5; index_X(2,165)= 2; index_X(3,165)= 5; d_X(165)= -2.7607401534485181E+00_ark
   index_X(1,166)= 5; index_X(2,166)= 2; index_X(3,166)= 6; d_X(166)=  1.1336488829973064E+00_ark
   index_X(1,167)= 5; index_X(2,167)= 2; index_X(3,167)= 7; d_X(167)= -2.6106089391979026E-01_ark
   index_X(1,168)= 5; index_X(2,168)= 2; index_X(3,168)= 8; d_X(168)=  2.5381814023432980E-02_ark
   index_X(1,169)= 6; index_X(2,169)= 0; index_X(3,169)= 1; d_X(169)= -1.3048766184702743E+00_ark
   index_X(1,170)= 6; index_X(2,170)= 0; index_X(3,170)= 2; d_X(170)=  1.7297283283916400E+01_ark
   index_X(1,171)= 6; index_X(2,171)= 0; index_X(3,171)= 3; d_X(171)= -5.6539929686259853E+01_ark
   index_X(1,172)= 6; index_X(2,172)= 0; index_X(3,172)= 4; d_X(172)=  8.5477538346249077E+01_ark
   index_X(1,173)= 6; index_X(2,173)= 0; index_X(3,173)= 5; d_X(173)= -7.0085131970198404E+01_ark
   index_X(1,174)= 6; index_X(2,174)= 0; index_X(3,174)= 6; d_X(174)=  3.2153164619975051E+01_ark
   index_X(1,175)= 6; index_X(2,175)= 0; index_X(3,175)= 7; d_X(175)= -7.7584716951285841E+00_ark
   index_X(1,176)= 6; index_X(2,176)= 0; index_X(3,176)= 8; d_X(176)=  7.6610206344729903E-01_ark
   index_X(1,177)= 6; index_X(2,177)= 2; index_X(3,177)= 1; d_X(177)=  1.6135086065826698E-02_ark
   index_X(1,178)= 6; index_X(2,178)= 2; index_X(3,178)= 2; d_X(178)= -6.5411995854831595E-02_ark
   index_X(1,179)= 6; index_X(2,179)= 2; index_X(3,179)= 3; d_X(179)=  1.1938369266735727E-01_ark
   index_X(1,180)= 6; index_X(2,180)= 2; index_X(3,180)= 4; d_X(180)= -1.2748840721739496E-01_ark
   index_X(1,181)= 6; index_X(2,181)= 2; index_X(3,181)= 5; d_X(181)=  8.5544011523523089E-02_ark
   index_X(1,182)= 6; index_X(2,182)= 2; index_X(3,182)= 6; d_X(182)= -3.5446945807093622E-02_ark
   index_X(1,183)= 6; index_X(2,183)= 2; index_X(3,183)= 7; d_X(183)=  8.1901856569501774E-03_ark
   index_X(1,184)= 6; index_X(2,184)= 2; index_X(3,184)= 8; d_X(184)= -7.9392288078405926E-04_ark
   index_X(1,185)= 7; index_X(2,185)= 0; index_X(3,185)= 1; d_X(185)=  9.1297847047918879E-02_ark
   index_X(1,186)= 7; index_X(2,186)= 0; index_X(3,186)= 2; d_X(186)= -1.1085362214654140E+00_ark
   index_X(1,187)= 7; index_X(2,187)= 0; index_X(3,187)= 3; d_X(187)=  3.5518472576224696E+00_ark
   index_X(1,188)= 7; index_X(2,188)= 0; index_X(3,188)= 4; d_X(188)= -5.3240810510830148E+00_ark
   index_X(1,189)= 7; index_X(2,189)= 0; index_X(3,189)= 5; d_X(189)=  4.3456428198497656E+00_ark
   index_X(1,190)= 7; index_X(2,190)= 0; index_X(3,190)= 6; d_X(190)= -1.9882595308045687E+00_ark
   index_X(1,191)= 7; index_X(2,191)= 0; index_X(3,191)= 7; d_X(191)=  4.7887660449723068E-01_ark
   index_X(1,192)= 7; index_X(2,192)= 0; index_X(3,192)= 8; d_X(192)= -4.7215822080697301E-02_ark
   index_X(1,193)= 8; index_X(2,193)= 0; index_X(3,193)= 1; d_X(193)= -2.6704121255463759E-03_ark
   index_X(1,194)= 8; index_X(2,194)= 0; index_X(3,194)= 2; d_X(194)=  3.0321455872724945E-02_ark
   index_X(1,195)= 8; index_X(2,195)= 0; index_X(3,195)= 3; d_X(195)= -9.5548340704489088E-02_ark
   index_X(1,196)= 8; index_X(2,196)= 0; index_X(3,196)= 4; d_X(196)=  1.4219126816072344E-01_ark
   index_X(1,197)= 8; index_X(2,197)= 0; index_X(3,197)= 5; d_X(197)= -1.1561371674879437E-01_ark
   index_X(1,198)= 8; index_X(2,198)= 0; index_X(3,198)= 6; d_X(198)=  5.2773082341791476E-02_ark
   index_X(1,199)= 8; index_X(2,199)= 0; index_X(3,199)= 7; d_X(199)= -1.2689504798380772E-02_ark
   index_X(1,200)= 8; index_X(2,200)= 0; index_X(3,200)= 8; d_X(200)=  1.2493518660616036E-03_ark
!********************************************************************************************
   n_r_Y = 9; n_theta_Y = 8; ncoeffs_Y = 225
   index_Y(1,001)= 0; index_Y(2,001)= 1; index_Y(3,001)= 0; d_Y(001)= -6.0072518557566873E+03_ark
   index_Y(1,002)= 0; index_Y(2,002)= 1; index_Y(3,002)= 1; d_Y(002)=  3.9310547921326004E+04_ark
   index_Y(1,003)= 0; index_Y(2,003)= 1; index_Y(3,003)= 2; d_Y(003)= -1.0584983864639507E+05_ark
   index_Y(1,004)= 0; index_Y(2,004)= 1; index_Y(3,004)= 3; d_Y(004)=  1.5066930955435100E+05_ark
   index_Y(1,005)= 0; index_Y(2,005)= 1; index_Y(3,005)= 4; d_Y(005)= -1.2352272503669836E+05_ark
   index_Y(1,006)= 0; index_Y(2,006)= 1; index_Y(3,006)= 5; d_Y(006)=  5.9671252931204042E+04_ark
   index_Y(1,007)= 0; index_Y(2,007)= 1; index_Y(3,007)= 6; d_Y(007)= -1.6497346837913385E+04_ark
   index_Y(1,008)= 0; index_Y(2,008)= 1; index_Y(3,008)= 7; d_Y(008)=  2.3457420732531464E+03_ark
   index_Y(1,009)= 0; index_Y(2,009)= 1; index_Y(3,009)= 8; d_Y(009)= -1.2533957796311006E+02_ark
   index_Y(1,010)= 0; index_Y(2,010)= 3; index_Y(3,010)= 0; d_Y(010)= -1.0292625077195407E+02_ark
   index_Y(1,011)= 0; index_Y(2,011)= 3; index_Y(3,011)= 1; d_Y(011)= -2.7318944343247858E+03_ark
   index_Y(1,012)= 0; index_Y(2,012)= 3; index_Y(3,012)= 2; d_Y(012)=  1.6253719547989313E+04_ark
   index_Y(1,013)= 0; index_Y(2,013)= 3; index_Y(3,013)= 3; d_Y(013)= -3.6169333503210859E+04_ark
   index_Y(1,014)= 0; index_Y(2,014)= 3; index_Y(3,014)= 4; d_Y(014)=  4.1475166034037131E+04_ark
   index_Y(1,015)= 0; index_Y(2,015)= 3; index_Y(3,015)= 5; d_Y(015)= -2.6807891692309029E+04_ark
   index_Y(1,016)= 0; index_Y(2,016)= 3; index_Y(3,016)= 6; d_Y(016)=  9.8146511252467753E+03_ark
   index_Y(1,017)= 0; index_Y(2,017)= 3; index_Y(3,017)= 7; d_Y(017)= -1.8867412841079058E+03_ark
   index_Y(1,018)= 0; index_Y(2,018)= 3; index_Y(3,018)= 8; d_Y(018)=  1.4655442674201913E+02_ark
   index_Y(1,019)= 0; index_Y(2,019)= 5; index_Y(3,019)= 0; d_Y(019)= -2.1804423586163921E+01_ark
   index_Y(1,020)= 0; index_Y(2,020)= 5; index_Y(3,020)= 1; d_Y(020)= -2.3838938174914801E+01_ark
   index_Y(1,021)= 0; index_Y(2,021)= 5; index_Y(3,021)= 2; d_Y(021)=  3.8445859857864696E+02_ark
   index_Y(1,022)= 0; index_Y(2,022)= 5; index_Y(3,022)= 3; d_Y(022)= -9.1243031912392325E+02_ark
   index_Y(1,023)= 0; index_Y(2,023)= 5; index_Y(3,023)= 4; d_Y(023)=  1.0227834007132551E+03_ark
   index_Y(1,024)= 0; index_Y(2,024)= 5; index_Y(3,024)= 5; d_Y(024)= -6.3261296192170994E+02_ark
   index_Y(1,025)= 0; index_Y(2,025)= 5; index_Y(3,025)= 6; d_Y(025)=  2.2052165623220208E+02_ark
   index_Y(1,026)= 0; index_Y(2,026)= 5; index_Y(3,026)= 7; d_Y(026)= -4.0438182465048158E+01_ark
   index_Y(1,027)= 0; index_Y(2,027)= 5; index_Y(3,027)= 8; d_Y(027)=  3.0172063735117263E+00_ark
   index_Y(1,028)= 0; index_Y(2,028)= 7; index_Y(3,028)= 0; d_Y(028)= -6.9830142543514739E-01_ark
   index_Y(1,029)= 0; index_Y(2,029)= 7; index_Y(3,029)= 1; d_Y(029)=  6.6435617076140261E+00_ark
   index_Y(1,030)= 0; index_Y(2,030)= 7; index_Y(3,030)= 2; d_Y(030)= -2.8092537597858609E+01_ark
   index_Y(1,031)= 0; index_Y(2,031)= 7; index_Y(3,031)= 3; d_Y(031)=  6.2249367190927387E+01_ark
   index_Y(1,032)= 0; index_Y(2,032)= 7; index_Y(3,032)= 4; d_Y(032)= -7.8378036412153051E+01_ark
   index_Y(1,033)= 0; index_Y(2,033)= 7; index_Y(3,033)= 5; d_Y(033)=  5.8124605370828704E+01_ark
   index_Y(1,034)= 0; index_Y(2,034)= 7; index_Y(3,034)= 6; d_Y(034)= -2.5129782974150203E+01_ark
   index_Y(1,035)= 0; index_Y(2,035)= 7; index_Y(3,035)= 7; d_Y(035)=  5.8569345320720458E+00_ark
   index_Y(1,036)= 0; index_Y(2,036)= 7; index_Y(3,036)= 8; d_Y(036)= -5.6850842548010405E-01_ark
   index_Y(1,037)= 0; index_Y(2,037)= 9; index_Y(3,037)= 0; d_Y(037)=  1.0294489539688367E-03_ark
   index_Y(1,038)= 0; index_Y(2,038)= 9; index_Y(3,038)= 1; d_Y(038)= -1.0864782603732692E-02_ark
   index_Y(1,039)= 0; index_Y(2,039)= 9; index_Y(3,039)= 2; d_Y(039)=  4.7504901712784431E-02_ark
   index_Y(1,040)= 0; index_Y(2,040)= 9; index_Y(3,040)= 3; d_Y(040)= -1.0617931096371080E-01_ark
   index_Y(1,041)= 0; index_Y(2,041)= 9; index_Y(3,041)= 4; d_Y(041)=  1.3108431779573948E-01_ark
   index_Y(1,042)= 0; index_Y(2,042)= 9; index_Y(3,042)= 5; d_Y(042)= -9.2839924004692875E-02_ark
   index_Y(1,043)= 0; index_Y(2,043)= 9; index_Y(3,043)= 6; d_Y(043)=  3.7460961480974220E-02_ark
   index_Y(1,044)= 0; index_Y(2,044)= 9; index_Y(3,044)= 7; d_Y(044)= -7.9819588909231243E-03_ark
   index_Y(1,045)= 0; index_Y(2,045)= 9; index_Y(3,045)= 8; d_Y(045)=  6.9450026444428659E-04_ark
   index_Y(1,046)= 1; index_Y(2,046)= 1; index_Y(3,046)= 0; d_Y(046)=  1.1116268738017519E+04_ark
   index_Y(1,047)= 1; index_Y(2,047)= 1; index_Y(3,047)= 1; d_Y(047)= -7.3865012597599925E+04_ark
   index_Y(1,048)= 1; index_Y(2,048)= 1; index_Y(3,048)= 2; d_Y(048)=  2.0266998769250905E+05_ark
   index_Y(1,049)= 1; index_Y(2,049)= 1; index_Y(3,049)= 3; d_Y(049)= -2.9510932195499167E+05_ark
   index_Y(1,050)= 1; index_Y(2,050)= 1; index_Y(3,050)= 4; d_Y(050)=  2.4887093516196194E+05_ark
   index_Y(1,051)= 1; index_Y(2,051)= 1; index_Y(3,051)= 5; d_Y(051)= -1.2472484049404738E+05_ark
   index_Y(1,052)= 1; index_Y(2,052)= 1; index_Y(3,052)= 6; d_Y(052)=  3.6278230853491346E+04_ark
   index_Y(1,053)= 1; index_Y(2,053)= 1; index_Y(3,053)= 7; d_Y(053)= -5.5684424631429138E+03_ark
   index_Y(1,054)= 1; index_Y(2,054)= 1; index_Y(3,054)= 8; d_Y(054)=  3.3976655071845744E+02_ark
   index_Y(1,055)= 1; index_Y(2,055)= 3; index_Y(3,055)= 0; d_Y(055)=  1.1029842191474745E+02_ark
   index_Y(1,056)= 1; index_Y(2,056)= 3; index_Y(3,056)= 1; d_Y(056)=  3.4823528419403010E+03_ark
   index_Y(1,057)= 1; index_Y(2,057)= 3; index_Y(3,057)= 2; d_Y(057)= -2.0378200872703572E+04_ark
   index_Y(1,058)= 1; index_Y(2,058)= 3; index_Y(3,058)= 3; d_Y(058)=  4.5226819432465360E+04_ark
   index_Y(1,059)= 1; index_Y(2,059)= 3; index_Y(3,059)= 4; d_Y(059)= -5.1928274336436531E+04_ark
   index_Y(1,060)= 1; index_Y(2,060)= 3; index_Y(3,060)= 5; d_Y(060)=  3.3688027555214707E+04_ark
   index_Y(1,061)= 1; index_Y(2,061)= 3; index_Y(3,061)= 6; d_Y(061)= -1.2406012416177662E+04_ark
   index_Y(1,062)= 1; index_Y(2,062)= 3; index_Y(3,062)= 7; d_Y(062)=  2.4050833530513337E+03_ark
   index_Y(1,063)= 1; index_Y(2,063)= 3; index_Y(3,063)= 8; d_Y(063)= -1.8907444483361178E+02_ark
   index_Y(1,064)= 1; index_Y(2,064)= 5; index_Y(3,064)= 0; d_Y(064)=  1.5781728775919873E+01_ark
   index_Y(1,065)= 1; index_Y(2,065)= 5; index_Y(3,065)= 1; d_Y(065)=  3.2182354955502888E+01_ark
   index_Y(1,066)= 1; index_Y(2,066)= 5; index_Y(3,066)= 2; d_Y(066)= -3.5761332050095916E+02_ark
   index_Y(1,067)= 1; index_Y(2,067)= 5; index_Y(3,067)= 3; d_Y(067)=  8.3533919283618161E+02_ark
   index_Y(1,068)= 1; index_Y(2,068)= 5; index_Y(3,068)= 4; d_Y(068)= -9.4713217996976891E+02_ark
   index_Y(1,069)= 1; index_Y(2,069)= 5; index_Y(3,069)= 5; d_Y(069)=  5.9977651685834280E+02_ark
   index_Y(1,070)= 1; index_Y(2,070)= 5; index_Y(3,070)= 6; d_Y(070)= -2.1624390952468821E+02_ark
   index_Y(1,071)= 1; index_Y(2,071)= 5; index_Y(3,071)= 7; d_Y(071)=  4.1485267560290595E+01_ark
   index_Y(1,072)= 1; index_Y(2,072)= 5; index_Y(3,072)= 8; d_Y(072)= -3.2863011558947619E+00_ark
   index_Y(1,073)= 1; index_Y(2,073)= 7; index_Y(3,073)= 0; d_Y(073)=  2.8542530666891253E-01_ark
   index_Y(1,074)= 1; index_Y(2,074)= 7; index_Y(3,074)= 1; d_Y(074)= -2.6837165902397828E+00_ark
   index_Y(1,075)= 1; index_Y(2,075)= 7; index_Y(3,075)= 2; d_Y(075)=  1.1159312553892960E+01_ark
   index_Y(1,076)= 1; index_Y(2,076)= 7; index_Y(3,076)= 3; d_Y(076)= -2.4433433178578525E+01_ark
   index_Y(1,077)= 1; index_Y(2,077)= 7; index_Y(3,077)= 4; d_Y(077)=  3.0500799611896127E+01_ark
   index_Y(1,078)= 1; index_Y(2,078)= 7; index_Y(3,078)= 5; d_Y(078)= -2.2463309823910095E+01_ark
   index_Y(1,079)= 1; index_Y(2,079)= 7; index_Y(3,079)= 6; d_Y(079)=  9.6524535689859476E+00_ark
   index_Y(1,080)= 1; index_Y(2,080)= 7; index_Y(3,080)= 7; d_Y(080)= -2.2367088140963460E+00_ark
   index_Y(1,081)= 1; index_Y(2,081)= 7; index_Y(3,081)= 8; d_Y(081)=  2.1589325017703231E-01_ark
   index_Y(1,082)= 2; index_Y(2,082)= 1; index_Y(3,082)= 0; d_Y(082)= -8.9027557048708950E+03_ark
   index_Y(1,083)= 2; index_Y(2,083)= 1; index_Y(3,083)= 1; d_Y(083)=  6.0084634519222847E+04_ark
   index_Y(1,084)= 2; index_Y(2,084)= 1; index_Y(3,084)= 2; d_Y(084)= -1.6786760345288322E+05_ark
   index_Y(1,085)= 2; index_Y(2,085)= 1; index_Y(3,085)= 3; d_Y(085)=  2.4957224041247845E+05_ark
   index_Y(1,086)= 2; index_Y(2,086)= 1; index_Y(3,086)= 4; d_Y(086)= -2.1572198659986624E+05_ark
   index_Y(1,087)= 2; index_Y(2,087)= 1; index_Y(3,087)= 5; d_Y(087)=  1.1143215726769256E+05_ark
   index_Y(1,088)= 2; index_Y(2,088)= 1; index_Y(3,088)= 6; d_Y(088)= -3.3691396399045712E+04_ark
   index_Y(1,089)= 2; index_Y(2,089)= 1; index_Y(3,089)= 7; d_Y(089)=  5.4497697207188758E+03_ark
   index_Y(1,090)= 2; index_Y(2,090)= 1; index_Y(3,090)= 8; d_Y(090)= -3.5914003819326172E+02_ark
   index_Y(1,091)= 2; index_Y(2,091)= 3; index_Y(3,091)= 0; d_Y(091)= -4.6876434105492081E+01_ark
   index_Y(1,092)= 2; index_Y(2,092)= 3; index_Y(3,092)= 1; d_Y(092)= -1.7968672728624661E+03_ark
   index_Y(1,093)= 2; index_Y(2,093)= 3; index_Y(3,093)= 2; d_Y(093)=  1.0351073956804525E+04_ark
   index_Y(1,094)= 2; index_Y(2,094)= 3; index_Y(3,094)= 3; d_Y(094)= -2.2904726270319370E+04_ark
   index_Y(1,095)= 2; index_Y(2,095)= 3; index_Y(3,095)= 4; d_Y(095)=  2.6311761069566099E+04_ark
   index_Y(1,096)= 2; index_Y(2,096)= 3; index_Y(3,096)= 5; d_Y(096)= -1.7110649912898836E+04_ark
   index_Y(1,097)= 2; index_Y(2,097)= 3; index_Y(3,097)= 6; d_Y(097)=  6.3263009915406728E+03_ark
   index_Y(1,098)= 2; index_Y(2,098)= 3; index_Y(3,098)= 7; d_Y(098)= -1.2334192738493075E+03_ark
   index_Y(1,099)= 2; index_Y(2,099)= 3; index_Y(3,099)= 8; d_Y(099)=  9.7731813974067336E+01_ark
   index_Y(1,100)= 2; index_Y(2,100)= 5; index_Y(3,100)= 0; d_Y(100)= -3.9060158322909047E+00_ark
   index_Y(1,101)= 2; index_Y(2,101)= 5; index_Y(3,101)= 1; d_Y(101)= -1.5790796532925015E+01_ark
   index_Y(1,102)= 2; index_Y(2,102)= 5; index_Y(3,102)= 2; d_Y(102)=  1.3000494940493445E+02_ark
   index_Y(1,103)= 2; index_Y(2,103)= 5; index_Y(3,103)= 3; d_Y(103)= -2.9815466145099981E+02_ark
   index_Y(1,104)= 2; index_Y(2,104)= 5; index_Y(3,104)= 4; d_Y(104)=  3.4251759952258726E+02_ark
   index_Y(1,105)= 2; index_Y(2,105)= 5; index_Y(3,105)= 5; d_Y(105)= -2.2265129389798881E+02_ark
   index_Y(1,106)= 2; index_Y(2,106)= 5; index_Y(3,106)= 6; d_Y(106)=  8.3167440640643690E+01_ark
   index_Y(1,107)= 2; index_Y(2,107)= 5; index_Y(3,107)= 7; d_Y(107)= -1.6671316934145921E+01_ark
   index_Y(1,108)= 2; index_Y(2,108)= 5; index_Y(3,108)= 8; d_Y(108)=  1.3921097961947453E+00_ark
   index_Y(1,109)= 2; index_Y(2,109)= 7; index_Y(3,109)= 0; d_Y(109)= -2.8559268304860552E-02_ark
   index_Y(1,110)= 2; index_Y(2,110)= 7; index_Y(3,110)= 1; d_Y(110)=  2.6702770424401479E-01_ark
   index_Y(1,111)= 2; index_Y(2,111)= 7; index_Y(3,111)= 2; d_Y(111)= -1.0989831652502744E+00_ark
   index_Y(1,112)= 2; index_Y(2,112)= 7; index_Y(3,112)= 3; d_Y(112)=  2.3877906080517732E+00_ark
   index_Y(1,113)= 2; index_Y(2,113)= 7; index_Y(3,113)= 4; d_Y(113)= -2.9624728072617472E+00_ark
   index_Y(1,114)= 2; index_Y(2,114)= 7; index_Y(3,114)= 5; d_Y(114)=  2.1695059595949715E+00_ark
   index_Y(1,115)= 2; index_Y(2,115)= 7; index_Y(3,115)= 6; d_Y(115)= -9.2699011528384290E-01_ark
   index_Y(1,116)= 2; index_Y(2,116)= 7; index_Y(3,116)= 7; d_Y(116)=  2.1356723511325981E-01_ark
   index_Y(1,117)= 2; index_Y(2,117)= 7; index_Y(3,117)= 8; d_Y(117)= -2.0491739016506472E-02_ark
   index_Y(1,118)= 3; index_Y(2,118)= 1; index_Y(3,118)= 0; d_Y(118)=  4.0278120122464061E+03_ark
   index_Y(1,119)= 3; index_Y(2,119)= 1; index_Y(3,119)= 1; d_Y(119)= -2.7609585315189965E+04_ark
   index_Y(1,120)= 3; index_Y(2,120)= 1; index_Y(3,120)= 2; d_Y(120)=  7.8463341023883346E+04_ark
   index_Y(1,121)= 3; index_Y(2,121)= 1; index_Y(3,121)= 3; d_Y(121)= -1.1886151475099110E+05_ark
   index_Y(1,122)= 3; index_Y(2,122)= 1; index_Y(3,122)= 4; d_Y(122)=  1.0494990515341470E+05_ark
   index_Y(1,123)= 3; index_Y(2,123)= 1; index_Y(3,123)= 5; d_Y(123)= -5.5577892654059484E+04_ark
   index_Y(1,124)= 3; index_Y(2,124)= 1; index_Y(3,124)= 6; d_Y(124)=  1.7315899373600783E+04_ark
   index_Y(1,125)= 3; index_Y(2,125)= 1; index_Y(3,125)= 7; d_Y(125)= -2.9084135180046287E+03_ark
   index_Y(1,126)= 3; index_Y(2,126)= 1; index_Y(3,126)= 8; d_Y(126)=  2.0144249819262586E+02_ark
   index_Y(1,127)= 3; index_Y(2,127)= 3; index_Y(3,127)= 0; d_Y(127)=  1.0312283677809319E+01_ark
   index_Y(1,128)= 3; index_Y(2,128)= 3; index_Y(3,128)= 1; d_Y(128)=  4.7637855601825868E+02_ark
   index_Y(1,129)= 3; index_Y(2,129)= 3; index_Y(3,129)= 2; d_Y(129)= -2.7091909382408630E+03_ark
   index_Y(1,130)= 3; index_Y(2,130)= 3; index_Y(3,130)= 3; d_Y(130)=  5.9776687977355905E+03_ark
   index_Y(1,131)= 3; index_Y(2,131)= 3; index_Y(3,131)= 4; d_Y(131)= -6.8647042637799750E+03_ark
   index_Y(1,132)= 3; index_Y(2,132)= 3; index_Y(3,132)= 5; d_Y(132)=  4.4682003891880777E+03_ark
   index_Y(1,133)= 3; index_Y(2,133)= 3; index_Y(3,133)= 6; d_Y(133)= -1.6548831479721775E+03_ark
   index_Y(1,134)= 3; index_Y(2,134)= 3; index_Y(3,134)= 7; d_Y(134)=  3.2343937513413039E+02_ark
   index_Y(1,135)= 3; index_Y(2,135)= 3; index_Y(3,135)= 8; d_Y(135)= -2.5711377731706307E+01_ark
   index_Y(1,136)= 3; index_Y(2,136)= 5; index_Y(3,136)= 0; d_Y(136)=  3.6611331468486696E-01_ark
   index_Y(1,137)= 3; index_Y(2,137)= 5; index_Y(3,137)= 1; d_Y(137)=  3.2236766187124886E+00_ark
   index_Y(1,138)= 3; index_Y(2,138)= 5; index_Y(3,138)= 2; d_Y(138)= -2.1380208602150105E+01_ark
   index_Y(1,139)= 3; index_Y(2,139)= 5; index_Y(3,139)= 3; d_Y(139)=  4.8111149580382971E+01_ark
   index_Y(1,140)= 3; index_Y(2,140)= 5; index_Y(3,140)= 4; d_Y(140)= -5.5858518789572827E+01_ark
   index_Y(1,141)= 3; index_Y(2,141)= 5; index_Y(3,141)= 5; d_Y(141)=  3.7113821546199915E+01_ark
   index_Y(1,142)= 3; index_Y(2,142)= 5; index_Y(3,142)= 6; d_Y(142)= -1.4260954405167013E+01_ark
   index_Y(1,143)= 3; index_Y(2,143)= 5; index_Y(3,143)= 7; d_Y(143)=  2.9538316429084261E+00_ark
   index_Y(1,144)= 3; index_Y(2,144)= 5; index_Y(3,144)= 8; d_Y(144)= -2.5572456334612070E-01_ark
   index_Y(1,145)= 4; index_Y(2,145)= 1; index_Y(3,145)= 0; d_Y(145)= -1.1251319340900045E+03_ark
   index_Y(1,146)= 4; index_Y(2,146)= 1; index_Y(3,146)= 1; d_Y(146)=  7.8301037060816379E+03_ark
   index_Y(1,147)= 4; index_Y(2,147)= 1; index_Y(3,147)= 2; d_Y(147)= -2.2604102905523239E+04_ark
   index_Y(1,148)= 4; index_Y(2,148)= 1; index_Y(3,148)= 3; d_Y(148)=  3.4813299305138469E+04_ark
   index_Y(1,149)= 4; index_Y(2,149)= 1; index_Y(3,149)= 4; d_Y(149)= -3.1298927010465628E+04_ark
   index_Y(1,150)= 4; index_Y(2,150)= 1; index_Y(3,150)= 5; d_Y(150)=  1.6914216783573960E+04_ark
   index_Y(1,151)= 4; index_Y(2,151)= 1; index_Y(3,151)= 6; d_Y(151)= -5.3942012049888890E+03_ark
   index_Y(1,152)= 4; index_Y(2,152)= 1; index_Y(3,152)= 7; d_Y(152)=  9.3141800289566345E+02_ark
   index_Y(1,153)= 4; index_Y(2,153)= 1; index_Y(3,153)= 8; d_Y(153)= -6.6741896029820055E+01_ark
   index_Y(1,154)= 4; index_Y(2,154)= 3; index_Y(3,154)= 0; d_Y(154)= -1.3294614192263907E+00_ark
   index_Y(1,155)= 4; index_Y(2,155)= 3; index_Y(3,155)= 1; d_Y(155)= -6.7462924491915146E+01_ark
   index_Y(1,156)= 4; index_Y(2,156)= 3; index_Y(3,156)= 2; d_Y(156)=  3.8125911251763318E+02_ark
   index_Y(1,157)= 4; index_Y(2,157)= 3; index_Y(3,157)= 3; d_Y(157)= -8.3950079924139391E+02_ark
   index_Y(1,158)= 4; index_Y(2,158)= 3; index_Y(3,158)= 4; d_Y(158)=  9.6290029694249915E+02_ark
   index_Y(1,159)= 4; index_Y(2,159)= 3; index_Y(3,159)= 5; d_Y(159)= -6.2604360516728320E+02_ark
   index_Y(1,160)= 4; index_Y(2,160)= 3; index_Y(3,160)= 6; d_Y(160)=  2.3154394062224856E+02_ark
   index_Y(1,161)= 4; index_Y(2,161)= 3; index_Y(3,161)= 7; d_Y(161)= -4.5162950730127704E+01_ark
   index_Y(1,162)= 4; index_Y(2,162)= 3; index_Y(3,162)= 8; d_Y(162)=  3.5789584815520357E+00_ark
   index_Y(1,163)= 4; index_Y(2,163)= 5; index_Y(3,163)= 0; d_Y(163)= -8.7214118116953898E-03_ark
   index_Y(1,164)= 4; index_Y(2,164)= 5; index_Y(3,164)= 1; d_Y(164)= -2.2967971234683660E-01_ark
   index_Y(1,165)= 4; index_Y(2,165)= 5; index_Y(3,165)= 2; d_Y(165)=  1.3132737764881881E+00_ark
   index_Y(1,166)= 4; index_Y(2,166)= 5; index_Y(3,166)= 3; d_Y(166)= -2.9052000128460094E+00_ark
   index_Y(1,167)= 4; index_Y(2,167)= 5; index_Y(3,167)= 4; d_Y(167)=  3.3983417490039471E+00_ark
   index_Y(1,168)= 4; index_Y(2,168)= 5; index_Y(3,168)= 5; d_Y(168)= -2.2949001213207794E+00_ark
   index_Y(1,169)= 4; index_Y(2,169)= 5; index_Y(3,169)= 6; d_Y(169)=  8.9983892195283488E-01_ark
   index_Y(1,170)= 4; index_Y(2,170)= 5; index_Y(3,170)= 7; d_Y(170)= -1.9055901773506889E-01_ark
   index_Y(1,171)= 4; index_Y(2,171)= 5; index_Y(3,171)= 8; d_Y(171)=  1.6879944665291191E-02_ark
   index_Y(1,172)= 5; index_Y(2,172)= 1; index_Y(3,172)= 0; d_Y(172)=  1.9854198727827600E+02_ark
   index_Y(1,173)= 5; index_Y(2,173)= 1; index_Y(3,173)= 1; d_Y(173)= -1.4015858436711642E+03_ark
   index_Y(1,174)= 5; index_Y(2,174)= 1; index_Y(3,174)= 2; d_Y(174)=  4.1031130630746902E+03_ark
   index_Y(1,175)= 5; index_Y(2,175)= 1; index_Y(3,175)= 3; d_Y(175)= -6.4094168431446515E+03_ark
   index_Y(1,176)= 5; index_Y(2,176)= 1; index_Y(3,176)= 4; d_Y(176)=  5.8490352933328522E+03_ark
   index_Y(1,177)= 5; index_Y(2,177)= 1; index_Y(3,177)= 5; d_Y(177)= -3.2123983128496311E+03_ark
   index_Y(1,178)= 5; index_Y(2,178)= 1; index_Y(3,178)= 6; d_Y(178)=  1.0429908791067919E+03_ark
   index_Y(1,179)= 5; index_Y(2,179)= 1; index_Y(3,179)= 7; d_Y(179)= -1.8378368896904669E+02_ark
   index_Y(1,180)= 5; index_Y(2,180)= 1; index_Y(3,180)= 8; d_Y(180)=  1.3484101145528030E+01_ark
   index_Y(1,181)= 5; index_Y(2,181)= 3; index_Y(3,181)= 0; d_Y(181)=  1.0926515275625093E-01_ark
   index_Y(1,182)= 5; index_Y(2,182)= 3; index_Y(3,182)= 1; d_Y(182)=  4.7108770306180077E+00_ark
   index_Y(1,183)= 5; index_Y(2,183)= 3; index_Y(3,183)= 2; d_Y(183)= -2.6836738626227202E+01_ark
   index_Y(1,184)= 5; index_Y(2,184)= 3; index_Y(3,184)= 3; d_Y(184)=  5.9103520229404012E+01_ark
   index_Y(1,185)= 5; index_Y(2,185)= 3; index_Y(3,185)= 4; d_Y(185)= -6.7628530104625895E+01_ark
   index_Y(1,186)= 5; index_Y(2,186)= 3; index_Y(3,186)= 5; d_Y(186)=  4.3777997594856487E+01_ark
   index_Y(1,187)= 5; index_Y(2,187)= 3; index_Y(3,187)= 6; d_Y(187)= -1.6085223821280920E+01_ark
   index_Y(1,188)= 5; index_Y(2,188)= 3; index_Y(3,188)= 7; d_Y(188)=  3.1075077037893948E+00_ark
   index_Y(1,189)= 5; index_Y(2,189)= 3; index_Y(3,189)= 8; d_Y(189)= -2.4279535236144056E-01_ark
   index_Y(1,190)= 6; index_Y(2,190)= 1; index_Y(3,190)= 0; d_Y(190)= -2.1590405575511369E+01_ark
   index_Y(1,191)= 6; index_Y(2,191)= 1; index_Y(3,191)= 1; d_Y(191)=  1.5439736774080194E+02_ark
   index_Y(1,192)= 6; index_Y(2,192)= 1; index_Y(3,192)= 2; d_Y(192)= -4.5741720728322628E+02_ark
   index_Y(1,193)= 6; index_Y(2,193)= 1; index_Y(3,193)= 3; d_Y(193)=  7.2284776661163482E+02_ark
   index_Y(1,194)= 6; index_Y(2,194)= 1; index_Y(3,194)= 4; d_Y(194)= -6.6748336843829225E+02_ark
   index_Y(1,195)= 6; index_Y(2,195)= 1; index_Y(3,195)= 5; d_Y(195)=  3.7116165904689103E+02_ark
   index_Y(1,196)= 6; index_Y(2,196)= 1; index_Y(3,196)= 6; d_Y(196)= -1.2211303243164140E+02_ark
   index_Y(1,197)= 6; index_Y(2,197)= 1; index_Y(3,197)= 7; d_Y(197)=  2.1829555446015608E+01_ark
   index_Y(1,198)= 6; index_Y(2,198)= 1; index_Y(3,198)= 8; d_Y(198)= -1.6274867776484427E+00_ark
   index_Y(1,199)= 6; index_Y(2,199)= 3; index_Y(3,199)= 0; d_Y(199)= -4.8316213686971921E-03_ark
   index_Y(1,200)= 6; index_Y(2,200)= 3; index_Y(3,200)= 1; d_Y(200)= -1.1945850803870606E-01_ark
   index_Y(1,201)= 6; index_Y(2,201)= 3; index_Y(3,201)= 2; d_Y(201)=  7.1020395500935152E-01_ark
   index_Y(1,202)= 6; index_Y(2,202)= 3; index_Y(3,202)= 3; d_Y(202)= -1.5731098589709838E+00_ark
   index_Y(1,203)= 6; index_Y(2,203)= 3; index_Y(3,203)= 4; d_Y(203)=  1.7919405288298380E+00_ark
   index_Y(1,204)= 6; index_Y(2,204)= 3; index_Y(3,204)= 5; d_Y(204)= -1.1473089309379709E+00_ark
   index_Y(1,205)= 6; index_Y(2,205)= 3; index_Y(3,205)= 6; d_Y(205)=  4.1419385495150784E-01_ark
   index_Y(1,206)= 6; index_Y(2,206)= 3; index_Y(3,206)= 7; d_Y(206)= -7.7926667627686186E-02_ark
   index_Y(1,207)= 6; index_Y(2,207)= 3; index_Y(3,207)= 8; d_Y(207)=  5.8468889732239404E-03_ark
   index_Y(1,208)= 7; index_Y(2,208)= 1; index_Y(3,208)= 0; d_Y(208)=  1.3211435942230620E+00_ark
   index_Y(1,209)= 7; index_Y(2,209)= 1; index_Y(3,209)= 1; d_Y(209)= -9.5519935440447430E+00_ark
   index_Y(1,210)= 7; index_Y(2,210)= 1; index_Y(3,210)= 2; d_Y(210)=  2.8568224211879119E+01_ark
   index_Y(1,211)= 7; index_Y(2,211)= 1; index_Y(3,211)= 3; d_Y(211)= -4.5544096308605127E+01_ark
   index_Y(1,212)= 7; index_Y(2,212)= 1; index_Y(3,212)= 4; d_Y(212)=  4.2419900743339653E+01_ark
   index_Y(1,213)= 7; index_Y(2,213)= 1; index_Y(3,213)= 5; d_Y(213)= -2.3794376588111639E+01_ark
   index_Y(1,214)= 7; index_Y(2,214)= 1; index_Y(3,214)= 6; d_Y(214)=  7.8986365960530893E+00_ark
   index_Y(1,215)= 7; index_Y(2,215)= 1; index_Y(3,215)= 7; d_Y(215)= -1.4251448800616546E+00_ark
   index_Y(1,216)= 7; index_Y(2,216)= 1; index_Y(3,216)= 8; d_Y(216)=  1.0728935899324099E-01_ark
   index_Y(1,217)= 8; index_Y(2,217)= 1; index_Y(3,217)= 0; d_Y(217)= -3.4771940960526423E-02_ark
   index_Y(1,218)= 8; index_Y(2,218)= 1; index_Y(3,218)= 1; d_Y(218)=  2.5351582752788820E-01_ark
   index_Y(1,219)= 8; index_Y(2,219)= 1; index_Y(3,219)= 2; d_Y(219)= -7.6323975153311896E-01_ark
   index_Y(1,220)= 8; index_Y(2,220)= 1; index_Y(3,220)= 3; d_Y(220)=  1.2237046278640684E+00_ark
   index_Y(1,221)= 8; index_Y(2,221)= 1; index_Y(3,221)= 4; d_Y(221)= -1.1457244338162937E+00_ark
   index_Y(1,222)= 8; index_Y(2,222)= 1; index_Y(3,222)= 5; d_Y(222)=  6.4582562302644675E-01_ark
   index_Y(1,223)= 8; index_Y(2,223)= 1; index_Y(3,223)= 6; d_Y(223)= -2.1537317799830588E-01_ark
   index_Y(1,224)= 8; index_Y(2,224)= 1; index_Y(3,224)= 7; d_Y(224)=  3.9024421399244842E-02_ark
   index_Y(1,225)= 8; index_Y(2,225)= 1; index_Y(3,225)= 8; d_Y(225)= -2.9488510705136636E-03_ark

 elseif(dipole_surface == "LTP2011S") then
   n_r_X = 6; n_theta_X = 7; ncoeffs_X = 112
   index_X(1,001)= 0; index_X(2,001)= 0; index_X(3,001)= 1; d_X(001)= -2.4567218078481318E+01_ark
   index_X(1,002)= 0; index_X(2,002)= 0; index_X(3,002)= 2; d_X(002)= -6.8081841551683851E+01_ark
   index_X(1,003)= 0; index_X(2,003)= 0; index_X(3,003)= 3; d_X(003)=  3.2349775989087141E+02_ark
   index_X(1,004)= 0; index_X(2,004)= 0; index_X(3,004)= 4; d_X(004)= -4.1717036379506840E+02_ark
   index_X(1,005)= 0; index_X(2,005)= 0; index_X(3,005)= 5; d_X(005)=  2.4749541716930176E+02_ark
   index_X(1,006)= 0; index_X(2,006)= 0; index_X(3,006)= 6; d_X(006)= -7.0157419291643492E+01_ark
   index_X(1,007)= 0; index_X(2,007)= 0; index_X(3,007)= 7; d_X(007)=  7.6565447225719936E+00_ark
   index_X(1,008)= 0; index_X(2,008)= 2; index_X(3,008)= 1; d_X(008)= -2.2383455956573215E+01_ark
   index_X(1,009)= 0; index_X(2,009)= 2; index_X(3,009)= 2; d_X(009)=  7.3398704401935220E+01_ark
   index_X(1,010)= 0; index_X(2,010)= 2; index_X(3,010)= 3; d_X(010)= -1.0397217255066653E+02_ark
   index_X(1,011)= 0; index_X(2,011)= 2; index_X(3,011)= 4; d_X(011)=  7.8471081504467293E+01_ark
   index_X(1,012)= 0; index_X(2,012)= 2; index_X(3,012)= 5; d_X(012)= -3.2568972928225548E+01_ark
   index_X(1,013)= 0; index_X(2,013)= 2; index_X(3,013)= 6; d_X(013)=  6.9649748009962877E+00_ark
   index_X(1,014)= 0; index_X(2,014)= 2; index_X(3,014)= 7; d_X(014)= -5.9254409882299797E-01_ark
   index_X(1,015)= 0; index_X(2,015)= 4; index_X(3,015)= 1; d_X(015)=  5.4120721679885975E-01_ark
   index_X(1,016)= 0; index_X(2,016)= 4; index_X(3,016)= 2; d_X(016)= -3.9475806816103756E+00_ark
   index_X(1,017)= 0; index_X(2,017)= 4; index_X(3,017)= 3; d_X(017)=  8.7231958444806992E+00_ark
   index_X(1,018)= 0; index_X(2,018)= 4; index_X(3,018)= 4; d_X(018)= -9.0560761182926655E+00_ark
   index_X(1,019)= 0; index_X(2,019)= 4; index_X(3,019)= 5; d_X(019)=  4.9084244070241425E+00_ark
   index_X(1,020)= 0; index_X(2,020)= 4; index_X(3,020)= 6; d_X(020)= -1.3391890317410002E+00_ark
   index_X(1,021)= 0; index_X(2,021)= 4; index_X(3,021)= 7; d_X(021)=  1.4415545079216940E-01_ark
   index_X(1,022)= 0; index_X(2,022)= 6; index_X(3,022)= 1; d_X(022)=  5.9123475365137068E-03_ark
   index_X(1,023)= 0; index_X(2,023)= 6; index_X(3,023)= 2; d_X(023)= -3.5427183432823028E-02_ark
   index_X(1,024)= 0; index_X(2,024)= 6; index_X(3,024)= 3; d_X(024)=  8.0105101389342792E-02_ark
   index_X(1,025)= 0; index_X(2,025)= 6; index_X(3,025)= 4; d_X(025)= -8.8466983988197967E-02_ark
   index_X(1,026)= 0; index_X(2,026)= 6; index_X(3,026)= 5; d_X(026)=  5.1166426434406276E-02_ark
   index_X(1,027)= 0; index_X(2,027)= 6; index_X(3,027)= 6; d_X(027)= -1.4904184489040517E-02_ark
   index_X(1,028)= 0; index_X(2,028)= 6; index_X(3,028)= 7; d_X(028)=  1.7264235782761261E-03_ark
   index_X(1,029)= 1; index_X(2,029)= 0; index_X(3,029)= 1; d_X(029)=  2.6152238815727287E+01_ark
   index_X(1,030)= 1; index_X(2,030)= 0; index_X(3,030)= 2; d_X(030)=  1.4662299819969337E+02_ark
   index_X(1,031)= 1; index_X(2,031)= 0; index_X(3,031)= 3; d_X(031)= -5.5387752692025151E+02_ark
   index_X(1,032)= 1; index_X(2,032)= 0; index_X(3,032)= 4; d_X(032)=  6.8093025782439634E+02_ark
   index_X(1,033)= 1; index_X(2,033)= 0; index_X(3,033)= 5; d_X(033)= -3.9609466592676318E+02_ark
   index_X(1,034)= 1; index_X(2,034)= 0; index_X(3,034)= 6; d_X(034)=  1.1124491238514679E+02_ark
   index_X(1,035)= 1; index_X(2,035)= 0; index_X(3,035)= 7; d_X(035)= -1.2099952398229107E+01_ark
   index_X(1,036)= 1; index_X(2,036)= 2; index_X(3,036)= 1; d_X(036)=  2.0312935167604167E+01_ark
   index_X(1,037)= 1; index_X(2,037)= 2; index_X(3,037)= 2; d_X(037)= -6.8322111585084940E+01_ark
   index_X(1,038)= 1; index_X(2,038)= 2; index_X(3,038)= 3; d_X(038)=  9.9578930642296200E+01_ark
   index_X(1,039)= 1; index_X(2,039)= 2; index_X(3,039)= 4; d_X(039)= -7.7369241246095498E+01_ark
   index_X(1,040)= 1; index_X(2,040)= 2; index_X(3,040)= 5; d_X(040)=  3.3128894984283079E+01_ark
   index_X(1,041)= 1; index_X(2,041)= 2; index_X(3,041)= 6; d_X(041)= -7.3425373864523067E+00_ark
   index_X(1,042)= 1; index_X(2,042)= 2; index_X(3,042)= 7; d_X(042)=  6.5283291983126901E-01_ark
   index_X(1,043)= 1; index_X(2,043)= 4; index_X(3,043)= 1; d_X(043)= -1.9233997419253512E-01_ark
   index_X(1,044)= 1; index_X(2,044)= 4; index_X(3,044)= 2; d_X(044)=  1.3872924195991558E+00_ark
   index_X(1,045)= 1; index_X(2,045)= 4; index_X(3,045)= 3; d_X(045)= -3.0289054385799794E+00_ark
   index_X(1,046)= 1; index_X(2,046)= 4; index_X(3,046)= 4; d_X(046)=  3.1100011257216522E+00_ark
   index_X(1,047)= 1; index_X(2,047)= 4; index_X(3,047)= 5; d_X(047)= -1.6711459732326830E+00_ark
   index_X(1,048)= 1; index_X(2,048)= 4; index_X(3,048)= 6; d_X(048)=  4.5340447306500664E-01_ark
   index_X(1,049)= 1; index_X(2,049)= 4; index_X(3,049)= 7; d_X(049)= -4.8686366608285070E-02_ark
   index_X(1,050)= 2; index_X(2,050)= 0; index_X(3,050)= 1; d_X(050)= -8.7202636323057447E+00_ark
   index_X(1,051)= 2; index_X(2,051)= 0; index_X(3,051)= 2; d_X(051)= -1.1819340065481248E+02_ark
   index_X(1,052)= 2; index_X(2,052)= 0; index_X(3,052)= 3; d_X(052)=  3.8330228083227621E+02_ark
   index_X(1,053)= 2; index_X(2,053)= 0; index_X(3,053)= 4; d_X(053)= -4.5282827301551060E+02_ark
   index_X(1,054)= 2; index_X(2,054)= 0; index_X(3,054)= 5; d_X(054)=  2.5881959146085364E+02_ark
   index_X(1,055)= 2; index_X(2,055)= 0; index_X(3,055)= 6; d_X(055)= -7.2063689588218608E+01_ark
   index_X(1,056)= 2; index_X(2,056)= 0; index_X(3,056)= 7; d_X(056)=  7.8113661205319431E+00_ark
   index_X(1,057)= 2; index_X(2,057)= 2; index_X(3,057)= 1; d_X(057)= -6.7444812065591577E+00_ark
   index_X(1,058)= 2; index_X(2,058)= 2; index_X(3,058)= 2; d_X(058)=  2.3441940344599686E+01_ark
   index_X(1,059)= 2; index_X(2,059)= 2; index_X(3,059)= 3; d_X(059)= -3.5152069673157484E+01_ark
   index_X(1,060)= 2; index_X(2,060)= 2; index_X(3,060)= 4; d_X(060)=  2.7999799762657290E+01_ark
   index_X(1,061)= 2; index_X(2,061)= 2; index_X(3,061)= 5; d_X(061)= -1.2279871957605877E+01_ark
   index_X(1,062)= 2; index_X(2,062)= 2; index_X(3,062)= 6; d_X(062)=  2.7903642955715071E+00_ark
   index_X(1,063)= 2; index_X(2,063)= 2; index_X(3,063)= 7; d_X(063)= -2.5503652427731183E-01_ark
   index_X(1,064)= 2; index_X(2,064)= 4; index_X(3,064)= 1; d_X(064)=  1.3456878405657235E-02_ark
   index_X(1,065)= 2; index_X(2,065)= 4; index_X(3,065)= 2; d_X(065)= -1.0144749965029121E-01_ark
   index_X(1,066)= 2; index_X(2,066)= 4; index_X(3,066)= 3; d_X(066)=  2.1954806647110203E-01_ark
   index_X(1,067)= 2; index_X(2,067)= 4; index_X(3,067)= 4; d_X(067)= -2.2134321834120485E-01_ark
   index_X(1,068)= 2; index_X(2,068)= 4; index_X(3,068)= 5; d_X(068)=  1.1673648035752571E-01_ark
   index_X(1,069)= 2; index_X(2,069)= 4; index_X(3,069)= 6; d_X(069)= -3.1135480950513994E-02_ark
   index_X(1,070)= 2; index_X(2,070)= 4; index_X(3,070)= 7; d_X(070)=  3.2882383180670161E-03_ark
   index_X(1,071)= 3; index_X(2,071)= 0; index_X(3,071)= 1; d_X(071)=  3.5308229330414775E-01_ark
   index_X(1,072)= 3; index_X(2,072)= 0; index_X(3,072)= 2; d_X(072)=  4.7330977812492193E+01_ark
   index_X(1,073)= 3; index_X(2,073)= 0; index_X(3,073)= 3; d_X(073)= -1.3739283971621572E+02_ark
   index_X(1,074)= 3; index_X(2,074)= 0; index_X(3,074)= 4; d_X(074)=  1.5683143503968424E+02_ark
   index_X(1,075)= 3; index_X(2,075)= 0; index_X(3,075)= 5; d_X(075)= -8.8202079408237267E+01_ark
   index_X(1,076)= 3; index_X(2,076)= 0; index_X(3,076)= 6; d_X(076)=  2.4353783310326335E+01_ark
   index_X(1,077)= 3; index_X(2,077)= 0; index_X(3,077)= 7; d_X(077)= -2.6301154562547993E+00_ark
   index_X(1,078)= 3; index_X(2,078)= 2; index_X(3,078)= 1; d_X(078)=  9.6738762743321161E-01_ark
   index_X(1,079)= 3; index_X(2,079)= 2; index_X(3,079)= 2; d_X(079)= -3.4761672506523098E+00_ark
   index_X(1,080)= 3; index_X(2,080)= 2; index_X(3,080)= 3; d_X(080)=  5.3339286424944277E+00_ark
   index_X(1,081)= 3; index_X(2,081)= 2; index_X(3,081)= 4; d_X(081)= -4.3208285272395273E+00_ark
   index_X(1,082)= 3; index_X(2,082)= 2; index_X(3,082)= 5; d_X(082)=  1.9222953234855531E+00_ark
   index_X(1,083)= 3; index_X(2,083)= 2; index_X(3,083)= 6; d_X(083)= -4.4287395095560100E-01_ark
   index_X(1,084)= 3; index_X(2,084)= 2; index_X(3,084)= 7; d_X(084)=  4.1065113411043797E-02_ark
   index_X(1,085)= 4; index_X(2,085)= 0; index_X(3,085)= 1; d_X(085)=  3.9374566094527985E-01_ark
   index_X(1,086)= 4; index_X(2,086)= 0; index_X(3,086)= 2; d_X(086)= -1.0085949605580609E+01_ark
   index_X(1,087)= 4; index_X(2,087)= 0; index_X(3,087)= 3; d_X(087)=  2.6881267312433984E+01_ark
   index_X(1,088)= 4; index_X(2,088)= 0; index_X(3,088)= 4; d_X(088)= -2.9770383888124947E+01_ark
   index_X(1,089)= 4; index_X(2,089)= 0; index_X(3,089)= 5; d_X(089)=  1.6492123680104811E+01_ark
   index_X(1,090)= 4; index_X(2,090)= 0; index_X(3,090)= 6; d_X(090)= -4.5163955307452088E+00_ark
   index_X(1,091)= 4; index_X(2,091)= 0; index_X(3,091)= 7; d_X(091)=  4.8577911409212682E-01_ark
   index_X(1,092)= 4; index_X(2,092)= 2; index_X(3,092)= 1; d_X(092)= -5.0001737097788101E-02_ark
   index_X(1,093)= 4; index_X(2,093)= 2; index_X(3,093)= 2; d_X(093)=  1.8444314169333298E-01_ark
   index_X(1,094)= 4; index_X(2,094)= 2; index_X(3,094)= 3; d_X(094)= -2.8650654462613900E-01_ark
   index_X(1,095)= 4; index_X(2,095)= 2; index_X(3,095)= 4; d_X(095)=  2.3313551940647348E-01_ark
   index_X(1,096)= 4; index_X(2,096)= 2; index_X(3,096)= 5; d_X(096)= -1.0380031391543509E-01_ark
   index_X(1,097)= 4; index_X(2,097)= 2; index_X(3,097)= 6; d_X(097)=  2.3890228357732579E-02_ark
   index_X(1,098)= 4; index_X(2,098)= 2; index_X(3,098)= 7; d_X(098)= -2.2100705497907214E-03_ark
   index_X(1,099)= 5; index_X(2,099)= 0; index_X(3,099)= 1; d_X(099)= -8.0434126044985677E-02_ark
   index_X(1,100)= 5; index_X(2,100)= 0; index_X(3,100)= 2; d_X(100)=  1.0912947046464767E+00_ark
   index_X(1,101)= 5; index_X(2,101)= 0; index_X(3,101)= 3; d_X(101)= -2.7167160129831514E+00_ark
   index_X(1,102)= 5; index_X(2,102)= 0; index_X(3,102)= 4; d_X(102)=  2.9287490775294729E+00_ark
   index_X(1,103)= 5; index_X(2,103)= 0; index_X(3,103)= 5; d_X(103)= -1.5995393633554107E+00_ark
   index_X(1,104)= 5; index_X(2,104)= 0; index_X(3,104)= 6; d_X(104)=  4.3446245134876044E-01_ark
   index_X(1,105)= 5; index_X(2,105)= 0; index_X(3,105)= 7; d_X(105)= -4.6520766695406002E-02_ark
   index_X(1,106)= 6; index_X(2,106)= 0; index_X(3,106)= 1; d_X(106)=  4.6914695230123593E-03_ark
   index_X(1,107)= 6; index_X(2,107)= 0; index_X(3,107)= 2; d_X(107)= -4.6861310006580892E-02_ark
   index_X(1,108)= 6; index_X(2,108)= 0; index_X(3,108)= 3; d_X(108)=  1.1038537113234781E-01_ark
   index_X(1,109)= 6; index_X(2,109)= 0; index_X(3,109)= 4; d_X(109)= -1.1616548824510167E-01_ark
   index_X(1,110)= 6; index_X(2,110)= 0; index_X(3,110)= 5; d_X(110)=  6.2588666551305502E-02_ark
   index_X(1,111)= 6; index_X(2,111)= 0; index_X(3,111)= 6; d_X(111)= -1.6859109782274725E-02_ark
   index_X(1,112)= 6; index_X(2,112)= 0; index_X(3,112)= 7; d_X(112)=  1.7960415093437735E-03_ark
!********************************************************************************************
   n_r_Y = 7; n_theta_Y = 7; ncoeffs_Y = 128
   index_Y(1,001)= 0; index_Y(2,001)= 1; index_Y(3,001)= 0; d_Y(001)=  2.7602888242402298E+02_ark
   index_Y(1,002)= 0; index_Y(2,002)= 1; index_Y(3,002)= 1; d_Y(002)= -9.4429499978529589E+02_ark
   index_Y(1,003)= 0; index_Y(2,003)= 1; index_Y(3,003)= 2; d_Y(003)=  7.5785885374558939E+02_ark
   index_Y(1,004)= 0; index_Y(2,004)= 1; index_Y(3,004)= 3; d_Y(004)=  8.8156305097691984E+02_ark
   index_Y(1,005)= 0; index_Y(2,005)= 1; index_Y(3,005)= 4; d_Y(005)= -1.9838570828139186E+03_ark
   index_Y(1,006)= 0; index_Y(2,006)= 1; index_Y(3,006)= 5; d_Y(006)=  1.4054991689346207E+03_ark
   index_Y(1,007)= 0; index_Y(2,007)= 1; index_Y(3,007)= 6; d_Y(007)= -4.4749368889416564E+02_ark
   index_Y(1,008)= 0; index_Y(2,008)= 1; index_Y(3,008)= 7; d_Y(008)=  5.4172767011290489E+01_ark
   index_Y(1,009)= 0; index_Y(2,009)= 3; index_Y(3,009)= 0; d_Y(009)=  2.2329060544254503E+01_ark
   index_Y(1,010)= 0; index_Y(2,010)= 3; index_Y(3,010)= 1; d_Y(010)= -7.8109258549579408E+01_ark
   index_Y(1,011)= 0; index_Y(2,011)= 3; index_Y(3,011)= 2; d_Y(011)=  1.1514715303757930E+02_ark
   index_Y(1,012)= 0; index_Y(2,012)= 3; index_Y(3,012)= 3; d_Y(012)= -9.7032672324743544E+01_ark
   index_Y(1,013)= 0; index_Y(2,013)= 3; index_Y(3,013)= 4; d_Y(013)=  6.0183351071978450E+01_ark
   index_Y(1,014)= 0; index_Y(2,014)= 3; index_Y(3,014)= 5; d_Y(014)= -3.0681444120528568E+01_ark
   index_Y(1,015)= 0; index_Y(2,015)= 3; index_Y(3,015)= 6; d_Y(015)=  1.0253415663006308E+01_ark
   index_Y(1,016)= 0; index_Y(2,016)= 3; index_Y(3,016)= 7; d_Y(016)= -1.4642327032256617E+00_ark
   index_Y(1,017)= 0; index_Y(2,017)= 5; index_Y(3,017)= 0; d_Y(017)=  8.1014552535220474E-01_ark
   index_Y(1,018)= 0; index_Y(2,018)= 5; index_Y(3,018)= 1; d_Y(018)= -4.2000929742390980E-01_ark
   index_Y(1,019)= 0; index_Y(2,019)= 5; index_Y(3,019)= 2; d_Y(019)= -5.6471175480941724E+00_ark
   index_Y(1,020)= 0; index_Y(2,020)= 5; index_Y(3,020)= 3; d_Y(020)=  1.2919965178073880E+01_ark
   index_Y(1,021)= 0; index_Y(2,021)= 5; index_Y(3,021)= 4; d_Y(021)= -1.2308209235325876E+01_ark
   index_Y(1,022)= 0; index_Y(2,022)= 5; index_Y(3,022)= 5; d_Y(022)=  5.9820589387793461E+00_ark
   index_Y(1,023)= 0; index_Y(2,023)= 5; index_Y(3,023)= 6; d_Y(023)= -1.4555087882786211E+00_ark
   index_Y(1,024)= 0; index_Y(2,024)= 5; index_Y(3,024)= 7; d_Y(024)=  1.3984548365829141E-01_ark
   index_Y(1,025)= 0; index_Y(2,025)= 7; index_Y(3,025)= 0; d_Y(025)=  3.2050668156529127E-03_ark
   index_Y(1,026)= 0; index_Y(2,026)= 7; index_Y(3,026)= 1; d_Y(026)= -4.0742501596042757E-02_ark
   index_Y(1,027)= 0; index_Y(2,027)= 7; index_Y(3,027)= 2; d_Y(027)=  1.6314399467479035E-01_ark
   index_Y(1,028)= 0; index_Y(2,028)= 7; index_Y(3,028)= 3; d_Y(028)= -2.9705557000488980E-01_ark
   index_Y(1,029)= 0; index_Y(2,029)= 7; index_Y(3,029)= 4; d_Y(029)=  2.8313139602997239E-01_ark
   index_Y(1,030)= 0; index_Y(2,030)= 7; index_Y(3,030)= 5; d_Y(030)= -1.4638994350268320E-01_ark
   index_Y(1,031)= 0; index_Y(2,031)= 7; index_Y(3,031)= 6; d_Y(031)=  3.8860165050678575E-02_ark
   index_Y(1,032)= 0; index_Y(2,032)= 7; index_Y(3,032)= 7; d_Y(032)= -4.1464398149173576E-03_ark
   index_Y(1,033)= 1; index_Y(2,033)= 1; index_Y(3,033)= 0; d_Y(033)= -3.8432139472588131E+02_ark
   index_Y(1,034)= 1; index_Y(2,034)= 1; index_Y(3,034)= 1; d_Y(034)=  1.3252809675262451E+03_ark
   index_Y(1,035)= 1; index_Y(2,035)= 1; index_Y(3,035)= 2; d_Y(035)= -1.1522307006858700E+03_ark
   index_Y(1,036)= 1; index_Y(2,036)= 1; index_Y(3,036)= 3; d_Y(036)= -9.8871909194310047E+02_ark
   index_Y(1,037)= 1; index_Y(2,037)= 1; index_Y(3,037)= 4; d_Y(037)=  2.4938766083963237E+03_ark
   index_Y(1,038)= 1; index_Y(2,038)= 1; index_Y(3,038)= 5; d_Y(038)= -1.8015979357545530E+03_ark
   index_Y(1,039)= 1; index_Y(2,039)= 1; index_Y(3,039)= 6; d_Y(039)=  5.7806771717135962E+02_ark
   index_Y(1,040)= 1; index_Y(2,040)= 1; index_Y(3,040)= 7; d_Y(040)= -7.0274281915151732E+01_ark
   index_Y(1,041)= 1; index_Y(2,041)= 3; index_Y(3,041)= 0; d_Y(041)= -1.8047557525087541E+01_ark
   index_Y(1,042)= 1; index_Y(2,042)= 3; index_Y(3,042)= 1; d_Y(042)=  6.1739597237700309E+01_ark
   index_Y(1,043)= 1; index_Y(2,043)= 3; index_Y(3,043)= 2; d_Y(043)= -9.0640390693787595E+01_ark
   index_Y(1,044)= 1; index_Y(2,044)= 3; index_Y(3,044)= 3; d_Y(044)=  7.7728417004245784E+01_ark
   index_Y(1,045)= 1; index_Y(2,045)= 3; index_Y(3,045)= 4; d_Y(045)= -4.9694202835600208E+01_ark
   index_Y(1,046)= 1; index_Y(2,046)= 3; index_Y(3,046)= 5; d_Y(046)=  2.5587862360683744E+01_ark
   index_Y(1,047)= 1; index_Y(2,047)= 3; index_Y(3,047)= 6; d_Y(047)= -8.4413526546265985E+00_ark
   index_Y(1,048)= 1; index_Y(2,048)= 3; index_Y(3,048)= 7; d_Y(048)=  1.1838772024136688E+00_ark
   index_Y(1,049)= 1; index_Y(2,049)= 5; index_Y(3,049)= 0; d_Y(049)= -3.2520439831553460E-01_ark
   index_Y(1,050)= 1; index_Y(2,050)= 5; index_Y(3,050)= 1; d_Y(050)=  3.4128648276201190E-01_ark
   index_Y(1,051)= 1; index_Y(2,051)= 5; index_Y(3,051)= 2; d_Y(051)=  1.6143398924409098E+00_ark
   index_Y(1,052)= 1; index_Y(2,052)= 5; index_Y(3,052)= 3; d_Y(052)= -4.2052943666292322E+00_ark
   index_Y(1,053)= 1; index_Y(2,053)= 5; index_Y(3,053)= 4; d_Y(053)=  4.1886212296860776E+00_ark
   index_Y(1,054)= 1; index_Y(2,054)= 5; index_Y(3,054)= 5; d_Y(054)= -2.0902968408237825E+00_ark
   index_Y(1,055)= 1; index_Y(2,055)= 5; index_Y(3,055)= 6; d_Y(055)=  5.1846653005009280E-01_ark
   index_Y(1,056)= 1; index_Y(2,056)= 5; index_Y(3,056)= 7; d_Y(056)= -5.0611720629774481E-02_ark
   index_Y(1,057)= 2; index_Y(2,057)= 1; index_Y(3,057)= 0; d_Y(057)=  2.1918743496194492E+02_ark
   index_Y(1,058)= 2; index_Y(2,058)= 1; index_Y(3,058)= 1; d_Y(058)= -7.6450372323722149E+02_ark
   index_Y(1,059)= 2; index_Y(2,059)= 1; index_Y(3,059)= 2; d_Y(059)=  7.2080653174528652E+02_ark
   index_Y(1,060)= 2; index_Y(2,060)= 1; index_Y(3,060)= 3; d_Y(060)=  4.1797734325992633E+02_ark
   index_Y(1,061)= 2; index_Y(2,061)= 1; index_Y(3,061)= 4; d_Y(061)= -1.2646969514690613E+03_ark
   index_Y(1,062)= 2; index_Y(2,062)= 1; index_Y(3,062)= 5; d_Y(062)=  9.3869378688125698E+02_ark
   index_Y(1,063)= 2; index_Y(2,063)= 1; index_Y(3,063)= 6; d_Y(063)= -3.0447047410884761E+02_ark
   index_Y(1,064)= 2; index_Y(2,064)= 1; index_Y(3,064)= 7; d_Y(064)=  3.7234428779570862E+01_ark
   index_Y(1,065)= 2; index_Y(2,065)= 3; index_Y(3,065)= 0; d_Y(065)=  5.2409269022925287E+00_ark
   index_Y(1,066)= 2; index_Y(2,066)= 3; index_Y(3,066)= 1; d_Y(066)= -1.7804244461312010E+01_ark
   index_Y(1,067)= 2; index_Y(2,067)= 3; index_Y(3,067)= 2; d_Y(067)=  2.6975177591892134E+01_ark
   index_Y(1,068)= 2; index_Y(2,068)= 3; index_Y(3,068)= 3; d_Y(068)= -2.4969534319308423E+01_ark
   index_Y(1,069)= 2; index_Y(2,069)= 3; index_Y(3,069)= 4; d_Y(069)=  1.7293110873518572E+01_ark
   index_Y(1,070)= 2; index_Y(2,070)= 3; index_Y(3,070)= 5; d_Y(070)= -9.0293501453666067E+00_ark
   index_Y(1,071)= 2; index_Y(2,071)= 3; index_Y(3,071)= 6; d_Y(071)=  2.8753674384406622E+00_ark
   index_Y(1,072)= 2; index_Y(2,072)= 3; index_Y(3,072)= 7; d_Y(072)= -3.8635910001269735E-01_ark
   index_Y(1,073)= 2; index_Y(2,073)= 5; index_Y(3,073)= 0; d_Y(073)=  2.9865528248732787E-02_ark
   index_Y(1,074)= 2; index_Y(2,074)= 5; index_Y(3,074)= 1; d_Y(074)= -2.5769277141414193E-02_ark
   index_Y(1,075)= 2; index_Y(2,075)= 5; index_Y(3,075)= 2; d_Y(075)= -1.8537093127383741E-01_ark
   index_Y(1,076)= 2; index_Y(2,076)= 5; index_Y(3,076)= 3; d_Y(076)=  4.7097890494002570E-01_ark
   index_Y(1,077)= 2; index_Y(2,077)= 5; index_Y(3,077)= 4; d_Y(077)= -4.7650923576138871E-01_ark
   index_Y(1,078)= 2; index_Y(2,078)= 5; index_Y(3,078)= 5; d_Y(078)=  2.4343213320304002E-01_ark
   index_Y(1,079)= 2; index_Y(2,079)= 5; index_Y(3,079)= 6; d_Y(079)= -6.2034057520264696E-02_ark
   index_Y(1,080)= 2; index_Y(2,080)= 5; index_Y(3,080)= 7; d_Y(080)=  6.2472338479864220E-03_ark
   index_Y(1,081)= 3; index_Y(2,081)= 1; index_Y(3,081)= 0; d_Y(081)= -6.5551885257176508E+01_ark
   index_Y(1,082)= 3; index_Y(2,082)= 1; index_Y(3,082)= 1; d_Y(082)=  2.3187666974386809E+02_ark
   index_Y(1,083)= 3; index_Y(2,083)= 1; index_Y(3,083)= 2; d_Y(083)= -2.3709872493424234E+02_ark
   index_Y(1,084)= 3; index_Y(2,084)= 1; index_Y(3,084)= 3; d_Y(084)= -7.8063830361784539E+01_ark
   index_Y(1,085)= 3; index_Y(2,085)= 1; index_Y(3,085)= 4; d_Y(085)=  3.2918255217036744E+02_ark
   index_Y(1,086)= 3; index_Y(2,086)= 1; index_Y(3,086)= 5; d_Y(086)= -2.5375797594657718E+02_ark
   index_Y(1,087)= 3; index_Y(2,087)= 1; index_Y(3,087)= 6; d_Y(087)=  8.3533330692747825E+01_ark
   index_Y(1,088)= 3; index_Y(2,088)= 1; index_Y(3,088)= 7; d_Y(088)= -1.0297944723905573E+01_ark
   index_Y(1,089)= 3; index_Y(2,089)= 3; index_Y(3,089)= 0; d_Y(089)= -6.4468105491265248E-01_ark
   index_Y(1,090)= 3; index_Y(2,090)= 3; index_Y(3,090)= 1; d_Y(090)=  2.2283502330610077E+00_ark
   index_Y(1,091)= 3; index_Y(2,091)= 3; index_Y(3,091)= 2; d_Y(091)= -3.6511865959138845E+00_ark
   index_Y(1,092)= 3; index_Y(2,092)= 3; index_Y(3,092)= 3; d_Y(092)=  3.8476781836540681E+00_ark
   index_Y(1,093)= 3; index_Y(2,093)= 3; index_Y(3,093)= 4; d_Y(093)= -2.9520174250588411E+00_ark
   index_Y(1,094)= 3; index_Y(2,094)= 3; index_Y(3,094)= 5; d_Y(094)=  1.5593226687926816E+00_ark
   index_Y(1,095)= 3; index_Y(2,095)= 3; index_Y(3,095)= 6; d_Y(095)= -4.7499900867517653E-01_ark
   index_Y(1,096)= 3; index_Y(2,096)= 3; index_Y(3,096)= 7; d_Y(096)=  6.0407181906811047E-02_ark
   index_Y(1,097)= 4; index_Y(2,097)= 1; index_Y(3,097)= 0; d_Y(097)=  1.0825011045884462E+01_ark
   index_Y(1,098)= 4; index_Y(2,098)= 1; index_Y(3,098)= 1; d_Y(098)= -3.8957767142692830E+01_ark
   index_Y(1,099)= 4; index_Y(2,099)= 1; index_Y(3,099)= 2; d_Y(099)=  4.3089977185643846E+01_ark
   index_Y(1,100)= 4; index_Y(2,100)= 1; index_Y(3,100)= 3; d_Y(100)=  4.7803864478065918E+00_ark
   index_Y(1,101)= 4; index_Y(2,101)= 1; index_Y(3,101)= 4; d_Y(101)= -4.6183571977701945E+01_ark
   index_Y(1,102)= 4; index_Y(2,102)= 1; index_Y(3,102)= 5; d_Y(102)=  3.7519332897891758E+01_ark
   index_Y(1,103)= 4; index_Y(2,103)= 1; index_Y(3,103)= 6; d_Y(103)= -1.2593171155996748E+01_ark
   index_Y(1,104)= 4; index_Y(2,104)= 1; index_Y(3,104)= 7; d_Y(104)=  1.5685499143757138E+00_ark
   index_Y(1,105)= 4; index_Y(2,105)= 3; index_Y(3,105)= 0; d_Y(105)=  2.8438540722877548E-02_ark
   index_Y(1,106)= 4; index_Y(2,106)= 3; index_Y(3,106)= 1; d_Y(106)= -1.0501497222346945E-01_ark
   index_Y(1,107)= 4; index_Y(2,107)= 3; index_Y(3,107)= 2; d_Y(107)=  1.9951569143079340E-01_ark
   index_Y(1,108)= 4; index_Y(2,108)= 3; index_Y(3,108)= 3; d_Y(108)= -2.4888868889716687E-01_ark
   index_Y(1,109)= 4; index_Y(2,109)= 3; index_Y(3,109)= 4; d_Y(109)=  2.0999196974309203E-01_ark
   index_Y(1,110)= 4; index_Y(2,110)= 3; index_Y(3,110)= 5; d_Y(110)= -1.1118952591081355E-01_ark
   index_Y(1,111)= 4; index_Y(2,111)= 3; index_Y(3,111)= 6; d_Y(111)=  3.2351835799717010E-02_ark
   index_Y(1,112)= 4; index_Y(2,112)= 3; index_Y(3,112)= 7; d_Y(112)= -3.8833615745925698E-03_ark
   index_Y(1,113)= 5; index_Y(2,113)= 1; index_Y(3,113)= 0; d_Y(113)= -9.3484381438870567E-01_ark
   index_Y(1,114)= 5; index_Y(2,114)= 1; index_Y(3,114)= 1; d_Y(114)=  3.4317881852124827E+00_ark
   index_Y(1,115)= 5; index_Y(2,115)= 1; index_Y(3,115)= 2; d_Y(115)= -4.0764022208982631E+00_ark
   index_Y(1,116)= 5; index_Y(2,116)= 1; index_Y(3,116)= 3; d_Y(116)=  2.7838139813782803E-01_ark
   index_Y(1,117)= 5; index_Y(2,117)= 1; index_Y(3,117)= 4; d_Y(117)=  3.3163816975473530E+00_ark
   index_Y(1,118)= 5; index_Y(2,118)= 1; index_Y(3,118)= 5; d_Y(118)= -2.8886400811141786E+00_ark
   index_Y(1,119)= 5; index_Y(2,119)= 1; index_Y(3,119)= 6; d_Y(119)=  9.9293453246556818E-01_ark
   index_Y(1,120)= 5; index_Y(2,120)= 1; index_Y(3,120)= 7; d_Y(120)= -1.2518679040823610E-01_ark
   index_Y(1,121)= 6; index_Y(2,121)= 1; index_Y(3,121)= 0; d_Y(121)=  3.2976694357696651E-02_ark
   index_Y(1,122)= 6; index_Y(2,122)= 1; index_Y(3,122)= 1; d_Y(122)= -1.2349687556161909E-01_ark
   index_Y(1,123)= 6; index_Y(2,123)= 1; index_Y(3,123)= 2; d_Y(123)=  1.5546881716498062E-01_ark
   index_Y(1,124)= 6; index_Y(2,124)= 1; index_Y(3,124)= 3; d_Y(124)= -3.1329479442763097E-02_ark
   index_Y(1,125)= 6; index_Y(2,125)= 1; index_Y(3,125)= 4; d_Y(125)= -9.6818120967887611E-02_ark
   index_Y(1,126)= 6; index_Y(2,126)= 1; index_Y(3,126)= 5; d_Y(126)=  9.1606357682491368E-02_ark
   index_Y(1,127)= 6; index_Y(2,127)= 1; index_Y(3,127)= 6; d_Y(127)= -3.2301323954391933E-02_ark
   index_Y(1,128)= 6; index_Y(2,128)= 1; index_Y(3,128)= 7; d_Y(128)=  4.1230694338540395E-03_ark
 endif

  !zLoaded = .true.
endif

!***********
! Check that input makes sense
if( abs(cos_theta) > 1._ark) write(*,*) 'WARNING: |cos(theta)| > 1 in DIPS'
if( r1 < 0._ark)             write(*,*) 'WARNING: r1 < 0 in DIPS'
if( r2 < 0._ark)             write(*,*) 'WARNING: r2 < 0 in DIPS'


!***************
! FITTING VARIABLES
theta = acos(cos_theta)
!***************
x1 = r1 + r2
x2 = r2 - r1
x3 = PI - theta
!***************

!***************
! PRE-COMPUTE POWERS OF THE VARIABLES
powers_x1(0) = 1._ark
powers_x2(0) = 1._ark
powers_x3(0) = 1._ark

do ii=1, max(n_r_X, n_r_Y)
 powers_x1(ii) =  x1*powers_x1(ii-1)
 powers_x2(ii) =  x2*powers_x2(ii-1)
enddo

do ii=1, max(n_theta_X, n_theta_Y)
 powers_x3(ii) =  x3*powers_x3(ii-1)
enddo

!***************
muX = 0._ark
do ii = 1, ncoeffs_X
  muX = muX + d_X(ii) * powers_x1(index_X(1,ii)) * powers_x2(index_X(2,ii)) * powers_x3(index_X(3,ii))
enddo
!***************
muY = 0._ark
do ii = 1, ncoeffs_Y
  muY = muY + d_Y(ii) * powers_x1(index_Y(1,ii)) * powers_x2(index_Y(2,ii)) * powers_x3(index_Y(3,ii))
enddo
!***************

!add damping/asymptotic terms
damp1 = damp(r1)
damp2 = damp(r2)
value1 = muOH(r1)*( 1._ark - damp2 )
value2 = muOH(r2)*( 1._ark - damp1 )

cos_theta_half = sqrt( (1._ark + cos_theta)/2._ark  )
sin_theta_half = sqrt( (1._ark - cos_theta)/2._ark  )

muX = muX*damp1*damp2 + cos_theta_half*( value1 + value2  )
muY = muY*damp1*damp2 + sin_theta_half*( value1 - value2  )

end subroutine DIPS_LPT2011

!==================================
function damp(r)
! Damping function.
! It is exactly =0 for r< rMIN;
! It is exactly =1 for r> rMAX;
IMPLICIT NONE
!integer, parameter :: fp =kind(1.d0)
real(ark), parameter :: PI=3.14159265358979323846264338328_ark
real(ark), intent(in) :: r
real(ark) :: damp
real(ark), parameter :: rMIN=3.5_ark, rMAX=5.5_ark
real(ark) :: t

t = ( ( (r-rMIN)/(rMAX-rMIN) ) -0.5_ark )

if( t <= -0.5_ark) then
   damp = 1._ark
   return
elseif(t >= 0.5_ark) then
   damp = 0._ark
   return
endif

damp = 1._ark / ( EXP(4._ark*tan(PI*t)) + 1._ark )

end function damp
!==================================
function muOH(r)
IMPLICIT NONE
! Lorenzo Lodi 7 April 2011
! It's a fit in even-tempered Gaussian functions of 59 AQCC[7,5]/XP aug-cc-pV5Z dipoles for the ground X state of OH
! between 1 and 4 bohrs. Dipoles should be accurate to about 2.e-3. Relativistic & core-corrections are NOT included.
! NB: There is a NASTY cancellation error with this fit which eats up about 5 digits of accuracy.
! double precision reals are compulsory. RMS of the fit is 5e-5 a.u. The largest residuals are about 0.7e-3 a.u.
! Note: Kahan summation algorithm (compensated summation) does not help here. The problem is ill-conditioned
!       as sum_i |x_i| / |sum_i x_i| = ~ 3e+7
! It might be solved by using as basis orthogonalised Gaussians, but I didn't deem it necessary at this time.
!
! Fit has reasonable behaviour for r-> +oo and also not too bad for r->0.
!integer, parameter :: fp = ark
real(ark), intent(in) :: r
real(ark) :: muOH
integer, parameter :: n_coeffs = 8
real(ark), parameter :: a = 0.934_ark
real(ark) :: coeffs(n_coeffs)
real(ark) :: powers_a(n_coeffs)
real(ark) :: r2
integer :: i
logical :: zFirstRun  = .true.
!save zFirstRun, coeffs, powers_a

if(zFirstRun .eqv. .true.) then !Load fit coefficients

 powers_a(1) = a
 do i=2, n_coeffs
   powers_a(i) = powers_a(i-1)*a
 enddo

 coeffs( 1) =   -51632.29191732509_ark
 coeffs( 2) =   336439.76799518_ark
 coeffs( 3) =  -948384.0949923932_ark
 coeffs( 4) =  1499393.4799029354_ark
 coeffs( 5) = -1436572.608839631_ark
 coeffs( 6) =   834824.6196200893_ark
 coeffs( 7) =  -272815.3271033254_ark
 coeffs( 8) =    38746.11983311285_ark

 !zFirstRun=.false.

endif

r2 = r**2
muOH = 0._ark
do i=1, n_coeffs
  muOH = muOH + coeffs(i)*exp( -powers_a(i)*r2 )
enddo

muOH = muOH*r2

end function muOH

recursive subroutine MLdipole_h2o_lpt2011(rank,ncoords,natoms,local,xyz,f)

    integer(ik),intent(in) ::  rank,ncoords,natoms
    real(ark),intent(in)   ::  local(ncoords),xyz(natoms,3)
    real(ark),intent(out)  ::  f(rank)
    !
    integer(ik)           :: k,imu,iterm
    real(ark)             :: y(3), mu(3),u1(3),u2(3),u3(3),tmat(3,3),n1(3),n2(3),x(2,3),r1,r2,alpha,re,ae
    real(ark)             :: xyz0(natoms,3),xi(3),mu_t,cos_theta
    !
    ! xyz are undefined for the local case
    if (all(abs(xyz)<small_)) then 
      !
      select case(trim(molec%coords_transform))
      case default
         write (out,"('MLdipole_h2o_lpt2011: coord. type ',a,' unknown')") trim(molec%coords_transform)
         stop 'MLdipole_h2o_lpt2011 - bad coord. type'
      case('R-RHO-Z')
         !
         xyz0 = MLloc2pqr_xy2(local)
         !
         x(1,:) = xyz0(2,:) - xyz0(1,:)
         x(2,:) = xyz0(3,:) - xyz0(1,:)
         !
      end select
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
    n1 = x(1,:) / r1
    n2 = x(2,:) / r2
    !
    alpha = acos(sum(n1*n2))
    !
    cos_theta = sum(n1*n2)
    !
    call DIPS_LPT2011(mu(2), mu(1), r1, r2, cos_theta)
    !
    mu(3) = 0
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
    f(1:3) = matmul(mu,tmat)
    !
    f = f/0.393430_ark ! to debye
    !
 end subroutine MLdipole_h2o_lpt2011



end module prop_xy2