!
!  This unit defines all specific routines for a triatomic molecule of XY2 type
!
module prop_xy2
  use accuracy
  use moltype
  use timer

  implicit none

  public prop_xy2_qmom_sym

  private
 
  integer(ik), parameter :: verbose     = 4                          ! Verbosity level


  contains


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

  !tmat_inv = tmat
  !call matrix_pseudoinverse_ark( 3, 3, tmat(1:3,1:3), tmat_inv(1:3,1:3) )
  !
  call MLinvmatark(tmat,tmat_inv,3,ierr)
  !
  if (ierr/=0) then
    write(out,"('xy2_dipole_sym error: failed inverse',i0)") ierr
    stop 'xy2_dipole_sym error: failed inverse'
  endif
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

  stop 'matrix_pseudoinverse_ark needs to be not activated'
  !!call matrix_pseudoinverse_ark( 3, 3, tmat(1:3,1:3), tmat_inv(1:3,1:3) )

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

  integer(ik),intent(in) ::  rank, ncoords, natoms
  real(ark),intent(in)   ::  local(ncoords), xyz(natoms,3)
  real(ark),intent(out)  ::  f(rank)

  integer(ik) :: iatom
  real(ark) :: xyz0(3), xyz_(natoms,3), r1, r2, alpha1, e1(3), e2(3), n1(3), n2(3), n3(3), tmat(3,3), &!
               coords(3), qmom_mb(3,3), qmom_xyz(3,3), qmom_xyz_(3,3), tmat_inv(3,3)

  if (rank/=6) then
    write(out, '(/a,1x,i3,1x,a)') 'xy2_qmom_sym error: rank of the dipole moment vector =', rank, ', expected 6'
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

  qmom_mb = 0
  qmom_mb(1,1) = fit_xy2_dipole_A1(extF%nterms(1), extF%coef(1:extF%nterms(1),1), coords)
  qmom_mb(2,2) = fit_xy2_dipole_A1(extF%nterms(2), extF%coef(1:extF%nterms(2),2), coords)
  qmom_mb(1,2) = fit_xy2_dipole_B2(extF%nterms(3), extF%coef(1:extF%nterms(3),3), coords)

  qmom_mb(3,3) = -( qmom_mb(1,1) + qmom_mb(2,2) )
  qmom_mb(2,1) = qmom_mb(1,2)

  stop 'matrix_pseudoinverse_ark needs to be not activated'
  !call matrix_pseudoinverse_ark( 3, 3, tmat(1:3,1:3), tmat_inv(1:3,1:3) )

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
    call xy2_qmom_sym(rank, ncoords, natoms, local, xyz_, mu_xyz(1:6))
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

  stop 'matrix_pseudoinverse_ark needs to be not activated'
  !call matrix_pseudoinverse_ark( 3, 3, tmat(1:3,1:3), tmat_inv(1:3,1:3) )

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


end module prop_xy2