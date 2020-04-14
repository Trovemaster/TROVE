module extfield
#define  _EMFIELD2_DEBUG_ 1

use accuracy
use timer
use rotme_cart_tens
use richmol_data
use fields, only: job
use moltype, only: intensity, molec, extF
use tran, only: TReigenvec_unit, bset_contr, read_contrind, read_eigenval, index_correlation, eigen, &
    Neigenlevels
use symmetry, only: sym
implicit none


real(rk), allocatable :: extf_vib_me(:,:,:)


contains


subroutine emf2_matelem

  type(rotme_cart_tens_type) :: tens
  integer(ik) :: jmin, jmax, dj, nJ, j, jind, info, oper_ielem
  integer(ik), allocatable :: Jval(:)
  real(rk) :: coef_tol, print_tol, linestr_tol, intens_tol, leading_coef_tol
  character(cl) :: oper, sielem

  oper = intensity%tens_oper
  oper_ielem = intensity%tens_oper_ielem
  jmin = minval(intensity%J(1:2))
  jmax = maxval(intensity%J(1:2))

  coef_tol = intensity%threshold%coeff
  linestr_tol = intensity%threshold%linestrength
  intens_tol = intensity%threshold%intensity
  print_tol = linestr_tol
  leading_coef_tol = intensity%threshold%leading_coeff

  if (sign(1.0_rk,coef_tol)<0) coef_tol = 1.0d-14
  if (sign(1.0_rk,linestr_tol)<0) linestr_tol = 1.0d-14
  if (sign(1.0_rk,intens_tol)<0) intens_tol = 1.0d-14
  if (sign(1.0_rk,print_tol)<0) print_tol = 1.0d-14
  if (sign(1.0_rk,leading_coef_tol)<0) leading_coef_tol = 1.0d-03


  Jmin = minval(intensity%J(1:2))
  Jmax = maxval(intensity%J(1:2))
  nJ   = Jmax - Jmin + 1

  if (Jmin>0) nJ = nJ + 1
  allocate(Jval(nJ), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'emfield2/emf2_matelem error: failed to allocate Jval(nJ)', 'nJ =', nJ
    stop 'STOP, error in emfield2/emf2_matelem'
  endif

  Jval = 0
  jind = 1
  do j=max(Jmin,1), Jmax
    jind = jind + 1
    Jval(jind) = j
  end do


  call read_contrind(nJ, Jval(1:nJ))
  call index_correlation(nJ, Jval(1:nJ))
  call read_eigenval(nJ, Jval(1:nJ))


  select case(trim(oper))

  case('VZZ')
    tens%func => rotme_vzz_trace0
    dj = 2
    call tens%init(jmin, jmax, dj, verbose=.true.)

    call read_extf_vib_me(tens%nelem)

    call rovib_me_storeall(tens, nJ, Jval, coef_tol, print_tol, leading_coef_tol)

  case('ALPHA')
    tens%func => rotme_alpha
    dj = 2
    call tens%init(jmin, jmax, dj, verbose=.true.)

    call read_extf_vib_me(tens%nelem)

    call rovib_me_storeall(tens, nJ, Jval, coef_tol, print_tol, leading_coef_tol)

  case('BETA')
    tens%func => rotme_beta
    dj = 3
    call tens%init(jmin, jmax, dj, verbose=.true.)

    call read_extf_vib_me(tens%nelem)

    call rovib_me_storeall(tens, nJ, Jval, coef_tol, print_tol, leading_coef_tol)

  case('MU')
    tens%func => rotme_mu
    dj = 1
    call tens%init(jmin, jmax, dj, verbose=.true.)

    call read_extf_vib_me(tens%nelem)

    call rovib_me_storeall(tens, nJ, Jval, coef_tol, print_tol, leading_coef_tol)

  case('COSTHETA')
    tens%func => rotme_costheta
    dj = 1
    call tens%init(jmin, jmax, dj, verbose=.true.)

    call init_extf_vib_me_overlap(tens%nelem)

    call rovib_me_storeall(tens, nJ, Jval, coef_tol, print_tol, leading_coef_tol)

  case('J')
    tens%func => rotme_j
    dj = 10
    call tens%init(jmin, jmax, dj, verbose=.true.)

    call init_extf_vib_me_overlap(tens%nelem)

    call rovib_me_storeall(tens, nJ, Jval, coef_tol, print_tol, leading_coef_tol)

  case('COS2THETA')
    tens%func => rotme_3cos2theta_min1
    dj = 2
    call tens%init(jmin, jmax, dj, verbose=.true.)

    call init_extf_vib_me_overlap(tens%nelem)

    call rovib_me_storeall(tens, nJ, Jval, coef_tol, print_tol, leading_coef_tol)

  case('RICHMOL_LEVELS_FILE')

    call store_richmol_enr(tens, nJ, Jval, coef_tol, print_tol, leading_coef_tol)

  case('MF_TENSOR')
    tens%func => rotme_mf
    dj = 0
    call tens%init(jmin, jmax, dj, verbose=.true.)

    write(sielem,*) oper_ielem
    tens%name = 'MF'//trim(adjustl(sielem))

    call read_extf_vib_me_ielem(oper_ielem)

    call rovib_me_storeall(tens, nJ, Jval, coef_tol, print_tol, leading_coef_tol)

  case default
    write(out, '(/3a)') 'emfield2/emf_matelem error: unexpected operator type = "', trim(oper), '"'
    stop 'STOP'

  end select

  deallocate(Jval)

end subroutine emf2_matelem


!###################################################################################################################################

! Creates Richmol files with matrix elements of a selected Cartesian tensor
! operator, for different pairs of bra and ket J quantum nunbers.

subroutine rovib_me_storeall(tens, nJ, Jval, coef_tol, print_tol, leading_coef_tol_)

  class(rotme_cart_tens_type), intent(in) :: tens
  integer(ik), intent(in) :: nJ, Jval(nJ)
  real(rk), intent(in) :: coef_tol, print_tol
  real(rk), intent(in), optional :: leading_coef_tol_

  integer(ik) :: nlevels(1000), ilevel, jind, info, maxnlevels, jind1, jind2
  integer(ik), allocatable :: level_ind(:,:)
  real(rk) :: leading_coef_tol

  leading_coef_tol = 1.0d-3
  if (present(leading_coef_tol_)) leading_coef_tol = leading_coef_tol_

  ! select levels that pass the energy and quanta filters, specified in "INTENSITY..END" input structure

  nlevels(:) = 0

  do ilevel=1, Neigenlevels
    if (enr_filter_intens( eigen(ilevel)%jval, eigen(ilevel)%quanta(1:), eigen(ilevel)%normal(0:), &
                           eigen(ilevel)%energy, 1 )) then
      jind = eigen(ilevel)%jind
      nlevels(jind) = nlevels(jind) + 1
      if (.not.enr_filter_intens( eigen(ilevel)%jval, eigen(ilevel)%quanta(1:), eigen(ilevel)%normal(0:), &
                                  eigen(ilevel)%energy, 2 )) then
        write(out, '(/a)') 'emfield2/rovib_me_storeall error: rovibrational state filters in "INTENSITY..END" &
            block are different for lower and upper states (must be the same)'
        stop 'STOP, error in emfield2/rovib_me_storeall'
      endif
    endif
  enddo

  maxnlevels = maxval(nlevels)
  allocate(level_ind(maxnlevels,nJ), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'emfield2/rovib_me_storeall error: failed to allocate &
        level_ind(maxnlevels)', 'maxnlevels =', maxnlevels
    stop 'STOP, error in emfield2/rovib_me_storeall'
  endif

  nlevels(:) = 0

  do ilevel=1, Neigenlevels
    if (enr_filter_intens( eigen(ilevel)%jval, eigen(ilevel)%quanta(1:), eigen(ilevel)%normal(0:), &
                           eigen(ilevel)%energy, 1 )) then
      jind = eigen(ilevel)%jind
      nlevels(jind) = nlevels(jind) + 1
      level_ind(nlevels(jind),jind) = ilevel
    endif
  enddo

  ! store rovibrational matrix elements of a tensor for different pairs of J-quanta

  do jind1=1, nJ
    do jind2=1, nJ

      !if (abs(jval(jind1)-jval(jind2))>tens%dj) cycle

      if (jval(jind1)/=intensity%J(1) .or. jval(jind2)/=intensity%J(2)) cycle

      call rovib_me_jpair( tens, nJ, Jval, jind1, jind2, nlevels(jind1), level_ind(1:nlevels(jind1),jind1), &
                           nlevels(jind2), level_ind(1:nlevels(jind2),jind2), coef_tol, print_tol )

    enddo
  enddo

  deallocate(level_ind)

end subroutine rovib_me_storeall


!###################################################################################################################################

! Creates Richmol-energy and Richmol-coefficients files.

subroutine store_richmol_enr(tens, nJ, Jval, coef_tol, print_tol, leading_coef_tol_)

  class(rotme_cart_tens_type), intent(in) :: tens
  integer(ik), intent(in) :: nJ, Jval(nJ)
  real(rk), intent(in) :: coef_tol, print_tol
  real(rk), intent(in), optional :: leading_coef_tol_

  integer(ik) :: nlevels(1000), ilevel, jind, info, maxnlevels, jind1, jind2
  integer(ik), allocatable :: level_ind(:,:)
  real(rk) :: leading_coef_tol

  leading_coef_tol = 1.0d-3
  if (present(leading_coef_tol_)) leading_coef_tol = leading_coef_tol_

  ! select levels that pass the energy and quanta filters, specified in "INTENSITY..END" input structure

  nlevels(:) = 0

  do ilevel=1, Neigenlevels
    if (enr_filter_intens( eigen(ilevel)%jval, eigen(ilevel)%quanta(1:), eigen(ilevel)%normal(0:), &
                           eigen(ilevel)%energy, 1 )) then
      jind = eigen(ilevel)%jind
      nlevels(jind) = nlevels(jind) + 1
      if (.not.enr_filter_intens( eigen(ilevel)%jval, eigen(ilevel)%quanta(1:), eigen(ilevel)%normal(0:), &
                                  eigen(ilevel)%energy, 2 )) then
        write(out, '(/a)') 'emfield2/rovib_me_storeall error: rovibrational state filters in "INTENSITY..END" &
            block are different for lower and upper states (must be the same)'
        stop 'STOP, error in emfield2/rovib_me_storeall'
      endif
    endif
  enddo

  maxnlevels = maxval(nlevels)
  allocate(level_ind(maxnlevels,nJ), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'emfield2/rovib_me_storeall error: failed to allocate &
        level_ind(maxnlevels)', 'maxnlevels =', maxnlevels
    stop 'STOP, error in emfield2/rovib_me_storeall'
  endif

  nlevels(:) = 0

  do ilevel=1, Neigenlevels
    if (enr_filter_intens( eigen(ilevel)%jval, eigen(ilevel)%quanta(1:), eigen(ilevel)%normal(0:), &
                           eigen(ilevel)%energy, 1 )) then
      jind = eigen(ilevel)%jind
      nlevels(jind) = nlevels(jind) + 1
      level_ind(nlevels(jind),jind) = ilevel
    endif
  enddo

  ! store energies in RichMol format

  call store_energies(nJ, Jval, nlevels, level_ind)
  call store_wf_leading(nJ, Jval, nlevels, level_ind, leading_coef_tol)

  deallocate(level_ind)

end subroutine store_richmol_enr


!###################################################################################################################################

! Creates Richmol matrix elements file for a given pair of the bra and ket J quantum numbers.

subroutine rovib_me_jpair( tens, nJ, Jval, jind1, jind2, nlevels1, level_ind1, nlevels2, level_ind2, &
                           coef_tol, print_tol )

  class(rotme_cart_tens_type), intent(in) :: tens
  integer(ik), intent(in) :: jind1, jind2, nlevels1, nlevels2, level_ind1(nlevels1), level_ind2(nlevels2), &
      nJ, Jval(nJ)
  real(rk), intent(in) :: coef_tol, print_tol

  integer(ik) :: nirrep, jval1, jval2, dimen1, dimen2, jind, ilevel_, ilevel, jlevel, jlevel_, isym1, &
      isym2, nsize1, nsize2, eigunit1, eigunit2, irec, ndeg1, ndeg2, irrep, ideg, jdeg, info, nsize, &
      maxdeg, maxdimen, iounit_me, ielem, m1, m2, isym, Jeigenvec_unit(size(Jval),sym%Nrepresen), &
      isign, icmplx, num_threads, n, ithread
  integer(ik), allocatable :: ind_sparse(:,:), nelem_sparse(:)
  integer(ik), external :: omp_get_max_threads, omp_get_thread_num
  real(rk) :: energy1, energy2, nu
  real(rk), allocatable :: vec_sym1(:,:), vec_sym2(:,:), vec_sparse(:,:,:), half_me(:,:,:), me(:,:,:,:)
  real(rk), external :: ddoti
  character(cl) :: fname_me, sj1, sj2
  logical :: tran_filter


  jval1 = Jval(jind1)
  jval2 = Jval(jind2)

  write(out, '(/a,1x,i3,1x,i3,1x,a)') 'Compute and store rovibrational matrix elements for J pair = ', &
      jval1, jval2, '(rovib_me_jpair)'

  if (jval1<tens%jmin.or.jval1>tens%jmax) then
    write(out, '(/a,1x,i3,1x,a,a,a)') 'emfield2/rovib_me_jpair error: initial state J value =', jval1, &
        'runs out of bounds for tensor "', trim(tens%name), '"'
    !stop 'STOP, error in emfield2/rovib_me_jpair'
    return
  endif

  if (jval2<tens%jmin.or.jval2>tens%jmax) then
    write(out, '(/a,1x,i3,1x,a,a,a)') 'emfield2/rovib_me_jpair error: finale state J value =', jval2, &
        'runs out of bounds for tensor "', trim(tens%name), '"'
    !stop 'STOP, error in emfield2/rovib_me_jpair'
    return
  endif


  nirrep = tens%nirrep

  dimen1 = bset_contr(jind1)%Maxcontracts
  dimen2 = bset_contr(jind2)%Maxcontracts

  num_threads = omp_get_max_threads()

  write(out, '(/1x,a,1x,i3)') 'number of parallel threads:', num_threads


  ! allocate workspace arrays

  nsize = 0
  maxdeg = 0
  do ilevel_=1, nlevels1
    ilevel = level_ind1(ilevel_)
    isym = eigen(ilevel)%igamma
    nsize = max(nsize,bset_contr(jind1)%nsize(isym))
    maxdeg = max(maxdeg,eigen(ilevel)%ndeg)
  enddo
  do ilevel_=1, nlevels2
    ilevel = level_ind2(ilevel_)
    isym = eigen(ilevel)%igamma
    nsize = max(nsize,bset_contr(jind2)%nsize(isym))
    maxdeg = max(maxdeg,eigen(ilevel)%ndeg)
  enddo
  maxdimen = max(dimen1,dimen2)

  write(out, '(/1x,a,1x,i6)') 'max sum-of-products dimension:', maxdimen
  write(out, '(1x,a,1x,i6)') 'max symmetrised dimension:', nsize
  write(out, '(1x,a,1x,i6)') 'max number of degenerate components:', maxdeg

  allocate( vec_sym1(nsize,nlevels1), vec_sym2(nsize,nlevels2), vec_sparse(maxdimen,maxdeg,0:num_threads-1), &
      ind_sparse(maxdimen,0:num_threads-1), half_me(dimen2,maxdeg,nirrep), &
      me(nirrep,maxdeg,maxdeg,0:num_threads-1), nelem_sparse(0:num_threads-1), stat=info )
  if (info/=0) then
    write(out, '(/a/a/a,10(1x,i8))') 'emfield2/rovib_me_jpair error: failed to allocate vec_sym1(nsize,nlevels1), &
        vec_sym2(nsize,nlevels2),', 'vec_sparse(maxdimen,maxdeg,0:num_threads-1), ind_sparse(maxdimen,0:num_threads-1), &
        half_me(dimen2,maxdeg,nirrep), me(nirrep,maxdeg,maxdeg,0:num_threads-1)', 'nsize, nlevels1, &
        nlevels2, maxdimen, maxdeg, dimen2, nirrep, num_threads =', nsize, nlevels1, nlevels2, maxdimen, &
        maxdeg, dimen2, nirrep, num_threads
    stop 'STOP, error in emfield2/rovib_me_jpair'
  endif


  ! list of file units with eigenvectors

  do jind=1, nJ
    do isym=1, sym%Nrepresen
      Jeigenvec_unit(jind,isym) = TReigenvec_unit(jind,Jval,isym)
    enddo
  enddo


  ! read eigenvectors for all initial and final states (symmetrized representation)

  write(out, '(/1x,a,1x,i3,1x,a,1x,i6)') 'read eigenvectors for initial states, J =', jval1, &
      ', no.levels =', nlevels1

  do ilevel_=1, nlevels1
    ilevel   = level_ind1(ilevel_)
    isym1    = eigen(ilevel)%igamma
    nsize1   = bset_contr(jind1)%nsize(isym1)
    eigunit1 = Jeigenvec_unit(jind1,isym1)
    irec     = eigen(ilevel)%irec(1)
    read(eigunit1, rec=irec) vec_sym1(1:nsize1,ilevel_)
  enddo

  write(out, '(/1x,a,1x,i3,1x,a,1x,i6)') 'read eigenvectors for final states, J =', jval2, &
      ', no.levels =', nlevels2

  do ilevel_=1, nlevels2
    ilevel   = level_ind2(ilevel_)
    isym2    = eigen(ilevel)%igamma
    nsize2   = bset_contr(jind2)%nsize(isym2)
    eigunit2 = Jeigenvec_unit(jind2,isym2)
    irec     = eigen(ilevel)%irec(1)
    read(eigunit2, rec=irec) vec_sym2(1:nsize2,ilevel_)
  enddo


  ! open file to store matrix elements

  write(sj1,'(i4)') jval1
  write(sj2,'(i4)') jval2
  fname_me = 'matelem_'//trim(tens%name)//'_j'//trim(adjustl(sj1))//'_j'//trim(adjustl(sj2))//'.rchm'
  call IOStart(fname_me, iounit_me)
  write(out, '(1x,a,a,a,1x,i3,1x,i3,1x,a,1x,i5)') 'open file "', trim(fname_me), &
      '" to store matrix elements for j1/j2 = (', jval1, jval2, '), I/O unit =', iounit_me
  open(iounit_me, form='formatted', action='write', position='rewind', status='unknown', file=fname_me, iostat=info)
  if (info/=0) then
    write(out, '(/a,1x,a)') 'emfield2/rovib_me_jpair error while opening file', trim(fname_me)
    stop 'STOP, error in emfield2/rovib_me_jpair'
  endif

  rewind(iounit_me)
  write(iounit_me,'(a)') 'Start richmol format'
  write(iounit_me,'(a,i4,1x,i4)') trim(tens%name), tens%nirrep, tens%nelem

  write(iounit_me, '(a)') 'M-tensor'
  do ielem=1, tens%nelem
    if (tens%mmat_cmplx(ielem)==-1.and.tens%kmat_cmplx==-1) then
      icmplx = 0
      isign = -1
    elseif (tens%mmat_cmplx(ielem)==-1.and.tens%kmat_cmplx==0) then
      icmplx = -1
      isign = 1
    elseif (tens%mmat_cmplx(ielem)==0.and.tens%kmat_cmplx==-1) then
      icmplx = -1
      isign = 1
    elseif (tens%mmat_cmplx(ielem)==0.and.tens%kmat_cmplx==0) then
      icmplx = 0
      isign = 1
    else
      write(out, '(/a,1x,i2,1x,a,1x,i2,1x,a,1x,i3)') 'emfield2/rovib_me_jpair error: invalid combination &
          of tens%mmat_cmplx(ielem) =', tens%mmat_cmplx(ielem), 'and tens%kmat_cmplx(ielem) =', &
          tens%kmat_cmplx, 'ielem =', ielem
      stop 'STOP, error in emfield2/rovib_me_jpair'
    endif
    write(iounit_me, '(a,1x,i4,1x,i2,1x,a)') 'alpha', ielem, icmplx, trim(tens%selem(ielem))
    do m1=-jval1, jval1
      do m2=-jval2, jval2
        if (abs(m1-m2)>tens%dm) cycle
        if (any(abs(tens%mmat(jval1,jval2)%me(1:nirrep,ielem,m1,m2))>print_tol)) then
          write(iounit_me,'(i4,1x,i4,100(1x,f))') m1, m2, tens%mmat(jval1,jval2)%me(1:nirrep,ielem,m1,m2) * isign
        endif
      enddo
    enddo
  enddo

  write(iounit_me, '(a)') 'K-tensor'


  ! start calculations of matrix elements

  do ilevel_=1, nlevels1

    ithread = 0

    ilevel = level_ind1(ilevel_)

    jval1    = eigen(ilevel)%jval
    energy1  = eigen(ilevel)%energy
    isym1    = eigen(ilevel)%igamma
    ndeg1    = eigen(ilevel)%ndeg
    nsize1   = bset_contr(jind1)%nsize(isym1)

    ! transform symmetrized eigenvector "vec_sym1" to sum-of-products of rotational and vibrational functions

    call desym_eigvec( jind1, isym1, ndeg1, vec_sym1(1:nsize1,ilevel_), coef_tol, nelem_sparse(ithread), &
        ind_sparse(:,ithread), vec_sparse(:,:,ithread) )

    ! half transform matrix elements into eigenfunction representation for initial state

    n = nelem_sparse(ithread)

    call half1_rovib_me(tens, jind1, ndeg1, dimen1, n, ind_sparse(1:n,ithread), vec_sparse(1:n,1:ndeg1,ithread), &
                        jind2, dimen2, half_me(1:dimen2,1:ndeg1,1:nirrep))

    !$omp parallel do private(jlevel_,ithread,jlevel,jval2,energy2,isym2,ndeg2,nsize2,nu,tran_filter,irrep,jdeg,ideg,n) schedule(dynamic)
    do jlevel_=1, nlevels2

      ithread = omp_get_thread_num()

      jlevel = level_ind2(jlevel_)

      jval2    = eigen(jlevel)%jval
      energy2  = eigen(jlevel)%energy
      isym2    = eigen(jlevel)%igamma
      ndeg2    = eigen(jlevel)%ndeg
      nsize2   = bset_contr(jind2)%nsize(isym2)

      ! transform symmetrized eigenvector "vec_sym2" to sum-of-products of rotational and vibrational functions

      call desym_eigvec( jind2, isym2, ndeg2, vec_sym2(1:nsize2,jlevel_), coef_tol, nelem_sparse(ithread), &
          ind_sparse(:,ithread), vec_sparse(:,:,ithread) )

      nu = energy2 - energy1

      ! check if current initial/final state pair passes the transition filter, specified in "INTENSITY..END" input structure

      tran_filter = tran_filter_intens(jval1, isym1, energy1, jval2, isym2, energy2, nu, calc_intens=.false.)
      if (.not.tran_filter) cycle

      ! full-transform to eigen-basis

      do irrep=1, nirrep
        do jdeg=1, ndeg2
          do ideg=1, ndeg1
            n = nelem_sparse(ithread)
            me(irrep,ideg,jdeg,ithread) = ddoti(n, vec_sparse(1:n,jdeg,ithread), ind_sparse(1:n,ithread), &
                half_me(1:dimen2,ideg,irrep))
          enddo
        enddo
      enddo

      do jdeg=1, ndeg2
        do ideg=1, ndeg1
          if (any(abs(me(1:nirrep,ideg,jdeg,ithread))>print_tol)) then
            !$omp critical
            write(iounit_me,'(i8,1x,i8,1x,i4,1x,i4,100(1x,f))') ilevel_, jlevel_, ideg, jdeg, &
                me(1:nirrep,ideg,jdeg,ithread)
            !$omp end critical
          endif
        enddo
      enddo

    enddo ! jlevel_
    !$omp end parallel do

  enddo ! ilevel_

  write(iounit_me,'(a)') 'End richmol format'

  close(iounit_me)
  call IOStop(fname_me)

  deallocate(ind_sparse,vec_sym1,vec_sym2,vec_sparse,half_me,me,nelem_sparse)

  write(out, '(a)') 'done (rovib_me_jpair)'

end subroutine rovib_me_jpair


!###################################################################################################################################


subroutine desym_eigvec(jind, isym, ndeg, vec_sym, coef_tol, nelem_sparse, ind_sparse, vec_sparse)

  integer(ik), intent(in) :: jind, isym, ndeg
  real(rk), intent(in) :: vec_sym(:), coef_tol
  integer(ik), intent(out) :: nelem_sparse, ind_sparse(:)
  real(rk), intent(out) :: vec_sparse(:,:)

  integer(ik) :: dimen, idimen, irow, ib, iterm, nelem, ielem, isroot, ideg, nsymcoefs, nrepresen, info, irep, Nterms
  real(rk) :: coef0(ndeg)
  real(rk), allocatable :: kmat(:,:)

  dimen = bset_contr(jind)%Maxcontracts
  nsymcoefs = bset_contr(jind)%Maxsymcoeffs
  nrepresen = sym%Nrepresen

  allocate(kmat(nsymcoefs,nrepresen), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'emfield2/desym_eigvec error: failed to allocate kmat(nsymcoefs,nrepresen)', &
    'nsymcoefs, nrepresen =', nsymcoefs, nrepresen
    stop 'STOP, error in emfield2/desym_eigvec'
  endif
  kmat = 0

  do irep=1, nrepresen
    Nterms = 0
    do iterm=1, nsymcoefs
      kmat(iterm,irep) = Nterms
      Nterms = Nterms + bset_contr(jind)%irr(irep)%N(iterm)
    enddo
  enddo

  nelem_sparse = 0
  do idimen=1, dimen
    irow = bset_contr(jind)%icontr2icase(idimen,1)
    ib   = bset_contr(jind)%icontr2icase(idimen,2)
    iterm = kmat(irow,isym)
    coef0 = 0
    nelem = bset_contr(jind)%irr(isym)%N(irow)
    do ielem=1, nelem
      isroot = iterm + ielem
      do ideg=1, ndeg
        coef0(ideg) = coef0(ideg) + vec_sym(isroot)*bset_contr(jind)%irr(isym)%repres(isroot,ideg,ib)
      enddo
    enddo
    if (any(abs(coef0(1:ndeg))>coef_tol)) then
      nelem_sparse = nelem_sparse + 1
      ind_sparse(nelem_sparse) = idimen
      vec_sparse(nelem_sparse,1:ndeg) = coef0(1:ndeg)
    endif
  enddo

  deallocate(kmat)

end subroutine desym_eigvec



!###################################################################################################################################



subroutine half1_rovib_me(tens, jind1, ndeg1, dimen1, nelem1, ind1, coefs1, jind2, dimen2, half_me)

  class(rotme_cart_tens_type), intent(in) :: tens
  integer(ik), intent(in) :: jind1, jind2, dimen1, dimen2, ndeg1, nelem1, ind1(nelem1)
  real(rk), intent(in) :: coefs1(nelem1,ndeg1)
  real(rk), intent(out) :: half_me(dimen2,ndeg1,tens%nirrep)

  integer(ik) :: jdimen, ithread, nirrep, info, num_threads
  integer(ik), external :: omp_get_max_threads, omp_get_thread_num
  real(rk), allocatable :: tvec(:,:,:)

  nirrep = tens%nirrep
  num_threads = omp_get_max_threads()

  allocate(tvec(nelem1,nirrep,0:num_threads-1), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'emfield2/half1_rovib_me error: failed to allocate tvec(nelem1,nirrep,0:num_threads-1)', &
    'nelem1, nirrep, num_threads =', nelem1, nirrep, num_threads
    stop 'STOP, error in emfield2/half1_rovib_me'
  endif

  !$omp parallel do private(jdimen,ithread) schedule(dynamic)
  do jdimen=1, dimen2

    ithread = omp_get_thread_num()

    call prim_me(tens, jind1, jind2, jdimen, nelem1, ind1(1:nelem1), tvec(1:nelem1,1:nirrep,ithread))

    call dgemm('T', 'N', ndeg1, nirrep, nelem1, 1.0d0, coefs1(1:nelem1,1:ndeg1), nelem1, &!
               tvec(1:nelem1,1:nirrep,ithread), nelem1, 0.0d0, &!
               half_me(jdimen,1:ndeg1,1:nirrep), ndeg1)
  enddo
 !$omp end parallel do

  deallocate(tvec)

end subroutine half1_rovib_me



!###################################################################################################################################



subroutine prim_me(tens, jind1, jind2, idimen2, nelem, ind, res_vec)

  class(rotme_cart_tens_type), intent(in) :: tens
  integer(ik), intent(in) :: jind1, jind2, idimen2, nelem, ind(nelem)
  real(rk), intent(out) :: res_vec(:,:)

  integer(ik) :: j1, j2, dimen1, dimen2, irrep, idimen, icontr1, icontr2, k1, k2, tau1, tau2, ktau1, &
      ktau2, ktau1_, ktau2_, nirrep, ncart, ielem, dk
  real(rk) :: rot_me(tens%nelem), vib_me(tens%nelem)

  ncart = tens%nelem
  nirrep = tens%nirrep
  j1 = bset_contr(jind1)%jval
  j2 = bset_contr(jind2)%jval

#if defined(_EMFIELD2_DEBUG_)
  dimen2 = bset_contr(jind2)%Maxcontracts
  if (idimen2>dimen2) then
    write(out, '(/a,1x,i6,1x,a,1x,i6)') &!
    'prim_me error: primitive function index for final state =', idimen2,' exceeds dimension of the basis =', dimen2
    stop 'STOP, error in emfield2/prim_me'
  endif
#endif

  icontr2 = bset_contr(jind2)%iroot_correlat_j0(idimen2)
  ktau2   = bset_contr(jind2)%ktau(idimen2)

#if defined(_EMFIELD2_DEBUG_)
  k2     = bset_contr(jind2)%k(idimen2)
  tau2   = mod(ktau2,2)
  ktau2_ = tens%ktau_ind(k2,tau2)
  if (ktau2_<=0) then
    write(out, '(/a,1x,i6)') 'prim_me error: tens%ktau_ind(k2,tau2) =', ktau2_
    stop 'STOP, error in emfield2/prim_me'
  endif
#else
  ktau2_ = ktau2
#endif

  res_vec = 0

  do irrep=1, tens%nirrep

    do ielem=1, nelem

      idimen = ind(ielem)

#if defined(_EMFIELD2_DEBUG_)
      dimen1 = bset_contr(jind1)%Maxcontracts
      if (idimen>dimen1) then
        write(out, '(/a,1x,i6,1x,a,1x,i6)') &!
        'prim_me error: primitive function index for initial state =', idimen,' exceeds dimension of the basis =', dimen1
        stop 'STOP, error in emfield2/prim_me'
      endif
#endif

      icontr1 = bset_contr(jind1)%iroot_correlat_j0(idimen)
      ktau1   = bset_contr(jind1)%ktau(idimen)

#if defined(_EMFIELD2_DEBUG_)
      k1     = bset_contr(jind1)%k(idimen)
      tau1   = mod(ktau1,2)
      ktau1_ = tens%ktau_ind(k1,tau1)
      if (ktau1_==0) then
        write(out, '(/a,1x,i6)') 'prim_me error: tens%ktau_ind(k1,tau1) =', ktau1_
        stop 'STOP, error in emfield2/prim_me'
      endif
#else
      ktau1_ = ktau1
#endif

      dk = abs(ktau1_-ktau2_)

      if (dk<=tens%dk) then

        rot_me(1:ncart) = tens%kmat(j1,j2)%me(1:ncart,irrep,ktau1_,ktau2_)
        vib_me(1:ncart) = extf_vib_me(1:ncart,icontr1,icontr2)

        res_vec(ielem,irrep) = sum(rot_me(1:ncart) * vib_me(1:ncart))

      endif

    enddo

  enddo

end subroutine prim_me



!###################################################################################################################################



subroutine read_extf_vib_me(rank)

  integer(ik), intent(in) :: rank

  integer(ik) :: ncontr_t, irank, irank_t, info, chkptIO, i, j
  character(len=cl) :: job_is
  character(len=20) :: buf20

  write(out, '(/a,a,a)') 'read_extf_vib_me: read vibrational contracted matrix elements from file "', &
      trim(job%extFmat_file), '"'

  if (rank/=extF%rank) then
    write(out, '(/a,1x,i4,1x,a,1x,i4)') 'emfield2/read_extf_vib_me error: rank of Cartesian tensor =', &
        rank, 'does not agree with the rank of TROVE extF tensor =', extF%rank
    stop 'STOP, error in emfield2/read_extf_vib_me'
  endif

  job_is ='extf contracted matrix elements'
  call IOStart(trim(job_is),chkptIO)
  open(chkptIO, form='unformatted', action='read', position='rewind', status='old', file=job%extFmat_file)

  read(chkptIO) buf20
  if (buf20/='Start external field') then
    write (out, '(/a,a,a,a,a)') 'emfield2/read_extf_vib_me error: file "', trim(job%extFmat_file), &
        '" has bogus header = "', buf20, '"'
    stop 'STOP, error in emfield2/read_extf_vib_me'
  endif

  read(chkptIO) ncontr_t

  if (bset_contr(1)%Maxcontracts/=ncontr_t) then
    write (out, '(/a,1x,i6,1x,a,1x,i6,1x,a)') 'emfield2/read_extf_vib_me error: actual size of basis &
        set =',  bset_contr(1)%Maxcontracts, 'and stored one =', ncontr_t, 'do not agree'
    stop 'STOP, error in emfield2/read_extf_vib_me'
  endif

  if (rank<=0) then
    write(out, '(/a,1x,i3)') 'emfield2/read_extf_vib_me error: rank of external function =', rank
    stop 'STOP, error in emfield2/read_extf_vib_me'
  endif

  if (allocated(extf_vib_me)) deallocate(extf_vib_me)
  allocate(extf_vib_me(rank,ncontr_t,ncontr_t), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'emfield2/read_extf_vib_me error: failed to allocate &
        extf_vib_me(rank,ncontr_t,ncontr_t)', 'ncontr_t, rank =', ncontr_t, rank
    stop 'STOP, error in emfield2/read_extf_vib_me'
  endif
  extf_vib_me = 0.0

  do irank=1, rank

    read(chkptIO) irank_t
    if (irank_t/=irank) then
      write (out, '(/a,a,a,1x,i3,1x,a,1x,i3)') 'emfield2/read_extf_vib_me error: file "', &
          trim(job%extFmat_file), '" has bogus irank = ', irank_t, ', expected irank =', irank
      stop
    endif

    read(chkptIO) extf_vib_me(irank,:,:)

  enddo

  read(chkptIO) buf20(1:18)
  if (buf20(1:18)/='End external field') then
    write (out, '(/a,a,a,a,a)') 'emfield2/read_extf_vib_me error: file "', trim(job%extFmat_file), &
        '" has bogus footer = "', buf20(1:18), '"'
    stop 'STOP, error in emfield2/read_extf_vib_me'
  endif

  ! print vibrational matrix elements
  !do i=1, ncontr_t
  !  do j=1, i
  !    write(out, '(1x,i6,1x,i6,100(1x,f))') i,j, extf_vib_me(:,i,j)
  !  enddo
  !enddo

  close(chkptIO)
  call IOStop(job_is)

end subroutine read_extf_vib_me



subroutine init_extf_vib_me_overlap(rank)

  integer(ik), intent(in) :: rank

  integer(ik) :: ncontr_t, irank, info, i

  ncontr_t = bset_contr(1)%Maxcontracts

  if (allocated(extf_vib_me)) deallocate(extf_vib_me)
  allocate(extf_vib_me(rank,ncontr_t,ncontr_t), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'emfield2/init_extf_vib_me_overlap error: failed to allocate &
        extf_vib_me(rank,ncontr_t,ncontr_t)', 'ncontr_t, rank =', ncontr_t, rank
    stop 'STOP, error in emfield2/init_extf_vib_me_overlap'
  endif

  extf_vib_me = 0.0

  do irank=1, rank

    do i=1, ncontr_t
      extf_vib_me(irank,i,i) = 1.0_ark
    enddo

  enddo

end subroutine init_extf_vib_me_overlap



subroutine read_extf_vib_me_ielem(ielem)

  integer(ik), intent(in) :: ielem

  integer(ik) :: ncontr_t, irank, irank_t, info, chkptIO, i, j, rank
  character(len=cl) :: job_is
  character(len=20) :: buf20

  rank = 1

  write(out, '(/a,a,a)') 'read_extf_vib_me_ielem: read vibrational contracted matrix elements from &
        file "', trim(job%extFmat_file), '"'

  if (ielem>extF%rank) then
    write(out, '(/a,1x,i4,1x,a,1x,i4)') 'emfield2/read_extf_vib_me_ielem error: index of Cartesian &
        tensor =', ielem, 'is larger than rank of TROVE extF tensor =', extF%rank
    stop 'STOP, error in emfield2/read_extf_vib_me_ielem'
  endif

  job_is ='extf contracted matrix elements'
  call IOStart(trim(job_is),chkptIO)
  open(chkptIO, form='unformatted', action='read', position='rewind', status='old', file=job%extFmat_file)

  read(chkptIO) buf20
  if (buf20/='Start external field') then
    write (out, '(/a,a,a,a,a)') 'emfield2/read_extf_vib_me_ielem error: file "', trim(job%extFmat_file), &
        '" has bogus header = "', buf20, '"'
    stop 'STOP, error in emfield2/read_extf_vib_me_ielem'
  endif

  read(chkptIO) ncontr_t

  if (bset_contr(1)%Maxcontracts/=ncontr_t) then
    write (out, '(/a,1x,i6,1x,a,1x,i6,1x,a)') 'emfield2/read_extf_vib_me_ielem error: actual size of &
        basis set =',  bset_contr(1)%Maxcontracts, 'and stored one =', ncontr_t, 'do not agree'
    stop 'STOP, error in emfield2/read_extf_vib_me_ielem'
  endif

  if (rank<=0) then
    write(out, '(/a,1x,i3)') 'emfield2/read_extf_vib_me_ielem error: rank of external function =', rank
    stop 'STOP, error in emfield2/read_extf_vib_me_ielem'
  endif

  if (allocated(extf_vib_me)) deallocate(extf_vib_me)
  allocate(extf_vib_me(rank,ncontr_t,ncontr_t), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'emfield2/read_extf_vib_me_ielem error: failed to allocate &
        extf_vib_me(rank,ncontr_t,ncontr_t)', 'ncontr_t, rank =', ncontr_t, rank
    stop 'STOP, error in emfield2/read_extf_vib_me_ielem'
  endif
  extf_vib_me = 0.0

  do irank=1, extF%rank

    read(chkptIO) irank_t
    if (irank_t/=irank) then
      write (out, '(/a,a,a,1x,i3,1x,a,1x,i3)') 'emfield2/read_extf_vib_me_ielem error: file "', &
          trim(job%extFmat_file), '" has bogus irank = ', irank_t, ', expected irank =', irank
      stop
    endif

    read(chkptIO) extf_vib_me(1,:,:)
    if (irank==ielem) exit

  enddo

  read(chkptIO) buf20(1:18)
  if (buf20(1:18)/='End external field') then
    write (out, '(/a,a,a,a,a)') 'emfield2/read_extf_vib_me_ielem error: file "', trim(job%extFmat_file), &
        '" has bogus footer = "', buf20(1:18), '"'
    stop 'STOP, error in emfield2/read_extf_vib_me_ielem'
  endif

  ! print vibrational matrix elements
  !do i=1, ncontr_t
  !  do j=i, i
  !    write(out, '(1x,i6,1x,i6,100(1x,f))') i,j, extf_vib_me(:,i,j)
  !  enddo
  !enddo

  close(chkptIO)
  call IOStop(job_is)

end subroutine read_extf_vib_me_ielem


!###################################################################################################################################

! Stores the rovibrational energies of selected state into a Richmol-energy file.

subroutine store_energies(nJ, Jval, nlevels, level_ind)

  integer(ik), intent(in) :: nJ, Jval(nJ), nlevels(nJ), level_ind(:,:)

  integer(ik) :: iounit, ilevel_, ilevel, jind, Jrot, isym, ndeg, nmodes, info, nclasses, nclasses_, &
      nmodes_
  real(rk) :: energy
  character(cl) :: sj1, sj2, fname

  nmodes = molec%nmodes
  nclasses = size(eigen(1)%cgamma)-1

  if (nmodes==0) then
    write(out, '(/a)') 'emfield2/store_energies error: molec%nmodes = 0'
    stop 'STOP, error in emfield2/store_energies'
  endif

  write(sj1,'(i4)') Jval(1)
  write(sj2,'(i4)') Jval(nJ)
  fname = 'energies'//'_j'//trim(adjustl(sj1))//'_j'//trim(adjustl(sj2))//'.rchm'
  call IOStart(trim(fname), iounit)
  open(iounit, form='formatted', action='write', position='rewind', status='unknown', file=fname,iostat=info)
  if (info/=0) then
    write(out, '(/a,1x,a)') 'emfield2/store_energies error while opening file', trim(fname)
    stop 'STOP, error in emfield2/store_energies'
  endif
  write(out, '(/1x,a,1x,a)') 'store rovibrational energies (RichMol format) in file:', trim(fname)
  write(out, '(1x,a,1x,f)') 'zero-point energy:', intensity%ZPE

  write(out, '(/1x,a,1x,i3,1x,a,100(1x,i3))') 'J quanta (', nJ, '):', Jval(1:nJ)
  write(out, '(1x,a,1x,100(1x,i6))') '.. and respective number of energy levels:', nlevels(1:nJ)

  nclasses_ = nclasses + 1
  nmodes_ = nmodes + 1
  do jind=1, nJ
    do ilevel_=1, nlevels(jind)
      ilevel = level_ind(ilevel_,jind)
      Jrot   = eigen(ilevel)%jval
      energy = eigen(ilevel)%energy
      isym   = eigen(ilevel)%igamma
      ndeg   = eigen(ilevel)%ndeg
      write(iounit, '(i4,1x,i8,1x,a5,1x,i4,1x,f20.12,1x,i4,<nmodes>(1x,i4),3x,i8,2x,<nclasses_>(1x,a5),2x,<nmodes_>(1x,i4),3x,es16.8)') &
          Jrot, ilevel_, sym%label(isym), ndeg, energy-intensity%ZPE, eigen(ilevel)%krot, &
          eigen(ilevel)%quanta(1:nmodes), eigen(ilevel)%icoeff, eigen(ilevel)%cgamma(0:nclasses), &
          eigen(ilevel)%normal(0:nmodes), eigen(ilevel)%largest_coeff
    enddo
  enddo

  close(iounit)

end subroutine store_energies


!###################################################################################################################################

! Stores the leading contributions to the rovibraitonal wavefunctions of
! selected states into a Richmol-coefficients file.

subroutine store_wf_leading(nJ, Jval, nlevels, level_ind, leading_coef_tol)

  integer(ik), intent(in) :: nJ, Jval(nJ), nlevels(nJ), level_ind(:,:)
  real(rk), intent(in) :: leading_coef_tol

  integer(ik), parameter :: maxnelem = 1000000
  integer(ik) :: ilevel, jind, Jrot, isym, ndeg, dimen, nsize, eigunit, ideg, irec, icontr, k, tau, &
      ktau, irow, ib, iterm, idimen, nelem, isroot, Jeigenvec_unit(size(Jval),sym%Nrepresen), Nterms, &
      info, maxnsize, maxdimen, maxdeg, ielem, iounit, rot_ncoefs, ilevel_, ksign(2), icoef
  integer(ik), allocatable :: v0(:), img(:), k0(:)
  real(rk) :: energy, coef_tol
  real(rk), allocatable :: vec_sym(:), vec(:,:), coef0(:), vec0(:)
  complex(rk) :: rot_coefs(2), vec0_
  character(cl) :: sj1, sj2, fname
  type(rotme_cart_tens_type) :: tens0
  type DkmatT
    integer(ik),pointer   :: kmat(:,:)
  end type DkmatT
  type(DkmatT), allocatable :: ijterm(:)

  coef_tol = 1.0d-14

  ! open file to store coefficients

  write(sj1,'(i4)') Jval(1)
  write(sj2,'(i4)') Jval(nJ)
  fname = 'coefficients'//'_j'//trim(adjustl(sj1))//'_j'//trim(adjustl(sj2))//'.rchm'
  call IOStart(trim(fname),iounit)
  open(iounit, form='formatted', action='write', position='rewind', status='unknown', file=fname,iostat=info)
  if (info/=0) then
    write(out, '(/a,a,a)') 'emfield2/store_wf_leading error while opening file "', trim(fname), '"'
    stop 'STOP, error in emfield2/store_wf_leading'
  endif
  write(out, '(/1x,a,a,a)') 'store leading contributions in file "', trim(fname), '"'


  ! list of file units with eigenvectors

  do jind=1, nJ
    do isym=1, sym%Nrepresen
      Jeigenvec_unit(jind,isym) = TReigenvec_unit(jind,Jval,isym)
    enddo
  enddo


  ! maximal dimensions

  maxnsize = 0
  maxdimen = 0
  maxdeg = 0
  do jind=1, nJ
    do ilevel_=1, nlevels(jind)
      ilevel = level_ind(ilevel_,jind)
      isym = eigen(ilevel)%igamma
      maxnsize = max(maxnsize,bset_contr(jind)%nsize(isym))
      maxdimen = max(maxdimen,bset_contr(jind)%Maxcontracts)
      maxdeg = max(maxdeg,eigen(ilevel)%ndeg)
    enddo
  enddo

  write(out, '(/1x,a,1x,i3,1x,a,100(1x,i3))') 'J quanta (', nJ, '):', Jval(:)
  write(out, '(1x,a,1x,100(1x,i6))') 'number of energy levels:', nlevels
  write(out, '(1x,a,1x,i6,1x,i6)') 'max dimensions of eigenvector:', maxnsize, maxdimen
  write(out, '(1x,a,1x,i6)') 'max degeneracy:', maxdeg


  ! allocate workspace arrays

  allocate(vec_sym(maxnsize), vec(maxdeg,maxdimen), coef0(maxdeg), ijterm(nJ), vec0(maxnelem), &
      v0(maxnelem), img(maxnelem), k0(maxnelem), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i8))') 'emfield2/store_wf_leading error: failed to allocate vec_sym(maxnsize), &
        vec(maxdeg,maxdimen), coef0(maxdeg), ijterm(nJ), vec0(maxnelem), v0(maxnelem), img(maxnelem), &
        k0(maxnelem)', 'maxnsize, maxdeg, maxdimen, nJ, maxnelem =', maxnsize, maxdeg, maxdimen, nJ, &
        maxnelem
    stop 'STOP, error in emfield2/store_wf_leading'
  endif


  do jind=1, nJ
    dimen = bset_contr(jind)%Maxcontracts
    allocate (ijterm(jind)%kmat(bset_contr(jind)%Maxsymcoeffs,sym%Nrepresen), stat=info)
    if (info/=0) then
      write(out, '(/a/a,10(1x,i6))') 'emfield2/store_wf_leading error: failed to allocate &
          ijterm(jind)%kmat(bset_contr(jind)%Maxsymcoeffs,sym%Nrepresen)', &
          'jind, bset_contr(jind)%Maxsymcoeffs, sym%Nrepresen = ', jind, bset_contr(jind)%Maxsymcoeffs, &
          sym%Nrepresen
      stop 'STOP, error in emfield2/store_wf_leading'
    endif
    do isym=1, sym%Nrepresen
      Nterms = 0
      do iterm=1, bset_contr(jind)%Maxsymcoeffs
        ijterm(jind)%kmat(iterm,isym) = Nterms
        Nterms = Nterms + bset_contr(jind)%irr(isym)%N(iterm)
      enddo
    enddo
  enddo


  do jind=1, nJ

    Jrot = Jval(jind)
    dimen = bset_contr(jind)%Maxcontracts

    ! read eigenvectors (in symmetrized representation)

    do ilevel_=1, nlevels(jind)

      ilevel  = level_ind(ilevel_,jind)

      isym    = eigen(ilevel)%igamma
      Jrot    = eigen(ilevel)%jval
      energy  = eigen(ilevel)%energy
      ndeg    = eigen(ilevel)%ndeg
      dimen   = bset_contr(jind)%Maxcontracts
      nsize   = bset_contr(jind)%nsize(isym)
      eigunit = Jeigenvec_unit(jind,isym)
      irec    = eigen(ilevel)%irec(1)
      read(eigunit, rec=irec) vec_sym(1:nsize)

      ! transform eigenvector to non-symmetrized representation

      do idimen=1, dimen
        irow = bset_contr(jind)%icontr2icase(idimen,1)
        ib   = bset_contr(jind)%icontr2icase(idimen,2)
        iterm = ijterm(jind)%kmat(irow,isym)
        coef0(:) = 0
        nelem = bset_contr(jind)%irr(isym)%N(irow)
        do ielem=1, nelem
          isroot = iterm + ielem
          do ideg=1, ndeg
            coef0(ideg) = coef0(ideg) + vec_sym(isroot)*bset_contr(jind)%irr(isym)%repres(isroot,ideg,ib)
          enddo
        enddo
        vec(1:ndeg,idimen) = coef0(1:ndeg)
      enddo

      ! print into file the leading contributions to the wave function

      do ideg=1, ndeg

        nelem = 0

        do idimen=1, dimen

          if (vec(ideg,idimen)**2>=leading_coef_tol) then

            icontr = bset_contr(jind)%iroot_correlat_j0(idimen)
            k      = bset_contr(jind)%k(idimen)
            ktau   = bset_contr(jind)%ktau(idimen)
            tau    = mod(ktau,2)

            !call symrot_coefs(Jrot, k, tau, rot_ncoefs, rot_coefs)
            call wang_jktau_coefs(tens0, Jrot, k, tau, rot_ncoefs, rot_coefs)

            ksign(1:2) = (/1,-1/)

            do icoef=1, rot_ncoefs

              nelem = nelem + 1

              if (nelem>maxnelem) then
                write(out, '(/a,1x,i6,1x,a,1x,es16.8)') 'emfield2/store_wf_leading error: number &
                    of leading contributions to the wave function exceeds maximum =', maxnelem, ', &
                    coefficient thresh =', leading_coef_tol
                stop 'STOP, error in emfield2/store_wf_leading'
              endif

              vec0_ = vec(ideg,idimen) * rot_coefs(icoef)
              k0(nelem) = k*ksign(icoef)
              v0(nelem) = icontr

              if (abs(real(vec0_,rk))>=coef_tol .and. abs(aimag(vec0_))<coef_tol) then
                img(nelem) = 0
                vec0(nelem) = real(vec0_)
              elseif (abs(real(vec0_,rk))<coef_tol .and. abs(aimag(vec0_))>=coef_tol) then
                img(nelem) = 1
                vec0(nelem) = aimag(vec0_)
              else
                write(out, '(/a,2(1x,f))') 'emfield2/store_wf_leading error: coefficient of the wave &
                    function, re-expressed in terms of symmetric-top functions, is neither pure real &
                        or pure imaginary:', vec0_
                stop 'STOP, error in emfield2/store_wf_leading'
              endif

            enddo !icoef

          endif

        enddo ! idimen

        write(iounit, '(1x,i4,1x,i8,1x,a5,1x,i4,1x,f,1x,i6,<maxnelem>(1x,f,1x,i1,1x,i6,1x,i4))') &
            Jrot, ilevel_, sym%label(isym), ideg, energy-intensity%ZPE, nelem, &
            (vec0(ielem), img(ielem), v0(ielem), k0(ielem), ielem=1, nelem)

      enddo ! ideg

    enddo ! ilevel_
  enddo ! jind1

  close(iounit)

end subroutine store_wf_leading


!###################################################################################################################################

! Filters rovibrational states that will be used for computing the Richmol data.

function enr_filter_intens(jval, vib_quanta, normal, energy, uplow) result(f)

  integer(ik), intent(in) :: jval, vib_quanta(:), uplow, normal(0:)
  real(rk), intent(in) :: energy
  logical :: f

  integer(ik) :: nmodes, i
  logical :: passed_qv, passed_enr, passed_j, passed_vibmom

  nmodes = molec%nmodes

  if (nmodes==0) then
    write(out, '(/a)') 'enr_filter_intens error: molec%nmodes = 0'
    stop 'STOP, error in emfield2/enr_filter_intens'
  endif

  ! J-quanta filter

  passed_j = .false.

  if ( jval>=minval(intensity%J(1:2)) .and. jval<=maxval(intensity%J(1:2)) ) then
    passed_j = .true.
  endif

  ! vibrational quanta filter

  passed_qv = .false.

  if (uplow==1) then ! lower energy level

    if (intensity%nvib_quanta_low>0) then ! use keyword "LOW_STATE"

      do i=1, intensity%nvib_quanta_low
        if (all(vib_quanta(1:nmodes)==intensity%vib_quanta_low(1:nmodes,i))) then
          passed_qv = .true.
          exit
        endif
      enddo

    else ! use keyword "LOWER"/"LOW"/"L"

      if (all(vib_quanta(1:nmodes)>=intensity%v_low(1:nmodes,1)) .and. all(vib_quanta(1:nmodes)<=intensity%v_low(1:nmodes,2)) ) then
        passed_qv = .true.
      endif

    endif

  else ! upper energy level

    if (intensity%nvib_quanta_upp>0) then  ! use keyword "UPP_STATE"

      do i=1, intensity%nvib_quanta_upp
        if (all(vib_quanta(1:nmodes)==intensity%vib_quanta_upp(1:nmodes,i))) then
          passed_qv = .true.
          exit
        endif
      enddo

    else ! use keyword "UPPER"/"UPP"/"U"

      if (all(vib_quanta(1:nmodes)>=intensity%v_upp(1:nmodes,1)) .and. all(vib_quanta(1:nmodes)<=intensity%v_upp(1:nmodes,2)) ) then
        passed_qv = .true.
      endif

    endif

  endif

  ! energy filter

  passed_enr = .false.

  if (uplow==1) then ! lower energy level
    if ( energy-intensity%ZPE>=intensity%erange_low(1) .and. energy-intensity%ZPE<=intensity%erange_low(2) ) then
      passed_enr = .true.
    endif
  else ! upper energy level
    if ( energy-intensity%ZPE>=intensity%erange_upp(1) .and. energy-intensity%ZPE<=intensity%erange_upp(2) ) then
      passed_enr = .true.
    endif
  endif

  ! vibrational angular momentum filter, it is needed to get rid of unphysical
  ! states in linear molecules

  passed_vibmom = .true.
  if (job%triatom_sing_resolve) then
    if (jval==0.and.normal(0)/=0) then
        passed_vibmom = .false.
    endif
  endif 

  f = passed_j .and. passed_qv .and. passed_enr .and. passed_vibmom

end function enr_filter_intens


!###################################################################################################################################



function tran_filter_intens(jval1, isym1, energy1, jval2, isym2, energy2, nu, calc_intens) result(f)

  integer(ik), intent(in) :: jval1, jval2, isym1, isym2
  real(rk), intent(in) :: energy1, energy2, nu
  logical, intent(in) :: calc_intens
  logical :: f

  real(rk) :: gns, small, gns1, gns2
  logical :: passed_e, passed_freq

  small = 1.0d-12
  gns1 =intensity%gns(isym1)
  gns2 =intensity%gns(isym2)

  f = .false.

  passed_e = .true.
  passed_freq = .true.

  if (calc_intens) then
    if (jval1==jval2 .and. energy2<energy1) then
      passed_e = .false.
    endif
    if (abs(nu)<small) then
      passed_freq = .false.
    endif
  endif

  if ( abs(nu)>=intensity%freq_window(1) .and. &!
       abs(nu)<=intensity%freq_window(2) .and. &!
       gns1>small .and. gns2>small .and. &!
       passed_e .and. passed_freq ) f = .true.

end function tran_filter_intens


end module extfield
