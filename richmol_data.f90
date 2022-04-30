module richmol_data

use accuracy
use timer
implicit none


type states_type
  integer(ik)                 ::  nlevels
  integer(ik), allocatable    ::  id(:)
  real(rk), allocatable       ::  energy(:)
  integer(ik), allocatable    ::  ndeg(:)
  integer(ik), allocatable    ::  k_quanta(:)
  integer(ik), allocatable    ::  v_quanta(:,:)
  character(cl), allocatable  ::  sym(:)
  integer(ik), allocatable    ::  isym(:)
  character(cl), allocatable  ::  sym_rot(:)
  character(cl), allocatable  ::  sym_vib(:)
  integer(ik), allocatable    ::  istate(:,:)
end type states_type


type state_filter_type
  integer(ik)  ::  jmin = 0
  integer(ik)  ::  jmax = 0
  integer(ik)  ::  mmin = 0
  integer(ik)  ::  mmax = 0
  real(rk)     ::  emin = 0
  real(rk)     ::  emax = 10000.0d0
  integer(ik)  ::  nmodes = 0
end type state_filter_type


type mmat_type
  integer(ik)               ::  nmpairs = 0
  integer(ik), allocatable  ::  mpair(:,:)     ! mpair(1:2,impair) = [m1,m2]
  integer(ik), allocatable  ::  mpair_ind(:,:) ! mpair_ind(m1,m2) = impair
  real(rk), allocatable     ::  me(:,:)        ! me(1:nirrep,impair)
end type mmat_type


type trans_type
  logical                       ::  null = .false.
  character(cl)                 ::  tens_name = 'N/A'
  integer(ik)                   ::  nelem = 0
  integer(ik)                   ::  nirrep = 0
  integer(ik)                   ::  icmplx = 0
  complex(rk)                   ::  cmplx_fac = 1.0_rk
  integer(ik)                   ::  nkpairs = 0
  integer(ik), allocatable      ::  kpair(:,:)          ! kpair(1:2,ikpair) = [state1,state2]
  real(rk), allocatable         ::  kme(:,:)            ! kme(1:nirrep,ikpair)
  type(mmat_type), allocatable  ::  mmat(:)             ! mmat(1:nelem)
end type trans_type



contains


! Reads rovibrational transition matrix elements from the RichMol database file.
! j2, j1 - bra and ket quantum numbers of the rotational angular momentum;
! tens_name - name of the transition (Cartesian tensor) operator;
! states - array of states_type objects containing the information about the rovibrational states, 
!          i.e. state ID numbers, energies, quantum numbers, etc.,
!          initialized by read_states subroutine;
! me_thresh - threshold for neglecting transition matrix elements;
! trans - trans_type object containing the rovibrational matrix elements;

subroutine read_tran(j1, j2, jmin, jmax, tens_name, states, me_thresh, trans, verbose_)

  integer(ik), intent(in) :: j1, j2, jmin, jmax
  character(*), intent(in) :: tens_name
  type(states_type), intent(in) :: states(jmin:jmax)
  real(rk), intent(in) :: me_thresh
  type(trans_type), intent(out) :: trans
  logical, optional, intent(in) :: verbose_

  integer(ik), parameter :: nmpairs0=10, nmpairs_incr=10, nkpairs0=1000, nkpairs_incr=1000, nirrep_max=100
  integer(ik) :: m1min, m1max, m2min, m2max, iline, iounit, info, nirrep, ncart, nmpairs_max, ielem, ielem_, &
                 impair, m1, m2, m1_, m2_,nkpairs_max, ikpair, istate1, istate2, istate, id1, id2, ideg1, ideg2, &
                 id1_, id2_, ideg1_, ideg2_, icmplx
  real(rk) :: me(nirrep_max)
  character(cl) :: sj1, sj2, fname, fname_transp, tens_name_, sbuf, cart_label, ioname
  character(1000) :: sbuf1000
  logical :: transp, verbose=.false.

  if (present(verbose_)) verbose=verbose_

  if (verbose) write(out, '(100(1h.),1x,a)') 'richmol_data/read_tran start'

  if (verbose) then
    write(out, '(1x,a)') 'Read RichMol transitions file'
    write(out, '(1x,a,1x,a)') 'operator:', trim(tens_name)
    write(out, '(1x,a,1x,i3,1x,i3)') 'J quanta:', j1, j2
    write(out, '(1x,a,1x,es16.8)') 'threshold for neglecting matrix elements:', me_thresh
  endif

  write(sj1,'(i4)') j1
  write(sj2,'(i4)') j2

  fname = 'matelem_'//trim(tens_name)//'_j'//trim(adjustl(sj1))//'_j'//trim(adjustl(sj2))//'.rchm'
  fname_transp = 'matelem_'//trim(tens_name)//'_j'//trim(adjustl(sj2))//'_j'//trim(adjustl(sj1))//'.rchm'

  ioname = fname
  call IOStart(ioname, iounit)
  if (verbose) write(out, '(1x,a,1x,a)') 'read file:', trim(fname)
  open(iounit, form='formatted', action='read', position='rewind', status='old', file=fname, iostat=info)
  if (info/=0) then
    if (verbose) write(out, '(1x,a,1x,a,1x,a)') '..file not found, try reding file:', trim(fname_transp), '(Hermitian transpose)'
    open(iounit, form='formatted', action='read', position='rewind', status='old', file=fname_transp, iostat=info)
    if (info/=0) then
      if (verbose) write(out, '(6x,a,a)') 'both files are not found, perhaps all transitions are zero ',&
                                          'or forbidden by Delta J or Delta M selection rules'
      !write(out, '(/a,1x,a,1x,a,1x,a)') 'richmol_data/read_tran error while opening file', trim(fname), 'or', trim(fname_transp)
      !stop 'STOP, error in richmol_data/read_tran'
      trans%null = .true.
      return
    endif
  	if (verbose) write(out, '(1x,a)') '....OK'
    transp = .true.
    fname = fname_transp
  else
  	if (verbose) write(out, '(1x,a)') '..OK'
    transp = .false.
  endif

  rewind(iounit)
  iline = 0

  read(iounit,'(a)',iostat=info) sbuf
  iline = iline + 1
  call check_iostat
  if (trim(adjustl(sbuf))/='Start richmol format') then
    write(out, '(/a,1x,a,1x,a,1x,a,1x,a)') &!
    'richmol_data/read_tran error: file', trim(fname), 'has bogus header = "', trim(sbuf), '" (expected "Start richmol format")'
    stop 'STOP, error in richmol_data/read_tran'
  endif

  read(iounit,*,iostat=info) tens_name_, nirrep, ncart
  iline = iline + 1
  call check_iostat

  if (trim(tens_name_)/=trim(tens_name)) then
    write(out, '(/a,1x,a,1x,a,1x,a)') &!
    'richmol_data/read_tran error: input operator name and name in file do not agree:', trim(tens_name), '.ne.', trim(tens_name_)
    stop 'STOP, error in richmol_data/read_tran'
  endif

  trans%tens_name = trim(tens_name)
  trans%nirrep = nirrep
  trans%nelem = ncart

  if (verbose) then
    write(out, '(1x,a,1x,i3)') 'no. irreps:', trans%nirrep
    write(out, '(1x,a,1x,i3)') 'no. elements:', trans%nelem
  endif

  ! read M-tensor

  if (verbose) write(out, '(1x,a)') 'read M-tensor'

  if (allocated(trans%mmat)) deallocate(trans%mmat)
  allocate(trans%mmat(ncart), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'richmol_data/read_tran error: failed to allocate trans%mmat(ncart)', 'ncart =', ncart
    stop 'STOP, error in richmol_data/read_tran'
  endif

  read(iounit,'(a)',iostat=info) sbuf
  iline = iline + 1
  call check_iostat
  if (trim(adjustl(sbuf))/='M-tensor') then
    write(out, '(/a,1x,a,1x,a,1x,a,1x,a,1x,i8)') &!
    'richmol_data/read_tran error: file', trim(fname), 'has bogus M-block header = "', trim(sbuf), &
    '" (expected "M-tensor"), line =', iline
    stop 'STOP, error in richmol_data/read_tran'
  endif

  do ielem=1, ncart

    read(iounit,*,iostat=info) sbuf, ielem_, icmplx, cart_label
    iline = iline + 1
    call check_iostat

    if (trim(adjustl(sbuf))/='alpha') then
      write(out, '(/a,a,a/a,1x,a,1x,a,1x,i8)') 'richmol_data/read_tran error: sbuf = "', trim(sbuf), '", while expected "alpha"', &!
      'file =', trim(fname), 'line =', iline
      stop 'STOP, error in richmol_data/read_tran'
    endif

    if (ielem_>ncart) then
      write(out, '(/a,1x,i3,1x,a,1x,i3/a,1x,a,1x,a,1x,i8)') &!
      'richmol_data/read_tran error: Cartesian component "ielem" = ', ielem_, '> total number of Cartesian components "ncart" =',&
       ncart,'file =', trim(fname), 'line =', iline
      stop 'STOP, error in richmol_data/read_tran'
    endif

    trans%icmplx = icmplx
    if (icmplx==-1) then
      trans%cmplx_fac = cmplx(0.0_rk,1.0_rk)
    elseif(icmplx==0) then
      trans%cmplx_fac = cmplx(1.0_rk,0.0_rk)
    else
      write(out,'(/a,1x,i3/a,1x,a,1x,a,1x,i8)') 'richmol_data/read_tran error: illegal imaginary/real number index "icmplx" =',&
           icmplx,'file =', trim(fname), 'line =', iline
      stop 'STOP, error in richmol_data/read_tran'
    endif

    nmpairs_max = nmpairs0
    allocate( trans%mmat(ielem_)%mpair(2,nmpairs_max), trans%mmat(ielem_)%mpair_ind(-j1:j1,-j2:j2),&
              trans%mmat(ielem_)%me(nirrep,nmpairs_max), stat=info)
    if (info/=0) then
      write(out, '(/a,a/a,10(1x,i6))') &!
      'richmol_data/read_tran error: failed to allocate trans%mmat(ielem_)%mpair(2,nmpairs_max),',&
      ' trans%mmat(ielem_)%mpair_ind(-j1:j1,-j2:j2), trans%mmat(ielem_)%me(nirrep,nmpairs_max)', &!
      'ielem_, nmpairs_max, nirrep, j1, j2 =', ielem_, nmpairs_max, nirrep, j1, j2
      stop 'STOP, error in richmol_data/read_tran'
    endif
    trans%mmat(ielem_)%mpair = 0
    trans%mmat(ielem_)%mpair_ind = 0
    trans%mmat(ielem_)%me = 0

    impair = 0
    do while(.true.)
      read(iounit,*,iostat=info) sbuf
      iline = iline + 1
      call check_iostat
      if (trim(sbuf)=='alpha'.or.trim(sbuf)=='K-tensor') then
        backspace(iounit)
        iline = iline - 1
        exit
      endif
      backspace(iounit)
      read(iounit,*,iostat=info) m1_, m2_, me(1:nirrep)
      call check_iostat
      if (transp) then
        m1 = m2_
        m2 = m1_
        if (icmplx==-1) me(1:nirrep) = -me(1:nirrep)
      else
        m1 = m1_
        m2 = m2_
      endif
      impair = impair + 1
      if (impair>nmpairs_max) call extend_mmat
      trans%mmat(ielem_)%mpair(1:2,impair) = (/m1,m2/)
      trans%mmat(ielem_)%mpair_ind(m1,m2) = impair
      trans%mmat(ielem_)%me(1:nirrep,impair) = me(1:nirrep)
    enddo

    trans%mmat(ielem_)%nmpairs = impair

    if (verbose) then
      write(out, '(1x,i3,1x,a10,5x,a,1x,i4,5x,a,2(1x,f6.3),1x,a)') &!
      ielem_, trim(cart_label), 'no. m-pairs:', trans%mmat(ielem_)%nmpairs, 'coef = (', trans%cmplx_fac, 'i )'
    endif

  enddo ! ielem

  ! read K-tensor

  if (verbose) write(out, '(1x,a)') 'read K-tensor'

  nkpairs_max = nkpairs0
  if (allocated(trans%kpair)) deallocate(trans%kpair)
  if (allocated(trans%kme)) deallocate(trans%kme)
  allocate( trans%kpair(2,nkpairs_max), trans%kme(nirrep,nkpairs_max), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') &!
    'richmol_data/read_tran error: failed to allocate trans%kpair(2,nkpairs_max), trans%kme(nirrep,nkpairs_max)', &!
    'nkpairs_max, nirrep =', nkpairs_max, nirrep
    stop 'STOP, error in richmol_data/read_tran'
  endif
  trans%kpair = 0
  trans%kme = 0


  read(iounit,'(a)') sbuf
  iline = iline + 1
  call check_iostat
  if (trim(adjustl(sbuf))/='K-tensor') then
    write(out, '(/a,1x,a,1x,a,1x,a,1x,a,1x,i8)') &!
    'richmol_data/read_tran error: file', trim(fname), 'has bogus K-block header = "', trim(sbuf), &
    '" (expected "K-tensor"), line =', iline
    stop 'STOP, error in richmol_data/read_tran'
  endif

  ikpair = 0

  do while(.true.)

    read(iounit,'(a)',iostat=info) sbuf1000
    iline = iline + 1
    call check_iostat

    if (trim(adjustl(sbuf1000))=='End richmol format') exit

    read(sbuf1000,*,iostat=info) id1_, id2_, ideg1_, ideg2_, me(1:nirrep)
    call check_iostat

    if (transp) then
      id1 = id2_
      id2 = id1_
      ideg1 = ideg2_
      ideg2 = ideg1_
    else
      id1 = id1_
      id2 = id2_
      ideg1 = ideg1_
      ideg2 = ideg2_
    endif

    istate1 = 0
    do istate=1, states(j1)%nlevels
      if (states(j1)%id(istate)==id1.and.states(j1)%ndeg(istate)>=ideg1) then
        istate1 = istate
        exit
      endif
    enddo

    istate2 = 0
    do istate=1, states(j2)%nlevels
      if (states(j2)%id(istate)==id2.and.states(j2)%ndeg(istate)>=ideg2) then
        istate2 = istate
        exit
      endif
    enddo

    if (istate1>0.and.istate2>0) then
      if (any(abs(me(1:nirrep))>=me_thresh)) then
        ikpair = ikpair + 1
        if (ikpair>nkpairs_max) call extend_kmat
        trans%kpair(1:2,ikpair) = (/states(j1)%istate(istate1,ideg1),states(j2)%istate(istate2,ideg2)/)
        trans%kme(1:nirrep,ikpair) = me(1:nirrep)
      endif
    endif

  enddo

  close(iounit)
  call IOStop(ioname)


  trans%nkpairs = ikpair

  if (verbose) write(out, '(1x,a,1x,i8)') 'number of rovibrational state pairs:', trans%nkpairs

  if (trans%nkpairs==0)  trans%null = .true.

  if (verbose) write(out, '(100(1h.),1x,a)') 'richmol_data/read_tran done'


  contains

  subroutine check_iostat
    if (info/=0) then
      write(out, '(/a,1x,a,1x,a,1x,i8)') 'richmol_data/read_tran error while reading file', trim(fname), 'line =', iline
      stop 'STOP, error in richmol_data/read_tran'
    endif
  end subroutine check_iostat

  subroutine extend_mmat
    integer(ik) :: nmpairs_max_
    integer(ik), allocatable :: itmp(:,:)
    real(rk), allocatable :: rtmp(:,:)
    call move_alloc(trans%mmat(ielem_)%mpair,itmp)
    call move_alloc(trans%mmat(ielem_)%me,rtmp)
    nmpairs_max_ = nmpairs_max
    nmpairs_max = nmpairs_max +  nmpairs_incr
    allocate( trans%mmat(ielem_)%mpair(2,nmpairs_max), trans%mmat(ielem_)%me(nirrep,nmpairs_max), stat=info)
    if (info/=0) then
      write(out, '(/a/a,10(1x,i6))') 'richmol_data/read_tran error: failed to allocate trans%mmat(ielem_)%mpair(2,nmpairs_max), &
                 trans%mmat(ielem_)%me(nirrep,nmpairs_max)','ielem_, nmpairs_max, nirrep =', ielem_, nmpairs_max, nirrep
      stop 'STOP, error in richmol_data/read_tran'
    endif
    trans%mmat(ielem_)%mpair = 0
    trans%mmat(ielem_)%me = 0
    trans%mmat(ielem_)%mpair(1:2,1:nmpairs_max_) = itmp
    trans%mmat(ielem_)%me(1:nirrep,1:nmpairs_max_) = rtmp
    deallocate(itmp,rtmp)
  end subroutine extend_mmat

  subroutine extend_kmat
    integer(ik) :: nkpairs_max_
    integer(ik), allocatable :: itmp(:,:)
    real(rk), allocatable :: rtmp(:,:)
    call move_alloc(trans%kpair,itmp)
    call move_alloc(trans%kme,rtmp)
    nkpairs_max_ = nkpairs_max
    nkpairs_max = nkpairs_max +  nkpairs_incr
    allocate( trans%kpair(2,nkpairs_max), trans%kme(nirrep,nkpairs_max), stat=info)
    if (info/=0) then
      write(out, '(/a,a/a,10(1x,i6))') 'richmol_data/read_tran error: failed to allocate trans%kpair(2,nkpairs_max),',&
                 ' trans%kme(nirrep,nkpairs_max)', &!
                 'nkpairs_max, nirrep =', nkpairs_max, nirrep
      stop 'STOP, error in richmol_data/read_tran'
    endif
    trans%kpair = 0
    trans%kme = 0
    trans%kpair(1:2,1:nkpairs_max_) = itmp
    trans%kme(1:nirrep,1:nkpairs_max_) = rtmp
    deallocate(itmp,rtmp)
  end subroutine extend_kmat

end subroutine read_tran



!###################################################################################################################################


! Reads rovibrational states from the RichMol database
! fname - database file name
! filter - state_filter_type object describing the rovibrational state filters
! states - states_type object containing the information about the rovibrational states, i.e. state ID numbers, energies, 
!          quantum numbers, etc.

subroutine read_states(fname, filter, states)

  character(*), intent(in) :: fname
  type(state_filter_type), intent(in) :: filter
  type(states_type), allocatable, intent(out) :: states(:)

  integer(ik), parameter :: maxnsym = 100
  integer(ik) :: jmin, jmax, info, iounit, iline, j, id, ndeg, k_quanta, v_quanta(1000), nmodes, nlevels, ilevel, maxndeg, ideg, &!
                 istate(filter%jmin:filter%jmax), isym
  real(rk) :: energy
  character(cl) :: sym, sym_vib, sym_rot, sym_arr(maxnsym)=''

  write(out, '(100(1h.),1x,a)') 'richmol_data/read_states start'

  nmodes = filter%nmodes
  if (nmodes<=0) then
    write(out, '(/a)') 'richmol_data/read_states error: filter%nmodes <= 0'
    stop 'STOP, error in richmol_data/read_states'
  endif

  jmin = filter%jmin
  jmax = filter%jmax

  ! allocate object "states"

  if (allocated(states)) deallocate(states)
  allocate(states(jmin:jmax), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'richmol_data/read_states error: failed to allocate states(jmin:jmax)','jmin, jmax =',jmin,jmax
    stop 'STOP, error in richmol_data/read_states'
  endif

  ! open RichMol states file

  call IOstart(fname, iounit)

  write(out, '(1x,a,1x,a)') 'read RichMol states file:', trim(fname)

  open(unit=iounit,form='formatted',action='read',position='rewind',status='old',file=trim(fname),iostat=info)
  if (info/=0) then
    write(out, '(/a,1x,a)') 'richmol_data/read_states error: failed to open file', trim(fname)
    stop 'STOP, error in richmol_data/read_states'
  endif

  ! estimate total number of states for each J quantum number

  states(:)%nlevels = 0
  iline = 0
  maxndeg = 0
  rewind(iounit)
  do
    read(iounit,*,iostat=info) j, id, sym, isym, ndeg, energy
    iline = iline + 1
    if (info>0) then ! reading error
      write(out, '(/a,1x,a,1x,a,1x,i4)') 'richmol_data/read_states error while reading first five columns in file', trim(fname),&
                   ', line =', iline
      stop 'STOP, error in richmol_data/read_states'
    elseif (info<0) then ! EOF
      exit
    endif
    if (.not.state_filter(filter, j=j, energy=energy)) cycle ! apply states filter
    states(j)%nlevels = states(j)%nlevels + 1
    maxndeg = max(maxndeg,ndeg)
  enddo

  write(out, '(/5x,a,2x,a)') 'J', 'no.states'

  do j=jmin, jmax
    nlevels = states(j)%nlevels
    write(out, '(3x,i3,5x,i6)') j, nlevels
    if (nlevels==0) then
      write(out, '(/a,1x,i4)') 'richmol_data/read_states error: no states found for J =', j
      stop 'STOP, error in richmol_data/read_states'
    endif
  enddo

   ! allocate "states(j)%id(:)", "states(j)%ndeg(:), etc. arrays

  do j=jmin, jmax
    nlevels = states(j)%nlevels
    allocate( states(j)%id(nlevels), states(j)%ndeg(nlevels), states(j)%sym(nlevels), states(j)%isym(nlevels), &
              states(j)%energy(nlevels), &!
              states(j)%k_quanta(nlevels), states(j)%v_quanta(nmodes,nlevels), states(j)%sym_rot(nlevels), &!
              states(j)%sym_vib(nlevels), states(j)%istate(nlevels,maxndeg), stat=info )
    if (info/=0) then
      write(out, '(/a/a,10(1x,i6))') &!
      'richmol_data/read_states error: failed to allocate states(j)%id(nlevels), states(j)%ndeg(nlevels), ...', &!
      'j, nlevels, nmodes =', j, nlevels, nmodes
      stop 'STOP, error in richmol_data/read_states'
    endif
    states(j)%id = 0
    states(j)%ndeg = 0
    states(j)%sym = ''
    states(j)%isym = 0
    states(j)%energy = 0.0
    states(j)%k_quanta = 0
    states(j)%v_quanta = 0
    states(j)%sym_rot = ''
    states(j)%sym_vib = ''
    states(j)%istate(:,:) = 0
  enddo

  ! read energies, IDs and quantum numbers

  states(:)%nlevels = 0
  iline = 0
  istate(:) = 0
  rewind(iounit)
  do
    read(iounit,*,iostat=info) j, id, sym, isym, ndeg, energy, k_quanta, v_quanta(1:nmodes), sym_rot, sym_vib
    iline = iline + 1
    if (info>0) then ! reading error
      write(out, '(/a,1x,a,1x,a,1x,i4)') 'richmol_data/read_states error while reading file', trim(fname), ', line =', iline
      stop 'STOP, error in richmol_data/read_states'
    elseif (info<0) then ! EOF
      exit
    endif
    if (.not.state_filter(filter, j=j, energy=energy)) cycle ! apply states filter
    states(j)%nlevels = states(j)%nlevels + 1
    ilevel = states(j)%nlevels
    states(j)%id(ilevel) = id
    states(j)%ndeg(ilevel) = ndeg
    states(j)%sym(ilevel) = sym
    states(j)%isym(ilevel) = isym
    states(j)%energy(ilevel) = energy
    states(j)%k_quanta(ilevel) = k_quanta
    states(j)%v_quanta(1:nmodes,ilevel) = v_quanta(1:nmodes)
    states(j)%sym_rot(ilevel) = sym_rot
    states(j)%sym_vib(ilevel) = sym_vib
    do ideg=1, ndeg
      istate(j) = istate(j) + 1
      states(j)%istate(ilevel,ideg) = istate(j)
    enddo
    !
  enddo

   write(out, '(/1x,a,1x,i6/1x,a,1x,i6)') 'total number of states:', sum(states(:)%nlevels), &!
   'those including degenerate:', sum( (/(sum(states(j)%ndeg(:)), j=jmin, jmax)/) )

  close(iounit)
  call IOstop(fname)

  write(out, '(100(1h.),1x,a)') 'richmol_data/read_states done'

end subroutine read_states



!###################################################################################################################################



function state_filter(filter, j, sym, energy, k, sym_rot, qvib, sym_vib) result(f)

  type(state_filter_type), intent(in) :: filter
  integer(ik), intent(in), optional :: j, k, qvib(:)
  real(rk), intent(in), optional :: energy
  character(cl), intent(in), optional :: sym, sym_rot, sym_vib
  logical :: f

  logical :: passed_j, passed_e

  passed_j = .false.
  passed_e = .false.

  if (present(j)) then
    if (j>=filter%jmin .and. j<=filter%jmax) passed_j = .true.
  endif

  if (present(energy)) then
    if (energy>=filter%emin .and. energy<=filter%emax) passed_e = .true.
  endif

  ! other filters ... to be added

  f = passed_j .and. passed_e

end function state_filter


end module richmol_data
