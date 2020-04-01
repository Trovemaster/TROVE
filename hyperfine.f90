module hyperfine

#define fwig  0
!
use accuracy
use timer
use richmol_data
#if (fwig>0)
  use fwigxjpf
#endif
use lapack
implicit none


real(rk), parameter :: VQ_mbau_cm = 7.83758020955394d-006 ! to convert from V[au]*Q[mb] to V*Q[cm^-1]
                                                          !   where V[au] is electric field gradient in atomic units
                                                          !   and Q[mb] is nuclear quadrupole constant in units millibarn (1mb=1E-31meter^2)


type hyperfine_type
  character(cl)               ::  task = 'N/A'
  real(rk)                    ::  ftot = 0
  integer(ik)                 ::  nspin = 0
  integer(ik)                 ::  nsym
  real(rk), allocatable       ::  spin(:)
  real(rk), allocatable       ::  gns(:)
  real(rk), allocatable       ::  quad_const(:)
  character(cl), allocatable  ::  efg_tens_name(:)
  character(cl)               ::  dipole_tens_name = 'MU'
  type(state_filter_type)     ::  filter
  character(cl)               ::  states_file = 'N/A'
  real(rk)                    ::  me_thresh = 1.0d-14
  real(rk)                    ::  eigenval_symtol_khz = 0.1d0
  real(rk)                    ::  eigencoef_symtol = 1.0d-08
  real(rk)                    ::  eigencoef_thresh = 1.0d-14
  character(cl)               ::  eigenfile_quanta = 'hyper_eigen_quanta'
  character(cl)               ::  eigenfile_descr = 'hyper_eigen_descr'
  character(cl)               ::  eigenfile_vectors = 'hyper_eigen_vectors'
  real(rk)                    ::  zpe, part_func, temperature
  character(cl)               ::  specfile = 'hyper_intens'
end type hyperfine_type


type, extends(states_type) :: hyper_states_type
  integer(ik), allocatable :: jval(:), ideg(:), ilevel(:)
  real(rk), allocatable    :: enr_bas(:)
end type hyper_states_type


contains



subroutine hyper_hmat(hyper)

  type(hyperfine_type), intent(in) :: hyper

  integer(ik), parameter :: max_coupl_dimen=1000, max_nirrep=100, maxncoef_print=10
  integer(ik) :: nmom, coupl_dimen, ielem, nelem, i, max_ang, nspin, jmin, jmax, coupl_qrot(max_coupl_dimen), j, idimen, nlevels, &!
                 ilevel, istate, ideg, nquad, j1, j2, jval, iquad, ialpha, impair, nirrep, ikpair, istate1, istate2, jelem, jdimen, &!
                 maxnapp, maxnstates, info, hmat_dimen, rank, icoupl, jcoupl, ielem_print(maxncoef_print), ncoef_print, nmodes, id, &!
                 iounit, nnz, nelem0, job(6), idimen0, isym1, isym2, j0, ilevel0, ideg0, isym0, ndeg, ndeg0
  integer(ik), allocatable :: jmap_nelem(:,:), jmap_ind(:,:,:), imap_ind(:), jmap_ind_inv(:,:), hmat_csr_ja(:), hmat_csr_ia(:)
  real(rk) :: small, qmom(hyper%nspin+1), fval, coupl(hyper%nspin+1), qmom_coupl(hyper%nspin+1,max_coupl_dimen), print_coef_thresh, &!
              coupl_qspin(hyper%nspin,max_coupl_dimen), coupl_qf(max_coupl_dimen), mmat(max_nirrep), me, small_, efg_coef
  real(rk), allocatable :: rme_efg_coef(:,:), rme_quad(:,:,:), hmat(:,:), hmat_diag_khz(:), enr(:), tmp_vec(:), hmat_acsr(:)
  character(cl) :: si(hyper%nspin), oper, fname, sfval, sym1, sym2, sym0
  type(trans_type) :: trans
  type(states_type), allocatable :: states(:)

  write(out, '(100(1h.),1x,a)') 'hyperfine/hyper_hmat start'

  write(out, '(///10x,a)') '---------------------------------------'
  write(out, '(10x,a)')    'Build hyperfine interaction Hamiltonian'
  write(out, '(10x,a)')    '---------------------------------------'


  small = 1.0d-06                  ! threshold used for checking whether two half integer numbers agree (e.g. spin and total angular momentum quantum numbers)
  small_ = epsilon(1.0d0)*100.0d0  ! threshold for zero
  print_coef_thresh = 0.2_rk       ! threshold for basis function coefficient^2 to print the corresponding basis quantum numbers in the state assignment line
  nmodes = hyper%filter%nmodes     ! number of vibrational coordinates


  fval = hyper%ftot         ! value of the total (rotations plus spin) angular momentum F
  jmin = hyper%filter%jmin  ! min value of the rotational angular momentum J
  jmax = hyper%filter%jmax  ! max value of J

  nspin = hyper%nspin  ! number of nuclear spins
  nmom = nspin + 1     ! number of coupled angular momenta = number of nuclear spins + rotation

  qmom(1:nmom-1) = (/hyper%spin(1:nspin)/)  ! nuclear spin quanta (I={I_1,I_2,..I_nspin})


  ! for given quantum number of the total angular momentum (fval)
  ! generate all possible valid combinations of the nuclear spin and rotational quanta (jval=jmin:jmax)

  coupl_dimen = 0 ! number of valid I+J combinations

  do jval=jmin, jmax

    qmom(nmom) = real(jval,rk)

    call coupl_ang_quanta(1, nmom, qmom(1:nmom), coupl(1:nmom), nelem, qmom_coupl)

    do ielem=1, nelem
      if (abs(qmom_coupl(nmom,ielem)-fval)<=small) then
        coupl_dimen = coupl_dimen + 1
        if (coupl_dimen>max_coupl_dimen) then
          write(out, '(/a,1x,i6)') 'hyperfine/hyper_hmat error: "coupl_dimen" exceeds "max_coupl_dimen" =', max_coupl_dimen
          stop 'STOP, error in hyperfine/hyper_hmat'
        endif
        coupl_qspin(1:nspin,coupl_dimen) = qmom_coupl(1:nmom-1,ielem)
        coupl_qrot(coupl_dimen) = jval
        coupl_qf(coupl_dimen) = qmom_coupl(nmom,ielem)
      endif
    enddo

  enddo

  ! print valid I+J combinations

  write(out, '(/1x,a,1x,i3,1x,a)') 'coupled nuclear spin (I) and rotational (J) quantum numbers (', coupl_dimen, '):'

  do i=1, nspin
    write(si(i),*) i
  enddo

  write(out, '(2(/8x,a,<nspin>(3x,a),6x,a))') 'I_1', ('I_1'//trim(adjustl(si(i))), i=2, nspin), '   J', 'F',  '---', ('----', i=2, nspin), '   -', '-'
  do icoupl=1, coupl_dimen
    write(out, '(1x,i3,100(1x,f6.1))') icoupl, coupl_qspin(1:nspin,icoupl), real(coupl_qrot(icoupl),rk), coupl_qf(icoupl)
  enddo



  ! read rovibrational states file (RichMol database format)

  call read_states(hyper%states_file, hyper%filter, states)



  ! establish mapping between the global spin-rovibrational basis indices and the pure rovibrational state indices
  ! estimate the total dimension of the Hamiltonian matrix

  jmin = minval(coupl_qrot(1:coupl_dimen))
  jmax = maxval(coupl_qrot(1:coupl_dimen))

  write(out, '(/1x,a,1x,i3,1x,i3)') 'min and max spin-coupled J quanta:', jmin, jmax

  maxnstates = maxval( (/(sum(states(j)%ndeg(:)), j=jmin, jmax)/) )
  maxnapp = maxval( (/(count(coupl_qrot(1:coupl_dimen)==j), j=jmin, jmax)/) )

  write(out, '(/1x,a,1x,i6)') 'max number of rovibrational states (within J-block):', maxnstates
  write(out, '(1x,a,1x,i3)') 'max number of appearences of same J in different spin-rovibration coupling blocks:', maxnapp

  allocate( jmap_nelem(maxnstates,jmin:jmax), jmap_ind(maxnapp,maxnstates,jmin:jmax), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') &!
    'hyperfine/hyper_hmat error: failed to allocate jmap_nelem(maxnstates,jmin:jmax), jmap_ind(maxnapp,maxnstates,jmin:jmax)', &!
    'maxnstates, maxnapp, jmin, jmax =', maxnstates, maxnapp, jmin, jmax
    stop 'STOP, error in hyperfine/hyper_hmat'
  endif

  jmap_nelem = 0
  jmap_ind = 0

  idimen = 0
  do icoupl=1, coupl_dimen
    j = coupl_qrot(icoupl)
    do ilevel=1, states(j)%nlevels
      do ideg=1, states(j)%ndeg(ilevel)
      	istate = states(j)%istate(ilevel,ideg)
        idimen = idimen + 1
        jmap_nelem(istate,j) = jmap_nelem(istate,j) + 1
        nelem = jmap_nelem(istate,j)
        jmap_ind(nelem,istate,j) = idimen
      enddo
    enddo
  enddo

  hmat_dimen = idimen

  write(out, '(/1x,a,1x,i6)') 'total dimension of the hyperfine Hamiltonian matrix:', hmat_dimen



  ! establish mapping between global spin-rovibration basis indices and nuclear spin basis indices

  allocate( imap_ind(hmat_dimen), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'hyperfine/hyper_hmat error: failed to allocate imap_ind(hmat_dimen)', 'hmat_dimen =', hmat_dimen
    stop 'STOP, error in hyperfine/hyper_hmat'
  endif

  imap_ind = 0

  idimen = 0
  do icoupl=1, coupl_dimen
    j = coupl_qrot(icoupl)
    do ilevel=1, states(j)%nlevels
      do ideg=1, states(j)%ndeg(ilevel)
        idimen = idimen + 1
        imap_ind(idimen) = icoupl
      enddo
    enddo
  enddo



  ! establish mapping between global spin-rovibration basis indices and rovibrational basis indices

  allocate( jmap_ind_inv(3,hmat_dimen), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'hyperfine/hyper_hmat error: failed to allocate jmap_ind_inv(3,hmat_dimen)', 'hmat_dimen =', hmat_dimen
    stop 'STOP, error in hyperfine/hyper_hmat'
  endif

  jmap_ind_inv = 0

  idimen = 0
  do icoupl=1, coupl_dimen
    j = coupl_qrot(icoupl)
    do ilevel=1, states(j)%nlevels
      do ideg=1, states(j)%ndeg(ilevel)
        idimen = idimen + 1
        jmap_ind_inv(1:3,idimen) = (/j,ilevel,ideg/)
      enddo
    enddo
  enddo



  ! init module "fwigxjpf" for computing 3j, 6j, and 9j symbols

  max_ang = nint( maxval(coupl_qspin(1:nspin,1:coupl_dimen)) + maxval(coupl_qrot(1:coupl_dimen)) )

  write(out, '(/1x,a,1x,i2,a,1x,i3)') 'initialize "fwigxjpf" for max', 6, '-j symbols and max angular momentum =', max_ang


#if (fwig>0)
  call fwig_table_init( int(2*max_ang,kind=4), int(6,kind=4) )
  call fwig_temp_init( int(2*max_ang,kind=4) )
#endif



  ! precompute reduced matrix elemens of the nuclear quadrupole tensor operator <I_12',..I_1N'||Q(iquad)||I_12,...I_1N>
  ! and multiply them with the coefficient (-1)**(J+I_1N'+F)*{ F  I_1N'  J'  }
  !                                                          { 2   J    I_1N }

  write(out, '(/1x,a)') 'compute reduced matrix elements of nuclear quadrupole'

  allocate( rme_quad(coupl_dimen,coupl_dimen,nspin), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'hyperfine/hyper_hmat error: failed to allocate rme_quad(coupl_dimen,coupl_dimen,nspin)', &!
    'coupl_dimen, nspin =', coupl_dimen, nspin
    stop 'STOP, error in hyperfine/hyper_hmat'
  endif

  rme_quad = 0
  rank = 2
  oper = 'QUAD'

  write(out, '(2(/1x,a,1x,a,1x,a,13x,a))') 'iquad', 'ibra', 'iket', '<eQ>', '-----', '----', '----', '----'

  do iquad=1, nspin
    do ielem=1, coupl_dimen
      do jelem=1, coupl_dimen

        rme_quad(jelem,ielem,iquad) = nucspin_reduced_me(hyper, nspin, iquad, coupl_qspin(1:nspin,jelem), coupl_qspin(1:nspin,ielem), rank, oper)

        ! print reduced matrix element

        if (jelem>=ielem) then
          write(out, '(3x,i3,2x,i3,2x,i3,1x,es16.8,3x,a,<nspin>(1x,f4.1),a,i2,a,<nspin>(1x,f4.1),a)') &!
          iquad, jelem,ielem, rme_quad(jelem,ielem,iquad), &!
          '= <', coupl_qspin(1:nspin,jelem), ' || Q(', iquad, ') || ', coupl_qspin(1:nspin,ielem), ' >'
        endif

        ! multiply reduced matrix element with 6j coefficient

        rme_quad(jelem,ielem,iquad) = rme_quad(jelem,ielem,iquad) * (-1.0_rk)**(fval+real(coupl_qrot(ielem),rk)+coupl_qspin(nspin,jelem)) &!
                                    * sixj_symbol( fval,   coupl_qspin(nspin,jelem),   real(coupl_qrot(jelem),rk), &!
                                                   2.0_rk, real(coupl_qrot(ielem),rk), coupl_qspin(nspin,ielem) )

      enddo
    enddo
  enddo



  ! precompute coefficients coef(J',J) = 2*(-1)**(J'-J)*( J'  2  J )
  !                                                     ( -J  0  J )
  ! these are needed to compute reduced matrix elements of electric field gradient (EFG) tensor as:
  ! <J'||V(iquad)||J> = <J',m=J|V_ZZ(iquad)|J,m=J> / coef(J',J)

  allocate(rme_efg_coef(jmin:jmax,jmin:jmax), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'hyperfine/hyper_hmat error: failed to allocate rme_efg_coef(jmin:jmax,jmin:jmax)', &!
    'jmin, jmax =', jmin, jmax
    stop 'STOP, error in hyperfine/hyper_hmat'
  endif

  rme_efg_coef = 0

  do j1=jmin, jmax
    do j2=jmin, jmax
      rme_efg_coef(j1,j2) = threej_symbol(real(j2,rk),real(-j1,rk), 2.0_rk,0.0_rk, real(j1,rk),real(j1,rk)) / (0.5_rk*(-1)**(j1-j2))
    enddo
  enddo


  ! start building the total hyperfine Hamiltonian matrix

  allocate(hmat(hmat_dimen,hmat_dimen), hmat_diag_khz(hmat_dimen), enr(hmat_dimen), tmp_vec(hmat_dimen), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') &!
    'hyperfine/hyper_hmat error: failed to allocate hmat(hmat_dimen,hmat_dimen), hmat_diag_khz(hmat_dimen), enr(hmat_dimen)', &!
    'hmat_dimen =', hmat_dimen
    stop 'STOP, error in hyperfine/hyper_hmat'
  endif

  hmat = 0

  ialpha = 6 ! ZZ component of the second-rank electric field gradient (EFG) tensor

  do j1=jmin, jmax
  	do j2=jmin, jmax

      efg_coef = rme_efg_coef(j1,j2)
      if (abs(efg_coef)<=small_) cycle ! due to selection rules for the corresponding three-j symbol

      do iquad=1, nspin

        if (abs(hyper%quad_const(iquad))<=hyper%me_thresh) cycle ! skip for zero-quadrupole nuclei

        ! read j1/j2 rovibrational matrix elements (of the EFG tensor)
        !   these were previously computed and stored in the RichMol database format files
        !   named as matelem_<hyper%efg_tens_name(iquad)>_j<j1>_j<j2>.rchm

        call read_tran(j1, j2, hyper%filter%jmin, hyper%filter%jmax, trim(hyper%efg_tens_name(iquad)), states, hyper%me_thresh, trans)

        if (trans%null) cycle ! trans%null=.true. means that all transition matrix elements for given j1/j2 pair are zero (by rigorous selection rules)
                              ! or simply the file with matrix elements does not exist

        ! double check if we are reading matrix elements of a second-order tensor operator (EFG tensor)
        !  i.e. the number of Cartesian components is equal to six

        if (trans%nelem/=6) then
          write(out, '(/a,1x,i3,1x,i3,1x,a,1x,a/a,1x,i3,1x,a)') &!
          'hyperfine/hyper_hmat error: file with j1/j2=', j1, j2, 'rovibrational matrix elements of tensor operator', trim(trans%tens_name), &!
          'returns number of tensor components =', trans%nelem, ', which is not equal to six as expected for the EFG tensor'
          stop 'STOP, error in hyperfine/hyper_hmat'
        endif

        nirrep = trans%nirrep

        ! for the EFG matrix elements we need only <J',m'=J|Vzz|J,m=J> if J'>J and <J',m'=J'|Vzz|J,m=J'> if J'<J
        if (j2>=j1) then
          impair = trans%mmat(ialpha)%mpair_ind(j1,j1)
        else
          impair = trans%mmat(ialpha)%mpair_ind(j2,j2)
        endif
        mmat(1:nirrep) = trans%mmat(ialpha)%me(1:nirrep,impair)

        do ikpair=1, trans%nkpairs

          me = sum( trans%kme(1:nirrep,ikpair) * mmat(1:nirrep) )
          me = me / efg_coef

          istate1 = trans%kpair(1,ikpair)
          istate2 = trans%kpair(2,ikpair)

          do ielem=1, jmap_nelem(istate1,j1)
            idimen = jmap_ind(ielem,istate1,j1)
            icoupl = imap_ind(idimen)
            do jelem=1, jmap_nelem(istate2,j2)
              jdimen = jmap_ind(jelem,istate2,j2)
              jcoupl = imap_ind(jdimen)
              if (idimen>0.and.jdimen>0) then
                hmat(idimen,jdimen) = hmat(idimen,jdimen) + me * rme_quad(jcoupl,icoupl,iquad)
              endif
            enddo
          enddo
        enddo ! ikpair

      enddo ! iquad
    enddo ! j2
  enddo ! j1

  ! save diagonal (largest) hyperfine contributions

  do j=jmin, jmax
    do ilevel=1, states(j)%nlevels
      do ideg=1, states(j)%ndeg(ilevel)
        istate = states(j)%istate(ilevel,ideg)
        do ielem=1, jmap_nelem(istate,j)
          idimen = jmap_ind(ielem,istate,j)
          hmat_diag_khz(idimen) = hmat(idimen,idimen)*vellgt/1000.0d0 ! in kHz
        enddo
      enddo
    enddo
  enddo


  ! add pure rovibrational part to Hamiltonian

  do j=jmin, jmax
    do ilevel=1, states(j)%nlevels
      do ideg=1, states(j)%ndeg(ilevel)
        istate = states(j)%istate(ilevel,ideg)
        do ielem=1, jmap_nelem(istate,j)
          do jelem=1, jmap_nelem(istate,j)
            idimen = jmap_ind(ielem,istate,j)
            jdimen = jmap_ind(jelem,istate,j)
            hmat(idimen,jdimen) = hmat(idimen,jdimen) + states(j)%energy(ilevel)
          enddo
        enddo
      enddo
    enddo
  enddo


  ! diagonalise matrix

  write(out, '(/1x,a)') 'Diagonalize spin-rovibrational Hamiltonian matrix...'

  call  lapack_dsyev(hmat, enr)

  write(out, '(1x,a)') '...done'


  ! print energies and quantum numbers

  write(out, '(/1x,a)') 'Energies and quantum numbers'

  write(out, '(6x,a,6x,a,2x,a,11x,a,8x,a,3x,a,2x,a,1x,a,6x,a,2x,a,4x,a,3x,a,3x,a/30x,a,34x,a,10x,a,8x,a)') &!
  'F', '#', 'sym', 'energy', 'coef', 'J', 'sym', 'ideg', 'enr(bas)', 'energy-enr(bas)', 'H(diag)', 'k', 'vib quanta', &!
  '[cm-1]', '[cm-1]', '[kHz]', '[kHz]'

  do idimen=1, hmat_dimen

    tmp_vec = abs(hmat(1:hmat_dimen,idimen))
    ielem_print = 0
    ncoef_print = 0
    do i=1, maxncoef_print
      ielem_print(i) = maxloc(tmp_vec(1:hmat_dimen),dim=1)
      if (tmp_vec(ielem_print(i))**2<print_coef_thresh) then
        exit
      endif
      tmp_vec(ielem_print(i)) = 0
      ncoef_print = ncoef_print + 1
    enddo

    idimen0 = ielem_print(1)
    j0 = jmap_ind_inv(1,idimen0)
    ilevel0 = jmap_ind_inv(2,idimen0)
    sym0 = states(j0)%sym(ilevel0)

    write(out, '(1x,f6.1,1x,i6,1x,a4,1x,f16.8,1x, 100(3x,a,1x,f6.3,1x,i3,1x,a4,1x,i2,f16.8,2x,2(1x,f12.3),1x,i3,<nmodes>(1x,i3),1x,a) )') &!
    fval, idimen, trim(sym0), enr(idimen), &                                                           ! print assignments for leading basis contributions
    ('(', hmat(ielem_print(i),idimen), &                                                                                    ! coefficient
          jmap_ind_inv(1,ielem_print(i)), &                                                                                 ! J
          trim(states( jmap_ind_inv(1,ielem_print(i)) )%sym( jmap_ind_inv(2,ielem_print(i)) )), &                           ! total symmetry
          jmap_ind_inv(3,ielem_print(i)), &                                                                                 ! ideg
          states( jmap_ind_inv(1,ielem_print(i)) )%energy( jmap_ind_inv(2,ielem_print(i)) ), &                              ! energy
          (enr(idimen)-states( jmap_ind_inv(1,ielem_print(i)) )%energy( jmap_ind_inv(2,ielem_print(i)) ))*vellgt/1.0d+3, &! ! energy-energy(basis function)
          hmat_diag_khz(ielem_print(i)), &!                                                                                 ! diagonal contribution from basis function
          states( jmap_ind_inv(1,ielem_print(i)) )%k_quanta( jmap_ind_inv(2,ielem_print(i)) ), &                            ! k-quantum number
          states( jmap_ind_inv(1,ielem_print(i)) )%v_quanta(:, jmap_ind_inv(2,ielem_print(i)) ), &                          ! vibrational quanta
    ')', i=1, ncoef_print )

  enddo


  !---------------------------------------------------------------------!
  ! since we don't enforce symmetry for spin-rovibrational calculations !
  ! here we check if molecular symmetry is still preserved              !
  !---------------------------------------------------------------------!

  write(out, '(/1x,a)') 'check if symmetry of spin-rovibrational states is preserved ...'
  write(out, '(1x,a,1x,es16.8)') 'symmetry tolerace for eigenvector coefficient:', hyper%eigencoef_symtol

  do istate=1, hmat_dimen
    idimen0 = maxloc(abs(hmat(1:hmat_dimen,istate)),1)
    if (idimen0==0) then
      write(out, '(/a,1x,i8,1x,a)') 'hyperfine/hyper_hmat error: eigenvector for state no.', istate, 'is zero-valued'
      stop 'STOP, error in hyperfine/hyper_hmat'
    endif
    j = jmap_ind_inv(1,idimen0)
    ilevel = jmap_ind_inv(2,idimen0)
    isym1 = states(j)%isym(ilevel)
    sym1 = trim(states(j)%sym(ilevel))
    do idimen=1, hmat_dimen
      if (idimen==idimen0) cycle
      if (abs(hmat(idimen,istate))>hyper%eigencoef_symtol) then
        j = jmap_ind_inv(1,idimen)
        ilevel = jmap_ind_inv(2,idimen)
        isym2 = states(j)%isym(ilevel)
        sym2 = trim(states(j)%sym(ilevel))
        if (isym1/=isym2) then
          write(out, '(/a,1x,i8/2(3x,a,1x,a,1x,a,1x,es16.8))') 'hyperfine/hyper_hmat error: eigenvector is not symmetry-adapted for state no.', istate, &!
          'symmetry1 =', trim(sym1), 'coef1 =', hmat(idimen0,istate), 'symmetry2 =', trim(sym2), 'coef2 =', hmat(idimen,istate)
          !stop 'STOP, error in hyperfine/hyper_hmat'
        endif
      endif
    enddo
  enddo

  write(out, '(1x,a)') '...done'


  !---------------------------!
  ! store eigenvalues in file !
  !---------------------------!

  write(sfval,'(f6.1)') fval
  fname = trim(adjustl(hyper%eigenfile_descr))//'_'//trim(adjustl(sfval))//'.chk'
  call IOStart(fname,iounit)

  write(out, '(/1x,a,1x,a)') 'store eigenvalues in file:', trim(fname)

  open(iounit,file=fname,form='formatted',action='write',status='unknown',position='rewind')
  write(iounit,'(a)') 'Start hyperfine energies'
  write(iounit,*) hmat_dimen

  do idimen=1, hmat_dimen
    idimen0 = maxloc(abs(hmat(1:hmat_dimen,idimen)),1)
    j0 = jmap_ind_inv(1,idimen0)
    ilevel0 = jmap_ind_inv(2,idimen0)
    sym0 = states(j0)%sym(ilevel0)
    isym0 = states(j0)%isym(ilevel0)
    ideg0 = jmap_ind_inv(3,idimen0)
    ndeg0 = states(j0)%ndeg(ilevel0)

    write(iounit, '(1x,f6.1,1x,i6,1x,a4,1x,i2,1x,f20.12,5x,i3,1x,i5,1x,i2,1x,i2,1x,a4,1x,f16.8,2x,f12.3,3x,i3,2x,<nmodes>(1x,i3),2x,a4,1x,a4)') &!
    fval, idimen, trim(sym0), isym0, enr(idimen), j0, ilevel0, ideg0, ndeg0, trim(sym0), states(j0)%energy(ilevel0), (enr(idimen)-states(j0)%energy(ilevel0))*vellgt/1.0d+3, &!
    states(j0)%k_quanta(ilevel0), states(j0)%v_quanta(:,ilevel0), trim(states(j0)%sym_rot(ilevel0)), trim(states(j0)%sym_vib(ilevel0))
  enddo

  write(iounit,'(a)') 'End hyperfine energies'
  close(iounit)
  write(out, '(1x,a)') '...done'
  call IOStop(fname)


  !-----------------------------------------------!
  ! store spin-rovibrational eigenvector indexing !
  !-----------------------------------------------!

  write(sfval,'(f6.1)') fval
  fname = trim(adjustl(hyper%eigenfile_quanta))//'_'//trim(adjustl(sfval))//'.chk'
  call IOStart(fname,iounit)

  write(out, '(/1x,a,1x,a)') 'store spin-rovibrational eigenvector indexig in file:', trim(fname)

  open(iounit,file=fname,form='formatted',action='write',position='rewind',status='unknown')

  write(iounit,'(a)') 'Start hyperfine basis indexes'
  write(iounit,'(1x,f6.1)') fval
  write(iounit,'(1x,i4,1x,i4,1x,i4)') jmin, jmax, maxnapp
  write(iounit,'(1x,i8)') hmat_dimen
  do idimen=1, hmat_dimen
    icoupl = imap_ind(idimen)
    j = jmap_ind_inv(1,idimen)
    ilevel = jmap_ind_inv(2,idimen)
    ideg = jmap_ind_inv(3,idimen)
    istate = states(j)%istate(ilevel,ideg)
    write(iounit,'(1x,f6.1,3x,i4,3x,i6)') coupl_qspin(nspin,icoupl), j, istate ! I_1N, J, ISTATE
  enddo
  write(iounit,'(a)') 'End hyperfine basis indexes'
  close(iounit)
  write(out, '(1x,a)') '...done'
  call IOStop(fname)


  !----------------------------------------------------------------!
  ! convert spin-rovibrational eigenvectors into sparse CSR format !
  !----------------------------------------------------------------!

  write(out, '(/1x,a)') 'convert eigenvalues into CSR sparse format ...'

  nelem0 = 0
  do idimen=1, hmat_dimen
    do jdimen=1, hmat_dimen
      if (abs(hmat(jdimen,idimen))<=hyper%eigencoef_thresh) then
        hmat(jdimen,idimen) = 0
        nelem0 = nelem0 + 1
      endif
    enddo
  enddo
  nnz = hmat_dimen * hmat_dimen - nelem0

  write(out, '(1x,a,1x,i8,1x,a,1x,es16.8,1x,a)') 'number of nonzero elements in eigenvector matrix:', nnz, '(>', hyper%eigencoef_thresh, ')'

  allocate( hmat_acsr(nnz), hmat_csr_ja(nnz), hmat_csr_ia(hmat_dimen+1), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') &!
    'hyperfine/hyper_hmat error: failed to allocate  hmat_acsr(nnz), hmat_csr_ja(nnz), hmat_csr_ia(hmat_dimen+1)', &!
    'nnz, hmat_dimen =', nnz, hmat_dimen
    stop 'STOP, error in hyperfine/hyper_hmat'
  endif

  job(1:6) = (/0, 1, 1, 2, nnz, 1/)

  call mkl_ddnscsr(job, hmat_dimen, hmat_dimen, hmat(1:hmat_dimen,1:hmat_dimen), hmat_dimen, hmat_acsr, hmat_csr_ja, hmat_csr_ia, info)

  if (info/=0) then
    write(out, '(/a)') 'hyperfine/hyper_hmat error: dens-to-csr format conversion with mkl_ddnscsr has failed'
    stop 'STOP, error in hyperfine/hyper_hmat'
  endif

  write(out, '(1x,a)') '...done'


  !----------------------------!
  ! store eigenvectors in file !
  !----------------------------!

  fname = trim(adjustl(hyper%eigenfile_vectors))//'_'//trim(adjustl(sfval))//'.chk'
  call IOStart(fname,iounit)

  write(out, '(/1x,a,1x,a)') 'store CSR-compressed eigenvectors in file:', trim(fname)

  open(iounit,file=fname,form='unformatted',action='write',status='unknown',position='rewind')
  write(iounit) 'Start hyperfine vectors'
  write(iounit) hmat_dimen, nnz
  write(iounit) hmat_acsr
  write(iounit) hmat_csr_ja
  write(iounit) hmat_csr_ia
  write(iounit) 'End hyperfine vectors'
  close(iounit)
  write(out, '(1x,a)') '...done'
  call IOStop(fname)



  deallocate(jmap_nelem, jmap_ind, rme_efg_coef, rme_quad, imap_ind, hmat, hmat_diag_khz, enr, tmp_vec, jmap_ind_inv)
  deallocate(hmat_acsr, hmat_csr_ja, hmat_csr_ia)

  write(out, '(100(1h.),1x,a)') 'hyperfine/hyper_hmat done'

end subroutine hyper_hmat



!###################################################################################################################################



subroutine hyper_read_quanta(hyper, fval, jmin, jmax, dimen, ind, jmap_nelem, jmap_ind, verbose_)

  type(hyperfine_type), intent(in) :: hyper
  real(rk), intent(in) :: fval
  integer(ik), intent(out) :: jmin, jmax, dimen
  integer(ik), allocatable, intent(out) :: ind(:,:), jmap_nelem(:,:), jmap_ind(:,:,:)
  logical, optional, intent(in) :: verbose_

  integer(ik) :: iounit, idimen, info, maxnapp, maxnstates, j, istate, nelem
  real(rk) :: small, fval_
  character(cl) :: fname, sfval, sbuf
  logical :: verbose=.false.

  if (present(verbose_)) verbose=verbose_

  if (verbose) write(out, '(100(1h.),1x,a)') 'hyperfine/hyper_read_quanta start'

  small = 1.0d-06 ! threshold for two half integer numbers to be equal

  write(sfval,'(f6.1)') fval
  fname = trim(adjustl(hyper%eigenfile_quanta))//'_'//trim(adjustl(sfval))//'.chk'
  call IOStart(fname,iounit)

  if (verbose) write(out, '(1x,a,1x,a)') 'read checkpoint file:', trim(fname)

  open(iounit,file=fname,form='formatted',action='read',position='rewind',status='old',iostat=info)
  if (info/=0) then
    write(out, '(/a,1x,a)') 'hyperfine/hyper_read_quanta error while opening file:', trim(fname)
    stop 'STOP, error in hyperfine/hyper_read_quanta'
  endif

  read(iounit,'(a)') sbuf

  if (trim(sbuf)/='Start hyperfine basis indexes') then
    write(out, '(/a,1x,a,1x,a,a,a)') &!
    'hyperfine/hyper_read_quanta error: file', trim(fname), 'has bogus header = "', trim(sbuf), '", (expected "Start hyperfine basis indexes")'
    stop 'STOP, error in hyperfine/hyper_read_quanta'
  endif

  read(iounit,*) fval_

  if (abs(fval_-fval)>small) then
    write(out, '(/a,1x,f6.1,1x,a,1x,f6.1,1x,a)') &!
    'hyperfine/hyper_read_quanta error: value of the quantum number F =', fval_, ', read from file, does not agree with F =', fval, 'used in the file naming'
    stop 'STOP, error in hyperfine/hyper_read_quanta'
  endif

  read(iounit,*) jmin, jmax, maxnapp
  read(iounit,*) dimen

  if (allocated(ind)) deallocate(ind)
  allocate( ind(3,dimen), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'hyperfine/hyper_read_quanta error: failed to allocate ind(3,dimen)', 'dimen =', dimen
    stop 'STOP, error in hyperfine/hyper_read_quanta'
  endif

  ind = 0
  do idimen=1, dimen
    read(iounit,*) ind(1:3,idimen)
  enddo

  read(iounit,'(a)') sbuf
  if (trim(sbuf)/='End hyperfine basis indexes') then
    write(out, '(/a,1x,a,1x,a,a,a)') &!
    'hyperfine/hyper_read_quanta error: file', trim(fname), 'has bogus footer = "', trim(sbuf), '", (expected "End hyperfine basis indexes")'
    stop 'STOP, error in hyperfine/hyper_read_quanta'
  endif

  close(iounit)
  call IOStop(fname)

  maxnstates = maxval(ind(3,:))

  if (verbose) then
    write(out, '(1x,a,1x,i3,1x,i3)') 'min and max spin-coupled J quanta:', jmin, jmax
    write(out, '(1x,a,1x,i6)') 'max number of rovibrational states (within J-block):', maxnstates
    write(out, '(1x,a,1x,i3)') 'max number of appearences of same J in different spin-rovibration coupling blocks:', maxnapp
  endif

  if (allocated(jmap_nelem)) deallocate(jmap_nelem)
  if (allocated(jmap_ind)) deallocate(jmap_ind)
  allocate( jmap_nelem(maxnstates,jmin:jmax), jmap_ind(maxnapp,maxnstates,jmin:jmax), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') &!
    'hyperfine/hyper_read_quanta error: failed to allocate  jmap_nelem(maxnstates,jmin:jmax), jmap_ind(maxnapp,maxnstates,jmin:jmax)', &!
    'maxnstates, jmin, jmax, maxnapp =', maxnstates, jmin, jmax, maxnapp
    stop 'STOP, error in hyperfine/hyper_read_quanta'
  endif

  jmap_nelem = 0
  jmap_ind = 0

  do idimen=1, dimen
    istate = ind(3,idimen)
    j = ind(2,idimen)
    jmap_nelem(istate,j) = jmap_nelem(istate,j) + 1
    nelem = jmap_nelem(istate,j)
    jmap_ind(nelem,istate,j) = idimen
  enddo

  if (verbose) write(out, '(100(1h.),1x,a)') 'hyperfine/hyper_read_quanta done'

end subroutine hyper_read_quanta



!###################################################################################################################################



subroutine hyper_read_vectors(hyper, fval, hmat_dimen, nnz, hmat_acsr, hmat_csr_ja, hmat_csr_ia, verbose_)

  type(hyperfine_type), intent(in) :: hyper
  real(rk), intent(in) :: fval
  integer(ik), intent(out) :: hmat_dimen, nnz
  real(rk), allocatable, intent(out) :: hmat_acsr(:)
  integer(ik), allocatable, intent(out) :: hmat_csr_ja(:), hmat_csr_ia(:)
  logical, optional, intent(in) :: verbose_

  integer(ik) :: iounit, info
  character(cl) :: sfval, fname
  character(23) :: sbuf23
  character(21) :: sbuf21
  logical :: verbose=.false.

  if (present(verbose_)) verbose=verbose_

  if (verbose) write(out, '(100(1h.),1x,a)') 'hyperfine/hyper_read_vectors start'

  write(sfval,'(f6.1)') fval
  fname = trim(adjustl(hyper%eigenfile_vectors))//'_'//trim(adjustl(sfval))//'.chk'
  call IOStart(fname,iounit)

  if (verbose) write(out, '(1x,a,1x,a)') 'read CSR-compressed eigenvectors from file:', trim(fname)

  open(iounit,file=fname,form='unformatted',action='read',status='old',position='rewind',iostat=info)
  if (info/=0) then
    write(out, '(/a,1x,a)') 'hyperfine/hyper_read_vectors error while opening file:', trim(fname)
    stop 'STOP, error in hyperfine/hyper_read_vectors'
  endif

  read(iounit) sbuf23

  if (trim(sbuf23)/='Start hyperfine vectors') then
    write(out, '(/a,1x,a,1x,a,a,a)') &!
    'hyperfine/hyper_read_vectors error: file', trim(fname), 'has bogus header = "', trim(sbuf23), '", (expected "Start hyperfine vectors")'
    stop 'STOP, error in hyperfine/hyper_read_vectors'
  endif

  read(iounit) hmat_dimen, nnz

  if (verbose) then
    write(out, '(1x,a,1x,i6)') 'dimension of basis:', hmat_dimen
    write(out, '(1x,a,1x,i8)') 'number of nonzero elements:', nnz
  endif

  if (allocated(hmat_acsr)) deallocate(hmat_acsr)
  if (allocated(hmat_csr_ja)) deallocate(hmat_csr_ja)
  if (allocated(hmat_csr_ia)) deallocate(hmat_csr_ia)
  allocate( hmat_acsr(nnz), hmat_csr_ja(nnz), hmat_csr_ia(hmat_dimen+1), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') &!
    'hyperfine/hyper_read_vectors error: failed to allocate hmat_acsr(nnz), hmat_csr_ja(nnz), hmat_csr_ia(hmat_dimen+1)', &!
    'nnz, hmat_dimen =', nnz, hmat_dimen
    stop 'STOP, error in hyperfine/hyper_read_vectors'
  endif

  read(iounit) hmat_acsr
  read(iounit) hmat_csr_ja
  read(iounit) hmat_csr_ia

  read(iounit) sbuf21

  if (trim(sbuf21)/='End hyperfine vectors') then
    write(out, '(/a,1x,a,1x,a,a,a)') &!
    'hyperfine/hyper_read_vectors error: file', trim(fname), 'has bogus footer = "', trim(sbuf21), '", (expected "End hyperfine vectors")'
    stop 'STOP, error in hyperfine/hyper_read_vectors'
  endif

  close(iounit)
  call IOStop(fname)

  if (verbose) write(out, '(100(1h.),1x,a)') 'hyperfine/hyper_read_vectors done'

end subroutine hyper_read_vectors



!###################################################################################################################################



subroutine hyper_read_descr(hyper, fval, levels, verbose_)

  type(hyperfine_type), intent(in) :: hyper
  real(rk), intent(in) :: fval
  type(hyper_states_type), allocatable, intent(out) :: levels
  logical, optional, intent(in) :: verbose_

  integer(ik) :: iounit, info, nmodes, nlevels, ilevel, ilevel_, ideg0_, iline
  real(rk) :: fval_, de_
  character(cl) :: fname, sfval, sbuf, sym0_
  logical :: verbose=.false.

  if (present(verbose_)) verbose=verbose_

  if (verbose) write(out, '(100(1h.),1x,a)') 'hyperfine/hyper_read_descr start'

  nmodes = hyper%filter%nmodes

  write(sfval,'(f6.1)') fval
  fname = trim(adjustl(hyper%eigenfile_descr))//'_'//trim(adjustl(sfval))//'.chk'
  call IOStart(fname,iounit)

  if (verbose) write(out, '(1x,a,1x,a)') 'read eigenvalues from file:', trim(fname)

  open(iounit,file=fname,form='formatted',action='read',status='old',position='rewind')
  if (info/=0) then
    write(out, '(/a,1x,a)') 'hyperfine/hyper_read_descr error while opening file:', trim(fname)
    stop 'STOP, error in hyperfine/hyper_read_descr'
  endif

  read(iounit,'(a)') sbuf
  iline = 1
  if (trim(sbuf)/='Start hyperfine energies') then
    write(out, '(/a,1x,a,1x,a,a,a)') &!
    'hyperfine/hyper_read_descr error: file', trim(fname), 'has bogus header = "', trim(sbuf), '", (expected "Start hyperfine energies")'
    stop 'STOP, error in hyperfine/hyper_read_descr'
  endif

  read(iounit,*) nlevels
  iline = iline + 1

  if (verbose) write(out, '(1x,a,1x,i6)') 'number of energy levels stored:', nlevels

  if (allocated(levels)) deallocate(levels)
  allocate( levels, stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'hyperfine/hyper_read_descr error: failed to allocate levels'
    stop 'STOP, error in hyperfine/hyper_read_descr'
  endif

  levels%nlevels = nlevels

  allocate( levels%sym(nlevels), levels%isym(nlevels), levels%energy(nlevels), levels%jval(nlevels), levels%enr_bas(nlevels), &!
            levels%k_quanta(nlevels), levels%v_quanta(nmodes,nlevels), levels%sym_rot(nlevels), levels%sym_vib(nlevels), &!
            levels%ideg(nlevels), levels%ndeg(nlevels), levels%ilevel(nlevels), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') &!
    'hyperfine/hyper_read_descr error: failed to allocate levels%sym(nlevels), levels%energy(nlevels), levels%jval(nlevels), ...', &!
    'nlevels, nmodes =', nlevels, nmodes
    stop 'STOP, error in hyperfine/hyper_read_descr'
  endif
  levels%sym = ''
  levels%isym = 0
  levels%energy = 0.0
  levels%k_quanta = 0
  levels%v_quanta = 0
  levels%sym_rot = ''
  levels%sym_vib = ''
  levels%jval = 0
  levels%ilevel = 0
  levels%ideg = 0
  levels%ndeg = 0
  levels%enr_bas = 0

  do ilevel=1, nlevels
    read(iounit,*,iostat=info) &!
    fval_, ilevel_, levels%sym(ilevel), levels%isym(ilevel), levels%energy(ilevel), levels%jval(ilevel), levels%ilevel(ilevel), &!
    levels%ideg(ilevel), levels%ndeg(ilevel), sym0_, levels%enr_bas(ilevel), de_, levels%k_quanta(ilevel), levels%v_quanta(:,ilevel), &!
    levels%sym_rot(ilevel), &!
    levels%sym_vib(ilevel)
    iline = iline + 1
    if (info/=0) then
      write(out, '(/a,1x,a,1x,a,1x,i6)') 'hyperfine/hyper_read_descr error while reading file =', trim(fname), ', line =', iline
      stop 'STOP, error in hyperfine/hyper_read_descr'
    endif
  enddo

  read(iounit,'(a)') sbuf
  iline = iline + 1
  if (trim(sbuf)/='End hyperfine energies') then
    write(out, '(/a,1x,a,1x,a,a,a)') &!
    'hyperfine/hyper_read_descr error: file', trim(fname), 'has bogus footer = "', trim(sbuf), '", (expected "End hyperfine energies")'
    stop 'STOP, error in hyperfine/hyper_read_descr'
  endif

  close(iounit)
  call IOStop(fname)

  if (verbose) write(out, '(100(1h.),1x,a)') 'hyperfine/hyper_read_descr done'

end subroutine hyper_read_descr



!###################################################################################################################################



subroutine hyper_int_fpair(hyper, f1, f2, linestr_tol, intens_tol)

  type(hyperfine_type), intent(in) :: hyper
  real(rk), intent(in) :: f1, f2, linestr_tol, intens_tol

  integer(ik) :: jmin, jmax, j1, j2, info, max_ang, dipole_nnz, vec1_dimen, vec2_dimen, vec_nnz, dimen, tmp_nnz, irequest, &!
                 ilevel, jlevel, ideg, jdeg, isym, isym1, iounit, nmodes, imode, ideg1, ideg2, ndeg1, ndeg2, nlevels1, nlevels2, &!
                 maxndeg, ilevel_, jlevel_
  integer(ik), allocatable :: dipole_csr_ja(:), dipole_csr_ia(:), vec_csr_ja(:), vec_csr_ia(:), tmp_csr_ia(:), tmp_csr_ja(:)
  integer(ik), allocatable :: level1_ind(:,:), level2_ind(:,:)
  real(rk) :: boltz_fac, boltz_beta, intens_cm_mol, linestr, intens, energy1, energy2, nu, nu_mhz, enr_degen_tol
  real(rk), allocatable :: dipole_acsr(:), vec_acsr(:), tmp_acsr(:), dipole(:,:)
  type(hyper_states_type), allocatable :: levels1, levels2
  logical :: f1_eq_f2
  character(cl) :: sf1, sf2, fname
  character(cl), allocatable :: squanta1(:), squanta2(:)

  jmin = hyper%filter%jmin
  jmax = hyper%filter%jmax

  nmodes = hyper%filter%nmodes

  ! init module "fwigxjpf" for computing 3j, 6j, and 9j symbols

  max_ang = jmax

  write(out, '(/1x,a,1x,i2,a,1x,i3)') 'initialize "fwigxjpf" for max', 6, '-j symbols and max angular momentum =', max_ang


#if (fwig>0)
  call fwig_table_init( int(2*max_ang,kind=4), int(6,kind=4) )
  call fwig_temp_init( int(2*max_ang,kind=4) )
#endif


  ! read energies and quantum numbers

  call hyper_read_descr(hyper, f1, levels1, verbose_=.true.)
  call hyper_read_descr(hyper, f2, levels2, verbose_=.true.)

  ! prepare text strings with assignments of initial and final states

  allocate(squanta1(levels1%nlevels),squanta2(levels2%nlevels), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'hyperfine/hyper_int_fpair error: failed to allocate squanta1(levels1%nlevels), squanta2(levels2%nlevels)', &!
    'levels1%nlevels, levels2%nlevels =', levels1%nlevels, levels2%nlevels
    stop 'STOP; error in hyperfine/hyper_int_fpair'
  endif

  do ilevel=1, levels1%nlevels
    write(squanta1(ilevel),'( 1(1h(),1x,a3,1x,i3,1x,i3,1x,1(1h)), 1x, 1(1h(),1x,a3,1x,<nmodes>(1x,i3),1x,1(1h))  )') &!
    trim(levels1%sym_rot(ilevel)), levels1%jval(ilevel), levels1%k_quanta(ilevel), &!
    trim(levels1%sym_vib(ilevel)), (levels1%v_quanta(imode,ilevel), imode=1, nmodes)
  enddo
  do ilevel=1, levels2%nlevels
    write(squanta2(ilevel),'( 1(1h(),1x,a3,1x,i3,1x,i3,1x,1(1h)), 1x, 1(1h(),1x,a3,1x,<nmodes>(1x,i3),1x,1(1h))  )') &!
    trim(levels2%sym_rot(ilevel)), levels2%jval(ilevel), levels2%k_quanta(ilevel), &!
    trim(levels2%sym_vib(ilevel)), (levels2%v_quanta(imode,ilevel), imode=1, nmodes)
  enddo



  ! read dipole moment in CSR format

  call hyper_dipole_me(hyper, f1, f2, dimen, dipole_nnz, dipole_acsr, dipole_csr_ja, dipole_csr_ia)

  ! read eigenvectors for final states in CSR format

  call hyper_read_vectors(hyper, f2, vec2_dimen, vec_nnz, vec_acsr, vec_csr_ja, vec_csr_ia, verbose_=.true.)


  ! half-transform dipole moment to hyperfine eigen-basis

  !   1. estimate the numbe rof non-zero elements

  allocate( tmp_csr_ia(dimen+1), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'hyperfine/hyper_int_fpair error: failed to allocate tmp_csr_ia(dimen+1)', 'dimen =', dimen
    stop 'STOP, error in hyperfine/hyper_int_fpair'
  endif

  write(out, '(/1x,a)') 'estimate number of non-zero elements in dipole*vec2 product matrix (nnz)'

  irequest = 1

  call mkl_dcsrmultcsr( 'N', irequest, 7, dimen, dimen, vec2_dimen, dipole_acsr, dipole_csr_ja, dipole_csr_ia, &!
                         vec_acsr, vec_csr_ja, vec_csr_ia, tmp_acsr, tmp_csr_ja, tmp_csr_ia, tmp_nnz, info )

  if (info/=0) then
    write(out, '(/a)') 'hyperfine/hyper_int_fpair error: failed to estimate nnz with mkl_dcsrmultcsr'
    stop 'STOP, error in hyperfine/hyper_int_fpair'
  endif

  tmp_nnz = tmp_csr_ia(dimen+1) - 1

  write(out, '(1x,a,1x,i8)') 'nnz =', tmp_nnz

  !   2. do the actual matrix multiplication

  allocate( tmp_acsr(tmp_nnz), tmp_csr_ja(tmp_nnz), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'hyperfine/hyper_int_fpair error: failed to allocate tmp_acsr(tmp_nnz), tmp_csr_ja(tmp_nnz)', 'tmp_nnz =', tmp_nnz
    stop 'STOP, error in hyperfine/hyper_int_fpair'
  endif

  write(out, '(/1x,a)') 'compute dipole*vec2 product matrix ...'

  irequest = 0

  call mkl_dcsrmultcsr( 'N', irequest, 7, dimen, dimen, vec2_dimen, dipole_acsr, dipole_csr_ja, dipole_csr_ia, &!
                         vec_acsr, vec_csr_ja, vec_csr_ia, tmp_acsr, tmp_csr_ja, tmp_csr_ia, tmp_nnz, info )

  if (info/=0) then
    write(out, '(/a)') 'hyperfine/hyper_int_fpair error: csr-csr matrix-matrix multiplication with mkl_dcsrmultcsr has failed'
    stop 'STOP, error in hyperfine/hyper_int_fpair'
  endif

  write(out, '(1x,a)') '... done'


  ! read eigenvectors for initial states in CSR format

  call hyper_read_vectors(hyper, f1, vec1_dimen, vec_nnz, vec_acsr, vec_csr_ja, vec_csr_ia, verbose_=.true.)

  ! full-transform dipole moment to hyperfine eigen-basis

  write(out, '(/1x,a)') 'compute vec1*dipole*vec2 product matrix ...'

  allocate( dipole(vec1_dimen,vec2_dimen), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'hyperfine/hyper_int_fpair error: failed to allocate dipole(vec1_dimen,vec2_dimen)', &!
    'vec1_dimen, vec2_dimen =', vec1_dimen, vec2_dimen
    stop 'STOP, error in hyperfine/hyper_int_fpair'
  endif

  call mkl_dcsrmultd( 'T', vec1_dimen, vec1_dimen, vec2_dimen, vec_acsr, vec_csr_ja, vec_csr_ia, &!
                      tmp_acsr, tmp_csr_ja, tmp_csr_ia, dipole, vec1_dimen )

  write(out, '(1x,a)') '... done'

  ! print intensities

  ! Bolzman beta constant
  boltz_beta = planck * vellgt / (boltz * hyper%temperature)
  ! convert dipole intensity into cm/mol if the linestrength is in Debye^2 and frequency is in cm^-1
  !intens_cm_mol = 8.0d-36 * real(pi,rk)**3 * avogno / (3.0_rk * planck * vellgt)
  intens_cm_mol = 250697.094895449972004014245288425_rk
  enr_degen_tol = 0.000001_rk

  maxndeg = max( maxval(levels1%ndeg(1:levels1%nlevels)), maxval(levels2%ndeg(1:levels2%nlevels)) )
  allocate(level1_ind(maxndeg,vec1_dimen),level2_ind(maxndeg,vec2_dimen), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'hyperfine/hyper_int_fpair error: failed to allocate level1_ind(maxndeg,vec1_dimen),level2_ind(maxndeg,vec2_dimen)', &!
    'maxndeg, vec1_dimen, vec2_dimen =', maxndeg, vec1_dimen, vec2_dimen
    stop 'STOP, error in hyperfine/hyper_int_fpair'
  endif

  call sym_degen_levels(vec1_dimen, levels1, enr_degen_tol, maxndeg, nlevels1, level1_ind)
  call sym_degen_levels(vec2_dimen, levels2, enr_degen_tol, maxndeg, nlevels2, level2_ind)

  write(sf1,'(f6.1)') f1
  write(sf2,'(f6.1)') f2
  fname = trim(adjustl(hyper%specfile))//'_f'//trim(adjustl(sf1))//'_f'//trim(adjustl(sf2))
  call IOStart(fname,iounit)

  open(iounit,file=fname,form='formatted',action='write',status='unknown',position='rewind')
  if (info/=0) then
    write(out, '(/a,1x,a)') 'hyperfine/hyper_int_fpair error while opening file:', trim(fname)
    stop 'STOP, error in hyperfine/hyper_int_fpair'
  endif

  f1_eq_f2 = .false.
  if (abs(f1-f2)<1.0d-03) f1_eq_f2 = .true.


  !do ilevel=1, vec1_dimen
  do ilevel_=1, nlevels1

    ilevel = level1_ind(1,ilevel_)

    isym1 =  levels1%isym(ilevel)
    energy1 = levels1%energy(ilevel)
    ideg1 = levels1%ideg(ilevel)
    ndeg1 = levels1%ndeg(ilevel)

    !do jlevel=1, vec2_dimen
    do jlevel_=1, nlevels2

      jlevel = level2_ind(1,jlevel_)
      ndeg2 = levels2%ndeg(jlevel)

      energy2 = levels2%energy(jlevel)
      nu = energy2 - energy1
      nu_mhz = nu * vellgt/1.0d+06 ! in MHz

      if (abs(nu)<1.0d-08) cycle
      if (f1_eq_f2.and.energy2<energy1) cycle

      linestr = 0
      do ideg=1, ndeg1
        ilevel = level1_ind(ideg,ilevel_)
        do jdeg=1, ndeg2
          jlevel = level2_ind(jdeg,jlevel_)
          linestr = linestr + dipole(ilevel,jlevel)**2
        enddo
      enddo

      linestr = hyper%gns(isym1) * linestr * (2.0_rk*f1+1.0_rk)*(2.0_rk*f2+1.0_rk) / real(ndeg1,rk)

      !linestr = hyper%gns(isym1)*dipole(ilevel,jlevel)**2
      boltz_fac = exp(-(energy1-hyper%ZPE) * boltz_beta) / hyper%part_func
      intens = linestr * boltz_fac * abs(nu) * (1.0_rk-exp(-abs(nu)*boltz_beta)) * intens_cm_mol

      if (intens>intens_tol.and.linestr>linestr_tol) then
        write(iounit, '( 1x,f4.1,1x,a4,1x,1(2h<-),1x,f4.1,1x,a4,3x,f12.6,1x,1(2h<-),1x,f12.6,3x,f12.6,2x,f16.4,3x,a,1x,1(2h<-),1x,a,3x,1(1x,es16.8),3x,1(1x,es16.8) )') &!
        f2, levels2%sym(jlevel), f1, levels1%sym(ilevel), &!
        energy2-hyper%ZPE, energy1-hyper%ZPE, abs(nu), abs(nu_mhz), trim(squanta2(jlevel)), trim(squanta1(ilevel)), linestr, intens
      endif

    enddo
  enddo

  close(iounit)


contains

  subroutine sym_degen_levels(nlevels, levels, enr_degen_tol, maxndeg, nlevels_deg, level_ind)

    integer(ik), intent(in) :: nlevels, maxndeg
    type(hyper_states_type), intent(in) :: levels
    real(rk), intent(in) :: enr_degen_tol
    integer(ik), intent(out) :: nlevels_deg
    integer(ik), intent(out) :: level_ind(maxndeg,nlevels)

    integer(ik) :: ilevel, jlevel, isym1, ideg1, ndeg1, ndeg, isym2, ideg2, ndeg2
    real(rk) :: energy1, energy2

    nlevels_deg = 0
    ilevel = 0

    do
      ilevel = ilevel + 1
      isym1 = levels%isym(ilevel)
      energy1 = levels%energy(ilevel)
      ideg1 = levels%ideg(ilevel)
      ndeg1 = levels%ndeg(ilevel)
      nlevels_deg = nlevels_deg + 1
      level_ind(ideg1,nlevels_deg) = ilevel
      ndeg = 1
      if (ndeg1>1) then
        if (ilevel==nlevels) then
          write(out, '(/a)') 'hyperfine/hyper_int_fpair/sym_degen_levels error: failed to determine degenerate components (i)'
          stop 'STOP, error in hyperfine/hyper_int_fpair/sym_degen_levels'
        endif
        do jlevel=ilevel+1, nlevels
          isym2 = levels%isym(jlevel)
          energy2 = levels%energy(jlevel)
          ideg2 = levels%ideg(jlevel)
          ndeg2 = levels%ndeg(jlevel)
          if (ndeg2==ndeg1.and.isym2==isym1.and.abs(energy2-energy1)<enr_degen_tol) then
            level_ind(ideg2,nlevels_deg) = jlevel
            ndeg = ndeg + 1
          else
            write(out, '(/a,1x,es16.8)') 'hyperfine/hyper_int_fpair/sym_degen_levels error: failed to determine degenerate components (ii)', abs(energy2-energy1)
            stop 'STOP, error in hyperfine/hyper_int_fpair/sym_degen_levels'
          endif
          if (ndeg==ndeg1) exit
        enddo
        ilevel = jlevel
      endif
      if (ilevel==nlevels) exit
    enddo

  end subroutine sym_degen_levels

end subroutine hyper_int_fpair



!###################################################################################################################################



subroutine hyper_dipole_me(hyper, f1, f2, dimen, dipole_nnz, dipole_acsr, dipole_csr_ja, dipole_csr_ia)

  type(hyperfine_type), intent(in) :: hyper
  real(rk), intent(in) :: f1, f2
  integer(ik), intent(out) :: dimen, dipole_nnz
  real(rk), allocatable, intent(out) :: dipole_acsr(:)
  integer(ik), allocatable, intent(out) :: dipole_csr_ja(:), dipole_csr_ia(:)

  integer(ik), parameter :: dipole_nnz_incr = 1000
  integer(ik) :: jmin1, jmax1, jmin2, jmax2, dimen1, dimen2, j1, j2, dipole_nnz_max, istate1, istate2, ikpair, ielem, jelem, &!
                 idimen, jdimen, info, job(8), jmin, jmax
  integer(ik), allocatable :: ind1(:,:), ind2(:,:), jmap_nelem1(:,:), jmap_nelem2(:,:), jmap_ind1(:,:,:), jmap_ind2(:,:,:), dipole_ind(:,:)
  real(rk) :: me, I1N_1, I1N_2, fac, fac1, fac2
  real(rk), allocatable :: dipole_coo(:)
  type(states_type), allocatable :: states(:)
  type(trans_type) :: trans

  write(out, '(100(1h.),1x,a)') 'hyperfine/hyper_dipole_me start'
  write(out, '(a,1x,f6.1,1x,f6.1)') 'Initialize matrix elements of the dipole moment for F-quanta pair =', f1, f2

  ! read rovibrational states file (RichMol database format)

  jmin = hyper%filter%jmin
  jmax = hyper%filter%jmax

  call read_states(hyper%states_file, hyper%filter, states)

  ! read hyperfine basis set indices

  write(out, '(/1x,a)') 'read spin-rovibrational basis set indices ...'

  call hyper_read_quanta(hyper, f1, jmin1, jmax1, dimen1, ind1, jmap_nelem1, jmap_ind1, verbose_=.false.)
  call hyper_read_quanta(hyper, f2, jmin2, jmax2, dimen2, ind2, jmap_nelem2, jmap_ind2, verbose_=.false.)

  write(out, '(1x,a)') '.. done'

  write(out, '(/1x,a,1x,f6.1/1x,a,1x,i6/1x,a,1x,i4,1x,i4)') &!
  'F1:', f1, 'dimension of the spin-rovibrational wavefunction:', dimen1, 'coupled min and max J quanta:', jmin1, jmax1

  write(out, '(/1x,a,1x,f6.1/1x,a,1x,i6/1x,a,1x,i4,1x,i4)') &!
  'F2:', f2, 'dimension of the spin-rovibrational wavefunction:', dimen2, 'coupled min and max J quanta:', jmin2, jmax2

  if (jmin1<jmin .or. jmin1>jmax .or. jmax1<jmin .or. jmax1>jmax) then
    write(out, '(/a,1x,i3,1x,i3,1x,a,1x,i3,1x,a,1x,i3,1x,a,1x,i3,1x,a)') &!
    'hyperfine/hyper_dipole_me error: jmin and jmax =', jmin1, jmax1, 'for F1 =', f1, 'run out of bounds j = [', jmin, ':', jmax, ']'
    stop 'STOP, error in hyperfine/hyper_dipole_me'
  endif

  if (jmin2<jmin .or. jmin2>jmax .or. jmax2<jmin .or. jmax2>jmax) then
    write(out, '(/a,1x,i3,1x,i3,1x,a,1x,i3,1x,a,1x,i3,1x,a,1x,i3,1x,a)') &!
    'hyperfine/hyper_dipole_me error: jmin and jmax =', jmin2, jmax2, 'for F2 =', f2, 'run out of bounds j = [', jmin, ':', jmax, ']'
    stop 'STOP, error in hyperfine/hyper_dipole_me'
  endif

  ! allocate array to keep matrix elements of dipole moment in sparse coordinate format

  dipole_nnz_max = dipole_nnz_incr
  allocate( dipole_ind(2,dipole_nnz_max), dipole_coo(dipole_nnz_max), stat=info )
  if (info/=0) then
    write(out, '(a/a,10(1x,i6))') &!
    'hyperfine/hyper_dipole_me error: failed to allocate dipole_ind(2,dipole_nnz_max), dipole_coo(dipole_nnz_max)', &!
    'dipole_nnz_max =', dipole_nnz_max
    stop 'STOP, error in hyperfine/hyper_dipole_me'
  endif
  dipole_ind = 0
  dipole_coo = 0

  write(out, '(/1x,a,1x,a)') 'start reading matrix elements of dipole moment from files', 'matelem_'//trim(hyper%dipole_tens_name)//'_j<j1>_j<j2>.rchm'

  dipole_nnz = 0

  do j1=jmin1, jmax1
    do j2=jmin2, jmax2

      call read_tran(j1, j2, jmin, jmax, trim(hyper%dipole_tens_name), states, hyper%me_thresh, trans, verbose_=.false.)

      if (trans%null) cycle ! trans%null=.true. means that all transition matrix elements for given j1/j2 pair are zero (by rigorous selection rules)
                            ! or simply the file with matrix elements does not exist

      ! double check if we are reading matrix elements of a first-rank tensor operator (dipole moment tensor)
      !  i.e. the number of Cartesian components is equal to 3

      if (trans%nelem/=3) then
        write(out, '(/a,1x,i3,1x,i3,1x,a,1x,a/a,1x,i3,1x,a)') &!
        'hyperfine/hyper_dipole_me error: file with j1/j2=', j1, j2, 'rovibrational matrix elements of tensor operator', trim(trans%tens_name), &!
        'returns number of tensor components =', trans%nelem, ', which is not equal to three as expected for the dipole moment tensor'
        stop 'STOP, error in hyperfine/hyper_dipole_me'
      endif

      fac = sqrt( (2.0_rk*j1+1.0_rk) * (2.0_rk*j2+1.0_rk) ) * (-1.0_rk)**(j1+j2)

      do ikpair=1, trans%nkpairs

        me = trans%kme(1,ikpair)

        istate1 = trans%kpair(1,ikpair)
        istate2 = trans%kpair(2,ikpair)

        do ielem=1, jmap_nelem1(istate1,j1)

          idimen = jmap_ind1(ielem,istate1,j1)
          I1N_1 = ind1(1,idimen)

          fac1 = fac * (-1.0_rk)**I1N_1

          do jelem=1, jmap_nelem2(istate2,j2)

            jdimen = jmap_ind2(jelem,istate2,j2)
            I1N_2 = ind2(1,jdimen)

            if (idimen>0.and.jdimen>0) then

              ! compute factor: (-1)^I_1N * sqrt((2J+1)*(2J'+1)) { J'  F'  I_1N }
              !                                                  { F   J     1  }

              fac2 = fac1 * sixj_symbol( real(j2,rk), f2, I1N_1, f1, real(j1,rk), 1.0_rk )
              if (abs(fac2)<hyper%me_thresh) cycle ! 6-j symbol nonvanishing rule

              dipole_nnz = dipole_nnz + 1
              if (dipole_nnz>dipole_nnz_max) call extend_dipole_me
              dipole_ind(1:2,dipole_nnz) = (/idimen,jdimen/)
              dipole_coo(dipole_nnz) = me * fac2

            endif
          enddo
        enddo

      enddo ! ikpair

    enddo ! j2
  enddo ! j1

  write(out, '(1x,a)') '.. done'

  write(out, '(/1x,a,1x,i8,1x,a,1x,es16.8,1x,a)') 'number of non-zero matrix elements:', dipole_nnz, '(>', hyper%me_thresh, ')'

  ! convert sparse matrix from coordinate to csr sparse format

  write(out, '(/1x,a)') 'convert dipole matrix from coordinate to CSR format ...'

  job(1:6) = (/1, 1, 1, 0, dipole_nnz, 0/)

  dimen = max(dimen1,dimen2)

  if (allocated(dipole_acsr)) deallocate(dipole_acsr)
  if (allocated(dipole_csr_ja)) deallocate(dipole_csr_ja)
  if (allocated(dipole_csr_ia)) deallocate(dipole_csr_ia)
  allocate( dipole_acsr(dipole_nnz), dipole_csr_ja(dipole_nnz), dipole_csr_ia(dimen+1), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') &!
    'hyperfine/hyper_dipole_me error: failed to allocate dipole_acsr(dipole_nnz), dipole_csr_ja(dipole_nnz), dipole_csr_ia(dimen+1)', &!
    'dipole_nnz, dimen =', dipole_nnz, dimen
    stop 'STOP, error in hyperfine/hyper_dipole_me'
  endif

  call mkl_dcsrcoo( job, dimen, dipole_acsr(1:dipole_nnz), dipole_csr_ja(1:dipole_nnz), dipole_csr_ia(1:dimen+1), &!
                    dipole_nnz, dipole_coo(1:dipole_nnz), dipole_ind(1,1:dipole_nnz), dipole_ind(2,1:dipole_nnz), info)

  if (info/=0) then
    write(out, '(/a)') 'hyperfine/hyper_dipole_me error: coo-to-csr format conversion with mkl_dcsrcoo has failed'
    stop 'STOP, error in hyperfine/hyper_dipole_me'
  endif

  write(out, '(1x,a)') '... done'

  deallocate(ind1, ind2, jmap_nelem1, jmap_nelem2, jmap_ind1, jmap_ind2, dipole_ind, dipole_coo)

  write(out, '(100(1h.),1x,a)') 'hyperfine/hyper_dipole_me done'

  contains

  subroutine extend_dipole_me
    integer(ik) :: dipole_nnz_max_
    integer(ik), allocatable :: itmp(:,:)
    real(rk), allocatable :: rtmp(:)
    call move_alloc(dipole_ind,itmp)
    call move_alloc(dipole_coo,rtmp)
    dipole_nnz_max_ = dipole_nnz_max
    dipole_nnz_max = dipole_nnz_max +  dipole_nnz_incr
    allocate( dipole_ind(2,dipole_nnz_max), dipole_coo(dipole_nnz_max), stat=info)
    if (info/=0) then
      write(out, '(/a/a,10(1x,i6))') &!
      'hyperfine/hyper_dipole_me error: failed to allocate dipole_ind(2,dipole_nnz_max), dipole_coo(dipole_nnz_max)', &!
      'dipole_nnz_max =', dipole_nnz_max
      stop 'STOP, error in hyperfine/hyper_dipole_me'
    endif
    dipole_ind = 0
    dipole_coo = 0
    dipole_ind(1:2,1:dipole_nnz_max_) = itmp
    dipole_coo(1:dipole_nnz_max_) = rtmp
    deallocate(itmp,rtmp)
  end subroutine extend_dipole_me

end subroutine hyper_dipole_me



!###################################################################################################################################



recursive subroutine coupl_ang_quanta(iang, nang, ang, ang_coupl, dimen, ang_coupl_list)

  integer(ik), intent(in) :: iang, nang
  integer(ik), intent(inout) :: dimen
  real(rk), intent(in) :: ang(nang)
  real(rk), intent(inout) :: ang_coupl(nang)
  real(rk), intent(inout), optional :: ang_coupl_list(:,:)

  real(rk) :: tol, ang12, ang1, ang2
  logical :: last

  tol = 1.0d-04

  if (iang==nang-1) then
    last = .true.
  else
    last = .false.
  endif

  if (iang==1) then
    dimen = 0
    ang_coupl = 0
    ang1 = ang(iang)
  else
    ang1 = ang_coupl(iang-1)
  endif

  ang2 = ang(iang+1)

  ang12 = abs(ang1-ang2)
  do while (ang12-(ang1+ang2)<=tol)
    ang_coupl(iang) = ang12
    if (last) then
      dimen = dimen + 1
      if (present(ang_coupl_list)) then
        if (dimen>size(ang_coupl_list,dim=2)) then
          write(out, '(/a,1x,i6)') 'hyperfine/coupl_ang_quanta error: "dimen" exceeds size of array "ang_coupl_list" =', size(ang_coupl_list,dim=2)
          stop 'STOP, error in hyperfine/coupl_ang_quanta'
        endif
        ang_coupl_list(1,dimen) = ang(1)
        ang_coupl_list(2:nang,dimen) = ang_coupl(1:nang-1)
      endif
    else
      call coupl_ang_quanta(iang+1, nang, ang, ang_coupl, dimen, ang_coupl_list)
    endif
    ang12 = ang12 + 1.0_rk
  enddo

end subroutine coupl_ang_quanta



!###################################################################################################################################



recursive function nucspin_reduced_me(hyper, n, iang, ang_coupl_bra, ang_coupl_ket, rank, oper) result(f)

  type(hyperfine_type), intent(in) :: hyper
  integer(ik), intent(in) :: n, iang, rank
  real(rk), intent(in) :: ang_coupl_bra(hyper%nspin), ang_coupl_ket(hyper%nspin)
  character(*), intent(in) :: oper
  real(rk) :: f

  integer(ik) :: ifac, nang
  real(rk) :: coef1, coef2, res, tol, fac, ang_bra(hyper%nspin), ang_ket(hyper%nspin), equad

  tol = 1.0d-06 ! tolerance factor to check if two half integer quantum numbers are equal

  nang = hyper%nspin
  ang_bra(1:nang) = hyper%spin(1:nang)
  ang_ket(1:nang) = hyper%spin(1:nang)
  equad = hyper%quad_const(iang)

  if (n==iang) then

    if (iang==1) then

      coef1 = 1.0_rk

    else

      fac = ang_coupl_ket(n-1) + ang_ket(n) + ang_coupl_bra(n) + real(rank,rk)
      ifac = nint(fac)
      if (abs(real(ifac,rk)-fac)>tol) then
        write(out, '(/a,1x,f20.12)') 'hyperfine/nucspin_reduced_me error: non-integer power factor in (-1)**f, f =', fac
        stop 'STOP, error in hyperfine/nucspin_reduced_me'
      endif
      fac = (-1.0_rk)**ifac

      coef1 = fac * sqrt( (2.0_rk*ang_coupl_bra(n)+1.0_rk) * (2.0_rk*ang_coupl_ket(n)+1.0_rk) ) &!
            * sixj_symbol(       ang_bra(n), ang_coupl_bra(n), ang_coupl_ket(n-1), &!
                           ang_coupl_ket(n),       ang_ket(n),      real(rank,rk) )
    endif

    select case(trim(oper))
    case('QUAD','quad')
      coef2 = reduced_me_quad(equad, ang_bra(iang), ang_ket(iang))
    case('SPIN','spin')
      coef2 = reduced_me_nsang(iang, ang_bra(iang), ang_ket(iang))
    case default
      write(out, '(/a,1x,a)') 'hyperfine/nucspin_reduced_me error: illegal nuclear spin operator type =', trim(oper)
      stop 'STOP, error in hyperfine/nucspin_reduced_me'
    end select

  elseif (n>iang) then

    fac = ang_coupl_bra(n-1) + ang_ket(n) + ang_coupl_ket(n) + real(rank,rk)
    ifac = nint(fac)
    if (abs(real(ifac,rk)-fac)>tol) then
      write(out, '(/a,1x,f20.12)') 'hyperfine/nucspin_reduced_me error: non-integer power factor in (-1)**f, f =', fac
      stop 'STOP, error in hyperfine/nucspin_reduced_me'
    endif
    fac = (-1.0_rk)**ifac

    coef1 = fac * sqrt( (2.0_rk*ang_coupl_bra(n)+1.0_rk) * (2.0_rk*ang_coupl_ket(n)+1.0_rk) ) &!
          * sixj_symbol( ang_coupl_bra(n-1),   ang_coupl_bra(n),    ang_ket(n), &!
                           ang_coupl_ket(n), ang_coupl_ket(n-1), real(rank,rk) )

    coef2 = nucspin_reduced_me(hyper, n-1, iang, ang_coupl_bra, ang_coupl_ket, rank, oper)

  else

    write(out, '(/a,1x,i4,1x,a,1x,i4)') 'hyperfine/nucspin_reduced_me error: n =', n, '< iang =', iang
    stop 'STOP, error in hyperfine/nucspin_reduced_me'

  endif

  f = coef1 * coef2

end function nucspin_reduced_me



!###################################################################################################################################



! Computes reduced matrix element of the nuclear quadruple moment < I_bra || Q^(2) || I_ket >

function reduced_me_quad(equad, I_bra, I_ket) result(f)

  real(rk), intent(in) :: equad, I_bra, I_ket
  real(rk) :: f

  real(rk) :: small, small_, threej

  small = 1.0d-06 ! threshold for two half integer (spin) quantum numbers to be equal
  small_ = epsilon(1.0d0)*100.0d0 ! threshold for zero

  if (abs(I_bra-I_ket)<=small) then

    threej = threej_symbol(I_bra,-I_bra, 2.0_rk,0.0_rk, I_bra,I_bra)

    if (abs(threej)<=small_.and.abs(equad)<=small_) then
      f = 0
    elseif (abs(threej)<=small_.and.abs(equad)>small_) then
      write(out, '(/a,1x,f20.12,1x,a,1x,f20.12)') 'hyperfine/reduced_me_quad error: illegal division eQ/3j =', equad, '/', threej
      stop 'STOP, error in hyperfine/reduced_me_quad'
    else
      f = 0.5_rk * equad / threej * VQ_mbau_cm
    endif

  else

    f = 0

  endif

end function reduced_me_quad



!###################################################################################################################################



! Computes reduced matrix element of nuclear spin angular momentum < I_bra || I^(1) || I_ket >

function reduced_me_nsang(iang, I_bra, I_ket) result(f)

  integer(ik), intent(in) :: iang
  real(rk), intent(in) :: I_bra, I_ket
  real(rk) :: f

  real(rk) :: small, small_, threej

  small = 1.0d-06 ! threshold for two half integer (spin) quantum numbers to be equal
  small_ = epsilon(1.0d0)*100.0d0 ! threshold for zero

  if (abs(I_bra-I_ket)<=small) then

    threej = threej_symbol(I_bra,-I_bra, 1.0_rk,0.0_rk, I_bra,I_bra)

    if (abs(threej)<=small_.and.abs(I_bra)<=small_) then
      f = 0
    elseif (abs(threej)<=small_.and.abs(I_bra)>small_) then
      write(out, '(/a,1x,f20.12,1x,a,1x,f20.12)') 'hyperfine/reduced_me_nsang error: illegal division I/3j =', I_bra, '/', threej
      stop 'STOP, error in hyperfine/reduced_me_nsang'
    else
      f = I_bra / threej
    endif

  else
    f = 0
  endif

end function reduced_me_nsang



!###################################################################################################################################



subroutine hyper_readinp(fname, hyper)

  character(*), intent(in) :: fname
  type(hyperfine_type), intent(out) :: hyper

  integer(ik), parameter :: maxnw=100
  integer(ik) info, nw, iw, iw_, i1, i2
  character(cl) w, w2, wrd(maxnw), s, ww
  character(500) sbuf

  integer(ik), parameter :: maxnspins=100, maxnsym=100
  integer(ik) :: iline, inp_unit, i, ind(0:maxnspins), nspin, nquad, ntens_efg, nsym
  real(rk) :: spin(maxnspins), quad(maxnspins), gns(maxnsym), ftot
  character(cl) :: tens_efg(maxnspins)

  write(out, '(///10x,a)') '---------------------------------'
  write(out, '(10x,a)')    'Read input for hyperfine coupling'
  write(out, '(10x,a)')    '---------------------------------'

  inp_unit = inp
  !call IOstart(fname,inp_unit)
  !open(unit=inp_unit,file=fname,form='formatted',status='old',iostat=info)
  !if (info/=0) then
  !  write(out, '(/a,1x,a)') 'hyperfine/hyper_readinp error: failed to open file', trim(fname)
  !  stop 'STOP, error in hyperfine/hyper_readinp'
  !endif

  nspin = 0
  nquad = 0
  ntens_efg = 0
  nsym = 0

  write(out, '(/a)') '------- echo of the input -------'
  rewind(inp_unit)
  do
    read(inp_unit,'(a)',iostat=info) sbuf
    if (info/=0) exit
    call tokenize(sbuf, nw, wrd)

    if (trim(adjustl(to_upper(wrd(1))))=='HYPERFINE') then

      write(out, '(a)') trim(sbuf)

      iline = 0
      do
        read(inp_unit,'(a)',iostat=info) sbuf   ;   write(out, '(a)') trim(sbuf)
        if (info/=0) then
          write(out, '(/a,1x,i3,1x,a)') 'hyperfine/hyper_readinp error while reading line', iline+1, 'in "HYPERFINE...END" input block'
          stop 'STOP, error in hyperfine/hyper_readinp'
        endif
        iline = iline + 1

        call tokenize(sbuf, nw, wrd)

        if (nw==0) then
          cycle
        elseif (trim(to_upper(wrd(1)))=='END') then
          exit
        endif

        ! skip lines beginning with '#'

        w = trim(adjustl(to_upper(wrd(1))))
        if (w(1:1)=='#') cycle

        iw = 0
        do iw_=1, nw
          iw = iw + 1  ;  if (iw>nw) exit
          w = to_upper(wrd(iw))
          w2 = wrd(iw)
          !
          if (w(1:5)=='TASK=') then
            hyper%task = trim(w(6:))
          elseif (w(1:2)=='F=') then
            read(w(3:),*,iostat=info) hyper%ftot
          elseif (w(1:5)=='JMIN=') then
            read(w(6:),*,iostat=info) hyper%filter%jmin
          elseif (w(1:5)=='JMAX=') then
            read(w(6:),*,iostat=info) hyper%filter%jmax
          elseif (w(1:5)=='EMIN=') then
            read(w(6:),*,iostat=info) hyper%filter%emin
          elseif (w(1:5)=='EMAX=') then
            read(w(6:),*,iostat=info) hyper%filter%emax
          elseif (w(1:7)=='NMODES=') then
            read(w(8:),*,iostat=info) hyper%filter%nmodes
          elseif (w(1:10)=='ME_THRESH=') then
            read(w(11:),*,iostat=info) hyper%me_thresh
          elseif (w(1:7)=='STATES=') then
            hyper%states_file = w2(8:)
          elseif (w(1:6)=='SPINS=') then
            ww = trim(w2(7:))
            nspin = 0
            i = 0
            ind(0) = 0
            do while(.true.)
              i = i + 1
            	ind(i) = index(ww(ind(i-1)+1:),',')
              if (ind(i)==0) then
            	  read(ww(ind(i-1)+1:),*,iostat=info) spin(i)
              else
 	            	ind(i) = ind(i) + ind(i-1)
              	read(ww(ind(i-1)+1:ind(i)-1),*,iostat=info) spin(i)
              endif
              if (info/=0) then
                write(out, '(/a,1x,i3,1x,a)') 'hyperfine/hyper_readinp error while reading line', iline, 'in "HYPERFINE...END" input block'
                stop 'STOP, error in hyperfine/hyper_readinp'
              endif
          	  if (ind(i)==0) exit
            enddo
            nspin = i
          elseif (w(1:5)=='QUAD=') then
            ww = trim(w2(6:))
            nquad = 0
            i = 0
            ind(0) = 0
            do while(.true.)
              i = i + 1
            	ind(i) = index(ww(ind(i-1)+1:),',')
              if (ind(i)==0) then
            	  read(ww(ind(i-1)+1:),*,iostat=info) quad(i)
              else
 	            	ind(i) = ind(i) + ind(i-1)
              	read(ww(ind(i-1)+1:ind(i)-1),*,iostat=info) quad(i)
              endif
              if (info/=0) then
                write(out, '(/a,1x,i3,1x,a)') 'hyperfine/hyper_readinp error while reading line', iline, 'in "HYPERFINE...END" input block'
                stop 'STOP, error in hyperfine/hyper_readinp'
              endif
          	  if (ind(i)==0) exit
            enddo
            nquad = i
          elseif (w(1:4)=='EFG=') then
            ww = trim(w2(5:))
            ntens_efg = 0
            i = 0
            ind(0) = 0
            do while(.true.)
              i = i + 1
            	ind(i) = index(ww(ind(i-1)+1:),',')
              if (ind(i)==0) then
                tens_efg(i) = trim(adjustl(ww(ind(i-1)+1:)))
              else
 	            	ind(i) = ind(i) + ind(i-1)
 	            	tens_efg(i) = trim(adjustl(ww(ind(i-1)+1:ind(i)-1)))
              endif
          	  if (ind(i)==0) exit
            enddo
            ntens_efg = i
          elseif (w(1:4)=='GNS=') then
            ww = trim(w2(5:))
            nsym = 0
            i = 0
            ind(0) = 0
            do while(.true.)
              i = i + 1
            	ind(i) = index(ww(ind(i-1)+1:),',')
              if (ind(i)==0) then
            	  read(ww(ind(i-1)+1:),*,iostat=info) gns(i)
              else
 	            	ind(i) = ind(i) + ind(i-1)
              	read(ww(ind(i-1)+1:ind(i)-1),*,iostat=info) gns(i)
              endif
              if (info/=0) then
                write(out, '(/a,1x,i3,1x,a)') 'hyperfine/hyper_readinp error while reading line', iline, 'in "HYPERFINE...END" input block'
                stop 'STOP, error in hyperfine/hyper_readinp'
              endif
          	  if (ind(i)==0) exit
            enddo
            nsym = i
          elseif (w(1:5)=='TEMP=') then
            read(w(6:),*,iostat=info) hyper%temperature
          elseif (w(1:4)=='ZPE=') then
            read(w(5:),*,iostat=info) hyper%zpe
          elseif (w(1:2)=='Q=') then
            read(w(3:),*,iostat=info) hyper%part_func
          else
            write(out, '(/a,a,a)') 'hyperfine/hyper_readinp error: illegal keyword = "', trim(w2), '" in "HYPERFINE...END" input block'
            stop 'STOP, error in hyperfine/hyper_readinp'
          endif
          !
          if (info/=0) then
            write(out, '(/a,1x,i3,1x,a)') 'hyperfine/hyper_readinp error while reading line', iline, 'in "HYPERFINE...END" input block'
            stop 'STOP, error in hyperfine/hyper_readinp'
          endif
          !
        enddo

      enddo
    endif

  enddo
  write(out, '(a)') '---------------------------------'

  !close(inp_unit)
  !call IOstop(fname)

  if (trim(hyper%task)=='ENERGIES') then
    write(out, '(/1x,a)') 'hyperfine-resolved energies will be computed'
  elseif (trim(hyper%task)=='SPECTRUM') then
    write(out, '(/1x,a)') 'hyperfine-resolved spectrum will be computed'
  else
    write(out, '(/a,1x,a)') 'hyperfine/hyper_readinp error: illegal value for TASK keyword =', trim(hyper%task)
    stop 'STOP, error in hyperfine/hyper_readinp'
  endif

  write(out, '(/1x,a,1x,f4.1)') 'quantum number of total F=J+I angular momentum:', hyper%ftot

  write(out, '(1x,a,1x,i3,1x,a,100(1x,f4.1))') 'nuclear spins (', nspin, '):', spin(1:nspin)

  if (nspin==0) then
    write(out, '(/a)') 'hyperfine/hyper_readinp error: number of input nuclear spins = 0 (check "SPINS=<spin1>,<spin2>,...,<spinN>")'
    stop 'STOP, error in hyperfine/hyper_readinp'
  endif

  write(out, '(/1x,a,1x,i3,1x,a)') 'nuclear quadrupole constants (', nquad, '):'

  if (nquad/=nspin) then
    write(out, '(/a,1x,i3,1x,a,1x,i3/a)') &!
    'hyperfine/hyper_readinp error: number of input nuclear spins =', nspin, '.ne. number of input nuclear quadrupole constants =', nquad, &!
    '(check "SPINS=<spin1>,<spin2>,...,<spinN>" and "QUAD=<const1>,<const2>,...,<constN>")'
    stop 'STOP, error in hyperfine/hyper_readinp'
  endif

  if (nquad/=ntens_efg) then
    write(out, '(/a,1x,i3,1x,a,1x,i3/a)') &!
    'hyperfine/hyper_readinp error: number of input quadrupolar nuclei =', nquad, '.ne. number of input electric field gradient tensor labels =', ntens_efg, &!
    '(check "QUAD=<const1>,<const2>,...,<constN>" and "EFG=<label1>,<label2>,...,<labelN>")'
    stop 'STOP, error in hyperfine/hyper_readinp'
  endif

  write(out, '(2(/9x,a,12x,a,3x,a))') 'spin', 'eQ [mb]', 'EFG tensor', '----', '-------', '----------'

  do i=1, nquad
    write(out, '(3x,i3,3x,f4.1,3x,es16.8,3x,a)') i, spin(i), quad(i), trim(tens_efg(i))
  enddo

  if (nsym==0) then
    write(out, '(/a)') 'hyperfine/hyper_readinp error: number of symmetries = 0 (check "GNS=<gns1>,<gns2>,...,<gnsN>")'
    stop 'STOP, error in hyperfine/hyper_readinp'
  endif

  write(out, '(/1x,a,1x,i3,1x,a,100(1x,f6.1))') 'Gns factors (', nsym, '):', gns(1:nsym)


  hyper%nspin = nspin
  hyper%nsym = nsym

  if (allocated(hyper%spin)) deallocate(hyper%spin)
  if (allocated(hyper%quad_const)) deallocate(hyper%quad_const)
  if (allocated(hyper%efg_tens_name)) deallocate(hyper%efg_tens_name)
  if (allocated(hyper%gns)) deallocate(hyper%gns)
  allocate( hyper%spin(nspin), hyper%quad_const(nspin), hyper%efg_tens_name(nspin), hyper%gns(nsym), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') &!
    'hyperfine/hyper_readinp error: failed to allocate hyper%spin(nspin), hyper%quad_const(nspin), hyper%efg_tens_name(nspin), hyper%gns(nsym)', &!
    'nspin, nsym =', nspin, nsym
    stop 'STOP, error in hyperfine/hyper_readinp'
  endif
  hyper%spin(1:nspin) = spin(1:nspin)
  hyper%quad_const(1:nspin) = quad(1:nspin)
  hyper%efg_tens_name(1:nspin) = tens_efg(1:nspin)
  hyper%gns(1:nsym) = gns(1:nsym)

end subroutine hyper_readinp



!###################################################################################################################################



function threej_symbol(j1,m1, j2,m2, j3,m3) result(f)

  real(rk), intent(in) :: j1,j2,j3,m1,m2,m3
  real(rk) :: f

  integer(4) :: two_j1, two_j2, two_j3, two_m1, two_m2, two_m3

  two_j1 = nint(2.0_rk*j1)
  two_j2 = nint(2.0_rk*j2)
  two_j3 = nint(2.0_rk*j3)

  two_m1 = nint(2.0_rk*m1)
  two_m2 = nint(2.0_rk*m2)
  two_m3 = nint(2.0_rk*m3)


#if (fwig>0)
  f = fwig3jj(two_j1, two_j2, two_j3, two_m1, two_m2, two_m3)
#else
  f = 0 
#endif


end function threej_symbol



!###################################################################################################################################



function sixj_symbol(j1,j2,j3, j4,j5,j6) result(f)

  real(rk), intent(in) :: j1,j2,j3,j4,j5,j6
  real(rk) :: f

  integer(4) :: two_j1, two_j2, two_j3, two_j4, two_j5, two_j6

  two_j1 = nint(2.0_rk*j1)
  two_j2 = nint(2.0_rk*j2)
  two_j3 = nint(2.0_rk*j3)

  two_j4 = nint(2.0_rk*j4)
  two_j5 = nint(2.0_rk*j5)
  two_j6 = nint(2.0_rk*j6)


#if (fwig>0)
  f = fwig6jj(two_j1, two_j2, two_j3, two_j4, two_j5, two_j6)
#else
  f = 0 
#endif


end function sixj_symbol


end module hyperfine
