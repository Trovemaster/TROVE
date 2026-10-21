module FITPES_ch3br
  !
  ! Standalone least-squares fitting program for the 9D symmetry-adapted
  ! Potential Energy Function (PEF) of ZXY3 molecules (CH3Br / CH3Cl)
  ! in C3v(M) symmetry (TROVE POT_TYPE = poten_zxy3_sym).
  !
  ! Adapted from FITPES_ch3cl (C2H6.f90) and TROVE pot_zxy3.f90.
  !
  implicit none
  public :: ik, rk, do_fitting
  private

  integer, parameter :: ik        = selected_int_kind(8)
  integer, parameter :: rk        = selected_real_kind(12, 25)
  integer, parameter :: f_inp     = 5, f_out = 6, f_res = 9, f_trove = 11
  integer, parameter :: nparam_eq = 5
  integer(ik), parameter :: verbose = 4

contains

  subroutine do_fitting
    integer, parameter          :: enermax = 150000, nparam = 1000
    integer, parameter          :: findif = 3
    double precision, parameter :: fact_corr = 0.4d0, too_big = 1.0d9
    double precision, parameter :: a_wats = 0.1d0
    double precision, parameter :: factordeltax = 0.001d0
    integer, parameter          :: ireplace_by = 0
    double precision, parameter :: bohr = 0.52917715d0, hartre = 2.194746354d5
    integer, parameter          :: ishift = 10
    character(len=40), parameter :: fit_type = 'linur'

    double precision  :: maxdeviat
    integer           :: i, npts, nused, itmax, numpar, ncol, iter, nrow, l, irow, icolumn
    character(len=80) :: title(4), longlabel
    character(len=1)  :: mark
    character(len=10) :: label
    double precision  :: wtsum, ssq
    double precision  :: stadev, epsil, stadev_old, v, pi
    double precision  :: stability, corr, ve
    integer           :: ierror, i1, i2, last_i, robustfit
    double precision  :: stab_best, conf_int
    integer           :: maxivar, minivar, tivar
    double precision  :: last_stadev, factor_back, factor_out
    double precision  :: rms, ZPE, local(18), stadev_best
    integer           :: last_ivar, ndigits, alloc
    integer           :: ipar, i0, l0, parmax, fivar_t, nmax, taumax
    logical           :: yes, still_run, ifopen, find_best_params
    integer           :: ivartmp, imin, imax, ilarge_sterr, ilow_sterr, Nout, parmax_t
    double precision  :: conv_length, conv_energy, paramtmp, relat_sterr, relat_sterr_min, fparam_t
    character(len=68) :: fmt1
    character(len=2)  :: fmt0
    character(len=120):: char_tmp
    integer           :: ipower_(18), iparam

    double precision, allocatable :: rjacob(:,:), last_param(:), param(:)
    double precision, allocatable :: param_tmp(:)
    double precision, allocatable :: crr(:), coord(:,:), coord_tmp(:,:)
    double precision, allocatable :: energy(:), energy_tmp(:)
    double precision, allocatable :: eps(:), wt(:), wt_in(:), parold(:), sigma(:), wt_tmp(:), parinitial(:)
    double precision, allocatable :: al(:,:), bl(:), dx(:)
    double precision, allocatable :: ai(:,:), sterr(:), dV(:,:)
    double precision, allocatable :: sens(:)
    character(len=10), allocatable:: parnam(:), parnam_tmp(:)
    integer, allocatable          :: ivar(:), icorr(:), ivar0(:), ivar_tmp(:), ipower_tmp(:,:), ipower(:,:)

    double precision, allocatable :: Tsing(:)
    integer                       :: rank0, info
    integer                       :: lwork
    double precision, allocatable :: wspace(:)

    last_stadev = 1.0d6
    last_i      = 0
    maxivar     = 8
    maxdeviat   = 2000.0d0

    pi = 4.0d0 * datan2(1.0d0, 1.0d0)

    call skiplines(f_inp, 1)

    do i = 1, 4
      read (f_inp, "(a80)") title(i)
    enddo

    write (f_out, "(3(132('*')/' '),16('*'),100x,16('*'))")
    write (f_out, "(4(' ',16('*'),10x,a80,10x,16('*')/))") (title(i), i = 1, 4)
    write (f_out, "(' ',16('*'),100x,16('*')/' ',3(132('*')/' '))")

    call skiplines(f_inp, 4)

    read (f_inp, *) itmax
    read (f_inp, *) stab_best, stadev_best
    read (f_inp, *) conv_length
    read (f_inp, *) conv_energy
    read (f_inp, *) factor_out
    read (f_inp, *) factor_back
    read (f_inp, *) robustfit

    write (f_out, "(i10,' <= number of iterations in the fit')")                  itmax
    write (f_out, "(f10.6,' <= The fit is over if the stability is reached ')")   stab_best
    write (f_out, "(f10.6,' <= The fit is over if the st.dev. best is reached ')")stadev_best
    write (f_out, "(f10.6,' <= sterr parameter out ')")                           factor_out
    write (f_out, "(f10.6,' <= stdev parameter out ')")                           factor_back
    write (f_out, "(i10,  ' <= watson robust on (1) or off (0) ')")               robustfit

    if (robustfit < 0) then
      find_best_params = .false.
    else
      find_best_params = .true.
    endif

    parmax_t = nparam
    lwork    = parmax_t * 16

    allocate(param_tmp(parmax_t), parnam_tmp(parmax_t), ivar_tmp(parmax_t), &
             ipower_tmp(1:9, parmax_t), stat=alloc)
    if (alloc /= 0) stop 'param_tmp - out of memory'

    ipower_tmp = 0
    parnam_tmp = 'xxxxxx'
    param_tmp  = 0.0d0
    ivar_tmp   = 0

    call skiplines(f_inp, 3)
    do i = 1, nparam_eq
      read (f_inp, "(a80)") longlabel
      read (longlabel, *) label, ivartmp, paramtmp
      parnam_tmp(i) = label
      ivar_tmp(i)   = ivartmp
      param_tmp(i)  = paramtmp
    enddo

    write (f_out, "(' parameters of structure ')")
    do i = 1, nparam_eq
      write (f_out, "(a10,'=',f16.8)") parnam_tmp(i), param_tmp(i)
    enddo

    yes    = .true.
    i      = nparam_eq
    nmax   = 0
    taumax = 0
    read_loop: do while (yes)
      read (f_inp, "(a120)") char_tmp
      if (char_tmp(1:3) /= 'end' .and. char_tmp(1:3) /= 'END') then
        i = i + 1
        if (i > parmax_t) stop 'Exceeded nparam! Increase nparam in FITPES_ch3br.'
        read (char_tmp, *) parnam_tmp(i), ipower_(1:9), fivar_t, fparam_t
        write (f_out, "(a2,9i3,' =',f18.8)") trim(parnam_tmp(i)), ipower_(1:9), fparam_t
        param_tmp(i)       = fparam_t
        ivar_tmp(i)        = fivar_t
        ipower_tmp(1:9, i) = ipower_(1:9)
      else
        yes = .false.
      endif
    enddo read_loop

    parmax = i

    allocate(param(parmax), parnam(parmax), ivar(parmax), icorr(parmax), &
             ivar0(parmax), ipower(1:18, parmax), stat=alloc)
    if (alloc /= 0) stop 'param - out of memory'

    ipower               = 0
    param(1:parmax)      = param_tmp(1:parmax)
    ivar(1:parmax)       = ivar_tmp(1:parmax)
    ipower(1:9, 1:parmax)= ipower_tmp(1:9, 1:parmax)
    parnam(1:parmax)     = parnam_tmp(1:parmax)

    deallocate(param_tmp, parnam_tmp, ipower_tmp, ivar_tmp)

    allocate(coord_tmp(10, enermax), energy_tmp(enermax), wt_in(enermax), stat=alloc)
    if (alloc /= 0) stop 'coord_tmp - out of memory'

    !---input energy values------!
    npts = 0
    yes  = .true.
    Nout = 1
    do while (yes)
      npts = npts + 1
      if (npts > enermax) then
        write(f_out, "('Too many ab initio points, increase enermax:',I8,' vs ',I8)") npts, enermax
        stop 'enermax exceeded'
      endif
      read (f_inp, *, end=199) coord_tmp(1:9, npts), energy_tmp(npts), wt_in(npts)
      if (wt_in(npts) < 0.0d0) wt_in(npts) = abs(wt_in(npts)) * factor_back
      if (coord_tmp(1, npts) < 0.0d0) yes = .false.
    enddo
199 continue
    npts = npts - 1
    write (f_out, "(' Number of ab initio points read: ', i8)") npts

    allocate(parold(parmax), parinitial(parmax), last_param(parmax), sterr(parmax), &
             Tsing(parmax), wspace(lwork), crr(parmax), al(parmax, parmax), bl(parmax), &
             dx(parmax), ai(parmax, parmax), rjacob(npts, parmax), coord(10, npts), &
             energy(npts), eps(npts), wt(npts), wt_tmp(npts), sigma(npts), sens(parmax), &
             dV(npts, parmax), stat=alloc)
    if (alloc /= 0) stop 'fit arrays - out of memory'

    coord(1:9, 1:npts) = coord_tmp(1:9, 1:npts)
    energy(1:npts)     = energy_tmp(1:npts)
    wt(1:npts)         = wt_in(1:npts)
    deallocate(coord_tmp, energy_tmp, wt_in)

    parinitial = param

    ! "normalising" the weight factors
    nused = 0
    wtsum = 0.0d0
    do i = 1, npts
      if (wt(i) > 0.0d0) nused = nused + 1
      wtsum = wtsum + wt(i)
    enddo
    wtsum = wtsum / dble(nused)
    wt(:) = wt(:) / wtsum

    ZPE       = 0.0d0
    energy(:) = (energy(:) - ZPE) * conv_energy

    do nrow = 1, npts
      if (wt(nrow) > 0.0d0) then
        sigma(nrow) = 1.0d0 / wt(nrow)
      else
        sigma(nrow) = 1000.0d0
      endif
    enddo

    still_run = .true.
    outer_loop: do while (still_run)

      if (itmax /= 0) then
        write(f_out, "(/3x,66('-'))")
        write(f_out, "('   |   iter  |  points |    deviat     |     rms       |  stability  |')")
        write(f_out, "(3x,66('-'))")
      endif

      stadev_old = 1.0d10
      stability  = 1.0d10
      stadev     = 1.0d10

      numpar  = 0
      maxivar = 0
      minivar = 1000000
      do i = 1, parmax
        if (ivar(i) > 0) numpar = numpar + 1
        if (ivar(i) > maxivar) maxivar = ivar(i)
        if (ivar(i) < minivar .and. ivar(i) /= 0) minivar = ivar(i)
      enddo

      if (nused < numpar) then
        write(f_out, "('warning! nused < numpar:',2I8)") nused, numpar
        stop 'nused < numpar'
      endif

      rjacob = 0.0d0
      iter   = 0
      sens   = 0.0d0
      do while (iter <= itmax .and. stadev > stadev_best .and. (stability >= stab_best))
        iter   = iter + 1
        ssq    = 0.0d0
        rms    = 0.0d0
        parold = param
        dV     = 0.0d0

        do nrow = 1, npts
          local(1:4) = coord(1:4, nrow) * conv_length
          local(5:9) = coord(5:9, nrow) * pi / 180.0d0
          call compute_all_dF(parmax, ipower, local, param(1:nparam_eq), dV(nrow, :))
        enddo

        do nrow = 1, npts
          local(1:4) = coord(1:4, nrow) * conv_length
          local(5:9) = coord(5:9, nrow) * pi / 180.0d0
          v = potential(parmax, dV(nrow, :), local, param, ipower, 0)

          if (itmax /= 0 .and. wt(nrow) > 0.0d0) then
            ncol = 0
            do i = 1, parmax
              if (ivar(i) /= 0) then
                ncol = ncol + 1
                rjacob(nrow, ncol) = potential(parmax, dV(nrow, :), local, param, ipower, i)
              endif
            enddo
          endif

          eps(nrow) = energy(nrow) - v
          ssq       = ssq + eps(nrow) * eps(nrow) * wt(nrow)
          rms       = rms + eps(nrow) * eps(nrow)
        enddo

        if (itmax /= 0) then
          do irow = 1, numpar
            do icolumn = 1, irow
              al(irow, icolumn) = sum(rjacob(:, icolumn) * rjacob(:, irow) * wt(:))
              al(icolumn, irow) = al(irow, icolumn)
            enddo
          enddo
          do irow = 1, numpar
            bl(irow) = sum(eps(:) * rjacob(:, irow) * wt(:))
          enddo

          select case (trim(fit_type))
          case default
            write (6, "('fit_type ',a,' unknown')") trim(fit_type)
            stop 'fit_type unknown'
          case ('linur')
            call linur(numpar, parmax, al, bl, dx, ierror)
            if (ierror /= 0) then
              ncol = 0
              do i = 1, parmax
                if (ivar(i) /= 0) then
                  ncol = ncol + 1
                  if (ncol == ierror) then
                    ivar(i)  = 0
                    param(i) = 0.0d0
                    if (ireplace_by /= 0) param(i) = parinitial(i)
                    write(f_out, "(i6,'-th is out (linur singular) - ',a10,9i3)") i, parnam(i), ipower(1:9, i)
                  endif
                endif
              enddo
              cycle outer_loop
            endif
          case ('dgelss')
            ai = al
            call dgelss(numpar, numpar, 1, ai(1:numpar, 1:numpar), numpar, &
                        bl(1:numpar), numpar, Tsing, 1.0d-12, rank0, wspace, lwork, info)
            if (info /= 0) then
              write(6, "('DGELSS: error ',i6)") info
              stop 'dgelss'
            endif
            dx = bl
          end select

          if (nused /= numpar) then
            stadev = dsqrt(ssq / dble(nused - numpar))
          else
            stadev = dsqrt(ssq / dble(nused))
          endif

          ncol = 0
          do i = 1, parmax
            if (ivar(i) /= 0) then
              ncol     = ncol + 1
              param(i) = param(i) + dx(ncol)
            endif
          enddo

          call invmat(al, ai, numpar, parmax)

          ncol = 0
          do i = 1, parmax
            if (ivar(i) /= 0) then
              ncol = ncol + 1
              if (nused == numpar) then
                sterr(ncol) = 0.0d0
              else
                sterr(ncol) = sqrt(abs(ai(ncol, ncol))) * stadev
              endif
            endif
          enddo

          stability  = abs((stadev - stadev_old) / stadev)
          stadev_old = stadev

          if (robustfit > 0) then
            do nrow = 1, npts
              if (wt(nrow) > 0.0d0) then
                wt(nrow) = 1.0d0 / (sigma(nrow)**2 + a_wats * eps(nrow)**2)
              endif
            enddo
            nused = 0
            wtsum = 0.0d0
            do i = 1, npts
              if (wt(i) > 0.0d0) nused = nused + 1
              wtsum = wtsum + wt(i)
            enddo
            wtsum = wtsum / dble(nused)
            wt(:) = wt(:) / wtsum
            ssq   = sum(eps(:)**2 * wt(:))
            if (nused == numpar) then
              stadev = dsqrt(ssq / dble(nused))
            else
              stadev = dsqrt(ssq / dble(nused - numpar))
            endif
          endif

        else
          stadev_old = stadev
          stadev     = dsqrt(ssq / dble(nused))
        endif

        rms = dsqrt(rms / dble(npts))

        if (itmax /= 0 .and. stadev < maxdeviat) then
          write (f_out, "('   |  ',i4,'   | ',i6,'  |  ',e12.5,' | ',e12.5,'  |  ',e10.3,' |')") &
                 iter, nused, stadev, rms, stability
        endif

      enddo ! --- iter

      if (robustfit > 0) write(f_out, "(/'Watson Alpha parameter',f8.4)") a_wats

      if (robustfit == 2) then
        robustfit = 0
      else if (robustfit == 3) then
        do i = 1, npts
          imin      = max(1, i - npts / ishift)
          imax      = min(npts, i + npts / ishift)
          ssq       = sum(wt(imin:imax)) / dble(imax - imin)
          wt_tmp(i) = ssq
        enddo
        wtsum     = sum(wt_tmp(:)) / dble(nused)
        wt(:)     = wt_tmp(:) / wtsum
        robustfit = 0
        cycle outer_loop
      endif

      ve = param(nparam_eq + 1)

      if (itmax /= 0 .and. find_best_params) then
        write (f_out, "(//'Potential parameters rounded in accord. with their standard errors'/)")
        l = 0
        do i = 1, parmax
          if (ivar(i) /= 0) then
            l        = l + 1
            ndigits  = 0
            conf_int = sterr(l)
            do while (conf_int <= 10.0d0 .and. ndigits < 12)
              ndigits  = ndigits + 1
              conf_int = conf_int * 10.0d0
            enddo
            write(fmt0, "(i2)") ndigits
            fmt0 = adjustl(fmt0)
            if (i <= nparam_eq) then
              fmt1 = '(a10,21x,I4,2x,f18.' // trim(fmt0) // ',''    ('',i8,'')'')'
              write (f_out, FMT=fmt1) trim(parnam(i)), ivar(i), param(i), nint(conf_int)
            else
              fmt1 = '(a2,9i3,I4,2x,f18.' // trim(fmt0) // ',''    ('',i8,'')'')'
              write (f_out, FMT=fmt1) trim(parnam(i)), ipower(1:9, i), ivar(i), param(i), nint(conf_int)
            endif
          else
            ndigits = 2
            if (param(i) /= 0.0d0) ndigits = 8
            write(fmt0, "(i2)") ndigits
            fmt0 = adjustl(fmt0)
            if (i <= nparam_eq) then
              fmt1 = '(a10,21x,I4,2x,f18.' // trim(fmt0) // ')'
              write (f_out, FMT=fmt1) trim(parnam(i)), ivar(i), param(i)
            else
              fmt1 = '(a2,9i3,I4,2x,f18.' // trim(fmt0) // ')'
              write (f_out, FMT=fmt1) trim(parnam(i)), ipower(1:9, i), ivar(i), param(i)
            endif
          endif
        enddo

        write (f_out, "('------------------------------')")
        write (f_out, "(a80)") (title(i), i = 1, 4)
        do i = 1, nparam_eq
          write (f_out, "(a10,21x,i4,2x,f18.8)") trim(parnam(i)), ivar(i), param(i)
        enddo
        do i = nparam_eq + 1, parmax
          write (f_out, "(a2,9i3,i4,2x,f18.8)") trim(parnam(i)), ipower(1:9, i), ivar(i), param(i)
        enddo
      endif

      if (itmax /= 0 .and. find_best_params) then
        if ((stadev - last_stadev) / stadev > Nout * factor_back .and. last_i /= 0) then
          ivar          = ivar0
          ivar(last_i)  = ivar(last_i) + 1
          param         = last_param
          last_stadev   = 1.0d7
          write(f_out, "(I6,'-th is back - ',A10)") last_i, parnam(last_i)
          cycle outer_loop
        endif

        write (f_out, "(/'Correlation coefficients larger than 0.9: ')")
        do i1 = 1, numpar
          do i2 = i1 + 1, numpar
            corr = ai(i1, i2) / sqrt(abs(ai(i1, i1) * ai(i2, i2)))
            if (abs(corr) > 0.9d0) write (f_out, "(10x,'corr(',i3,',',i3,') = ',f10.7)") i1, i2, abs(corr)
          enddo
        enddo

        ncol = 0
        do i = 1, parmax
          if (ivar(i) /= 0) then
            ncol = ncol + 1
            if (nused /= numpar) then
              ssq        = sum(rjacob(:, ncol)**2 * wt(:))
              sens(ncol) = stadev * 0.1d0 / dble(numpar) / dsqrt(ssq / dble(nused))
            endif
          endif
        enddo

        write (f_out, "(/15x,80('-')/15x,':     old parm   :     new parm   :  delta parm  :   std. err.  : sensitivity'/15x,80('-'))")
        l = 0
        do i = 1, parmax
          if (ivar(i) /= 0) then
            l      = l + 1
            epsil  = param(i) - parold(i)
            crr(l) = (sterr(l) * al(l, l) / stadev)**2
            mark   = ' '
            if (abs(sens(l)) >= abs(epsil)) mark = '!'
            if (abs(param(i)) > abs(sterr(l))) then
              write (f_out, "(i3,2x,a10,1x,2(':',f14.2,2x),':',3(f14.3,':'),a1)") &
                     l, parnam(i), parold(i), param(i), epsil, sterr(l), sens(l), mark
            else
              write (f_out, "(i3,2('#'),a10,1x,2(':',f14.2,2x),':',3(f14.3,':'),a1,'#')") &
                     l, parnam(i), parold(i), param(i), epsil, sterr(l), sens(l), mark
            endif
          endif
        enddo

        tivar      = minivar
        ivar0      = ivar
        icorr      = 0
        last_param = param
        Nout       = 0
        do while (tivar < maxivar .and. factor_out > 0.0d0)
          relat_sterr_min = 1.0d20
          ilow_sterr      = -1
          relat_sterr     = 0.0d0
          ilarge_sterr    = -1

          l = numpar + 1
          do i = parmax, 1, -1
            if (ivar(i) /= 0) then
              l = l - 1
              if (i > nparam_eq .and. abs(param(i)) > too_big) then
                last_param  = param
                last_stadev = stadev
                last_i      = i
                ivar(i)     = 0
                param(i)    = 0.0d0
                write(f_out, "(' parameter ',i6,' is out (too_big) - ',A10)") i, parnam(i)
                cycle outer_loop
              endif

              if (param(i) /= 0.0d0 .and. ivar(i) == tivar) then
                if (icorr(i) == 0 .and. factor_out < abs(sterr(l) / param(i))) then
                  icorr(i)     = 1
                  l0           = 0
                  relat_sterr  = abs(sterr(l) / param(i))
                  ilarge_sterr = i
                  do i0 = 1, parmax
                    if (ivar(i0) /= 0) then
                      l0   = l0 + 1
                      corr = abs(ai(l, l0)) / sqrt(abs(ai(l, l) * ai(l0, l0)))
                      if (i0 /= i .and. param(i0) /= 0.0d0 .and. ivar(i0) == tivar) then
                        if (icorr(i0) == 0 .and. corr > fact_corr) then
                          icorr(i0) = 1
                          if (relat_sterr < abs(sterr(l0) / param(i0))) then
                            relat_sterr  = abs(sterr(l0) / param(i0))
                            ilarge_sterr = i0
                          endif
                        endif
                      endif
                    endif
                  enddo
                  ivar(ilarge_sterr)  = 0
                  param(ilarge_sterr) = 0.0d0
                  if (ireplace_by /= 0) param(ilarge_sterr) = parinitial(ilarge_sterr)
                  Nout = Nout + 1
                  if (relat_sterr_min > relat_sterr) then
                    relat_sterr_min = relat_sterr
                    ilow_sterr      = ilarge_sterr
                  endif
                endif
              endif
            endif
          enddo

          if (relat_sterr_min > factor_out .and. ilow_sterr /= -1) then
            last_ivar   = ilow_sterr
            last_stadev = stadev
            last_i      = ilow_sterr
            write(f_out, "(I6,' parameters are out')") Nout
            cycle outer_loop
          endif
          tivar = tivar + 1
        enddo

        still_run = .false.
      endif
      still_run = .false.
    enddo outer_loop

    !--- Write 00.res ---!
    inquire(f_res, opened=ifopen)
    if (ifopen) then
      rewind(f_res)
    else
      open(f_res, file='00.res', status='replace')
    endif

    if (itmax /= 0) then
      write (f_res, "(a80)") (title(i), i = 1, 4)
      do i = 1, nparam_eq
        write (f_res, "(a10,21x,i4,2x,f18.8)") parnam(i), ivar(i), param(i)
      enddo
      do i = nparam_eq + 1, parmax
        write (f_res, "(a2,9i3,i4,2x,f18.8)") trim(parnam(i)), ipower(1:9, i), ivar(i), param(i)
      enddo
    endif

    write (f_res, "('------------------------------')")
    write (f_res, "(a80)") (title(i), i = 1, 4)
    do i = 1, nparam_eq
      write (f_res, "(a10,27x,i4,2x,f18.8)") parnam(i), ivar(i), param(i)
    enddo
    do i = nparam_eq + 1, parmax
      write (f_res, "(a2,9i3,i4,2x,f18.8)") trim(parnam(i)), ipower(1:9, i), ivar(i), param(i)
    enddo

    do nrow = 1, npts
      v = energy(nrow) - eps(nrow)
      write (f_res, "(4f15.6,1x,5f11.4,3f18.7,2x,e9.2)") &
            coord(1:9, nrow), energy(nrow) - ve, v - ve, eps(nrow), wt(nrow)
    enddo

    write (f_res, "(/3x,59('-')/'   |  param  |  points |   deviat    |    rms      | stability |')")
    write (f_res, "(3x,59('-'))")
    write (f_res, "('   |  ',i4,'   | ',i6,'  |  ',e12.5,' | ',e12.5,'  |  ',e11.3,' |')") &
          numpar, nused, stadev, rms, stability
    write (f_res, "(/92('#'))")

    !--- Write TROVE POTEN block (ch3br_poten.trove) ---!
    open(f_trove, file='ch3br_poten.trove', status='replace')
    write(f_trove, "('POTEN')")
    write(f_trove, "('NPARAM   ',i0)") parmax
    write(f_trove, "('POT_TYPE  poten_zxy3_sym')")
    write(f_trove, "('COEFF  powers (powers or list)')")
    do i = 1, nparam_eq
      write(f_trove, "(a5,11i3,4x,f16.8)") parnam(i)(1:5), 0,0,0,0,0,0,0,0,0,0, ivar(i), param(i)
    enddo
    do i = nparam_eq + 1, parmax
      write(f_trove, "(a2,9i3,i2,i3,2x,f18.8)") trim(parnam(i)), ipower(1:9, i), 0, ivar(i), param(i)
    enddo
    write(f_trove, "('end')")
    close(f_trove)

    write (f_out, "(/3x,59('-')/'   |  param  |  points |   deviat    |    rms      | stability |')")
    write (f_out, "(3x,59('-'))")
    write (f_out, "('   |  ',i4,'   | ',i6,'  |  ',e12.5,' | ',e12.5,'  |  ',e11.3,' |')") &
          numpar, nused, stadev, rms, stability
    write (f_out, "(/92('#'))")

  end subroutine do_fitting

  recursive function potential(parmax, dF, local, param, ipower, ipar) result (V)
    integer, intent(in)          :: parmax
    double precision, intent(in) :: dF(parmax), local(18), param(parmax)
    integer, intent(in)          :: ipar, ipower(18, parmax)

    double precision             :: param_(parmax), dF_(parmax)
    double precision             :: deltax, potright, potleft, V
    integer                      :: i, iparam
    double precision, parameter  :: factordeltax = 0.001d0

    V = 0.0d0
    if (ipar == 0) then
      do i = nparam_eq + 1, parmax
        V = V + param(i) * dF(i)
      enddo
    elseif (ipar > nparam_eq) then
      V = dF(ipar)
    else
      deltax = factordeltax * abs(param(ipar))
      if (deltax <= 1.0d-15) deltax = 1.0d-6

      param_       = param
      param_(ipar) = param(ipar) + deltax
      call compute_all_dF(parmax, ipower, local, param_(1:nparam_eq), dF_)
      potright = potential(parmax, dF_, local, param_, ipower, 0)

      param_(ipar) = param(ipar) - deltax
      call compute_all_dF(parmax, ipower, local, param_(1:nparam_eq), dF_)
      potleft = potential(parmax, dF_, local, param_, ipower, 0)

      V = (potright - potleft) / (2.0d0 * deltax)
    endif
  end function potential

  subroutine compute_all_dF(parmax, ipower, local, param, dF_row)
    implicit none
    integer, intent(in)           :: parmax, ipower(18, parmax)
    double precision, intent(in)  :: local(18), param(5)
    double precision, intent(out) :: dF_row(parmax)
    double precision              :: r1, r2, r3, r4, deg
    double precision              :: r1e, r2e, a, b, betae, pi
    double precision              :: xi(9), chi(9, 6), tau12, tau13, tau23, term
    integer                       :: ioper, iparam

    pi  = 4.0d0 * datan2(1.0d0, 1.0d0)
    deg = pi / 180.0d0

    r1 = local(1)
    r2 = local(2)
    r3 = local(3)
    r4 = local(4)

    r1e   = param(1)
    r2e   = param(2)
    betae = param(3) * deg
    a     = param(4)
    b     = param(5)

    xi(1) = 1.0d0 - exp(-a * (r1 - r1e))
    xi(2) = 1.0d0 - exp(-b * (r2 - r2e))
    xi(3) = 1.0d0 - exp(-b * (r3 - r2e))
    xi(4) = 1.0d0 - exp(-b * (r4 - r2e))

    xi(5) = local(5) - betae
    xi(6) = local(6) - betae
    xi(7) = local(7) - betae
    tau12 = local(8)
    tau13 = local(9)
    tau23 = 2.0d0 * pi - tau12 - tau13

    xi(8) = (2.0d0 * tau23 - tau13 - tau12) / sqrt(6.0d0)
    xi(9) = (tau12 - tau13) / sqrt(2.0d0)

    do ioper = 1, 6
      call ML_symmetry_transformation_XY3_II(ioper, xi, chi(:, ioper))
    enddo

    dF_row(1:nparam_eq) = 0.0d0
    do iparam = nparam_eq + 1, parmax
      term = 0.0d0
      do ioper = 1, 6
        term = term + product(chi(1:9, ioper)**ipower(1:9, iparam))
      enddo
      dF_row(iparam) = term
    enddo
  end subroutine compute_all_dF

  double precision function diff_V(ipower, local, param)
    implicit none
    integer, intent(in)          :: ipower(18)
    double precision, intent(in) :: local(18), param(5)
    double precision             :: r1, r2, r3, r4, deg
    double precision             :: r1e, r2e, a, b, betae, pi
    double precision             :: xi(9), chi(9), tau12, tau13, tau23, term
    integer                      :: ioper

    pi  = 4.0d0 * datan2(1.0d0, 1.0d0)
    deg = pi / 180.0d0

    r1 = local(1)
    r2 = local(2)
    r3 = local(3)
    r4 = local(4)

    r1e   = param(1)
    r2e   = param(2)
    betae = param(3) * deg
    a     = param(4)
    b     = param(5)

    xi(1) = 1.0d0 - exp(-a * (r1 - r1e))
    xi(2) = 1.0d0 - exp(-b * (r2 - r2e))
    xi(3) = 1.0d0 - exp(-b * (r3 - r2e))
    xi(4) = 1.0d0 - exp(-b * (r4 - r2e))

    xi(5) = local(5) - betae
    xi(6) = local(6) - betae
    xi(7) = local(7) - betae
    tau12 = local(8)
    tau13 = local(9)
    tau23 = 2.0d0 * pi - tau12 - tau13

    xi(8) = (2.0d0 * tau23 - tau13 - tau12) / sqrt(6.0d0)
    xi(9) = (tau12 - tau13) / sqrt(2.0d0)

    term = 0.0d0
    do ioper = 1, 6
      call ML_symmetry_transformation_XY3_II(ioper, xi, chi)
      term = term + product(chi(1:9)**ipower(1:9))
    enddo
    diff_V = term
  end function diff_V

  subroutine skiplines(inpunit, n)
    integer, intent(in) :: n, inpunit
    character(len=80)   :: label
    integer             :: i0
    do i0 = 1, n
      read (inpunit, "(a80)") label
    enddo
  end subroutine skiplines

  subroutine linur(dimen, npar, coeff, constant, solution, error)
    implicit none
    integer, intent(in)           :: dimen, npar
    integer, intent(out)          :: error
    double precision, intent(in)  :: coeff(npar, npar), constant(npar)
    double precision, intent(out) :: solution(npar)
    double precision              :: a0(npar, npar), c
    integer                       :: i1, i2, i, k8, k, l, k9

    do i1 = 1, dimen
      do i2 = 1, dimen
        a0(i1, i2) = coeff(i1, i2)
      enddo
    enddo
    do i = 1, dimen
      solution(i) = constant(i)
    enddo
    error = 0
    do i = 1, dimen
      c  = 0.0d0
      k8 = i - 1
      do k = 1, k8
        c = c + a0(k, i) * a0(k, i)
      enddo
      if (c >= a0(i, i)) then
        error = i
        return
      endif
      a0(i, i) = sqrt(a0(i, i) - c)
      if (a0(i, i) == 0.0d0) then
        error = i
        return
      endif
      k8 = i + 1
      do l = k8, dimen
        k9 = i - 1
        c  = 0.0d0
        do k = 1, k9
          c = c + a0(k, i) * a0(k, l)
        enddo
        a0(i, l) = (a0(l, i) - c) / a0(i, i)
      enddo
    enddo
    do i = 1, dimen
      k8 = i - 1
      c  = 0.0d0
      do k = 1, k8
        c = c + a0(k, i) * solution(k)
      enddo
      solution(i) = (solution(i) - c) / a0(i, i)
    enddo
    do i1 = 1, dimen
      i  = 1 + dimen - i1
      k8 = i + 1
      c  = 0.0d0
      do k = k8, dimen
        c = c + a0(i, k) * solution(k)
      enddo
      solution(i) = (solution(i) - c) / a0(i, i)
    enddo
  end subroutine linur

  subroutine invmat(al, ai, dimen, npar)
    implicit none
    integer, intent(in)           :: npar, dimen
    double precision, intent(in)  :: al(npar, npar)
    double precision, intent(out) :: ai(npar, npar)
    double precision              :: h(npar), p, q
    integer                       :: i1, k, i, j, k8, k9

    ai(1:dimen, 1:dimen) = al(1:dimen, 1:dimen)
    do i1 = 1, dimen
      k = dimen - i1 + 1
      p = ai(1, 1)
      do i = 2, dimen
        q    = ai(i, 1)
        h(i) = q / p
        if (i <= k) h(i) = -q / p
        do j = 2, i
          k8         = i - 1
          k9         = j - 1
          ai(k8, k9) = ai(i, j) + q * h(j)
        enddo
      enddo
      ai(dimen, dimen) = 1.0d0 / p
      do i = 2, dimen
        k8            = i - 1
        ai(dimen, k8) = h(i)
      enddo
    enddo
    do i = 1, dimen
      k8 = i - 1
      do j = 1, k8
        ai(j, i) = ai(i, j)
      enddo
    enddo
  end subroutine invmat

  subroutine ML_symmetry_transformation_XY3_II(ioper, src, dst)
    implicit none
    integer, intent(in)           :: ioper
    double precision, intent(in)  :: src(1:9)
    double precision, intent(out) :: dst(1:9)
    double precision              :: repres(6, 9, 9), a, b

    a = 0.5d0
    b = 0.5d0 * sqrt(3.0d0)
    repres = 0.0d0

    repres(1,1,1) = 1.0d0; repres(1,2,2) = 1.0d0; repres(1,3,3) = 1.0d0
    repres(1,4,4) = 1.0d0; repres(1,5,5) = 1.0d0; repres(1,6,6) = 1.0d0
    repres(1,7,7) = 1.0d0; repres(1,8,8) = 1.0d0; repres(1,9,9) = 1.0d0

    repres(2,1,1) = 1.0d0
    repres(2,2,4) = 1.0d0; repres(2,3,2) = 1.0d0; repres(2,4,3) = 1.0d0
    repres(2,5,7) = 1.0d0; repres(2,6,5) = 1.0d0; repres(2,7,6) = 1.0d0
    repres(2,8,8) = -a;    repres(2,8,9) =  b
    repres(2,9,8) = -b;    repres(2,9,9) = -a

    repres(3,1,1) = 1.0d0
    repres(3,2,3) = 1.0d0; repres(3,3,4) = 1.0d0; repres(3,4,2) = 1.0d0
    repres(3,5,6) = 1.0d0; repres(3,6,7) = 1.0d0; repres(3,7,5) = 1.0d0
    repres(3,8,8) = -a;    repres(3,8,9) = -b
    repres(3,9,8) =  b;    repres(3,9,9) = -a

    repres(4,1,1) = 1.0d0
    repres(4,2,2) = 1.0d0; repres(4,3,4) = 1.0d0; repres(4,4,3) = 1.0d0
    repres(4,5,5) = 1.0d0; repres(4,6,7) = 1.0d0; repres(4,7,6) = 1.0d0
    repres(4,8,8) = 1.0d0; repres(4,9,9) = -1.0d0

    repres(5,1,1) = 1.0d0
    repres(5,2,4) = 1.0d0; repres(5,3,3) = 1.0d0; repres(5,4,2) = 1.0d0
    repres(5,5,7) = 1.0d0; repres(5,6,6) = 1.0d0; repres(5,7,5) = 1.0d0
    repres(5,8,8) = -a;    repres(5,8,9) =  b
    repres(5,9,8) =  b;    repres(5,9,9) =  a

    repres(6,1,1) = 1.0d0
    repres(6,2,3) = 1.0d0; repres(6,3,2) = 1.0d0; repres(6,4,4) = 1.0d0
    repres(6,5,6) = 1.0d0; repres(6,6,5) = 1.0d0; repres(6,7,7) = 1.0d0
    repres(6,8,8) = -a;    repres(6,8,9) = -b
    repres(6,9,8) = -b;    repres(6,9,9) =  a

    if (ioper < 1 .or. ioper > 6) stop 'ML_symmetry_transformation_XY3_II - bad operation'
    dst = matmul(repres(ioper, :, :), src)
  end subroutine ML_symmetry_transformation_XY3_II

end module FITPES_ch3br

program driver
  use FITPES_ch3br
  implicit none
  call do_fitting
end program driver
