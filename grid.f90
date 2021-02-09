module grid

use accuracy
use fields
use splines
implicit none


type grid_1d_type
  character(cl)            :: qtype
  integer(ik)              :: polyad_number
  integer(ik)              :: maxorder
  integer(ik), allocatable :: nqpoints(:)
  real(rk), allocatable    :: qpoints(:,:)
  real(rk), allocatable    :: qweights(:,:)
  real(rk), allocatable    :: xme(:,:)
  real(rk)                 :: xborders(2)
  integer(ik)              :: phi_npoints
  real(rk), allocatable    :: phi_val(:,:)
  real(rk), allocatable    :: phi_xval(:)
end type grid_1d_type

integer(ik)                     :: grid_nmodes
type(grid_1d_type), allocatable :: grid_1d(:)


contains


subroutine init_grid_1d

  integer(ik) :: nmodes, info, imode, iorder, maxorder, maxnpoints, quad_order_nd, npoints_nd, ibas, jbas, nbas, &
                 phi_npoints, lwork, liwork, ipoint, npoints, nelem1, nelem2, ideg, ibastype, nbas_
  integer(ik), allocatable :: iwork(:)
  real(rk) :: x, step
  real(rk), allocatable :: work(:), xval(:), xvec(:,:), phival(:,:), phival_dvr(:,:)
  real(rk), allocatable :: points_nd(:,:), weights_nd(:)

  write(out, '(/a)') 'init_grid_1d: Initialize 1D quadrature formulas'

  nmodes = trove%Nmodes
  grid_nmodes = nmodes

  allocate(grid_1d(nmodes), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'init_grid_1d error: failed to allocate grid_1d(nmodes)', 'nmodes =', nmodes
    stop
  endif


  do imode=1, nmodes

    ibastype = job%bset(imode)%species                 ! index of 1D basis set
    nbas = bset%bs1D(ibastype)%Size + 1                ! number of 1D basis functions
    phi_npoints = size(bset%bs1D(ibastype)%phi,dim=1)  ! number of grid points 1D basis functions computed at

    maxorder = nbas   ! max order of 1D quadrature

    allocate(grid_1d(imode)%nqpoints(maxorder), &
             grid_1d(imode)%xme(maxorder,maxorder), &
             grid_1d(imode)%qpoints(maxorder,maxorder), &
             grid_1d(imode)%qweights(maxorder,maxorder), &
             grid_1d(imode)%phi_val(phi_npoints,nbas), &
             grid_1d(imode)%phi_xval(phi_npoints), stat=info)

    if (info/=0) then
      write(out, '(/a/a/a,10(1x,i6))') &
      'init_grid_1d error: failed to allocate grid_1d(imode)%nqpoints(maxorder), grid_1d(imode)%qpoints(maxorder,maxorder), grid_1d(imode)%qweights(maxorder,maxorder),', &
      'grid_1d(imode)%xme(maxorder,maxorder), grid_1d(imode)%phi_val(phi_npoints,nbas), grid_1d(imode)%phi_xval(phi_npoints)', &
      'imode, maxorder, nbas, phi_npoints =', imode, maxorder, nbas, phi_npoints
      stop
    endif

    grid_1d(imode)%maxorder = maxorder
    grid_1d(imode)%nqpoints(1:maxorder) = (/( iorder, iorder=1,maxorder )/)   ! number of quadrature points for each order=1..maxorder
    grid_1d(imode)%polyad_number = job%bset(imode)%res_coeffs

    ! matrix elements of X
    ideg = 1
    do ibas=1, nbas
      do jbas=1, nbas
        grid_1d(imode)%xme(jbas,ibas) = bset%bs1D(ibastype)%matelements(-1,ideg,jbas-1,ibas-1)   ! basis set indices in "matelements" begin from zero
      enddo                                                                                      ! "matelements(-1,.." returns matrix elements of coordinates chosen for expansion of KEO
    enddo                                                                                        ! "matelements(-1,ideg..." returns matrix elements of coordinates to the "ideg"-th power

    ! values of coordinate and basis functions on a grid
    grid_1d(imode)%xborders(1:2) = job%bset(imode)%borders(1:2)
    grid_1d(imode)%phi_npoints = phi_npoints
    step = (grid_1d(imode)%xborders(2)-grid_1d(imode)%xborders(1))/real(phi_npoints-1,kind=rk)
    do ipoint=1, phi_npoints
      x = grid_1d(imode)%xborders(1)+real(ipoint-1,kind=rk)*step
      grid_1d(imode)%phi_xval(ipoint) = x
      grid_1d(imode)%phi_val(ipoint,1:nbas) = bset%bs1D(ibastype)%phi(ipoint-1,0:nbas-1)
    enddo

    ! generate quadrature points and weights

    write(out, '(/1x,a,1x,i3/1x,a)') 'quadrature points (x) and weights(w) for imode =', imode, 'iorder'

    ! allocate arrays to keep eigenvalues and eigenvectors of X, values of basis functions on grid, as well as work arrays for LAPACK diagonalizer DSPEVD
    maxnpoints = maxval(grid_1d(imode)%nqpoints(1:maxorder))
    nelem1 = 2*maxnpoints*maxnpoints+6*maxnpoints+1               ! recommended sizes of work and iwork arrays
    nelem2 = 5*maxnpoints+3                                       ! for DSYEVD procedure
    allocate(xval(maxnpoints), xvec(maxnpoints,maxnpoints), phival(maxnpoints,nbas), phival_dvr(maxnpoints,maxnpoints), work(nelem1), iwork(nelem2), stat=info)
    if (info/=0) then
      write(out, '(/a/a,10(1x,i8))') &
      'init_grid_1d error: failed to allocate xval(maxnpoints), xvec(maxnpoints,maxnpoints), phival(maxnpoints,nbas), phival_dvr(maxnpoints,maxnpoints), work(nelem1), iwork(nelem2)', &
      'maxnpoints, nbas, nelem1, nelem2 =', maxnpoints, nbas, nelem1, nelem2
      stop
    endif

    do iorder=1, maxorder

      npoints = grid_1d(imode)%nqpoints(iorder)

      ! diagonalize X to get quadrature points
      liwork = size(iwork)
      lwork = size(work)
      xvec(1:npoints,1:npoints) = grid_1d(imode)%xme(1:npoints,1:npoints)
      call dsyevd('V', 'L', npoints, xvec(1:npoints,1:npoints), npoints, xval(1:npoints), work, lwork, iwork, liwork, info)
      if (info/=0) then
        write(out, '(/a,1x,i6,1x,a,1x,i3,1x,a,1x,i3)') &
        'init_grid_1d error: DSYEVD failed to diagonalize X, npoints =', npoints, ', imode =', imode, ', info =', info
        stop
      endif

      nbas_ = npoints

      ! compute values of basis functions at quadrature points
      do ibas=1, nbas_
        phival(1:npoints,ibas) = spline3(grid_1d(imode)%phi_xval(1:phi_npoints), grid_1d(imode)%phi_val(1:phi_npoints,ibas), xval(1:npoints))
      enddo

      ! transform basis to DVR
      phival_dvr(1:npoints,1:nbas_) = matmul(phival(1:npoints,1:nbas_), xvec(1:nbas_,1:nbas_))

      ! init quadrature points and weights
      grid_1d(imode)%qpoints(:,iorder) = 0.0
      grid_1d(imode)%qweights(:,iorder) = 0.0
      grid_1d(imode)%qpoints(1:npoints,iorder) = xval(1:npoints)
      grid_1d(imode)%qweights(1:npoints,iorder) = (/( 1.0_rk/phival_dvr(ipoint,ipoint)**2, ipoint=1, npoints )/)

      write(out, '(3x,i4,5x,a,1x,100(1x,es16.8))') iorder, 'x:', grid_1d(imode)%qpoints(1:npoints,iorder)
      write(out, '(12x,a,1x,100(1x,es16.8))') 'w:', grid_1d(imode)%qweights(1:npoints,iorder)

    enddo ! iorder

    deallocate(xval, xvec, phival, phival_dvr, work, iwork)

  enddo ! imode


  quad_order_nd = 20
  call init_grid_nd(quad_order_nd, npoints_nd, points_nd, weights_nd)

  write(out, '(/a)') 'init_grid_1d: Done'

end subroutine init_grid_1d



subroutine init_grid_nd(quad_order_nd, npoints_nd, points_nd, weights_nd)

  integer(ik), intent(in) :: quad_order_nd
  integer(ik), intent(out) :: npoints_nd
  real(rk), allocatable, intent(out) :: points_nd(:,:), weights_nd(:)

  integer(ik), parameter :: npoints_incr = 10000
  integer(ik) :: imode, nmodes, i, smolyak_nprod, iprod, ipoint, ipoint_, prod_nelem, info, iorder, maxorder, order, &
                 maxnpoints, maxnpoints_nd, np, mode_weights(grid_nmodes), ngrids_1d(grid_nmodes), ii(grid_nmodes), npoints(grid_nmodes)
  integer(ik), allocatable :: orders_1d(:,:), smolyak_ind(:,:), ind_points(:,:), prod_ind(:,:), smolyak_coef(:)
  real(rk), allocatable :: tmp_points(:,:), tmp_weights(:)

  write(out, '(/a)') 'init_grid_nd: Generate N-dimensional sparse grid'

  write(out, '(/1x,a,1x,i4)') 'quadrature order:', quad_order_nd
  write(out, '(1x,a,1x,i4)') 'min order:', sum(grid_1d(:)%polyad_number)

  if (.not.allocated(grid_1d)) then
    write(out, '(/a)') 'init_grid_nd error: 1D quadrature formulas are not initialized'
    stop
  endif

  nmodes = grid_nmodes
  maxorder = maxval((/(grid_1d(imode)%maxorder, imode=1, nmodes)/))
  maxnpoints = maxval((/(maxval(grid_1d(imode)%nqpoints(:)), imode=1, nmodes)/))

  allocate(orders_1d(maxorder,nmodes), ind_points(maxnpoints,nmodes), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i8))') 'init_grid_nd error: failed to allocate orders_1d(maxorder,nmodes) ind_points(maxnpoints,nmodes)', &
   'maxorder, maxnpoints, nmodes =', maxorder, maxnpoints, nmodes
    stop
  endif

  do imode=1, nmodes
    maxorder = grid_1d(imode)%maxorder
    ngrids_1d(imode) = maxorder
    orders_1d(1:maxorder,imode) = (/(iorder, iorder=1, maxorder)/)
  enddo

  do imode=1, nmodes
    mode_weights(imode) = grid_1d(imode)%polyad_number
  enddo

  ! estimate number of elements in Smolyak sum

  call smolyak_product(nmodes, ngrids_1d(1:nmodes), orders_1d(1:maxval(ngrids_1d),1:nmodes), mode_weights(1:nmodes), quad_order_nd, 1, ii, smolyak_nprod)

  ! allocate arrays to keep Smolyak coefficients and indices of 1D quadratures

  allocate(smolyak_ind(nmodes,smolyak_nprod), smolyak_coef(smolyak_nprod), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i8))') 'init_grid_nd error: failed to allocate smolyak_ind(nmodes,smolyak_nprod), smolyak_coef(smolyak_nprod)', &
    'smolyak_nprod, nmodes =', smolyak_nprod, nmodes
    stop
  endif

  ! compute Smolyak coefficients and indices of 1D quadratures

  call smolyak_product(nmodes, ngrids_1d(1:nmodes), orders_1d(1:maxval(ngrids_1d),1:nmodes), mode_weights(1:nmodes), quad_order_nd, 1, ii, smolyak_nprod, smolyak_ind, smolyak_coef)

  write(out, '(/1x,a,1x,i6,1x,a)') 'tensor products in Smolyak sum (', smolyak_nprod, '):'
  if (smolyak_nprod<=1000) then
    write(out, '(2x,a,7x,a,8x,a)') 'iprod', 'coef', 'orders of 1D quadratures'
    do iprod=1, smolyak_nprod
      write(out, '(1x,i6,5x,i6,5x,<nmodes>(1x,i3))') iprod, smolyak_coef(iprod), smolyak_ind(1:nmodes,iprod)
    enddo
  else
    write(out, '(1x,a)') 'too many elements in Smolyak sum, skip printing'
  endif

  ! generate sparse grid

  npoints_nd = 0
  maxnpoints_nd = 0

  if (allocated(points_nd)) deallocate(points_nd)
  if (allocated(weights_nd)) deallocate(weights_nd)

  do iprod=1, smolyak_nprod

    ! for each mode generate set of 1D-quadratures

    do imode=1, nmodes

      ! order of 1D quadrature for mode=imode
      order = smolyak_ind(imode,iprod)

      ! find index of 1D quadrature for given order
      iorder = 0
      do i=1, ngrids_1d(imode)
        if (orders_1d(i,imode)==order) then
          iorder = i
          exit
        endif
      enddo
      if (iorder<=0) then
        write(out, '(/a,1x,i3,1x,a,1x,i3)') 'init_grid_nd error: could not locate 1D quadrature of', order, 'order for mode #', imode
        write(out, '(a,1x,i3,1x,a,100(1x,i3))') 'list of available quadratures for mode #', imode, ', orders =', orders_1d(1:ngrids_1d(imode),imode)
        stop
      endif

      ! number and indices of points in "iorder"-th 1D quadrature
      npoints(imode) = grid_1d(imode)%nqpoints(iorder)
      ind_points(1:npoints(imode),imode) = (/(ipoint, ipoint=1, npoints(imode))/)

    enddo ! imode

    ! generate tensor product of 1D quadratures

    call direct_product(nmodes, npoints, ind_points(1:maxval(npoints),1:nmodes), 1, ii, prod_nelem)

    allocate(prod_ind(nmodes,prod_nelem), stat=info)
    if (info/=0) then
      write(out, '(/a/a,10(1x,i8))') 'init_grid_nd error: failed to allocate prod_ind(nmodes,prod_nelem)', 'nmodes, prod_nelem =', nmodes, prod_nelem
      stop
    endif

    call direct_product(nmodes, npoints, ind_points(1:maxval(npoints),1:nmodes), 1, ii, prod_nelem, prod_ind)

    ! add points and weights

    do ipoint=1, prod_nelem

      ipoint_ = npoints_nd + ipoint

      ! increase sizes of output arrays "npoints_nd" and "weights_nd"
      if (ipoint_>maxnpoints_nd) then
        if (allocated(points_nd)) then
          call move_alloc(points_nd, tmp_points)
          call move_alloc(weights_nd, tmp_weights)
        endif
        maxnpoints_nd = maxnpoints_nd + max(npoints_incr,prod_nelem)
        allocate(points_nd(nmodes,maxnpoints_nd), weights_nd(maxnpoints_nd), stat=info)
        if (info/=0) then
          write(out, '(/a/a,10(1x,i8))') 'init_grid_nd error: failed to allocate points_nd(nmodes,maxnpoints_nd), weights_nd(maxnpoints_nd)', &
          'nmodes, maxnpoints_nd =', nmodes, maxnpoints_nd
          stop
        endif
        points_nd = 0
        weights_nd = 0.0
        if (allocated(tmp_points)) then
          np = size(tmp_weights)
          points_nd(1:nmodes,1:np) = tmp_points(1:nmodes,1:np)
          weights_nd(1:np) = tmp_weights(1:np)
          deallocate(tmp_points,tmp_weights)
        endif
      endif

      points_nd(1:nmodes,ipoint_) = (/( grid_1d(imode)%qpoints(prod_ind(imode,ipoint),iorder), imode=1, nmodes )/)
      weights_nd(ipoint_) = product( (/( grid_1d(imode)%qweights(prod_ind(imode,ipoint),iorder), imode=1, nmodes )/) ) * real(smolyak_coef(iprod),rk)

      !write(out, '(1x,i8,3x,es16.8,3x,<nmodes>(1xes16.8))') ipoint_, weights_nd(ipoint_), points_nd(1:nmodes,ipoint_)

    enddo

    npoints_nd = npoints_nd + prod_nelem

    deallocate(prod_ind)

  enddo ! iprod

  write(out, '(/1x,a,1x,i8)') 'number of points in sparse grid:', npoints_nd

  deallocate(smolyak_ind, smolyak_coef, orders_1d, ind_points)


  ! adjust sizes (cut tails) of output arrays "points_nd" and "weights_nd"
  if (maxnpoints_nd/=npoints_nd) then
    call move_alloc(points_nd, tmp_points)
    call move_alloc(weights_nd, tmp_weights)
    allocate(points_nd(nmodes,npoints_nd), weights_nd(npoints_nd), stat=info)
    if (info/=0) then
      write(out, '(/a/a,10(1x,i8))') 'init_grid_nd error: failed to allocate points_nd(nmodes,npoints_nd), weights_nd(npoints_nd)', &
      'nmodes, npoints_nd =', nmodes, npoints_nd
      stop
    endif
    points_nd(1:nmodes,1:npoints_nd) = tmp_points(1:nmodes,1:npoints_nd)
    weights_nd(1:npoints_nd) = tmp_weights(1:npoints_nd)
    deallocate(tmp_points,tmp_weights)
  endif

  write(out, '(/a)') 'init_grid_nd: Done'

end subroutine init_grid_nd



! Returns Cartesian product of "nvec" integer vectors
!
recursive subroutine direct_product(nvec, nelem, elem, ivec, ind, prod_nelem, prod_ind)

  integer(ik), intent(in)              :: nvec
  integer(ik), intent(in)              :: nelem(nvec)
  integer(ik), intent(in)              :: elem(maxval(nelem),nvec)
  integer(ik), intent(in)              :: ivec
  integer(ik), intent(inout)           :: ind(nvec)
  integer(ik), intent(inout)           :: prod_nelem
  integer(ik), intent(inout), optional :: prod_ind(:,:)

  integer(ik) :: ielem

  if (ivec==1) prod_nelem = 0

  do ielem=1, nelem(ivec)
    ind(ivec) = elem(ielem,ivec)
    if (ivec==nvec) then
      prod_nelem = prod_nelem + 1
      if (present(prod_ind)) prod_ind(:,prod_nelem) = ind
    else
      call direct_product(nvec, nelem, elem, ivec+1, ind, prod_nelem, prod_ind)
    endif
  enddo

end subroutine direct_product



! Returns set of coefficients "smolyak_coef" and 1D-quadrature orders "smolyak_ord" for Smolyak sparse grid of order="smolyak_order"
!
recursive subroutine smolyak_product(nmodes, nsets, orders, weights, smolyak_order, imode, ord, smolyak_nprod, smolyak_ord, smolyak_coef)

  integer(ik), intent(in)              :: nmodes
  integer(ik), intent(in)              :: nsets(nmodes)
  integer(ik), intent(in)              :: orders(maxval(nsets),nmodes)
  integer(ik), intent(in)              :: weights(nmodes)
  integer(ik), intent(in)              :: smolyak_order
  integer(ik), intent(in)              :: imode
  integer(ik), intent(inout)           :: ord(nmodes)
  integer(ik), intent(inout)           :: smolyak_nprod
  integer(ik), intent(inout), optional :: smolyak_ord(:,:)
  integer(ik), intent(inout), optional :: smolyak_coef(:)

  integer(ik) :: iset, wsum_min, wsum_max, wsum

  if (imode==1) smolyak_nprod = 0

  wsum_min = max(smolyak_order-nmodes+1, nmodes)
  wsum_max = smolyak_order

  do iset=1, nsets(imode)
    ord(imode) = orders(iset,imode)
    wsum = sum(ord(1:imode)*weights(1:imode))
    if (wsum>wsum_max) cycle
    if (imode==nmodes) then
      if (wsum<wsum_min) cycle
      smolyak_nprod = smolyak_nprod + 1
      if (present(smolyak_ord)) smolyak_ord(:,smolyak_nprod) = ord
      if (present(smolyak_coef)) then
        smolyak_coef(smolyak_nprod) = bico(nmodes-1,smolyak_order-wsum)
        if (mod(smolyak_order-wsum,2)/=0) smolyak_coef(smolyak_nprod) = -smolyak_coef(smolyak_nprod)
      endif
    else
      call smolyak_product(nmodes, nsets, orders, weights, smolyak_order, imode+1, ord, smolyak_nprod, smolyak_ord, smolyak_coef)
    endif
  enddo

end subroutine smolyak_product



! Computes binomial coefficient n choose k
!
function bico(n, k)

  integer(ik), intent(in) :: k, n
  real(rk) :: bico

  bico = nint(exp(factln(n)-factln(k)-factln(n-k)))

end function bico



! Returns ln(n!)
!
function factln(n)

  integer(ik), intent(in) :: n
  real(rk) :: factln

  if (n<0) then
    write(out, '(/a,1x,i4)') 'factln error: negative factorial =', n
    stop
  endif

  factln = gammln(real(n+1,kind=rk))

end function factln



! Returns ln(gamma(x))); call C lgamma function
!
function gammln(x)

  use iso_c_binding
  real(rk), intent(in) :: x
  real(rk) :: gammln

  interface
    real(c_double) function lgamma (y) bind(c)
    use iso_c_binding
    real(c_double), value :: y
    end function lgamma
  end interface

  gammln = lgamma(x)

end function gammln



end module grid
