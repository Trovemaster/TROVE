module rotme_cart_tens

use accuracy
use timer
use fwigxjpf
implicit none


integer(ik), parameter :: maxnelem = 30


type kmat_type
  real(rk), allocatable :: me(:,:,:,:)
end type kmat_type


type rotme_cart_tens_type
  integer(ik)                   ::  nelem = 0
  integer(ik)                   ::  nirrep = 0
  character(cl)                 ::  name = 'N/A'
  character(cl), allocatable    ::  selem(:)
  character(cl), allocatable    ::  sirrep(:)
  complex(rk), allocatable      ::  tmat_s(:,:)
  complex(rk), allocatable      ::  tmat_x(:,:)
  logical                       ::  sym_a = .true.
  integer(ik)                   ::  jmin = 0
  integer(ik)                   ::  jmax = 0
  integer(ik), allocatable      ::  nktau(:)
  integer(ik), allocatable      ::  ktau(:,:,:)
  integer(ik), allocatable      ::  ktau_ind(:,:)
  integer(ik)                   ::  wang_ksign(2) = (/-1,1/)
  type(kmat_type), allocatable  ::  mmat(:,:)
  type(kmat_type), allocatable  ::  kmat(:,:)
  integer(ik)                   ::  dj = 1000
  integer(ik)                   ::  dm = 1000
  integer(ik)                   ::  dk = 1000
  integer(ik)                   ::  kmat_cmplx = 0
  integer(ik), allocatable      ::  mmat_cmplx(:)
  real(rk)                      ::  zero_tol = 1.0d-14
  procedure(prim_me_func), pointer, nopass :: func => null()
  contains
  procedure, public  ::  init => init_rotme_cart_tens_type
  procedure, public  ::  wang_coefs => wang_jktau_coefs
end type rotme_cart_tens_type


abstract interface
  subroutine prim_me_func(q1,q2,name,nelem,nirrep,mf,lf,sirrep,selem)
    use accuracy
    integer(ik), intent(in) :: q1(:), q2(:)
    character(cl), intent(out) :: name
    integer(ik), intent(inout) :: nirrep, nelem
    complex(rk), intent(out), optional :: mf(:,:), lf(:,:)
    character(cl), intent(out), optional :: sirrep(:), selem(:)
  end subroutine prim_me_func
end interface


contains


#include 'rotme_vzz.f90'
#include 'rotme_alpha.f90'
#include 'rotme_beta.f90'
#include 'rotme_mu.f90'
#include 'rotme_costheta.f90'
#include 'rotme_3cos2theta_min1.f90'
#include 'rotme_mf.f90'
#include 'rotme_j.f90'



!###################################################################################################################################



subroutine init_rotme_cart_tens_type(tens, jmin, jmax, dj, verbose)

  class(rotme_cart_tens_type) :: tens
  integer(ik), intent(in) :: jmin, jmax, dj
  logical, optional, intent(in) :: verbose

  integer(ik) :: info, j, ktau, tau_range(2), tau, nimg_j, j1, j2, nktau1, nktau2, nimg_k, ktau1, ktau2, ncoefs1, ncoefs2, icoef, jcoef, nimg_elem, &!
                 irrep, ielem, dk_max, dj_max, k, k1, k2, tau1, tau2, nelem, nirrep, dm_max, m1, m2, m1min, m1max, m2min, m2max, mmin, mmax, &!
                 nimg_elem_lf(maxnelem), nimg_k_lf(maxnelem), nimg_j_lf(maxnelem)
  real(rk) :: zero_tol
  complex(rk) :: coefs1(2), coefs2(2), me(maxnelem,maxnelem), dens, res(maxnelem,maxnelem)
  logical :: verbose_, img_j((jmax+1)**2), img_k((2*jmax+1)**2), img_elem(maxnelem**2), img_elem_lf(maxnelem**2,maxnelem), &!
             img_k_lf((2*jmax+1)**2,maxnelem), img_j_lf((jmax+1)**2,maxnelem)
  character(cl) :: name

  write(out, '(///a)') 'Compute rotational matrix elements of Cartesian tensor operator (rotme_cart_tens/init_rotme_cart_tens_type)'

  verbose_ = .false.
  if (present(verbose)) verbose_ = verbose

  zero_tol = tens%zero_tol
  tens%jmin = jmin
  tens%jmax = jmax

  write(out, '(/1x,a,2(1x,i4)/1x,a,1x,i4)') 'min and max J quanta:', jmin, jmax, "user-defined |j-j'| selection rules:", dj

  call tens%func( (/0,0,0/), (/0,0,0/), name,nelem,nirrep)
  tens%name = name
  tens%nelem = nelem
  tens%nirrep = nirrep

  write(out, '(/1x,a,1x,a/1x,a,1x,i4/1x,a,1x,i4)') 'tensor:', tens%name, 'no. Cartesian components:', tens%nelem, 'no. irreps:', tens%nirrep

  if (allocated(tens%sirrep)) deallocate(tens%sirrep)
  if (allocated(tens%selem)) deallocate(tens%selem)
  allocate(tens%sirrep(tens%nirrep),tens%selem(tens%nelem), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'rotme_cart_tens/init_rotme_cart_tens_type error: failed to allocate tens%sirrep(tens%nirrep),tens%selem(tens%nelem)', &!
    'tens%nelem, tens%nirrep =', tens%nelem, tens%nirrep
    stop 'STOP'
  endif

  call tens%func( (/0,0,0/), (/0,0,0/), name,nelem,nirrep, sirrep=tens%sirrep(1:tens%nirrep), selem=tens%selem(1:tens%nelem))

  write(out, '(1x,a,100(3x,a))') 'Cartesian components:', (trim(tens%selem(ielem)), ielem=1, tens%nelem)
  write(out, '(1x,a,100(3x,a))') 'irreps:', (trim(tens%sirrep(irrep)), irrep=1, tens%nirrep)


  ! initialize rotational quanta K and tau

  if (allocated(tens%nktau)) deallocate(tens%nktau)
  if (allocated(tens%ktau)) deallocate(tens%ktau)
  if (allocated(tens%ktau_ind)) deallocate(tens%ktau_ind)
  allocate( tens%nktau(jmin:jmax), tens%ktau(2,2*jmax+1,jmin:jmax), tens%ktau_ind(0:jmax,0:1), stat=info )
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'rotme_cart_tens/init_rotme_cart_tens_type error: failed to allocate tens%nktau(jmin:jmax), tens%ktau(2,2*jmax+1,jmin:jmax), , tens%ktau_ind(0:jmax,0:1)', &!
    'jmin, jmax =', jmin, jmax
    stop 'STOP'
  endif
  tens%nktau = 0
  tens%ktau_ind = 0

  do j=jmin, jmax
    ktau = 0
    do k=0, j
      tau_range = (/0,1/)
      if (k==0) tau_range = mod(j,2)
      do tau=tau_range(1), tau_range(2)
        ktau = ktau + 1
        tens%ktau(1:2,ktau,j) = (/k,tau/)
        tens%ktau_ind(k,tau) = ktau
      enddo
    enddo
    tens%nktau(j) = ktau
  enddo
  tens%wang_ksign = (/1,-1/)


  ! init module "fwigxjpf" for computing 3j symbols

  write(out, '(/1x,a,1x,i3)') 'initialize "fwigxjpf" for computing 3-j symbols with max angular momentum =', jmax
  call fwig_table_init(int(2*(jmax+1),kind=4), int(3,kind=4))
  call fwig_temp_init(int(2*(jmax+1),kind=4))


  dj_max = 0
  dk_max = 0
  dm_max = 0


  ! compute rotaitonal matrix elements, molecule-fixed K-tensor

  if (allocated(tens%kmat)) deallocate(tens%kmat)
  allocate(tens%kmat(jmin:jmax,jmin:jmax), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'rotme_cart_tens/init_rotme_cart_tens_type error: failed to allocate tens%kmat(jmin:jmax,jmin:jmax)', &!
    'jmin, jmax =', jmin, jmax
    stop 'STOP'
  endif

  if (verbose_) then
    write(out, '(/1x,a/2x,a,2x,a,4x,a,1x,a,4x,a,1x,a,9x,a,8x,a)') &!
    'molecule-fixed frame matrix elements', 'j1', 'j2', 'k1', 'tau1', 'k2', 'tau2', 'elem', 'irrep'
  endif

  nimg_j = 0
  img_j(:) = .false.

  do j1=jmin, jmax
    do j2=jmin, jmax

      if (abs(j1-j2)>dj) cycle

      nktau1 = tens%nktau(j1)
      nktau2 = tens%nktau(j2)

      if (allocated(tens%kmat(j1,j2)%me)) deallocate(tens%kmat(j1,j2)%me)
      allocate(tens%kmat(j1,j2)%me(tens%nelem,tens%nirrep,nktau1,nktau2), stat=info)
      if (info/=0) then
        write(out, '(/a/a,10(1x,i6))') 'rotme_cart_tens/init_rotme_cart_tens_type error: failed to allocate tens%kmat(j1,j2)%me(tens%nelem,tens%nirrep,nktau1,nktau2)', &!
        'j1, j2, tens%nelem, tens%nirrep, nktau1, nktau2 =', j1, j2, tens%nelem, tens%nirrep, nktau1, nktau2
        stop 'STOP'
      endif

      tens%kmat(j1,j2)%me = 0

      nimg_k = 0
      img_k(:) = .false.

      do ktau1=1, nktau1

        k1 = tens%ktau(1,ktau1,j1)
        tau1 = tens%ktau(2,ktau1,j1)
        call tens%wang_coefs(j1, k1, tau1, ncoefs1, coefs1)

        do ktau2=1, nktau2

          k2 = tens%ktau(1,ktau2,j2)
          tau2 = tens%ktau(2,ktau2,j2)
          call tens%wang_coefs(j2, k2, tau2, ncoefs2, coefs2)

          me(:,:) = 0
          do icoef=1, ncoefs1
            do jcoef=1, ncoefs2
              dens = coefs1(icoef) * conjg(coefs2(jcoef))
              call tens%func( (/j1,k1*tens%wang_ksign(icoef),0/), (/j2,k2*tens%wang_ksign(jcoef),0/), name,nelem,nirrep, mf=res )
              me(1:nelem,1:nirrep) = me(1:nelem,1:nirrep) + res(1:nelem,1:nirrep) * dens
            enddo
          enddo

          nimg_elem = 0
          img_elem(:) = .false.

          do irrep=1, tens%nirrep
            do ielem=1, tens%nelem

              if ( abs(real(me(ielem,irrep),rk))>zero_tol .and. abs(aimag(me(ielem,irrep)))<=zero_tol ) then
                nimg_elem = nimg_elem + 1
                img_elem(nimg_elem) = .false.
                tens%kmat(j1,j2)%me(ielem,irrep,ktau1,ktau2) = real(me(ielem,irrep),rk)
              elseif ( abs(real(me(ielem,irrep),rk))<=zero_tol .and. abs(aimag(me(ielem,irrep)))>zero_tol ) then
                nimg_elem = nimg_elem + 1
                img_elem(nimg_elem) = .true.
                tens%kmat(j1,j2)%me(ielem,irrep,ktau1,ktau2) = aimag(me(ielem,irrep))
              elseif ( abs(real(me(ielem,irrep),rk))<=zero_tol .and. abs(aimag(me(ielem,irrep)))<=zero_tol ) then
                cycle
              else
                write(out, '(/a)') 'rotme_cart_tens/init_rotme_cart_tens_type error: element of K-tensor is a complex-valued number'
                write(out, '(1x,a,2(1x,i3),1x,a,2(1x,es16.8),1x,a)') 'ielem/irrep/value:', ielem,irrep,'(', me(ielem,irrep), 'i)'
                stop 'STOP'
              endif

              if (verbose_) then
                if (img_elem(nimg_elem)) then
                  write(out, '(1x,i3,1x,i3,3x,i3,2x,i3,3x,i3,2x,i3,3x,a10,3x,a10,3x,a,1x,es16.8,1x,es16.8,1x,a)') &!
                  j1, j2, k1, tau1, k2, tau2, trim(tens%selem(ielem)), trim(tens%sirrep(irrep)), &!
                  '(', 0.0, tens%kmat(j1,j2)%me(ielem,irrep,ktau1,ktau2), ' i )'
                else
                  write(out, '(1x,i3,1x,i3,3x,i3,2x,i3,3x,i3,2x,i3,3x,a10,3x,a10,3x,a,1x,es16.8,1x,es16.8,1x,a)') &!
                  j1, j2, k1, tau1, k2, tau2, trim(tens%selem(ielem)), trim(tens%sirrep(irrep)), &!
                  '(', tens%kmat(j1,j2)%me(ielem,irrep,ktau1,ktau2), 0.0, ' i )'
                endif
              endif

            enddo
          enddo

          if (nimg_elem>0) then
            nimg_k = nimg_k + 1
            if (all(img_elem(1:nimg_elem))) then
              img_k(nimg_k) = .true.
            elseif (all(.not.img_elem(1:nimg_elem))) then
              img_k(nimg_k) = .false.
            else
              write(out, '(/a/a,6(1x,i3))') 'rotme_cart_tens/init_rotme_cart_tens_type error: K-tensor has mixed real- and imaginary-valued elements', &!
              'j1,k1,tau1/j2,k2,tau2 =', j1,k1,tau1, j2,k2,tau2
              stop 'STOP'
            endif
          endif

          if (any(abs(tens%kmat(j1,j2)%me(:,:,ktau1,ktau2))>zero_tol)) then
            dk_max = max(dk_max,abs(ktau1-ktau2))
            dj_max = max(dj_max,abs(j1-j2))
          endif

        enddo ! ktau2
      enddo ! ktau1

      if (nimg_k>0) then
        nimg_j = nimg_j + 1
        if (all(img_k(1:nimg_k))) then
          img_j(nimg_j) = .true.
        elseif (all(.not.img_k(1:nimg_k))) then
          img_j(nimg_j) = .false.
        else
          write(out, '(/a/a,2(1x,i3))') &!
          'rotme_cart_tens/init_rotme_cart_tens_type error: K-tensor has mixed real- and imaginary-valued elements for different pairs of k,tau-quanta', &!
          ' j1/j2 =', j1, j2
          stop 'STOP'
        endif
      endif

    enddo ! j2
  enddo ! j1

  if (all(img_j(1:nimg_j))) then
    tens%kmat_cmplx = -1 ! imaginary numbers
  elseif (all(.not.img_j(1:nimg_j))) then
    tens%kmat_cmplx = 0 ! real numbers
  else
    write(out, '(/a)') &!
    'rotme_cart_tens/init_rotme_cart_tens_type error: K-tensor has mixed real- and imaginary-valued elements for different pairs of J-quanta'
    stop 'STOP'
  endif


  ! compute rotaitonal matrix elements, laboratory-fixed M-tensor

  mmin = -jmax
  mmax = jmax

  if (allocated(tens%mmat)) deallocate(tens%mmat)
  if (allocated(tens%mmat_cmplx)) deallocate(tens%mmat_cmplx)
  allocate(tens%mmat(jmin:jmax,jmin:jmax), tens%mmat_cmplx(tens%nelem), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') 'rotme_cart_tens/init_rotme_cart_tens_type error: failed to allocate tens%mmat(jmin:jmax,jmin:jmax), tens%mmat_cmplx(tens%nelem)', &!
    'jmin, jmax, tens%nelem =', jmin, jmax, tens%nelem
    stop 'STOP'
  endif

  if (verbose_) then
    write(out, '(/1x,a/2x,a,2x,a,4x,a,2x,a,8x,a,9x,a)') &!
    'laboratory-fixed frame matrix elements', 'j1', 'j2', 'm1', 'm2', 'irrep', 'elem'
  endif

  nimg_j_lf(:) = 0
  img_j_lf(:,:) = .false.

  do j1=jmin, jmax
    do j2=jmin, jmax

      if (abs(j1-j2)>dj) cycle

      m1min = max(mmin,-j1)
      m1max = min(mmax,j1)
      m2min = max(mmin,-j2)
      m2max = min(mmax,j2)

      if (allocated(tens%mmat(j1,j2)%me)) deallocate(tens%mmat(j1,j2)%me)
      allocate(tens%mmat(j1,j2)%me(tens%nirrep,tens%nelem,m1min:m1max,m2min:m2max), stat=info)
      if (info/=0) then
        write(out, '(/a/a,10(1x,i6))') 'rotme_cart_tens/init_rotme_cart_tens_type error: failed to allocate tens%mmat(j1,j2)%me(tens%nirrep,tens%nelem,m1min:m1max,m2min:m2max)', &!
        'j1, j2, tens%nirrep, tens%nelem, m1min, m1max, m2min, m2max =', j1, j2, tens%nirrep, tens%nelem, m1min, m1max, m2min, m2max
        stop 'STOP'
      endif

      tens%mmat(j1,j2)%me = 0

      nimg_k_lf(:) = 0
      img_k_lf(:,:) = .false.

      do m1=m1min, m1max
        do m2=m2min, m2max

          call tens%func( (/j1,0,m1/), (/j2,0,m2/), name,nelem,nirrep, lf=res)
          me(1:nelem,1:nirrep) = res(1:nelem,1:nirrep)


          nimg_elem_lf(:) = 0
          img_elem_lf(:,:) = .false.

          do ielem=1, tens%nelem
            do irrep=1, tens%nirrep

              if ( abs(real(me(ielem,irrep),rk))>zero_tol .and. abs(aimag(me(ielem,irrep)))<=zero_tol ) then
                nimg_elem_lf(ielem) = nimg_elem_lf(ielem) + 1
                img_elem_lf(nimg_elem_lf(ielem),ielem) = .false.
                tens%mmat(j1,j2)%me(irrep,ielem,m1,m2) = real(me(ielem,irrep),rk)
              elseif ( abs(real(me(ielem,irrep),rk))<=zero_tol .and. abs(aimag(me(ielem,irrep)))>zero_tol ) then
                nimg_elem_lf(ielem) = nimg_elem_lf(ielem) + 1
                img_elem_lf(nimg_elem_lf(ielem),ielem) = .true.
                tens%mmat(j1,j2)%me(irrep,ielem,m1,m2) = aimag(me(ielem,irrep))
              elseif ( abs(real(me(ielem,irrep),rk))<=zero_tol .and. abs(aimag(me(ielem,irrep)))<=zero_tol ) then
                cycle
              else
                write(out, '(/a)') 'rotme_cart_tens/init_rotme_cart_tens_type error: element of M-tensor is a complex-valued number'
                write(out, '(1x,a,2(1x,i3),1x,a,2(1x,es16.8),1x,a)') 'ielem/irrep/value:', ielem,irrep,'(', me(ielem,irrep), 'i)'
                stop 'STOP'
              endif

              if (verbose_) then
                if (img_elem_lf(nimg_elem_lf(ielem),ielem)) then
                  write(out, '(1x,i3,1x,i3,3x,i3,1x,i3,3x,a10,3x,a10,3x,a,1x,es16.8,1x,es16.8,1x,a)') &!
                  j1, j2, m1, m2, trim(tens%sirrep(irrep)), trim(tens%selem(ielem)), &!
                  '(', 0.0, tens%mmat(j1,j2)%me(irrep,ielem,m1,m2), ' i )'
                else
                  write(out, '(1x,i3,1x,i3,3x,i3,1x,i3,3x,a10,3x,a10,3x,a,1x,es16.8,1x,es16.8,1x,a)') &!
                  j1, j2, m1, m2, trim(tens%sirrep(irrep)), trim(tens%selem(ielem)), &!
                  '(', tens%mmat(j1,j2)%me(irrep,ielem,m1,m2), 0.0, ' i )'
                endif
              endif

            enddo ! irrep

            if (nimg_elem_lf(ielem)>0) then
              nimg_k_lf(ielem) = nimg_k_lf(ielem) + 1
              if (all(img_elem_lf(1:nimg_elem_lf(ielem),ielem))) then
                img_k_lf(nimg_k_lf(ielem),ielem) = .true.
              elseif (all(.not.img_elem_lf(1:nimg_elem_lf(ielem),ielem))) then
                img_k_lf(nimg_k_lf(ielem),ielem) = .false.
              else
                write(out, '(/a/a,4(1x,i3))') 'rotme_cart_tens/init_rotme_cart_tens_type error: M-tensor has mixed real- and imaginary-valued elements', &!
                'j1,m1/j2,m2 =', j1,m1, j2,m2
                stop 'STOP'
              endif
            endif

          enddo ! ielem

          if (any(abs(tens%mmat(j1,j2)%me(:,:,m1,m2))>zero_tol)) then
            dm_max = max(dm_max,abs(m1-m2))
            dj_max = max(dj_max,abs(j1-j2))
          endif

        enddo ! m2
      enddo ! m1

      do ielem=1, nelem
        if (nimg_k_lf(ielem)>0) then
          nimg_j_lf(ielem) = nimg_j_lf(ielem) + 1
          if (all(img_k_lf(1:nimg_k_lf(ielem),ielem))) then
            img_j_lf(nimg_j_lf(ielem),ielem) = .true.
          elseif (all(.not.img_k_lf(1:nimg_k_lf(ielem),ielem))) then
            img_j_lf(nimg_j_lf(ielem),ielem) = .false.
          else
            write(out, '(/a/a,2(1x,i3))') &!
            'rotme_cart_tens/init_rotme_cart_tens_type error: M-tensor has mixed real- and imaginary-valued elements for different pairs of m-quanta', &!
            'j1/j2 =', j1, j2
            stop 'STOP'
          endif
        endif
      enddo

    enddo ! j2
  enddo ! j1

  do ielem=1, nelem
    if (all(img_j_lf(1:nimg_j_lf(ielem),ielem))) then
      tens%mmat_cmplx(ielem) = -1
    elseif (all(.not.img_j_lf(1:nimg_j_lf(ielem),ielem))) then
      tens%mmat_cmplx(ielem) = 0
    else
      write(out, '(/a)') &!
      'rotme_cart_tens/init_rotme_cart_tens_type error: M-tensor has mixed real- and imaginary-valued elements for different pairs of J-quanta'
      stop 'STOP'
    endif
  enddo


  tens%dj = dj_max
  tens%dm = dm_max
  tens%dk = dk_max

  write(out, '(/1x,a,1x,i4)') "|ktau-ktau'| selection rules:", tens%dk
  write(out, '(1x,a,1x,i4)') "|m-m'| selection rules:", tens%dm
  write(out, '(1x,a,1x,i4)') "|j-j'| selection rules:", tens%dj

  write(out,'(a)') 'done (rotme_cart_tens/init_rotme_cart_tens_type)'

end subroutine init_rotme_cart_tens_type



!###################################################################################################################################



subroutine wang_jktau_coefs(tens, j, k, tau, ncoefs, coefs)

  class(rotme_cart_tens_type), intent(in) :: tens
  integer(ik), intent(in) :: j, k, tau
  integer(ik), intent(out) :: ncoefs
  complex(rk), intent(out), optional :: coefs(2)

  integer(ik) :: sigma
  real(rk) :: fac1, fac2

  sigma = mod(k,3)*tau
  fac1 = real((-1)**sigma,rk)/sqrt(2.0_rk)
  fac2 = fac1 * (-1)**(j+k)

  if (tau==0) then

    if (k==0) then
      ncoefs = 1
      if (present(coefs)) coefs = (/ cmplx(1.0_rk,0.0_rk,kind=rk), cmplx(0.0_rk,0.0_rk,kind=rk) /)
    elseif (k>0) then
      ncoefs = 2
      if (present(coefs)) coefs = (/ cmplx(fac1,0.0_rk,kind=rk), cmplx(fac2,0.0_rk,kind=rk) /)
    else
      write(out, '(/a,1x,i3,1x,a)') 'rotme_cart_tens/wang_jktau_coefs error: rotational quantum number "k" =', k, 'is less than zero (expected >=0)'
      stop 'STOP'
    endif

  elseif (tau==1) then

    if (k==0) then
      ncoefs = 1
      if (present(coefs)) coefs = (/ cmplx(0.0_rk,1.0_rk,kind=rk), cmplx(0.0_rk,0.0_rk,kind=rk) /)
    elseif (k>0) then
      ncoefs = 2
      if (present(coefs)) coefs = (/ cmplx(0.0_rk,fac1,kind=rk), cmplx(0.0_rk,-fac2,kind=rk) /)
    else
      write(out, '(/a,1x,i3,1x,a)') 'rotme_cart_tens/wang_jktau_coefs error: rotational quantum number "k" =', k, 'is less than zero (expected >=0)'
      stop 'STOP'
    endif

  else

    write(out, '(/a,1x,i3,1x,a)') 'rotme_cart_tens/wang_jktau_coefs error: rotational quantum number "tau" =', tau, '(expected 0 or 1)'
    stop 'STOP'

  endif

end subroutine wang_jktau_coefs



!###################################################################################################################################



subroutine pseudoinverse(nrows, ncols, mat, invmat, tol_)

  integer(ik), intent(in) :: nrows, ncols
  complex(rk), intent(in) :: mat(nrows,ncols)
  complex(rk), intent(out) :: invmat(ncols,nrows)
  real(rk), intent(in), optional :: tol_

  integer(ik) lwork, info, i
  complex(rk) :: work1(1)
  complex(rk), allocatable :: work(:), matu(:,:), matvt(:,:), mat_d(:,:), matt(:,:)
  double precision :: tol, rwork(5*max(nrows,ncols))
  double precision, allocatable :: sv(:)

  if (present(tol_)) then
    tol = tol_
  else
    tol = 1.0d-14 !epsilon(1.0_rk)
  endif

  allocate(matt(nrows,ncols), sv(min(ncols,nrows)), matu(nrows,nrows), matvt(ncols,ncols), mat_d(ncols,nrows), stat=info)
  if (info/=0) then
    write(out, '(/a/a,10(1x,i6))') &!
    'rotbas/pseudoinverse error: failed to allocate matt(nrows,ncols), sv(min(ncols,nrows)), matu(nrows,nrows), matvt(ncols,ncols), mat_d(ncols,nrows)', &!
    'ncols, nrows =', ncols, nrows
    stop 'STOP'
  endif
  matt = 0
  sv = 0
  matu = 0
  matvt = 0
  mat_d = 0

  matt = mat

  ! estimate size of the workspace array (by setting lwork<0 to ZGESVD)
  lwork = -1
  call zgesvd('A', 'A', nrows, ncols, matt, nrows, sv, matu, nrows, matvt, ncols, work1, lwork, rwork, info)
  lwork = int(work1(1), kind=ik)

  ! allocate workspace array
  allocate(work(lwork), stat=info)
  if (info/=0) then
    write(out, '(/a,1x,i6)') 'rotbas/pseudoinverse error: failed to allocate workspace array required for SVD, size =', lwork
    stop 'STOP'
  endif

  ! perform actual SVD
  call zgesvd('A', 'A', nrows, ncols, matt, nrows, sv, matu, nrows, matvt, ncols, work, lwork, rwork, info)

  if (info/=0) then
    write(out, '(/a,1x,i6)') 'rotbas/pseudoinverse error: SVD failed, info =', info
    stop 'STOP'
  endif

  ! inverse diagonal matrix
  mat_d = 0.0
  do i=1, min(nrows,ncols)
    if (sv(i)>=tol) then
      mat_d(i,i) = 1.0_rk/sv(i)
    !else
    !  write(out, '(/a)') 'rotbas/pseudoinverse error: matrix is singular'
    !  write(out, '(a)') 'singular elements:'
    !  do j=1, min(nrows,ncols)
    !    write(out, '(1x,i3,1x,f)') j, sv(j)
    !  enddo
    !  stop 'STOP'
    endif
  enddo

  ! compute pseudoinverse
  mat_d = matmul(mat_d, conjg(transpose(matu)))
  invmat(1:ncols,1:nrows) = matmul(conjg(transpose(matvt)), mat_d)

  ! clear workspace
  deallocate(work,matt,mat_d,sv,matu,matvt)

end subroutine pseudoinverse



!###################################################################################################################################



function threej_symbol(j1,j2,j3, m1,m2,m3) result(f)

  integer(ik), intent(in) :: j1,j2,j3,m1,m2,m3
  real(rk) :: f

  integer(4) :: two, two_j1, two_j2, two_j3, two_m1, two_m2, two_m3

  two = 2_ik

  two_j1 = two*j1
  two_j2 = two*j2
  two_j3 = two*j3

  two_m1 = two*m1
  two_m2 = two*m2
  two_m3 = two*m3

  f = fwig3jj(two_j1, two_j2, two_j3, two_m1, two_m2, two_m3)

end function threej_symbol


end module rotme_cart_tens
