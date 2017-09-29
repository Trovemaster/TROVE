module diag
 
!dec$ define arpack_ = 0
!dec$ define blacs_  = 0
!dec$ define mpi_    = 0
!dec$ define omparpack_ = 0
!dec$ define propack_ = 0

!
!  Simplistic type-agnostic LAPACK interface
!
  use accuracy
  use timer
  implicit none
  private verbose
  !
  !
  public diag_tridiag,diag_tridiag_pack,diag_dsyev_i8,diag_dgelss,diag_syev_ilp,diag_syev_i8,diag_dseupd,diag_dseupd_p
  public dseupd_omp_arpack,diag_propack,daprod
  !
  integer,parameter:: verbose = 4
  integer(hik),allocatable,save   :: hikparm(:)
  real(rk),allocatable,save       :: hmat_(:,:)
  !
  
  contains

  subroutine diag_tridiag(h,e,solver,jobz,rng,iroots,vrange,irange,tol)
    double precision, intent(inout) :: h(:,:)
    double precision, intent(out)   :: e(:)
    character(len=cl),intent(in)    :: solver
    character(len=1),intent(in),optional   :: rng,jobz
    double precision, intent(in),optional  :: vrange(2)  ! In:  energy range to be computed
    integer(ik)     , intent(in),optional  :: irange(2)  ! In:  index to be computed
    integer(ik)     , intent(out),optional :: iroots     ! Out: Number of roots found
    double precision, intent(in),optional  :: tol    ! In: abs. tolerance 
    !
    double precision,allocatable    :: work(:),offd(:),d(:),tau(:)
    double precision,allocatable :: a(:,:)
    integer,allocatable :: iwork(:),ifail(:)
    integer,allocatable :: isuppz(:)
    !
    integer          :: n,lwork,niwork,m
    integer          :: nh1,il, ir, nb,ilaenv,iinfo
    !
    integer(ik)      :: msize,alloc,k,iksize,info
    !
    double precision :: vl,vu,abstol
    character(len=1) :: rng_,jobz_
    integer(hik)     ::  matsize
    !
    if (verbose>=2) call TimerStart('diag_tridiag')
    !
    if (verbose>=3) write(out,"('Converting to a tridiagonal form ...')") 
    !
    n = size(h,dim=1)
    !
    jobz_ = 'V'
    !
    if (present(jobz)) then
       jobz_ = jobz
    endif
    !
    nb = ilaenv( 1, 'DSYTRD','U',n, -1, -1, -1 )
    !
    if (verbose>=3) write(out,"('The optional block size  = ',i0)") nb
    !
    lwork = (nb+10)*n
    !
    allocate(work(lwork),d(n),offd(n-1),tau(n-1),stat=alloc)
    !
    iksize = lwork+3*n-2
    !
    call ArrayStart('diag_tridiag-arrays',alloc,iksize,rk)
    !
    if (verbose>=2) call TimerStart('diag_dsytrd')
    !
    call dsytrd('U',n,h,n,d,offd,tau,work,-1,info)
    !
    if (info/=0) then
      write (out,"(' dsytrd-1 returned ',i0)") info
      stop 'dsytrd  failed'
    end if
    !
    if (int(work(1))>size(work)) then 
      !
      lwork = int(work(1))
      !
      if (verbose>=3) write(out,"('work-size = ',i0)") lwork
      !
      deallocate(work)
      !
      !ArrayMinus('diag_tridiag-arrays')
      !
      allocate(work(lwork),stat=alloc)
      !
      call ArrayStart('diag_tridiag-arrays',alloc,int(lwork,ik),rk)
      !
    endif 
    !
    !call dsytrd('U',n,h,n,d,offd,tau,work,lwork,info)
    !
    call diag_dsytrd('U',n,h,d,offd,tau,info)
    !
    if (info/=0) then
      write (out,"(' dsytrd returned ',i0)") info
      stop 'dsytrd  failed'
    end if
    !
    if (verbose>=2) call TimerStop('diag_dsytrd')
    !
    if (verbose>=3) write(out,"('...done!')") 
    !
    if (verbose>=3) write(out,"('Diagonalizing the tridiagonal matrix...')") 
    !
    if (verbose>=2) call TimerStart('diag_dstev')
    !
    rng_ = 'A'
    il = 1   ; ir = 100000
    vl = -.001 ; vu = 100000.0
    !
    if (present(rng)) then
       rng_ = rng
    endif
    !
    if (present(vrange).and.rng_=='V') then
       vl = vrange(1) ; vu = vrange(2)
    endif
    !
    if (present(irange).and.rng_=='I') then
       il = irange(1) ; ir = irange(2)
       if (irange(2)>=size(h,dim=1)) rng_ = 'A' 
    endif
    !
    abstol =(small_)
    if (present(tol)) then
       abstol = tol
    endif
    !
    msize = n
    if (rng=='I'.and.trim(solver)/='DSTEV'.and.trim(solver)/='DSTEVD')  msize = ir-il+1
    !
    m = msize
    !
    allocate(ifail(n),stat=alloc)
    !
    iksize = size(ifail)
    call ArrayStart('diag_tridiag-arrays',alloc,iksize,ik)
    !
    allocate (a(n,msize),stat=alloc)
    matsize = int(n,hik)*int(msize,hik)
    call ArrayStart('diag_tridiag-arrays',alloc,1,rk,matsize )
    !
    !
    if (verbose>=3) call MemoryReport
    !
    if (verbose>=4) write(out,"(/'abstol = ',e18.8)") abstol
    !
    select case (trim(solver))
      !
    case('DSTEV')
      !
      lwork = max(1,2*N-2)
      !
      deallocate(work)
      !
      allocate(work(lwork),stat=alloc)
      !
      iksize = size(work)
      call ArrayStart('diag_tridiag-arrays',alloc,iksize,rk)
      !
      call dstev(jobz_,n,d,offd,a,n,work,info)
      !
      e = d
      !
    case('DSTEVX')
      !
      niwork = 5*n
      lwork = 5*n
      deallocate(work)
      !
      allocate(work(lwork),stat=alloc)
      call ArrayStart('diag_tridiag-arrays',alloc,int(lwork,ik),rk)
      !
      allocate(iwork(niwork),stat=alloc)
      call ArrayStart('diag_tridiag-arrays',alloc,int(niwork,ik),ik)
      !
      call dstevx(jobz_,rng,n,d,offd,vl,vu,il,ir,abstol,m,e,a,n,work,iwork,ifail,info)
      !
    case('DSTEVR')
      !
      niwork = 20*n
      !
      allocate(iwork(niwork),stat=alloc)
      call ArrayStart('diag_tridiag-arrays',alloc,int(niwork,ik),ik)
      !
      allocate(isuppz(2*m),stat=alloc)
      iksize = 2*m
      call ArrayStart('diag_tridiag-arrays',alloc,iksize,ik)
      !
      call dstevr(jobz_,rng,n,d,offd,vl,vu,il,ir,abstol,m,e,a,n,isuppz,work,-1,iwork,-1,info)
      !
      if (int(work(1))>size(work)) then 
        !
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork),stat=alloc)
        !
      endif
      !
      if (int(iwork(1))>size(work)) then 
        !
        niwork = int(iwork(1))
        deallocate(iwork)
        allocate(iwork(niwork),stat=alloc)
        !
      endif
      !
      call dstevr(jobz_,rng,n,d,offd,vl,vu,il,ir,abstol,m,e,a,n,isuppz,work,lwork,iwork,niwork,info)
      !
    case('DSTEVD')
      !
      niwork = 5*n+3
      !
      allocate(iwork(niwork),stat=alloc)
      call ArrayStart('diag_tridiag-arrays',alloc,int(niwork,ik),ik)
      !
      call dstevd(jobz_,n,d,offd,a,n,work,-1,iwork,-1,info)
      !
      !
      if (int(work(1))>size(work)) then 
        !
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork),stat=alloc)
        !
      endif
      !
      if (int(iwork(1))>size(work)) then 
        !
        niwork = int(iwork(1))
        deallocate(iwork)
        allocate(iwork(niwork),stat=alloc)
        !
      endif
      !
      call dstevd(jobz_,n,d,offd,a,n,work,lwork,iwork,niwork,info)
      !
      e = d
      !
    end select
    !
    !
    if (verbose>=2) call TimerStop('diag_dstev')
    !
    if (verbose>=3) write(out,"('...done!')") 
    !
    if (info/=0) then
      write (out,"(' dstevx returned ',i8)") info
      stop 'diag - dstevxfailed'
    end if
    !
    if (verbose>=2) call TimerStart('diag_dormtr')
    !
    if (verbose>=3) write(out,"('Unitary transformation to the original representation...')") 
    !
    call dormtr('L','U','N',n,m,h,n,tau,a,n,work,-1,info )
    !
    if (int(work(1))>size(work)) then 
      !
      lwork = int(work(1))
      !
      deallocate(work)
      !
      allocate(work(lwork),stat=alloc)
      !
      call ArrayStart('diag_tridiag-arrays',alloc,int(lwork,ik),rk)
      !
    endif 
    !
    call dormtr('L','U','N',n,m,h,n,tau,a(1:n,1:msize),n,work,lwork,info )
    !
    if (verbose>=3) write(out,"('...done!')") 
    !
    if (verbose>=2) call TimerStop('diag_dormtr')
    !
    if (verbose>=2) call TimerStart('diag_copy')
    !
    !$omp parallel do private(k) shared(h) schedule(dynamic)
    do k=1,m
       !h(:,k) = a(:,k)
       call dcopy(n,a(:,k),1,h(:,k),1)
    enddo
    !$omp end parallel do
    !
    if (verbose>=2) call TimerStop('diag_copy')
    !
    if (allocated(isuppz)) deallocate(isuppz)
    if (allocated(ifail)) deallocate(ifail)
    if (allocated(iwork)) deallocate(iwork)
    !
    deallocate(a)
    !
    if (present(iroots))  iroots = m
    !
    deallocate(work,offd,tau,d)
    !
    call ArrayStop('diag_tridiag-arrays')
    !
    if (verbose>=2) call TimerStop('diag_tridiag')
    !
  end subroutine diag_tridiag


  subroutine diag_tridiag_pack(h,e,solver,jobz,rng,iroots,vrange,irange,tol)
    double precision, intent(inout) :: h(:,:)
    double precision, intent(out)   :: e(:)
    character(len=cl),intent(in)    :: solver
    character(len=1),intent(in),optional   :: rng,jobz
    double precision, intent(in),optional  :: vrange(2)  ! In:  energy range to be computed
    integer(ik)     , intent(in),optional  :: irange(2)  ! In:  index to be computed
    integer(ik)     , intent(out),optional :: iroots     ! Out: Number of roots found
    double precision, intent(in),optional  :: tol    ! In: abs. tolerance 
    !
    double precision,allocatable    :: work(:),offd(:),d(:),tau(:)
    double precision,allocatable :: ap(:)
    integer,allocatable :: iwork(:),ifail(:)
    integer,allocatable :: isuppz(:)
    !
    integer(ik)      :: alloc,i,j,iksize
    !
    integer           :: info
    integer           :: n,ilaenv
    integer(ik)       :: lwork,niwork
    !
    double precision  :: vl,vu,abstol
    integer           :: nh1,il, ir,m,nb
    character(len=1)  :: rng_,jobz_
    integer           :: msize
    integer(hik)      :: matsize,k
    !
    if (verbose>=2) call TimerStart('diag_tridiag')
    !
    n = size(h,dim=1)
    !
    jobz_ = 'V'
    !
    if (present(jobz)) then
       jobz_ = jobz
    endif
    !    
    allocate(d(n),offd(n-1),tau(n-1),stat=alloc)
    !
    iksize = n
    call ArrayStart('diag_tridiag-arrays',alloc,iksize,rk)
    call ArrayStart('diag_tridiag-arrays',alloc,iksize,rk)
    call ArrayStart('diag_tridiag-arrays',alloc,iksize,rk)
    !
    matsize = n*(n+1)/2_hik
    !
    allocate(ap(matsize),stat=alloc)
    !
    call ArrayStart('diag_tridiag-arrays',alloc,1,rk,matsize)
    !
    if (verbose>=2) call TimerStart('diag_copy_compact')
    !
    !$omp parallel do private(i,j,k) shared(ap) schedule(dynamic)
    do j=1,n
      do i=1,j
        !
        k = int(i,hik)+int(j,hik)*int(j-1,hik)/2_hik
        !
        ap(k) = h(i,j)
        !
      enddo
    enddo
    !$omp end parallel do
    !
    if (verbose>=2) call TimerStop('diag_copy_compact')
    !
    if (verbose>=3) write(out,"('Converting to a tridiagonal form ...')") 
    !
    if (verbose>=2) call TimerStart('diag_dsptrd')
    !
    nb = ilaenv( 1, 'DSYTRD','U',n, -1, -1, -1 )
    !
    if (verbose>=3) write(out,"('The optional block size  = ',i0)") nb
    !
    !call dsptrd('U',n,ap,d,offd,tau,info)
    !
    call diag_dsptrd('U',n,ap,d,offd,tau,info)
    !
    !call diag_dsptrd_omp('U',n,ap,d,offd,tau,info)
    !
    if (info/=0) then
      write (out,"(' dsptrd returned ',i0)") info
      stop 'dsytrd  failed'
    end if
    !
    if (verbose>=2) call TimerStop('diag_dsptrd')
    !
    if (verbose>=3) write(out,"('...done!')") 
    !
    if (verbose>=3) write(out,"('Diagonalizing the tridiagonal matrix...')") 
    !
    if (verbose>=2) call TimerStart('diag_dstev')
    !
    rng_ = 'A'
    il = 1   ; ir = 100000
    vl = -.001 ; vu = 100000.0
    !
    if (present(rng)) then
       rng_ = rng
    endif
    !
    if (present(vrange).and.rng_=='V') then
       vl = vrange(1) ; vu = vrange(2)
    endif
    !
    if (present(irange).and.rng_=='I') then
       il = irange(1) ; ir = irange(2)
       if (irange(2)>=size(h,dim=1)) rng_ = 'A' 
    endif
    !
    abstol =(small_)
    if (present(tol)) then
       abstol = tol
    endif
    !
    msize = n
    if (rng=='I'.and.trim(solver)/='DSTEV-P'.and.trim(solver)/='DSTEVD-P')  msize = ir-il+1
    !
    m = msize
    !
    allocate(ifail(n),stat=alloc)
    call ArrayStart('diag_tridiag-arrays',alloc,int(n,ik),ik)
    !
    if (verbose>=3) call MemoryReport
    !
    if (verbose>=4) write(out,"(/'abstol = ',e18.8)") abstol
    !
    select case (trim(solver))
      !
    case('DSTEV-P')
      !
      lwork = max(1,2*N-2)
      !
      allocate(work(lwork),stat=alloc)
      call ArrayStart('diag_tridiag-arrays',alloc,lwork,rk)
      !
      call dstev(jobz_,n,d,offd,h(1:n,1:msize),n,work,info)
      !
      e = d
      !
    case('DSTEVX-P')
      !
      niwork = 5*n
      lwork = 5*n
      !
      allocate(work(lwork),stat=alloc)
      call ArrayStart('diag_tridiag-arrays',alloc,lwork,rk)
      !
      allocate(iwork(niwork),stat=alloc)
      call ArrayStart('diag_tridiag-arrays',alloc,niwork,ik)
      !
      call dstevx(jobz_,rng,n,d,offd,vl,vu,il,ir,abstol,m,e,h(1:n,1:msize),n,work,iwork,ifail,info)
      !
    case('DSTEVR-P')
      !
      lwork = 30*n
      niwork = 30*n
      !
      allocate(work(lwork),stat=alloc)
      call ArrayStart('diag_tridiag-arrays',alloc,lwork,rk)
      !
      allocate(iwork(niwork),stat=alloc)
      call ArrayStart('diag_tridiag-arrays',alloc,niwork,ik)
      !
      allocate(isuppz(2*m),stat=alloc)
      call ArrayStart('diag_tridiag-arrays',alloc,int(2*m,ik),ik)
      !
      call dstevr(jobz_,rng,n,d,offd,vl,vu,il,ir,abstol,m,e,h(1:n,1:m),n,isuppz,work,-1,iwork,-1,info)
      !
      if (int(work(1))>size(work)) then 
        !
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork),stat=alloc)
        !
      endif
      !
      if (int(iwork(1))>size(work)) then 
        !
        niwork = int(iwork(1))
        deallocate(iwork)
        allocate(iwork(niwork),stat=alloc)
        !
      endif
      !
      call dstevr(jobz_,rng,n,d,offd,vl,vu,il,ir,abstol,m,e,h(1:n,1:m),n,isuppz,work,lwork,iwork,niwork,info)
      !
    case('DSTEVD-P')
      !
      niwork = 5*n+3
      lwork = 5 + 5*n+n**2
      !
      allocate(work(lwork),stat=alloc)
      call ArrayStart('diag_tridiag-arrays',alloc,lwork,rk)
      !
      allocate(iwork(niwork),stat=alloc)
      call ArrayStart('diag_tridiag-arrays',alloc,niwork,ik)
      !
      call dstevd(jobz_,n,d,offd,h(1:n,1:n),n,work,-1,iwork,-1,info)
      !
      if (int(work(1))>size(work)) then 
        !
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork),stat=alloc)
        !
      endif
      !
      if (int(iwork(1))>size(work)) then 
        !
        niwork = int(iwork(1))
        deallocate(iwork)
        allocate(iwork(niwork),stat=alloc)
        !
      endif
      !
      call dstevd(jobz_,n,d,offd,h(1:n,1:n),n,work,lwork,iwork,niwork,info)
      !
      e = d
      !
    end select
    !
    !
    if (verbose>=2) call TimerStop('diag_dstev')
    !
    if (verbose>=3) write(out,"('...done!')") 
    !
    if (info/=0) then
      write (out,"(' dstev returned ',i0)") info
      stop 'diag - dstev failed'
    end if
    !
    if (verbose>=2) call TimerStart('diag_dopmtr')
    !
    if (verbose>=3) write(out,"('Unitary transformation to the original representation...')") 
    !
    lwork = n
    !
    deallocate(work)
    ! 
    allocate(work(lwork),stat=alloc)
    !
    call ArrayStart('diag_tridiag-arrays',alloc,lwork,rk)
    !
    call dopmtr('L','U','N',n,m,ap,tau,h(1:n,1:msize),n,work,info)
    !
    if (info/=0) then
      write (out,"(' dopmtr returned ',i0)") info
      stop 'diag - dopmtr failed'
    end if
    !
    if (verbose>=3) write(out,"('...done!')") 
    !
    if (verbose>=2) call TimerStop('diag_dopmtr')
    !
    if (allocated(isuppz)) deallocate(isuppz)
    if (allocated(ifail)) deallocate(ifail)
    if (allocated(iwork)) deallocate(iwork)
    !
    deallocate(ap)
    !
    if (present(iroots))  iroots = m
    !
    deallocate(work,offd,tau,d)
    !
    call ArrayStop('diag_tridiag-arrays')
    !
    if (verbose>=2) call TimerStop('diag_tridiag')
    !
  end subroutine diag_tridiag_pack




  subroutine diag_dsyev_i8(h,e,jobz)
    double precision, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    double precision, intent(out)   :: e(:)    ! Out: Eigenvalues
    character(len=1),intent(in),optional :: jobz
    character(len=1) :: jobz_

    double precision,allocatable    :: work(:)
    integer(hik)      :: info
    integer(hik)      :: nh1, nh2
    integer(hik)      :: lwork
    !
    !if (verbose>=2) call TimerStart('lapack_dsyev: diagonalization')
    !
    jobz_ = 'V'
    !
    if (present(jobz)) then
       jobz_ = jobz
    endif
    !
    print("('kind(hik) = ',i0)"), kind(info)
    !
    lwork = 50*size(h,dim=1)
    !
    allocate(work(lwork))
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call dsyev(jobz_,'U',nh1,h(1:nh1,1:nh2),nh1,e,work,-1,info)
    !
    if (int(work(1))>size(work)) then
      !
      lwork = int(work(1))
      !
      if (verbose>=3) write(out,"('lwork = ',i0)")  lwork
      !
      deallocate(work)
      !
      allocate(work(lwork))
      !
      call dsyev(jobz_,'U',nh1,h(1:nh1,1:nh2),nh1,e,work,lwork,info)
      !
    else
      !
      call dsyev(jobz_,'U',nh1,h(1:nh1,1:nh2),nh1,e,work,lwork,info)
      !
    endif
    !
    deallocate(work)
    !
    if (info/=0) then
      write (out,"(' dsyev returned ',i8)") info
      stop 'lapack_dsyev - dsyev failed'
    end if
    !
    !if (verbose>=2) call TimerStop('lapack_dsyev: diagonalization')
    !
  end subroutine diag_dsyev_i8


  subroutine diag_syev_i8(h,e,jobz)
    real(4), intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    real(4), intent(out)   :: e(:)    ! Out: Eigenvalues
    character(len=1),intent(in),optional :: jobz
    character(len=1) :: jobz_

    real(4),allocatable    :: work(:)
    integer(hik)      :: info
    integer(hik)      :: nh1, nh2
    integer(hik)      :: lwork
    !
    !if (verbose>=2) call TimerStart('lapack_dsyev: diagonalization')
    !
    jobz_ = 'V'
    !
    if (present(jobz)) then
       jobz_ = jobz
    endif
    !
    print("('kind(hik) = ',i0)"), kind(info)
    !
    lwork = 50*size(h,dim=1)
    !
    allocate(work(lwork))
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call ssyev(jobz_,'U',nh1,h(1:nh1,1:nh2),nh1,e,work,-1,info)
    !
    if (int(work(1))>size(work)) then
      !
      lwork = int(work(1))
      !
      if (verbose>=3) write(out,"('lwork = ',i0)")  lwork
      !
      deallocate(work)
      !
      allocate(work(lwork))
      !
      call ssyev(jobz_,'U',nh1,h(1:nh1,1:nh2),nh1,e,work,lwork,info)
      !
    else
      !
      call ssyev(jobz_,'U',nh1,h(1:nh1,1:nh2),nh1,e,work,lwork,info)
      !
    endif
    !
    deallocate(work)
    !
    if (info/=0) then
      write (out,"(' ssyev returned ',i8)") info
      stop 'lapack_ssyev - dsyev failed'
    end if
    !
    !if (verbose>=2) call TimerStop('lapack_dsyev: diagonalization')
    !
  end subroutine diag_syev_i8



  subroutine diag_dgelss(a,b,ierror)
    double precision, intent(inout) :: a(:,:)
    double precision, intent(inout) :: b(:,:)
    integer,intent(inout),optional  :: ierror

    external dgelss
    double precision    :: s    (   min(size(a,dim=1),size(a,dim=2))),tol
    double precision    :: work (70*max(size(a,dim=1),size(a,dim=2),size(b,dim=2)))
    double precision,allocatable  :: work_t (:),h(:,:)

    integer(ik)           :: info1,iw,info2
    integer(hik)          :: rank,info
    integer(hik)          :: na1, na2, nb1, nb2
    !   
    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    !
    allocate(h(na1,na2),stat=info1)
    !
    h = a
    !
    tol = -1.0d-12
    !
    call dgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s, tol, rank, work, -1, info)
    !
    iw = int(work(1))
    allocate(work_t(iw),stat=info2)
    !
    call dgelss(na1,na2,nb2,h(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s, tol, rank, work_t, iw, info)
    !
    a = h 
    !
    deallocate(work_t,h)
    !
    if (info/=0.or.info1/=0.or.info2/=0) then
      write (out,"(' dgelss returned ',i8)") info
      stop 'diag_dgelss - dgelss failed'
    end if
    !
    if (rank/=na2.and.present(ierror)) then 
      !
      ierror = 1
      ! 
    endif 
    !
  end subroutine diag_dgelss


  subroutine diag_syev_ilp(h,e,solver,jobz,rng,iroots,vrange,irange,tol)
    double precision, intent(inout) :: h(:,:)
    double precision, intent(out)   :: e(:)
    character(len=cl),intent(in)    :: solver
    character(len=1),intent(in),optional   :: rng,jobz
    double precision, intent(in),optional  :: vrange(2)  ! In:  energy range to be computed
    integer(ik)     , intent(in),optional  :: irange(2)  ! In:  index to be computed
    integer(ik)     , intent(out),optional :: iroots     ! Out: Number of roots found
    double precision, intent(in),optional  :: tol    ! In: abs. tolerance 
    !
    double precision,allocatable    :: work(:)
    double precision,allocatable :: a(:,:)
    integer,allocatable :: iwork(:),ifail(:)
    integer,allocatable :: isuppz(:)
    !
    integer(hik)      :: info,n,lwork,il,ir,m,niwork
    integer(ik)       :: iksize,alloc,k,msize
    !
    double precision :: vl,vu,abstol
    character(len=1) :: rng_,jobz_
    integer(hik) ::  matsize
    !
    if (verbose>=2) call TimerStart('diag_syev_i8')
    !
    if (verbose>=3) write(out,"('Diagonalization with ILP64 ...')") 
    !
    n = size(h,dim=1)
    !
    jobz_ = 'V'
    !
    if (present(jobz)) then
       jobz_ = jobz
    endif
    
    lwork = 50*size(h,dim=1)
    !
    allocate(work(lwork),stat=alloc)
    !
    iksize = lwork
    !
    call ArrayStart('diag_tridiag-arrays',alloc,iksize,rk)
    !
    if (verbose>=2) call TimerStart('diag_dsyev')
    !
    rng_ = 'A'
    il = 1   ; ir = 100000
    vl = -.001 ; vu = 100000.0
    !
    if (present(rng)) then
       rng_ = rng
    endif
    !
    if (present(vrange).and.rng_=='V') then
       vl = vrange(1) ; vu = vrange(2)
    endif
    !
    if (present(irange).and.rng_=='I') then
       il = irange(1) ; ir = irange(2)
       if (irange(2)>=size(h,dim=1)) rng_ = 'A' 
    endif
    !
    abstol =(small_)
    if (present(tol)) then
       abstol = tol
    endif
    !
    msize = n
    !
    if (trim(solver)=='DSYEVR-ILP'.or.trim(solver)=='DSYEVX-ILP') then 
       !
       if (rng=='I') msize = ir-il+1
       !
       allocate (a(n,msize),stat=alloc)
       matsize = int(n,hik)*int(msize,hik)
       call ArrayStart('diag_tridiag-arrays',alloc,1,rk,matsize )
       !
       if (verbose>=3) call MemoryReport
       !
    endif
    !
    m = msize
    !
    allocate(ifail(n),stat=alloc)
    !
    iksize = size(ifail)
    call ArrayStart('diag_tridiag-arrays',alloc,iksize,ik)
    !
    if (verbose>=4) write(out,"(/'abstol = ',e18.8)") abstol
    !
    select case (trim(solver))
      !
    case('DSYEV-ILP')
      !
      call dsyev(jobz_,'U',n,h,n,e,work,-1,info)
      !
      if (int(work(1),hik)>lwork) then 
        !
        lwork = int(work(1),hik)
        deallocate(work)
        allocate(work(lwork),stat=alloc)
        !
        if (verbose>=3) write(out,"(/'lwork = ',i0)") lwork
        !
      endif
      !
      call dsyev(jobz_,'U',n,h,n,e,work,lwork,info)
      !
    case('DSYEVX-ILP')
      !
      niwork = 5*n
      !
      allocate(iwork(niwork),stat=alloc)
      call ArrayStart('diag_tridiag-arrays',alloc,int(niwork,ik),ik)
      !
      call dsyevx('V',rng,'U',n,h,n,vl,vu,il,ir,abstol,m,e,a,n,work,-1,iwork,ifail,info) 
      !
      if (int(work(1))>size(work)) then 
        !
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork),stat=alloc)
        !
      endif
      !
      if (int(iwork(1))>size(work)) then 
        !
        niwork = int(iwork(1))
        deallocate(iwork)
        allocate(iwork(niwork),stat=alloc)
        !
      endif
      !
      call dsyevx('V',rng,'U',n,h,n,vl,vu,il,ir,abstol,m,e,a,n,work,lwork,iwork,ifail,info) 
      !
      if (verbose>=2) call TimerStart('diag_copy')
      !
      !$omp parallel do private(k) shared(h) schedule(dynamic)
      do k=1,m
         call dcopy(n,a(:,k),1,h(:,k),1)
      enddo
      !$omp end parallel do
      !
      if (verbose>=2) call TimerStop('diag_copy')
      !
      deallocate(a)
      !
    case('DSYEVR-ILP')
      !
      niwork = 50*n
      !
      allocate(iwork(niwork),stat=alloc)
      call ArrayStart('diag_tridiag-arrays',alloc,int(niwork,ik),ik)
      !
      allocate(isuppz(2*m),stat=alloc)
      iksize = 2*m
      call ArrayStart('diag_tridiag-arrays',alloc,iksize,ik)
      !
      call dsyevr(jobz_,rng_,'U',n,h,n,vl,vu,il,ir,abstol,m,e,h,n,isuppz,work,-1,iwork,-1,info)
      !
      if (int(work(1))>size(work)) then 
        !
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork),stat=alloc)
        !
      endif
      !
      if (int(iwork(1))>size(work)) then 
        !
        niwork = int(iwork(1))
        deallocate(iwork)
        allocate(iwork(niwork),stat=alloc)
        !
      endif
      !
      call dsyevr(jobz_,rng_,'U',n,h,n,vl,vu,il,ir,abstol,m,e,a,n,isuppz,work,lwork,iwork,niwork,info)
      !
      if (verbose>=2) call TimerStart('diag_copy')
      !
      !$omp parallel do private(k) shared(h) schedule(dynamic)
      do k=1,m
         call dcopy(n,a(:,k),1,h(:,k),1)
      enddo
      !$omp end parallel do
      !
      if (verbose>=2) call TimerStop('diag_copy')
      !
      deallocate(a)
      !
    case('DSYEVD-ILP')
      !
      niwork = 50*n
      !
      allocate(iwork(niwork),stat=alloc)
      call ArrayStart('diag_tridiag-arrays',alloc,int(niwork,ik),ik)
      !
      call dsyevd('V','U',n,h,n,e,work,-1,iwork,-1,info)
      !
      if (int(work(1))>size(work)) then 
        !
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork),stat=alloc)
        !
      endif
      !
      if (int(iwork(1))>size(work)) then 
        !
        niwork = int(iwork(1))
        deallocate(iwork)
        allocate(iwork(niwork),stat=alloc)
        !
      endif
      !
      call dsyevd('V','U',n,h,n,e,work,lwork,iwork,niwork,info)
      !
    end select
    !
    if (verbose>=2) call TimerStop('diag_dsyev')
    !
    if (verbose>=3) write(out,"('...done!')") 
    !
    if (info/=0) then
      write (out,"(' dsyev returned ',i0)") info
      stop 'diag - dsyev failed'
    end if
    !
    if (allocated(isuppz)) deallocate(isuppz)
    if (allocated(ifail)) deallocate(ifail)
    if (allocated(iwork)) deallocate(iwork)
    !
    if (present(iroots))  iroots = m
    !
    deallocate(work)
    !
    call ArrayStop('diag_tridiag-arrays')
    !
    if (verbose>=2) call TimerStop('diag_syev_i8')
    !
  end subroutine diag_syev_ilp



  subroutine diag_dsytrd( uplo, n, a, d, e, tau, info )
      !
      !     .. scalar arguments ..
      character(len=1) ::  uplo
      integer          ::  info, n
      !     ..
      !     .. array arguments ..
      !
      double precision ::  a(:,:),d(:), e(:), tau(:)

      double precision :: one, zero, half
      parameter          ( one = 1.0d0, zero = 0.0d0,half = 1.0d0 / 2.0d0 )
         !     ..
         !     .. local scalars ..
      integer            i, i1, i1i1, ii, k, j, ij1, il1, l, kk, kp1, jl, jj, m, ii1, lda, ldw, nb, nx
      double precision   alpha, taui, temp, temp1, temp2
      !
      integer(ik)  :: alloc
      integer(hik) :: matsize  
      !
      double precision,allocatable ::  w(:,:)
      !     ..
      !     .. external subroutines ..
            integer :: ilaenv
      !     external           daxpy, dlarfg, dspmv, dspr2, xerbla
      !     ..
      !     .. external functions ..
      !     logical            lsame
      !     double precision   :: ddot
         !     external           lsame, ddot
         !
         info = 0
         !
         if (verbose>=2) call TimerStart('diag_dsytrd_in')
         !
         if (verbose>=3) write(out,"('Tridiagonal form with diag_dsytrd...')") 
         !
         if (n<=0) then 
            info = -2
            return
         endif

         lda = n
         ldw = n
         !
         nb = ilaenv( 1, 'dsytrd', uplo, n, -1, -1, -1 )
         nx = max( nb, ilaenv( 3, 'dsytrd', uplo, n, -1, -1, -1 ) )
         nx = min(nx,n)
         !
         if (verbose>=3) write(out,"('the optimal block size nb,nx =  ',2i0)") nb,nx
         !
         allocate(w(ldw,nb),stat=alloc)
         !
         matsize = nb*ldw
         !
         call ArrayStart('diag_tridiag-arrays',alloc,1,ik,matsize)
         !
         kk = n - ( ( n-nx+nb-1 ) / nb )*nb
         !
         do i = n - nb + 1, kk + 1, -nb
            !
            if (verbose>=3) write(out,"('  i = ',i0)") i
            !
            !reduce columns i:i+nb-1 to tridiagonal form and form the
            !matrix w which is needed to update the unreduced part of
            !the matrix
            !
            if (verbose>=2) call TimerStart('dlatrd')
            !
            if (verbose>=3) write(out,"('  dlatrd...')") 
            !
            !call dlatrd('U', i+nb-1, nb, a, lda, e, tau, w, ldw )
            !
            call diag_dlatrd('U', i+nb-1, nb, a, lda, e, tau, w, ldw )
            !
            if (verbose>=3) write(out,"('  ...done!')") 
            !
            if (verbose>=2) call TimerStop('dlatrd')
            !
            !update the unreduced submatrix a(1:i-1,1:i-1), using an
            !update of the form:  a := a - v*w' - w*v'
            !
            if (verbose>=2) call TimerStart('dsyr2k')
            !
            if (verbose>=3) write(out,"('  dsyr2k...')") 
            !
            call dsyr2k('U','N',i-1,nb,-one,a(1,i),lda,w,ldw,one,a,lda )
            !
            if (verbose>=3) write(out,"('  ...done!')") 
            !
            if (verbose>=2) call TimerStop('dsyr2k')
            !
            !copy superdiagonal elements back into a, and diagonal
            !elements into d
            !
            !$omp parallel do private(j) shared(a,d) schedule(guided)
            do j = i, i + nb - 1
               a( j-1, j ) = e( j-1 )
               d( j ) = a( j, j )
            enddo
            !$omp end parallel do
         enddo
         !
         ! Use unblocked code to reduce the last or only block
         !
         !
         if (verbose>=2) call TimerStart('dsyr2k')
         !
         if (verbose>=3) write(out,"('  dsyr2k...')") 
         !
         call dsytd2('U',kk,a,lda,d,e,tau,info)
         !
         !
         if (verbose>=3) write(out,"('  ...done!')") 
         !
         if (verbose>=2) call TimerStop('dsyr2k')
         !
         deallocate(w)
         !
         call ArrayStop('diag_tridiag-arrays')
         !
         if (verbose>=3) write(out,"('...done!')") 
         !
         if (verbose>=2) call TimerStop('diag_dsytrd_in')
         !
  end subroutine diag_dsytrd



  subroutine diag_dlatrd(uplo,n,nb,a,lda,e,tau,w,ldw)
      !
      character   ::       uplo
      integer     ::       lda, ldw, n, nb
      !     ..
      !     .. array arguments ..
      !
      double precision  :: a( lda, * ), e( * ), tau( * ), w( ldw, * )
      !
      !     .. parameters ..
      double precision   zero, one, half
      parameter          ( zero = 0.0d+0, one = 1.0d+0, half = 0.5d+0 )
      !
      !     .. local scalars ..
      integer            i, iw
      double precision   alpha
      !     ..
      !     .. external subroutines ..
      external           daxpy, dgemv, dlarfg, dscal, dsymv
      !     ..
      !     .. external functions ..
      logical            lsame
      double precision   ddot
      external           lsame, ddot
      !     ..
      !     .. intrinsic functions ..
      intrinsic          min

      !
      do i = n, n - nb + 1, -1
         iw = i - n + nb
         if( i.lt.n ) then
            !
            !update a(1:i,i)
            !
            if (verbose>=2) call TimerStart('dgemv')
            !
            if (verbose>=3) write(out,"('    dgemv-2...')")
            !
            call dgemv( 'no transpose', i, n-i, -one, a( 1, i+1 ),lda, w( i, iw+1 ), ldw, one, a( 1, i ), 1 )
            call dgemv( 'no transpose', i, n-i, -one, w( 1, iw+1 ),ldw, a( i, i+1 ), lda, one, a( 1, i ), 1 )
            !
            if (verbose>=3) write(out,"('    ...done!')")
            !
            if (verbose>=2) call TimerStop('dgemv')
            !
         end if
         !
         if( i.gt.1 ) then
            !
            !generate elementary reflector h(i) to annihilate
            !a(1:i-2,i)
            !
            if (verbose>=2) call TimerStart('dlarfg')
            !
            if (verbose>=3) write(out,"('    dlarfg...')")
            !
            call dlarfg( i-1, a( i-1, i ), a( 1, i ), 1, tau( i-1 ) )
            !
            if (verbose>=3) write(out,"('    ...done')")
            !
            if (verbose>=2) call TimerStop('dlarfg')
            !
            e( i-1 ) = a( i-1, i )
            a( i-1, i ) = one
            !
            !compute w(1:i-1,i)
            !
            if (verbose>=2) call TimerStart('dsymv')
            !
            if (verbose>=3) write(out,"('    dsymv...')")
            !
            call diag_dsymv('U', i-1, one, a, lda, a( 1, i ), 1,zero, w( 1, iw ), 1 )
            !
            if (verbose>=3) write(out,"('    ...done!')")
            !
            if (verbose>=2) call TimerStop('dsymv')
            !
            if( i.lt.n ) then
               !
               if (verbose>=2) call TimerStart('dgemv')
               !
               if (verbose>=3) write(out,"('    dgemv-4...')")
               !
               call diag_dgemv( 'transpose', i-1, n-i, one, w( 1, iw+1 ),ldw, a( 1, i ), 1, zero, w( i+1, iw ), 1 )
               call diag_dgemv( 'no transpose', i-1, n-i, -one,a( 1, i+1 ), lda, w( i+1, iw ), 1, one,w( 1, iw ), 1 )
               call diag_dgemv( 'transpose', i-1, n-i, one, a( 1, i+1 ),lda, a( 1, i ), 1, zero, w( i+1, iw ), 1 )
               call diag_dgemv( 'no transpose', i-1, n-i, -one,w( 1, iw+1 ), ldw, w( i+1, iw ), 1, one,w( 1, iw ), 1 )
               !
               if (verbose>=3) write(out,"('    ...done')")
               !
               if (verbose>=2) call TimerStop('dgemv')
               !
            end if
            !
            if (verbose>=2) call TimerStart('dscal-daxpy')
            !
            call dscal( i-1, tau( i-1 ), w( 1, iw ), 1 )
            alpha = -half*tau( i-1 )*ddot( i-1, w( 1, iw ), 1, a( 1, i ), 1 )
            call daxpy( i-1, alpha, a( 1, i ), 1, w( 1, iw ), 1 )
            !
            if (verbose>=2) call TimerStop('dscal-daxpy')
            !
         end if
         !
      enddo
      !
  end subroutine diag_dlatrd



  subroutine diag_dsymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
      !.. scalar arguments ..
      double precision alpha,beta
      integer incx,incy,lda,n
      character uplo
      !..
      !.. array arguments ..
      double precision a(lda,*),x(*),y(*)
      !..
      !
      !  purpose
      !   =======
      ! 
      !   dsymv  performs the matrix-vector  operation
      ! 
      !y := alpha*a*x + beta*y,
      ! 
      !   where alpha and beta are scalars, x and y are n element vectors and
      !   a is an n by n symmetric matrix.
      !.. parameters ..
      double precision one,zero
      parameter (one=1.0d+0,zero=0.0d+0)
      !..
      !.. local scalars ..
      double precision temp1,temp2
      integer i,info,ix,iy,j
      !..
      !.. external functions ..
      logical lsame
      external lsame
      !..
      !.. external subroutines ..
      external xerbla
      !..
      !.. intrinsic functions ..
      intrinsic max
      !..
      !
      !test the input parameters.
      !
      info = 0
      if (.not.lsame(uplo,'u') .and. .not.lsame(uplo,'l')) then
          info = 1
      else if (n.lt.0) then
          info = 2
      else if (lda.lt.max(1,n)) then
          info = 5
      else if (incx/=1) then
          info = 7
      else if (incy/=1) then
          info = 10
      end if
      if (info.ne.0) then
          call xerbla('dsymv ',info)
          return
      end if
      !
      !quick return if possible.
      !
      if ((n.eq.0) .or. ((alpha.eq.zero).and. (beta.eq.one))) return
      !
      !start the operations. in this version the elements of a are
      !accessed sequentially with one pass through the triangular part
      !of a.
      !
      !first form  y := beta*y.
      !
      if (beta.ne.one) then
          if (beta.eq.zero) then
              !
              !$omp parallel do private(i) shared(y) schedule(guided)
              do i = 1,n
                  y(i) = zero
              end do
              !$omp end parallel do
          else
              !$omp parallel do private(i) shared(y) schedule(guided)
              do i = 1,n
                  y(i) = beta*y(i)
              end do
              !$omp end parallel do
          end if
      end if
      !
      if (alpha.eq.zero) return
      !
      if (lsame(uplo,'u')) then
          !
          !   form  y  when a is stored in upper triangle.
          !
          !$omp parallel do private(j,temp2,i) shared(y) schedule(guided)
          do j = 1,n
              !
              temp2 = zero
              !
              do i = 1,j
                  !
                  !y(i) = y(i) + temp1*a(i,j)
                  temp2 = temp2 + a(i,j)*x(i)
                  !
              end do
              !
              do i = j + 1,n
                  temp2 = temp2 + a(j,i)*x(i)
              end do
              !
              y(j) = y(j) + alpha*temp2
              !
          end do
          !$omp end parallel do
          !
      else
          !
          !   form  y  when a is stored in lower triangle.
          !
          !$omp parallel do private(j,temp2,i) shared(y) schedule(guided)
          do j = 1,n
              !
              temp2 = zero
              !
              !y(j) = y(j) + temp1*a(j,j)
              !
              do i = j ,n
                  !y(i) = y(i) + temp1*a(i,j)
                  !
                  temp2 = temp2 + a(i,j)*x(i)
                  !
              end do
              !
              do i = 1,j - 1 
                  temp2 = temp2 + a(j,i)*x(i)
              end do
              !
              y(j) = y(j) + alpha*temp2
              !
          end do
          !$omp end parallel do
          !
      end if
      !
      !end of dsymv .
      !
  end subroutine diag_dsymv




  subroutine diag_dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
      !.. scalar arguments ..
      double precision alpha,beta
      integer lda,m,n,incx,incy
      character trans
      !..
      !.. array arguments ..
      double precision a(lda,*),x(*),y(*)
      !..
      !
      !  purpose
      !  =======
      !
      !  dgemv  performs one of the matrix-vector operations
      !
      !y := alpha*a*x + beta*y,   or   y := alpha*a'*x + beta*y,
      !
      !  where alpha and beta are scalars, x and y are vectors and a is an
      !  m by n matrix.
      !
      !.. parameters ..
      double precision one,zero
      parameter (one=1.0d+0,zero=0.0d+0)
      !..
      !.. local scalars ..
      double precision temp
      integer i,info,ix,iy,j, lenx, leny
      !..
      !.. external functions ..
      logical lsame
      external lsame
      !..
      !.. external subroutines ..
      external xerbla
      !..
      !.. intrinsic functions ..
      intrinsic max
      !..
      !
      !test the input parameters.
      !

      !
      info = 0
      if (.not.lsame(trans,'n') .and. .not.lsame(trans,'t') .and..not.lsame(trans,'c')) then
          info = 1
      else if (m.lt.0) then
          info = 2
      else if (n.lt.0) then
          info = 3
      else if (lda.lt.max(1,m)) then
          info = 6
      else if (incx/=1) then
          info = 8
      else if (incy/=1) then
          info = 11
      end if
      if (info.ne.0) then
          call xerbla('dgemv ',info)
          return
      end if
      !
      !quick return if possible.
      !
      if ((m.eq.0) .or. (n.eq.0) .or.((alpha.eq.zero).and. (beta.eq.one))) return
      !
      if (lsame(trans,'n')) then
          lenx = n
          leny = m
      else
          lenx = m
          leny = n
      end if
      !
      !start the operations. in this version the elements of a are
      !accessed sequentially with one pass through a.
      !
      !first form  y := beta*y.
      !
      if (beta.ne.one) then
          if (beta.eq.zero) then
              !
              !$omp parallel do private(i) shared(y) schedule(guided)
              do i = 1,leny
                  y(i) = zero
              enddo
              !$omp end parallel do
              !
          else
              !
              !$omp parallel do private(i) shared(y) schedule(guided)
              do i = 1,leny
                  y(i) = beta*y(i)
              enddo
              !$omp end parallel do
              !
          end if
      end if
      if (alpha.eq.zero) return
      if (lsame(trans,'n')) then
      !
      !   form  y := alpha*a*x + y.
      !
         !$omp parallel do private(j,temp,i) shared(y) schedule(guided)
         do j = 1,m
              temp = zero
              do i = 1,n
                  temp = temp + a(j,i)*x(i)
              end do
              y(j) = y(j) + alpha*temp
         end do
         !$omp end parallel do
      else
      !
      !   form  y := alpha*a'*x + y.
      !
          !$omp parallel do private(j,temp,i) shared(y) schedule(guided)
          do j = 1,n
              temp = zero
              do i = 1,m
                  temp = temp + a(i,j)*x(i)
              end do
              y(j) = y(j) + alpha*temp
          end do
          !$omp end parallel do
      end if
      !
      !end of dgemv .
      !
  end subroutine diag_dgemv




  subroutine diag_dsptrd( uplo, n, ap, d, e, tau, info )
      !
      !     .. scalar arguments ..
      character(len=1) ::  uplo
      integer          ::  info, n
      !     ..
      !     .. array arguments ..
      !
      double precision ::  ap(:),d(:), e(:), tau(:)

      double precision :: one, zero, half
      parameter          ( one = 1.0d0, zero = 0.0d0,half = 1.0d0 / 2.0d0 )
         !     ..
         !     .. local scalars ..
      integer            i, i1, i1i1, ii, k, j, ij1, il1, l, kk, kp1, jl, jj, m, ii1
      double precision   alpha, taui, temp, temp1, temp2
      !     ..
      !     .. external subroutines ..
      !     external           daxpy, dlarfg, dspmv, dspr2, xerbla
      !     ..
      !     .. external functions ..
      !     logical            lsame
      double precision   :: ddot
         !     external           lsame, ddot
         !
         info = 0
         !
         if (verbose>=2) call TimerStart('diag_dsptrd_in')
         !
         if (verbose>=3) write(out,"('Tridiagonal form with diag_dsptrd...')") 
         !
         if (n<=0) then 
            info = -2
            return
         endif
         ! reduce the upper triangle of a.
         ! i1 is the index in ap of a(1,i+1).
         !
         i1 = n*( n-1 ) / 2 + 1
         !
         if (verbose>=3) write(out,"('i1 = ',i0)") i1
         !
         do i = n - 1, 1, -1
            !
            if (verbose>=3) write(out,"('   i = ',i0)") i
            !
            ! generate elementary reflector h(i) = i - tau * v * v'
            ! to annihilate a(1:i-1,i+1)
            !
            if (verbose>=2) call TimerStart('dlarfg')
            !
            if (verbose>=3) write(out,"('  dlarfg...')") 
            !
            call dlarfg( i, ap( i1+i-1 ), ap( i1 ), 1, taui )
            !
            if (verbose>=3) write(out,"('  ...done!')") 
            !
            if (verbose>=2) call TimerStop('dlarfg')
            !
            e( i ) = ap( i1+i-1 )
            !
            if( taui.ne.zero ) then
               !
               ! apply h(i) from both sides to a(1:i,1:i)
               !
               ap( i1+i-1 ) = one
               !
               ! compute  y := tau * a * v  storing y in tau(1:i)
               !
               if (verbose>=2) call TimerStart('dspmv')
               !
               if (verbose>=3) write(out,"('  dspmv...')") 
               !
               !call dspmv(uplo,i,taui,ap,ap(i1),1,zero,tau,1 )
               !
               ! AP(i + (j-1)*j/2) = A(i,j)
               !
               !$omp parallel do private(j) shared(tau) schedule(guided)
               do j = 1,i
                  tau(j) = zero
               enddo
               !$omp end parallel do 
               !
               !kk = 1
               !
               !$omp parallel do private(j,temp,k,m,l) shared(tau) schedule(guided)
               do j = 1,i
                  !
                  temp = zero
                  !
                  k = i1
                  m = (j-1)*j/2 + 1
                  !
                  do l = 1,j
                      !
                      temp = temp + ap(k)*Ap(m)
                      !
                      k = k + 1
                      m = m + 1
                      !
                  enddo
                  !
                  do l = j + 1,i 
                      !
                      m = j + (l-1)*l/2
                      !
                      temp = temp + Ap(m)*ap(k)
                      !
                      k = k + 1
                      !
                  enddo
                  !
                  tau(j) = taui*temp
                  !
               enddo
               !$omp end parallel do
               !
               if (verbose>=3) write(out,"('  ...done!')") 
               !
               if (verbose>=2) call TimerStop('dspmv')
               !
               if (verbose>=2) call TimerStart('daxpy')
               !
               if (verbose>=3) write(out,"('  daxpy...')") 
               !
               ! compute  w := y - 1/2 * tau * (y'*v) * v
               !
               alpha = -half*taui*ddot( i, tau, 1, ap( i1 ), 1 )
               !
               !call daxpy( i, alpha, ap( i1 ), 1, tau, 1 )
               !
               k = mod(i,4)
               if (k/=0) then
                 !
                 !$omp parallel do private(j,ij1) shared(tau) schedule(guided)
                 do j = 1,k
                   !
                   ij1 = i1 + j - 1
                   !
                   tau(j) = tau(j) + alpha*ap(ij1)
                 enddo
                 !$omp end parallel do
                 !
               endif
               !
               if (i>=4) then 
                 !
                 kp1 = k + 1
                 ii1 = i1 - 1
                 !$omp parallel do private(j,ij1) shared(tau) schedule(guided)
                 do j = kp1,i,4
                   !
                   ij1 = ii1 + j
                   !
                   tau(j)   = tau(j  ) + alpha*ap(ij1  )
                   tau(j+1) = tau(j+1) + alpha*ap(ij1+1)
                   tau(j+2) = tau(j+2) + alpha*ap(ij1+2)
                   tau(j+3) = tau(j+3) + alpha*ap(ij1+3)
                   !
                 enddo 
                 !$omp end parallel do
                 !
               endif
               !
               if (verbose>=3) write(out,"('  ...done!')") 
               !
               if (verbose>=2) call TimerStop('daxpy')
               !
               ! apply the transformation as a rank-2 update:
               !    a := a - v * w' - w * v'
               !
               if (verbose>=2) call TimerStart('dspr2')
               !
               if (verbose>=3) write(out,"('  dspr2...')") 
               !
               !call dspr2( uplo, i, -one, ap( i1 ), 1, tau, 1, ap )
               !
               alpha = -one
               !
               !$omp parallel do private(j,ij1,temp1,temp2,k,m,l) shared(ap) schedule(guided)
               do j = 1,i
                  !
                  ij1 = i1 + j - 1
                  !
                  if ((ap(ij1).ne.zero) .or. (tau(j).ne.zero)) then
                      !
                      temp1 = alpha*tau(j)
                      temp2 = alpha*ap(ij1)
                      !
                      k = (j-1)*j/2 + 1
                      m = i1
                      !
                      do l = 1,j
                          !
                          ap(k) = ap(k) + ap(m)*temp1 + tau(l)*temp2
                          !
                          k = k + 1
                          m = m + 1
                          !
                      enddo
                      !
                  end if
                  !
               enddo
               !$omp end parallel do
               !
               if (verbose>=3) write(out,"('  ...done!')") 
               !
               if (verbose>=2) call TimerStop('dspr2')
               !
               ap( i1+i-1 ) = e( i )
               !
            end if
            !
            d( i+1 ) = ap( i1+i )
            tau( i ) = taui
            i1 = i1 - i
            !
         enddo
         !
         d( 1 ) = ap( 1 )
         !
         if (verbose>=3) write(out,"('...done!')") 
         !
         if (verbose>=2) call TimerStop('diag_dsptrd_in')
         !
  end subroutine diag_dsptrd




  subroutine diag_dsptrd_omp( uplo, n, ap, d, e, tau, info )
      !
      !     .. scalar arguments ..
      character(len=1) ::  uplo
      integer          ::  info, n
      !     ..
      !     .. array arguments ..
      !
      double precision ::  ap(:),d(:), e(:), tau(:)

      double precision :: one, zero, half
      parameter          ( one = 1.0d0, zero = 0.0d0,half = 1.0d0 / 2.0d0 )
         !     ..
         !     .. local scalars ..
      integer            i, i1, i1i1, ii, np, ip, j, alloc, i0
      double precision   alpha, taui
      !     ..
      double precision,allocatable ::  tau_(:,:),taui_(:)
      !
      !     .. external subroutines ..
      !     external           daxpy, dlarfg, dspmv, dspr2, xerbla
      !     ..
      !     .. external functions ..
      !     logical            lsame
      double precision   :: ddot
         !     external           lsame, ddot
         !
         info = 0
         !
         if (verbose>=2) call TimerStart('diag_dsptrd_omp')
         !
         if (verbose>=3) write(out,"('Tridiagonal form with diag_dsptrd_omp...')")
         !
         np = 2
         allocate (tau_(n-1,np),taui_(np),stat=alloc)     
         !
         if (n<=0) then 
            info = -2
            return
         endif
         ! reduce the upper triangle of a.
         ! i1 is the index in ap of a(1,i+1).
         !
         i1 = n*( n-1 ) / 2 + 1
         !
         if (verbose>=3) write(out,"('i1 = ',i0)") i1
         !
         do i = n - 1, np+1, -np
            !
            i0 = max(i-np+1,np+1)
            !
            !$omp parallel do private(j,taui,ip,alpha) shared(ap,e,tau_) schedule(guided)
            do j = i,i0,-1
              !
              ip = i-j+1
              !
              if (verbose>=3) write(out,"('   i,j = ',2i0)") i,j
              !
              ! generate elementary reflector h(i) = i - tau * v * v'
              ! to annihilate a(1:i-1,i+1)
              !
              call dlarfg( j, ap( i1+j-1 ), ap( i1 ), 1, taui )
              !
              taui_(ip) = taui
              !
              e( j ) = ap( i1+j-1 )
              !
              if( taui.ne.zero ) then
                 !
                 ! apply h(i) from both sides to a(1:i,1:i)
                 !
                 ap( i1+j-1 ) = one
                 !
                 ! compute  y := tau * a * v  storing y in tau(1:i)
                 !
                 call dspmv(uplo,j,taui,ap,ap(i1),1,zero,tau_(:,ip),1 )
                 !
                 ! compute  w := y - 1/2 * tau * (y'*v) * v
                 !
                 alpha = -half*taui*ddot( j, tau_(:,ip), 1, ap( i1 ), 1 )
                 !
                 call daxpy( j, alpha, ap( i1 ), 1, tau_(:,ip), 1 )
                 !
              end if
              !
            enddo
            !$omp end parallel do 
            !
            if (verbose>=2) call TimerStart('dspr2-comb')
            !
            if (verbose>=3) write(out,"('  dspr2-comb...')") 
            !
            do j = i,i0,-1
              !
              ip = i-j+1
              !
              if (verbose>=3) write(out,"('   j = ',i0)") j
              !
              taui = taui_(ip)
              !
              if( taui.ne.zero ) then
                 !
                 ! apply the transformation as a rank-2 update:
                 !    a := a - v * w' - w * v'
                 !
                 call dspr2( uplo, j, -one, ap( i1 ), 1, tau_(:,np), 1, ap )
                 !
                 tau(:) = tau_(:,np)
                 !
                 ap( i1+j-1 ) = e( j )
                 !
              end if
              !
              d( j+1 ) = ap( i1+j )
              tau( j ) = taui
              !
            enddo
            !
            if (verbose>=3) write(out,"('  ...done!')") 
            !
            if (verbose>=2) call TimerStop('dspr2-comb')
            !
            !i1 = i1 - i
            !
            i1 = i0*( i0-1 ) / 2 + 1
            !
         enddo
         !
         do i = min(np,n-1), 1, -1
            !
            if (verbose>=3) write(out,"('   i = ',i0)") i
            !
            ! generate elementary reflector h(i) = i - tau * v * v'
            ! to annihilate a(1:i-1,i+1)
            !
            if (verbose>=2) call TimerStart('dlarfg')
            !
            if (verbose>=3) write(out,"('  dlarfg...')") 
            !
            call dlarfg( i, ap( i1+i-1 ), ap( i1 ), 1, taui )
            !
            if (verbose>=3) write(out,"('  ...done!')") 
            !
            if (verbose>=2) call TimerStop('dlarfg')
            !
            e( i ) = ap( i1+i-1 )
            !
            if( taui.ne.zero ) then
               !
               ! apply h(i) from both sides to a(1:i,1:i)
               !
               ap( i1+i-1 ) = one
               !
               ! compute  y := tau * a * v  storing y in tau(1:i)
               !
               if (verbose>=2) call TimerStart('daxpy')
               !
               if (verbose>=3) write(out,"('  daxpy...')") 
               !
               call dspmv(uplo,i,taui,ap,ap(i1),1,zero,tau,1 )
               !
               ! compute  w := y - 1/2 * tau * (y'*v) * v
               !
               alpha = -half*taui*ddot( i, tau, 1, ap( i1 ), 1 )
               !
               call daxpy( i, alpha, ap( i1 ), 1, tau, 1 )
               !
               if (verbose>=3) write(out,"('  ...done!')") 
               !
               if (verbose>=2) call TimerStop('daxpy')
               !
               ! apply the transformation as a rank-2 update:
               !    a := a - v * w' - w * v'
               !
               if (verbose>=2) call TimerStart('dspr2')
               !
               if (verbose>=3) write(out,"('  dspr2...')") 
               !
               call dspr2( uplo, i, -one, ap( i1 ), 1, tau, 1, ap )
               !
               if (verbose>=3) write(out,"('  ...done!')") 
               !
               if (verbose>=2) call TimerStop('dspr2')
               !
               ap( i1+i-1 ) = e( i )
               !
            end if
            !
            d( i+1 ) = ap( i1+i )
            tau( i ) = taui
            i1 = i1 - i
            !
         enddo
         !
         d( 1 ) = ap( 1 )
         !
         deallocate(tau_,taui_)
         !
         if (verbose>=3) write(out,"('...done!')") 
         !
         if (verbose>=2) call TimerStop('diag_dsptrd_omp')
         !
  end subroutine diag_dsptrd_omp



  subroutine diag_dseupd(n,bterm,nroots,factor,maxitr_,tol,h,e)

    integer,intent(in)     :: n,bterm(n,2),maxitr_
    double precision,intent(in)    :: factor
    integer,intent(inout)  :: nroots
    double precision,intent(in)    :: tol
    double precision,intent(inout) :: h(:,:)
    double precision,intent(out)   :: e(nroots)

    double precision,allocatable :: v(:,:),workl(:),workd(:),d(:,:),resid(:)
    !
    logical,allocatable   ::   select(:)
    integer(ik) ::   iparam(11), ipntr(11),ldv,iter
    integer(hik) ::  msize
    !

    character(len=1) ::  bmat
    character(len=2) ::  which
    !
    integer(ik)   :: ido,  nev, ncv, lworkl, info, ierr, j, &
                     mode, ishfts, alloc, maxitr,k
                     !
    logical          rvec
    double precision  ::  sigma
    !
    !double precision,external :: dnrm2
    !
      integer  logfil, ndigit, mgetv0,&
              msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,&
              mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,&
              mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ &
              logfil, ndigit, mgetv0,&
              msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,&
              mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,&
              mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd


      if (verbose>=2) call TimerStart('diag_dseupd: diagonalization')
      !
      ndigit = -3
      logfil = out
      msgets = 0
      msaitr = 0
      msapps = 0
      msaupd = 0
      msaup2 = 0
      mseigt = 0
      mseupd = 0
      !
      if (verbose>=3) msaupd = 1
      !
      nev = nroots
      !
      if (nroots==0) nev=max(100,n)
      !
      ncv = factor*nev
      !
      if (ncv==nev) ncv = int(21*nev/10)+1
      !
      ncv = min(n,ncv)
      !
      bmat = 'I'
      which = 'SM'
      !
      ! The work array WORKL is used in DSAUPD as workspace. 
      ! The parameter TOL determines the stopping criterion.  
      ! If TOL<=0, machine precision is used.  The variable IDO is used for 
      ! reverse communication and is initially set to 0. 
      ! Setting INFO=0 indicates that a random vector is generated in DSAUPD 
      ! to start the Arnoldi  iteration.  
      !
 
      lworkl = ncv*(ncv+10)
      ldv = n
      !
      !tol = tol
      info = 0
      ido = 0
      !
      allocate(v(ldv,ncv), workl(lworkl),workd(3*n),d(ncv,2),resid(n),select(ncv),stat=alloc)
      !
      msize = int(ldv,hik)*int(ncv,hik)
      !
      if (verbose>=3) write(out,"(/'Arpack: msize = ',i0)") msize
      !
      call ArrayStart('diag_dseupd:v',alloc,1,kind(v),msize)
      call ArrayStart('diag_dseupd:workl',alloc,size(workl),kind(workl))
      call ArrayStart('diag_dseupd:workd',alloc,size(workd),kind(workd))
      call ArrayStart('diag_dseupd:d',alloc,size(d),kind(d))
      call ArrayStart('diag_dseupd:resid',alloc,size(resid),kind(resid))
      call ArrayStart('diag_dseupd:select',alloc,size(select),kind(select))
      !
      !
      iparam = 0 
      ipntr = 0 
      !
      ! This routone uses exact shifts with respect to    
      ! the current Hessenberg matrix (IPARAM(1) = 1).    
      ! IPARAM(3) specifies the maximum number of Arnoldi 
      ! iterations allowed.  Mode 1 of DSAUPD is used     
      ! (IPARAM(7) = 1).  
      !
      ishfts = 1
      maxitr = maxitr_*nev
      mode   = 1
      !      
      iparam(1) = ishfts 
      iparam(3) = maxitr 
      iparam(7) = mode 
      !
      !
      ! M A I N   L O O P (Reverse communication) 
      !
      !if (verbose>=4) write(out,"(/'Arpack: eigenvalues computed:')") 
      !
      if (verbose>=3) write(out,"(/'Arpack: n = ',i8,' nev = ',i8,' ncv = ',i8,' maxitr = ',i8,' tol = ',e12.5)") n,nev,ncv,maxitr,tol
      !
      iter = 0 
      !
      do while(ido<=1)
        !
        iter = iter + 1
        !
        if (verbose>=3.and.mod(iter,500)==0) write(out,"(' iter  = ',i8)") iter 
        !
        ! Repeatedly call the routine DSAUPD and take 
        ! actions indicated by parameter IDO until    
        ! either convergence is indicated or maxitr   
        ! has been exceeded.                          
        !

        !dec$ if (arpack_ > 0)
           !
           call dsaupd ( ido, bmat, n, which, nev, tol, resid, &
                         ncv, v, ldv, iparam, ipntr, workd, workl, &
                         lworkl, info )
        !dec$ else 
            !
            write(out,"(/'Arpack was not activated yet. Please uncomment dsaupd and  dseupd bellow')")
            stop 'Arpack was not activated'
            !
        !dec$ end if
        !
        !if (verbose>=4.and.iparam(5)>0) then 
        !  !
        !  write(out,"(i8)") iparam(5)
        !  !
        !endif
        !
        if (ido==-1 .or. ido==1) then 
          !
          ! Perform matrix vector multiplication 
          !        y <--- OP*x 
          !
          !call matvec(n,bterm,h,workd(ipntr(1)), workd(ipntr(2)))
          !
          call matvec_band(n,h,workd(ipntr(1)), workd(ipntr(2)))
          !
        endif 
        !
      enddo
      !
      !  Either we have convergence or there is  an error. 
      !
      if ( info < 0 ) then
         write(out,"(/'Error with _saupd, info = ',i8)") info
         stop 'Error with _saupd'
      endif
      !
      !
      !  No fatal errors occurred. Post-Process using DSEUPD. 
      !  Computed eigenvalues may be extracted. 
      !
      !  Eigenvectors may also be computed now if desired.  (indicated by rvec = .true.)  
      !             
      rvec = .true.
      !
      !dec$ if (arpack_ > 0)
        !
        if (verbose>=5) write(out,"(/'Arpack: dseupd')") 
        !
        if (verbose>=2) call TimerStart('diag_dseupd: Eigenvectors')
        !
        call dseupd ( rvec, 'All', select, d, v, ldv, sigma, &
                      bmat, n, which, nev, tol, resid, ncv, v, ldv, &
                      iparam, ipntr, workd, workl, lworkl, ierr )
        !
        if (verbose>=2) call TimerStop('diag_dseupd: Eigenvectors')
        !
        if (verbose>=5) write(out,"(/'Arpack: done!')") 
        !
      !dec$ else 
          !
          write(out,"(/'Arpack was not activated yet. Please uncomment dsaupd and  dseupd bellow')")
          stop 'Arpack was not activated'
          !
      !dec$ end if
      !
      if ( ierr < 0 ) then
         write(out,"(/'Error with_seupd, info = ',i8)") ierr
         stop 'Error with _saupd'
      endif
      !
      k = size(h,dim=1)
      !
      nroots =  iparam(5)
      !
      nroots = min(nroots,k)
      !
      e(1:nroots) = d(1:nroots,1)
      !
      !$omp parallel do private(j) shared(h) schedule(static)
      do j=1,n
         !
         call dcopy(nroots,v(j,1:nroots),1,h(1:nroots,j),1)
         !
         !h(j,1:nroots) = v(j,1:nroots)
         !
      enddo
      !$omp end parallel do
      ! 
      !  Eigenvalues are returned in the first column 
      !  of the two dimensional array D and the       
      !  corresponding eigenvectors are returned in   
      !  the first NEV columns of the two dimensional 
      !  array V if requested.  Otherwise, an         
      !  orthogonal basis for the invariant subspace  
      !  corresponding to the eigenvalues in D is     
      !  returned in V.                               
      !
      !do j=1, nconv
        !
        ! Compute the residual norm ||  A*x - lambda*x ||  
        ! for the NCONV accurately computed eigenvalues and 
        ! eigenvectors.  (iparam(5) indicates how many are   
        ! accurate to the requested tolerance)               
        !
        !call matvec(n,bterm,v(1,j), ax)
        !
        !call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
        !d(j,2) = dnrm2(n, ax, 1)
        !d(j,2) = d(j,2) / abs(d(j,1))
        !
      !enddo     

      deallocate(v,workl,workd,d,resid,select)
      !
      call ArrayStop('diag_dseupd:v')
      call ArrayStop('diag_dseupd:workl')
      call ArrayStop('diag_dseupd:workd')
      call ArrayStop('diag_dseupd:d')
      call ArrayStop('diag_dseupd:resid')
      call ArrayStop('diag_dseupd:select')
      !
      if (verbose>=2) call TimerStop('diag_dseupd: diagonalization')
      !
  end subroutine diag_dseupd
! 
! ------------------------------------------------------------------
!     matrix vector subroutine
!
  subroutine matvec(n,bterm,h,z,w)
  
     integer,intent(in)  :: n,bterm(n,2)
     double precision,intent(in)  :: h(:,:),z(n)
     double precision,intent(out)  :: w(n)

     ! matvec performs w = hamil * z for SEUPD 
     ! written to take advantage of the block structure
     ! note:
     ! hamil contains arrays diag & offdg relying on them being
     ! adjacent in the dynamic store allocation
  
      double precision,parameter :: alpha = 1.0d0,beta=0.0d0
      integer               :: k,istart,iend,dimen,m
      double precision,external    :: ddot
      !
      !$omp parallel do private(k,istart,iend,dimen) shared(w) schedule(dynamic)
      do k=1,n
         !
         istart = bterm(k,1)
         iend   = bterm(k,2)
         !
         dimen = iend-istart+1
         !
         w(k) = ddot(dimen,h(k,1:dimen),1,z(istart:iend),1)
         !
      enddo 
      !$omp end parallel do 
      !
      !enddo
      !
   return
  end subroutine matvec 



  subroutine matvec_band(n,h,z,w)
  
     integer,intent(in)  :: n
     double precision,intent(in)  :: h(:,:),z(n)
     double precision,intent(out)  :: w(n)

     ! matvec performs w = hamil * z for SEUPD 
     ! written to take advantage of the block structure
     ! note:
     ! hamil contains arrays diag & offdg relying on them being
     ! adjacent in the dynamic store allocation
  
      double precision,parameter :: alpha = 1.0d0,beta=0.0d0
      integer               :: k,istart,iend,dimen,m,lda
      double precision,external    :: ddot
      !
      k = size(h,dim=1)-1
      lda = k+1
      !
      call dsbmv('L',n,k,alpha,h,lda,z,1,beta,w,1)
      !
  end subroutine matvec_band 




  subroutine diag_dsyev_lin(h,e,jobz)
    double precision, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    double precision, intent(out)   :: e(:)    ! Out: Eigenvalues
    character(len=1),intent(in),optional :: jobz
    character(len=1) :: jobz_

    double precision,allocatable    :: work(:),a(:,:),b(:,:),s(:)
    integer(ik)      :: info
    integer(ik)      :: n,i,j,k,rank
    integer(ik)      :: lwork
    !
    integer(hik)      :: asize
    double precision  :: tol,cross_prod,factor,ddot
    !
    if (verbose>=2) call TimerStart('diag_dsyev_lin: diagonalization')
    !
    jobz_ = 'N'
    !
    lwork = 50*size(h,dim=1)
    !
    allocate(work(lwork))
    !
    n = size(h,dim=1) 
    call dsyev(jobz_,'U',n,h(1:n,1:n),n,e,work,-1,info)
    !
    asize = int(n,hik)*int(n,hik)
    !
    allocate (a(n,n),stat=info)
    call ArrayStart('lapack_dsyevr-a',info,1,kind(a),asize)
    !
    !$omp parallel do private(k) shared(a) schedule(dynamic)
    do k=1,n
       a(:,k) = h(:,k)
    enddo
    !$omp end parallel do
    !
    if (int(work(1))>size(work)) then
      !
      lwork = int(work(1))
      !
      if (verbose>=3) write(out,"('lwork = ',i0)")  lwork
      !
      deallocate(work)
      !
      allocate(work(lwork))
      !
      call dsyev(jobz_,'U',n,h(1:n,1:n),n,e,work,lwork,info)
      !
    else
      !
      call dsyev(jobz_,'U',n,h(1:n,1:n),n,e,work,lwork,info)
      !
    endif
    !
    deallocate(work)
    !
    lwork = 10*n
    !
    allocate(work(lwork),b(n,n),s(n))
    !
    asize = asize + int(n,hik)
    !
    call ArrayStart('lapack_dsyevr-a',info,1,kind(b),asize)
    !
    tol = -1.0d-12
    !
    do i = 1,n
      !
      do j = 1,n
        !
        h(:,j) = a(:,j)
        !
      enddo
      !
      do j = 1,n
        !
        h(j,j) = a(j,j) - e(i)
        !
      enddo
      !
      call dgelss(n,n,1,h(1:n,1:n),n,b(1:n,i),n,s,tol,rank,work,-1,info)
      !
      if (int(work(1))>size(work)) then
        !
        deallocate(work)
        !
        lwork = int(work(1))
        !
        allocate(work(lwork))
        ! 
      endif
      !
      call dgelss(n,n,1,h(1:n,1:n),n,b(1:n,i),n,s,tol,rank,work,lwork,info)
      !
      cross_prod = ddot(n,b(:,i),1,b(:,i),1)
      !
      factor = 1.0_rk/sqrt(cross_prod)
      !
      call dscal(n,factor,b(:,i),1)
      !
      !omp parallel do private(jelem,cross_prod) shared(mat) schedule(dynamic)
      do j = 1,i-1
        !
        cross_prod = ddot(n,b(:,j),1,b(:,j),1)
        !
        cross_prod = -cross_prod
        !
        !mat(:,jelem) = mat(:,jelem)-cross_prod*mat(:,ielem)
        !
        call daxpy(n,cross_prod,b(:,j),1,b(:,i),1)
        ! 
      enddo 
      !omp end parallel do
      !
      !
    enddo
    !
    do j = 1,n
      !
      h(:,j) = b(:,j)
      !
    enddo
    !
    deallocate (a,b)
    call ArrayStop('lapack_dsyevr-a')


    if (info/=0) then
      write (out,"(' dsyev returned ',i8)") info
      stop 'lapack_dsyev - dsyev failed'
    end if
    !
    if (verbose>=2) call TimerStop('diag_dsyev_lin: diagonalization')
    !
  end subroutine diag_dsyev_lin


  subroutine matvec_p(comm,n,nloc,nprocs,kstart,dx,bterm,h,mv_buf,z,w)
  
     integer,intent(in)  :: n,nloc,nprocs,bterm(n,2),kstart(0:nprocs),dx
     double precision,intent(in)  :: h(:,:),z(n),mv_buf(n)
     double precision,intent(out)  :: w(n)
     integer(ik) ::    comm, nprow, npcol, myprow, mypcol, myid, ierr, nprocs_

     ! matvec performs w = hamil * z for SEUPD 
     ! written to take advantage of the block structure
     ! note:
     ! hamil contains arrays diag & offdg relying on them being
     ! adjacent in the dynamic store allocation
  
      double precision,parameter :: alpha = 1.0d0,beta=0.0d0
      integer               :: k,istart,iend,dimen,m,kend,iprev,inext,nx
      double precision,external    :: ddot
      integer,parameter  :: MPI_DOUBLE_PRECISION = 17
	  !
	  myid = 1
	  nprow = 1
      !
      !dec$ if (blacs_ > 0)
        call BLACS_GRIDINFO( comm, nprow, npcol, myprow, mypcol )
        myid= myprow
      !dec$ elseif (mpi_ > 0)
        call MPI_COMM_RANK( comm, myid, ierr )
        call MPI_COMM_SIZE( comm, nprocs_, ierr )
        if (nprocs_/=nprocs) then
          write(out,"('matvec_p: inconsistent number of  nprocs  = ',2i0)") nprocs_,nprocs
          stop 'matvec_p: inconsistent number of  nprocs s'
        endif
      !dec$ end if
      !
      kend = kstart(myid) + nloc-1
      !
      !nprow = nprocs
      !

      do k=kstart(myid),kend
         !
         ! main part that belongs to the same block 
         !
         istart = max(bterm(k,1),kstart(myid))
         iend   = min(bterm(k,2),kend)
         !
         dimen = iend-istart+1
         !
         w(k) = ddot(dimen,h(k,istart:iend),1,z(istart:iend),1)
         !
         ! receiving and sending parts from the previos blocks
         !
         do iprev = max(myid-dx,0),myid-1
           !
           !if (bterm(k,1)>=kstart(iprev+1)) cycle
           !
           istart = kstart(iprev)
           !
           iend = kstart(iprev+1)-1
           !
           nx = iend-istart+1
           !
           !dec$ if (blacs_ > 0)
             call dgerv2d( comm, nx, 1, mv_buf(istart:iend), nx, iprev, mypcol )
           !dec$ end if
           !dec$ if (mpi_ > 0)
             call mpi_recv(mv_buf(istart),nx,MPI_DOUBLE_PRECISION,iprev,myid,comm,ierr)
           !dec$ end if
           !
           istart = max(istart,bterm(k,1))
           !
           nx = iend-istart+1
           !
           w(k) = w(k) + ddot(nx,h(k,istart:iend),1,mv_buf(istart:iend),1)
           !
           istart = kstart(myid)
           !
           iend = kstart(myid+1)-1
           !
           nx = iend-istart+1
           !
           !dec$ if (blacs_ > 0)
             call dgesd2d( comm, nx, 1, z(istart:iend), nx, iprev, mypcol )
           !dec$ end if
           !dec$ if (mpi_ > 0)
             call mpi_send(mv_buf(istart),nx,MPI_DOUBLE_PRECISION,iprev,iprev,comm,ierr)
           !dec$ end if
           !
         enddo
         !
         ! receiving and sending parts from the next blocks
         !
         !
         do inext = myid+1,min(nprow-1,myid+dx)
           !
           !if (bterm(k,2)<kstart(inext)) exit
           !
           istart = kstart(inext)
           !
           iend = kstart(inext+1)-1
           !
           nx = iend-istart+1
           !
           !dec$ if (blacs_ > 0)
             call dgerv2d( comm, nx, 1, mv_buf(istart:iend), nx, inext, mypcol )
           !dec$ end if
           !dec$ if (mpi_ > 0)
             call mpi_recv(mv_buf(istart),nx,MPI_DOUBLE_PRECISION,inext,myid,comm,ierr)
           !dec$ end if
           !
           iend = min(bterm(k,2),iend)
           !
           nx = iend-istart+1
           !
           w(k) = w(k) + ddot(nx,h(k,istart:iend),1,mv_buf(istart:iend),1)
           !
           istart = kstart(myid)
           !
           iend = kstart(myid+1)-1
           !
           nx = iend-istart+1
           !
           !dec$ if (blacs_ > 0)
             call dgesd2d( comm, nx, 1, z(istart:iend), nx, inext, mypcol )
           !dec$ end if
           !dec$ if (mpi_ > 0)
             call mpi_send(mv_buf(istart),nx,MPI_DOUBLE_PRECISION,inext,inext,comm,ierr)
           !dec$ end if
           !
         enddo
         !
      enddo 
      !
      !enddo
      !
   return
  end subroutine matvec_p



  subroutine diag_dseupd_p(n,bterm,nroots,factor,maxitr_,tol,h,e)

    integer,intent(in)     :: n,bterm(n,2),maxitr_
    integer,intent(inout)  :: nroots
    double precision,intent(in)    :: factor
    double precision,intent(in)    :: tol
    double precision,intent(inout) :: h(:,:)
    double precision,intent(out)   :: e(nroots)

    double precision,allocatable :: v(:,:),workl(:),workd(:),d(:,:),resid(:)
    !
    logical,allocatable   ::   select(:)
    integer(ik) ::   iparam(11), ipntr(11),ldv,iter
    !
    integer,parameter ::  maxnprocs=1024
    !
    integer(ik)       ::  kstart(0:maxnprocs),iproc,dx,i,jproc
    !
    character(len=1) ::  bmat
    character(len=2) ::  which
    character(len=cl) ::  blacs_or_mpi = 'NONE'
    !
    integer(ik)   :: ido,  nev, ncv, lworkl, info, ierr, j, &
                     mode, ishfts, alloc, maxitr
                     !
    logical          rvec
    double precision  ::  sigma

    !
    !double precision,external :: dnrm2
    !

!
!-----------------------------------------------------------------------
!
      integer  logfil, ndigit, mgetv0,&
              msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,&
              mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,&
              mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ &
              logfil, ndigit, mgetv0,&
              msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,&
              mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,&
              mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd

      real       t0, t1, t2, t3, t4, t5
      save       t0, t1, t2, t3, t4, t5
!
      integer    nopx, nbx, nrorth, nitref, nrstrt
      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,&
                 tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,&
                 tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,&
                 tmvopx, tmvbx, tgetv0, titref, trvec
      common /timing/ &
                 nopx, nbx, nrorth, nitref, nrstrt,&
                 tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,&
                 tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,&
                 tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,&
                 tmvopx, tmvbx, tgetv0, titref, trvec


 
!     %-----------------%
!     | BLACS INTERFACE |
!     %-----------------%
!
      integer           comm, iam, nprocs, nloc, &
                        nprow, npcol, myprow, mypcol
!
      !external          BLACS_PINFO, BLACS_SETUP, BLACS_GET, &
      !                  BLACS_GRIDINIT, BLACS_GRIDINFO
      !integer,parameter :: MPI_COMM_WORLD=0
!
 
!     %---------------%
!     | MPI INTERFACE |
!     %---------------%
 
      integer           myid, rc
!     include 'mpif.h'


!
!     %----------------------------------------------%
!     | Local Buffers needed for BLACS communication |
!     %----------------------------------------------%
!
      Double precision,allocatable   :: mv_buf(:)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Double precision   zero
      parameter        (zero = 0.0)
!  
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Double precision pdnorm2 
      external         pdnorm2, daxpy
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic         abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!

      write(out,"('Start PARPACK-diagonalization')")


      !dec$ if (blacs_ > 0)
        write(out,"('BLAS-PINFO-start')")
        call BLACS_PINFO( iam, nprocs )
        print *,nprocs
        blacs_or_mpi = 'BLACS'
      !dec$ end if
      !
      write(out,"('BLAS-PINFO-done')")
      !
      !call BLACS_PINFO( iam, nprocs )
      !print *,nprocs

      !dec$ if (mpi_ > 0)
        call MPI_INIT( ierr )
        comm = MPI_COMM_WORLD
        call MPI_COMM_RANK( comm, myid, ierr )
        call MPI_COMM_SIZE( comm, nprocs, ierr )
        !
        print *,comm,myid,nprocs
        !
        if (trim(blacs_or_mpi)=='BLACS') then
          write(out,"('diag_dseupd_p: IT IS ILLEGAL TO USE MPI AND BLACS AT THE SAME TIME')")
          stop 'diag_dseupd_p: IT IS ILLEGAL TO USE MPI AND BLACS AT THE SAME TIME'
        endif
        !
        blacs_or_mpi = 'MPI'
        !
      !dec$ end if
      !


!
!     If in PVM, create virtual machine if it doesn't exist
!
      if (nprocs .lt. 1) then
         nprocs = 1
         !dec$ if (blacs_ > 0)
           call BLACS_SETUP( iam, nprocs )
         !dec$ end if
	 !
	 print *,nprocs
	 !
      endif
      if (nprocs >maxnprocs) stop 'nprocs > maxnprocs'
      !
      !     Set up processors in 1D Grid
      !
      nprow = nprocs
      npcol = 1
      !
      do iproc = 0,nprocs
        !
        kstart(iproc) = iproc*(n/nprocs+1)+1
        if (iproc>=mod(n,nprocs)) kstart(iproc) = iproc*(n/nprocs)+1+mod(n,nprocs)
        print *,iproc,kstart(iproc)
        !
      enddo
      !
      if (verbose>=2) call TimerStart('diag_dseupd_p: diagonalization')
      !
      nev = nroots
      !
      if (nroots==0) nev=max(100,n)
      !
      ncv = factor*nev
      !
      if (ncv==nev) ncv = int(21*nev/10)+1
      !
      ncv = min(n,ncv)
      !
      bmat = 'I'
      which = 'SM'
      !
      ! The work array WORKL is used in DSAUPD as workspace. 
      ! The parameter TOL determines the stopping criterion.  
      ! If TOL<=0, machine precision is used.  The variable IDO is used for 
      ! reverse communication and is initially set to 0. 
      ! Setting INFO=0 indicates that a random vector is generated in DSAUPD 
      ! to start the Arnoldi  iteration.  
      !
 
      lworkl = ncv*(ncv+10)
      ldv = n
      !
      info = 0
      ido = 0
      !
      allocate(v(ldv,ncv), workl(lworkl),workd(3*n),d(ncv,2),resid(n),select(ncv),mv_buf(n),stat=alloc)
      call ArrayStart('diag_dseupd_p',alloc,size(v),kind(v))
      call ArrayStart('diag_dseupd_p',alloc,size(workl),kind(workl))
      call ArrayStart('diag_dseupd_p',alloc,size(workd),kind(workd))
      call ArrayStart('diag_dseupd_p',alloc,size(d),kind(d))
      call ArrayStart('diag_dseupd_p',alloc,size(resid),kind(resid))
      call ArrayStart('diag_dseupd_p',alloc,size(select),kind(select))
      call ArrayStart('diag_dseupd_p',alloc,size(mv_buf),kind(mv_buf))

      !
      ! Estimate block overlapping:
      !
      dx = 0
      !
      do iproc = 0,nprocs-1
        !
        do i =  kstart(iproc),kstart(iproc+1)-1
          !
          do jproc = 0,iproc-1
            !
            if (bterm(i,1)<kstart(jproc+1)-1) dx = max(dx,iproc-jproc)
            !
          enddo 
          !
          do jproc = iproc+1,nprocs-1
            !
            if (bterm(i,2)>kstart(jproc)) dx = max(dx,jproc-iproc)
            !
          enddo 
          !
        enddo
        !
      enddo 
      !
!
!     Get default system context, and define grid
!
      !
      myprow = 1 ; mypcol = 1 ; myid = 1
      !
      !dec$ if (blacs_ > 0)
        call BLACS_GET( 0, 0, comm )
        call BLACS_GRIDINIT( comm, 'Row', nprow, npcol )
        call BLACS_GRIDINFO( comm, nprow, npcol, myprow, mypcol )
        !
        write(out,"(' myprow, nprow, mypcol, npcol = ',4i8)") myprow, nprow, mypcol, npcol
        !
      !dec$ end if
      !
      if (verbose>=2.and.trim(blacs_or_mpi)=='BLACS') write(out,"('myprow, nprow, mypcol, npcol, nprocs = ',5i8)") myprow, nprow, mypcol, npcol, nprocs
      if (verbose>=2.and.trim(blacs_or_mpi)=='MPI') write(out,"('myid,nproc = ',2i8)") myid,nprocs
!     !
!     If I'm not in grid, go to end of program
!
      if ( trim(blacs_or_mpi)=='BLACS'.and.(myprow .ge. nprow) .or. (mypcol .ge. npcol) ) then 
        write(out,"('not in grid, myprow, nprow, mypcol, npcol = ',4i8)") myprow, nprow, mypcol, npcol
        return 
      endif 
!
      ndigit = -3
      logfil = 6
      msaupd = 1

 
      !--------------------------------------!
      ! Set up distribution of data to nodes !
      !--------------------------------------!
 
      nloc = (n/nprocs)
      if ( mod(n,nprocs)>myprow ) nloc = nloc + 1
      !
      !
      !allocate(mv_buf(n),stat=alloc)
      !call ArrayStart('diag_dseupd_p',alloc,size(mv_buf),kind(mv_buf))
      !
      iparam = 0 
      ipntr = 0 
      !
      ! This routone uses exact shifts with respect to    
      ! the current Hessenberg matrix (IPARAM(1) = 1).    
      ! IPARAM(3) specifies the maximum number of Arnoldi 
      ! iterations allowed.  Mode 1 of DSAUPD is used     
      ! (IPARAM(7) = 1).  
      !
      ishfts = 1
      maxitr = maxitr_*nev
      mode   = 1
      !      
      iparam(1) = ishfts 
      iparam(3) = maxitr 
      iparam(7) = mode 
      !
      !
      ! M A I N   L O O P (Reverse communication) 
      !
      !if (verbose>=4) write(out,"(/'Arpack: eigenvalues computed:')") 
      !
      if (verbose>=3) write(out,"(/'Arpack: n = ',i8,' nev = ',i8,' maxitr = ',i8,' tol = ',e12.5)") n,nev,maxitr,tol
      !
      iter = 0 
      !
      do while(ido<=1)
        !
        iter = iter + 1
        !
        if (verbose>=3.and.mod(iter,10)==0) write(out,"(' iter  = ',i8)") iter 
        !
        ! Repeatedly call the routine DSAUPD and take 
        ! actions indicated by parameter IDO until    
        ! either convergence is indicated or maxitr   
        ! has been exceeded.                          
        !

        !dec$ if (blacs_ > 0.or.mpi_ > 0)
           !
           call pdsaupd ( comm, ido, bmat, nloc, which, nev, tol, resid, &
                          ncv, v, ldv, iparam, ipntr, workd, workl, &
                          lworkl, info )
        !
        !dec$ else 
            !
           call dsaupd ( ido, bmat, n, which, nev, tol, resid, &
                          ncv, v, ldv, iparam, ipntr, workd, workl, &
                          lworkl, info )
            !
        !dec$ end if
        !
        !if (verbose>=4.and.iparam(5)>0) then 
        !  !
        !  write(out,"(i8)") iparam(5)
        !  !
        !endif
        !
        if (ido==-1 .or. ido==1) then 
          !
          ! Perform matrix vector multiplication 
          !        y <--- OP*x 
          !
          call matvec_p(comm,n,nloc,nprocs,kstart(0:nprocs),dx,bterm,h,mv_buf,workd(ipntr(1)),workd(ipntr(2)))
          !
        endif 
        !
      enddo
      !
      !  Either we have convergence or there is  an error. 
      !
      if ( info < 0 ) then
         write(out,"(/'Error with _saupd, info = ',i8)") info
         stop 'Error with _saupd'
      endif
      !
      !
      !  No fatal errors occurred. Post-Process using DSEUPD. 
      !  Computed eigenvalues may be extracted. 
      !
      !  Eigenvectors may also be computed now if desired.  (indicated by rvec = .true.)  
      !             
      rvec = .true.
      !
      !dec$ if (blacs_ > 0)
        !
        if (verbose>=5) write(out,"(/'Arpack: dseupd')") 
        !
        call pdseupd ( comm, rvec, 'All', select, d, v, ldv, sigma, &
                       bmat, n, which, nev, tol, resid, ncv, v, ldv, &
                       iparam, ipntr, workd, workl, lworkl, ierr )


        if (verbose>=5) write(out,"(/'Arpack: done!')") 

        !                    
      !dec$ else 
          !
          write(out,"(/'Arpack was not activated yet. Please uncomment dsaupd and  dseupd bellow')")
          stop 'Arpack was not activated'
          !
      !dec$ end if
      !
      if ( ierr < 0 ) then
         write(out,"(/'Error with_seupd, info = ',i8)") ierr
         stop 'Error with _saupd'
      endif
      !
      nroots =  iparam(5)
      !
      e(1:nroots) = d(1:nroots,1)
      !
      !$omp parallel do private(j) shared(h) schedule(dynamic)
      do j=1,n
         h(j,1:nroots) = v(j,1:nroots)
      enddo
      !$omp end parallel do
      ! 
      !  Eigenvalues are returned in the first column 
      !  of the two dimensional array D and the       
      !  corresponding eigenvectors are returned in   
      !  the first NEV columns of the two dimensional 
      !  array V if requested.  Otherwise, an         
      !  orthogonal basis for the invariant subspace  
      !  corresponding to the eigenvalues in D is     
      !  returned in V.                               
      !
      !do j=1, nconv
        !
        ! Compute the residual norm ||  A*x - lambda*x ||  
        ! for the NCONV accurately computed eigenvalues and 
        ! eigenvectors.  (iparam(5) indicates how many are   
        ! accurate to the requested tolerance)               
        !
        !call matvec(n,bterm,v(1,j), ax)
        !
        !call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
        !d(j,2) = dnrm2(n, ax, 1)
        !d(j,2) = d(j,2) / abs(d(j,1))
        !
      !enddo     

!
!     %---------------------------%
!     | Done with program pdsdrv1.|
!     %---------------------------%
!
 9000 continue
!
      !dec$ if (blacs_ > 0)
        call BLACS_GRIDEXIT ( comm )
        call BLACS_EXIT(0)
      !dec$ end if
      !dec$ if (mpi_ > 0)
         call MPI_FINALIZE(rc)
      !dec$ end if


      deallocate(v,workl,workd,d,resid,select,mv_buf)
      !
      call ArrayStop('diag_dseupd_p')
      !
      if (verbose>=2) call TimerStop('diag_dseupd_p: diagonalization')
      !
  end subroutine diag_dseupd_p



  subroutine dseupd_omp_arpack(n,bterm,nroots,factor,maxitr_,tol,h,e)

    integer,intent(in)     :: n,bterm(n,2),maxitr_
    double precision,intent(in)    :: factor
    integer,intent(inout)  :: nroots
    double precision,intent(in)    :: tol
    double precision,intent(inout) :: h(:,:)
    double precision,intent(out)   :: e(nroots)

    double precision,allocatable :: v(:,:),workl(:),workd(:),d(:,:),resid(:)
    !
    logical,allocatable   ::   select(:)
    integer(ik) ::   iparam(11), ipntr(11),ldv,iter
    integer(hik) ::  msize
    !

    character(len=1) ::  bmat
    character(len=2) ::  which
    !
    integer(ik)   :: ido,  nev, ncv, lworkl, info, ierr, j, &
                     mode, ishfts, alloc, maxitr,k
                     !
    logical          rvec
    double precision  ::  sigma
    !
    !double precision,external :: dnrm2
    !
      integer  logfil, ndigit, mgetv0,&
              msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,&
              mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,&
              mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ &
              logfil, ndigit, mgetv0,&
              msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,&
              mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,&
              mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd


      if (verbose>=2) call TimerStart('dseupd_omp_arpack: diagonalization')
      !
      ndigit = -3
      logfil = out
      msgets = 0
      msaitr = 0
      msapps = 0
      msaupd = 0
      msaup2 = 0
      mseigt = 0
      mseupd = 0
      !
      if (verbose>=3) msaupd = 1
      !
      nev = nroots
      !
      if (nroots==0) nev=max(100,n)
      !
      ncv = factor*nev
      !
      if (ncv==nev) ncv = int(21*nev/10)+1
      !
      ncv = min(n,ncv)
      !
      bmat = 'I'
      which = 'SM'
      !
      ! The work array WORKL is used in DSAUPD as workspace. 
      ! The parameter TOL determines the stopping criterion.  
      ! If TOL<=0, machine precision is used.  The variable IDO is used for 
      ! reverse communication and is initially set to 0. 
      ! Setting INFO=0 indicates that a random vector is generated in DSAUPD 
      ! to start the Arnoldi  iteration.  
      !
 
      lworkl = ncv*(ncv+10)
      ldv = n
      !
      !tol = tol
      info = 0
      ido = 0
      !
      allocate(v(ldv,ncv), workl(lworkl),workd(3*n),d(ncv,2),resid(n),select(ncv),stat=alloc)
      !
      msize = int(ldv,hik)*int(ncv,hik)
      !
      if (verbose>=3) write(out,"(/'Arpack: msize = ',i0)") msize
      !
      call ArrayStart('dseupd_omp_arpack:v',alloc,1,kind(v),msize)
      call ArrayStart('dseupd_omp_arpack:workl',alloc,size(workl),kind(workl))
      call ArrayStart('dseupd_omp_arpack:workd',alloc,size(workd),kind(workd))
      call ArrayStart('dseupd_omp_arpack:d',alloc,size(d),kind(d))
      call ArrayStart('dseupd_omp_arpack:resid',alloc,size(resid),kind(resid))
      call ArrayStart('dseupd_omp_arpack:select',alloc,size(select),kind(select))
      !
      !
      iparam = 0 
      ipntr = 0 
      !
      ! This routone uses exact shifts with respect to    
      ! the current Hessenberg matrix (IPARAM(1) = 1).    
      ! IPARAM(3) specifies the maximum number of Arnoldi 
      ! iterations allowed.  Mode 1 of DSAUPD is used     
      ! (IPARAM(7) = 1).  
      !
      ishfts = 1
      maxitr = maxitr_*nev
      mode   = 1
      !      
      iparam(1) = ishfts 
      iparam(3) = maxitr 
      iparam(7) = mode 
      !
      !
      ! M A I N   L O O P (Reverse communication) 
      !
      !if (verbose>=4) write(out,"(/'Arpack: eigenvalues computed:')") 
      !
      if (verbose>=3) write(out,"(/'Arpack: n = ',i8,' nev = ',i8,' ncv = ',i8,' maxitr = ',i8,' tol = ',e12.5)") n,nev,ncv,maxitr,tol
      !
      iter = 0 
      !
      do while(ido<=1)
        !
        iter = iter + 1
        !
        if (verbose>=3.and.mod(iter,500)==0) write(out,"(' iter  = ',i8)") iter 
        !
        ! Repeatedly call the routine DSAUPD and take 
        ! actions indicated by parameter IDO until    
        ! either convergence is indicated or maxitr   
        ! has been exceeded.                          
        !

        !dec$ if (omparpack_ > 0)
           !
           call dsaupd_ ( ido, bmat, n, which, nev, tol, resid, &
                         ncv, v, ldv, iparam, ipntr, workd, workl, &
                         lworkl, info )
        !dec$ else 
            !
            write(out,"(/'Arpack was not activated yet. Please uncomment dsaupd and  dseupd bellow')")
            stop 'Arpack was not activated'
            !
        !dec$ end if
        !
        !if (verbose>=4.and.iparam(5)>0) then 
        !  !
        !  write(out,"(i8)") iparam(5)
        !  !
        !endif
        !
        if (ido==-1 .or. ido==1) then 
          !
          ! Perform matrix vector multiplication 
          !        y <--- OP*x 
          !
          !call matvec(n,bterm,h,workd(ipntr(1)), workd(ipntr(2)))
          !
          call matvec_band(n,h,workd(ipntr(1)), workd(ipntr(2)))
          !
        endif 
        !
      enddo
      !
      !  Either we have convergence or there is  an error. 
      !
      if ( info < 0 ) then
         write(out,"(/'Error with _saupd, info = ',i8)") info
         stop 'Error with _saupd'
      endif
      !
      !
      !  No fatal errors occurred. Post-Process using DSEUPD. 
      !  Computed eigenvalues may be extracted. 
      !
      !  Eigenvectors may also be computed now if desired.  (indicated by rvec = .true.)  
      !             
      rvec = .true.
      !
      !dec$ if (omparpack_ > 0)
        !
        if (verbose>=5) write(out,"(/'Arpack: dseupd')") 
        !
        if (verbose>=2) call TimerStart('dseupd_omp_arpack: Eigenvectors')
        !
        call dseupd_ ( rvec, 'All', select, d, v, ldv, sigma, &
                      bmat, n, which, nev, tol, resid, ncv, v, ldv, &
                      iparam, ipntr, workd, workl, lworkl, ierr )
        !
        if (verbose>=2) call TimerStop('dseupd_omp_arpack: Eigenvectors')
        !
        if (verbose>=5) write(out,"(/'Arpack: done!')") 
        !
      !dec$ else 
          !
          write(out,"(/'Arpack was not activated yet. Please uncomment dsaupd and  dseupd bellow')")
          stop 'Arpack was not activated'
          !
      !dec$ end if
      !
      if ( ierr < 0 ) then
         write(out,"(/'Error with_seupd, info = ',i8)") ierr
         stop 'Error with _saupd'
      endif
      !
      k = size(h,dim=1)
      !
      nroots =  iparam(5)
      !
      nroots = min(nroots,k)
      !
      e(1:nroots) = d(1:nroots,1)
      !
      !$omp parallel do private(j) shared(h) schedule(static)
      do j=1,n
         !
         call dcopy(nroots,v(j,1:nroots),1,h(1:nroots,j),1)
         !
         !h(j,1:nroots) = v(j,1:nroots)
         !
      enddo
      !$omp end parallel do
      ! 
      !  Eigenvalues are returned in the first column 
      !  of the two dimensional array D and the       
      !  corresponding eigenvectors are returned in   
      !  the first NEV columns of the two dimensional 
      !  array V if requested.  Otherwise, an         
      !  orthogonal basis for the invariant subspace  
      !  corresponding to the eigenvalues in D is     
      !  returned in V.                               
      !
      !do j=1, nconv
        !
        ! Compute the residual norm ||  A*x - lambda*x ||  
        ! for the NCONV accurately computed eigenvalues and 
        ! eigenvectors.  (iparam(5) indicates how many are   
        ! accurate to the requested tolerance)               
        !
        !call matvec(n,bterm,v(1,j), ax)
        !
        !call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
        !d(j,2) = dnrm2(n, ax, 1)
        !d(j,2) = d(j,2) / abs(d(j,1))
        !
      !enddo     

      deallocate(v,workl,workd,d,resid,select)
      !
      call ArrayStop('dseupd_omp_arpack:v')
      call ArrayStop('dseupd_omp_arpack:workl')
      call ArrayStop('dseupd_omp_arpack:workd')
      call ArrayStop('dseupd_omp_arpack:d')
      call ArrayStop('dseupd_omp_arpack:resid')
      call ArrayStop('dseupd_omp_arpack:select')
      !
      if (verbose>=2) call TimerStop('dseupd_omp_arpack: diagonalization')
      !
  end subroutine dseupd_omp_arpack





  subroutine diag_propack(n,bterm,nroots,factor,maxiter,iverbose,tol,h,e)
    !
    integer,intent(in)     :: n,bterm(n,2),maxiter,iverbose
    double precision,intent(in)    :: factor
    integer,intent(inout)  :: nroots
    double precision,intent(in)    :: tol
    double precision,intent(inout) :: h(:,:)
    double precision,intent(out)   :: e(nroots)
    !
    double precision,allocatable :: u(:,:),v(:,:),B(:,:),work(:),dparm(:),S(:),uv(:,:),energies(:)
    !
    integer(ik) ::   ldu,ldv,ldb,k,k0,ioption(2),nev,ncv,alloc,lwork,info,i,istart,iend,j,dimen,nb,dim,iter,n2
    !
    integer(ik),allocatable ::   iparm(:),iwork(:)
    !
    integer(hik)     ::  matsize
    !
    double precision :: doption(3),rnorm,conv_nroots
    !
    character(len=1) :: which, jobu,jobv
    !
    !external :: daprod


      if (iverbose>=2) call TimerStart('diag_propack: diagonalization')
      !
      !
      !dec$ if (propack_ < 1)
        !
        write(out,'("PROPACK is not activated!")')
        stop 'PROPACK is not activated!'
        e = 0 
        nroots = 0
        return
        !
      !dec$ endif
      !
      nev = nroots
      !
      if (nroots==0) nev=max(100,n)
      !
      ncv = factor*nev
      !
      if (ncv==nev) ncv = int(21*nev/10)+1
      !
      ncv = min(int(n*0.9),ncv)
      !
      ! K: INTEGER. On entry: The desired dimension of the Lanczos 
      !    bidiagonalization.
      !
      k = ncv 
      !
      doption(1) = 1e-32
      doption(2) = 1e-32
      !
      doption(3) = 1
      !
      ioption(1) = 1
      !
      ioption(2) = 1
      !
      matsize = 0 
      do i=1,n
         istart = bterm(i,1)
         iend   = bterm(i,2)
         matsize = matsize + iend-istart+1
      enddo
      !
      matsize = 1
      !
      allocate(dparm(matsize),iparm(n),hikparm(n),energies(n),stat=alloc)
      call ArrayStart('diag_propack:dparm',alloc,1,kind(dparm),matsize)
      call ArrayStart('diag_propack:iparm',alloc,size(iparm),kind(iparm))
      call ArrayStart('diag_propack:hikparm',alloc,size(hikparm),kind(hikparm))
      call ArrayStart('diag_propack:energies',alloc,size(energies),kind(energies))
      !
      matsize = size(h)
      !
      n2 = size(h,dim=2)
      !
      allocate(hmat_(n2,n),stat=alloc)
      call ArrayStart('diag_propack:hmat',alloc,size(energies),kind(energies))
      !
      matsize = 0 
      do i=1,n
         istart = bterm(i,1)
         iend   = bterm(i,2)
         !
         hmat_(:,i) = h(i,:)
         !
         iparm(i) = max(bterm(i,2)-i,i-bterm(i,1))
         !
         hikparm(i) = matsize
         !
         dimen = iend-istart+1
         !
         !do j = 1,dimen
         !  !
         !  matsize = matsize + 1
         !  !
         !  dparm(matsize) = h(i,j)
         !  !
         !enddo
         !
      enddo
      !
      k0 = 0
      k = ncv
      !
      ldu = n ; ldv = n 
      !
      allocate (U(ldu,k+1),V(ldv,k),stat=alloc)
      !
      matsize  = int(ldu*(k+1)*2,hik)
      call ArrayStart('diag_propack:UV',alloc,1,kind(U),matsize)
      !
      U = 0 ; V = 0
      !
      allocate (B(n,2),stat=alloc)
      call ArrayStart('diag_propack:B',alloc,size(B),kind(B))
      !
      e = 0
      !
      B = 0 
      !
      ldb = n
      !
      do iter = 1,maxiter
          !
          lwork = 2.5*(n+n+k+1)
          !
          allocate(work(lwork),iwork(2*k+1),stat=alloc)
          !
          matsize = lwork
          call ArrayStart('diag_propack:work',alloc,1,kind(work),matsize)
          call ArrayStart('diag_propack:work',alloc,size(iwork),kind(iwork))
          !
          rnorm = 0
          !
          b = 0
          !
          u = 0 ; v = 0 
          k0 = 0
          !
          call dlanbpro( n, n, k0, k, daprod, U, ldu, V, ldv, B, ldb,rnorm, doption, ioption, work, iwork, h, iparm, info)
          !
          !  Either we have convergence or there is  an error. 
          !
          if ( info < 0 ) then
             write(out,"(/'Error from dlanbpr:  An invariant subspace of dimension -J was found, info = ',i8)") info
             stop 'Error with dlanbpro'
          endif
          !
          if ( info > 0 ) then
             write(out,"(/'Error from dlanbpr:  The computation succeeded, but the algorithm came close to computing an invariant subspace, info = ',i8)") info
             !stop 'Error with dlanbpro'
          endif
          !
          deallocate(work,iwork)
          call ArrayStop('diag_propack:work')
          !
          dim = k
          !
          !k = nev
          !
          NB = 32
          !
          lwork = (3*dim**2 + max(3*dim**2+4*dim+4,nb*n))*1.2
          !
          allocate(work(lwork),iwork(8*dim),stat=alloc)
          !
          matsize = lwork
          call ArrayStart('diag_propack:work',alloc,1,kind(work),matsize)
          call ArrayStart('diag_propack:work',alloc,size(iwork),kind(iwork))
          !
          !allocate (U(ldu,dim+1),V(ldu,dim),stat=alloc)
          !
          !matsize  = int(ldu*(2*dim+2),hik)
          !call ArrayStart('diag_propack:UV',alloc,1,kind(U),matsize)
          !
          which = 's'
          jobu = 'y'
          jobv = 'y'
          !
          allocate (S(k),stat=alloc)
          call ArrayStart('diag_propack:S',alloc,size(S),kind(S))
          !
          S = 0
          !
          call dritzvec(which,jobu,jobv,n,n,nev,dim,B(1:dim,1),B(1:dim,2),S,U,ldu,V,ldv,work,lwork,iwork)
          !
          nroots = min(nroots,nev)
          !
          !do i = 1,nroots
          !  e(i) = b(dim-i+1,1)
          !enddo
          !
          energies(1:dim) = b(dim:1:-1,1)
          !
          conv_nroots = abs( e(nroots) - energies(nroots) )
          !
          if (iverbose>=5) then 
            !
            write(out,'(/"Energies (1/cm): ")')
            !
            do i = 1,dim
              write(out,'(f18.8)') energies(i)
            enddo
            !
            write(out,'("Convergence of the ",i0,"th root = ",g18.8)') nroots,conv_nroots
            !
          endif 
          !
          e(1:nroots) = energies(1:nroots)
          !
          k0 = k
          !
          k = k + ncv
          !
          if (k<n-1.and.iter+1/=maxiter.and.conv_nroots>=tol) then 
            !
            !allocate (UV(ldu,k0+1),stat=alloc)
            !
            !matsize  = int(ldu*(k0+1),hik)
            !call ArrayStart('diag_propack:UV',alloc,1,kind(UV),matsize)
            !
            !UV = U
            !
            deallocate(U)
            allocate (U(ldu,k+1),stat=alloc)
            matsize  = int(ldu*(k+1),hik)
            call ArrayStart('diag_propack:UV',alloc,1,kind(U),matsize)
            !
            !U(:,1:k0+1) = UV(:,1:k0+1)
            !
            !UV(:,1:k0) = V(:,1:k0)
            !
            deallocate(V)
            allocate (V(ldv,k),stat=alloc)
            matsize  = int(ldv*(k+1),hik)
            call ArrayStart('diag_propack:UV',alloc,1,kind(V),matsize)
            !
            !V(:,1:k0) = UV(:,1:k0)
            !
            !deallocate (UV)
            !call ArrayStop('diag_propack:UV')
            !
            !matsize  = int(ldu*(k+1)*2,hik)
            !call ArrayStart('diag_propack:UV',alloc,1,kind(U),matsize)
            !
          endif
          !
          deallocate(S)
          call ArrayStop('diag_propack:S')
          !
          deallocate(work,iwork)
          call ArrayStop('diag_propack:work')
          !
          if (k>=n.or.conv_nroots<tol) exit
          !
      enddo
      !
      if (conv_nroots>=tol.and.k<n) then 
        write(out,'("Convergence is not reached after ",i0," iterations, ",g19.8)') iter,conv_nroots
        stop 'diag_propack: No convergence'
      endif
      !
      !omp parallel do private(i) shared(h) schedule(dynamic)
      do i = 1,nroots
          h(:,i) = v(:,dim-i+1)
          !
          write(out,'(i4,10g18.8)') i,v(1:min(10,n),dim-i+1)
          !
      enddo
      !omp end parallel do
      !
      !
      !
      deallocate(b)
      call ArrayStop('diag_propack:B')
      !
      deallocate(u,v)
      call ArrayStop('diag_propack:UV')
      !
      !e(1:nroots) = b(n:n-nroots+1:-1,1)
      !
      !omp parallel do private(j) shared(h) schedule(static)
      !do j=1,n
      !   !
      !   call dcopy(nroots,v(j,1:nroots),1,h(1:nroots,j),1)
      !   !
      !enddo
      !omp end parallel do
      !
      deallocate(hmat_)
      call ArrayStop('diag_propack:hmat')
      !
      deallocate(dparm,iparm,energies,hikparm)
      call ArrayStop('diag_propack:dparm')
      call ArrayStop('diag_propack:iparm')
      call ArrayStop('diag_propack:energies')
      call ArrayStop('diag_propack:hikparm')
      !
      if (iverbose>=2) call TimerStop('diag_propack: diagonalization')
      !
  end subroutine diag_propack

  subroutine daprod1(transa,m,n,x,y,dparm,iparm)
    !
      character(len=1) :: transa
      integer(ik)      ::  m,n,iparm(*),bterm(n)
      double precision :: x(*),y(*),dparm(*),f_t
      !
      double precision,parameter :: alpha = 1.0d0,beta=0.0d0
      integer                   :: lda,k,iend,istart,dimen
      integer(hik)              :: i_t
      double precision,external :: ddot
      !
      f_t = dparm(1)
      i_t = iparm(1)
      !
      !$omp parallel do private(k,istart,iend,dimen,i_t) shared(y) schedule(dynamic)
      do k=1,n
         !
         istart = max(1,k-iparm(k))
         iend   = min(n,k+iparm(k))
         !
         dimen = iend-istart+1
         !
         i_t = hikparm(k)
         !
         y(k) = ddot(dimen,dparm(i_t+1:i_t+dimen),1,x(istart:iend),1)
         !
         !i_t = i_t + dimen
         !
      enddo 
      !$omp end parallel do 
      !
      !k = size(h,dim=1)-1
      !lda = k+1
      !
      !call dsbmv('L',n,k,alpha,h,lda,x,1,beta,y,1)
    !
  end subroutine daprod1


  subroutine daprod(transa,m,n,x,y,dparm,iparm)
    !
      character(len=1) :: transa
      integer(ik)      ::  m,n,iparm(*),bterm(n)
      double precision :: x(*),y(*),dparm(*),f_t
      !
      double precision,parameter :: alpha = 1.0d0,beta=0.0d0
      integer                   :: lda,k,iend,istart,dimen
      integer(hik)              :: i_t
      double precision,external :: ddot
      !
      !f_t = dparm(1)
      !i_t = iparm(1)
      !
      !$omp parallel do private(k,istart,iend,dimen,i_t) shared(y) schedule(dynamic)
      do k=1,n
         !
         istart = max(1,k-iparm(k))
         iend   = min(n,k+iparm(k))
         !
         dimen = iend-istart+1
         !
         !i_t = hikparm(k)
         !
         !y(k) = ddot(dimen,dparm(i_t+1:i_t+dimen),1,x(istart:iend),1)
         !
         y(k) = ddot(dimen,hmat_(1:dimen,k),1,x(istart:iend),1)
         !
         !i_t = i_t + dimen
         !
      enddo 
      !$omp end parallel do 
      !
      !k = size(h,dim=1)-1
      !lda = k+1
      !
      !call dsbmv('L',n,k,alpha,h,lda,x,1,beta,y,1)
    !
  end subroutine daprod 


end module diag

