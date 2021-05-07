module lapack
 
#define arpack_  0
#define blacs_   0
#define mpi_     0
#define omparpack_  0

!
!  Simplistic type-agnostic LAPACK interface
!
  use accuracy
  use timer
  implicit none
  private verbose

  interface lapack_gelss
    module procedure lapack_cgelss
    module procedure lapack_zgelss
    module procedure lapack_sgelss
    module procedure lapack_dgelss
  end interface ! lapack_gelss

  interface lapack_stev
    module procedure lapack_sstev
    module procedure lapack_dstev
  end interface ! lapack_stev

  interface lapack_heev
    module procedure lapack_cheev
    module procedure lapack_zheev
  end interface ! lapack_heev

  interface lapack_syev
    module procedure lapack_dsyev
    module procedure lapack_ssyev
  end interface ! lapack_syev

  interface lapack_syevd
    module procedure lapack_dsyevd
    module procedure lapack_ssyevd
  end interface ! lapack_syev

  interface lapack_syevr
    module procedure lapack_dsyevr
    module procedure lapack_ssyevr
  end interface ! lapack_syevr

  interface lapack_gesvd
    module procedure lapack_dgesvd
    module procedure lapack_sgesvd
  end interface ! lapack_gesvd

  interface lapack_syevx
    module procedure lapack_dsyevx
  end interface ! lapack_syevx

  interface lapack_ginverse
    module procedure lapack_ginverse_real
    module procedure lapack_ginverse_double
  end interface ! lapack_ginverse

  integer,parameter:: verbose = 3
  real(rk),parameter :: singtol = -1.0d-12 ! 100.0d0*spacing(1.0d0)
  !
  contains

  subroutine lapack_cgelss(a,b)
    complex, intent(inout) :: a(:,:)
    complex, intent(inout) :: b(:,:)

    external cgelss
    real                   :: s    (   min(size(a,dim=1),size(a,dim=2)))
    complex                :: work (50*max(size(a,dim=1),size(a,dim=2)))
    real                   :: rwork( 5*min(size(a,dim=1),size(a,dim=2)))
    integer                :: rank, info
    integer                :: na1, na2, nb1, nb2
    
    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call cgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s, 100.0*spacing(1.0), rank, work, size(work), rwork, info)

    if (info/=0) then
      write (out,"(' cgelss returned ',i8)") info
      stop 'lapack_cgelss - cgelss failed'
    end if
  end subroutine lapack_cgelss

  subroutine lapack_zgelss(a,b)
    double complex, intent(inout) :: a(:,:)
    double complex, intent(inout) :: b(:,:)

    external zgelss
    double precision       :: s    (   min(size(a,dim=1),size(a,dim=2)))
    double complex         :: work (50*max(size(a,dim=1),size(a,dim=2)))
    double precision       :: rwork( 5*min(size(a,dim=1),size(a,dim=2)))
    integer                :: rank, info
    integer                :: na1, na2, nb1, nb2
    
    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call zgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s,singtol, rank, work, size(work), rwork, info)

    if (info/=0) then
      write (out,"(' cgelss returned ',i8)") info
      stop 'lapack_cgelss - cgelss failed'
    end if
  end subroutine lapack_zgelss

  subroutine lapack_sgelss(a,b)
    real, intent(inout) :: a(:,:)
    real, intent(inout) :: b(:,:)

    external sgelss
    real                :: s    (   min(size(a,dim=1),size(a,dim=2)))
    real                :: work (50*max(size(a,dim=1),size(a,dim=2),size(b,dim=2)))
    integer             :: rank, info
    integer             :: na1, na2, nb1, nb2
    
    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call sgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s, 100.0*spacing(1.0), rank, work, size(work), info)

    if (info/=0) then
      write (out,"(' sgelss returned ',i8)") info
      stop 'lapack_sgelss - sgelss failed'
    end if
  end subroutine lapack_sgelss

  subroutine lapack_dgelss(a,b,ierror)
    double precision, intent(inout) :: a(:,:)
    double precision, intent(inout) :: b(:,:)
    integer,intent(inout),optional  :: ierror

    external dgelss
    double precision    :: s    (   min(size(a,dim=1),size(a,dim=2))),tol
    double precision    :: work (70*max(size(a,dim=1),size(a,dim=2),size(b,dim=2)))
    double precision,allocatable  :: work_t (:),h(:,:)
    integer(ik)          :: rank,info1,info2,info
    integer(ik)          :: na1, na2, nb1, nb2, iw

   
    if (verbose>=4) write(out,"('lapack_dgelss...')") 
    !
    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    !
    if (verbose>=5) write(out,"('lapack_dgelss...na1,nb1,na2,nb2 = ',4i8)") na1,nb1,na2,nb2
    !
    allocate(h(na1,na2),stat=info1)
    !
    if (verbose>=5) write(out,"('lapack_dgelss...allocated')") 
    !
    if (info1/=0) then
      !
      write (out,"(' allocation (h) error ',i8,' size = ',2i8)") info1,na1,na2
      stop 'lapack_dgelss - allocation of h failed'
      !
    end if
    !
    if (verbose>=5) write(out,"('lapack_dgelss: a=>h ')")
    !
    h = a
    !
     if (verbose>=5) write(out,"('lapack_dgelss: a=>h ... done')")
    !
    tol = singtol
    !
    call dgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s,singtol, rank, work, -1, info)
    !
    iw = int(work(1))
    !
    if (verbose>=5) write(out,"('  size(work) = ',i8)") iw 
    !
    allocate(work_t(iw),stat=info2)
    !
    if (verbose>=5) write(out,"('  work_t allocated? ')")
    !
    if (info2/=0) then
      !
      write (out,"(' allocation (of work) error ',i8,' iwork = ',i8)") info2,iw
      stop 'lapack_dgelss - allocation of work_t failed'
      !
    end if
    !
    call dgelss(na1,na2,nb2,h(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s,singtol, rank, work_t, iw, info)
    !
    a = h 
    !
    deallocate(work_t,h)
    !
    if (info/=0.or.info1/=0.or.info2/=0) then
      !
      write (out,"(' dgelss returned ',i8)") info
      write (out,"(' dgelss: na1 = ',i8,' nb1 =  ',i8,' na2 = ',i8,' nb2 =  ',i8,' iwork = ',i8)") na1,nb1,na2,nb2,iw
      !
      if (present(ierror).and.info/=0) then 
        !
        ierror = info
        !
      else
        !
        stop 'lapack_dgelss - dgelss failed'
        !
      endif
      !
    end if
    !
    if (rank/=na2.and.present(ierror)) then 
      !
      ierror = 1
      ! 
    endif
    !
    if (verbose>=4) write(out,"('...lapack_dgelss done!')") 
    !
  end subroutine lapack_dgelss


  subroutine lapack_sstev(d,e,z)
    real, intent(inout) :: d(:)   ! In:  Diagonal elements of the matrix 
                                  ! Out: Eigenvalues, ascending order
    real, intent(inout) :: e(:)   ! In:  Sub-/super-diagonal elements of the matrix
                                  ! Out: Destroyed
    real, intent(out)   :: z(:,:) ! Out: Eigenvectors

    real    :: work(max(1,2*size(d)-2))
    integer :: info
    integer :: nz1, nz2
    
    nz1 = size(z,dim=1) ; nz2 = size(z,dim=2)
    call sstev('V',size(d),d,e,z(1:nz1,1:nz2),nz1,work,info)

    if (info/=0) then
      write (out,"(' sstev returned ',i8)") info
      stop 'lapack_sstev - sstev failed'
    end if
  end subroutine lapack_sstev

  subroutine lapack_dstev(d,e,z)
    double precision, intent(inout) :: d(:)   ! In:  Diagonal elements of the matrix 
                                              ! Out: Eigenvalues, ascending order
    double precision, intent(inout) :: e(:)   ! In:  Sub-/super-diagonal elements of the matrix
                                              ! Out: Destroyed
    double precision, intent(out)   :: z(:,:) ! Out: Eigenvectors

    double precision :: work(max(1,2*size(d)-2))
    integer          :: info
    integer          :: nz1, nz2
    
    nz1 = size(z,dim=1) ; nz2 = size(z,dim=2)
    call dstev('V',size(d),d,e,z(1:nz1,1:nz2),nz1,work,info)

    if (info/=0) then
      write (out,"(' dstev returned ',i8)") info
      stop 'lapack_dstev - dstev failed'
    end if
  end subroutine lapack_dstev

  subroutine lapack_cheev(h,e)
    complex, intent(inout) :: h(:,:)  ! In:  Hermitian matrix to be diagonalized
                                      ! Out: Eigenvectors
    real, intent(out)   :: e(:)       ! Out: Eigenvalues

    complex :: work(50*size(h,dim=1))
    real    :: rwork(3*size(h,dim=1))
    integer :: info
    integer :: nh1, nh2
    
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call cheev('V','L',nh1,h(1:nh1,1:nh2),nh1,e(:),work,size(work),rwork,info)
    if (info/=0) then
      write (out,"(' cheev returned ',i8)") info
      stop 'lapack_cheev - cheev failed'
    end if
  end subroutine lapack_cheev

  subroutine lapack_zheev(h,e)
    double complex, intent(inout) :: h(:,:)  ! In:  Hermitian matrix to be diagonalized
                                             ! Out: Eigenvectors
    double precision, intent(out)   :: e(:)  ! Out: Eigenvalues

    double complex   :: work(50*size(h,dim=1))
    double precision :: rwork(3*size(h,dim=1))
    integer          :: info
    integer          :: nh1, nh2
    
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call zheev('V','L',nh1,h(1:nh1,1:nh2),nh1,e(:),work,size(work),rwork,info)
    if (info/=0) then
      write (out,"(' zheev returned ',i8)") info
      stop 'lapack_zheev - zheev failed'
    end if
  end subroutine lapack_zheev


  subroutine lapack_dgesvd(h)
    double precision, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    character(len=1) :: jobu,jobvt

    double precision,allocatable    :: work(:),u(:,:),vt(:,:),s(:)
    integer          :: info
    integer          :: nh1, nh2
    integer          :: lwork
    !
    if (verbose>=2) call TimerStart('lapack_dgesvd: diagonalization')
    !
    jobu = 'O'
    !
    jobvt = 'N'
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    !
    lwork = 50*nh1
    !
    allocate(work(lwork),u(1,nh1),vt(1,nh2),s(min(nh1,nh2)),stat=info)
    !
    call ArrayStart('lapack_gesvd-arrays-work',info,size(u),kind(u))
    call ArrayStart('lapack_gesvd-arrays-work',info,size(vt),kind(vt))
    call ArrayStart('lapack_gesvd-arrays-work',info,size(s),kind(s))
    !
    call dgesvd(jobu,jobvt,nh1,nh2,h(1:nh1,1:nh2),nh1,s,u,nh1,vt,nh2,work,-1,info)
    !
    if (int(work(1))>size(work)) then 
      !
      lwork = int(work(1))
      !
      deallocate(work)
      !
      allocate(work(lwork))
      !
    endif 
    !
    call dgesvd(jobu,jobvt,nh1,nh2,h(1:nh1,1:nh2),nh1,s,u,nh1,vt,nh2,work,lwork,info)
    !
    deallocate(work,u,vt,s)
    !
    call ArrayStop('lapack_gesvd-arrays-work')
    !
    if (info/=0) then
      write (out,"(' dgesvd returned ',i8)") info
      stop 'lapack_dgesvd - dgesvd failed'
    end if
    !
    if (verbose>=2) call TimerStop('lapack_dgesvd: diagonalization')
    !
  end subroutine lapack_dgesvd



  subroutine lapack_sgesvd(h)
    real, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    character(len=1) :: jobu,jobvt

    real,allocatable    :: work(:),u(:,:),vt(:,:),s(:)
    integer          :: info
    integer          :: nh1, nh2
    integer          :: lwork
    !
    if (verbose>=2) call TimerStart('lapack_sgesvd: diagonalization')
    !
    jobu = 'O'
    !
    jobvt = 'N'
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    !
    lwork = 50*nh1
    !
    allocate(work(lwork),u(nh1,nh1),vt(nh2,nh2),s(min(nh1,nh2)),stat=info)
    !
    call ArrayStart('lapack_gesvd-arrays-work',info,size(u),kind(u))
    call ArrayStart('lapack_gesvd-arrays-work',info,size(vt),kind(vt))
    call ArrayStart('lapack_gesvd-arrays-work',info,size(s),kind(s))
    !
    call sgesvd(jobu,jobvt,nh1,nh2,h(1:nh1,1:nh2),nh1,s,u,nh1,vt,nh2,work,-1,info)
    !
    if (int(work(1))>size(work)) then 
      !
      lwork = int(work(1))
      !
      deallocate(work)
      !
      allocate(work(lwork))
      !
    endif 
    !
    call sgesvd(jobu,jobvt,nh1,nh2,h(1:nh1,1:nh2),nh1,s,u,nh1,vt,nh2,work,lwork,info)
    !
    deallocate(work,u,vt,s)
    !
    call ArrayStop('lapack_gesvd-arrays-work')
    !
    if (info/=0) then
      write (out,"(' dgesvd returned ',i8)") info
      stop 'lapack_dgesvd - sgesvd failed'
    end if
    !
    !
    if (verbose>=2) call TimerStop('lapack_sgesvd: diagonalization')
    !
  end subroutine lapack_sgesvd




  subroutine lapack_dsyev(h,e,jobz)
    double precision, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    double precision, intent(out)   :: e(:)    ! Out: Eigenvalues
    character(len=1),intent(in),optional :: jobz
    character(len=1) :: jobz_

    double precision,allocatable    :: work(:)
    integer(ik)      :: info
    integer(ik)      :: nh1, nh2
    integer(ik)      :: lwork
    !
    !if (verbose>=2) call TimerStart('lapack_dsyev: diagonalization')
    !
    jobz_ = 'V'
    !
    if (present(jobz)) then
       jobz_ = jobz
    endif
    
    lwork = 50*size(h,dim=1)
    !
    allocate(work(lwork))
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call dsyev(jobz_,'L',nh1,h(1:nh1,1:nh2),nh1,e,work,-1,info)
    !
    if (int(work(1))>size(work)) then 
      !
      lwork = int(work(1))
      !
      deallocate(work)
      !
      allocate(work(lwork))
      !
      call dsyev(jobz_,'L',nh1,h(1:nh1,1:nh2),nh1,e,work,lwork,info)
      !
    else
      !
      call dsyev(jobz_,'L',nh1,h(1:nh1,1:nh2),nh1,e,work,size(work),info)
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
  end subroutine lapack_dsyev


  subroutine lapack_ssyev(h,e,jobz)
    real, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                   ! Out: Eigenvectors
    real, intent(out)   :: e(:)    ! Out: Eigenvalues
    character(len=1),intent(in),optional :: jobz
    character(len=1) :: jobz_

    real,allocatable    :: work(:)
    integer          :: info
    integer          :: nh1, nh2,i
    integer          :: lwork
    !
    if (verbose>=2) call TimerStart('lapack_ssyev: diagonalization')
    !
    jobz_ = 'V'
    !
    if (present(jobz)) then
       jobz_ = jobz
    endif
    !
    lwork = 50*size(h,dim=1)
    !
    allocate(work(lwork))
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call ssyev(jobz_,'L',nh1,h(1:nh1,1:nh2),nh1,e,work,-1,info)
    !
    if (int(work(1))>size(work)) then 
      !
      lwork = int(work(1))
      !
      deallocate(work)
      !
      allocate(work(lwork))
      !
      call ssyev(jobz_,'L',nh1,h(1:nh1,1:nh2),nh1,e,work,lwork,info)
      !
    else
      !
      call ssyev(jobz_,'L',nh1,h(1:nh1,1:nh2),nh1,e,work,size(work),info)
      !
    endif 
    !
    deallocate(work)
    !
    if (jobz == 'N') then  
      h = 0 
      forall(i=1:min(nh1,nh2)) h(i,i) = 1.0d0
    endif
    !
    if (info/=0) then
      write (out,"(' ssyev returned ',i8)") info
      stop 'lapack_ssyev - dsyev failed'
    end if
    !
    if (verbose>=2) call TimerStop('lapack_ssyev: diagonalization')
    !
  end subroutine lapack_ssyev



  subroutine lapack_dsyevd(h,e)
    double precision, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    double precision, intent(out)   :: e(:)    ! Out: Eigenvalues

    integer          :: info
    integer          :: nh1, nh2,liwork,lwork
    double precision,allocatable :: work(:)
    integer,allocatable :: iwork(:)
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    !
    lwork  = 2*nh1**2+6*nh1+10
    liwork = 3*nh1+10
    !
    allocate(work(lwork),stat=info)
    call ArrayStart('lapack_ssyevd-arrays-work',info,size(work),kind(work))
    allocate(iwork(liwork),stat=info)
    call ArrayStart('lapack_ssyevd-arrays-work',info,size(iwork),kind(iwork))
    !
    if (info/=0) then
      write (out,"(' dsyevd returned allocation work failed',i8)") info
      stop 'lapack_dsyevd - allocation work failed'
    end if
    call dsyevd('V','L',nh1,h(1:nh1,1:nh2),nh1,e,work,-1,iwork,-1,info)
    !
    if (int(work(1))>size(work).or.int(iwork(1))>size(iwork)) then 
      !
      lwork = int(work(1))
      liwork = int(iwork(1))
      !
      deallocate(work,iwork)
      !
      allocate(work(lwork),iwork(liwork))
      !
      call dsyevd('V','L',nh1,h(1:nh1,1:nh2),nh1,e,work,lwork,iwork,liwork,info)
      !
    else
      !
      call dsyevd('V','L',nh1,h(1:nh1,1:nh2),nh1,e,work,lwork,iwork,liwork,info)
      !
    endif 
    !
    if (info/=0) then
      write (out,"(' dsyevd returned ',i8)") info
      stop 'lapack_dsyevd - dsyev failed'
    end if
    !
    deallocate(work,iwork)
    call ArrayStop('lapack_ssyevd-arrays-work')
    !
  end subroutine lapack_dsyevd



  subroutine lapack_ssyevd(h,e)
    real, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                   ! Out: Eigenvectors
    real, intent(out)   :: e(:)    ! Out: Eigenvalues

    real     :: work(10*size(h,dim=1))
    integer  :: iwork(10*size(h,dim=1))
    integer  :: info
    integer  :: nh1, nh2,liwork,lwork
    real,allocatable :: twork(:)
    integer,allocatable :: tiwork(:)
    
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    !
    lwork = size(work)
    iwork = size(iwork)
    !
    call ssyevd('V','L',nh1,h(1:nh1,1:nh2),nh1,e,work,-1,iwork,-1,info)
    !
    lwork = int(work(1))
    liwork = int(iwork(1))
    !
    allocate(twork(lwork),tiwork(liwork))
    !
    call ssyevd('V','L',nh1,h(1:nh1,1:nh2),nh1,e,work,lwork,iwork,liwork,info)
    !
    if (info/=0) then
      write (out,"(' ssyev returned ',i8)") info
      stop 'lapack_ssyevd - ssyevd failed'
    end if
    !
    deallocate(twork,tiwork)
    !
  end subroutine lapack_ssyevd


  subroutine lapack_ssyevr(h,e,rng,jobz,iroots,vrange,irange)
    real, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    real, intent(out)   :: e(:)    ! Out: Eigenvalues
    character(len=1),intent(in),optional   :: rng,jobz
    real, intent(in),optional  :: vrange(2)  ! In:  energy range to be computed
    integer         , intent(in),optional  :: irange(2)  ! In:  index to be computed
    integer         , intent(out),optional :: iroots     ! Out: Number of roots found
    !
    real :: vl,vu,abstol
    integer          :: info
    integer          :: nh1, nh2, il, ir,n,k, alloc
    character(len=1) :: rng_,jobz_
    integer          :: ib,jb,niwork,nwork,ldz,msize
    real,allocatable :: a(:,:),work(:),twork(:)
    integer,allocatable :: iwork(:),isuppz(:),tiwork(:)
    !
    if (verbose>=2) call TimerStart('lapack_ssyevr: diagonalization')
    !
    rng_ = 'A'
    jobz_ = 'V'
    il = 1   ; ir = 100000
    vl = -.001 ; vu = 100000.0
    !
    if (present(rng)) then
       rng_ = rng
    endif
    !
    if (present(jobz)) then
       jobz_ = jobz
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
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    !
    msize = nh2
    if (rng_=='I') msize = ir-il+1
    !
    if (nh1/=size(e)) then
      write (out,"(' ssyevr inconsistent sizes of h and e :',2i8)") nh1/=size(e)
      stop 'lapack_ssyevr inconsistent sizes of h and e'
    end if
    !
    abstol =(small_)
    !
    nwork = 50*nh1
    niwork = 20*nh1
    ldz = 2*nh1
    !
    if (verbose>=4) then 
      write(out,"(//'MATRIX:')") 
      do ib=1,nh1
        do jb=ib,nh1
          write(out,"(2i8,f19.8)") ib,jb,h(jb,ib)
        enddo
      enddo
    endif 
    !
    allocate(work(nwork),stat=alloc)
    call ArrayStart('lapack_ssyevr-arrays-work',alloc,size(work),kind(work))
    allocate(iwork(niwork),stat=alloc)
    call ArrayStart('lapack_ssyevr-arrays',alloc,size(iwork),kind(iwork))
    allocate(isuppz(2*msize),stat=alloc)
    call ArrayStart('lapack_ssyevr-arrays',alloc,size(isuppz),kind(isuppz))
    !
    !if (verbose>=3) call MemoryReport
    !if (verbose>=4) write(out,"(/'range = ',a)") rng_
    !
    if (verbose>=4) write(out,"(/'abstol = ',e18.8)") abstol
    call ssyevr(jobz_,rng_,'L',nh1,h(1:nh1,1:nh2),nh1,vl,vu,il,ir,abstol,n,e(1:nh1),h(1:nh1,1:nh2),&
                                                          nh1,isuppz,work,-1,iwork,-1,info)
    !
    if (verbose>=4) write(out,"(/'nwork = ',i8,' niwork = ',i8,' changed to ',2i8)") nwork,niwork,int(work(1)),int(iwork(1))
    !
    nwork = int(work(1))
    niwork = int(iwork(1))
    !
    !if (verbose>=3) write(out,"(/'nwork = ',i8,' niwork = ',i8)") nwork,niwork
    !
    deallocate(work,iwork)
    call ArrayStop('lapack_ssyevr-arrays')
    call ArrayStop('lapack_ssyevr-arrays-work')
    !
    allocate(twork(nwork),stat=alloc)
    call ArrayStart('lapack_ssyevr-arrays-twork',alloc,size(twork),kind(twork))
    allocate(tiwork(niwork),stat=alloc)
    call ArrayStart('lapack_ssyevr-arrays',alloc,size(tiwork),kind(tiwork))
    !
    allocate (a(nh1,msize),stat=alloc)
    call ArrayStart('lapack_ssyevr-a',alloc,size(a),kind(a))
    !
    if (verbose>=3) call MemoryReport
    !
    call ssyevr(jobz_,rng_,'L',nh1,h(1:nh1,1:nh2),nh1,vl,vu,il,ir,abstol,n,e(1:nh1),a(1:nh1,1:msize),&
                nh1,isuppz,twork(1:nwork),nwork,tiwork(1:niwork),niwork,info)
    !
    if (verbose>=3) call TimerStop('lapack_ssyevr: diagonalization finished')
    !
    !$omp parallel do private(k) shared(h) schedule(dynamic)
    do k=1,msize
       h(:,k) = a(:,k)
    enddo
    !$omp end parallel do

    !h(1:nh1,1:nh2) = a(1:nh1,1:nh2)
    !
    deallocate(a,twork,tiwork,isuppz)
    call ArrayStop('lapack_ssyevr-a')
    call ArrayStop('lapack_ssyevr-arrays')
    call ArrayStop('lapack_ssyevr-arrays-twork')
    !
    if (verbose>=5) then 
      write(out,"(//'ssyevr solution:')") 
      do ib=1,n
         write(out,"(i8,f19.8)") ib,e(ib)
         do jb=1,nh1
           write(out,"(2i8,f19.8)") ib,jb,h(jb,ib)
         enddo
      enddo
      !
      call MemoryReport
      !
    endif 
    if (present(iroots))  iroots = n
    !
    if (info/=0) then
      write (out,"(' ssyevr returned ',i8)") info
      stop 'lapack_ssyevr - ssyevr failed'
    end if
    !
    if (verbose>=2) call TimerStop('lapack_ssyevr: diagonalization')
    !
  end subroutine lapack_ssyevr


  subroutine lapack_dsyevr(h,e,rng,jobz,iroots,vrange,irange)
    double precision, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    double precision, intent(out)   :: e(:)    ! Out: Eigenvalues
    character(len=1),intent(in),optional   :: rng
    character(len=1),intent(in),optional   :: jobz
    double precision, intent(in),optional  :: vrange(2)  ! In:  energy range to be computed
    integer         , intent(in),optional  :: irange(2)  ! In:  index to be computed
    integer         , intent(out),optional :: iroots     ! Out: Number of roots found

    double precision :: vl,vu,abstol
    integer          :: info
    integer          :: nh1, nh2, il, ir,n,k, alloc
    character(len=1) :: rng_,jobz_
    integer          :: ib,jb,niwork,nwork,ldz,msize
    integer(hik)     :: hsize
    double precision,allocatable :: a(:,:),work(:),twork(:)
    integer,allocatable :: iwork(:),isuppz(:),tiwork(:)
    !
    if (verbose>=4) call TimerStart('lapack_dsyevr: diagonalization')
    !
    rng_ = 'A'
    jobz_ = 'V'
    il = 1   ; ir = 100000
    vl = -.001 ; vu = 100000.0
    !
    if (present(rng)) then
       rng_ = rng
    endif
    !
    if (present(jobz)) then
       jobz_ = jobz_
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
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    !
    msize = nh2
    if (rng_=='I') msize = ir-il+1
    !
    if (nh1/=size(e)) then
      write (out,"(' dsyevr inconsistent sizes of h and e :',2i8)") nh1/=size(e)
      stop 'lapack_dsyevr inconsistent sizes of h and e'
    end if
    !
    abstol =(small_)
    !
    nwork = 50*nh1
    niwork = 20*nh1
    ldz = 2*nh1
    !
    if (verbose>=4) write(out,'("jobz,rng = ",2a," vl,vu,il,iu =  ",2f12.2,2i8)') jobz_,rng_,vl,vu,il,ir
    !
    if (verbose>=6) then 
      write(out,"(//'MATRIX:')")
      write(out,"(i8)") nh1
      do ib=1,nh1
        do jb=ib,nh1
          write(out,"(2i8,f19.8)") ib,jb,h(jb,ib)
        enddo
      enddo
    endif 
    !
    if (verbose>=4) write(out,"(/'nwork = ',i8,' niwork = ',i8)") nwork,niwork
    !
    allocate(work(nwork),stat=alloc)
    call ArrayStart('lapack_dsyevr-arrays-work',alloc,size(work),kind(work))
    allocate(iwork(niwork),stat=alloc)
    call ArrayStart('lapack_dsyevr-arrays',alloc,size(iwork),kind(iwork))
    allocate(isuppz(2*msize),stat=alloc)
    call ArrayStart('lapack_dsyevr-arrays',alloc,size(isuppz),kind(isuppz))
    !
    if (verbose>=4) call MemoryReport
    !if (verbose>=4) write(out,"(/'range = ',a)") rng_
    !
    hsize = int(nh1,hik)*int(msize,hik)
    !
    allocate (a(nh1,msize),stat=alloc)
    call ArrayStart('lapack_dsyevr-a',alloc,1,kind(a),hsize)
    !
    if (verbose>=4) write(out,"(/'abstol = ',e18.8)") abstol
    call dsyevr(jobz_,rng_,'L',nh1,h,nh1,vl,vu,il,ir,abstol,n,e,a,&
                                                          nh1,isuppz,work,-1,iwork,-1,info)
    !
    if (verbose>=4) write(out,"(/'nwork = ',i8,' niwork = ',i8,' changed to ',2i8)") nwork,niwork,int(work(1)),int(iwork(1))
    !
    nwork = int(work(1))
    niwork = int(iwork(1))
    !
    deallocate(work,iwork)
    call ArrayStop('lapack_dsyevr-arrays')
    call ArrayStop('lapack_dsyevr-arrays-work')
    !
    allocate(twork(nwork),stat=alloc)
    call ArrayStart('lapack_dsyevr-arrays-twork',alloc,size(twork),kind(twork))
    allocate(tiwork(niwork),stat=alloc)
    call ArrayStart('lapack_dsyevr-arrays',alloc,size(tiwork),kind(tiwork))
    !
    if (verbose>=4) call MemoryReport
    !
    call dsyevr(jobz_,rng_,'L',nh1,h,nh1,vl,vu,il,ir,abstol,n,e,a,&
                nh1,isuppz,twork,nwork,tiwork,niwork,info)
    !
    !$omp parallel do private(k) shared(h) schedule(dynamic)
    do k=1,msize
       h(:,k) = a(:,k)
    enddo
    !$omp end parallel do

    !h(1:nh1,1:nh2) = a(1:nh1,1:nh2)
    !
    deallocate(a,twork,tiwork,isuppz)
    call ArrayStop('lapack_dsyevr-a')
    call ArrayStop('lapack_dsyevr-arrays')
    call ArrayStop('lapack_dsyevr-arrays-twork')
    !
    if (verbose>=6) then 
      write(out,"(//'dsyevr solution:')") 
      do ib=1,n
         write(out,"(i8,f19.8)") ib,e(ib)
         do jb=1,nh1
           write(out,"(2i8,f19.8)") ib,jb,h(jb,ib)
         enddo
      enddo
      !
      call MemoryReport
      !
    endif 
    if (present(iroots))  iroots = n
    !
    if (info/=0) then
      write (out,"(' dsyevr returned ',i8)") info
      stop 'lapack_dsyevr - dsyevr failed'
    end if
    !
    if (verbose>=4) call TimerStop('lapack_dsyevr: diagonalization')
    !
  end subroutine lapack_dsyevr



  subroutine lapack_dsyevx(h,e,iroots,vrange,irange,tol)
    double precision, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    double precision, intent(out)   :: e(:)    ! Out: Eigenvalues
    double precision, intent(in),optional  :: vrange(2)  ! In:  energy range to be computed
    integer         , intent(in),optional  :: irange(2)  ! In:  index to be computed
    integer         , intent(out),optional :: iroots     ! Out: Number of roots found
    double precision, intent(in),optional  :: tol    ! In: abs. tolerance 

    double precision :: vl,vu,abstol
    integer          :: info
    integer          :: nh1, nh2, il, ir,n,k, alloc
    character(len=1) :: rng,jobz
    integer          :: ib,jb,niwork,nwork,ldz,msize
    double precision,allocatable :: a(:,:),work(:),twork(:)
    integer,allocatable :: iwork(:),ifail(:)
    !
    if (verbose>=2) call TimerStart('lapack_dsyevx: diagonalization')
    !
    rng = 'A'
    jobz = 'V'
    il = 1   ; ir = 100000
    vl = -.001 ; vu = 100000.0
    !
    if (present(vrange)) then
       rng = 'V'
       vl = vrange(1) ; vu = vrange(2)
    endif
    !
    abstol = sqrt(small_)
    if (present(tol))  abstol = tol
    !
    if (present(irange)) then
       rng = 'I'
       il = irange(1) ; ir = irange(2)
       if (irange(2)>=size(h,dim=1)) rng = 'A' 
    endif
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    !
    msize = nh2
    if (rng=='I') msize = ir-il+1
    !
    if (nh1/=size(e)) then
      write (out,"(' dsyevr inconsistent sizes of h and e :',2i8)") nh1/=size(e)
      stop 'lapack_dsyevx inconsistent sizes of h and e'
    end if
    !
    nwork = 50*nh1
    niwork = 5*nh1
    !
    if (verbose>=4) then 
      write(out,"(//'MATRIX:')") 
      do ib=1,nh1
        do jb=ib,nh1
          write(out,"(2i8,f19.8)") ib,jb,h(jb,ib)
        enddo
      enddo
    endif 
    !
    if (verbose>=4) write(out,"(/'nwork = ',i8,' niwork = ',i8)") nwork,niwork
    !
    allocate(work(nwork),stat=alloc)
    call ArrayStart('lapack_dsyevx-arrays-work',alloc,size(work),kind(work))
    allocate(iwork(niwork),stat=alloc)
    call ArrayStart('lapack_dsyevx-arrays',alloc,size(iwork),kind(iwork))
    allocate(ifail(nh1),stat=alloc)
    call ArrayStart('lapack_dsyevx-arrays',alloc,size(ifail),kind(ifail))
    !
    if (verbose>=3) call MemoryReport
    !if (verbose>=3) write(out,"(/'range = ',a)") rng
    !
    if (verbose>=4) write(out,"(/'abstol = ',e18.8)") abstol
    call dsyevx('V',rng,'L',nh1,h(1:nh1,1:nh2),nh1,vl,vu,il,ir,abstol,n,e(1:nh1),h(1:nh1,1:nh2),&
                                                          nh1,work,-1,iwork,ifail,info) 
    !
    nwork = int(work(1))
    !
    if (verbose>=4) write(out,"(/'nwork = ',i8)") nwork
    !
    deallocate(work,iwork)
    call ArrayStop('lapack_dsyevx-arrays')
    call ArrayStop('lapack_dsyevx-arrays-work')
    !
    allocate(twork(nwork),stat=alloc)
    call ArrayStart('lapack_dsyevx-arrays-twork',alloc,size(twork),kind(twork))
    !
    allocate (a(nh1,msize),stat=alloc)
    call ArrayStart('lapack_dsyevx-a',alloc,size(a),kind(a))
    !
    if (verbose>=3) call MemoryReport
    !
    call dsyevx('V',rng,'L',nh1,h(1:nh1,1:nh2),nh1,vl,vu,il,ir,abstol,n,e(1:nh1),a(1:nh1,1:msize),&
                nh1,twork(1:nwork),nwork,iwork(1:nh1),ifail,info)
    !
    !$omp parallel do private(k) shared(h) schedule(dynamic)
    do k=1,msize
       h(:,k) = a(:,k)
    enddo
    !$omp end parallel do

    !h(1:nh1,1:nh2) = a(1:nh1,1:nh2)
    !
    deallocate(a,twork,ifail)
    call ArrayStop('lapack_dsyevx-a')
    call ArrayStop('lapack_dsyevx-arrays')
    call ArrayStop('lapack_dsyevx-arrays-twork')
    !
    if (verbose>=5) then 
      write(out,"(//'dsyevr solution:')") 
      do ib=1,n
         write(out,"(i8,f19.8)") ib,e(ib)
         do jb=1,nh1
           write(out,"(2i8,f19.8)") ib,jb,h(jb,ib)
         enddo
      enddo
      !
      call MemoryReport
      !
    endif 
    if (present(iroots))  iroots = n
    !
    if (info/=0) then
      write (out,"(' dsyevr returned ',i8)") info
      stop 'lapack_dsyevx - dsyevr failed'
    end if
    !
    if (verbose>=2) call TimerStop('lapack_dsyevx: diagonalization')
    !
  end subroutine lapack_dsyevx

 

  subroutine lapack_ginverse_real(amat)
    real, intent(inout) :: amat(:,:) ! In: matrix to invert
                                         ! Out: generalized inverse of the matrix

    real :: eval(size(amat,dim=1))
    real :: evec(size(amat,dim=1),size(amat,dim=1))
    real :: eps

    evec = amat
    call lapack_syev(evec(:,:),eval(:))
    eps = 100.0*spacing(maxval(eval))
    where (abs(eval)>eps) 
      eval = 1.0 / sqrt(eval)
    elsewhere  
      eval = 0.0
    end where
    evec = evec * spread(eval,dim=1,ncopies=size(evec,dim=1))
    amat = matmul(evec,transpose(evec))

  end subroutine lapack_ginverse_real

  subroutine lapack_ginverse_double(amat)
    double precision, intent(inout) :: amat(:,:) ! In: matrix to invert
                                         ! Out: generalized inverse of the matrix

    double precision :: eval(size(amat,dim=1))
    double precision :: evec(size(amat,dim=1),size(amat,dim=1))
    double precision :: eps

    evec = amat
    call lapack_syev(evec(:,:),eval(:))
    eps = 100.0d0*spacing(maxval(eval))
    !
    if (any(eval(:)<0)) then 
      !
      write(out,"('lapack_ginverse_double: negative evals')")
      stop 'lapack_ginverse_double: negative evals'
      !
    endif
    !
    where (abs(eval)>eps) 
      eval = 1.0d0 / sqrt(eval)
    elsewhere  
      eval = 0.0d0
    end where
    evec = evec * spread(eval,dim=1,ncopies=size(evec,dim=1))
    amat = matmul(evec,transpose(evec))

  end subroutine lapack_ginverse_double



  subroutine dseupd_arpack(n,bterm,nroots,factor,maxitr_,tol,h,e,IO)

    integer,intent(in)     :: n,bterm(n,2),maxitr_
    integer,intent(inout)  :: nroots
    double precision,intent(in)    :: factor
    double precision,intent(in)    :: tol
    double precision,intent(inout) :: h(:,:)
    double precision,intent(out)   :: e(nroots)
    character(len=cl),optional     :: IO

    double precision,allocatable :: v(:,:),workl(:),workd(:),d(:,:),resid(:),vect_chk(:)
    !
    logical,allocatable   ::   select(:)
    integer(ik) ::   iparam(11), ipntr(11),ldv,iter,iparam_chk(11),ipntr_chk(11)
    integer(hik) ::  msize, np, nconv, nev2
    !

    character(len=1) ::  bmat
    character(len=2) ::  which
    character(len=cl)::  ioname,IO_job
    character(len=10)::  buf
    !
    integer(ik)   :: ido,  nev, ncv, lworkl, info, ierr, j, &
                     mode, ishfts, alloc, maxitr, n_chk,nev_chk,iounit
                     !
    logical          rvec
    double precision  ::  sigma,rnorm
    !
    character(len=20) ::  valfmt
    !
    double precision  :: tol_chk
    integer(ik)       :: lworkl_chk,info_chk,ncv_chk,i,ifmt,len,nperli

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


      if (verbose>=2) call TimerStart('dseupd_arpack: diagonalization')
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

      if (verbose>=4) then
        msgets = 1
        msaitr = 0
        msapps = 0
        msaupd = 1
        msaup2 = 0
        mseigt = 1
        mseupd = 0
      endif 
      !
      if (verbose>=5) then
        msgets = 1
        msaitr = 2
        msapps = 0
        msaupd = 1
        msaup2 = 3
        mseigt = 0
        mseupd = 1
      endif 
      !
      if (present(IO)) IO_job = IO
      !
      if (verbose>=3) msaupd = 1
      !
      nev = nroots
      !
      if (nroots==0) nev=max(100,n)
      !
      !ncv = iwork
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
      msize = ldv*ncv
      call ArrayStart('dseupd_arpack:v',alloc,1,kind(v),msize)
      call ArrayStart('dseupd_arpack:workl',alloc,size(workl),kind(workl))
      call ArrayStart('dseupd_arpack:workd',alloc,size(workd),kind(workd))
      call ArrayStart('dseupd_arpack:d',alloc,size(d),kind(d))
      call ArrayStart('dseupd_arpack:resid',alloc,size(resid),kind(resid))
      call ArrayStart('dseupd_arpack:select',alloc,size(select),kind(select))
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
      if (verbose>=3) write(out,"(/'Arpack: n = ',i8,' nev = ',i8,' ncv = ',i8,' maxitr = ',i8,' tol = ',e12.5)") n,nev,& 
          ncv,maxitr,tol
      !
      iter = 0 
      !
      if (verbose>=3) call TimerStart('dseupd_arpack: Eigenvalues')
      !
      if (trim(IO_job)=='READ') then
         !
         ioname ='Eigenvectors in symm. representation'
         call IOStart(trim(ioname),iounit)
         !
         read(iounit,"(a10)") buf
         !
         if (buf /= 'next gamma') then 
           write(out,"('Arpack-read_symm_eigenvect error: wrong header ',a)") buf
           stop 'read_symm_eigenvect error: wrong header'
         endif
         !
         read(iounit,"(i2,a1,i14,a2,i14,d23.16,16x,/,12i5,12x,/,13i5,7x,/,i5,d23.16,i5,i5)") & 
                        ido, bmat, n_chk, which, nev_chk, tol_chk, ncv_chk,iparam_chk, &
                        ipntr_chk, lworkl_chk, info_chk, np, rnorm, nconv, nev2

        !
        if (verbose>=3) write(out,"(/'Arpack-stored: n = ',i8,' nev = ',i8)") n_chk,nev_chk
        !
        if (n_chk>n.or.nev_chk<ncv) then 
          write(out,"('dseupd_arpack error: stored vectors are inconsistent: either  n_chk>n = ',2i8,', or nev_chk<ncv = ',2i8)") &
                n_chk,n,nev_chk,ncv
          stop 'dseupd_arpack error: stored vectors are inconsistent'
        endif
        !
        valfmt = '(3e22.16)'
        resid = 0
        read(iounit,valfmt) (resid(i), i = 1, n_chk)
        !
        workd = 0
        read(iounit,valfmt) (workd(i), i =  1, n_chk)
        read(iounit,valfmt) (workd(i), i =  n_chk+1, 2*n_chk)
        read(iounit,valfmt) (workd(i), i =2*n_chk+1, 3*n_chk)
        !
        workl = 0 
        read(iounit,valfmt) (workl(i), i = 1, min(lworkl_chk,lworkl))
        !
        v = 0 
        do j=1,min(ncv_chk,ncv)
           read(iounit,valfmt) (v(i,j), i = 1, n_chk)
        enddo
        !
        !resid = 0
        !
        workd = 0
        !
        !workl = 0 
        !
        ido = -2
        !
        np = ncv - nev
        nev2 = nev
        rnorm = 1e-5
        !
#if (arpack_ > 0)
        !
        call dsaupd( ido, bmat, n, which, nev, tol, resid, &
                      ncv, v, ldv, iparam, ipntr, workd, workl, &
                      lworkl, info) !,np, rnorm, nconv, nev2 )
        !
#else 
            !
            write(out,"(/'Arpack was not activated yet. Please uncomment dsaupd and  dseupd bellow')")
            stop 'Arpack was not activated'
            !
#endif
        !
        ido = -1
        !
      endif 
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

#if (arpack_ > 0)
           !
           call dsaupd ( ido, bmat, n, which, nev, tol, resid, &
                         ncv, v, ldv, iparam, ipntr, workd, workl, &
                         lworkl, info) !,np, rnorm, nconv, nev2)
#else 
            !
            write(out,"(/'Arpack was not activated yet. Please uncomment dsaupd and  dseupd bellow')")
            stop 'Arpack was not activated'
            !
#endif
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
          call matvec(n,bterm,h,workd(ipntr(1)), workd(ipntr(2)))
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

      if (trim(IO)=='SAVE') then
         !
         ioname ='Eigenvectors in symm. representation'
         call IOStart(trim(ioname),iounit)
         !
         write(iounit,"(a10)") 'next gamma'
         !
         !    After maxitr iterations without convergence,
         !    output the computed quantities to the file state.
         !
         write(iounit,"(i2,a1,i14,a2,i14,d23.16,16x,/,12i5,12x,/,13i5,7x,/,i5,d23.16,i5,i5)") &
                        ido, bmat, n, which, nev, tol, ncv, iparam, &
                        ipntr, lworkl, info, np, rnorm, nconv, nev2 
         ifmt = 16
         len  = ifmt + 6
         nperli = 3
         write(valfmt,"(1h(,i1,1he,i2,1h.,i2,1h))") nperli,len,ifmt
         write(iounit,valfmt) (resid(i), i = 1, n)
         write(iounit,valfmt) (workd(i), i = 1, 3*n)
         write(iounit,valfmt) (workl(i), i = 1, lworkl)
         do j=1,ncv
            write(iounit,valfmt) (v(i,j), i = 1, n)
         enddo
         !
      endif
      !
      if (verbose>=3) call TimerStop('dseupd_arpack: Eigenvalues')
      !
      !  No fatal errors occurred. Post-Process using DSEUPD. 
      !  Computed eigenvalues may be extracted. 
      !
      !  Eigenvectors may also be computed now if desired.  (indicated by rvec = .true.)  
      !             
      rvec = .true.
      !
#if (arpack_ > 0)
        !
        if (verbose>=5) write(out,"(/'Arpack: dseupd')") 
        !
        if (verbose>=2) call TimerStart('dseupd_arpack: Eigenvectors')
        !
        call dseupd ( rvec, 'All', select, d, v, ldv, sigma, &
                      bmat, n, which, nev, tol, resid, ncv, v, ldv, &
                      iparam, ipntr, workd, workl, lworkl, ierr )
        !
        if (verbose>=2) call TimerStop('dseupd_arpack: Eigenvectors')
        !
        if (verbose>=5) write(out,"(/'Arpack: done!')") 
        !
#else 
          !
          write(out,"(/'Arpack was not activated yet. Please uncomment dsaupd and  dseupd bellow')")
          stop 'Arpack was not activated'
          !
#endif
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
      !$omp parallel do private(j) shared(h) schedule(static)
      do j=1,n
         !
         call dcopy(nroots,v(j,1:nroots),1,h(j,1:nroots),1)
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
      if (verbose>=3) write(out,"('Total number of iteratons  = ',i8/)") iter 
      !
      call ArrayStop('dseupd_arpack:v')
      call ArrayStop('dseupd_arpack:workl')
      call ArrayStop('dseupd_arpack:workd')
      call ArrayStop('dseupd_arpack:d')
      call ArrayStop('dseupd_arpack:resid')
      call ArrayStop('dseupd_arpack:select')
      !
      if (verbose>=2) call TimerStop('dseupd_arpack: diagonalization')
      !
  end subroutine dseupd_arpack
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
#if (blacs_ > 0)
        call BLACS_GRIDINFO( comm, nprow, npcol, myprow, mypcol )
        myid= myprow
#elif (mpi_ > 0)
        call MPI_COMM_RANK( comm, myid, ierr )
        call MPI_COMM_SIZE( comm, nprocs_, ierr )
        if (nprocs_/=nprocs) then
          write(out,"('matvec_p: inconsistent number of  nprocs  = ',2i8)") nprocs_,nprocs
          stop 'matvec_p: inconsistent number of  nprocs s'
        endif
#endif
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
#if (blacs_ > 0)
             call dgerv2d( comm, nx, 1, mv_buf(istart:iend), nx, iprev, mypcol )
#endif
#if (mpi_ > 0)
             call mpi_recv(mv_buf(istart),nx,MPI_DOUBLE_PRECISION,iprev,myid,comm,ierr)
#endif
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
#if (blacs_ > 0)
             call dgesd2d( comm, nx, 1, z(istart:iend), nx, iprev, mypcol )
#endif
#if (mpi_ > 0)
             call mpi_send(mv_buf(istart),nx,MPI_DOUBLE_PRECISION,iprev,iprev,comm,ierr)
#endif
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
#if (blacs_ > 0)
             call dgerv2d( comm, nx, 1, mv_buf(istart:iend), nx, inext, mypcol )
#endif
#if (mpi_ > 0)
             call mpi_recv(mv_buf(istart),nx,MPI_DOUBLE_PRECISION,inext,myid,comm,ierr)
#endif
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
#if (blacs_ > 0)
             call dgesd2d( comm, nx, 1, z(istart:iend), nx, inext, mypcol )
#endif
#if (mpi_ > 0)
             call mpi_send(mv_buf(istart),nx,MPI_DOUBLE_PRECISION,inext,inext,comm,ierr)
#endif
           !
         enddo
         !
      enddo 
      !
      !enddo
      !
   return
  end subroutine matvec_p



  subroutine dseupd_p_arpack(n,bterm,nroots,factor,maxitr_,tol,h,e)

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
    integer,parameter ::  maxnprocs=256
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
      integer,parameter :: MPI_COMM_WORLD=0
!
 
!     %---------------%
!     | MPI INTERFACE |
!     %---------------%
 
      integer           myid, rc


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

#if (blacs_ > 0)
        call BLACS_PINFO( iam, nprocs )
        blacs_or_mpi = 'BLACS'
#endif

#if (mpi_ > 0)
        call MPI_INIT( ierr )
        comm = MPI_COMM_WORLD
        call MPI_COMM_RANK( comm, myid, ierr )
        call MPI_COMM_SIZE( comm, nprocs, ierr )
        !
        print *,comm,myid,nprocs
        !
        if (trim(blacs_or_mpi)=='BLACS') then
          write(out,"('dseupd_p_arpack: IT IS ILLEGAL TO USE MPI AND BLACS AT THE SAME TIME')")
          stop 'dseupd_p_arpack: IT IS ILLEGAL TO USE MPI AND BLACS AT THE SAME TIME'
        endif
        !
        blacs_or_mpi = 'MPI'
        !
#endif
      !


!
!     If in PVM, create virtual machine if it doesn't exist
!
      if (nprocs .lt. 1) then
         nprocs = 1
#if (blacs_ > 0)
           call BLACS_SETUP( iam, nprocs )
#endif
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
      if (verbose>=2) call TimerStart('dseupd_p_arpack: diagonalization')
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
      call ArrayStart('dseupd_p_arpack',alloc,size(v),kind(v))
      call ArrayStart('dseupd_p_arpack',alloc,size(workl),kind(workl))
      call ArrayStart('dseupd_p_arpack',alloc,size(workd),kind(workd))
      call ArrayStart('dseupd_p_arpack',alloc,size(d),kind(d))
      call ArrayStart('dseupd_p_arpack',alloc,size(resid),kind(resid))
      call ArrayStart('dseupd_p_arpack',alloc,size(select),kind(select))
      call ArrayStart('dseupd_p_arpack',alloc,size(mv_buf),kind(mv_buf))

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
#if (blacs_ > 0)
        call BLACS_GET( 0, 0, comm )
        call BLACS_GRIDINIT( comm, 'Row', nprow, npcol )
        call BLACS_GRIDINFO( comm, nprow, npcol, myprow, mypcol )
#endif
      !
      if (verbose>=2.and.trim(blacs_or_mpi)=='BLACS') write(out,"('myprow, nprow, mypcol, npcol, nprocs = ',5i8)") & 
          myprow, nprow, mypcol, npcol, nprocs
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
      !call ArrayStart('dseupd_p_arpack',alloc,size(mv_buf),kind(mv_buf))
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

#if (blacs_ > 0 || mpi_ > 0)
           !
           call pdsaupd ( comm, ido, bmat, nloc, which, nev, tol, resid, &
                          ncv, v, ldv, iparam, ipntr, workd, workl, &
                          lworkl, info )
        !
#elseif (arpack_>0)
            !
           call dsaupd ( ido, bmat, n, which, nev, tol, resid, &
                          ncv, v, ldv, iparam, ipntr, workd, workl, &
                          lworkl, info )
            !
#else
            write(out,"(/'Arpack was not activated yet. Please uncomment dsaupd and  dseupd bellow')")
            stop 'Arpack was not activated'

#endif
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
#if (blacs_ > 0)
        !
        if (verbose>=5) write(out,"(/'Arpack: dseupd')") 
        !
        call pdseupd ( comm, rvec, 'All', select, d, v, ldv, sigma, &
                       bmat, n, which, nev, tol, resid, ncv, v, ldv, &
                       iparam, ipntr, workd, workl, lworkl, ierr )


        if (verbose>=5) write(out,"(/'Arpack: done!')") 

        !                    
#else 
          !
          write(out,"(/'Arpack was not activated yet. Please uncomment dsaupd and  dseupd bellow')")
          stop 'Arpack was not activated'
          !
#endif
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
#if (blacs_ > 0)
        call BLACS_GRIDEXIT ( comm )
        call BLACS_EXIT(0)
#endif
#if (mpi_ > 0)
         call MPI_FINALIZE(rc)
#endif


      deallocate(v,workl,workd,d,resid,select,mv_buf)
      !
      call ArrayStop('dseupd_p_arpack')
      !
      if (verbose>=2) call TimerStop('dseupd_p_arpack: diagonalization')
      !
  end subroutine dseupd_p_arpack


!subroutine BLACS_GRIDINFO( comm, nprow, npcol, myprow, mypcol )
!
!  integer :: comm, nprow, npcol, myprow, mypcol
!
!  !
!  myprow = 1
!  mypcol = 1
!  !
!end subroutine BLACS_GRIDINFO




end module lapack
