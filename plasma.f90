module plasma
 
#define plasma_ 0

!
!  Simplistic type-agnostic PLASMA interface
!
  use accuracy
  use timer
  implicit none
  private verbose

  interface plasma_sytrdx
    module procedure plasma_diag_dsytrdx
  end interface ! plasma_sytrdx

  integer,parameter:: verbose = 5
  !
  contains


  subroutine plasma_diag_dsytrdx(n,a,w,nroots,Ethres_)
    !
#if (plasma_ > 0) 
      INCLUDE "plasmaf.h"
      integer(ik), parameter :: VEC = PlasmaVec
      integer(ik), parameter :: UPLO = PlasmaLower
      EXTERNAL PLASMA_DSYTRDX
      INTEGER PLASMA_DSYTRDX
#endif
    !
    integer         , intent(in)    :: n
    double precision, intent(inout) :: a(n,n)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    double precision, intent(out)   :: w(n)    ! Out: Eigenvalues
    integer         , intent(out)    :: nroots
    double precision,optional,intent(in)   :: Ethres_
    integer          :: info,nprocs_,VBLKSIZ
    integer          :: nh,alloc,k,OMP_GET_NUM_THREADS,tid
    integer          :: msize,nb,ib,lda,ldq,corea,coreb,corec,i,j,nvect
    integer(hik)     :: hsize,jhik
    !double precision,allocatable :: a(:,:),w(:),d(:)
    !
    double precision :: mat_elem,dot_product, Ethres
    real(rk) :: real_end,real_start,cpu_end,cpu_start,real_time,cpu_time
    !
    double precision q(1),d(1),e(1)
    pointer (pq,q), (pd,d),  (pe,e)
    !
    integer*8 lena, len
    !
    if (verbose>=4) call TimerStart('plasma_sytrdx: diagonalization')
    !
    w = 0 
    !
    Ethres = 1000000.0d0
    !
    nroots = 1
    !
    nvect = n
    !
    ! estimate number of roots using the diagonal elements
    nvect = 1
    do while (a(nvect,nvect)<Ethres_.and.nvect<n)
      nvect = nvect+1 
    enddo
    !
    ! Use the energy threshold if given in as input
    !
    if (present(Ethres_)) Ethres = Ethres_
    !
#if (plasma_ == 0) 
       write(out,"('Plasma is not activated, in plasma.f90 please set plasma_ to 1')")
#endif
    !
#if (plasma_ > 0)
      !
      !call getsize(N,LDA,nprocs_)
      !
#endif
    !  
    !$omp parallel private(tid)
      if (tid==0) then
         nprocs_ = OMP_GET_NUM_THREADS()
      endif
    !$omp end parallel
    !
    corea = nprocs_
    !
    if (n>=64000)  corea = max( min(int( ( 632.0+6.0*N/1000.0 )/17.0,ik ),nprocs_),1 )
    if (n>200000)  corea = min( 132,nprocs_ )
    !
    COREB = nprocs_
    COREC = (nprocs_/8) * 2
    VBLKSIZ = 64
    !
    ! block and tile sizes
    !
    LDA = N
    !
    LDQ = LDA
    NB = 720
    !
    !if (n>150000) NB = 384
    !if (n>200000) NB = 448
    !
    IB = 144
    !
    if (verbose>=4) write(out,"('nthreads = ',i7)") nprocs_
    !
    ! Check the parallel environment
    if (verbose>=4) then
      !
      write(out,"('COREA (nbrthrdA) = ',i7)") COREA
      call print_openmp
      !
    endif
    !
    ! Set OMP_SET_NUM_THREADS to 1: important for plasma!
    !
    call OMP_SET_NUM_THREADS(1)
    !
    !omp parallel private(tid)
    !  if (tid==0) then
    !     call OMP_SET_NUM_THREADS(1)
    !  endif
    !omp end parallel
    !
    if (verbose>=4) then
      !
      write(out,"('Check if OMP_SET_NUM_THREADS is 1:')") 
      !
      call print_openmp
      !
    endif
    !
    !call mkl_set_num_threads(1)
    !
    call estimate_plasma_memory(N,NB,IB,VBLKSIZ,nvect,hsize)
    !
    hsize = hsize - int(n*n+n,hik)
    !
    alloc = 0
    !
    call ArrayStart('plasma_sytrdx-all-plasma',alloc,1,kind(w),hsize)
    !
    len = 8
    !
    lena = len * N * LDQ
    pq = malloc(lena)
    if (pq .eq. 0) then
      print*,'cannot allocate matrix Q',N,N
      stop
    end if
    !
    lena = len * N
    pd = malloc(lena)
    if (pd .eq. 0) then
      print*,'cannot allocate matrix d',N
      stop
    end if
    !
    pe = malloc(lena)
    if (pe .eq. 0) then
      print*,'cannot allocate matrix e',N
      stop
    end if
    !
    alloc = 0
    !
    !if (verbose>=4) write(out,"(' allocation of a matrix of ',i,'x',i)") N, N
    !hsize = int(n,hik)*int(n,hik)
    !
    if (verbose>=4) call MemoryReport
    !
    if (verbose>=3) write(out,"(/   'Diagonalization with PLASMA_dsytrdX...')")
    !
    real_start = get_real_time()
    cpu_start  = get_cpu_time ()
    !
    ! call lapack_dsyev(mat,e,jobz)
    !
    info = -1
    ! 
    ! set up my PLASMA_DSYTRDX environmenta, then call the eigensolver
    !
#if (plasma_ > 0)
        !
        call resetcore(corea,corec)
        CALL PRINTARGS(VEC, UPLO, N, LDA, LDQ, COREA, COREB, COREC, NB, IB)
        CALL SETPLASMAENV(COREA)
        !
        print*,'check args again'
        !
        CALL PRINTARGS(VEC, UPLO, N, LDA, LDQ, COREA, COREB, COREC, NB, IB)
        !
        INFO = PLASMA_DSYTRDX( VEC, UPLO, N, A, LDA, D, E, W, Q, LDQ, &
                               COREA, COREB, COREC, NB, IB,Ethres )
        IF (INFO .NE. PLASMA_SUCCESS) THEN
          PRINT*,"PLASMA_DSYTRDX RETURNS ", INFO, " ABORTING!"
          STOP
        END IF
        CALL USETPLASMAENV()
#endif
    !
    real_end = get_real_time()
    cpu_end  = get_cpu_time ()
    !
    real_time = real_end - real_start
    cpu_time = cpu_end - cpu_start
    !
    write (out,"(t2,'real ',f9.1,'s ; cpu ',f9.1,'s')")  real_time, cpu_time
    !
    if (verbose>=3) write(out,"('   ...done!')") 
    !
    nroots = 0
    !$omp parallel do private(k,jhik) shared(A) reduction(+:nroots)  schedule(dynamic)
    do k=1,n
       !
       if (w(k)>Ethres) cycle 
       !
       jhik = int(n*(k-1)+1,hik)
       !
       A(:,k) = q(jhik:)
       !
       nroots  = nroots + 1
       !
    enddo
    !$omp end parallel do
    !
    !if (verbose>=5) then 
    !  write(out,"(//'Solution:')") 
    !  do i=1,n
    !     write(out,"(i8,f19.8)") i,e(i)
    !     do j=1,n
    !       write(out,"(2i8,f19.8)") i,j,a(j,i)
    !     enddo
    !  enddo
    !  !
    !  call MemoryReport
    !  !
    !endif 
    !
    !
    if (verbose>=5) then 
      write(out,"(/'Eigenvalues:')")
      !
      !omp parallel do private(i,j,mat_elem) schedule(dynamic)
      do i=1,n
         !
         write(out,"(i8,f19.8)") i,w(i)-w(1)
         !
         !do j=i,n
         !  !
         !  mat_elem = dot_product(h(:,i),h(:,j))
         !  !
         !  write(out,"(2i8,g19.8)") i,j,mat_elem
         !  !
         !enddo
      enddo
      !omp end parallel do
      !
      call MemoryReport
      !
    endif 
    !
    !deallocate(a,d,w)
    !
    call free(pq)
    call free(pd)
    call free(pe)
    !
    call ArrayStop('plasma_sytrdx-all-plasma')
    !
    if (verbose>=4) then
      call print_openmp
    endif
    !
    !$omp parallel private(tid)
      if (tid==0) then
         call OMP_SET_NUM_THREADS(nprocs_)
      endif
    !$omp end parallel
    !
    if (verbose>=4) call TimerStop('plasma_sytrdx: diagonalization')
    !
  end subroutine plasma_diag_dsytrdx

  subroutine print_openmp
    !
    implicit none
    integer(ik)      :: TID,PROCS,NTHREADS,MAXT
    logical          :: INPAR,DYNAMIC,NESTED,OMP_IN_PARALLEL,OMP_GET_DYNAMIC,OMP_GET_NESTED
    integer(ik)      :: OMP_GET_THREAD_NUM,OMP_GET_NUM_PROCS,OMP_GET_MAX_THREADS,OMP_GET_NUM_THREADS
    !
    !     Start parallel region
    !$OMP PARALLEL PRIVATE(NTHREADS, TID)
    !
    !  Obtain thread number
    TID = OMP_GET_THREAD_NUM()
    !
    !     Only master thread does this
    IF (TID .EQ. 0) THEN
        !
        PRINT *, 'Thread',tid,'getting environment information'
        !     Get environment information
        PROCS = OMP_GET_NUM_PROCS() 
        NTHREADS = OMP_GET_NUM_THREADS()
        MAXT = OMP_GET_MAX_THREADS()
        INPAR = OMP_IN_PARALLEL()
        DYNAMIC = OMP_GET_DYNAMIC()
        NESTED = OMP_GET_NESTED()
        !     Print environment information
        !
        PRINT *, 'Number of processors = ', PROCS
        PRINT *, 'Number of threads = ', NTHREADS
        PRINT *, 'Max threads = ', MAXT
        PRINT *, 'In parallel? = ', INPAR
        PRINT *, 'Dynamic threads enabled? = ', DYNAMIC
        PRINT *, 'Nested parallelism supported? = ', NESTED
        !
        !call OMP_SET_NUM_THREADS(1) 
        !
        !
      END IF
      !     Done
    !$OMP END PARALLEL
    !
  end subroutine print_openmp


  subroutine estimate_plasma_memory(N,NB,IB,VBLKSIZ,nevect,hsize)
   !
   integer,intent(in)  :: n,NB,ib,vblksiz,nevect
   integer(8),intent(out) :: hsize
   real(rk) :: m,bcnt,n_,NB_,ib_,vblksiz_,nevect_
   real(rk) :: total
     !
     N_ = real(N,rk)
     NB_= real(NB,rk)
     IB_= real(IB,rk)
     VBLKSIZ_ = real(VBLKSIZ,rk)
     nevect_ = real(nevect,rk)
     !
     if (verbose>=4) print*,'N,NB,IB,VBLKSIZ = ',N,NB,IB,VBLKSIZ
     total = 0
     total = total + n_*n_
     total = total + (n_+ib_/nb_+1.0_rk)*(n_+ib_/nb_+1.0_rk)
     total = total + (2.0_rk*nb_+6.0_rk)*n_
     if (verbose>=4) print*,'Memory required for computing plasma eigenvalues'
     if (verbose>=4) print*,'only is ',total+6.0_rk*n_,' double words.'
     m = n_/vblksiz_ + 1.0_rk
     bcnt = m* (n_/nb_ + 1.0_rk) - ( m*(m-1.0_rk)*0.5_rk )*vblksiz_/nb_
     total = total + vblksiz_*bcnt*vblksiz_
     total = total + bcnt*vblksiz_
     total = total + (nb_+vblksiz_-1.0_rk)*bcnt*vblksiz_
     total = total + nevect_*n_
     !
     hsize = int(total,8)
     !
     if (verbose>=4) print*,'Memory required for all eigenvalues and ',nevect,' eigenvectors'
     if (verbose>=4) print*,'only is ',total,' double words',hsize/1024.0**3,'Gb'
     !
  end subroutine estimate_plasma_memory


   !
   !  Get CPU time, whatever this means
   !
   function get_cpu_time() result(t)
     real(rk) :: t
     !
     call cpu_time(t)
   end function get_cpu_time


  function get_real_time() result(t)
    real(rk) :: t
    !
    integer         :: count, count_rate, count_max
    real(rk), save :: overflow   =  0
    integer, save   :: last_count = -1
    !
    call system_clock(count,count_rate,count_max)
    !
    ! Try to detect a rollover
    !
    if (count<last_count) then
      overflow = overflow + count_max
    end if
    last_count = count
    !
    ! Convert to seconds
    !
    t = (overflow+count)/count_rate
    !
  end function get_real_time


end module plasma
