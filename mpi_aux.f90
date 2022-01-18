module mpi_aux
#ifdef TROVE_USE_MPI_
  use mpi_f08
#endif
  use timer
  use accuracy
  implicit none

  public co_init_comms, co_finalize_comms, co_init_distr, co_distr_data, co_write_matrix_distr, co_read_matrix_distr
  public co_block_type_init, co_read_matrix_distr_ordered

  public send_or_recv, comm_size, mpi_rank
  public co_startdim, co_enddim

  public blacs_size, blacs_rank
  public nprow,npcol,myprow,mypcol
  public mpi_real_size, mpi_int_size

#ifndef TROVE_USE_MPI_
  public MPI_Datatype, MPI_File, MPI_Status, MPI_Request
#endif

  interface co_sum
    module procedure :: co_sum_double
  end interface

#ifndef TROVE_USE_MPI_
  !
  ! When not using MPI, create a dummy type to maintain the interface of this module,
  ! specifically the last argument of co_block_type_init.
  ! The value passed will not be accessed.
  !
  type MPI_Datatype
    integer :: dummy = 0
  end type MPI_Datatype

  type MPI_File
    integer :: dummy = 0
  end type MPI_File

  type MPI_Status
    integer :: dummy = 0
  end type MPI_Status

  type MPI_Request
    integer :: dummy = 0
  end type MPI_Request

  integer, parameter :: MPI_OFFSET_KIND=8
#endif

  integer,dimension(:),allocatable              :: send_or_recv
  integer                                       :: comm_size, mpi_rank
  integer                                       :: co_startdim, co_enddim, co_curr_dimen, &
    co_blocksize, co_localsize
  logical                                       :: comms_inited = .false., distr_inited=.false.
  type(MPI_Datatype)                            :: mpitype_column
  type(MPI_Datatype),dimension(:), allocatable  :: mpi_blocktype
  integer                                       :: mpi_real_size, mpi_int_size

  !blacs/pblas
  integer :: blacs_size, blacs_rank, blacs_ctxt
  integer :: nprow,npcol,myprow,mypcol
  integer,dimension(2)  :: blacs_dims

contains

  subroutine co_init_blacs()
    implicit none

    if (.not. comms_inited) stop "CO_INIT_BLACS COMMS NOT INITED"

#ifdef TROVE_USE_MPI_
    ! Must be initialised to zero - if stack contains garbage here MPI_Dims_create WILL fail
    blacs_dims = 0

    call blacs_pinfo(blacs_rank, blacs_size)
    if (blacs_rank .lt. 0) return

    call MPI_Dims_create(blacs_size, 2, blacs_dims)

    call blacs_get(-1, 0, blacs_ctxt)
    call blacs_gridinit(blacs_ctxt, 'R', blacs_dims(1), blacs_dims(2))
    call blacs_gridinfo(blacs_ctxt, nprow, npcol, myprow, mypcol)

    !write(*,"('BLACS: [',i2,',',i2'](',i4,i4,i4,i4',)')") mpi_rank,blacs_rank,nprow,npcol,myprow,mypcol
#else
    blacs_size = 1
    blacs_rank = 0
    nprow = 1
    npcol = 1
    myprow = 0
    mypcol = 0
#endif
  end subroutine co_init_blacs

  subroutine co_block_type_init(smat, dimx, dimy, descr, allocinfo, mpi_type)
    implicit none

    real(rk),intent(out),dimension(:,:),allocatable   :: smat

    ! Note: allocated matrix will have total dimensions (dimy x dimx)
    integer,intent(in)                                :: dimx, dimy
    integer,intent(out),dimension(9)                  :: descr
    integer,intent(out)                               :: allocinfo

    type(MPI_Datatype),intent(out),optional           :: mpi_type

#ifdef TROVE_USE_MPI_
    integer,dimension(2)                              :: global_size, distr, dargs
    integer                                           :: MB,NB,MLOC,NLOC,ierr

    integer,external  :: NUMROC

    if (.not. comms_inited) stop "CO_BLOCK_TYPE_INIT COMMS NOT INITED"

    MB = dimy/nprow
    NB = dimx/npcol
    MLOC = NUMROC( dimy, MB, myprow, 0, nprow )
    NLOC = NUMROC( dimx, NB, mypcol, 0, npcol )
    call DESCINIT(descr, dimy, dimx, MB, NB, 0, 0, blacs_ctxt, max(MLOC,1), ierr)

    allocate(smat(MLOC,NLOC), stat=allocinfo)
    if(allocinfo.ne.0) return

    if (present(mpi_type)) then
      global_size = (/dimy, dimx/)
      distr = (/MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC/)
      dargs = (/MB, NB/)
      call MPI_Type_create_darray(blacs_size, blacs_rank, 2, global_size, distr, dargs, blacs_dims, &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mpi_type, ierr)
      call MPI_Type_commit(mpi_type, ierr)
    endif
#else
    allocate(smat(dimy,dimx), stat=allocinfo)
    if(allocinfo.ne.0) return
#endif

  end subroutine co_block_type_init

  subroutine co_sum_double(x, root_process)

    real(rk), intent(inout), dimension(..) :: x
    integer, optional :: root_process

    if (comm_size.eq.1) return
#ifdef TROVE_USE_MPI_
    call TimerStart('co_sum_double')

    if (present(root_process)) then

      if (mpi_rank .eq. 0) then
        call mpi_reduce(mpi_in_place, x, size(x), mpi_double_precision, mpi_sum, 0, mpi_comm_world)
      else
        call mpi_reduce(x, x, size(x), mpi_double_precision, mpi_sum, 0, mpi_comm_world)
      endif
    else
      call mpi_allreduce(mpi_in_place, x, size(x), mpi_double_precision, mpi_sum, mpi_comm_world)
    end if

    call TimerStop('co_sum_double')
#endif
  end subroutine


  subroutine co_init_comms()
#ifdef TROVE_USE_MPI_
    integer :: ierr

    call mpi_init(ierr)
    if (ierr .gt. 0) stop "MPI_INIT"
    call mpi_comm_size(mpi_comm_world, comm_size, ierr)
    if (ierr .gt. 0) stop "MPI_COMM_SIZE"
    call mpi_comm_rank(mpi_comm_world, mpi_rank, ierr)
    if (ierr .gt. 0) stop "MPI_COMM_RANK"

    call MPI_Type_size(mpi_double_precision, mpi_real_size,ierr)
    call MPI_Type_size(mpi_integer, mpi_int_size,ierr)

    if (mpi_rank.ne.0) then
      open(newunit=out, file='/dev/null', status='replace', iostat=ierr, action="write")
    endif
#else
      comm_size = 1
      mpi_rank = 0
#endif
    comms_inited = .true.

    call co_init_blacs()




  end subroutine co_init_comms

  subroutine co_finalize_comms()
    integer :: ierr

    if (.not. comms_inited) stop "CO_FINALIZE_COMMS COMMS NOT INITED"

#ifdef TROVE_USE_MPI_
    call mpi_finalize(ierr)

    if (ierr .gt. 0) stop "MPI_FINALIZE"
#endif

  end subroutine co_finalize_comms


  subroutine co_init_distr(dimen, startdim, enddim, blocksize)
    integer,intent(in) :: dimen
    integer,intent(out) :: startdim, enddim, blocksize
    integer :: ierr

    integer :: proc_index, localsize
    integer :: i, to_calc, ioslice_width, ioslice_maxwidth

    if (.not. comms_inited) stop "COMMS NOT INITIALISED"
    !if (distr_inited) stop "DISTRIBUTION ALREADY INITIALISED"

    if (dimen < comm_size) then
      stop "co_init_distr: Cannot distribute matrix of dimension less than comm_size"
    endif

    if (.not. distr_inited) then
      allocate(send_or_recv(comm_size),stat=ierr)
      if (ierr .gt. 0) stop "CO_INIT_DISTR ALLOCATION FAILED"
    endif

    co_curr_dimen = dimen ! set which dimension array we're currently distributing

    if (comm_size .eq. 1) then
      co_startdim = 1
      co_enddim = dimen
      co_blocksize = dimen*dimen
      co_localsize = dimen
      send_or_recv(1) = 0
    else

      localsize = 1+dimen/comm_size
      co_startdim = mpi_rank*localsize + 1
      if (mpi_rank == comm_size-1) then
        ! Last process gets full dimension
        co_enddim = dimen
      else
        co_enddim = (mpi_rank+1)*localsize
      endif
      co_blocksize = localsize*(comm_size*localsize)

      co_localsize = localsize

      if(.not. distr_inited) then
        allocate(mpi_blocktype(comm_size))
      endif

      do i=1,comm_size
        if (mod(comm_size,2).eq.1) then
          to_calc = comm_size/2+1
        else
          if ((mod(comm_size, 4).eq.0) .and. ((mod(i,2).eq.1 .and. i.le.comm_size/2).or.(mod(i,2).eq.0 .and. i.gt.comm_size/2))) then
            to_calc = comm_size/2+1
          elseif (mod(i,2).eq.1 .and. (mod(comm_size, 4).gt.0)) then
            to_calc = comm_size/2+1
          else
            to_calc = comm_size/2
          endif
        endif

        proc_index = mpi_rank+1

        if (i.eq.proc_index) then
          send_or_recv(i) = 0
        elseif ( ((i.gt.(proc_index - to_calc) .and. i.lt.proc_index)) .or. &
            ((proc_index-to_calc).lt.1 .and. (i-comm_size).gt.(proc_index-to_calc))) then
          send_or_recv(i) = 1 ! send
#ifdef TROVE_USE_MPI_
          call co_create_type_subarray(co_localsize, co_blocksize, co_localsize, i, mpi_blocktype(i))
#endif
        else
          send_or_recv(i) = -1 ! recv
        endif
      end do
    endif

    startdim = co_startdim
    enddim = co_enddim
    blocksize = co_blocksize

    ioslice_width = co_enddim-co_startdim+1
#ifdef TROVE_USE_MPI_
    if (comm_size .eq. 1) then
      call co_create_type_column(dimen,dimen,dimen)
    else
      call co_create_type_column(ioslice_width, dimen, comm_size*co_localsize)
    endif
#endif

    distr_inited = .true.
  end subroutine co_init_distr

  subroutine co_validate_dimensions(dimen)
    integer, intent(in) :: dimen

    if (dimen .ne. co_curr_dimen) then
      stop "Tried to use an array of different dimension to the current distributed setup"
    endif
  end subroutine

  subroutine co_create_distr_array(arr, dimen)
    real(rk), pointer, intent(out) :: arr(:,:)
    integer, intent(in) :: dimen

    call co_validate_dimensions(dimen)

    allocate(arr(comm_size*co_localsize, co_startdim:co_enddim))
  end subroutine

  !
  ! Distribute the contents of an array among processes.
  ! If only using one process or not using MPI, do nothing.
  !
  subroutine co_distr_data(x, tmp, blocksize, lb, ub)

    real(rk),dimension(:,lb:),intent(inout) :: x
    real(rk),dimension(:,:,:),intent(inout) :: tmp
    integer,intent(in)                :: blocksize, lb, ub

#ifdef TROVE_USE_MPI_
    integer :: i, icoeff, jcoeff, offset, ierr, k
    type(MPI_Request)  :: reqs(comm_size)

    if (comm_size.eq.1) return

    call TimerStart('MPI_transpose')
    call TimerStart('MPI_transpose_sendrecv')

    do i=1,comm_size
      reqs(i)= MPI_REQUEST_NULL
    end do

    do i=1,comm_size
      if (send_or_recv(i).eq.1) then
        call mpi_isend(x,1,mpi_blocktype(i),i-1,0,mpi_comm_world,reqs(i),ierr)
      elseif (send_or_recv(i).eq.-1) then
        call mpi_irecv(tmp(:,:,i),blocksize*blocksize,mpi_double_precision,i-1,mpi_any_tag,mpi_comm_world,reqs(i),ierr)
      else
        reqs(i) = MPI_REQUEST_NULL
      endif
    enddo

    call mpi_waitall(comm_size,reqs,mpi_statuses_ignore,ierr)
    call TimerStop('MPI_transpose_sendrecv')
    call TimerStart('MPI_transpose_local')

    do i=1,comm_size
      if (send_or_recv(i).eq.-1) then
        offset = (i-1)*blocksize
        !$omp parallel do private(icoeff,jcoeff) shared(i,x,tmp,lb,ub,offset,blocksize) schedule(static)
        do icoeff=lb,ub
          do jcoeff=offset+1,offset+blocksize
            x(jcoeff,icoeff) = tmp(icoeff-lb+1,jcoeff-offset,i)
          enddo
        enddo
        !$omp end parallel do
      endif
    enddo
    call TimerStop('MPI_transpose_local')
    call TimerStop('MPI_transpose')
#endif

  end subroutine co_distr_data

  subroutine co_read_matrix_distr_ordered(x, longdim, lb, ub, infile)

    real(rk),dimension(:,lb:),intent(out) :: x
    integer,intent(in)                :: longdim, lb, ub
    type(MPI_File),intent(inout) :: infile

#ifdef TROVE_USE_MPI_
    type(MPI_Status) :: writestat
    integer(kind=MPI_Offset_kind) :: offset_start,offset_end
    integer :: readcount, ierr

    call mpi_barrier(mpi_comm_world, ierr)

    call TimerStart('MPI_read_matrix')

    call MPI_File_read_ordered(infile,x,1,mpitype_column,writestat,ierr)

    call TimerStop('MPI_read_matrix')
#endif

  end subroutine co_read_matrix_distr_ordered

  subroutine co_read_matrix_distr(x, longdim, lb, ub, infile)

    real(rk),dimension(:,lb:),intent(out) :: x
    integer,intent(in)                :: longdim, lb, ub
    type(MPI_File),intent(inout) :: infile

#ifdef TROVE_USE_MPI_
    type(MPI_Status) :: writestat
    integer(kind=MPI_Offset_kind) :: offset_start,offset_end
    integer :: readcount, ierr

    call mpi_barrier(mpi_comm_world, ierr)

    call TimerStart('MPI_read_matrix')

    if (comm_size.gt.1) then
      offset_start = (lb-1)*int(longdim, MPI_OFFSET_KIND)*mpi_real_size
      offset_end = (longdim-ub)*int(longdim, MPI_OFFSET_KIND)*mpi_real_size

      call MPI_File_seek(infile, offset_start, MPI_SEEK_CUR)
      call MPI_File_read_all(infile,x,1,mpitype_column,writestat,ierr)
      call MPI_File_seek(infile, offset_end, MPI_SEEK_CUR)
    else
      call MPI_File_read_all(infile,x,1,mpitype_column,writestat,ierr)
    endif

    call TimerStop('MPI_read_matrix')
#endif

  end subroutine co_read_matrix_distr

  subroutine co_write_matrix_distr(x, longdim, lb, ub, outfile)

    real(rk),dimension(:,lb:),intent(in) :: x
    integer,intent(in)                :: longdim, lb, ub
    type(MPI_File),intent(inout) :: outfile

#ifdef TROVE_USE_MPI_
    integer :: writecount, ierr
    integer(kind=MPI_Offset_kind) :: offset_start, offset_end
    type(MPI_Status) :: writestat

    call mpi_barrier(mpi_comm_world, ierr)

    call TimerStart('MPI_write_matrix')

    call MPI_File_write_ordered(outfile,x,1,mpitype_column,writestat,ierr)

    call TimerStop('MPI_write_matrix')
#endif

  end subroutine co_write_matrix_distr

#ifdef TROVE_USE_MPI_
  subroutine co_create_type_column(ncols, extent, blocksize)
    integer, intent(in) :: extent, blocksize, ncols
    integer :: ierr,writecount

    call MPI_Type_vector(ncols, extent, blocksize, mpi_double_precision, mpitype_column, ierr)
    call MPI_Type_commit(mpitype_column, ierr)

  end subroutine co_create_type_column

  subroutine co_create_type_subarray(extent, coldim, rowdim, blockid, mpi_newtype)
    integer,intent(in) :: extent, coldim, rowdim, blockid
    type(MPI_Datatype),intent(inout) :: mpi_newtype
    integer,dimension(2) :: array_of_sizes, array_of_subsizes, array_of_starts
    integer :: ierr

    array_of_sizes(1) = comm_size * extent!coldim
    array_of_sizes(2) = extent
    array_of_subsizes(:) = extent
    array_of_starts(1) = (blockid - 1) * extent
    array_of_starts(2) = 0


    call MPI_Type_create_subarray(2, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, mpi_double_precision, mpi_newtype, ierr)
    call MPI_Type_commit(mpi_newtype, ierr)

  end subroutine co_create_type_subarray
#endif

end module
