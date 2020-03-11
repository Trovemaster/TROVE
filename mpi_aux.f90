module mpi_aux
  use mpi_f08
  use timer
  use accuracy
  implicit none

  public co_init_comms, co_finalize_comms, co_init_distr, co_distr_data, co_write_matrix_distr, co_read_matrix_distr
  public co_block_type_init

  public send_or_recv, comm_size, mpi_rank
  public co_startdim, co_enddim

  public blacs_size, blacs_rank, blacs_ctxt
  public nprow,npcol,myprow,mypcol
  public mpi_real_size, mpi_int_size

  interface co_sum
    module procedure :: co_sum_double
  end interface

  integer,dimension(:),allocatable  :: proc_sizes, proc_offsets, send_or_recv
  integer                           :: comm_size, mpi_rank
  integer                           :: co_startdim, co_enddim
  logical                           :: comms_inited = .false., distr_inited=.false.
  type(MPI_Datatype) :: mpitype_column
  type(MPI_Datatype),dimension(:), allocatable :: mpi_blocktype
  integer                           :: mpi_real_size, mpi_int_size

  !blacs/pblas
  integer :: blacs_size, blacs_rank, blacs_ctxt
  integer :: nprow,npcol,myprow,mypcol
  integer,dimension(2)  :: blacs_dims

contains

  subroutine co_init_blacs()
    implicit none

    if (.not. comms_inited) stop "CO_INIT_BLACS COMMS NOT INITED"

    ! Must be initialised to zero - if stack contains garbage here MPI_Dims_create WILL fail
    blacs_dims = 0

    call blacs_pinfo(blacs_rank, blacs_size)
    if (blacs_rank .lt. 0) return

    call MPI_Dims_create(blacs_size, 2, blacs_dims)

    call blacs_get(-1, 0, blacs_ctxt)
    call blacs_gridinit(blacs_ctxt, 'R', blacs_dims(1), blacs_dims(2))
    call blacs_gridinfo(blacs_ctxt, nprow, npcol, myprow, mypcol)

    !write(*,"('BLACS: [',i2,',',i2'](',i4,i4,i4,i4',)')") mpi_rank,blacs_rank,nprow,npcol,myprow,mypcol
  end subroutine co_init_blacs

  subroutine co_block_type_init(smat, dimx, dimy, descr, allocinfo, mpi_type)
    implicit none

    real(rk),intent(out),dimension(:,:),allocatable   :: smat

    integer,intent(in)                                :: dimx, dimy
    integer,intent(out),dimension(9)                  :: descr
    integer,intent(out)                               :: allocinfo

    type(MPI_Datatype),intent(out),optional           :: mpi_type

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

  end subroutine co_block_type_init

  subroutine co_sum_double(x, root_process)

    real(rk), intent(inout), dimension(..) :: x
    integer, optional :: root_process

    if (comm_size.eq.1) return
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
  end subroutine


  subroutine co_init_comms()
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

    comms_inited = .true.

    call co_init_blacs()

  end subroutine co_init_comms

  subroutine co_finalize_comms()
    integer :: ierr

    if (.not. comms_inited) stop "CO_FINALIZE_COMMS COMMS NOT INITED"

    call mpi_finalize(ierr)

    if (ierr .gt. 0) stop "MPI_FINALIZE"

  end subroutine co_finalize_comms

  subroutine co_init_distr(dimen, startdim, enddim, blocksize)
    integer,intent(in) :: dimen
    integer,intent(out) :: startdim, enddim, blocksize
    integer,dimension(:),allocatable  :: starts, ends
    integer :: localsize, proc_index, localsize_
    integer :: i, ierr, to_calc

    if (.not. comms_inited) stop "COMMS NOT INITIALISED"
    !if (distr_inited) stop "DISTRIBUTION ALREADY INITIALISED"

    proc_index = mpi_rank+1

    if (.not. distr_inited) then
      allocate(proc_sizes(comm_size),proc_offsets(comm_size),send_or_recv(comm_size),starts(comm_size),ends(comm_size),stat=ierr)
      if (ierr .gt. 0) stop "CO_INIT_DISTR ALLOCATION FAILED"
    else
      allocate(starts(comm_size),ends(comm_size),stat=ierr)
    endif

    if (comm_size .eq. 1) then
      startdim = 1
      enddim = dimen
      co_startdim = 1
      co_enddim = dimen
      blocksize = dimen*dimen
      send_or_recv(1) = 0
    else

      if (mpi_rank .eq. 0) then !root

        localsize = dimen/comm_size
        localsize_ = int(1+real(dimen/comm_size))

        starts(1) = 1
        ends(1) = localsize_
        proc_sizes(1) = localsize_*(comm_size*localsize_)
        proc_offsets(1) = 0

        do i=2,comm_size-1
          starts(i) = (i-1)*localsize_+1
          ends(i) = i*localsize_
          proc_sizes(i) = localsize_ * (comm_size*localsize_)!dimen
          proc_offsets(i) = localsize_*(i-1)*(comm_size*localsize_)!dimen
        end do

        starts(comm_size) = (i-1) * localsize_ + 1
        ends(comm_size) = dimen!comm_size*localsize_!dimen
        proc_sizes(comm_size) = localsize_*comm_size*localsize_!dimen

        proc_offsets(comm_size) = (comm_size-1)*localsize_*(comm_size*localsize_)!dimen
      endif

      call mpi_bcast(starts, comm_size, mpi_integer, 0, mpi_comm_world)
      call mpi_bcast(ends, comm_size, mpi_integer, 0, mpi_comm_world)
      call mpi_bcast(proc_sizes, comm_size, mpi_integer, 0, mpi_comm_world)
      call mpi_bcast(proc_offsets, comm_size, mpi_integer, 0, mpi_comm_world)



      blocksize = proc_sizes(proc_index)
      startdim = starts(proc_index)
      enddim = ends(proc_index)

      co_startdim = startdim
      co_enddim = enddim

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


        if (i.eq.proc_index) then
          send_or_recv(i) = 0
        elseif ( ((i.gt.(proc_index - to_calc) .and. i.lt.proc_index)) .or. &
            ((proc_index-to_calc).lt.1 .and. (i-comm_size).gt.(proc_index-to_calc))) then
          send_or_recv(i) = 1 ! send
          call co_create_type_subarray(int(1+real(dimen/comm_size)), blocksize, int(1+real(dimen/comm_size)), i, mpi_blocktype(i))
        else
          send_or_recv(i) = -1 ! recv
        endif
      end do

    endif

    if (comm_size .eq. 1) then
      call co_create_type_column(dimen,dimen,dimen)
    else
      call co_create_type_column(dimen,comm_size*(int(1+real(dimen/comm_size))),enddim-startdim+1)
    endif

    deallocate(starts,ends)

    distr_inited = .true.
  end subroutine co_init_distr

  subroutine co_distr_data(x, tmp, blocksize, lb, ub)

    real(rk),dimension(:,lb:),intent(inout) :: x
    real(rk),dimension(:,:,:),intent(inout) :: tmp
    integer,intent(in)                :: blocksize, lb, ub

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

  end subroutine co_distr_data

  subroutine co_read_matrix_distr(x, longdim, lb, ub, infile)

    real(rk),dimension(:,lb:),intent(out) :: x
    integer,intent(in)                :: longdim, lb, ub

    type(MPI_File),intent(inout) :: infile
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
      call MPI_File_read(infile,x,1,mpitype_column,writestat,ierr)
    endif

    call TimerStop('MPI_read_matrix')

  end subroutine co_read_matrix_distr

  subroutine co_write_matrix_distr(x, longdim, lb, ub, outfile)

    real(rk),dimension(:,lb:),intent(in) :: x
    integer,intent(in)                :: longdim, lb, ub
    type(MPI_File),intent(inout) :: outfile
    integer :: ierr
    integer(kind=MPI_Offset_kind) :: offset_start, offset_end
    type(MPI_Status) :: writestat

    call mpi_barrier(mpi_comm_world, ierr)

    call TimerStart('MPI_write_matrix')

    if (comm_size.gt.1) then
      offset_start = (lb-1)*int(longdim,MPI_OFFSET_KIND)*mpi_real_size
      offset_end = 0
      call mpi_barrier(mpi_comm_world, ierr)

      call MPI_File_seek(outfile, offset_start, MPI_SEEK_END)
      call MPI_File_write_all(outfile,x,1,mpitype_column,writestat,ierr)
      call mpi_barrier(mpi_comm_world, ierr)
      call MPI_File_seek(outfile, offset_end, MPI_SEEK_END)
    else
      call MPI_File_write(outfile,x,1,mpitype_column,writestat,ierr)
    endif

    call TimerStop('MPI_write_matrix')

  end subroutine co_write_matrix_distr

  subroutine co_create_type_column(extent, blocksize, ncols)
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
    array_of_starts(1) = (blockid - 1) * extent + 0
    array_of_starts(2) = 0


    call MPI_Type_create_subarray(2, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, mpi_double_precision, mpi_newtype, ierr)
    call MPI_Type_commit(mpi_newtype, ierr)

  end subroutine co_create_type_subarray

end module
