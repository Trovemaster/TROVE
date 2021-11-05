#include "errors.fpp"

#define MPI_WRAPPER(function, handle, obj, size, mytype, status, err) \
select type(obj); \
    type is (integer(ik)); \
      call function(handle, obj, size, mytype, status, err); \
    type is (integer(hik)); \
      call function(handle, obj, size, mytype, status, err); \
    type is (real); \
      call function(handle, obj, size, mytype, status, err); \
    type is (real(rk)); \
      call function(handle, obj, size, mytype, status, err); \
    type is (real(ark)); \
      call function(handle, obj, size, mytype, status, err); \
    type is (complex); \
      call function(handle, obj, size, mytype, status, err); \
    type is (character(len=*)); \
      call function(handle, obj, size, mytype, status, err); \
    class default; \
      stop 'MPI: tried to handle unsupported type'; \
end select


module io_handler_mpi
  use mpi_f08
  use mpi_aux
  use io_handler_base
  use errors
  use accuracy, only : rk, ark, ik, hik

  implicit none

  type, extends(ioHandlerBase) :: ioHandlerMPI
    integer (kind=MPI_Offset_kind) :: bookendBytes = 0
    integer :: rank = -1
    type(MPI_File) :: fileh
    logical :: isOpen = .false.
    character (len=20) :: accessVal
  contains
    procedure :: writeScalar => writeScalarMPI
    procedure :: write1DArray => write1DArrayMPI
    procedure :: write2DArray => write2DArrayMPI
    procedure :: write2DArrayDistBlacs => write2DArrayDistBlacsMPI
    procedure :: write2DArrayDistColumn => write2DArrayDistColumnMPI
    procedure :: readScalar => readScalarMPI
    procedure :: read1DArray => read1DArrayMPI
    procedure :: read2DArray => read2DArrayMPI
    procedure :: read2DArrayDistBlacs => read2DArrayDistBlacsMPI
    procedure :: read2DArrayDistColumn => read2DArrayDistColumnMPI
    procedure :: open
    procedure :: close
    procedure :: seek
    final :: destroyIoHandlerMPI
  end type ioHandlerMPI

  private

  public :: ioHandlerMPI

  contains

    subroutine destroyIoHandlerMPI(this)
      type(ioHandlerMPI) :: this
      call this%close()
    end subroutine

    subroutine open(this, fname, err, action, position, status, form, access)
      ! writer MPI constructor
      class(ioHandlerMPI) :: this
      character (len = *), intent(in) :: fname
      type(ErrorType), intent(inout) :: err
      character (len = *), intent(in) :: action
      character (len = *), intent(in), optional :: position, status, form, access
      character (len = 20) :: positionVal, statusVal, formVal
      integer :: ierr

      if (this%isOpen) then
        RAISE_ERROR("ERROR: Tried to open second file", err)
      endif

      if (present(position)) then
        positionVal = position
      else
        positionVal = 'append'
      end if

      if (present(status)) then
        statusVal = status
      else
        statusVal = 'unknown'
      end if

      if (present(form)) then
        formVal = form
      else
        formVal = 'formatted'
      end if

      if (present(access)) then
        this%accessVal = access
      else
        this%accessVal = 'sequential'
      endif

      if (trim(this%accessVal) == "sequential") then
        this%bookendBytes = 4
      else
        this%bookendBytes = 0
      endif

      ! FIXME use above flags to change open behaviour

      call MPI_Comm_rank(MPI_COMM_WORLD, this%rank, ierr)

      if(this%rank == 0) then
        print *, "MPI: Opening ", trim(fname), " with ", \
          trim(action), " ", trim(positionVal), " ", trim(statusVal), " ", trim(formVal), " ", trim(this%accessVal)
      endif

      ! FIXME is there a better way to set MPI_MODE_* flags?
      if(trim(action) == 'write') then
        call MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, this%fileh, ierr)
      else if(trim(action) == 'read') then
        call MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, this%fileh, ierr)
      endif

      call MPI_File_set_errhandler(this%fileh, MPI_ERRORS_ARE_FATAL)
      ! FIXME handle error
    end subroutine open

    subroutine close(this)
      class(ioHandlerMPI) :: this
      integer :: ierr
      if (this%isOpen) then
        call MPI_File_close(this%fileh, ierr)
        this%isOpen = .false.
      endif
      ! FIXME handle error
    end subroutine close

    subroutine seek(this, offset)
      class(ioHandlerMPI) :: this
      integer, intent(in) :: offset
      integer(kind=MPI_OFFSET_KIND) :: total_offset

      if (trim(this%accessVal) == "sequential") then
        ! Add two bookend offsets
        total_offset = offset + 2*4
      endif
      call MPI_File_seek(this%fileh, total_offset, MPI_SEEK_CUR)
    end subroutine

    subroutine getMPIVarInfo(object, byteSize, mpiType)
      class(*), intent(in) :: object
      integer, intent(out) :: byteSize
      type(MPI_Datatype), intent(out) :: mpiType

      select type(object)
      type is (integer(kind=4))
        byteSize = sizeof(object)
        mpiType = MPI_INTEGER
      type is (integer(kind=8))
        byteSize = sizeof(object)
        mpiType = MPI_LONG
      type is (real(kind=4))
        byteSize = sizeof(object)
        mpiType = MPI_FLOAT
      type is (real(kind=8))
        byteSize = sizeof(object)
        mpiType = MPI_DOUBLE
      type is (complex(kind=4))
        byteSize = sizeof(object)
        mpiType = MPI_COMPLEX
      type is (complex(kind=8))
        byteSize = sizeof(object)
        mpiType = MPI_DOUBLE_COMPLEX
      type is (character(len=*))
        byteSize = len(object) * sizeof('a')
        mpiType = MPI_CHARACTER
      class default
        print *, "ERROR: Unknown type"
      end select
    end subroutine

    subroutine writeBookendBytes(this, bytes)
      ! write first and last bookends containing object size in bytes

      class(ioHandlerMPI) :: this
      integer :: bytes, ierr
      integer(kind=MPI_OFFSET_KIND) :: offset

      if (this%rank == 0 .and. trim(this%accessVal)=="sequential") then
        ! Get original position
        call MPI_File_get_position(this%fileh, offset, ierr)

        ! Write bookends
        call MPI_File_write(this%fileh, bytes, 1, MPI_INTEGER, &
                            MPI_STATUS_IGNORE, ierr)
        call MPI_File_seek(this%fileh, int(bytes,MPI_OFFSET_KIND), MPI_SEEK_CUR, ierr)
        call MPI_File_write(this%fileh, bytes, 1, MPI_INTEGER, &
                            MPI_STATUS_IGNORE, ierr)

        ! Reset position
        call MPI_File_seek(this%fileh, offset, MPI_SEEK_SET, ierr)
      endif
    end subroutine

    subroutine writeScalarMPI(this, object)
      class(ioHandlerMPI) :: this
      class(*), intent(in) :: object

      integer :: byteSize, ierr, length
      type(MPI_Datatype) :: mpiType
      integer(kind = MPI_OFFSET_KIND) :: offset

      call getMPIVarInfo(object, byteSize, mpiType)

      select type(object)
      type is (character(len=*))
        length = len(object)
      class default
        length = 1
      end select

      call writeBookendBytes(this, byteSize)

      call MPI_File_seek(this%fileh, this%bookendBytes, MPI_SEEK_CUR, ierr)
      if (this%rank == 0) then
        MPI_WRAPPER(MPI_File_write, this%fileh, object, length, mpiType, MPI_STATUS_IGNORE, ierr)
      else
        call MPI_File_seek(this%fileh, int(byteSize,kind=MPI_OFFSET_KIND), MPI_SEEK_CUR, ierr)
      end if
      call MPI_File_seek(this%fileh, this%bookendBytes, MPI_SEEK_CUR, ierr)
    end subroutine

    subroutine write1DArrayMPI(this, object)
      class(ioHandlerMPI) :: this
      class(*), intent(in) :: object(:)
      stop "ERROR: 1D array saving not currently supported by MPI writer"
    end subroutine

    subroutine write2DArrayMPI(this, object)
      class(ioHandlerMPI) :: this
      class(*), intent(in) :: object(:,:)

      type(MPI_Datatype) :: mpiType
      integer :: byteSize, globalSize, ierr, length, arrSizeBytes

      integer(kind = MPI_OFFSET_KIND) :: offset

      globalSize = size(object)

      call getMPIVarInfo(object(1,1), byteSize, mpiType)
      arrSizeBytes = globalSize*byteSize

      call writeBookendBytes(this, arrSizeBytes)

      call MPI_File_seek(this%fileh, this%bookendBytes, MPI_SEEK_CUR, ierr)
      if (this%rank == 0) then
        MPI_WRAPPER(MPI_File_write, this%fileh, object, globalSize, mpiType, MPI_STATUS_IGNORE, ierr)
      else
        call MPI_File_seek(this%fileh, int(arrSizeBytes,kind=MPI_OFFSET_KIND), MPI_SEEK_CUR, ierr)
      end if
      call MPI_File_seek(this%fileh, this%bookendBytes, MPI_SEEK_CUR, ierr)
    end subroutine

    subroutine write2DArrayDistBlacsMPI(this, object, descr, block_type)
      class(ioHandlerMPI) :: this
      class(*), intent(in) :: object(:,:)
      integer, intent(in) :: descr(9) ! Description array outputted from co_block_type_init
      type(MPI_Datatype), intent(in) :: block_type ! subarray type outputed from co_block_type_init

      type(MPI_Datatype) :: mpiType
      integer :: byteSize, arrSizeBytes, globalSize, ierr
      integer(kind = MPI_OFFSET_KIND) :: offset, disp

      integer :: dims(2)

      dims(:) = descr(3:4)
      globalSize = dims(1)*dims(2)

      call getMPIVarInfo(object(1,1), byteSize, mpiType)
      arrSizeBytes = globalSize*byteSize

      call writeBookendBytes(this, arrSizeBytes)

      call MPI_File_get_position(this%fileh, offset, ierr)
      ! Get initial displacement in file
      call MPI_File_get_byte_offset(this%fileh, offset, disp, ierr)

      ! Set file view including offsetting bookend
      call MPI_File_set_view(this%fileh, disp+this%bookendBytes, mpiType, block_type, &
                             'native', MPI_INFO_NULL, ierr)
      ! Write array in parallel
      MPI_WRAPPER(MPI_File_write_all, this%fileh, object, size(object), mpiType, MPI_STATUS_IGNORE, ierr)
      ! Reset file view
      call MPI_File_set_view(this%fileh, int(0,MPI_OFFSET_KIND), MPI_BYTE, MPI_BYTE, &
                             'native', MPI_INFO_NULL, ierr)
      call MPI_File_seek(this%fileh, disp+this%bookendBytes+arrSizeBytes+this%bookendBytes, MPI_SEEK_SET)
    end subroutine

    subroutine write2DArrayDistColumnMPI(this, object, mdimen)
      class(ioHandlerMPI) :: this
      class(*), intent(in) :: object(:,:)
      integer, intent(in) :: mdimen ! Dimension of entire distributed array

      type(MPI_Datatype) :: mpiType
      integer :: byteSize, globalSize, ierr, writestat, arrSizeBytes
      integer(kind = MPI_OFFSET_KIND) :: offset, disp

      globalSize = mdimen**2

      call getMPIVarInfo(object(1,1), byteSize, mpiType)
      arrSizeBytes = globalSize*byteSize

      call writeBookendBytes(this, arrSizeBytes)

      ! Get individual pointer offset
      call MPI_File_get_position(this%fileh, offset, ierr)
      ! Set shared pointer to individual pointer + bookend
      call MPI_File_seek_shared(this%fileh, offset+this%bookendBytes, MPI_SEEK_SET, ierr)
      ! Write array in parallel
      MPI_WRAPPER(MPI_File_write_ordered,this%fileh,object,1,mpitype_column,MPI_STATUS_IGNORE,ierr)
      ! Skip over last bookend
      call MPI_File_seek_shared(this%fileh, int(this%bookendBytes,MPI_OFFSET_KIND), MPI_SEEK_CUR, ierr)

      ! Set individual pointer to match shared
      call MPI_File_get_position_shared(this%fileh, offset, ierr)
      call MPI_File_seek(this%fileh, offset, MPI_SEEK_SET, ierr)
    end subroutine

    subroutine readScalarMPI(this, object)
      class(ioHandlerMPI) :: this
      class(*), intent(out) :: object

      integer :: byteSize, ierr, length
      type(MPI_Datatype) :: mpiType

      call getMPIVarInfo(object, byteSize, mpiType)

      select type(object)
      type is (character(len=*))
        length = len(object)
      class default
        length = 1
      end select

      call MPI_File_seek(this%fileh, int(this%bookendBytes,MPI_OFFSET_KIND), MPI_SEEK_CUR, ierr)
      MPI_WRAPPER(MPI_File_read_all, this%fileh, object, length, mpiType, MPI_STATUS_IGNORE, ierr)
      call MPI_File_seek(this%fileh, int(this%bookendBytes,MPI_OFFSET_KIND), MPI_SEEK_CUR, ierr)
    end subroutine

    subroutine read1DArrayMPI(this, object)
      class(ioHandlerMPI) :: this
      class(*), dimension(:), intent(out) :: object
      stop "reading 1D array with MPI IO not supported"
    end subroutine

    subroutine read2DArrayMPI(this, object)
      class(ioHandlerMPI) :: this
      class(*), intent(out) :: object(:,:)

      type(MPI_Datatype) :: mpiType
      integer :: byteSize, globalSize, ierr, length, arrSizeBytes

      integer(kind = MPI_OFFSET_KIND) :: offset

      globalSize = size(object)

      call getMPIVarInfo(object(1,1), byteSize, mpiType)
      arrSizeBytes = globalSize*byteSize

      call MPI_File_seek(this%fileh, this%bookendBytes, MPI_SEEK_CUR, ierr)
      MPI_WRAPPER(MPI_File_read_all, this%fileh, object, globalSize, mpiType, MPI_STATUS_IGNORE, ierr)
      call MPI_File_seek(this%fileh, this%bookendBytes, MPI_SEEK_CUR, ierr)
    end subroutine

    subroutine read2DArrayDistBlacsMPI(this, object, descr, block_type)
      class(ioHandlerMPI) :: this
      class(*), intent(out) :: object(:,:)
      integer, intent(in) :: descr(9) ! Description array outputted from co_block_type_init
      type(MPI_Datatype), intent(in) :: block_type ! subarray type outputed from co_block_type_init
      type(MPI_Datatype) :: mpiType
      integer :: byteSize, arrSizeBytes, globalSize, ierr
      integer(kind = MPI_OFFSET_KIND) :: offset, disp

      integer :: dims(2)

      dims(:) = descr(3:4)
      globalSize = dims(1)*dims(2)

      call getMPIVarInfo(object(1,1), byteSize, mpiType)
      arrSizeBytes = globalSize*byteSize

      call MPI_File_get_position(this%fileh, offset, ierr)
      ! Get initial displacement in file
      call MPI_File_get_byte_offset(this%fileh, offset, disp, ierr)
      ! Set file view including offsetting bookend
      call MPI_File_set_view(this%fileh, disp+this%bookendBytes, mpiType, block_type, &
                             'native', MPI_INFO_NULL, ierr)
      ! Read array in parallel
      MPI_WRAPPER(MPI_File_read_all, this%fileh, object, size(object), mpiType, MPI_STATUS_IGNORE, ierr)
      ! Reset file view back to regular ol bytes, including bookends and array we've just written
      call MPI_File_set_view(this%fileh, int(0,MPI_OFFSET_KIND), MPI_BYTE, MPI_BYTE, &
                             'native', MPI_INFO_NULL, ierr)
      call MPI_File_seek(this%fileh, disp+this%bookendBytes+arrSizeBytes+this%bookendBytes, MPI_SEEK_SET)
    end subroutine

    subroutine read2DArrayDistColumnMPI(this, object, dimen)
      class(ioHandlerMPI) :: this
      class(*), intent(out) :: object(:,:)
      integer, intent(in) :: dimen ! Dimension of entire distributed array

      integer :: ierr
      integer(kind = MPI_OFFSET_KIND) :: offset, disp

      ! Get individual pointer offset
      call MPI_File_get_position(this%fileh, offset, ierr)
      ! Set shared pointer to individual pointer + bookend
      call MPI_File_seek_shared(this%fileh, offset+this%bookendBytes, MPI_SEEK_SET, ierr)
      ! Read array in parallel
      MPI_WRAPPER(MPI_File_read_ordered,this%fileh,object,1,mpitype_column,MPI_STATUS_IGNORE,ierr)
      ! Skip over last bookend
      call MPI_File_seek_shared(this%fileh, int(this%bookendBytes,MPI_OFFSET_KIND), MPI_SEEK_CUR, ierr)

      ! Set individual pointer to match shared
      call MPI_File_get_position_shared(this%fileh, offset, ierr)
      call MPI_File_seek(this%fileh, offset, MPI_SEEK_SET, ierr)
    end subroutine

end module
