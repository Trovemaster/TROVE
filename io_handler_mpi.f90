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
      write(*,*) 'Not covered'; \
end select


module io_handler_mpi
  use mpi_f08
  use mpi_aux
  use io_handler_base
  use errors
  use accuracy, only : rk, ark, ik, hik

  implicit none

  type, extends(ioHandlerBase) :: ioHandlerMPI
    integer (kind=MPI_Offset_kind) :: offset = 0
    integer :: rank = 0
    type(MPI_File) :: fileh
    logical :: isOpen = .false.
  contains
    procedure :: writeScalar => writeScalarMPI
    procedure :: write1DArray => write1DArrayMPI
    procedure :: write2DArray => write2DArrayMPI
    procedure :: write2DArrayDistBlacs => write2DArrayDistBlacsMPI
    procedure :: write2DArrayDistColumn => write2DArrayDistColumnMPI
    procedure :: readScalar => readScalarMPI
    procedure :: read1DArray => read1DArrayMPI
    procedure :: read2DArray => read2DArrayMPI
    procedure :: open
    procedure :: close
    final :: destroyIoHandlerMPI
  end type ioHandlerMPI

  interface ioHandlerMPI
    procedure :: newIoHandlerMPI
  end interface ioHandlerMPI

  private

  public :: ioHandlerMPI

  contains

    type(ioHandlerMPI) function newIoHandlerMPI(fname, err, action, position, status, form, access) result(this)
      ! writer MPI constructor
      character (len = *), intent(in) :: fname
      type(ErrorType), intent(inout) :: err
      character (len = *), intent(in) :: action
      character (len = *), intent(in), optional :: position, status, form, access

      call this%open(fname, err, action, position, status, form, access)
    end function

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
      character (len = 20) :: positionVal, statusVal, formVal, accessVal
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
        accessVal = access
      else
        accessVal = 'sequential'
      end if

      print *, "MPI: Opening ", trim(fname), " with ", \
        trim(positionVal), " ", trim(statusVal), " ", trim(formVal), " ", trim(accessVal)

      ! FIXME use above flags to change open behaviour

      call MPI_Comm_rank(MPI_COMM_WORLD, this%rank, ierr)

      ! FIXME is there a better way to set MPI_MODE_* flags?
      if(trim(action) == 'write') then
        call MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, this%fileh, ierr)
      else
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

    subroutine writeScalarMPI(this, object)
      class(ioHandlerMPI) :: this
      class(*), intent(in) :: object

      integer :: byteSize, ierr, length
      type(MPI_Datatype) :: mpiType

      call getMPIVarInfo(object, byteSize, mpiType)
      this%offset = this%offset + 4+byteSize+4

      select type(object)
      type is (character(len=*))
        length = len(object)
      class default
        length = 1
      end select

      if (this%rank /= 0) then
        return
      end if

      call MPI_File_write(this%fileh, byteSize, 1, MPI_INTEGER, &
                          MPI_STATUS_IGNORE, ierr)
      MPI_WRAPPER(MPI_File_write, this%fileh, object, length, mpiType, MPI_STATUS_IGNORE, ierr)
      call MPI_File_write(this%fileh, byteSize, 1, MPI_INTEGER, &
                          MPI_STATUS_IGNORE, ierr)
    end subroutine

    subroutine write1DArrayMPI(this, object)
      class(ioHandlerMPI) :: this
      class(*), intent(in) :: object(:)
      print *, "ERROR: 1D array saving not currently supported"
    end subroutine

    subroutine write2DArrayMPI(this, object)
      class(ioHandlerMPI) :: this
      class(*), intent(in) :: object(:,:)

      type(MPI_Datatype) :: mpiType
      integer :: byteSize, globalSize, ierr, length

      integer(kind = MPI_OFFSET_KIND) :: arrSizeBytes

      globalSize = size(object)

      call getMPIVarInfo(object(1,1), byteSize, mpiType)
      arrSizeBytes = globalSize*byteSize

      this%offset = this%offset + 4 + arrSizeBytes + 4

      if (this%rank /= 0) then
        return
      end if

      call MPI_File_write(this%fileh, arrSizeBytes, 1, MPI_INTEGER, &
                          MPI_STATUS_IGNORE, ierr)
      MPI_WRAPPER(MPI_File_write, this%fileh, object, globalSize, mpiType, MPI_STATUS_IGNORE, ierr)
      call MPI_File_write(this%fileh, arrSizeBytes, 1, MPI_INTEGER, &
                          MPI_STATUS_IGNORE, ierr)
    end subroutine

    subroutine write2DArrayDistBlacsMPI(this, object, descr, block_type)
      class(ioHandlerMPI) :: this
      class(*), intent(in) :: object(:,:)
      integer, intent(in) :: descr(9) ! Description array outputted from co_block_type_init
      type(MPI_Datatype), intent(in) :: block_type ! subarray type outputed from co_block_type_init

      type(MPI_Datatype) :: mpiType
      integer :: byteSize, globalSize, ierr
      integer(kind = MPI_OFFSET_KIND) :: arrSizeBytes

      integer :: dims(2)

      dims(:) = descr(3:4)
      globalSize = dims(1)*dims(2)

      call getMPIVarInfo(object(1,1), byteSize, mpiType)
      arrSizeBytes = globalSize*byteSize

      if (this%rank == 0) then
        ! write first and last bookends containing array size in bytes
        call MPI_File_write(this%fileh, arrSizeBytes, 1, MPI_INTEGER, &
                            MPI_STATUS_IGNORE, ierr)
        call MPI_File_seek(this%fileh, arrSizeBytes, MPI_SEEK_CUR, ierr)
        call MPI_File_write(this%fileh, arrSizeBytes, 1, MPI_INTEGER, &
                            MPI_STATUS_IGNORE, ierr)
      endif
      ! Offset first bookend
      this%offset = this%offset + 4
      ! Set file view including offset
      call MPI_File_set_view(this%fileh, this%offset, mpiType, block_type, &
                             'native', MPI_INFO_NULL, ierr)
      ! Write array in parallel
      MPI_WRAPPER(MPI_File_write_all, this%fileh, object, size(object), mpiType, MPI_STATUS_IGNORE, ierr)
      ! Offset by size of array and end bookend integer
      this%offset = this%offset + arrSizeBytes + 4
      ! Reset file view back to regular ol bytes
      call MPI_File_set_view(this%fileh, this%offset, MPI_BYTE, MPI_BYTE, &
                             'native', MPI_INFO_NULL, ierr)
    end subroutine

    subroutine write2DArrayDistColumnMPI(this, object, mdimen)
      class(ioHandlerMPI) :: this
      class(*), intent(in) :: object(:,:)
      integer, intent(in) :: mdimen ! Dimension of entire distributed array

      type(MPI_Datatype) :: mpiType
      integer :: byteSize, globalSize, ierr, writestat
      integer(kind = MPI_OFFSET_KIND) :: arrSizeBytes

      globalSize = mdimen**2

      call getMPIVarInfo(object(1,1), byteSize, mpiType)
      arrSizeBytes = globalSize*byteSize

      ! TODO what if format isn't sequential??
      if (this%rank == 0) then
        ! write first and last bookends containing array size in bytes
        call MPI_File_write(this%fileh, arrSizeBytes, 1, MPI_INTEGER, &
                            MPI_STATUS_IGNORE, ierr)
        call MPI_File_seek(this%fileh, arrSizeBytes, MPI_SEEK_CUR, ierr)
        call MPI_File_write(this%fileh, arrSizeBytes, 1, MPI_INTEGER, &
                            MPI_STATUS_IGNORE, ierr)
      endif
      ! Offset first bookend
      this%offset = this%offset + 4
      ! Seek to byte after bookend
      call MPI_File_seek_shared(this%fileh, this%offset, MPI_SEEK_SET, ierr)
      ! Write array in parallel
      MPI_WRAPPER(MPI_File_write_ordered,this%fileh,object,1,mpitype_column,MPI_STATUS_IGNORE,ierr)
      ! Offset by size of array and end bookend integer
      this%offset = this%offset + arrSizeBytes + 4
      ! Ensure all file pointers point to end of array
      call MPI_File_seek(this%fileh, this%offset, MPI_SEEK_SET, ierr)
    end subroutine

    subroutine readScalarMPI(this, object)
      class(ioHandlerMPI) :: this
      class(*), intent(out) :: object
      print *, "reading object to MPI IO"
      ! Example object handling
      select type(object)
      type is (integer)
        print *, object
      type is (real)
        print *, object
      type is (complex)
        print *, object
      class default
        print *, "Unsupported type!"
      end select
    end subroutine

    subroutine read1DArrayMPI(this, object)
      class(ioHandlerMPI) :: this
      class(*), dimension(:), intent(out) :: object
      print *, "reading 1D array to MPI IO"
    end subroutine

    subroutine read2DArrayMPI(this, object)
      class(ioHandlerMPI) :: this
      class(*), dimension(:,:), intent(out) :: object
      print *, "reading 2D array to MPI IO"
    end subroutine

end module
