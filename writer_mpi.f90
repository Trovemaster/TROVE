#include "errors.fpp"

module writer_mpi
  use mpi_f08
  use mpi_aux
  use writer_base
  use errors

  implicit none

  type, extends(writerBase) :: writerMPI
    integer (kind=MPI_Offset_kind) :: offset = 0
    integer :: rank = 0
    type(MPI_File) :: fileh
    logical :: isOpen = .false.
  contains
    procedure :: writeScalar => writeScalarMPI
    procedure :: write1DArray => write1DArrayMPI
    procedure :: write2DArray => write2DArrayMPI
    procedure :: write2DArrayDist => write2DArrayMPIDist
    procedure :: open
    procedure :: close
    final :: destroyWriterMPI
  end type writerMPI

  interface writerMPI
    procedure :: newWriterMPI
  end interface writerMPI

  private

  public :: writerMPI

  contains

    type(writerMPI) function newWriterMPI(fname, err, position, status, form, access) result(this)
      ! writer MPI constructor
      type(ErrorType), intent(inout) :: err
      character (len = *), intent(in) :: fname
      character (len = *), intent(in), optional :: position, status, form, access

      call this%open(fname, err, position, status, form, access)
    end function

    subroutine destroyWriterMPI(this)
      type(writerMPI) :: this
      call this%close()
    end subroutine

    subroutine open(this, fname, err, position, status, form, access)
      ! writer MPI constructor
      class(writerMPI) :: this
      type(ErrorType), intent(inout) :: err
      character (len = *), intent(in) :: fname
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

      print *, "MPI: Opening ", fname, " with ", \
        positionVal, statusVal, formVal, accessVal

      ! FIXME use above flags to change open behaviour

      call MPI_Comm_rank(MPI_COMM_WORLD, this%rank, ierr)

      call MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, this%fileh, ierr)
      ! FIXME handle error
    end subroutine open

    subroutine close(this)
      class(writerMPI) :: this
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
        byteSize = 4
        mpiType = MPI_INTEGER
      type is (integer(kind=8))
        byteSize = 8
        mpiType = MPI_LONG
      type is (real(kind=4))
        byteSize = 4
        mpiType = MPI_FLOAT
      type is (real(kind=8))
        byteSize = 8
        mpiType = MPI_DOUBLE
      type is (complex(kind=4))
        byteSize = 8
        mpiType = MPI_COMPLEX
      type is (complex(kind=8))
        byteSize = 16
        mpiType = MPI_DOUBLE_COMPLEX
      class default
        print *, "ERROR: Unknown type"
      end select
    end subroutine

    subroutine writeScalarMPI(this, object)
      class(writerMPI) :: this
      class(*), intent(in) :: object

      integer :: byteSize, ierr
      type(MPI_Datatype) :: mpiType

      call getMPIVarInfo(object, byteSize, mpiType)
      this%offset = this%offset + 4+byteSize+4

      if (this%rank /= 0) then
        return
      end if

      call MPI_File_write(this%fileh, byteSize, 1, MPI_INTEGER, &
                          MPI_STATUS_IGNORE, ierr)
      call MPI_File_write(this%fileh, object, 1, mpiType, &
                          MPI_STATUS_IGNORE, ierr)
      call MPI_File_write(this%fileh, byteSize, 1, MPI_INTEGER, &
                          MPI_STATUS_IGNORE, ierr)
    end subroutine

    subroutine write1DArrayMPI(this, object)
      class(writerMPI) :: this
      class(*), intent(in) :: object(:)
      print *, "ERROR: 1D array saving not currently supported"
    end subroutine

    subroutine write2DArrayMPI(this, object)
      class(writerMPI) :: this
      class(*), intent(in) :: object(:,:)
      print *, "ERROR: Writing non-distributed array using MPI writer not supported."
    end subroutine

    subroutine write2DArrayMPIDist(this, object, descr, block_type)
      class(writerMPI) :: this
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

      !print *, globalSize, arrSizeBytes, byteSize, mpiType

      if (this%rank == 0) then
        ! write first and last bookends containing array byte size
        call MPI_File_write(this%fileh, arrSizeBytes, 1, MPI_INTEGER, &
                            MPI_STATUS_IGNORE, ierr)
        call MPI_File_seek(this%fileh, arrSizeBytes, MPI_SEEK_CUR, ierr)
        call MPI_File_write(this%fileh, arrSizeBytes, 1, MPI_INTEGER, &
                            MPI_STATUS_IGNORE, ierr)
      endif
      ! offset first bookend
      this%offset = this%offset + 4
      ! Set file view including offset
      call MPI_File_set_view(this%fileh, this%offset, mpiType, block_type, &
                             'native', MPI_INFO_NULL, ierr)
      ! Write array in parallel
      call MPI_File_write_all(this%fileh, object, size(object), mpiType, &
                              MPI_STATUS_IGNORE, ierr)
      ! Set offset and reset file view
      this%offset = this%offset + arrSizeBytes + 4
      call MPI_File_set_view(this%fileh, this%offset, MPI_BYTE, MPI_BYTE, &
                             'native', MPI_INFO_NULL, ierr)
    end subroutine

end module
