#include "errors.fpp"

module writer_mpi
  use mpi
  use writer_base
  use errors

  implicit none

  type, extends(writerBase) :: writerMPI
    integer (kind=MPI_Offset_kind) :: offset
    integer :: fileh, rank
    logical :: isOpen = .false.
  contains
    procedure :: writeScalar => writeScalarMPI
    procedure :: write1DArray => write1DArrayMPI
    procedure :: write2DArray => write2DArrayMPI
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

      this%isOpen = .false.
      this%offset = 0
      this%fileh = 0
      this%rank = 0

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
      integer, intent(out) :: byteSize, mpiType

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

      integer :: byteSize, mpiType, ierr

      if (this%rank /= 0) then
        return
      end if

      call getMPIVarInfo(object, byteSize, mpiType)

      this%offset = this%offset + 4+byteSize+4
      call MPI_File_write(this%fileh, byteSize, 1, MPI_INTEGER, &
                          MPI_STATUS_IGNORE, ierr)
      call MPI_File_write(this%fileh, object, 1, mpiType, &
                          MPI_STATUS_IGNORE, ierr)
      call MPI_File_write(this%fileh, byteSize, 1, MPI_INTEGER, &
                          MPI_STATUS_IGNORE, ierr)
    end subroutine

    subroutine write1DArrayMPI(this, object)
      class(writerMPI) :: this
      class(*), dimension(:), intent(in) :: object
      print *, "writing 1D array to MPI IO"
    end subroutine

    subroutine write2DArrayMPI(this, object)
      class(writerMPI) :: this
      class(*), dimension(:,:), intent(in) :: object
      print *, "writing 2D array to MPI IO"
    end subroutine

end module
