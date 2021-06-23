module writer_mpi
  use mpi
  use writer_type

  implicit none

  type, extends(writerType) :: writerMPI
    integer (kind=MPI_Offset_kind) :: offset
    integer :: fileh, rank
  contains
    procedure :: writeScalar => writeScalar_MPI
    procedure :: write1DArray => write1DArray_MPI
    procedure :: write2DArray => write2DArray_MPI
  end type writerMPI

  interface writerMPI
    procedure :: new_writerMPI
  end interface writerMPI

  private

  public :: writerMPI

  contains

    type(writerMPI) function new_writerMPI(fname, position, status, form, access)
      ! writer MPI constructor
      character (len = *), intent(in) :: fname
      character (len = *), intent(in), optional :: position, status, form, access
      character (len = 20) :: position_val, status_val, form_val, access_val
      integer :: ierr

      print *, "Creating new writerMPI!"

      if (present(position)) then
        position_val = position
      else
        position_val = 'append'
      end if

      if (present(status)) then
        status_val = status
      else
        status_val = 'unknown'
      end if

      if (present(form)) then
        form_val = form
      else
        form_val = 'formatted'
      end if

      if (present(access)) then
        access_val = access
      else
        access_val = 'sequential'
      end if

      print *, position_val, status_val, form_val, access_val

      ! FIXME use above flags to change open behaviour

      call MPI_Comm_rank(MPI_COMM_WORLD, new_writerMPI%rank, ierr)

      call MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, new_writerMPI%fileh, ierr)
      ! FIXME handle error
    end function new_writerMPI

    subroutine destroy_writerMPI(this)
      type(writerMPI) :: this
      integer :: ierr
      print *, "Closing file"
      call MPI_File_close(this%fileh, ierr)
      ! FIXME handle error
    end subroutine

    subroutine writeScalar_MPI(this, object)
      class(writerMPI) :: this
      class(*), intent(in) :: object

      integer :: byte_size, mpi_type, ierr

      if (this%rank /= 0) then
        return
      end if

      print *, "writing object to MPI IO"

      select type(object)
      type is (integer(kind=4))
        byte_size = 4
        mpi_type = MPI_INTEGER
      type is (integer(kind=8))
        byte_size = 8
        mpi_type = MPI_LONG
      type is (real(kind=4))
        byte_size = 4
        mpi_type = MPI_FLOAT
      type is (real(kind=8))
        byte_size = 8
        mpi_type = MPI_DOUBLE
      type is (complex(kind=4))
        byte_size = 8
        mpi_type = MPI_COMPLEX
      type is (complex(kind=8))
        byte_size = 16
        mpi_type = MPI_DOUBLE_COMPLEX
      class default
        print *, "ERROR: Tried to write unsupported type"
        return
      end select

      this%offset = this%offset + 4+byte_size+4
      call MPI_File_write(this%fileh, byte_size, 1, MPI_INTEGER, &
                          MPI_STATUS_IGNORE, ierr)
      call MPI_File_write(this%fileh, object, 1, mpi_type, &
                          MPI_STATUS_IGNORE, ierr)
      call MPI_File_write(this%fileh, byte_size, 1, MPI_INTEGER, &
                          MPI_STATUS_IGNORE, ierr)
    end subroutine

    subroutine write1DArray_MPI(this, object)
      class(writerMPI) :: this
      class(*), dimension(:), intent(in) :: object
      print *, "writing 1D array to MPI IO"
    end subroutine

    subroutine write2DArray_MPI(this, object)
      class(writerMPI) :: this
      class(*), dimension(:,:), intent(in) :: object
      print *, "writing 2D array to MPI IO"
    end subroutine

end module
