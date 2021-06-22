module reader_mpi
  use mpi
  use reader_type

  implicit none

  type, extends(readerType) :: readerMPI
  contains
    procedure :: readScalar => readScalar_MPI
    procedure :: read1DArray => read1DArray_MPI
    procedure :: read2DArray => read2DArray_MPI
  end type readerMPI

  private

  public :: readerMPI

  contains

    subroutine readScalar_MPI(this, object)
      class(readerMPI) :: this
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

    subroutine read1DArray_MPI(this, object)
      class(readerMPI) :: this
      class(*), dimension(:), intent(out) :: object
      print *, "reading 1D array to MPI IO"
    end subroutine

    subroutine read2DArray_MPI(this, object)
      class(readerMPI) :: this
      class(*), dimension(:,:), intent(out) :: object
      print *, "reading 2D array to MPI IO"
    end subroutine


end module
