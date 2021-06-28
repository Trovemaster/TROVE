module reader_mpi
  use mpi
  use reader_base

  implicit none

  type, extends(readerBase) :: readerMPI
  contains
    procedure :: readScalar => readScalarMPI
    procedure :: read1DArray => read1DArrayMPI
    procedure :: read2DArray => read2DArrayMPI
  end type readerMPI

  private

  public :: readerMPI

  contains

    subroutine readScalarMPI(this, object)
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

    subroutine read1DArrayMPI(this, object)
      class(readerMPI) :: this
      class(*), dimension(:), intent(out) :: object
      print *, "reading 1D array to MPI IO"
    end subroutine

    subroutine read2DArrayMPI(this, object)
      class(readerMPI) :: this
      class(*), dimension(:,:), intent(out) :: object
      print *, "reading 2D array to MPI IO"
    end subroutine


end module
