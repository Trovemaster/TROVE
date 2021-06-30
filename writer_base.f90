module writer_base
  use mpi_aux

  implicit none

  type, abstract :: writerBase
  contains
    generic :: write => writeScalar, write1DArray, write2DArray, write2DArrayDist
    procedure(writeScalar), deferred :: writeScalar
    procedure(write1DArray), deferred :: write1DArray
    procedure(write2DArray), deferred :: write2DArray
    procedure(write2DArrayDist), deferred :: write2DArrayDist
  end type writerBase

  abstract interface
    subroutine writeScalar(this, object)
      import writerBase
      class(writerBase) :: this
      class(*), intent(in) :: object
    end subroutine
    subroutine write1DArray(this, object)
      import writerBase
      class(writerBase) :: this
      class(*), dimension(:), intent(in) :: object
    end subroutine
    subroutine write2DArray(this, object)
      import writerBase
      class(writerBase) :: this
      class(*), dimension(:,:), intent(in) :: object
    end subroutine
    subroutine write2DArrayDist(this, object, descr, block_type)
      import writerBase
      import MPI_Datatype
      class(writerBase) :: this
      class(*), dimension(:,:), intent(in) :: object
      integer, dimension(9), intent(in) :: descr
      type(MPI_Datatype), intent(in) :: block_type
    end subroutine
  end interface

  private

  public :: writerBase

end module
