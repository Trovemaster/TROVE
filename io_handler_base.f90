module io_handler_base
  use mpi_aux

  implicit none

  type, abstract :: ioHandlerBase
  contains
    generic :: write => writeScalar, write1DArray, write2DArray, write2DArrayDist
    procedure(writeScalar), deferred :: writeScalar
    procedure(write1DArray), deferred :: write1DArray
    procedure(write2DArray), deferred :: write2DArray
    procedure(write2DArrayDist), deferred :: write2DArrayDist
    generic :: read => readScalar, read1DArray, read2DArray
    procedure(readScalar), deferred :: readScalar
    procedure(read1DArray), deferred :: read1DArray
    procedure(read2DArray), deferred :: read2DArray
  end type ioHandlerBase

  abstract interface
    subroutine writeScalar(this, object)
      import ioHandlerBase
      class(ioHandlerBase) :: this
      class(*), intent(in) :: object
    end subroutine
    subroutine write1DArray(this, object)
      import ioHandlerBase
      class(ioHandlerBase) :: this
      class(*), dimension(:), intent(in) :: object
    end subroutine
    subroutine write2DArray(this, object)
      import ioHandlerBase
      class(ioHandlerBase) :: this
      class(*), dimension(:,:), intent(in) :: object
    end subroutine
    subroutine write2DArrayDist(this, object, descr, block_type)
      import ioHandlerBase
      import MPI_Datatype
      class(ioHandlerBase) :: this
      class(*), dimension(:,:), intent(in) :: object
      integer, intent(in) :: descr(9)
      type(MPI_Datatype), intent(in) :: block_type
    end subroutine
    subroutine readScalar(this, object)
      import ioHandlerBase
      class(ioHandlerBase) :: this
      class(*), intent(out) :: object
    end subroutine
    subroutine read1DArray(this, object)
      import ioHandlerBase
      class(ioHandlerBase) :: this
      class(*), dimension(:), intent(out) :: object
    end subroutine
    subroutine read2DArray(this, object)
      import ioHandlerBase
      class(ioHandlerBase) :: this
      class(*), dimension(:,:), intent(out) :: object
    end subroutine
  end interface

  private

  public :: ioHandlerBase

end module
