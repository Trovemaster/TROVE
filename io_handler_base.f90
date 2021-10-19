module io_handler_base
  use mpi_aux

  implicit none

  type, abstract :: ioHandlerBase
  contains
    generic :: write => writeScalar, write1DArray, write2DArray, write2DArrayDistBlacs, write2DArrayDistColumn
    procedure(writeScalar), deferred :: writeScalar
    procedure(write1DArray), deferred :: write1DArray
    procedure(write2DArray), deferred :: write2DArray
    procedure(write2DArrayDistBlacs), deferred :: write2DArrayDistBlacs
    procedure(write2DArrayDistColumn), deferred :: write2DArrayDistColumn
    generic :: read => readScalar, read1DArray, read2DArray, read2DArrayDistBlacs
    procedure(readScalar), deferred :: readScalar
    procedure(read1DArray), deferred :: read1DArray
    procedure(read2DArray), deferred :: read2DArray
    procedure(read2DArrayDistBlacs), deferred :: read2DArrayDistBlacs
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
    subroutine write2DArrayDistBlacs(this, object, descr, block_type)
      import ioHandlerBase
      import MPI_Datatype
      class(ioHandlerBase) :: this
      class(*), dimension(:,:), intent(in) :: object
      integer, intent(in) :: descr(9)
      type(MPI_Datatype), intent(in) :: block_type
    end subroutine
    subroutine write2DArrayDistColumn(this, object, mdimen)
      import ioHandlerBase
      import MPI_Datatype
      class(ioHandlerBase) :: this
      class(*), dimension(:,:), intent(in) :: object
      integer, intent(in) :: mdimen
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
    subroutine read2DArrayDistBlacs(this, object, descr, block_type)
      import ioHandlerBase
      import MPI_Datatype
      class(ioHandlerBase) :: this
      class(*), dimension(:,:), intent(out) :: object
      integer, intent(in) :: descr(9) ! Description array outputted from co_block_type_init
      type(MPI_Datatype), intent(in) :: block_type
    end subroutine
  end interface

  private

  public :: ioHandlerBase

end module
