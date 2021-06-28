module reader_base
  implicit none

  type, abstract :: readerBase
  contains
    generic :: read => readScalar, read1DArray, read2DArray
    procedure(readScalar), deferred :: readScalar
    procedure(read1DArray), deferred :: read1DArray
    procedure(read2DArray), deferred :: read2DArray
  end type readerBase

  abstract interface
    subroutine readScalar(this, object)
      import readerBase
      class(readerBase) :: this
      class(*), intent(out) :: object
    end subroutine
    subroutine read1DArray(this, object)
      import readerBase
      class(readerBase) :: this
      class(*), dimension(:), intent(out) :: object
    end subroutine
    subroutine read2DArray(this, object)
      import readerBase
      class(readerBase) :: this
      class(*), dimension(:,:), intent(out) :: object
    end subroutine
  end interface

  private

  public :: readerBase

end module
