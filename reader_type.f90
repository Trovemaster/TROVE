module reader_type
  implicit none

  type, abstract :: readerType
    integer :: iounit
  contains
    generic :: read => readScalar, read1DArray, read2DArray
    procedure(readScalar), deferred :: readScalar
    procedure(read1DArray), deferred :: read1DArray
    procedure(read2DArray), deferred :: read2DArray
  end type readerType

  abstract interface
    subroutine readScalar(this, object)
      import readerType
      class(readerType) :: this
      class(*), intent(out) :: object
    end subroutine
    subroutine read1DArray(this, object)
      import readerType
      class(readerType) :: this
      class(*), dimension(:), intent(out) :: object
    end subroutine
    subroutine read2DArray(this, object)
      import readerType
      class(readerType) :: this
      class(*), dimension(:,:), intent(out) :: object
    end subroutine
  end interface

  private

  public :: readerType

end module
