module writer_base
  implicit none

  type, abstract :: writerBase
  contains
    generic :: write => writeScalar, write1DArray, write2DArray
    procedure(writeScalar), deferred :: writeScalar
    procedure(write1DArray), deferred :: write1DArray
    procedure(write2DArray), deferred :: write2DArray
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
  end interface

  private

  public :: writerBase

end module
