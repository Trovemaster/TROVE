module writer_type
  implicit none

  type, abstract :: writerType
  contains
    generic :: write => writeScalar, write1DArray, write2DArray
    procedure(writeScalar), deferred :: writeScalar
    procedure(write1DArray), deferred :: write1DArray
    procedure(write2DArray), deferred :: write2DArray
  end type writerType

  abstract interface
    subroutine writeScalar(this, object)
      import writerType
      class(writerType) :: this
      class(*), intent(in) :: object
    end subroutine
    subroutine write1DArray(this, object)
      import writerType
      class(writerType) :: this
      class(*), dimension(:), intent(in) :: object
    end subroutine
    subroutine write2DArray(this, object)
      import writerType
      class(writerType) :: this
      class(*), dimension(:,:), intent(in) :: object
    end subroutine
  end interface

  private

  public :: writerType

end module
