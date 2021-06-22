module reader_ftn
  use reader_type

  implicit none

  type, extends(readerType) :: readerFTN
  contains
    procedure :: readScalar => readScalar_FTN
    procedure :: read1DArray => read1DArray_FTN
    procedure :: read2DArray => read2DArray_FTN
    final :: destroy_readerFTN
  end type readerFTN

  interface readerFTN
    procedure :: new_readerFTN
  end interface readerFTN

  private

  public :: readerFTN

  contains

    type(readerFTN) function new_readerFTN(fname, form, access)
      ! reader FTN constructor
      character (len = *), intent(in) :: fname
      character (len = *), intent(in), optional :: form, access
      character (len = 20) :: form_val, access_val

      print *, "Creating new readerFTN!"

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

      print *, form_val, access_val

      open(newunit=new_readerFTN%iounit, action='read', access=access_val, form=form_val, file=fname)
    end function new_readerFTN

    subroutine destroy_readerFTN(this)
      type(readerFTN) :: this
      print *, "Closing file"
      close(this%iounit)
    end subroutine

    subroutine readScalar_FTN(this, object)
      class(readerFTN) :: this
      class(*), intent(out) :: object
      print *, "reading object with FTN IO"
      select type(object)
      type is (integer)
        read(this%iounit) object
      type is (real)
        read(this%iounit) object
      type is (complex)
        read(this%iounit) object
      class default
        print *, "Unsupported type!"
      end select
    end subroutine

    subroutine read1DArray_FTN(this, object)
      class(readerFTN) :: this
      class(*), dimension(:), intent(out) :: object
      print *, "reading 1D array with FTN IO"
      select type(object)
      type is (integer)
        read(this%iounit) object
      type is (real)
        read(this%iounit) object
      type is (complex)
        read(this%iounit) object
      class default
        print *, "Unsupported type!"
      end select
    end subroutine

    subroutine read2DArray_FTN(this, object)
      class(readerFTN) :: this
      class(*), dimension(:,:), intent(out) :: object
      print *, "reading 2D array with FTN IO"
      select type(object)
      type is (integer)
        read(this%iounit) object
      type is (real)
        read(this%iounit) object
      type is (complex)
        read(this%iounit) object
      class default
        print *, "Unsupported type!"
      end select
    end subroutine
end module
