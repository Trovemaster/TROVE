module writer_ftn
  use writer_type

  implicit none

  type, extends(writerType) :: writerFTN
    integer :: iounit
  contains
    procedure :: writeScalar => writeScalar_FTN
    procedure :: write1DArray => write1DArray_FTN
    procedure :: write2DArray => write2DArray_FTN
    final :: destroy_writerFTN
  end type writerFTN

  interface writerFTN
    procedure :: new_writerFTN
  end interface writerFTN

  private

  public :: writerFTN

  contains

    type(writerFTN) function new_writerFTN(fname, position, status, form, access)
      ! writer FTN constructor
      character (len = *), intent(in) :: fname
      character (len = *), intent(in), optional :: position, status, form, access
      character (len = 20) :: position_val, status_val, form_val, access_val

      print *, "Creating new writerFTN!"

      if (present(position)) then
        position_val = position
      else
        position_val = 'append'
      end if

      if (present(status)) then
        status_val = status
      else
        status_val = 'unknown'
      end if

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

      print *, position_val, status_val, form_val, access_val

      open(newunit=new_writerFTN%iounit, action='write', form=form_val, position=position_val, status=status_val, file=fname)
    end function new_writerFTN

    subroutine destroy_writerFTN(this)
      type(writerFTN) :: this
      print *, "Closing file"
      close(this%iounit)
    end subroutine

    subroutine writeScalar_FTN(this, object)
      class(writerFTN) :: this
      class(*), intent(in) :: object
      print *, "writing object with FTN IO"
      select type(object)
      type is (integer)
        write(this%iounit) object
      type is (real)
        write(this%iounit) object
      type is (complex)
        write(this%iounit) object
      class default
        print *, "Unsupported type!"
      end select
    end subroutine

    subroutine write1DArray_FTN(this, object)
      class(writerFTN) :: this
      class(*), dimension(:), intent(in) :: object
      print *, "writing 1D array with FTN IO"
      select type(object)
      type is (integer)
        write(this%iounit) object
      type is (real)
        write(this%iounit) object
      type is (complex)
        write(this%iounit) object
      class default
        print *, "Unsupported type!"
      end select
    end subroutine

    subroutine write2DArray_FTN(this, object)
      class(writerFTN) :: this
      class(*), dimension(:,:), intent(in) :: object
      print *, "writing 2D array with FTN IO"
      select type(object)
      type is (integer)
        write(this%iounit) object
      type is (real)
        write(this%iounit) object
      type is (complex)
        write(this%iounit) object
      class default
        print *, "Unsupported type!"
      end select
    end subroutine
end module
