#include "errors.fpp"

module writer_ftn
  use writer_base
  use errors

  implicit none

  type, extends(writerBase) :: writerFTN
    integer :: iounit = 0
    integer :: stat = 0
    logical :: isOpen = .false.
  contains
    procedure :: writeScalar => writeScalarFTN
    procedure :: write1DArray => write1DArrayFTN
    procedure :: write2DArray => write2DArrayFTN
    procedure :: open
    procedure :: close
    final :: destroyWriterFTN
  end type writerFTN

  ! Constructor
  interface writerFTN
    procedure :: newWriterFTN
  end interface writerFTN

  private

  public :: writerFTN

  contains

    type(writerFTN) function newWriterFTN(fname, err, position, status, form, access) result(this)
      ! writer FTN constructor
      type(ErrorType), intent(inout) :: err
      character (len = *), intent(in) :: fname
      character (len = *), intent(in), optional :: position, status, form, access

      this%isOpen = .false.
      this%stat = 0
      this%iounit = 0

      call this%open(fname, err, position, status, form, access)
    end function

    subroutine destroyWriterFTN(this)
      type(writerFTN) :: this
      call this%close()
    end subroutine

    subroutine open(this, fname, err, position, status, form, access)
      class(writerFTN) :: this
      type(ErrorType), intent(inout) :: err
      character (len = *), intent(in) :: fname
      character (len = *), intent(in), optional :: position, status, form, access
      character (len = 20) :: positionVal, statusVal, formVal, accessVal

      if (this%isOpen) then
        RAISE_ERROR("ERROR: Tried to open second file", err)
      endif

      if (present(position)) then
        positionVal = position
      else
        positionVal = 'asis'
      end if

      if (present(status)) then
        statusVal = status
      else
        statusVal = 'unknown'
      end if

      if (present(form)) then
        formVal = form
      else
        formVal = 'unformatted'
      end if

      if (present(access)) then
        accessVal = access
      else
        accessVal = 'sequential'
      end if

      print *, "Opening ", fname, " with ", \
        positionVal, statusVal, formVal, accessVal

      open(newunit=this%iounit, action='write',\
        form=formVal, position=positionVal, status=statusVal, file=fname,\
        iostat=this%stat)

      if (this%stat == 0) then
        this%isOpen = .true.
      else
        RAISE_ERROR("ERROR: Could not open file", err)
      endif
    end subroutine

    subroutine close(this)
      class(writerFTN) :: this
      if (this%isOpen) then
        close(this%iounit)
        this%isOpen = .false.
      endif
    end subroutine

    subroutine writeScalarFTN(this, object)
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
        print *, "ERROR: Tried to write unsupported type"
      end select
    end subroutine

    subroutine write1DArrayFTN(this, object)
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
        print *, "ERROR: Tried to write unsupported type"
      end select
    end subroutine

    subroutine write2DArrayFTN(this, object)
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
        print *, "ERROR: Tried to write unsupported type"
      end select
    end subroutine
end module
