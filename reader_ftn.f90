#include "errors.fpp"

module reader_ftn
  use reader_base
  use errors

  implicit none

  type, extends(readerBase) :: readerFTN
    integer :: iounit = 0
    integer :: stat = 0
    logical :: isOpen = .false.
  contains
    procedure :: readScalar => readScalarFTN
    procedure :: read1DArray => read1DArrayFTN
    procedure :: read2DArray => read2DArrayFTN
    procedure :: open
    procedure :: close
    final :: destroyReaderFTN
  end type readerFTN

  interface readerFTN
    procedure :: newReaderFTN
  end interface readerFTN

  private

  public :: readerFTN

  contains

    type(readerFTN) function newReaderFTN(fname, err, position, status, form, access) result(this)
      ! reader FTN constructor
      type(ErrorType), intent(inout) :: err
      character (len = *), intent(in) :: fname
      character (len = *), intent(in), optional :: position, status, form, access

      this%isOpen = .false.
      this%stat = 0
      this%iounit = 0

      call this%open(fname, err, position, status, form, access)
    end function

    subroutine destroyReaderFTN(this)
      type(readerFTN) :: this
      call this%close()
    end subroutine

    subroutine open(this, fname, err, position, status, form, access)
      class(readerFTN) :: this
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

      open(newunit=this%iounit, action='read',\
        form=formVal, position=positionVal, status=statusVal, file=fname,\
        iostat=this%stat)

      if (this%stat == 0) then
        this%isOpen = .true.
      else
        RAISE_ERROR("ERROR: Could not open file", err)
      endif
    end subroutine

    subroutine close(this)
      class(readerFTN) :: this
      if (this%isOpen) then
        close(this%iounit)
        this%isOpen = .false.
      endif
    end subroutine

    subroutine readScalarFTN(this, object)
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

    subroutine read1DArrayFTN(this, object)
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

    subroutine read2DArrayFTN(this, object)
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
