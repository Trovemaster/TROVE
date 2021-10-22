#include "errors.fpp"

module io_handler_ftn
  use mpi_aux
  use io_handler_base
  use errors
  use iso_fortran_env

  implicit none

  type, extends(ioHandlerBase) :: ioHandlerFTN
    integer :: iounit = 0
    integer :: stat = 0
    logical :: isOpen = .false.
  contains
    procedure :: writeScalar => writeScalarFTN
    procedure :: write1DArray => write1DArrayFTN
    procedure :: write2DArray => write2DArrayFTN
    procedure :: write2DArrayDistBlacs => write2DArrayDistBlacsFTN
    procedure :: write2DArrayDistColumn => write2DArrayDistColumnFTN
    procedure :: readScalar => readScalarFTN
    procedure :: read1DArray => read1DArrayFTN
    procedure :: read2DArray => read2DArrayFTN
    procedure :: read2DArrayDistBlacs => read2DArrayDistBlacsFTN
    procedure :: open
    procedure :: close
    final :: destroyIoHandlerFTN
  end type ioHandlerFTN

  private

  public :: ioHandlerFTN

  contains

    subroutine destroyIoHandlerFTN(this)
      type(ioHandlerFTN) :: this
      call this%close()
    end subroutine

    subroutine open(this, fname, err, action, position, status, form, access)
      class(ioHandlerFTN) :: this
      type(ErrorType), intent(inout) :: err
      character (len = *), intent(in) :: fname
      character (len = *), intent(in) :: action
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

      print *, "FTN: Opening ", trim(fname), " with ", &
        trim(action), " ", &
        trim(positionVal), " ", &
        trim(statusVal), " ", &
        trim(formVal), " ", &
        trim(accessVal)

      open(newunit=this%iounit, action=action, &
        form=formVal, position=positionVal, status=statusVal, & 
        access=accessVal, file=fname, iostat=this%stat)

      if (this%stat == 0) then
        this%isOpen = .true.
      else
        RAISE_ERROR("ERROR: Could not open file", err)
      endif
    end subroutine

    subroutine close(this)
      class(ioHandlerFTN) :: this
      if (this%isOpen) then
        close(this%iounit)
        this%isOpen = .false.
      endif
    end subroutine

    subroutine writeScalarFTN(this, object)
      class(ioHandlerFTN) :: this
      class(*), intent(in) :: object

      select type(object)
      type is (integer)
        write(this%iounit) object
      type is (real)
        write(this%iounit) object
      type is (complex)
        write(this%iounit) object
      type is (character(len=*))
        write(this%iounit) trim(object)
      class default
        stop "ioHandlerFTN%writeScalarFTN: Tried to write unsupported type"
      end select
    end subroutine

    subroutine write1DArrayFTN(this, object)
      class(ioHandlerFTN) :: this
      class(*), dimension(:), intent(in) :: object

      select type(object)
      type is (integer)
        write(this%iounit) object
      type is (real)
        write(this%iounit) object
      type is (complex)
        write(this%iounit) object
      class default
        stop "ioHandlerFTN%write1DArrayFTN: Tried to write unsupported type"
      end select
    end subroutine

    subroutine write2DArrayFTN(this, object)
      class(ioHandlerFTN) :: this
      class(*), dimension(:,:), intent(in) :: object

      select type(object)
      type is (integer(int32))
        write(this%iounit) object
      type is (integer(int64))
        write(this%iounit) object
      type is (real(real32))
        write(this%iounit) object
      type is (real(real64))
        write(this%iounit) object
      type is (complex(kind=4))
        write(this%iounit) object
      type is (complex(kind=8))
        write(this%iounit) object
      class default
        stop "ioHandlerFTN%write2DArrayFTN: Tried to write unsupported type"
      end select
    end subroutine

    subroutine write2DArrayDistBlacsFTN(this, object, descr, block_type)
      ! Write arrays distributed using blacs

      class(ioHandlerFTN) :: this
      class(*), dimension(:,:), intent(in) :: object
      integer, intent(in) :: descr(9)
      type(MPI_Datatype), intent(in) :: block_type

      ! Using the fortran io_handler means array isn't distributed, just write normally
      call this%write2DArray(object)
    end subroutine

    subroutine write2DArrayDistColumnFTN(this, object, mdimen)
      ! Write arrays distributed as columns using co_distr_data

      class(ioHandlerFTN) :: this
      class(*), intent(in) :: object(:,:)
      integer, intent(in) :: mdimen ! Dimension of entire distributed array

      ! Using the fortran io_handler means array isn't distributed, just write normally
      call this%write2DArray(object)
    end subroutine

    subroutine readScalarFTN(this, object)
      class(ioHandlerFTN) :: this
      class(*), intent(out) :: object

      select type(object)
      type is (integer(int32))
        read(this%iounit) object
      type is (integer(int64))
        read(this%iounit) object
      type is (real(real32))
        read(this%iounit) object
      type is (real(real64))
        read(this%iounit) object
      type is (complex(kind=4))
        read(this%iounit) object
      type is (complex(kind=8))
        read(this%iounit) object
      type is (character(len=*))
        !print *, object
        !stop "DEBUG readScalarFTN"
        read(this%iounit) object
      class default
        stop "ioHandlerFTN%readScalarFTN: Tried to read unsupported type"
      end select
    end subroutine

    subroutine read1DArrayFTN(this, object)
      class(ioHandlerFTN) :: this
      class(*), dimension(:), intent(out) :: object

      select type(object)
      type is (integer(int32))
        read(this%iounit) object
      type is (integer(int64))
        read(this%iounit) object
      type is (real(real32))
        read(this%iounit) object
      type is (real(real64))
        read(this%iounit) object
      type is (complex(kind=4))
        read(this%iounit) object
      type is (complex(kind=8))
        read(this%iounit) object
      class default
        stop "ioHandlerFTN%read1DArrayFTN: Tried to read unsupported type"
      end select
    end subroutine

    subroutine read2DArrayFTN(this, object)
      class(ioHandlerFTN) :: this
      class(*), dimension(:,:), intent(out) :: object

      select type(object)
      type is (integer(int32))
        read(this%iounit) object
      type is (integer(int64))
        read(this%iounit) object
      type is (real(real32))
        read(this%iounit) object
      type is (real(real64))
        read(this%iounit) object
      type is (complex(kind=4))
        read(this%iounit) object
      type is (complex(kind=8))
        read(this%iounit) object
      class default
        stop "ioHandlerFTN%read2DArrayFTN: Tried to read unsupported type"
      end select
    end subroutine

    subroutine read2DArrayDistBlacsFTN(this, object, descr, block_type)
      ! read arrays distributed as columns using co_distr_data

      class(ioHandlerFTN) :: this
      class(*), dimension(:,:), intent(out) :: object
      integer, intent(in) :: descr(9) ! Description array outputted from co_block_type_init
      type(MPI_Datatype), intent(in) :: block_type

      ! Using the fortran io_handler means array isn't distributed, just read normally
      call this%read2DArray(object)
    end subroutine
end module
