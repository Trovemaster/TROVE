module io_factory
  use mpi_aux
  use io_handler_base
  use io_handler_ftn
#ifdef TROVE_USE_MPI_
  use io_handler_mpi
#endif
  use errors
  implicit none

  contains
    subroutine allocateHandler(ioHandler)
      class(ioHandlerBase), allocatable, intent(out) :: ioHandler
#ifdef TROVE_USE_MPI_
      if (blacs_size > 1) then
        ! Only allocate MPI IO when compiled with MPI *and* 
        ! running with more than one MPI process
        allocate(ioHandlerMPI::ioHandler)
      else
        allocate(ioHandlerFTN::ioHandler)
      endif
#else
      allocate(ioHandlerFTN::ioHandler)
#endif
    end subroutine

    subroutine openFile(ioHandler, filename, err, action, position, status, form, access)
      class(ioHandlerBase), allocatable, intent(out) :: ioHandler
      character (len = *), intent(in) :: filename
      type(ErrorType), intent(inout) :: err
      character (len = *), intent(in) :: action
      character (len = *), intent(in), optional :: position, status, form, access

      call allocateHandler(ioHandler)
      call ioHandler%open(filename, err, action, position, status, form, access)
    end subroutine
end module
