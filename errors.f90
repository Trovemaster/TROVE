module Errors
  implicit none

  integer, parameter :: ERR_None = 0, &
                        ERR_Default = 1, &
                        ERR_FileNotFound = 2

  type:: ErrorType
    integer :: code
    character(len=256) :: message
  end type
end module Errors
