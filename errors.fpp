! Error.fpp

! From http://www.luckingtechnotes.com/fortran-error-handling-techniques on 23/06/2021

! Macros for error handling. 
! Enables user to store errors and exit the subroutine in single statement. 
! Fortran preprocessor must be enabled: -fpp.

! Raise Error
! Store the error code and info (only if the current code is zero).
! Return from the subroutine.
#define RAISE_ERROR(msg, err) if (err%Code == ERR_None) then; err = ErrorType(Code=ERR_Default, Message=msg); end if; return;
 
! Pass Error
! Returns if there's an error.
#define HANDLE_ERROR(err) if (err%Code /= ERR_None) then; print *, err%message; return; end if;
