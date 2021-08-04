BASE_DIR := ../..

include $(BASE_DIR)/lib/pFUnit/build/PFUNIT.mk

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I../..

test_io_TESTS := test_io.pf
test_io_OTHER_SRCS := $(addprefix $(BASE_DIR)/, io_handler_ftn.f90 io_handler_base.f90 )
$(eval $(call make_pfunit_test,test_io))

ifdef USE_MPI
FFLAGS += -DTROVE_USE_MPI_
test_mpi_io_TESTS := test_mpi_io.pf
test_mpi_io_OTHER_LIBRARIES := ${LAPACK} # passed in from base makefile
test_mpi_io_OTHER_SRCS := $(addprefix $(BASE_DIR)/, io_handler_mpi.f90 io_handler_base.f90 mpi_aux.f90 timer.f90 accuracy.f90)
$(eval $(call make_pfunit_test,test_mpi_io))
endif

clean:
	$(RM) *.o *.mod *.a *.inc
	$(RM) test_io.F90 test_io
	$(RM) test_mpi_io.F90 test_mpi_io
