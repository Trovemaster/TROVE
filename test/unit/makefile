BASE_DIR := ../..

include $(BASE_DIR)/lib/pFUnit/build/PFUNIT.mk

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I../..

LAPACK = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
ifdef USE_MPI
LAPACK += -lmkl_blacs_openmpi_lp64 -lmkl_scalapack_lp64
endif

test_io_TESTS := test_io.pf
test_io_OTHER_SRCS := $(addprefix $(BASE_DIR)/, writer_ftn.f90 writer_base.f90 reader_ftn.f90 reader_base.f90)
$(eval $(call make_pfunit_test,test_io))

ifdef USE_MPI
FFLAGS += -DTROVE_USE_MPI_
test_mpi_io_TESTS := test_mpi_io.pf
test_mpi_io_OTHER_LIBRARIES := ${LAPACK}
test_mpi_io_OTHER_SRCS := $(addprefix $(BASE_DIR)/, writer_mpi.f90 writer_base.f90 reader_mpi.f90 reader_base.f90 mpi_aux.f90 timer.f90 accuracy.f90)
$(eval $(call make_pfunit_test,test_mpi_io))
endif

clean:
	$(RM) *.o *.mod *.a *.inc
	$(RM) test_io.F90 test_io test.dat
	$(RM) test_mpi_io.F90 test_mpi_io test.dat
