BASE_DIR := ../..

include $(BASE_DIR)/lib/pFUnit/build/PFUNIT.mk

all: test_io

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I../..

test_io_TESTS := test_io.pf
test_io_OTHER_SRCS := $(addprefix $(BASE_DIR)/, writer_ftn.f90 writer_base.f90)
$(eval $(call make_pfunit_test,test_io))

clean:
	$(RM) *.o *.mod *.a *.inc
	$(RM) test_io.F90 test_io test.dat
