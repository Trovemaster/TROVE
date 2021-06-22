include ../../lib/pFUnit/build/PFUNIT.mk

all: test_simple

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I../..

test_simple_TESTS := test_simple.pf test_simple_fixture.pf
$(eval $(call make_pfunit_test,test_simple))

clean:
	$(RM) *.o *.mod *.a  *.inc
	$(RM) test_simple.F90
