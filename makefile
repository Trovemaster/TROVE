################################################################################
## USER OPTIONS
################################################################################

PLAT ?=
EXE=j-trove$(PLAT).x
pot_user = pot_H2O_Conway

################################################################################
## COMPILER OPTIONS
################################################################################

COMPILER ?= intel
MODE ?= release

# === MPI

#FOR = mpif90
##FFLAGS =  -qopenmp -xHost -O3 -ip -g3 -traceback -DTROVE_USE_MPI_
#FFLAGS = -fopenmp -ffree-line-length-512 -march=native -O0 -fcray-pointer -g -fallow-argument-mismatch -fbacktrace -DTROVE_USE_MPI_ 

##LAPACK = -mkl=parallel -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64
#LAPACK = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_blacs_intelmpi_lp64 -lmkl_scalapack_lp64 -lmkl_core -lgomp -lpthread -lm -ldl

#ARPACK = -larpack

#LIB = $(LAPACK) $(ARPACK)

# === MPI

# Intel
#######
ifeq ($(strip $(COMPILER)),intel)
	FOR = ifort
	FFLAGS = -cpp -qopenmp -module $(OBJDIR)
	LAPACK = -mkl

	ifeq ($(strip $(MODE)),debug)
		FFLAGS += -O0 -g
	else
		FFLAGS += -O3
	endif

	# Alternative flags:
	#FFLAGS = -ip -align -ansi-alias -traceback -qopenmp -mcmodel=medium -parallel -nostandard-realloc-lhs -module $(OBJDIR)
	#LAPACK = -mkl=parallel

# gfortran
##########
else ifeq ($(strip $(COMPILER)),gfortran)
	FOR = gfortran
	FFLAGS = -cpp -std=gnu -fopenmp -march=native -ffree-line-length-512 -fcray-pointer -I$(OBJDIR) -J$(OBJDIR)

	GCC_VERSION_GT_10 := $(shell expr `gcc -dumpversion | cut -f1 -d.` \>= 10)
	ifeq "${GCC_VERSION_GT_10}" "1"
		# gcc 10 complains about mismatched argument types
		FFLAGS += -fallow-argument-mismatch
	endif

	ifeq ($(strip $(MODE)),debug)
		FFLAGS += -O0 -g -Wall -Wextra -fbacktrace
	else ifeq ($(strip $(MODE)),ci)
		FFLAGS += -O0 -g
	else
		FFLAGS += -O3
	endif

	# Use non-MKL LAPACK:
	#LAPACK = -llapack -lblas

	# Use MKL LAPACK:
	LAPACK += -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
else
$(error Compiler option "$(COMPILER)" not defined.)
endif

CPPFLAGS = -D_EXTFIELD_DEBUG_

################################################################################
## LIBRARIES
################################################################################

WIGXJPF_DIR = wigxjpf-1.5
WIGXJPF_LIB = $(WIGXJPF_DIR)/lib/libwigxjpf.a
LIB     =   $(LAPACK) $(LIBS) $(WIGXJPF_LIB)

################################################################################
## SOURCES & DIRECTORIES
################################################################################

BINDIR=.
SRCDIR=.
OBJDIR=.
user_pot_dir=.

SRCS := timer.f90 accuracy.f90 diag.f90 dipole.f90 extfield.f90 fields.f90 fwigxjpf.f90 input.f90 kin_xy2.f90 lapack.f90 \
	me_bnd.f90 me_numer.f90 me_rot.f90 me_str.f90 \
	mol_abcd.f90 mol_c2h4.f90 mol_c2h6.f90 mol_c3h6.f90 mol_ch3oh.f90 mol_xy.f90 \
	mol_xy2.f90 mol_xy3.f90 mol_xy4.f90 mol_zxy2.f90 mol_zxy3.f90 \
	molecules.f90 moltype.f90 mpi_aux.f90 perturbation.f90 plasma.f90 \
	pot_abcd.f90 pot_c2h4.f90 pot_c2h6.f90 pot_c3h6.f90 pot_ch3oh.f90 \
	pot_xy2.f90 pot_xy3.f90 pot_xy4.f90 pot_zxy2.f90 pot_zxy3.f90 \
	prop_xy2.f90 prop_xy2_quad.f90 prop_xy2_spinrot.f90 prop_xy2_spinspin.f90 \
	refinement.f90 richmol_data.f90 rotme_cart_tens.f90 symmetry.f90 tran.f90 trove.f90 $(pot_user).f90
OBJS := ${SRCS:.f90=.o}

VPATH = $(SRCDIR):$(user_pot_dir):$(OBJDIR)

################################################################################
## TARGETS
################################################################################

.PHONY: all, clean, cleanall, tarball, checkin, test

all: $(BINDIR) $(OBJDIR) $(BINDIR)/$(EXE)

%.o : %.f90
	$(FOR) -c $(FFLAGS) $(CPPFLAGS) -o $(OBJDIR)/$@ $<

$(BINDIR)/$(EXE): $(OBJS) $(WIGXJPF_LIB)
	$(FOR) $(FFLAGS) -o $@ $(addprefix $(OBJDIR)/,$(OBJS)) $(LIB)

$(WIGXJPF_LIB):
	$(MAKE) -C $(WIGXJPF_DIR)

$(OBJDIR):
ifneq ($(OBJDIR),.)
	mkdir -p $(OBJDIR)
endif

$(BINDIR):
ifneq ($(BINDIR),.)
	mkdir -p $(BINDIR)
endif

clean:
	rm -rf $(BINDIR)/$(EXE) $(OBJDIR)/*.mod $(OBJDIR)/*.o

cleanall: clean
	$(MAKE) -C $(WIGXJPF_DIR) clean

tarball:
	tar cf trove.tar makefile *.f90

checkin:
	ci -l Makefile *.f90

test: $(BINDIR)/$(EXE)
	cd test; ./run_regression_tests.sh

################################################################################
## DEPENDENCIES
################################################################################

pot_user_deps=$(shell grep -io '^\s*use [a-zA-Z0-9_]*' ${user_pot_dir}/${pot_user}.f90 | awk '{print $$2".o"}' | tr '\n' ' ')
$(pot_user).o: $(pot_user_deps)

accuracy.o: accuracy.f90 
diag.o: diag.f90 accuracy.o timer.o
dipole.o: dipole.f90 accuracy.o fields.o timer.o molecules.o moltype.o symmetry.o tran.o
extfield.o: extfield.f90 accuracy.o timer.o rotme_cart_tens.o richmol_data.o fields.o moltype.o tran.o symmetry.o
fields.o: fields.f90 accuracy.o molecules.o lapack.o me_str.o me_bnd.o me_numer.o me_rot.o timer.o moltype.o symmetry.o input.o accuracy.o moltype.o accuracy.o moltype.o accuracy.o moltype.o
fwigxjpf.o: fwigxjpf.f90 $(WIGXJPF_LIB)
grid.o: grid.f90 accuracy.o fields.o splines.o iso_c_binding.o iso_c_binding.o
input.o: input.f90 accuracy.o
kin_xy2.o: kin_xy2.f90 accuracy.o moltype.o
lapack.o: lapack.f90 accuracy.o timer.o
ltp2011_water_dipole_surface.o: ltp2011_water_dipole_surface.f90 
me_bnd.o: me_bnd.f90 accuracy.o me_numer.o molecules.o timer.o lapack.o moltype.o
me_numer.o: me_numer.f90 accuracy.o molecules.o moltype.o lapack.o timer.o
me_rot.o: me_rot.f90 accuracy.o timer.o
me_str.o: me_str.f90 accuracy.o timer.o me_numer.o
mol_abcd.o: mol_abcd.f90 accuracy.o moltype.o lapack.o pot_abcd.o symmetry.o
mol_c2h4.o: mol_c2h4.f90 accuracy.o moltype.o
mol_c2h6.o: mol_c2h6.f90 accuracy.o moltype.o
mol_c3h6.o: mol_c3h6.f90 accuracy.o moltype.o
mol_ch3oh.o: mol_ch3oh.f90 accuracy.o moltype.o lapack.o pot_ch3oh.o
mol_ch4.o: mol_ch4.f90 accuracy.o moltype.o lapack.o symmetry.o
molecules.o: molecules.f90 accuracy.o lapack.o moltype.o mol_xy.o mol_xy2.o mol_xy3.o mol_xy4.o mol_zxy2.o mol_zxy3.o mol_ch3oh.o mol_abcd.o mol_c2h4.o mol_c2h6.o mol_c3h6.o pot_xy2.o pot_xy3.o pot_abcd.o pot_zxy2.o pot_zxy3.o pot_xy4.o pot_ch3oh.o pot_c2h4.o pot_c2h6.o pot_c3h6.o prop_xy2.o prop_xy2_quad.o prop_xy2_spinrot.o prop_xy2_spinspin.o kin_xy2.o symmetry.o $(pot_user).o
moltype.o: moltype.f90 accuracy.o lapack.o
mol_user.o: mol_user.f90 accuracy.o moltype.o lapack.o symmetry.o
mol_xy2.o: mol_xy2.f90 accuracy.o moltype.o symmetry.o
mol_xy3.o: mol_xy3.f90 accuracy.o moltype.o lapack.o
mol_xy4.o: mol_xy4.f90 accuracy.o moltype.o lapack.o symmetry.o pot_xy4.o
mol_xy.o: mol_xy.f90 accuracy.o moltype.o
mol_zxy2.o: mol_zxy2.f90 accuracy.o moltype.o
mol_zxy3.o: mol_zxy3.f90 accuracy.o moltype.o lapack.o
mpi_aux.o: mpi_aux.f90 accuracy.o timer.o
perturbation.o: perturbation.f90 accuracy.o molecules.o moltype.o lapack.o plasma.o fields.o timer.o symmetry.o me_numer.o diag.o mpi_aux.o
plasma.o: plasma.f90 accuracy.o timer.o
pot_abcd.o: pot_abcd.f90 accuracy.o moltype.o lapack.o
pot_c2h4.o: pot_c2h4.f90 accuracy.o moltype.o
pot_c2h6.o: pot_c2h6.f90 accuracy.o moltype.o mol_c2h6.o
pot_c3h6.o: pot_c3h6.f90 accuracy.o moltype.o
pot_ch3oh.o: pot_ch3oh.f90 accuracy.o moltype.o lapack.o
pot_xy2.o: pot_xy2.f90 accuracy.o moltype.o
pot_xy3.o: pot_xy3.f90 accuracy.o moltype.o lapack.o
pot_xy4.o: pot_xy4.f90 accuracy.o moltype.o lapack.o symmetry.o
pot_zxy2.o: pot_zxy2.f90 accuracy.o moltype.o
pot_zxy3.o: pot_zxy3.f90 accuracy.o moltype.o
prop_xy2.o: prop_xy2.f90 accuracy.o moltype.o timer.o pot_xy2.o
prop_xy2_quad.o: prop_xy2_quad.f90 accuracy.o moltype.o timer.o pot_xy2.o
prop_xy2_spinrot.o: prop_xy2_spinrot.f90 accuracy.o moltype.o pot_xy2.o timer.o
prop_xy2_spinspin.o: prop_xy2_spinspin.f90 accuracy.o moltype.o pot_xy2.o timer.o
refinement.o: refinement.f90 accuracy.o fields.o timer.o molecules.o moltype.o symmetry.o lapack.o tran.o
richmol_data.o: richmol_data.f90 accuracy.o timer.o
rotme_cart_tens.o: rotme_cart_tens.f90 accuracy.o timer.o fwigxjpf.o moltype.o accuracy.o
symmetry.o: symmetry.f90 accuracy.o timer.o
timer.o: timer.f90 accuracy.o
tran.o: tran.f90 accuracy.o timer.o me_numer.o molecules.o fields.o moltype.o symmetry.o perturbation.o mpi_aux.o
trove.o: trove.f90 accuracy.o fields.o perturbation.o symmetry.o timer.o moltype.o dipole.o refinement.o tran.o extfield.o
