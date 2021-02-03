goal:   trove.x

tarball:
	tar cf trove.tar makefile *.f90
        
checkin:
	ci -l Makefile *.f90



pot_user = pot_H2O_Conway

PLAT = _0209
###FOR  = ifort
FOR = ifort 
#FFLAGS = -ip -O3 -align -ansi-alias -g -traceback  -qopenmp -mcmodel=medium -parallel -cpp -nostandard-realloc-lhs
FFLAGS = -O3 -qopenmp -cpp
CPPFLAGS = -D_EXTFIELD_DEBUG_


#ARPACK =  ~/libraries/ARPACK/libarpack_omp_64.a

#LAPACK = -mkl
LAPACK = -mkl=parallel -qopenmp
#LIBS   =  ./libplasma.a ./libcoreblas.a ./libquark.a ./libmrrr.a  -lpthread  -lnuma -lm
#LIBS   =   -lpthread  -lnuma -lm

LIB     =   $(LAPACK) $(LIBS) wigxjpf-1.5/lib/libwigxjpf.a


%.o : %.f90
	$(FOR) -c $(FFLAGS) $(CPPFLAGS) $<


###############################################################################

user_pot_dir=.
pot_user_deps=$(shell grep -io '^\s*use [a-zA-Z0-9_]*' ${user_pot_dir}/${pot_user}.f90 | awk '{print $$2".o"}' | tr '\n' ' ')
$(pot_user).o: $(pot_user_deps)

trove.x:        trove.o accuracy.o perturbation.o fields.o symmetry.o molecules.o me_numer.o me_str.o me_bnd.o me_rot.o \
	            lapack.o plasma.o moltype.o refinement.o dipole.o refinement.o tran.o diag.o timer.o input.o \
                mol_xy.o mol_xy2.o mol_xy3.o mol_xy4.o mol_zxy2.o mol_zxy3.o mol_ch3oh.o mol_abcd.o mol_c2h4.o mol_c2h6.o mol_c3h6.o \
				pot_xy2.o pot_xy3.o pot_xy4.o pot_zxy2.o pot_zxy3.o pot_ch3oh.o pot_abcd.o pot_c2h4.o pot_c2h6.o pot_c3h6.o  $(pot_user).o \
				prop_xy2.o prop_xy2_quad.o prop_xy2_spinrot.o prop_xy2_spinspin.o rotme_cart_tens.o fwigxjpf.o extfield.o \
				richmol_data.o kin_xy2.o 
				$(FOR) $(FFLAGS) -o j-trove$(PLAT).x $^ $(LIB)

accuracy.o: accuracy.f90 
diag.o: diag.f90 accuracy.o timer.o
dipole.o: dipole.f90 accuracy.o fields.o timer.o molecules.o moltype.o symmetry.o tran.o
extfield.o: extfield.f90 accuracy.o timer.o rotme_cart_tens.o richmol_data.o fields.o moltype.o tran.o symmetry.o
fields.o: fields.f90 accuracy.o molecules.o lapack.o me_str.o me_bnd.o me_numer.o me_rot.o timer.o moltype.o symmetry.o input.o accuracy.o moltype.o accuracy.o moltype.o accuracy.o moltype.o
#fwigxjpf.o: fwigxjpf.f90 
grid.o: grid.f90 accuracy.o fields.o splines.o iso_c_binding.o iso_c_binding.o
input.o: input.f90 
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
moltype.o: moltype.f90 accuracy.o lapack.o accuracy.o accuracy.o accuracy.o accuracy.o accuracy.o accuracy.o accuracy.o accuracy.o
mol_user.o: mol_user.f90 accuracy.o moltype.o lapack.o symmetry.o
mol_xy2.o: mol_xy2.f90 accuracy.o moltype.o symmetry.o
mol_xy3.o: mol_xy3.f90 accuracy.o moltype.o lapack.o
mol_xy4.o: mol_xy4.f90 accuracy.o moltype.o lapack.o symmetry.o pot_xy4.o
mol_xy.o: mol_xy.f90 accuracy.o moltype.o
mol_zxy2.o: mol_zxy2.f90 accuracy.o moltype.o
mol_zxy3.o: mol_zxy3.f90 accuracy.o moltype.o lapack.o
perturbation.o: perturbation.f90 accuracy.o molecules.o moltype.o lapack.o plasma.o fields.o timer.o symmetry.o me_numer.o diag.o
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
tran.o: tran.f90 accuracy.o timer.o me_numer.o molecules.o fields.o moltype.o symmetry.o perturbation.o
trove.o: trove.f90 accuracy.o fields.o perturbation.o symmetry.o timer.o moltype.o dipole.o refinement.o tran.o extfield.o

fwigxjpf.o:
	$(FOR) -c $(FFLAGS) -cpp fwigxjpf.f90
	cd wigxjpf-1.5  && $(MAKE)

clean:
	rm -f *.mod *.o
