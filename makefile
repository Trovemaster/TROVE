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

trove.x:        trove.o accuracy.o perturbation.o fields.o symmetry.o molecules.o me_numer.o me_str.o me_bnd.o me_rot.o \
	            lapack.o plasma.o moltype.o refinement.o dipole.o refinement.o tran.o diag.o timer.o input.o \
                mol_xy.o mol_xy2.o mol_xy3.o mol_xy4.o mol_zxy2.o mol_zxy3.o mol_ch3oh.o mol_abcd.o mol_c2h4.o mol_c2h6.o mol_c3h6.o \
				pot_xy2.o pot_xy3.o pot_xy4.o pot_zxy2.o pot_zxy3.o pot_ch3oh.o pot_abcd.o pot_c2h4.o pot_c2h6.o pot_c3h6.o  $(pot_user).o \
				prop_xy2.o prop_xy2_quad.o prop_xy2_spinrot.o prop_xy2_spinspin.o rotme_cart_tens.o fwigxjpf.o extfield.o \
				richmol_data.o kin_xy2.o 
				$(FOR) $(FFLAGS) -o j-trove$(PLAT).x $^ $(LIB)

trove.o:        accuracy.o fields.o perturbation.o symmetry.o timer.o moltype.o dipole.o refinement.o tran.o extfield.o
perturbation.o: accuracy.o molecules.o lapack.o fields.o timer.o symmetry.o diag.o plasma.o
fields.o:       accuracy.o molecules.o lapack.o me_str.o timer.o me_numer.o input.o me_rot.o moltype.o symmetry.o me_bnd.o
symmetry.o:     accuracy.o
molecules.o:    accuracy.o moltype.o mol_xy.o mol_xy2.o mol_xy3.o mol_xy4.o mol_zxy2.o mol_zxy3.o mol_ch3oh.o mol_abcd.o mol_c2h4.o mol_c3h6.o  mol_c2h6.o \
				lapack.o pot_xy2.o pot_xy3.o mol_xy4.o pot_zxy2.o pot_zxy3.o pot_ch3oh.o pot_abcd.o pot_c2h4.o  pot_c3h6.o  pot_c2h6.o \
				symmetry.o $(pot_user).o prop_xy2.o prop_xy2_quad.o prop_xy2_spinrot.o prop_xy2_spinspin.o kin_xy2.o

me_numer.o:     accuracy.o molecules.o timer.o
me_str.o:       accuracy.o timer.o me_numer.o
me_bnd.o:       accuracy.o timer.o me_numer.o
me_rot.o:       accuracy.o timer.o

lapack.o:       accuracy.o timer.o
plasma.o:       accuracy.o timer.o
moltype.o:      accuracy.o lapack.o
dipole.o:       accuracy.o fields.o molecules.o timer.o tran.o
refinement.o:   accuracy.o fields.o molecules.o timer.o tran.o
tran.o:         accuracy.o fields.o molecules.o timer.o me_numer.o fields.o moltype.o symmetry.o perturbation.o
diag.o:         accuracy.o timer.o
timer.o:        accuracy.o

mol_xy.o:       accuracy.o moltype.o
mol_xy2.o:      accuracy.o moltype.o symmetry.o
mol_xy3.o:      accuracy.o moltype.o
mol_xy4.o:      accuracy.o moltype.o symmetry.o pot_xy4.o
mol_zxy2.o:     accuracy.o moltype.o
mol_zxy3.o:     accuracy.o moltype.o
mol_ch3oh.o:    accuracy.o moltype.o pot_ch3oh.o
mol_c2h4.o:	    accuracy.o moltype.o
mol_c3h6.o:	    accuracy.o moltype.o
mol_c2h6.o:     accuracy.o moltype.o
mol_abcd.o:     accuracy.o moltype.o pot_abcd.o

pot_ch4.o:      accuracy.o moltype.o
pot_xy2.o:      accuracy.o moltype.o
pot_xy3.o:      accuracy.o moltype.o
pot_xy4.o:      accuracy.o moltype.o symmetry.o
pot_zxy2.o:     accuracy.o moltype.o
pot_zxy3.o:     accuracy.o moltype.o
pot_c2h6.o:     accuracy.o moltype.o
pot_ch3oh.o:	  accuracy.o moltype.o
pot_c2h4.o:	    accuracy.o moltype.o
pot_c3h6.o:	    accuracy.o moltype.o
pot_c2h6.o:     accuracy.o moltype.o
pot_abcd.o:     accuracy.o moltype.o lapack.o
prop_xy2.o:     accuracy.o moltype.o

kin_xy2.o:      accuracy.o moltype.o


#H2O_POTV_2018.o:    


prop_xy2.o:     accuracy.o moltype.o timer.o pot_xy2.o
prop_xy2_quad.o:     accuracy.o moltype.o timer.o pot_xy2.o
prop_xy2_spinrot.o:  accuracy.o moltype.o timer.o pot_xy2.o
prop_xy2_spinspin.o:  accuracy.o moltype.o timer.o pot_xy2.o

extfield.o:         accuracy.o moltype.o fields.o symmetry.o tran.o timer.o rotme_cart_tens.o tran.o richmol_data.o 
rotme_cart_tens.o:  accuracy.o timer.o fwigxjpf.o
richmol_data.o:     accuracy.o timer.o
fwigxjpf.o:
	$(FOR) -c $(FFLAGS) -cpp fwigxjpf.f90
	cd wigxjpf-1.5  && $(MAKE)


clean:
	rm -f *.mod *.o
