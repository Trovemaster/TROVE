goal:   trove.x

tarball:
	tar cf trove.tar makefile *.f90
        
checkin:
	ci -l Makefile *.f90



pot_user = pot_H2O_Conway

PLAT = _0209
###FOR  = ifort
FOR = ifort 
FFLAGS = -ip   -qopenmp -O3  -static -cpp


#ARPACK =  ~/libraries/ARPACK/libarpack_omp_64.a

#LAPACK = -mkl
LAPACK = -mkl=parallel


LIB     =   $(LAPACK) 

EXTERNAL = H2O_POTV_2018.o

%.o : %.f90
	$(FOR) -c $(FFLAGS) $<


###############################################################################

trove.x:        trove.o accuracy.o perturbation.o fields.o symmetry.o molecules.o me_numer.o me_str.o me_bnd.o me_rot.o \
	                 lapack.o plasma.o moltype.o refinement.o dipole.o refinement.o tran.o diag.o timer.o input.o \
                   mol_xy.o mol_xy2.o mol_xy3.o mol_xy4.o mol_zxy2.o mol_zxy3.o mol_ch3oh.o mol_abcd.o mol_c2h4.o mol_c2h6.o mol_c3h6.o \
									          pot_xy2.o pot_xy3.o pot_xy4.o pot_zxy2.o pot_zxy3.o \
                   pot_ch3oh.o pot_abcd.o pot_c2h4.o pot_c3h6.o pot_c2h6.o $(pot_user).o prop_xy2.o kin_xy2.o  # $(EXTERNAL)

	$(FOR) $(FFLAGS) -o j-trove$(PLAT).x $^ $(LIB)

trove.o:        accuracy.o fields.o perturbation.o symmetry.o timer.o moltype.o dipole.o refinement.o tran.o
perturbation.o: accuracy.o molecules.o lapack.o fields.o timer.o symmetry.o diag.o plasma.o
fields.o:       accuracy.o molecules.o lapack.o me_str.o timer.o me_numer.o input.o me_rot.o moltype.o symmetry.o me_bnd.o
symmetry.o:     accuracy.o
molecules.o:    accuracy.o moltype.o mol_xy.o mol_xy2.o mol_xy3.o mol_xy4.o mol_zxy2.o mol_zxy3.o mol_ch3oh.o mol_abcd.o mol_c2h4.o mol_c2h6.o mol_c3h6.o \
													 lapack.o	          pot_xy2.o pot_xy3.o mol_xy4.o pot_zxy2.o pot_zxy3.o pot_ch3oh.o pot_abcd.o pot_c2h4.o pot_c2h6.o pot_c3h6.o \
													 symmetry.o mol_c3h6.o pot_c3h6.o  $(pot_user).o prop_xy2.o kin_xy2.o

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


clean:
	rm -f *.mod *.o
