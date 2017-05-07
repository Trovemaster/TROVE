goal:   trove.x

tarball:
	tar cf trove.tar makefile *.f90
        
checkin:
	ci -l Makefile *.f90


FPATH = 

pot_user = pot_ch4

PLAT = 0105_i16_II
###FOR  = ifort
FOR = ifort 
# -multiple-processes=4
###FOR = /home/syurchenko/opt/intel/intel_cluster/composerxe-2011.5.220/composerxe-2011.5.220/bin/intel64/ifort

#FFLAGS = -C -static
# FFLAGS = -O3 -ip  
FFLAGS = -ip   -openmp -O3 -static

#-guide-vec 
###FFLAGS = -O0 -fpe0  -fltconsistency -stand f03 -check all -warn all -traceback -fp-stack-check  # debugging options


#-guide-vec -parallel 
# -openmp-report2 
#FFLAGS = -C -check bounds -g  -gen-interfaces -warn interfaces  -check arg_temp_created -prof-value-profiling=all -warn all


#LAPACK =   -L  /usr/local/share/intel/mkl11.1/lib/em64t/ -lmkl_solver  -lmkl_lapack -lmkl_em64t -lmkl_intel_thread  -lguide -lpthread

#LDIR = /opt/intel/composerxe-2011.4.191/mkl/lib/intel64/
#LAPACK =  -Wl,--start-group $(LDIR)libmkl_intel_lp64.a   $(LDIR)libmkl_intel_thread.a $(LDIR)libmkl_core.a   -Wl,--end-group  -lpthread

ARPACK =  ~/libraries/ARPACK/libarpack_omp_64.a

#LAPACK = -mkl
LAPACK = -mkl=parallel

PROPACK = ~/libraries/PROPACK/libdpropack_LINUX_ICC_X86_64.a

#LDIR = /usr/local/intel/Compiler/11.1/059/mkl/lib/32/
#LAPACK =  -Wl,--start-group $(LDIR)libmkl_intel.a   \
# $(LDIR)libmkl_intel_thread.a $(LDIR)libmkl_core.a   -Wl,--end-group -lguide -lpthread
#ARPACK =  ~/libraries/ARPACK/libarpack_omp.a


LIB     =  $(ARPACK) $(LAPACK) 



###############################################################################

OBJ = perturbation.o molecules.o me_str.o symmetry.o \
      fields.o accuracy.o lapack.o timer.o me_numer.o  \
      input.o me_rot.o moltype.o dipole.o tran.o refinement.o \
      mol_xy.o mol_xy2.o mol_xy3.o mol_xy4.o mol_zxy2.o mol_abcd.o  mol_ch3oh.o mol_c2h4.o diag.o me_bnd.o mol_zxy3.o \
               pot_xy2.o pot_xy3.o pot_xy4.o pot_zxy2.o pot_abcd.o  pot_ch3oh.o pot_c2h4.o  $(pot_user).o plasma.o pot_zxy3.o

trove.x:       $(OBJ) trove.o
	$(FOR) -o j-trove$(PLAT).x $(OBJ) $(FFLAGS) trove.o $(LIB)

trove.o:       trove.f90 $(OBJ) 
	$(FOR) -c trove.f90 $(FFLAGS)

perturbation.o:  perturbation.f90 accuracy.o molecules.o lapack.o fields.o timer.o symmetry.o diag.o plasma.o
	$(FOR) -c perturbation.f90 -O3 -ip  $(FFLAGS)  

#-inline-level=2 -opt-matmul

fields.o:  fields.f90 accuracy.o molecules.o lapack.o me_str.o timer.o me_numer.o input.o me_rot.o moltype.o symmetry.o me_bnd.o
	$(FOR) -c fields.f90 $(FFLAGS)

symmetry.o:  symmetry.f90 accuracy.o
	$(FOR) -c symmetry.f90 $(FFLAGS)

me_numer.o:  me_numer.f90 accuracy.o molecules.o timer.o
	$(FOR) -c me_numer.f90 $(FFLAGS)

molecules.o:  molecules.f90 accuracy.o lapack.o mol_xy.o mol_xy2.o mol_xy3.o mol_xy4.o mol_zxy2.o mol_c2h4.o mol_abcd.o mol_ch3oh.o mol_zxy3.o \
              pot_xy2.o pot_xy3.o pot_xy4.o pot_zxy2.o pot_abcd.o pot_ch3oh.o pot_c2h4.o symmetry.o pot_zxy3.o   $(pot_user).o
	$(FOR) -c molecules.f90 $(FFLAGS)

me_str.o:  me_str.f90 accuracy.o
	$(FOR) -c me_str.f90 $(FFLAGS)

me_bnd.o:  me_bnd.f90 accuracy.o timer.o me_numer.o
	$(FOR) -c me_bnd.f90 $(FFLAGS)

me_rot.o:  me_rot.f90 accuracy.o
	$(FOR) -c me_rot.f90 $(FFLAGS)

accuracy.o:  accuracy.f90
	$(FOR) -c accuracy.f90 $(FFLAGS)

lapack.o:  lapack.f90 accuracy.o timer.o 
	$(FOR) -c lapack.f90 $(FFLAGS)

plasma.o:  plasma.f90 accuracy.o timer.o 
	$(FOR) -c plasma.f90 $(FFLAGS)

timer.o:  timer.f90
	$(FOR) -c timer.f90 $(FFLAGS)

mol_xy.o:  mol_xy.f90 moltype.o
	$(FOR) -c mol_xy.f90 $(FFLAGS)

mol_xy2.o:  mol_xy2.f90 moltype.o symmetry.o
	$(FOR) -c mol_xy2.f90 $(FFLAGS)

mol_c2h4.o:  mol_c2h4.f90 moltype.o
	$(FOR) -c mol_c2h4.f90 $(FFLAGS)

mol_xy3.o:  mol_xy3.f90 moltype.o
	$(FOR) -c mol_xy3.f90 $(FFLAGS)

mol_xy4.o:  mol_xy4.f90 moltype.o symmetry.o pot_xy4.o
	$(FOR) -c mol_xy4.f90 $(FFLAGS)

mol_zxy2.o:  mol_zxy2.f90 moltype.o accuracy.o
	$(FOR) -c mol_zxy2.f90 $(FFLAGS)

mol_zxy3.o:  mol_zxy3.f90 moltype.o accuracy.o
	$(FOR) -c mol_zxy3.f90 $(FFLAGS)

mol_ch3oh.o:	mol_ch3oh.f90 moltype.o pot_ch3oh.o
	$(FOR) -c mol_ch3oh.f90 $(FFLAGS)
        
mol_abcd.o: moltype.o mol_abcd.f90 pot_abcd.o
	$(FOR) -c mol_abcd.f90 $(FFLAGS)

pot_xy2.o:  pot_xy2.f90 moltype.o 
	$(FOR) -c pot_xy2.f90 $(FFLAGS)

pot_xy3.o:  pot_xy3.f90 moltype.o
	$(FOR) -c pot_xy3.f90 $(FFLAGS)

pot_xy4.o:  pot_xy4.f90 moltype.o symmetry.o
	$(FOR) -c pot_xy4.f90 $(FFLAGS)

pot_zxy2.o:  pot_zxy2.f90 moltype.o accuracy.o
	$(FOR) -c pot_zxy2.f90 $(FFLAGS)

pot_zxy3.o:  pot_zxy3.f90 moltype.o accuracy.o
	$(FOR) -c pot_zxy3.f90 $(FFLAGS)

pot_ch3oh.o:	pot_ch3oh.f90 moltype.o
	$(FOR) -c pot_ch3oh.f90 $(FFLAGS)

pot_c2h4.o:	pot_c2h4.f90 moltype.o
	$(FOR) -c pot_c2h4.f90 $(FFLAGS)

        
pot_abcd.o: pot_abcd.f90 
	$(FOR) -c pot_abcd.f90 $(FFLAGS)

moltype.o:  moltype.f90  accuracy.o lapack.o
	$(FOR) -c moltype.f90 $(FFLAGS)

dipole.o:  dipole.f90 accuracy.o fields.o molecules.o timer.o tran.o
	$(FOR) -c dipole.f90 $(FFLAGS)

refinement.o:  refinement.f90 accuracy.o fields.o molecules.o timer.o tran.o
	$(FOR) -c refinement.f90 $(FFLAGS)

tran.o:  tran.f90 accuracy.o fields.o molecules.o timer.o
	$(FOR) -c tran.f90 $(FFLAGS)

#pot_ch4.o:  pot_ch4.f90
#	$(FOR) -c pot_ch4.f90 $(FFLAGS)

$(pot_user).o:  $(pot_user).f90
	$(FOR) -c $(pot_user).f90 $(FFLAGS)

input.o:  input.f90
	$(FOR) -c input.f90 $(FFLAGS)

diag.o:  diag.f90
	$(FOR) -c diag.f90 $(FFLAGS)


clean:
	rm $(OBJ) *.mod trove.o

