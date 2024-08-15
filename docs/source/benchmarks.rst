Benchmarks and examples
***********************

A set of TROVE benchmarks for a number of typical systems can be found here

`Benchmarks <https://github.com/Trovemaster/TROVE-benchmarks>`__.


PH\ :sub:`3`
============


This PH\ :sub:`3` benchmark is to produce a small line list for phosphine (PH\ :sub:`3`) using a relatively  large primitive bais set.

It includes the following set of jobs:
::
    
   ./trove.x <ph3_P18_step1.inp > ph3_P18_step1.out
   ./trove.x <ph3_P18_step2.inp > ph3_P18_step3.out
   ......
   ./trove.x <ph3_P18_step10.inp > ph3_P18_step10.out
   ./trove.x <ph3_P18_intensity.inp > ph3_P18_intensity.out
    

covering the following steps:

- Step 1: Produce the initial smaller :math:`J=0` HAMILTONIAN matrix and diagonalize it using ``DSYEV`` (step1)
- Step 2: Convert to :math:`J=0` representation; prepare matrix elements for  :math:`J>0` calculations (step2)
- Step 3: Generate a :math:`J=0-5` 10 energies and eigenfunctions  (step3-step9) by diagonalising three double real, symmetric matrices for each J using  DSYEV.
         Each calculation is independent, but depends on the checkpoints produced at steps 1-2. It requires checkpoints from steps 1-2.
- Step 4: Generate a :math:``J=100` (double real, symmetric) matrix  and store to the disk without diagonalisation (step 10)
- Step 5: Intensity calculations between J=0-5 states (ph3_P18_step_intensity.inp).
        It will depend on the checkpoints produced at steps 1-3.

Memory and time requirements:
 
 -Step1:    84 Gb 27032.7
 -Step2:   9.9 Gb 1220.9 s
 -Step3:   0.8 Gb 2.6 s
 -Step4:   6.8 Gb 73 s
 -Step5:   7.8 Gb 198 s
 -Step6:   8.4 Gb  502 s
 -Step7:   9.0 Gb  830 s
 -Step8:   9.6 Gb  1449 s
 -Step9:  19.3 Gb  9032 s
 -Step10: 32.5 Gb  12096 s
 -Step_intensity: 89 Gb  1592 s

Disk space used: 300 Gb

Here :math:`J` is the rotational angular momentum. 

There are three  matrices to diagonalise for each symmetry :math:`\Gamma=A_1, A_2, E`. The dimensions of the matrices are  given approximately  by

.. math:: 
  
  N(A) = (2J+1) x N_0(A) 
   N(E) = (2J+1) x N_0(E)
   
      
where :math:`N_0` are the  dimensions of the :math:`\Gamma=A_1, A_2` and  :math:`E` matrices at :math:`J=0`, :math:`N_0=1000` and 2000, respectively  (approximately). 

The first step in this case is relatively expensive, 84 Gb and 27032.7 s on a 16 processor workstation.



