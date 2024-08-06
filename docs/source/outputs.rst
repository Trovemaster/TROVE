Outputs
*******

.. _outputs:

In this chapter important and useful parts of the TROVE main output and checkpoint files will be described. The TROVE output file is rather large, especially if the value of ``verbose`` is high. Here only the most important parts of the output will be discussed in detail.

The TROVE output will differ depending on the checkpoints read/saved. The order given below is from setting up a calculation
from the start.

Beginning Output
================

The first part of the TROVE output file is a repeat of the input followed by statements about setting up the calculation. Details on the geometry of the molecule are given such as the number of angles, connectivity and Cartesian coordinates. When setting up a new molecule (see Chapter `New molecule <https://spectrove.readthedocs.io/en/latest/newmolecules.html>`__) this section is important as errors are often highlighted. Details of this will be given in that chapter.

Details of the S-matrix calculation are also printed in the section (see Theory chapter).

Trove will then print all of the finite difference steps which will be carried out to calculate the potential and then print out the expansion parameters of the potential and pseudo-potential.


Numerov basis
-------------

The details of the Numerov-Cooley generated basis can be reached by searching the output file for ``Numerov``. This section contains a grid for each coordinate of the form
::

     Numerov matrix elements calculations
     vmax = 10
     maxorder = 8
     icoord = 1
     rho_b (x) = -0.5000  0.6000
     rhostep (x) = 0.0003
     grid values (i,rho,rho_kinet,rho_poten,poten, mu_rr, f):
            0  -0.500000  -0.500000  -1.561761  299084.  5.61921  -0.920017E-25
            2  -0.499450  -0.499450  -1.559112  296470.  5.61921  -0.203945E-25

`vmax` is the maximum vibrational quantum number to be generated for the one-dimensional basis, 10 in this example.

``Maxorder`` is the expansion order of the potential as set by ``PotOrder``.
``icoord`` is which coordinate the basis is being generated for. TROVE will only generate
one set of Numerov basis functions if coordinates have been grouped together in the Basis block and have the same
ranges.
``rho_b (x)`` and ``rhostep (x)`` specify the range and step size for the grid as specified in the Basis block.
The grid is then given explicitly. A plot of ``rho`` against ``poten`` will show the potential expanded by TROVE.
mu_rr is the reduced mass for this coordinate. ``f`` gives the Numerov basis function for the ground state,
this is only printed with a sufficiently high value of ``verbose``.

At the end of the grid the outer most points where the one-dimensional vibrational wavefunctions have a minimum set value is given. The energies of each basis function are given (adjusted so that :math:`\nu = 0` has zero energy) along with the absolute zero point energy. A check is then carried out to see if the basis functions are orthonormal (to within numerical tolerance).

Contracted basis
----------------

After the primitive basis functions haven been generated, often using the Numerov method, TROVE then builds contractions of these functions. This procedure has been discussed in detail in [17YuYaOv]_ and here in Chapter `Theory <https://spectrove.readthedocs.io/en/latest/theory.html>`__. TROVE diagonalises a reduced Hamiltonian and the energies and primitive functions are given in a list, for example for PF\ :sub:`3`
::

      Variational eigenvalues:
            i        value       quanta
            1        0.00000000  0  0  0  0  0  0  0
            2      863.60110491  0  1  0  0  0  0  0
            3      863.60110491  0  0  0  1  0  0  0
            4      885.56391158  0  0  0  1  0  0  0
            .         .


The symmetry of these eigenfunctions of the reduced Hamiltonian are then reported
::

      Symmetry of the contracted solution, class:   1
            i       ener         deg  symmetry  quanta:
            1        0.00000000   1   1  A1     0   0   0   0   0   0    0
            2      863.60110491   2   3  E      1   0   0   0   0   0    1
            3      885.56391158   1   1  A1     0   0   1   0   0   0    0
            4     1718.03977668   1   1  A1     2   0   0   0   0   0    0


TROVE will then make use of the symmetry of these functions to set up matrix elements of the full Hamiltonian.


TROVE energy output
-------------------

The final step of a vibrational or rotational calculation is an output of the rotational-vibrational energies. These are ordered by energy and separated into symmetry blocks.


Consider an example of H\ :sub:`2`\ S calculations. The the J=0 energy output for the total symmetry :math:`A_1` is given by
::

      ------ ------ ------------- ------- --- -- --- ---- ----- ---- --- ------ -------- ------ --- --- --- -------- -----
        1        2        3           4    5  6  7    8    9     10   11  12       13      14  15  16  17       18    19
      ------ ------ ------------- ------- --- -- --- ---- ----- ---- --- ------ -------- ------ --- --- --- -------- -----
      G_tot      n     Energy        G_r   J  K  t    Gv1 Gv2     v1  v2  v3       C_i      n1  n2  n3  l    ivib1   ivib2
      ------ ------ ------------- ------- --- -- --- ---- ----- ---- --- ------ -------- ------ --- --- --- -------- -----
       A1        1      0.000000   ( A1 ;  0  0  0 ) ( A1  A1 ;   0   0   0 )      1.00 (   0   0   0   0 ) (    1    1 )
       A1        2   1172.667646   ( A1 ;  0  0  0 ) ( A1  A1 ;   0   0   1 )      1.00 (   0   0   1   0 ) (    1    2 )
       A1        3   2335.297519   ( A1 ;  0  0  0 ) ( A1  A1 ;   0   0   2 )      0.98 (   0   0   2   0 ) (    1    3 )
       A1        4   2608.713940   ( A1 ;  0  0  0 ) ( A1  A1 ;   1   0   0 )      0.99 (   1   0   0   0 ) (    2    1 )
       A1        5   3503.042415   ( A1 ;  0  0  0 ) ( A1  A1 ;   0   0   3 )      0.95 (   0   0   3   0 ) (    1    4 )
       A1        6   3765.459944   ( A1 ;  0  0  0 ) ( A1  A1 ;   1   0   1 )      0.95 (   1   0   1   0 ) (    2    2 )
       A1        7   4675.006191   ( A1 ;  0  0  0 ) ( A1  A1 ;   0   0   4 )      0.92 (   0   0   4   0 ) (    1    5 )
       ...........
       ...........



where the designation of the columns is as follows

  - Col 1: ``G_tot`` is the total symmetry of a ro-vibrational state;
  - Col 2: ``n`` is the counting number of the energy;
  - Col 3: ``Energy`` term value of the state;
  - Col 4: ``G_r`` is the rotational symmetry;
  - Col 5: ``J`` is the total angular momentum rotational quantum number;
  - Col 6: ``k`` is a rotational quantum number (projection of :math:`J` on the molecular axis :math:`z` );
  - Col 7: ``t`` is a rotational index defining the state parity :math:`\tau`;
  - Col 8-9: ``Gv1`` are ``Gv2`` are the vibrational symmetries of the corresponding vibrational sub-classes;
  - Cols 10-12: ``v1``, ``v2``, ``v3`` are the TROVE (local mode) vibrational quantum numbers;
  - Col 13: ``C_i`` is the largest eigen-coefficient used in the assignment.
  - Cols 14-17: ``K, n1, n2, n3`` are placeholder for the user-defined quantum numbers to be propagated to the final ro-vibrational eigenstates.
  - Cols 18-19: ``ivib1``, ``ivib2`` are the counting indices of sub-classes in the representation of direct products of the symmetry adapted 'contracted' basis set.


It should be noted that for equivalent modes, such as the two stretches in the case of H\ :sub:`2`\ S, only their total quanta :math:`v_1+v_2` is meaningful, not the individual values. For example, the following TROVE stretching states  :math:`(v_1,v_2) = (2,0), (1,1)`\ , and :math:`(0,2)` (:math:`v_1+v+2 = 2`\ ) are equivalent and cannot be distinguished without some extra information (e.g. their symmetry).


Similarity, the :math:`B_2` symmetry TROVE output is given by
::

      Variational solution - irreducible representation
        Gamma     i       value             j  k  t   quanta
        B2        1   3280.145078   ( A1 ;  0  0  0 ) ( B2  A1 ;   0   1   0 )      1.00 (   1   0   0   0 ) (    3    1 )
        B2        2   4415.876421   ( A1 ;  0  0  0 ) ( B2  A1 ;   0   1   1 )      0.99 (   1   0   1   0 ) (    3    2 )
        B2        3   5556.806722   ( A1 ;  0  0  0 ) ( B2  A1 ;   0   1   2 )      0.97 (   1   0   2   0 ) (    3    3 )
        B2        4   5785.428853   ( A1 ;  0  0  0 ) ( B2  A1 ;   0   2   0 )      0.99 (   2   0   0   0 ) (    5    1 )
        B2        5   6717.570020   ( A1 ;  0  0  0 ) ( B2  A1 ;   0   1   3 )      0.96 (   1   0   3   0 ) (    3    4 )
        B2        6   6914.548146   ( A1 ;  0  0  0 ) ( B2  A1 ;   0   2   1 )      0.96 (   2   0   1   0 ) (    5    2 )
        B2        7   8041.707663   ( A1 ;  0  0  0 ) ( B2  A1 ;   0   2   2 )      0.98 (   2   0   2   0 ) (    5    3 )
        ....


The non-rigourous quantum numbers :math:`K` and :math:`v_i` are defined using the largest eigen-coefficient  approach and are approximate. They represent the measure of how the given wavefunction is similar to a single selected basis set function selected as the largest contribution the corresponding expansion.  The quality of the assignment can be judged based on the expansion eigen-coefficients
(column with numbers :math:`\le 1` and two decimal points): coefficients smaller than 0.7 indicate that the corresponding quantum number are less reliable. Due to this approximate nature of the TROVE quantum numbers, the TROVE assignment is usually not complete and unambiguous. It is common to find states with duplicate assignments as well as some missing combinations (see Quantum Numbers).


This output section can be reached by searching for ``Zero-point-energy`` (continuing past the basis set sections). This gives the zero-point energy for the vibrational ground state of the molecule, an important quantity. Below this the rotational-vibrational energies for each symmetry are given in order of 'reducing' symmetry.


Rotational-Vibrational energies and quantum numbers
===================================================

The vibrational energies of PF\ :sub:`3` will be given as an example.
::

      Variational solution - irreducible representation
      Gamma  i   value    j  k  t   quanta
      A1 1 0.000000   (A1; 0 0 0)(A1 A1; 0 0 0 0 0 0 ) 0.96 (0 0 0 0 0 0 0) (1 1)
      A1 2 487.299315 (A1; 0 0 0)(A1 A1; 0 0 0 1 0 0 ) 0.86 (0 0 0 1 0 0 0) (1 3)
      A1 3 692.280535 (A1; 0 0 0)(A1 A1; 0 0 0 0 0 2 ) 0.89 (0 0 0 0 0 2 0) (1 4)

In this example, ``Gamma`` is the symmetry, in this case the totally symmetric :math:`A_1` class. ``i`` is just an integer label of the states. ``value`` is the energy of the vibrational levels with respect to the zero point energy in wavenumbers. The rest of the information relates to the eigenfunction of the level.

``j  k  t   quanta`` are related to the rotational states and are discussed below. The next two brackets are the quantum numbers of the state in both normal coordinates and local coordinates used by TROVE. Unless the relations between these quantum numbers have been set up this will not be automatically correct.

The decimal before the second set of quantum numbers gives the certainty of that state consisting of the specified quantum numbers. This is related to the magnitude of the expansion coefficient of this state. For example, here the second row is a fundamental mode of PF\ :sub:`3` while the third row is an overtone with :math:`\nu = 2`. Often states need to be compared to experimental assignments. For vibrational states the total excitation number is usually reliable if not the actual states included.

An example from a :math:`J=2` calculation on PF\ :sub:`3` is shown below.
::

      Variational solution - irreducible representation
      Gamma     i    value       j  k  t   quanta
      E 1  1.157546  (E; 2 2 0) (A1; 0 0 0 0 0 0) 1.00 (0 0 0 0 0 0 0) (1)
      E 2  1.458987  (E; 2 1 0) (A1; 0 0 0 0 0 0) 1.00 (0 0 0 0 0 0 0) (1)
      E 3 347.957388 (E; 2 1 0) (E ; 0 0 0 0 0 1) 1.00 (0 0 0 0 0 1 0) (2)
      E 4 348.255477 (E; 2 2 0) (E ; 0 0 0 0 0 1) 0.73 (0 0 0 0 0 1 0) (2)

In this case the energies are from the doubly degenerate :math:`E` symmetry class. The first two rows are pure rotational states. The ``j k t`` section for these two states are ``2 2 0`` and ``2 1 0`` respectively. This means the total angular momentum is 2 and the projection of the angular momentum onto an axis (usually the :math:`z`-axis is chosen) is 2 and 1 respectively. The third and fourth row are ro-vibrational states with the same vibrational quantum numbers :math:`v_1, v_2,\ldots,`,  but different values of :math:`K`. For a more detailed explanation of the quantum number scheme in TROVE see 
Chapter `Quantum Numbers <https://spectrove.readthedocs.io/en/latest/quantumnumbers.html>`__).

Transition Moment output
========================

The output for a transition moment calculation (for :math:`J=0` only) is similar to the output for intensities discussed below. The section starts at the line
::

     Linestrength S(f<-i) [Debye**2], Transition moments [Debye], ...


A typical output has the following form 
::

        J' G'        J G     Type        E'            E            nu        Gr'   K'  Gv'   v1' v2' v3'      Gr   K    Gv     v1  v2  v3      mu             Int(Tref)         i     n1' n2' n3'         n1  n2  n3           mux             muy               muz                z
        0 A1    <-   0 A1    D       0.000000 <-      0.000000     0.000000  (A1 ;  0) (A1 ;   0   0   0)  <- (A1 ;  0) (A1 ;   0   0   0)  9.70832585E-01  0.00000000E+00       1  (   0   0   0)  <-  (   0   0   0)       0.00000000       0.00000000      -0.97083259
        0 A1    <-   0 A1    D    1172.667646 <-      0.000000  1172.667646  (A1 ;  0) (A1 ;   0   0   1)  <- (A1 ;  0) (A1 ;   0   0   0)  1.08010176E-02  1.92265321E-24       2  (   0   0   1)  <-  (   0   0   0)       0.00000000       0.00000000       0.01080102
        0 A1    <-   0 A1    D    2335.297519 <-      0.000000  2335.297519  (A1 ;  0) (A1 ;   0   0   2)  <- (A1 ;  0) (A1 ;   0   0   0)  3.33036698E-03  4.81671652E-25       3  (   0   0   2)  <-  (   0   0   0)       0.00000000       0.00000000      -0.00333037
        0 A1    <-   0 A1    D    2608.713940 <-      0.000000  2608.713940  (A1 ;  0) (A1 ;   1   0   0)  <- (A1 ;  0) (A1 ;   0   0   0)  4.12201261E-03  8.46925928E-25       4  (   1   0   0)  <-  (   0   0   0)       0.00000000       0.00000000      -0.00412201
        0 B2    <-   0 A1    D    3280.145078 <-      0.000000  3280.145078  (A1 ;  0) (B2 ;   0   1   0)  <- (A1 ;  0) (A1 ;   0   0   0)  1.14535592E-03  8.57021154E-26       5  (   1   0   0)  <-  (   0   0   0)      -0.00114536       0.00000000       0.00000000
        0 A1    <-   0 A1    D    3503.042415 <-      0.000000  3503.042415  (A1 ;  0) (A1 ;   0   0   3)  <- (A1 ;  0) (A1 ;   0   0   0)  8.28791839E-04  4.83387175E-26       6  (   0   0   3)  <-  (   0   0   0)       0.00000000       0.00000000      -0.00082879
        0 A1    <-   0 A1    D    3765.459944 <-      0.000000  3765.459944  (A1 ;  0) (A1 ;   1   0   1)  <- (A1 ;  0) (A1 ;   0   0   0)  7.11292429E-03  3.85777487E-24       7  (   1   0   1)  <-  (   0   0   0)       0.00000000       0.00000000      -0.00711292
        0 B2    <-   0 A1    D    4415.876421 <-      0.000000  4415.876421  (A1 ;  0) (B2 ;   0   1   1)  <- (A1 ;  0) (A1 ;   0   0   0)  1.50378709E-02  2.04819238E-23       8  (   1   0   1)  <-  (   0   0   0)      -0.01503787       0.00000000       0.00000000
        0 A1    <-   0 A1    D    4675.006191 <-      0.000000  4675.006191  (A1 ;  0) (A1 ;   0   0   4)  <- (A1 ;  0) (A1 ;   0   0   0)  1.42747914E-04  1.96021622E-27       9  (   0   0   4)  <-  (   0   0   0)       0.00000000       0.00000000      -0.00014275
        0 A1    <-   0 A1    D    4927.853585 <-      0.000000  4927.853585  (A1 ;  0) (A1 ;   1   0   2)  <- (A1 ;  0) (A1 ;   0   0   0)  5.50415352E-04  3.07955132E-26      10  (   1   0   2)  <-  (   0   0   0)       0.00000000       0.00000000      -0.00055042


and provides the total vibrational transition dipole moment :math:`\bar\my=u = \sqrt{\bar\mu_x^2+\bar\mu_y^2+\bar\mu_z^2}`, the individual components of the transition dipole :math:`\mu_\alpha` as well the vibrational band intensity computed for the reference temperature. 


A list of information on the transition moments between vibrational states is then given. Similar to the output of the rotational-vibrational energy levels, the symmetry and energy of the upper and lower vibrational states is given along with the corresponding vibrational quantum numbers and transition frequency between the states.

The transition moments are printed out along with the line strength. The end of the row shows the values of the transition moment for the x,y and z directions.


Intensity output
================

The intensity output section also starts after the line
::

      Linestrength S(f<-i) [Debye**2], Transition moments [Debye],...

A typical intensity output is given by 
::
     
      J' G'         J G     Type       E'            E          nu        Gr'   K'     Gv'   v1' v2' v3'       Gr    K      Gv    v1  v2  v3       S(f<-i)          A(if)            I(f<-i)             Ni        Nf        N               normal mode          normal mode    S(deg-component)
      1 B1     <-   1 B2     Q    1212.1915 <-     50.2853   1161.9063  ( B1    1 ) (  A1     0   0   1 ) <- ( B2    1 ) (  A1     0   0   0 )    6.87365073E-04   3.75716703E-02   1.06423096E-23        2 <-      1       30        2 (    0   0   1 ) <-  (    0   0   0 )    5.04558873E-03
      1 B1     <-   1 B2     Q    2375.0878 <-     50.2853   2324.8025  ( B1    1 ) (  A1     0   0   2 ) <- ( B2    1 ) (  A1     0   0   0 )    5.09854985E-05   2.23236209E-02   2.09790812E-24        3 <-      1       31        3 (    0   0   2 ) <-  (    0   0   0 )   -1.37417313E-03
      1 B1     <-   1 B2     Q    2647.3320 <-     50.2853   2597.0467  ( B1    1 ) (  A1     1   0   0 ) <- ( B2    1 ) (  A1     0   0   0 )    8.81230563E-05   5.37883956E-02   4.16280721E-24        4 <-      1       32        4 (    1   0   0 ) <-  (    0   0   0 )   -1.80660369E-03
      1 B1     <-   1 B2     Q    3297.7254 <-     50.2853   3247.4401  ( A2    0 ) (  B2     0   1   0 ) <- ( B2    1 ) (  A1     0   0   0 )    1.99024255E-03   2.37514001E+00   1.22488437E-22        5 <-      1       33        5 (    1   0   0 ) <-  (    0   0   0 )    8.58560930E-03
      1 B1     <-   1 B2     Q    3543.1999 <-     50.2853   3492.9147  ( B1    1 ) (  A1     0   0   3 ) <- ( B2    1 ) (  A1     0   0   0 )    2.99417558E-06   4.44632247E-03   2.00135595E-25        6 <-      1       34        6 (    0   0   3 ) <-  (    0   0   0 )   -3.33009597E-04
      1 B1     <-   1 B2     Q    3804.3651 <-     50.2853   3754.0798  ( B1    1 ) (  A1     1   0   1 ) <- ( B2    1 ) (  A1     0   0   0 )    2.22051406E-04   4.09377345E-01   1.60805195E-23        7 <-      1       35        7 (    1   0   1 ) <-  (    0   0   0 )    2.86777373E-03
      1 B1     <-   1 B2     Q    4433.7691 <-     50.2853   4383.4838  ( A2    0 ) (  B2     0   1   1 ) <- ( B2    1 ) (  A1     0   0   0 )    8.72066094E-04   2.55957224E+00   7.46795222E-23        8 <-      1       36        8 (    1   0   1 ) <-  (    0   0   0 )    5.68319841E-03
      1 B1     <-   1 B2     Q    4717.4255 <-     50.2853   4667.1402  ( B1    1 ) (  A1     0   0   4 ) <- ( B2    1 ) (  A1     0   0   0 )    4.98410218E-08   1.76562816E-04   4.56074049E-27        9 <-      1       37        9 (    0   0   4 ) <-  (    0   0   0 )   -4.29646805E-05
      1 B1     <-   1 B2     Q    4967.4402 <-     50.2853   4917.1550  ( B1    1 ) (  A1     1   0   2 ) <- ( B2    1 ) (  A1     0   0   0 )    1.64999991E-06   6.83573742E-03   1.59462746E-25       10 <-      1       38       10 (    1   0   2 ) <-  (    0   0   0 )   -2.47206610E-04
    

It contains state energies, quantum numbers, linestrengths (D\ :sup:`2`), both total an per degenerate component (last column), Einstein A coefficients (1/S), absorption intensities (cm/molecule) for the reference input temperature. 

This section is similar to the transition moment output. The symmetries, quantum numbers and energies of the lower and upper states are given along with the transition frequency. The intensity is given for the transitions along with the line strength and the Einstein A coefficient (see Chapter `Theory <https://spectrove.readthedocs.io/en/latest/theory.html>`__).


Checkpoint File Outputs
=======================

See a detailed description in (see Chapter `Checkpoints  https://spectrove.readthedocs.io/en/latest/checkpoints.html`__. 

As well as the main TROVE output file, useful information is also contained in the descr checkpoint files. These will be described here.

Contr Files
-----------

The contr files describe the details of the contracted functions formed by grouping basis with the same symmetry class.

The file contr-quanta.chk gives the vibrational quantum numbers for the primitive basis functions used for each class of contractions. This is just columns of integers corresponding to the primitive basis functions.

The file contr-descr.chk give the details of the contracted functions themselves. This file first gives some detail on the masses of the atoms and geometries and symmetry of the molecule. This is followed by a summary of how the primitive functions were generated, for example a summary of the Numerov parameters. Details are then given on the contraction. For each class. For example for PF\ :sub:`3` the first class is
::

     Class #       1
     120           120  <-  number of roots and dimension of basis
     1  1  1   1   1954.033595307337   0   0   0   0   0   0   0   0   0   0   0   0   0   0    0.99846636
     2  3  2   1   2817.634700213870   0   1   0   0   0   0   0   0   1   0   0   0   0   0   -0.76056863
     3  3  2   2   2817.634700213870   0   1   0   0   0   0   0   0   1   0   0   0   0   0    -0.76056863
     4  1  3   1   2839.597506890540   0   0   0   1   0   0   0   0   0   0   1   0   0   0    -0.57531184
     5  1  4   1   3672.073371984382   0   2   0   0   0   0   0   0   2   0   0   0   0   0     0.49580488
     6  3  5   1   3676.006458469679   0   2   0   0   0   0   0   0   2   0   0   0   0   0    -0.61014685

The number of roots is the total number of eigenfunctions (contracted basis functions) for this class. This is limited by polyad number or energy cut offs. The rows give details on each contracted function. The energies for the contracted function is then given along with the vibrational quantum numbers of the constituent primitive functions. The final column is the largest coefficient of the linear combination of primitives making up the contracted function.


Eigen files
-----------

The details of the eigenfunctions for the full Hamiltonian are given in the  eigen-descrn-m.chk files where n and m are the :math:`J` and symmetry numbers of the eigenfunctions respectively. This file is very similar to the contr-des files described
above. If the :math:`J=0` method is used then j0eigen-descrn-m.chk files are generated which have the same structure. The j0contr-descr.chk also contains similar information.


