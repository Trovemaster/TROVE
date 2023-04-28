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


Rotational-Vibrational energies
===============================

The final step of a vibrational or rotational calculation is an output of the rotational-vibrational energies. These are ordered by energy and separated into symmetry blocks.

This output section can be reached by searching for ``Zero-point-energy`` (continuing past the basis set sections). This gives the zero-point energy for the vibrational ground state of the molecule, an important quantity. Below this the rotational-vibrational energies for each symmetry are given in order of 'reducing' symmetry.

The vibrational energies of PF\ :sub:`3` will be given as an example.
::

      Variational solution - irreducible representation
      Gamma  i   value    j  k  t   quanta
      A1 1 0.000000   (A1; 0 0 0)(A1 A1; 0 0 0 0 0 0 ) 0.96 (0 0 0 0 0 0 0) (1 1)
      A1 2 487.299315 (A1; 0 0 0)(A1 A1; 0 0 0 1 0 0 ) 0.86 (0 0 0 1 0 0 0) (1 3)
      A1 3 692.280535 (A1; 0 0 0)(A1 A1; 0 0 0 0 0 2 ) 0.89 (0 0 0 0 0 2 0) (1 4)

In this example, ``Gamma`` is the symmetry, in this case the totally symmetric :math:`A_1` class. ``i`` is just an integer label of the states. ``value`` is the energy of the vibrational levels with respect to the zero point energy in wavenumbers. The rest of the information relates to the eigenfunction of the level.

``j  k  t   quanta`` are related to the rotational states and are discussed below. The next two brackets are the quantum numbers of the state in both normal coordinates and local coordinates used by TROVE. Unless the relations between these quantum numbers have been set up this will not be automatically correct.

The decimal before the second set of quantum numbers gives the certainty of that state consisting of the specified quantum numbers. This is related to the magnitude of the expansion coefficient of this state. For example, here the second row is a a fundamental mode of PF\ :sub:`3` while the third row is an overtone with :math:`\nu = 2`. Often states need to be compared to experimental assignments. For vibrational states the total excitation number is usually reliable if not the actual states included.

An example from a :math:`J=2` calculation on PF\ :sub:`3` is shown below.
::

      Variational solution - irreducible representation
      Gamma     i    value       j  k  t   quanta
      E 1  1.157546  (E; 2 2 0) (A1; 0 0 0 0 0 0) 1.00 (0 0 0 0 0 0 0) (1)
      E 2  1.458987  (E; 2 1 0) (A1; 0 0 0 0 0 0) 1.00 (0 0 0 0 0 0 0) (1)
      E 3 347.957388 (E; 2 1 0) (E ; 0 0 0 0 0 1) 1.00 (0 0 0 0 0 1 0) (2)
      E 4 348.255477 (E; 2 2 0) (E ; 0 0 0 0 0 1) 0.73 (0 0 0 0 0 1 0) (2)

In this case the energies are from the doubly degenerate E symmetry class. The first two rows are pure rotational states. The ``j k t`` section for these two states are ``2 2 0`` and ``2 1 0`` respectively. This means the total angular momentum is 2 and the projection of the angular momentum onto an axis (usually the z-axis is chosen) is 2 and 1 respectively. The third and fourth row are ro-vibrational states with the same vibrational quantum numbers but different values of :math:`k`.


Transition Moment output
========================


The output for a transition moment calculation is similar to the output for intensities discussed below. The section starts at the line
::

     Linestrength S(f<-i) [Debye**2], Transition moments [Debye], ...


A list of information on the transition moments between vibrational states is then given. Similar to the output of the rotational-vibrational energy levels, the symmetry and energy of the upper and lower vibrational states is given along with the corresponding vibrational quantum numbers and transition frequency between the states.

The transition moments are printed out along with the line strength. The end of the row shows the values of the transition moment for the x,y and z directions.


Intensity output
================

The intensity output section also starts after the line
::

      Linestrength S(f<-i) [Debye**2], Transition moments [Debye],...


This section is similar to the transition moment output. The symmetries, quantum numbers and energies of the lower and upper states are given along with the transition frequency. The intensity is given for the transitions along with the line strength and the Einstein A coefficient (see Chapter `Theory <https://spectrove.readthedocs.io/en/latest/theory.html>`__).


Checkpoint File Outputs
=======================


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


