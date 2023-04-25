Introduction
============
.. _sec-intro:

What is TROVE?
--------------


TROVE (Theoretical ROtational Vibrational Energies) is a suite of programs primarily designed for the
calculation of molecular infrared line lists [1]_.
It was (and is) developed by Sergey Yurchenko with contributions from others over the years.
It has been used to study a range of molecules and features continue to be added.

The main capabilities of TROVE are as follows: The calculation of rotational-vibrational energy levels of polyatomic molecules;
Refinement of *ab initio* potential energy surfaces (PES) using experimental data;
The calculation of infrared transition intensities and absorption cross sections.

The philosophy of TROVE is to enable the study of molecules of arbitrary structure.
This is possible since TROVE uses a numerically generated Hamiltonian unlike other programs which are hard coded with
specific analytical Hamiltonians.
So far TROVE has been used for  triatomics  [1]_, tetratomics,
[2]_-[8]_
pentatomics\cite [9]_,[10]_,[11]_ and hexatomic molecules [12]_ and there are plans to implement even larger molecules.
The structure of TROVE also allows very straightforward studies of new molecules of the same structure and symmetry to
molecules which have previously been implemented. For example, since phosphine, PH:sub:`3` is implemented, it is simple to then
carry out calculations on arsine, AsH:sub:`3`, provided a PES and dipole moment surface (DMS) is available.

TROVE calculates rotational-vibrational energy levels using the variational method.
In very simple terms this method requires a suitable choice of functions called a basis set,
the calculation of a Hamiltonian matrix using the basis set and the diagonalisation of the Hamiltonian.
The eigenvalues are the energy levels for the system while the eigenfunctions are the states/wavefunctions of the system.
TROVE carries out this procedure in a very sophisticated manner by building up basis sets from calculations on smaller
parts of the system and makes use of symmetry to simplify the calculation of matrix elements and reduce the size of
matrices to be diagonalised. Despite the complexity of the implementation in TROVE, the basic procedure should be
kept in mind.

What is this manual for?
------------------------

This aim of this manual is to describe the theory behind TROVE and especially how to operate TROVE in more detail than
would be possible in a journal publication. Since the original paper, many new features have been added to TROVE which have
been described as they were applied. In this manual these new features will be described in one place in more detail.

TROVE is very much a `living' program with new capabilities being added all the time. This may mean that some parts of
this manual are not be up to date or even superfluous for the reader. Nevertheless, it is anticipated that the core
parts of TROVE will not change significantly and that this document will still be of use, especially to new
TROVE users. The last chapter describes new features of TROVE and will be kept up to date with new developments.


The 'Quick start' chapter following this one is aimed at new users, especially new PhD students and postdocs, to get
 started doing calculations with TROVE assuming some set up has already been done for the molecule of interest.
Later chapters describe the underlying theory behind TROVE, a detailed description of outputs, the refinement procedure,
the calculation of line lists, a list of molecules already implemented in TROVE and summary of work on them,
how to set up TROVE for a new molecule and finally, new features in TROVE.


Citations
---------


.. [1] S. N. Yurchenko, W. Thiel, P. Jensen, J. Mol. Spectrosc. 245, 126 (2007), Theoretical ROVibrational Energies (TROVE): A robust numerical approach to the calculation of
rovibrational energies for polyatomic molecules. 

.. [2]_ S. N. Yurchenko, R. J. Barber, A. Yachmenev, W. Thiel, P. Jensen, J. Tennyson, J. Phys.
Chem. A 113, 11845 (2009), A variationally computed T=300 K line list for NH3.

.. [3]_ S. N. Yurchenko, R. J. Barber, J. Tennyson, Mon. Not. R. Astron. Soc. 413, 1828 (2011),
A variationally computed hot line list for NH3.

.. [4]_ D. S. Underwood, J. Tennyson, S. N. Yurchenko, Phys. Chem. Chem. Phys. 15, 10118
(2013), An ab initio variationally computed room-temperature line list for SO3.

.. [5]_ O. L. Polyansky, I. N. Kozin, P. Malyszek, J. Koput, J. Tennyson, S. N. Yurchenko, J.
Phys. Chem. A 117, 7367 (2013), Variational calculation of highly excited rovibrational
energy levels of H2O2.

.. [6]_ C. Sousa-Silva, S. N. Yurchenko, J. Tennyson, J. Mol. Spectrosc. 288, 28 (2013), A computed
room temperature line list for phosphine.

.. [7]_ A. F. Al-Refaie, S. N. Yurchenko, A. Yachmenev, J. Tennyson, Mon. Not. R. Astron.
Soc. 448, 1704 (2015), ExoMol line lists - VIII: A variationally computed line list for hot
formaldehyde.

.. [8]_ C. Sousa-Silva, A. F. Al-Refaie, J. Tennyson, S. N. Yurchenko, Mon. Not. R. Astron. Soc.
446, 2337 (2015), ExoMol line lists - VII. the rotation-vibration spectrum of phosphine
up to 1500 K.

.. [9]_ S. N. Yurchenko, J. Tennyson, Mon. Not. R. Astron. Soc. 440, 1649 (2014), ExoMol line
lists IV: The rotation-vibration spectrum of methane up to 1500 K.

.. [10]_ A. Owens, S. N. Yurchenko, A. Yachmenev, J. Tennyson, W. Thiel, J. Chem. Phys. 142,
244306 (2015), Accurate ab initio vibrational energies of methyl chloride.

.. [11]_ A. Owens, S. N. Yurchenko, A. Yachmenev, W. Thiel, J. Tennyson, Mon. Not. R. Astron.
Soc. 471, 5025 (2017), ExoMol molecular line lists XXII. The rotation-vibration spectrum
of silane up to 1200K.

.. [12]_ B. P. Mant, A. Yachmenev, J. Tennyson, S. N. Yurchenko, Mon. Not. R. Astron. Soc. 478,
3220 (2018), ExoMol molecular line lists - XXVII: spectra of C2H4.

