Introduction
************
.. _intro:

What is TROVE?
==============


TROVE (Theoretical ROtational Vibrational Energies) is a suite of programs primarily designed for the calculation of molecular infrared line lists [TROVE]_.
It was (and is) developed by Sergey Yurchenko with contributions from others over the years. It has been used to study a range of molecules and features continue to be added.

The main capabilities of TROVE are as follows: The calculation of rotational-vibrational energy levels of polyatomic molecules; Refinement of *ab initio* potential energy surfaces (PES) using experimental data; The calculation of infrared transition intensities and absorption cross sections.

The philosophy of TROVE is to enable the study of molecules of arbitrary structure. This is possible since TROVE uses a numerically generated Hamiltonian unlike other programs which are hard coded with specific analytical Hamiltonians. So far TROVE has been used for  triatomics  [TROVE]_, tetratomics,
[09YuBaYa]_, [11YuBaTe]_, [13SoYuTe]_, [13UnTeYu]_, [15AlOvPo]_, [15SoAlTe]_,  pentatomics\cite [14YuJo]_,[15OwYuTa]_,[17OwYuYa]_ and hexatomic molecules [18MaYaTe]_ and there are plans to implement even larger molecules. The structure of TROVE also allows very straightforward studies of new molecules of the same structure and symmetry to molecules which have previously been implemented. For example, since phosphine, PH:sub:`3` is implemented, it is simple to then carry out calculations on arsine, AsH:sub:`3`, provided a PES and dipole moment surface (DMS) is available.

TROVE calculates rotational-vibrational energy levels using the variational method. In very simple terms this method requires a suitable choice of functions called a basis set, the calculation of a Hamiltonian matrix using the basis set and the diagonalisation of the Hamiltonian. The eigenvalues are the energy levels for the system while the eigenfunctions are the states/wavefunctions of the system. TROVE carries out this procedure in a very sophisticated manner by building up basis sets from calculations on smaller parts of the system and makes use of symmetry to simplify the calculation of matrix elements and reduce the size of matrices to be diagonalised. Despite the complexity of the implementation in TROVE, the basic procedure should be kept in mind.

What is this manual for?
========================

This aim of this manual is to describe the theory behind TROVE and especially how to operate TROVE in more detail than would be possible in a journal publication. Since the original paper, many new features have been added to TROVE which have been described as they were applied. In this manual these new features will be described in one place in more detail.

TROVE is very much a 'living' program with new capabilities being added all the time. This may mean that some parts of this manual are not be up to date or even superfluous for the reader. Nevertheless, it is anticipated that the core parts of TROVE will not change significantly and that this document will still be of use, especially to new TROVE users. The last chapter describes new features of TROVE and will be kept up to date with new developments.


The 'Quick start' chapter following this one is aimed at new users, especially new PhD students and postdocs, to get  started doing calculations with TROVE assuming some set up has already been done for the molecule of interest. Later chapters describe the underlying theory behind TROVE, a detailed description of outputs, the refinement procedure, the calculation of line lists, a list of molecules already implemented in TROVE and summary of work on them, how to set up TROVE for a new molecule and finally, new features in TROVE.


Main Contributors 
=================

Sergey Yurchenko
Andrey Yachmenev
Thomas Mellor 
Oleksiy Smola
Barry Mant 
Roman Ovsyannikov
Per Jensen 
Walter Thiel 



