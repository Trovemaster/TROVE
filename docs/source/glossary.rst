.. _glossary:

Glossary
========

.. glossary::
   :sorted:

   KinOrder
      Expansion order of the KEO


   PotOrder
      Expansion order of the PEF.


   Natoms
      Number of atoms (nuclei).

   Nmodes
      Number of vibrational coordinates; normally 3N-6 in the nonlinear formulation.

   SYMGROUP
      Molecular symmetry group.

   verbose
      Verbosity level controlling the standard output.

   dstep
      Numerical difference step size (angstroms or radians).

   Coords
      Type of the coordinate, ``linear`` (``linearised``) or ``local`` (``curvilinear``). 

   Frame
      Molecular frame.

   Zmat
     Z-matrix block defining the Z-matrix coordinates and nuclear (atomic) masses.


   print
      Block with printing options

   Nassignments
     (alias ``N_EIGEN-CONTRIBUTIONS``) defines the number of the assignments to generate.
      

   exp_coeff_thresh
     Expansion coefficients of  the field checkpoints ``kinetic.chk``, ``potential.chk`` and ``external.chk`` that are smaller by magnitude than this threshold are not included in the corresponding checkpoint.

   Moltype
     The type of molecule (XYZ, XY2, XY3, XY4, ZXY3, etc).

   REFER-CONF
      reference configuration, ``RIGID`` or ``NON-RIGID``.

   PRIMITIVES
      block defining parameters of the primitive bases.

   Npolyads   
      Polyad cutoff in the PRIMITIVES or CONTRACTION block, respectively.

   CONTRACTION
      Block defining parameters of the contracted basis set.


   sample_points
      number of sampling points in the symmetrisation procedure.

   sample_attempts
      number of symmetrisation attempts.

   symm_toler
      Numerical tolerance used in symmetrisation.

   DIAGONALIZER
      Block defining the diagonaliser (eigensolver) as well as its options (number of roots, maximal energy etc).

   SYEV
      LAPACK Eigensolver type DSYEV.

   enermax
      Maximal energy (cm\ :sup:`-1`).


   control
      Control block (see :doc:`quickstart`).

   Basis
      Basis set block (See :doc:`basissets`).

   EQUILIBRIUM
     Equilibrium values of the molecule geometries in terms of the Z-matrix coordinates.

   SPECPARAM
     Special parameters used to define the coordinate to expand PEF, e.g. the Morse parameter :math:`a`.

   POTEN
      Potential block (see :doc:`potentials`).

   DIPOLE
      Dipole moment block (or ``external`` field block)
