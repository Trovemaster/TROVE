.. _glossary:

Glossary
========

.. glossary::

   KinOrder
      Expansion order of the KEO


   PotOrder
      Expansion order of the PEF.


   Natoms: Number of atoms (nuclei) :math:`N`.
   Nmodes: Number of modes or degrees of freedom :math:`M` (here :math:`M=3N-6`).
   SYMGROUP: Molecular symmetry group.
   verbose: Verbosity level controlling amount of information in the standard output.
    ``dstep``: numerical difference step size used in finite differences (Angstrom or radian).
      Stand-alone keywords

   For more information, refer to :doc:`quickstart`.


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
      Maximal number of polyads.

   CONTRACTION
      Block defining parameters of the contracted basis set.

   Npolyads
      Maximal number of polyads in the contracted basis.

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
      Control block (see **Quick start**).

   Basis
      Basis set block (See **Basis sets**).

   EQUILIBRIUM
     Equilibrium values of the molecule geometries in terms of the Z-matrix coordinates.

   SPECPARAM
     Special parameters used to define the coordinate to expand PEF, e.g. the Morse parameter :math:`a`.

   POTEN
      Potential block (see **Potential energy functions**).

   DIPOLE
      Dipole moment block (or ``external`` field block)


