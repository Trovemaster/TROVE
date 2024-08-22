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


    print:  Name of the block
    NASSIGNMENTS: (alias ``N_EIGEN-CONTRIBUTIONS``) defines the number of the assignments to generate.
      Block with printing options

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


    document name
      Since reStructuredText source files can have different extensions
      (some people like ``.txt``, some like ``.rst`` -- the extension can be
      configured with :confval:`source_suffix`)
      and different OSes have different path
      separators, Sphinx abstracts them: :dfn:`document names` are always
      relative to the :term:`source directory`, the extension is stripped, and
      path separators are converted to slashes.  All values, parameters and such
      referring to "documents" expect such document names.

      Examples for document names are ``index``, ``library/zipfile``, or
      ``reference/datamodel/types``.  Note that there is no leading or trailing
      slash.

   domain
      A domain is a collection of markup (reStructuredText :term:`directive`\ s
      and :term:`role`\ s) to describe and link to :term:`object`\ s belonging
      together, e.g. elements of a programming language.  Directive and role
      names in a domain have names like ``domain:name``, e.g. ``py:function``.

      Having domains means that there are no naming problems when one set of
      documentation wants to refer to e.g. C++ and Python classes.  It also
      means that extensions that support the documentation of whole new
      languages are much easier to write.

      For more information, refer to :doc:`/usage/domains/index`.

   environment
      A structure where information about all documents under the root is saved,
      and used for cross-referencing.  The environment is pickled after the
      parsing stage, so that successive runs only need to read and parse new and
      changed documents.

   extension
     A custom :term:`role`, :term:`directive` or other aspect of Sphinx that
     allows users to modify any aspect of the build process within Sphinx.

     For more information, refer to :doc:`/usage/extensions/index`.

   master document
      The document that contains the root :rst:dir:`toctree` directive.

   root document
      Same as :term:`master document`.

   object
      The basic building block of Sphinx documentation.  Every "object
      directive" (e.g. :rst:dir:`py:function` or :rst:dir:`object`) creates such
      a block; and most objects can be cross-referenced to.

   RemoveInSphinxXXXWarning
      The feature which is warned will be removed in Sphinx-XXX version.
      It usually caused from Sphinx extensions which is using deprecated.
      See also :ref:`when-deprecation-warnings-are-displayed`.

   role
      A reStructuredText markup element that allows marking a piece of text.
      Like directives, roles are extensible.  The basic syntax looks like this:
      ``:rolename:`content```.  See :ref:`rst-inline-markup` for details.

   source directory
      The directory which, including its subdirectories, contains all source
      files for one Sphinx project.

   reStructuredText
      An easy-to-read, what-you-see-is-what-you-get plaintext markup syntax and
      parser system.

