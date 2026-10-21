SPAR operator files and conventions
===================================

A ``BASIC-FUNCTION`` block defines the coordinate functions. The operator
file supplies their coefficients and index tuples. Both are needed to
reconstruct the model. TROVE has several checkpoint readers, selected by
the input options; a filename ending in ``.chk`` does not identify its format.

KEO components
--------------

The kinetic checkpoint contains four consecutive sections:

.. list-table:: Meaning of the component indices
   :class: wrap-table
   :header-rows: 1
   :widths: 22 24 54

   * - Section
     - Indices
     - Meaning
   * - ``Gvib``
     - ``k1, k2``
     - Two vibrational modes
   * - ``Grot``
     - ``k1, k2``
     - Two body-fixed rotational axes (1, 2, 3)
   * - ``Gcor``
     - ``k1, k2``
     - Vibrational mode followed by rotational axis
   * - Pseudopotential
     - ``k1, k2``
     - Scalar contribution; retain the file's placeholder fields

In the mode-indexed KEO reader, a data record has the layout::

   k1 k2 term grid coefficient l1 l2 ... lD

In the mass-factorised reader the layout is::

   k1 k2 term grid coefficient l0 l1 l2 ... lD

``term`` is local to the component. Keep component records together.
``grid`` is 0 for a rigid-reference file; non-rigid grid handling requires
the matching reference-configuration setup. A zero coordinate index selects
the constant function. In the mass-factorised format ``l0`` is a positive
atom index in the ``ZMAT`` mass order; it is not a coordinate exponent.

Each section has a header and ends with a sentinel record beginning with
``987654321``. Preserve the remaining sentinel fields as well: the reader
reads the whole record before checking the sentinel. Keep the final
threshold record and footer, where present. The abbreviated tables in the
paper omit some native checkpoint fields and are not complete files.

Units and mass factors
----------------------

The ordinary mode-indexed reader uses the supplied coefficients directly
in TROVE's internal energy convention. A checkpoint written by the matching
TROVE calculation is the safest starting point for this route.

The mass-factorised reader multiplies each coefficient by an internal
conversion factor divided by ``mass(l0)``. The factor in ``fields.f90`` is

.. math::

   C=\frac{h N_{\rm A}\,10^{16}}{4\pi^2 c},

with the physical constants in TROVE's own unit convention. The conversion
is consistent with distances in angstroms, angular coordinates in radians,
nuclear masses in atomic mass units, and output energies in
:math:`\mathrm{cm}^{-1}`. Do not multiply a saved internal-unit checkpoint
by this factor a second time. Compare a reconstructed operator at several
geometries before using converted coefficients for a spectrum.

Selecting a reader
-------------------

An external mode-indexed KEO commonly uses the following checkpoint option::

   CHECK_POINT
     kinetic read with_modes
   END

This is an excerpt, not a complete restart specification. Pair it with the
appropriate ``BASIC-FUNCTION`` block, ``KINETIC`` settings, coordinates,
frame, and basis. ``BASIC-FUNCTION MASS`` selects the additional mass column.
See :doc:`checkpoints` for the other checkpoint objects and the staged
TROVE workflow. Reusing wavefunctions or matrix elements requires compatible
dimensions and conventions as well as compatible files.

Do not retain ``KINETIC COMPACT`` when intending to read an external file:
at this revision that option unconditionally changes ``KINETIC READ`` to
``SAVE`` for the built-in construction route. In the H2S round-trip test,
the entire built-in ``KINETIC`` block is removed, the explicit function
library is retained, and ``kinetic read with_modes`` is specified after
the ``CONTROL`` block. The saved ``kinetic.chk`` is copied into the fresh
run directory. This exercises the external reader instead of rebuilding
the operator. See :doc:`basic_functions_example`.

Potential and dipole files
--------------------------

The potential and external-property readers have their own record formats
and function-index conventions. External-property records include a component
index; the meaning of that component depends on the chosen body-fixed frame
and property representation. The model must specify that mapping explicitly.
The compact tables in the SPAR supplementary material document the model
coefficients, while :doc:`spar_input` describes the current TROVE constructor.
Neither replaces a complete native checkpoint example.

Validation before reuse
-----------------------

1. Check atom, coordinate, and component order against the derivation.
2. Check every tuple against its function library and mass convention.
3. Evaluate the reconstructed functions at non-singular geometries and
   compare them with independent analytic expressions or reference data.
4. Run a small calculation from a clean directory. Inspect the program's
   completion message and selected levels; a zero process exit status alone
   does not prove that a Fortran program completed.
5. Check a nonzero-J case before claiming rotational/Coriolis validation,
   and a transition observable before claiming intensity validation.

The H2S example tests the explicit kinetic-function library and the ordinary
mode-indexed KEO checkpoint reader. It is separate from validation of the
mass-factorised and external PEF/DMF readers.
