SPAR: operator representations
==============================

SPAR (Summed Product Analytic Representation) describes molecular operators
through one-dimensional **basic functions**, numerical coefficients, and a
mapping between them [26SmMeYu]_. It provides a common representation for kinetic
energy operators (KEOs), potential energy functions (PEFs), and dipole moment
functions (DMFs). The representation can be exact, approximate, or hybrid.

The pages below describe the input reader in TROVE revision ``7831d55``
(7 September 2026). Published examples and older inputs can use different
constructor syntax; see :doc:`spar_input` before adapting them.

.. toctree::
   :maxdepth: 2

   spar_input
   spar_checkpoints
   basic_functions_example

From an operator to one-dimensional integrals
---------------------------------------------

For a scalar function, or one component of a tensor-valued function, write

.. math::
   :label: spar-sop

   F(q_1,\ldots,q_D)=\sum_{L=1}^{N_F} a_L
      \prod_{k=1}^{D} f^{(k)}_{\ell_k(L)}(q_k).

Here :math:`D` is the number of vibrational coordinates, :math:`L` labels a
product term, and :math:`\ell_k(L)` selects a basic function for mode
:math:`k`. A zero mode index represents the constant function 1 in the KEO
mapping. The coefficient and the complete index tuple belong together:
reordering the functions without updating the tuples changes the operator.

The operator basic functions :math:`f^{(k)}_\ell` are distinct from the
variational basis functions :math:`\phi^{(k)}_n`. For a multiplication operator
in a product basis, the matrix element factors as

.. math::

   \langle\boldsymbol n|F|\boldsymbol m\rangle
   =\sum_L a_L\prod_k
     \langle\phi^{(k)}_{n_k}|f^{(k)}_{\ell_k(L)}|
       \phi^{(k)}_{m_k}\rangle.

KEO matrix elements also contain derivatives and angular-momentum operators.
Their placement and the associated pseudopotential must follow the chosen
Hamiltonian convention; specifying a set of scalar functions alone does not
define a KEO. TROVE stores separate vibrational, rotational, Coriolis, and
pseudopotential contributions (see :doc:`spar_checkpoints`).

Accuracy and convergence
------------------------

An exact analytic KEO avoids an expansion error in that operator. It does not
remove errors from a fitted potential, numerical quadrature, a finite
variational basis, or contracted-state truncation. A calculation can combine
an exact KEO with an approximate PEF or DMF. Record the origin and approximation
of each component separately.

For a Taylor representation, a mode index may coincide with a power. For a
general SPAR representation it is a function label: increasing ``KinOrder``
does not automatically improve an exact finite function library. Convergence
still requires systematic changes to the basis, integration grids, and
contraction settings. See :doc:`basissets` and :doc:`quickstart`.

Mass dependence
---------------

``BASIC-FUNCTION MASS`` selects the KEO reader with an additional mass index
:math:`\ell_0`. Each term then contains an inverse nuclear mass selected from
the input ``ZMAT`` order. This permits coefficients with the mass factors
separated from their coordinate dependence. ``BASIC-FUNCTION KINETIC`` uses
coefficients with the masses already incorporated.

These two checkpoint conventions have different columns and unit conversion
behaviour. They cannot be interchanged by changing the block heading alone.
Furthermore, a frame or coordinate definition may itself depend on masses.
Check that definition before transferring a model to another isotopologue.

Singular coordinates
--------------------

Factors such as :math:`1/\sin\rho` require a compatible basis and integration
scheme at a linear geometry. TROVE recognises ``CSC_`` and ``COT_`` markers
for special matrix-element treatment, but they are not general replacements
for ``CSC`` and ``COT``. The basis, operator, and pseudopotential must be derived
and used together. The paper's supplementary material describes these cases.

The bent H2S tutorial stays away from linear configurations and uses ordinary
``CSC`` and ``COT`` functions. It does not validate singular-coordinate models.

Preparing a transferable model
------------------------------

Keep the following information with every operator file:

* Coordinate definitions, order, units, and reference geometry.
* Body-fixed axes, atom order, masses, and symmetry conventions.
* Function libraries and index mappings for each operator or property.
* Coefficient units, mass convention, differential ordering, integration
  measure, and pseudopotential convention.
* Potential and dipole parameters, their source, and their range of validity.
* TROVE revision, build choices, input files, and numerical comparison data.

An external operator does not supply a missing molecular coordinate
transformation, symmetry implementation, or dipole-frame transformation.
These remain prerequisites for the calculation.

Sources and examples
--------------------

The methodological reference is [26SmMeYu]_. The associated
`supplementary archive <https://doi.org/10.6084/m9.figshare.33197481>`__ contains
``InputsAndFrames.pdf``, derivations, function definitions, and operator files.
The public `SPAR repository <https://github.com/Trovemaster/SPAR>`__ contains
further examples and conversion tools. Record the repository commit when
using them; the input conventions must be reconciled with the TROVE revision.

The constructor description in this manual was checked against ``FLReadInput``
in ``fields.f90`` and the current HOPO input in ``inputs/HOPO``. The compact
XY2 example was checked against ``MLkinetic_compact_xy2_bisect_EKE_rigid`` in
``kin_xy2.f90``. The scientific interpretation of a new molecular model still
requires validation of its coordinates, operators, and observables.
