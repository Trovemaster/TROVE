SPAR input syntax and migration
===============================

Declare ``NMODES`` before any ``BASIC-FUNCTION`` block. Supply modes in order,
from 1 to ``NMODES``, and terminate the block with ``END``. Function names are
case insensitive. The KEO and PEF/DMF constructors have different grammars.

Kinetic functions
-----------------

The current grammar is::

   BASIC-FUNCTION KINETIC
   Mode k Nfunctions
   index Nfactors name1 power1 scale1 [name2 power2 scale2 ...]
   ...
   END

``KINETIC`` may be omitted. Use ``MASS`` instead only for a matching KEO
checkpoint with separated inverse mass factors. Each record defines

.. math::

   f^{(k)}_\ell(q_k)=\prod_{s=1}^{N_{\rm factors}}
       [h_s(b_s q_k)]^{p_s}.

The record has exactly two leading integers followed by three fields per
factor. ``index`` is the positive function index; give consecutive, unique
indices from 1 to ``Nfunctions``. Index 0 is implicit and is not listed.
``power`` is an integer (negative powers are allowed); ``scale`` is real.
An angular argument is in radians unless an input field explicitly converts
degrees. The scale multiplies the argument before the function is evaluated.

.. list-table:: Primitive KEO functions
   :class: wrap-table
   :header-rows: 1
   :widths: 24 42 34

   * - Name
     - Function :math:`h(x)`
     - Use
   * - ``I``, ``R``, ``A``
     - :math:`x`
     - Identity aliases
   * - ``SIN``, ``COS``, ``TAN``
     - :math:`\sin x`, :math:`\cos x`, :math:`\tan x`
     - Trigonometric functions
   * - ``CSC``, ``SEC``, ``COT``
     - :math:`1/\sin x`, :math:`1/\cos x`, :math:`1/\tan x`
     - Reciprocal functions
   * - ``CSC_``, ``COT_``
     - Singular forms with special reduced-function pointers
     - Only with the corresponding singularity treatment

For example, these three records define :math:`r^{-2}`,
:math:`\cos^2(\alpha/2)`, and :math:`\cot\alpha\,\csc\alpha`,
respectively. They are separate record examples, not a complete block::

   1 1 R -2 1.0
   1 1 COS 2 0.5
   1 2 COT 1 1.0 CSC 1 1.0

Select ``automatic`` as the kinetic expansion type in the appropriate
``BASIS`` rows when using this constructor. The basis family, coordinate
transformation, frame, KEO, and reference configuration must also agree.
See the complete :doc:`basic_functions_example`.

For ``automatic``, the evaluator adds the mode's reference value to its
displacement before applying the functions: :math:`q_k=\xi_k+q_{k,0}`.
This is why ``R -1 1`` represents an inverse bond length in the tutorial.
The separate ``AUTO`` expansion option evaluates the raw variable without
that shift; these names are not interchangeable coordinate conventions.

Use ``AUTO`` when the coordinate transformation itself returns absolute
values rather than displacements, as ``R-SYMPHI-BETA`` does for the bond
lengths and the umbrella angle (see :doc:`spar_xy3`). ``AUTO-HARMON``
evaluates the functions in the same way and is used for modes with the
``harmonic`` primitive basis. For such modes TROVE evaluates the KEO at the
reference geometry, including the harmonic mode itself, to obtain the
harmonic-oscillator parameters printed on the ``params`` line.

Potential and external-property functions
-----------------------------------------

Use ``BASIC-FUNCTION POTEN`` (alias ``POTENTIAL``) or
``BASIC-FUNCTION EXTERNAL`` (alias ``EXT``). The current grammar is::

   BASIC-FUNCTION POTEN
   Mode k Nfunctions
   Nfactors 0 : pmax primary_name [power2 name2 scale2 inner_power2 ...]
   ...
   END

Each range must start at 0. Put spaces around ``:`` or ``-``. There is no
explicit function-index column. A range from 0 to ``pmax`` contributes
``pmax + 1`` functions, stored consecutively; the total of all ranges must
equal ``Nfunctions``. Subsequent ranges append functions instead of resetting
their storage positions. This matters for an expansion augmented by, for
example, a sine factor.

.. list-table:: Primary function names accepted by the reader
   :class: wrap-table
   :header-rows: 1

   * - Name
     - Family
   * - ``X-XE``, ``(X-XE)``, ``X-X0``
     - Powers of a coordinate displacement
   * - ``MORSE(X-XE)``
     - Powers of a Morse coordinate
   * - ``COS(A)-COS(A0)``, ``COS(A0)-COS(A)``
     - Powers of a cosine difference, with the indicated sign
   * - ``COSNX``, ``SINNX``
     - Cosine or sine series indexed by the range value

The primary-function parameters come from the molecular reference and
expansion settings. Do not append equilibrium values using the portable
SPAR listing syntax. Additional factors accept ``SIN``, ``COS``, ``TAN``,
``CSC``, ``COT``, and ``SEC``; their four-field order differs from the KEO
three-field order.

For example, the current HOPO input defines nine stretching functions with::

   Mode 1 9
   1, 0 - 8, morse(x-xe)

and an 18-function angular library with::

   Mode 4 18
   1, 0 - 8, cos(a)-cos(a0)
   2, 0 - 8, cos(a)-cos(a0), 1 sin 1.0 1

These are constructor excerpts. They do not specify a complete potential,
property model, or runnable HOPO calculation. The checkpoint indices must
follow the corresponding POT/EXT reader's convention; do not substitute KEO
indices just because the function names look similar.

Migrating earlier inputs
------------------------

Older kinetic input records used four fields per factor::

   index Nfactors power name scale inner_power

For ``inner_power = 1`` the syntax translation is, for example::

   ! Earlier input
   2 1 -2 I 1.0 1
   ! Current input
   2 1 I -2 1.0

The current kinetic reader fixes the inner power to 1. A non-unit inner power
requires a mathematical reformulation or another supported representation;
deleting that number does not preserve the function. The old format remains
in some internal Fortran constructor strings handled by a separate helper.
Those strings do not specify the user-input grammar.

Migration also requires checking the mapping. The old H2S manual example used
six bending functions, while the current compact XY2 KEO uses seven in a
different order. Translating the fields alone would leave the operator
incorrect. Use the reconciled :doc:`basic_functions_example`.

The PEF/DMF listings in the paper's supplementary PDF use names such as
``Morse`` and ``Power``, place the name before the power range, and show
explicit structural parameters. They are not accepted verbatim by this
revision's reader. Treat them as model specifications and translate both the
constructor and its checkpoint mapping. The POT/EXT examples above document
the current reader; a full numerical migration of the published PEF/DMF
examples remains a separate validation task.
