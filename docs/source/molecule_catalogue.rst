Molecular models and example status
===================================

A molecular formula does not uniquely specify a TROVE model. Coordinates,
frame, kinetic representation, surface parameters, and basis settings must
be considered together. This catalogue starts with three models that help
explain those choices; the broader historical coverage is in :doc:`molecules`.

The entries below distinguish code/input inspection from numerical
validation. Production ExoMol settings can require substantially more
resources than a teaching example.

Hydrogen sulfide: H2S
---------------------

The :doc:`basic_functions_example` uses ``MOLTYPE XY2``, ``C2v(M)`` symmetry,
valence coordinates ``R-ALPHA``, the ``BISECT-Z`` frame, and a rigid reference
configuration. Its KEO is the built-in compact analytic
``KINETIC_XY2_EKE_BISECT_COMPACT_RIGID`` representation. The nuclear mass order
is S, H, H. This is a bent-coordinate example with integration bounds away
from the linear singularities.

The supplied input includes the ``POTEN_XY2_TYUTEREV`` parameter list and an
``XY2_PQ_COEFF`` dipole parameter list. The tutorial tests the kinetic
constructor and small-basis energy calculation; it does not establish an
intensity benchmark or reproduce a production ExoMol line list. See the
tutorial for the precise numerical validation scope and commands.

The associated code is in ``mol_xy2.f90``, ``kin_xy2.f90``, and
``pot_xy2.f90``. The downloadable input is
:download:`H2S_EKE_basic-functions_step1.inp <input/H2S_EKE_basic-functions_step1.inp>`.
Its small primitive and contracted cutoffs are teaching settings, not a
statement of spectroscopic convergence.

Silylene: SiH2 / CATS
---------------------

The repository's ``inputs/SiH2`` directory contains CATS inputs for steps
1 and 2, step 3 at several J values, and a step-4 intensity calculation.
The step-1 configuration uses ``MOLTYPE XY2``, ``C2v(M)``, ``COORDS LINEAR``,
``TRANSFORM R-RHO-Z``, and ``REFER-CONF NON-RIGID``. Here ``LINEAR`` denotes
linearised coordinates; it does not mean that the molecule is linear.
The mass order is Si, H, H.

The PEF is ``POTEN_XY2_MORSE_COS``. The input specifies ``KinOrder 6`` and
``PotOrder 8`` and uses Numerov primitive functions. This non-rigid,
linearised-coordinate model is different from the rigid-reference H2S
SPAR tutorial, even though both use the XY2 molecular implementation.

The `CATS publication <https://arxiv.org/abs/2004.01005>`__ describes the
production calculation. The manual also supplies
:download:`a SiH2 example <input/SiH2_XY2_MORSE_COS_step1.inp>`; its settings
must not be assumed identical to every CATS stage. These inputs have been
inspected for configuration, but have not been rerun as part of this
documentation update. Later stages require the checkpoints produced by the
earlier stages.

Thioformaldehyde: H2CS / MOTY
-----------------------------

The MOTY model [23MeOwTe]_ uses the ``ZXY2`` molecular implementation and
``R-THETA-TAU`` coordinates with a rigid reference configuration. Its six
coordinates comprise three bond lengths, two angles, and a torsion. See
:doc:`kinetic` for the function ordering and :doc:`frames` for the embedding.
Do not infer the nuclear mass order from the molecular formula; retain the
model's ``ZMAT`` order.

The `MOTY model archive <https://exomol.com/models/H2CS/1H2-12C-32S/MOTY/>`__
and the supplementary material of [23MeOwTe]_ provide the model provenance
and kinetic derivation. The inspected MOTY input requests an external
kinetic checkpoint. The model page also links the companion ``kinetic.chk``.
Obtain that file together with the input and check its function mapping and
derivation conventions before running it. The page identifies the original
executable as ``j-trove2509.x (2022)``; compatibility with a later revision
must be checked. This documentation update does not claim a reproduced
MOTY spectrum.

The repository also contains ``inputs/H2CS/H2CS_Andrey_PES2_rigid_01b_basic.inp``
with ``POTEN_H2CS_DAMP``, while the downloadable
:download:`H2CS_zxy2_morse_cos_step1.inp <input/H2CS_zxy2_morse_cos_step1.inp>`
uses ``POTEN_ZXY2_MORSE_COS``. These are different model specifications;
sharing a formula and coordinate family does not make the surfaces
interchangeable. Both remain reference inputs pending complete numerical
validation and dependency packaging.

Adding another entry
--------------------

Use the following information for each model:

* Formula, isotopologue, electronic-state scope, citation, and model name.
* Molecular family, symmetry, atom order, masses, coordinates, and frame.
* KEO approximation and source; potential and property functional forms.
* Full input and all external parameters, checkpoints, or compiled modules.
* Primitive and contracted basis settings and numerical integration bounds.
* Stage dependencies, representative resource use, and expected outputs.
* Validation state and evidence: inspected, run complete, numerically
  compared, or tutorial ready; specify the tested observables and tolerances.

Retain the original production specification when deriving a smaller
example, and state every change that affects the physical model or its
numerical approximation.
