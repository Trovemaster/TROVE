Manual validation and remaining work
===================================

This update establishes the first batch of the manual-improvement plan:
the documentation build, immediate scientific and syntax corrections,
the SPAR foundation, a reproducible H2S example, and initial catalogue entries.
It is not a validation of every molecular implementation or archived input.

Documentation
-------------

The unmodified baseline at ``7831d55f39662f4a4b2782e89a2011929cc25b01``
produced 132 Sphinx warnings, including malformed markup and broken local
references. The updated HTML builds with warnings treated as errors using
Sphinx 8.2.3 and sphinx-rtd-theme 3.0.2. The dependencies are pinned in
``docs/requirements.txt`` and the Read the Docs configuration now fails on
warnings. MathJax replaces the HTML build's LaTeX dependency.

Corrections include the rotational statistical factor ``2J+1``, the H2CS
mode-function descriptions, the secant keyword, four missing download paths,
and the SiH2 molecule heading. The duplicated XY2 walkthrough is consolidated
while retaining its old page URL. Tables, directives, citations, and internal
references were repaired as indicated by the build diagnostics.

The SPAR material distinguishes the current input grammar from the
supplementary publication's PEF/DMF listings. The supplementary archive and
public SPAR repository were inspected; their identifiers are recorded in
``spar-2026-10.json``. No third-party operator files were silently rewritten.

Numerical validation
--------------------

See ``test/spar_h2s/README.rst`` for the runnable checks, build choices, and
known baseline issues. The machine-readable results are in
``spar-2026-10.json``. The H2S tutorial includes the numerical comparison
and the primitive-basis regeneration workaround. The independent check
covers the vibrational metric; broader KEO, potential, and property
validation remains model specific.

Next batches
------------

* Validate a mass-factorised SPAR model, including isotopic changes and
  any mass-dependent coordinate/frame definitions.
* Reconcile and run the published PEF/DMF examples against the current
  constructor and index conventions.
* Package and validate the external H2CS/MOTY kinetic checkpoint linked
  from its ExoMol model page, then extend the molecule catalogue to the
  other first-release families.
* Expand the potential, kinetic, and property selector catalogue with
  parameter definitions, units, aliases, and validated examples.
* Extend convergence, refinement, line-list, and large-calculation workflows.
* Diagnose the separately recorded makefile omissions and primitive-basis
  checkpoint failure before changing program code.

The SiH2/CATS and H2CS entries currently report inspected configurations and
dependencies. They do not claim newly reproduced spectra. External websites
are checked separately from the local Sphinx build, since availability can
change and some publisher endpoints reject automated requests.

XY3 SPAR kinetic operators (October 2026)
-----------------------------------------

The page ``spar_xy3.rst`` documents the ``R-SYMPHI-BETA`` coordinates and the
1D-based and 2D-based external KEOs for PH3. Its code statements were
checked against commit ``a3902fd`` on the ``xy3-exact-2d-keo`` branch
(``mol_xy3.f90``, ``perturbation.f90``, ``fields.f90``, ``molecules.f90``).
The numerical results on that page were obtained on Linux with gfortran
(``-O3 -march=x86-64-v3 -fopenmp``) and OpenBLAS:

* 1D-based KEO, polyad 14, J = 0--4, with commit ``ae08647``. The
  downloadable step-1 input reproduced the J = 0 energies of the original
  run to 1e-6 cm-1 in a clean directory, and reruns of J = 1--4 agreed to
  5e-5 cm-1.
* 2D-based KEO, polyad 14, J = 0 and 1, with commit ``ae08647``. The J = 1
  run used the J = 0-contracted workflow with the table removed after
  step 1.
* 2D-based example of ``inputs/PH3``, polyad 6, J = 0 and 1, with commits
  ``ae08647`` (binary table) and ``a3902fd`` (ASCII table); the energies were
  identical and reproduce the reference values of ``inputs/PH3/README.md``.
* Taylor-mode 2D tables against TROVE's 1D x 1D integration, polyad 6:
  identical energies for 174 J = 0 and 518 J = 1 levels.

The comparison calculation with the built-in linearised KEO was run before
this update and was not repeated. Intensities, D3h(M) calculations, and
other isotopologues are outside the validated scope.
