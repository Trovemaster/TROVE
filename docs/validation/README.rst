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
