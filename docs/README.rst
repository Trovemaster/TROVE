Building the manual
===================

From the repository root, create a Python environment and install the
documentation dependencies::

   python -m venv .venv-docs
   # Activate the environment using the command appropriate to your shell.
   python -m pip install -r docs/requirements.txt
   python -m sphinx -b html -W --keep-going docs/source docs/_build/html

Open ``docs/_build/html/index.html``. HTML mathematics uses MathJax, loaded
from the Sphinx default CDN. The browser needs network access for that
script; an offline distribution must provide a local MathJax installation
and set ``mathjax_path`` accordingly. A LaTeX installation is no longer
required just to build HTML.

Keep ``-W --keep-going`` when checking edits. A clean Sphinx build checks
markup, internal references, and local downloads; it does not validate
the scientific content or run TROVE examples. Run the numerical checks
described on each example page separately.

Use ``literalinclude`` for input excerpts, keep full downloadable inputs
under ``docs/source/input``, and record the TROVE revision and validation
scope when updating an example. Avoid silently replacing a published model
with a smaller teaching calculation.
