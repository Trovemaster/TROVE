Getting and building TROVE
==========================

Obtain the source from `Trovemaster/TROVE on GitHub
<https://github.com/Trovemaster/TROVE>`__. Record the branch and commit used
for a calculation: input syntax and checkpoint formats can change.
The SPAR pages are checked against revision ``7831d55`` (7 September 2026).

Build requirements
------------------

The supplied ``makefile`` targets a Unix-like environment with GNU Make,
a Fortran compiler with OpenMP support, a C compiler, and BLAS/LAPACK.
It builds the bundled WIGXJPF library. Intel and GNU Fortran branches are
provided; inspect their flags and library settings on a new system.

The Intel branch defaults to ``ifort`` and MKL. The GNU branch also defaults
to MKL linkage. With a compatible system BLAS/LAPACK installation, the
corresponding configuration is::

   make COMPILER=gfortran LAPACK="-llapack -lblas"

Load the compiler and library environment before building. ``FOR=...`` can
override the compiler executable; verify the flags for that compiler.
The default executable is ``j-trove.x``. ``MODE=debug`` changes the compiler
options. Use separate builds when switching compilers, because Fortran
module and object files are compiler dependent.

.. note::

   At revision ``7831d55``, ``molecules.f90`` imports ``kin_zxy2`` and
   ``kin_abcd``, but the default makefile omits their source files.
   A clean build needs ``kin_zxy2.f90`` and ``kin_abcd.f90`` in its source
   list and module dependencies before compiling ``molecules.f90``.
   The local documentation-validation build included these modules.
   This is a build-configuration issue, separate from the manual edits.

Selecting a user potential
--------------------------

The makefile variable ``pot_user`` selects the source module implementing
the generic user potential/property interface; its default is
``pot_H2O_Conway``. A model using ``GENERAL`` may need another source
module and a rebuild. Record that choice with the executable. A built-in
``POT_TYPE`` is a separate input-level choice; see :doc:`potentials`.

Running a calculation
---------------------

Run each new example in a clean directory, because checkpoint paths are
normally relative to the working directory. A Unix-like shell command is::

   /path/to/j-trove.x < example.inp > example.out

For a Windows executable, use Command Prompt::

   C:\path\to\trove.exe < example.inp > example.out

PowerShell does not support the same ``<`` redirection; use ``cmd /c`` or
a runner that opens the input as standard input. Native Windows compilation
also needs platform-specific compiler and library settings; the Unix
makefile is not a Windows build recipe.

Check the completion message, expected files, and numerical results.
Some Fortran ``STOP`` paths return a successful process status without
completing the calculation. A successful run alone does not establish
basis convergence or validate a model.

Start with :doc:`quickstart` or :doc:`basic_functions_example`. Preserve the
input, build information, and reference results when reporting a problem
on `GitHub <https://github.com/Trovemaster/TROVE/issues>`__.
