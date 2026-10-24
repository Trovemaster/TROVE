PH3: 1D-based and 2D-based SPAR kinetic energy operators
========================================================

This page describes two kinetic energy operators (KEOs) for pyramidal XY3
molecules. Both are generated outside TROVE, read from ``kinetic.chk``, and
use the SPAR coordinates ``R-SYMPHI-BETA``. No built-in XY3 KEO is involved.
The two operators differ only in how they treat the bending pair
:math:`(s_4,s_5)`:

* **1D-based SPAR KEO.** Every term is a product of one-dimensional basic
  functions, as in Eq. :eq:`spar-sop` of :doc:`spar`. The stretching and
  umbrella factors are exact; the :math:`(s_4,s_5)` dependence is a Taylor
  expansion of order :math:`N`. TROVE integrates all factors with its
  standard 1D machinery.
* **2D-based SPAR KEO.** Every term is a product of 1D functions of
  :math:`r_1,r_2,r_3,\beta` and one two-dimensional function of
  :math:`(s_4,s_5)`. The 2D functions are not expanded in powers. Their
  matrix elements in the harmonic-oscillator basis of modes 4 and 5 are
  computed outside TROVE by 2D quadrature and read from ``keo_2d_s4s5.chk``.

The coordinates and the 2D integration were added on the
``xy3-exact-2d-keo`` branch. This page was checked against commit
``a3902fd`` (23 October 2026). The validation covers energies and
:math:`\mathcal{C}_{3\mathrm{v}}(\mathrm{M})` symmetry labels of
:sup:`31`\ PH\ :sub:`3` up to J = 4 for the 1D-based operator and up to J = 1
for the 2D-based operator. Intensities, the inversion (D3h(M)) problem, and
other isotopologues have not been validated.

Coordinates ``R-SYMPHI-BETA``
-----------------------------

Let :math:`\mathbf e_i` be the unit vector along the X--Y\ :sub:`i` bond and
:math:`\mathbf n` the trisector, the unit vector that makes the same angle
:math:`\beta` with all three bonds:

.. math::

   \mathbf e_1\cdot\mathbf n=\mathbf e_2\cdot\mathbf n
   =\mathbf e_3\cdot\mathbf n=\cos\beta .

Let :math:`\phi_i` be the dihedral angle between the two half-planes that
contain :math:`\mathbf n` and the bonds :math:`j` and :math:`k`
(:math:`\{i,j,k\}=\{1,2,3\}`), so that
:math:`\phi_1+\phi_2+\phi_3=2\pi`. The six vibrational coordinates are

.. list-table:: Coordinates of ``TRANSFORM R-SYMPHI-BETA``
   :class: wrap-table
   :header-rows: 1
   :widths: 12 40 48

   * - Mode
     - Coordinate
     - Note
   * - 1, 2, 3
     - :math:`\xi_i=r_i`
     - Bond lengths in Å
   * - 4
     - :math:`\xi_4=s_4=(2\phi_1-\phi_2-\phi_3)/\sqrt6`
     - First component of the :math:`E` bending pair
   * - 5
     - :math:`\xi_5=s_5=(\phi_2-\phi_3)/\sqrt2`
     - Second component of the :math:`E` bending pair
   * - 6
     - :math:`\xi_6=\beta`
     - Umbrella angle between each bond and :math:`\mathbf n`, in radians

The bond angle :math:`\alpha_i` between the bonds :math:`j` and :math:`k`
follows from :math:`\sin(\alpha_i/2)=\sin\beta\,\sin(\phi_i/2)`. At the
symmetric reference geometry :math:`s_4=s_5=0`, all three
:math:`r_i=r_{\rm e}`, and

.. math::

   \sin\beta_{\rm e}=\frac{2}{\sqrt3}\,\sin\frac{\alpha_{\rm e}}{2}.

The ``EQUILIBRIUM`` block keeps the usual XY3 local values
:math:`r_{\rm e}` and :math:`\alpha_{\rm e}`; TROVE converts them with the
coordinate transformation. For the PH3 inputs,
:math:`r_{\rm e}=1.4118221` Å and :math:`\alpha_{\rm e}=93.3685^\circ`,
which give :math:`\beta_{\rm e}=0.997545` rad (57.1551°).

The coordinates are *not* displacements. The bond lengths and :math:`\beta`
are used as absolute values, so the ``BASIS`` borders are absolute
(for example 0.91--3.51 Å for :math:`r_i` and 0.40--1.50 rad for
:math:`\beta`). Use the ``AUTO`` (or, for a harmonic basis,
``AUTO-HARMON``) kinetic expansion type, which evaluates the basic
functions at the raw coordinate. Do not use ``AUTOMATIC``: it adds the
reference value to the coordinate a second time (see :doc:`spar_input`).

In :math:`\mathcal{C}_{3\mathrm{v}}(\mathrm{M})` the bond lengths are
permuted, :math:`(s_4,s_5)` transform as the :math:`E` representation with
the same matrices as for ``R-SYMPHI-DELTA``, and :math:`\beta` is invariant.
In D3h(M) the operations that contain :math:`E^*` also map :math:`\beta` to
:math:`\pi-\beta`. The transformation and the symmetry operations are
implemented in ``ML_coordinate_transform_XY3`` and
``ML_symmetry_transformation_XY3`` in ``mol_xy3.f90``.

The two external operators below are derived in the frame with the
:math:`z` axis along :math:`\mathbf n`. The :math:`x` axis is fixed by
requiring that the azimuthal angles of the three bonds deviate from 0,
:math:`2\pi/3` and :math:`4\pi/3` by amounts that sum to zero. This frame is
a property of the KEO file; a dipole moment function used with these
operators must be expressed in the same frame.

Structure of the exact operator
-------------------------------

Each element of :math:`G_{\rm vib}`, :math:`G_{\rm rot}`, :math:`G_{\rm cor}`
and the pseudopotential :math:`U` is linear in the inverse masses of the
central atom and of the three equivalent Y atoms. In these coordinates and
this frame every term has the form

.. math::
   :label: xy3-keo-terms

   G(\boldsymbol\xi)=\sum_L \frac{1}{m_{\ell_0(L)}}
   \prod_{i=1}^{3} r_i^{-a_i(L)}\; B_{b(L)}(\beta)\; F_L(s_4,s_5),

with :math:`a_i\in\{0,1,2\}` and

.. math::

   B_b(\beta)\in\{\,1,\ \csc^2\beta,\ \cot\beta,\
   \sin\beta\cos\beta,\ \sin^2\beta\,\}.

The two generators described below stop if any other dependence on
:math:`\beta` appears. Both find these five functions only: the exact
operator used for the 2D-based KEO, and the Taylor operator up to
:math:`N=4`. Only :math:`F_L` couples :math:`s_4` and :math:`s_5`; all other
factors are already one-dimensional basic functions. The two
representations differ in how they treat :math:`F_L`.

Because of the mass factorisation, both ``kinetic.chk`` files use the
``BASIC-FUNCTION MASS`` format (see :doc:`spar_checkpoints`). The mass index
:math:`\ell_0=1` selects the mass of ``ZMAT`` atom 1 (P) and
:math:`\ell_0=2` the mass of atom 2 (H). The three Y masses must be equal.

1D-based SPAR KEO
-----------------

Representation
^^^^^^^^^^^^^^

The 2D factors are replaced by their Taylor expansions,

.. math::

   F_L(s_4,s_5)\approx\sum_{i_4+i_5\le N} c^{(L)}_{i_4 i_5}\,
   s_4^{\,i_4} s_5^{\,i_5},

so that every term is a product of the basic functions below. The example
uses :math:`N=4`.

.. list-table:: Function library of the 1D-based PH3 KEO
   :class: wrap-table
   :header-rows: 1
   :widths: 18 18 64

   * - Mode
     - Index
     - Function
   * - 1, 2, 3
     - 1, 2
     - :math:`r_i^{-1}`, :math:`r_i^{-2}`
   * - 4, 5
     - :math:`n=1,\ldots,N`
     - :math:`s^n`
   * - 6
     - 1, 2, 3, 4
     - :math:`\csc^2\beta`, :math:`\cot\beta`,
       :math:`\sin\beta\cos\beta`, :math:`\sin^2\beta`

In the input this library reads:

.. literalinclude:: input/PH3_SPAR_1D_P14_step1.inp
   :language: text
   :start-at: BASIC-FUNCTION
   :end-at: END

The second number of a factor is a power: ``Csc 2 1.0`` is
:math:`\csc^2\beta`, not :math:`\csc 2\beta`. Index 0 is the implicit
constant.

The operator file
^^^^^^^^^^^^^^^^^

The first records of ``kinetic.chk`` are::

     1  1  1   <-  g_vib
       1    1       1  0    1.00000000000000e+00   1  0  0  0  0  0  0
       1    1       2  0    1.00000000000000e+00   2  0  0  0  0  0  0
       1    2       1  0    1.00000000000000e+00   1  0  0  0  0  0  0
       1    2       2  0   -1.50000000000000e+00   1  0  0  0  0  0  4
       1    2       3  0    6.12372435695794e-01   1  0  0  0  0  1  4

The columns are ``k1 k2 term grid coefficient l0 l1 ... l6``. The first two
records give :math:`G^{\rm vib}_{11}=m_{\rm P}^{-1}+m_{\rm H}^{-1}`; the
next two start :math:`G^{\rm vib}_{12}=m_{\rm P}^{-1}
(1-\tfrac32\sin^2\beta+\ldots)`, which equals
:math:`m_{\rm P}^{-1}\cos\alpha_3` at :math:`s_4=s_5=0`. The file contains
2457 :math:`G_{\rm vib}`, 1271 :math:`G_{\rm rot}`, 1737 :math:`G_{\rm cor}`
and 170 pseudopotential records (5635 in total). The largest index is 2 for
the stretches and 4 for the other modes; :math:`i_4+i_5\le4`.

The generator
^^^^^^^^^^^^^

:download:`generate_xy3_spar_keo.py <input/generate_xy3_spar_keo.py>`
derives the operator with the :math:`t`--:math:`s` matrix formalism. It
uses exact arithmetic in the number field
:math:`\mathbb Q(\sqrt2,\sqrt3)`, so the coefficients carry no rounding
error before they are written. Run it in an otherwise empty directory::

   python3 generate_xy3_spar_keo.py

The script needs Python 3 and NumPy (NumPy is used only for the final
check) and takes about 10 s. It writes its files next to the script:

* ``kinetic.chk``: the mass-factorised operator for ``BASIC-FUNCTION MASS``.
  This is the file used below.
* ``kinetic_nomass.chk``: the same operator with the
  :sup:`31`\ PH\ :sub:`3` masses multiplied in. TROVE's energy-unit
  conversion is *not* applied, so this file cannot be used as a
  ``BASIC-FUNCTION KINETIC`` checkpoint without rescaling.
* ``combination.chk``: the function library in the earlier four-field
  constructor syntax. TROVE does not read it; the ``BASIC-FUNCTION`` block
  of the input is authoritative.

Finally, the script compares the operator with a direct numerical inversion
of the :math:`12\times12` Cartesian :math:`t` matrix, in the generator's
unscaled units (Da, Å). At :math:`s_4=s_5=0` with unequal bond lengths the
largest :math:`G` errors are below :math:`10^{-15}`. At
:math:`(s_4,s_5)=(0.02,-0.015)` they are at most :math:`2.3\times10^{-10}`,
consistent with the fifth-order truncation. The pseudopotential agrees
within the error of its finite-difference reference (a few times
:math:`10^{-6}`).

Basis set and staged run
^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: input/PH3_SPAR_1D_P14_step1.inp
   :language: text
   :start-at: BASIS
   :end-at: END

The coordinates split into three contracted classes: the stretches
(class 1), the bending pair :math:`(s_4,s_5)` (class 2), and :math:`\beta`
(class 3). The primitive basis therefore factorises as 3D + 2D + 1D. Modes
4 and 5 use harmonic-oscillator functions. TROVE integrates the basic
functions numerically on the grid given by ``points`` and ``borders``.
Choose the borders so that the highest harmonic function (here
:math:`v=14`) has decayed at the grid ends. The classical turning point of
:math:`v=14` lies at :math:`|s|\approx1.54`, so borders of :math:`\pm1.5`
are too narrow at polyad 14; the example uses :math:`\pm2.5`.
``reduced 3`` selects the reduced Hamiltonian used to contract the class.
``nocheck`` skips the test that the one-dimensional reference Hamiltonians
of the modes in the class are identical.

The full calculation uses three inputs:

* :download:`PH3_SPAR_1D_P14_step1.inp <input/PH3_SPAR_1D_P14_step1.inp>`
  (J = 0, saves the contracted matrix elements),
* :download:`PH3_SPAR_1D_P14_step2.inp <input/PH3_SPAR_1D_P14_step2.inp>`
  (``model j=0``; converts the matrix elements to the J = 0
  representation), and
* :download:`PH3_SPAR_1D_P14_step3_J1.inp <input/PH3_SPAR_1D_P14_step3_J1.inp>`
  (J = 1; change ``Jrot`` for other J).

Generate ``kinetic.chk`` in the working directory, then run the three
stages, for example::

   python3 generate_xy3_spar_keo.py
   /path/to/j-trove.x < PH3_SPAR_1D_P14_step1.inp    > step1.out
   /path/to/j-trove.x < PH3_SPAR_1D_P14_step2.inp    > step2.out
   /path/to/j-trove.x < PH3_SPAR_1D_P14_step3_J1.inp > step3_J1.out

Step 1 computes the potential expansion from the ``POTEN`` block. Steps 1
and 3 end with ``End of TROVE``, step 2 with ``...done!``. A J > 0 run
prints energies relative to the J = 0 zero-point energy if the file
``j0eigen_descr0_1.chk`` is present. That file is written by a step-3 run
with ``Jrot 0``. Without it, the energies are relative to the minimum of
the potential.

The step-3 input uses the ``syevr`` diagonaliser. For J = 1 its 11 142
eigenvalues were identical to those of ``syevd`` at the printed precision,
and the run was 21 times faster than with ``syev``.

A smaller self-contained example at polyad 6 with Taylor order :math:`N=8`
(17 646 records in ``kinetic.chk``) is provided in ``inputs/PH3/1D_SPAR``
together with ``make_spar_kinetic.py`` and ``run.sh``.

2D-based SPAR KEO
-----------------

Representation
^^^^^^^^^^^^^^

Here the 2D factors are kept exact. They are expanded in one set of 2D
functions :math:`\Phi_k`, shared by all terms of all components:

.. math::

   F_L(s_4,s_5)=F_L(0,0)+\sum_{k=1}^{N_r} c_{Lk}\,\Phi_k(s_4,s_5),
   \qquad \Phi_k(0,0)=0.

Each :math:`\Phi_k` is labelled by a pair :math:`(i_4,i_5)`, in order of
increasing :math:`i_4+i_5`: (1,0), (0,1), (2,0), (1,1), (0,2), and so on.
The label (0,0) denotes the constant. These pairs are *labels*, not powers.
In ``kinetic.chk`` they occupy the columns ``l4`` and ``l5``; the other
columns keep their 1D meaning. Such a ``kinetic.chk`` is meaningful only
together with the ``keo_2d_s4s5.chk`` generated with it.

Because every :math:`\Phi_k` vanishes at :math:`s_4=s_5=0`, as does every
monomial :math:`s_4^{i_4}s_5^{i_5}` with :math:`i_4+i_5>0`, the KEO at the
reference geometry is the same whether a label is read as :math:`\Phi_k` or
as a monomial. TROVE uses this value to estimate the harmonic parameters of
modes 4 and 5. Beyond the reference geometry, the 2D table must be used.

Construction of the 2D functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The generator ``exact_2d_keo.py`` (in ``inputs/PH3`` of the repository;
:download:`download <../../inputs/PH3/exact_2d_keo.py>`) evaluates the exact
operator on a polar grid with :math:`(s_4,s_5)=\rho(\cos\theta,\sin\theta)`.
It uses 80 Gauss--Legendre points in :math:`\rho\in[0,R]` and 96 equally
spaced values of :math:`\theta`. The dependence on :math:`r_i` and
:math:`\beta` stays analytic, as in Eq. :eq:`xy3-keo-terms`. For every
component, mass index, stretching powers and :math:`\beta` function, the
script collects :math:`F_L(s_4,s_5)-F_L(0,0)` on the grid. It then performs
a singular value decomposition of all these functions, weighted by
:math:`\exp(-\rho^2/4f_t^2)` and the quadrature weights. The leading
:math:`N_r` singular vectors define the :math:`\Phi_k`.

For PH3 the functions span a space of numerical rank 63. With
:math:`N_r=63` and :math:`R=2.45`, the 63rd and 64th singular values are
:math:`1.5\times10^{-9}` and :math:`2.3\times10^{-15}`, and the largest
reconstruction error, weighted with the Gaussian envelope, is
:math:`6\times10^{-16}`. In the region sampled by the basis functions the
SVD therefore adds no significant approximation. The 63 labels fit into
:math:`0\le i_4,i_5\le8`, which defines :math:`k_{\max}=8`.

The 2D integral table
^^^^^^^^^^^^^^^^^^^^^

For the harmonic-oscillator functions :math:`\psi_v` of modes 4 and 5, the
generator computes

.. math::

   I^{(o_4o_5)}_{k_4k_5}(v_4v_4';v_5v_5')=\iint
   D^{(o_4)}_{v_4v_4'}(s_4)\,D^{(o_5)}_{v_5v_5'}(s_5)\,
   \Phi_{k_4k_5}(s_4,s_5)\,ds_4\,ds_5,

with the one-dimensional operator densities

.. math::

   D^{(0)}_{vv'}=\psi_v\psi_{v'},\quad
   D^{(1)}_{vv'}=\psi_v\,\psi'_{v'},\quad
   D^{(2)}_{vv'}=-\psi'_v\,\psi_{v'},\quad
   D^{(3)}_{vv'}=-\psi'_v\,\psi'_{v'},

where :math:`\psi'=d\psi/ds`. These follow the derivative placement of
TROVE's 1D KEO matrix elements. For a term of :math:`G^{\rm vib}_{k_1k_2}`,
the label :math:`o_4` is 0 if neither :math:`k_1` nor :math:`k_2` equals 4,
1 if only :math:`k_2=4`, 2 if only :math:`k_1=4`, and 3 if both equal 4.
:math:`o_5` follows the same rule for mode 5. This needs the nine
combinations (0,0), (1,0), (2,0), (0,1), (0,2), (3,0), (0,3), (2,1) and
(1,2). :math:`G_{\rm rot}` and :math:`U` use (0,0). For the symmetrised
Coriolis term of mode 4, TROVE uses the sum of (1,0) and (2,0); for mode 5,
the sum of (0,1) and (0,2).

The table replaces the product of the 1D matrix elements of modes 4 and 5
in the routines that build the vibrational, rotational, Coriolis and
pseudopotential matrix elements. TROVE holds the whole array
:math:`(0{:}3,0{:}3,0{:}k_{\max},0{:}k_{\max},0{:}v_{\max},\ldots)` in
memory: :math:`16(k_{\max}+1)^2(v_{\max}+1)^4` double-precision numbers,
25 MB for :math:`v_{\max}=6` and 525 MB for :math:`v_{\max}=14`. By default
the generator writes only the non-zero elements, as ASCII text
(1 382 976 entries and about 83 MB for :math:`v_{\max}=6`); ``--binary``
writes the whole array as a binary stream file instead. Both formats are
described in :doc:`spar_checkpoints`.

Requirements and checks
^^^^^^^^^^^^^^^^^^^^^^^

There is no input keyword. TROVE reads ``keo_2d_s4s5.chk`` whenever the file
exists in the working directory and prints, for an ASCII table::

   Loaded exact 2D (s4,s5) KEO integrals from keo_2d_s4s5.chk (ASCII): kmax=  8, vmax=  6, entries=     1382976

(without ``(ASCII)`` and the entry count for a binary table). Look for this
line. TROVE stops with a message (``PTread_me2d_s4s5: ...``) if:

* there are fewer than five modes;
* modes 4 and 5 do not form a contracted class of their own;
* :math:`v_{\max}` of the table is smaller than the basis range of mode 4
  or 5;
* :math:`k_{\max}` is smaller than the largest index of modes 4 and 5 in
  :math:`G_{\rm vib}`, :math:`G_{\rm rot}` or :math:`G_{\rm cor}`;
* an ASCII entry lies outside the ranges given in the header, or the number
  of entries disagrees with the header;
* the file cannot be opened or read.

TROVE does not check, and the user must ensure, that:

* modes 4 and 5 use the ``harmonic`` primitive basis with the length
  parameter :math:`f_t` passed to the generator (``--ft``). TROVE prints
  :math:`f_t` as the first number of the ``params`` line of this basis, for
  example ``params = 0.28658334 1145.16544179 0.00000000``;
* the quadrature radius (``--rmax``) agrees with the ``borders`` of modes 4
  and 5;
* ``kinetic.chk`` and ``keo_2d_s4s5.chk`` come from the same generator run;
* the ``BASIC-FUNCTION`` block defines at least :math:`k_{\max}` functions
  for modes 4 and 5. The example lists ``I 1`` to ``I 8`` as placeholders
  so that the labels are valid indices.

The presence of the file alone changes the KEO integration, so remove it
before running unrelated calculations in the same directory. A Fortran
``stop`` returns exit status 0; check the output for the message and for
``End of TROVE``.

The file is needed whenever TROVE computes KEO matrix elements, which in
the examples happens in step 1. Two workflows were tested:

* **Without** ``model j=0`` (the ``inputs/PH3`` example): steps 1 and 2
  keep the classes. Keep the file in the directory; both runs load it.
* **With** ``model j=0`` (the staged workflow of the 1D-based example):
  the conversion and J > 0 runs merge the classes, and
  ``PTread_me2d_s4s5`` stops with "modes 4 and 5 must form a separate
  class". Rename or remove ``keo_2d_s4s5.chk`` after step 1. The contracted
  matrix elements saved in step 1 already contain the 2D integrals.

Running the example
^^^^^^^^^^^^^^^^^^^

The repository directories ``inputs/PH3`` and ``inputs/PH3/2D_exact``
contain the generator, a README, and two inputs at polyad 6:
:download:`PH3_exact2D_P6_J0.inp <../../inputs/PH3/PH3_exact2D_P6_J0.inp>`
and :download:`PH3_exact2D_P6_J1.inp <../../inputs/PH3/PH3_exact2D_P6_J1.inp>`
(``inputs/PH3/2D_exact`` also includes a precomputed ``kinetic.chk``, a
compressed ``keo_2d_s4s5.chk.xz``, and ``run.sh``).
Their basis rows for modes 4 and 5 are

.. literalinclude:: ../../inputs/PH3/PH3_exact2D_P6_J0.inp
   :language: text
   :start-at: 2,'harmonic'
   :end-before: 3,'numerov'

Run, in one directory::

   python3 exact_2d_keo.py --mode exact --vmax 6 --out-dir .
   /path/to/j-trove.x < PH3_exact2D_P6_J0.inp > PH3_exact2D_P6_J0.out
   /path/to/j-trove.x < PH3_exact2D_P6_J1.inp > PH3_exact2D_P6_J1.out

The generator needs NumPy and SciPy and takes a few seconds. Its options
must agree with the TROVE input:

.. list-table:: Options of ``exact_2d_keo.py``
   :class: wrap-table
   :header-rows: 1
   :widths: 22 78

   * - Option
     - Meaning and required value
   * - ``--mode``
     - ``exact`` writes ``kinetic.chk`` and ``keo_2d_s4s5.chk`` for the exact
       operator. ``taylor`` writes only the table for the monomials
       :math:`s_4^{i_4}s_5^{i_5}` with :math:`i_4+i_5\le` ``--kmax``.
   * - ``--vmax``
     - At least the basis range of modes 4 and 5 (``range 0, 6``).
   * - ``--ft``
     - :math:`f_t` from the ``params`` line. The default, 0.2865833376039929,
       belongs to the supplied PH3 inputs.
   * - ``--rmax``
     - Quadrature radius :math:`R`; use the ``borders`` of modes 4 and 5.
   * - ``--kmax``
     - Taylor mode only. The exact mode uses :math:`k_{\max}=8` and
       :math:`N_r=63`, set in the script.
   * - ``--binary``
     - Write the table as a binary stream file instead of ASCII.
   * - ``--out-dir``
     - Output directory.

Keep ``KinOrder 16``, the ``BASIC-FUNCTION`` block, and the three classes of
the ``BASIS`` block as in the example.

Checking the 2D integration
^^^^^^^^^^^^^^^^^^^^^^^^^^^

With ``--mode taylor --kmax N`` and a 1D-based ``kinetic.chk`` of order
:math:`N`, the 2D route must reproduce TROVE's 1D :math:`\times` 1D
integration. For PH3 with :math:`N=4` at polyad 6, the energies of all 174
J = 0 and 518 J = 1 levels were identical at the printed precision
(:math:`10^{-6}\,\mathrm{cm}^{-1}`). This holds with
``reduced 2, POST`` and, for J = 0, also with ``reduced 3, nocheck``.

Validation
----------

The SPAR calculations below used commit ``ae08647`` with a binary 2D table,
built with gfortran (``-O3 -march=x86-64-v3 -fopenmp``). The polyad-6
example was repeated with ``a3902fd`` and an ASCII table; all 174 J = 0 and
518 J = 1 energies were identical. The comparison calculation, labelled
R-ALPHA below, was run earlier with TROVE's built-in linearised KEO
(``COORDS LINEAR``, ``TRANSFORM R-ALPHA``, ``KinOrder 6``), the same
potential, and the same polyad cut-off. Its own KEO is a sixth-order
expansion and it uses a 3D bending basis. The differences from R-ALPHA
therefore compare two approximate routes; they are not errors of either
SPAR KEO.

Band centres at polyad 14
^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: J = 0 energies (cm\ :sup:`-1`) of :sup:`31`\ PH\ :sub:`3`, polyad 14, relative to the zero-point energy (ZPE)
   :class: wrap-table
   :header-rows: 1
   :widths: 16 8 19 19 19 19

   * - Band
     - Γ
     - 1D-based (Taylor, N = 4)
     - 2D-based (exact)
     - 1D − 2D
     - R-ALPHA
   * - ZPE
     - A1
     - 5214.0561
     - 5214.0559
     - 0.0002
     - 5214.0615
   * - ν2
     - A1
     - 992.1033
     - 992.1060
     - −0.0027
     - 992.1365
   * - ν4
     - E
     - 1118.3008
     - 1118.3066
     - −0.0058
     - 1118.3158
   * - 2ν2
     - A1
     - 1972.3979
     - 1972.4012
     - −0.0032
     - 1972.5815
   * - ν2 + ν4
     - E
     - 2108.0754
     - 2108.0933
     - −0.0180
     - 2108.1403
   * - 2ν4
     - A1
     - 2226.6854
     - 2226.7355
     - −0.0501
     - 2226.8041
   * - 2ν4
     - E
     - 2234.8512
     - 2234.8710
     - −0.0198
     - 2234.9097
   * - ν1
     - A1
     - 2321.1075
     - 2321.1091
     - −0.0017
     - 2321.1186
   * - ν3
     - E
     - 2326.8547
     - 2326.8553
     - −0.0007
     - 2326.8681
   * - 3ν2
     - A1
     - 2939.9619
     - 2939.9645
     - −0.0026
     - 2940.8581
   * - 3ν4
     - E
     - 3332.9974
     - 3333.1397
     - −0.1423
     - 3333.4300

The three calculations have the same number of J = 0 levels (5151).
Differences are computed from the unrounded energies. The Taylor
truncation in :math:`(s_4,s_5)` matters mainly for bending overtones. Over
all levels with the same symmetry and index, the largest difference between
the 1D- and 2D-based results is 0.14 cm\ :sup:`-1` below
4000 cm\ :sup:`-1`, 0.66 cm\ :sup:`-1` below 6000 cm\ :sup:`-1`, and
2.2 cm\ :sup:`-1` below 8000 cm\ :sup:`-1`.

Rotational levels
^^^^^^^^^^^^^^^^^

For J = 1--4 the 1D-based KEO and R-ALPHA give the same number of levels of
each total symmetry below 5000 cm\ :sup:`-1`. Taking the levels of each
symmetry in energy order, the rotational labels
:math:`(\Gamma_{\rm rot};J,K,\tau)` agree for every level below
3000 cm\ :sup:`-1` and for all but one below 4000 cm\ :sup:`-1`; at higher
energies, nearly degenerate levels are ordered differently in the two
calculations. The R-ALPHA J = 4 output does not contain the E block.

.. list-table:: Ground vibrational state, energies in cm\ :sup:`-1` relative to the ZPE
   :class: wrap-table
   :header-rows: 1
   :widths: 6 8 22 22 22 20

   * - J
     - Γ
     - :math:`(\Gamma_{\rm rot};J,K,\tau)`
     - 1D-based
     - R-ALPHA
     - 2D-based
   * - 1
     - E
     - (E; 1, 1, 0)
     - 8.371093
     - 8.371072
     - 8.371080
   * - 1
     - A2
     - (A2; 1, 0, 1)
     - 8.904314
     - 8.904217
     - 8.904203
   * - 2
     - E
     - (E; 2, 2, 0)
     - 24.577967
     - 24.577979
     -
   * - 2
     - E
     - (E; 2, 1, 0)
     - 26.177272
     - 26.177051
     -
   * - 2
     - A1
     - (A1; 2, 0, 0)
     - 26.709810
     - 26.709512
     -
   * - 3
     - A2
     - (A2; 3, 3, 1)
     - 48.616459
     - 48.616555
     -
   * - 3
     - A1
     - (A1; 3, 3, 0)
     - 48.616462
     - 48.616557
     -
   * - 3
     - E
     - (E; 3, 2, 0)
     - 51.282479
     - 51.282182
     -
   * - 3
     - E
     - (E; 3, 1, 0)
     - 52.878714
     - 52.878181
     -
   * - 3
     - A2
     - (A2; 3, 0, 1)
     - 53.410231
     - 53.409619
     -
   * - 4
     - A1
     - (A1; 4, 3, 0)
     - 84.214684
     - 84.214352
     -
   * - 4
     - A2
     - (A2; 4, 3, 1)
     - 84.214699
     - 84.214367
     -
   * - 4
     - A1
     - (A1; 4, 0, 0)
     - 88.996217
     - 88.995160
     -

The tiny :math:`K=3` splittings, 2--3\ :math:`\times10^{-6}`
cm\ :sup:`-1` at J = 3 and :math:`1.5\times10^{-5}` cm\ :sup:`-1` at J = 4,
agree with R-ALPHA at the printed precision. For the 2D-based KEO at
J = 1, the lowest :math:`\nu_2` and :math:`\nu_4` rovibrational levels
differ from the 1D-based ones by 0.0003 to 0.012 cm\ :sup:`-1`.

Polyad-6 example
^^^^^^^^^^^^^^^^

The ``inputs/PH3/2D_exact`` and ``inputs/PH3/1D_SPAR`` examples give the
following values (cm\ :sup:`-1`, relative to the respective J = 0 ZPE:
5214.173894 cm\ :sup:`-1` for the 2D-based KEO and 5214.173918 cm\ :sup:`-1`
for the 1D-based KEO at :math:`N=8`). The fifth and sixth columns show the
1D-based operator at Taylor orders :math:`N=8` (from ``inputs/PH3/1D_SPAR``)
and :math:`N=4` in the same basis.

.. list-table:: Polyad 6, ``reduced 2, POST``
   :class: wrap-table
   :header-rows: 1
   :widths: 6 8 24 16 16 16 14

   * - J
     - Γ
     - State
     - 2D-based (exact)
     - 1D-based (N = 8)
     - 1D-based (N = 4)
     - 1D(4) − 2D
   * - 0
     - A1
     - ν2
     - 992.4200
     - 992.4200
     - 992.4172
     - −0.0027
   * - 0
     - E
     - ν4
     - 1118.4931
     - 1118.4932
     - 1118.4864
     - −0.0067
   * - 0
     - A1
     - 2ν2
     - 1973.8687
     - 1973.8687
     - 1973.8649
     - −0.0038
   * - 0
     - E
     - ν2 + ν4
     - 2109.0113
     - 2109.0115
     - 2108.9914
     - −0.0200
   * - 0
     - A1
     - 2ν4
     - 2227.2674
     - 2227.2677
     - 2227.2149
     - −0.0524
   * - 0
     - E
     - 2ν4
     - 2235.4937
     - 2235.4942
     - 2235.4675
     - −0.0262
   * - 0
     - A1
     - ν1
     - 2321.6121
     - 2321.6121
     - 2321.6105
     - −0.0016
   * - 0
     - E
     - ν3
     - 2327.2029
     - 2327.2029
     - 2327.2022
     - −0.0007
   * - 1
     - E
     - ground state, K = 1
     - 8.3712
     - 8.3712
     - 8.3712
     - 0.0000
   * - 1
     - A2
     - ground state, K = 0
     - 8.9045
     - 8.9045
     - 8.9046
     - 0.0001
   * - 1
     - E
     - ν2, K = 1
     - 1000.6400
     - 1000.6400
     - 1000.6377
     - −0.0023
   * - 1
     - A2
     - ν2, K = 0
     - 1000.9660
     - 1000.9661
     - 1000.9642
     - −0.0018
   * - 1
     - E
     - ν4
     - 1123.3475
     - 1123.3477
     - 1123.3464
     - −0.0011
   * - 1
     - E
     - ν4
     - 1127.6001
     - 1127.6002
     - 1127.5935
     - −0.0066
   * - 1
     - A1
     - ν4
     - 1130.3632
     - 1130.3633
     - 1130.3513
     - −0.0120
   * - 1
     - A2
     - ν4
     - 1130.7113
     - 1130.7114
     - 1130.6983
     - −0.0130

Reproducibility checks
^^^^^^^^^^^^^^^^^^^^^^

* The downloadable step-1 input, run in a clean directory with only
  ``kinetic.chk``, reproduces the J = 0 energies of the original run to
  :math:`10^{-6}` cm\ :sup:`-1`. Reruns of the J = 1--4 stages agree with
  the original runs to :math:`5\times10^{-5}` cm\ :sup:`-1`, with identical
  rotational labels.
* ``generate_xy3_spar_keo.py`` reproduces the ``kinetic.chk`` used here byte
  for byte. ``exact_2d_keo.py`` is not bit-reproducible: rounding in the SVD
  can move a few coefficients across the :math:`10^{-14}` cut-off, so the
  number of records varies by a few (27 846 to 27 848 for PH3). This had no
  effect on the printed energies.
* The vibrational labels :math:`(v_4,v_5)` of the 2D harmonic class are
  products of 1D harmonic-oscillator labels. For degenerate overtones the
  printed assignment can switch between equivalent labels, for example
  (2, 0) and (0, 2), between builds or between the 1D and 2D routes. The
  energies and the symmetry labels do not change. Compare energies and
  symmetries, not :math:`(v_4,v_5)`.

Resources at polyad 14
^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Step 1 (J = 0), 16 OpenMP threads
   :class: wrap-table
   :header-rows: 1
   :widths: 34 22 22 22

   * - Quantity
     - R-ALPHA
     - 1D-based
     - 2D-based
   * - Primitive class structure
     - 3D + 3D
     - 3D + 2D + 1D
     - 3D + 2D + 1D
   * - Potential expansion terms
     - 3002
     - 1662
     - 1662
   * - KEO records in ``kinetic.chk``
     - (built in)
     - 5635
     - 27 846
   * - Peak memory (TROVE's ``Maximal memory``)
     - 11.9 GB
     - 1.3 GB
     - 1.9 GB
   * - Wall-clock time (TROVE timer)
     - 3286 s
     - 78 s
     - 157 s

The R-ALPHA run used a different build and BLAS library, and the
wall-clock times depend on the BLAS library. The 1D- and 2D-based runs used
multithreaded OpenBLAS; with the reference (single-threaded) BLAS, step 1
of the 1D-based calculation took 426 s. The memory reduction and the
smaller potential expansion are consequences of the 3D + 2D + 1D class
structure. The 2D-based run also holds the 525 MB table.

Choosing a representation
-------------------------

* The 1D-based KEO uses only the standard SPAR route and works with any
  basis for modes 4 and 5. Its error comes from the Taylor order in
  :math:`(s_4,s_5)`. At :math:`N=4` this is a few
  :math:`10^{-3}` cm\ :sup:`-1` for the fundamentals and grows for the
  bending overtones.
* The 2D-based KEO has no expansion error in :math:`(s_4,s_5)` beyond the
  quadrature and the (here negligible) SVD truncation. It requires a
  harmonic basis for modes 4 and 5 with a matching :math:`f_t`, a separate
  class for them, and a table that grows as :math:`(v_{\max}+1)^4`. The
  implementation is specific to modes 4 and 5 and to the file name
  ``keo_2d_s4s5.chk``.
* Both operators are exact in :math:`r_i` and :math:`\beta`. Both apply,
  through the mass index, to XY3 isotopologues with three equivalent Y
  nuclei, but only :sup:`31`\ PH\ :sub:`3` has been tested.
