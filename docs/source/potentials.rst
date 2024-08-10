Potential energy functions
**************************

TROVE provides a larger number of potential energy functions (PEFs) for different molecules already implemented. Most of these PEFs are in modules ``pot_*`` contained in  file ``pot_*.f90``.

 - ``pot_xy2.f90``
 - ``pot_xy3.f90``
 - ``pot_zxy2.f90``
 - ``pot_abcd.f90``
 - ``pot_xy4.f90``
 - ``pot_zxy3.f90``
 - ........

These are a part of the standard TROVE compilation set. Alternatively, a user-defined PEF can be included into the TROVE compilation as a generic 'user-defined' module ``pot_user``, see details below.





Potential Block
===============

The ``Potential`` (``Poten``) block used to specify a PEF, has the following generic structure (using XY\ :sub:`2` as an example)

::

      POTEN
      NPARAM  102
      POT_TYPE  POTEN_XY2_MORSE_COS
      COEFF  list  (powers or list)
       RE13         1  0.15144017558000E+01
       ALPHAE       1  0.92005073880000E+02
       AA           1  0.12705074620000E+01
       B1           1  0.50000000000000E+06
       B2           1  0.50000000000000E+05
       G1           1  0.15000000000000E+02
       G2           1  0.10000000000000E+02
       V0           0  0.00000000000000E+00
       F_0_0_1      1 -0.11243403302598E+02
       F_1_0_0      1 -0.94842865087918E+01
       F_0_0_2      1  0.17366522840412E+05
       F_1_0_1      1 -0.25278354456474E+04
       F_1_1_0      1  0.20295521820240E+03
       F_2_0_0      1  0.38448640879698E+05
       F_0_0_3      1  0.27058767090918E+04
       F_1_0_2      0 -0.47718397149800E+04
      ...
      end

For an example, see `SiH2_XY2_MORSE_COS_step1.inp <https://raw.githubusercontent.com/Trovemaster/TROVE/develop/docs/source/input/SiH2_XY2_MORSE_COS_step1.inp>`_  where this PES is used.  For an example, see `SiH2_XY2_MORSE_COS_step1.inp </input/SiH2_XY2_MORSE_COS_step1.inp>`_

Here ``NPARAM`` is used to specify the number of parameters used to define the PES. ``POT_TYPE`` is the name of the potential energy surface being used which is defined in the ``pot_*.f90 file``. The keywords ``COEFF`` indicates if the potential contains a list of parameter values (``LIST``) or  values with the corresponding expansion powers (``POWERS``), e.g. (for H\ :sub:`2`\ CO):
::

     POTEN
     NPARAM  114
     POT_TYPE  poten_zxy2_mep_r_alpha_rho_powers
     COEFF  powers
     f_0_0      0 0 0 0 0 0   1  0.00000000000000E+00
     f_0_1      1 0 0 0 0 0   1  0.00000000000000E+00
     f_0_1      0 1 0 0 0 0   1  0.00000000000000E+00
     f_0_2      0 0 0 1 0 0   1  0.00000000000000E+00
     f_1_1      0 0 0 0 0 1   1  0.13239727881219E+05
     f_2_1      0 0 0 0 0 2   1  0.46279621687684E+04
     f_3_1      0 0 0 0 0 3   1  0.14394787420943E+04
     f_4_1      0 0 0 0 0 4   1  0.10067584554678E+04
     f_5_1      0 0 0 2 0 0   1  0.31402651686299E+05
     .....
     end

which is from the variational calculations of H\ :sub:`2`\ CO, see the TROVE input `1H2-12C-16O__AYTY__TROVE.inp <https://raw.githubusercontent.com/Trovemaster/TROVE/develop/docs/source/input/1H2-12C-16O__AYTY__TROVE.inp>`_.


The potential parameters are listed after the keyword ``COEFF`` and terminated with the keyword ``END`` with exactly ``NPARAM``. For the ``COEFF  list`` option, the meaning of the columns is as follows:

  +---------+-----+-----------------------+
  | Label   |Index| Value                 |
  +=========+=====+=======================+
  |   b1    |   0 |  0.80000000000000E+06 |
  +---------+-----+-----------------------+
  |   b2    |   0 |  0.80000000000000E+05 |
  +---------+-----+-----------------------+
  |   g1    |   0 |  0.13000000000000E+02 |
  +---------+-----+-----------------------+
  |   g2    |   0 |  0.55000000000000E+01 |
  +---------+-----+-----------------------+
  |   f000  |   0 |  0.00000000000000E+00 |
  +---------+-----+-----------------------+
  |   f001  |   1 |  0.25298724728304E+01 |
  +---------+-----+-----------------------+
  |   f100  |   1 |  0.76001446034650E+01 |
  +---------+-----+-----------------------+


Here 'Labels' are the parameter names,  used only for printing purposes and not in any calculations. The 'Index' field can be used to as a switch to indicate if the corresponding parameter was fitted or can be fitted. Otherwise it has no impact on any evaluations of the PEF values. 'Values' are the actual potential parameters, listed in the order implemented in the corresponding PEF ``POT_TYPE``, for example

.. math::

   \begin{split}
   V(r_1,r_2,\alpha) &= f_{000} + f_{001} y_3 + f_{100} [ y_1 + y_2 ] + f_{100} [ y_1 + y_2 ] + f_{002} y_3^2 + \ldots +  \\
                     & + b_1 e^{-g_1 r_{\rm HH}} + b_2 e^{-g_2 r_{\rm HH}^2} \\
   \end{split}



where

.. math::

   \begin{split}
      y_1 & = 1-e^{-a (r_1 - r_{\rm e})}, \\
      y_2 & = 1-e^{-a (r_2 - r_{\rm e})}, \\
      y_3 &= \cos\alpha-\cos\alpha_{\rm e}, \\
      r_{HH}=\sqrt{r_1^2+r_2^2-2 r_1 r_2 \cos\alpha}.\\
   \end{split}



For the ``COEFF  Powers`` option, the meaning of the columns is as follows:

   +---------+-----+--+--+--+--+--+-----+-----------------------+
   | Label   |   n1|n2|n3|n4|n5|n6|Index| Value                 |
   +---------+-----+--+--+--+--+--+-----+-----------------------+
   |  f000000|    0| 0| 0| 0| 0| 0|    1|  0.00000000000000E+00 |
   +---------+-----+--+--+--+--+--+-----+-----------------------+
   |  f100000|    1| 0| 0| 0| 0| 0|    1|  0.00000000000000E+00 |
   +---------+-----+--+--+--+--+--+-----+-----------------------+
   |  f010000|    0| 1| 0| 0| 0| 0|    1|  0.00000000000000E+00 |
   +---------+-----+--+--+--+--+--+-----+-----------------------+
   |  f000100|    0| 0| 0| 1| 0| 0|    1|  0.00000000000000E+00 |
   +---------+-----+--+--+--+--+--+-----+-----------------------+
   |  f000001|    0| 0| 0| 0| 0| 1|    1|  0.13239727881219E+05 |
   +---------+-----+--+--+--+--+--+-----+-----------------------+
   |  f000002|    0| 0| 0| 0| 0| 2|    1|  0.46279621687684E+04 |
   +---------+-----+--+--+--+--+--+-----+-----------------------+
   |  f000003|    0| 0| 0| 0| 0| 3|    1|  0.14394787420943E+04 |
   +---------+-----+--+--+--+--+--+-----+-----------------------+

where

 - 'Labels' are the parameter name,  for printing purposes only;
 - 'n1', 'n2', 'n3', ... are the 'powers' of an expansion term, e.g.
   :math:`V(r_1,r_2,r_3,r_4,r_5, r_6) = \sum_{n_1,n_2,n_3,n_4,n_5,n_1} f_{n_1,n_2,n_3,n_4,n_5,n_1} \xi_1^{n_1} \xi_2^{n_2} \xi_3^{n_3} \xi_4^{n_4} \xi_5^{n_5} \xi_6^{n_6}`
  - 'Index' is a switch to indicate if the corresponding parameter was fitted or can be fitted, with no impact on any evaluations of the PEF values.
  - 'Values' are the actual potential parameters. Their order is not important for this implementation as long as the corresponding powers are defined.



In case the definition of PEF requires also structural parameters, such as equilibrium bond lengths :math:`r_{\rm e}`\ , equilibrium inter-bond angles :math:`\alpha_{\rm e}`, Morse exponents :math:`a` etc., in the ``COEFF  Powers`` form these parameters should be listed exactly in the order expected by the  implemented of the PEF (similar to the ``COEFF LIST`` form), but with dummy "powers" columns so that their 'values' appear in the right column. For example:
::

    POTEN
    NPARAM  58
    POT_TYPE  poten_C3_R_theta
    COEFF  powers  (powers or list)
    RE12          0      0      0      0     1.29397
    theta0        0      0      0      0     0.000000000000E+00
    f200          2      0      0      0        0.33240693
    f300          3      0      0      0       -0.35060064
    f400          4      0      0      0        0.22690209
    f500          5      0      0      0       -0.11822982
    .....


Here, ``RE12`` and ``theta0`` are two the equilibrium values and the three columns with ``0 0 0`` are given in order to parse their values using column 6.



TROVE configuration option used for PEF and KEO expansions
===========================================================

In the case of non-rigid reference configurations with periodicity used for representing the KEO and PEF as expansions, it is not uncommon that the numerically generated expansion parameter :math:`f_{ijk}(\rho_i)` on a grid of :math:`\rho_i` have discontinuities. They appear because of numerical problems at points like 180, 360, 720 :math:`^\circ`, i.e where trigonometric functions can change their values because of their periodicity constraints. In principle, these problems indicate that there are implementation in the code that could not be resolved. As ad hoc solution, TROVE can check and iron-out possible discontinuities by interpolating between points were the behaviour is nice and smooth. TROVE will always check and report if there are any discontinuity issues at the expansion stage, e.g.
:: 
   check_field_smoothness: potencheck_field_smoothness: potencheck_field_smoothness: poten; an outlier found for iterm =      18 at i =   400: -0.46876607222E+07 vs  0.28075523334E+05
   
In order to switch on the "ironing-out" procedure, the ``iron-out`` card must be added anywhere in the main body of the TROVE input:
::
   IRON-OUT
   




Implemented PEFs
================


XY\ :sub:`2` type
-----------------

There are several PEFs available for this molecule type.


``POTEN_XY2_MORBID``
^^^^^^^^^^^^^^^^^^^^

This form is given by

.. math::

   \begin{split}
   V(r_1,r_2,\alpha) &= f_{000} + f_{001} y_3 + f_{002} y_3^2 +  f_{003} y_4 + \ldots  \\
                   & + (f_{100} + f_{101} y_3 + f_{102} y_3^2 +  f_{103} y_4 + \ldots) [y_1 + y_2] \\
                   & + (f_{200} + f_{201} y_3 + f_{202} y_3^2 +  f_{203} y_4 + \ldots) [y_1^2 + y_2^2] \\
                   & + (f_{110} + f_{111} y_3 + f_{112} y_3^2 +  f_{113} y_4 + \ldots) y_1y_2 \\
                   & + \ldots \\
                     & + b_1 e^{-g_1 r_{\rm HH}} + b_2 e^{-g_2 r_{\rm HH}^2} \\
   \end{split}


where

.. math::

   \begin{split}
      y_1 & = 1-e^{-a (r_1 - r_{\rm e})}, \\
      y_2 & = 1-e^{-a (r_2 - r_{\rm e})}, \\
      y_3 &= \cos\alpha-\cos\alpha_{\rm e}, \\
      r_{HH}=\sqrt{r_1^2+r_2^2-2 r_1 r_2 \cos\alpha}.\\
   \end{split}


For an example, see `SO2_morbid_step1.inp <https://raw.githubusercontent.com/Trovemaster/TROVE/develop/docs/source/input/SO2_morbid_step1.inp>`_  where this PES is used.



``POTEN_XY2_MORSE_COS``
^^^^^^^^^^^^^^^^^^^^^^^

For description and example see 'Potential Block' above with the input file example for SiH\ :sub:`2`  `SiH2_XY2_MORSE_COS_step1.inp <https://raw.githubusercontent.com/Trovemaster/TROVE/develop/docs/source/input/SiH2_XY2_MORSE_COS_step1.inp>`_. This PES was used in [21ClYu]_ .



``POTEN_XY2_TYUTEREV``
^^^^^^^^^^^^^^^^^^^^^^

``POTEN_XY2_TYUTEREV`` is the essentially the same PES as ``POTEN_XY2_MORSE_COS``, with the difference that the ``POTEN`` input part does not contain the structural parameters:


::

      POTENTIAL
      NPARAM  99
      POT_TYPE  poten_xy2_tyuterev
      COEFF  list
      b1        0    0.80000000000000E+06
      b2        0    0.80000000000000E+05
      g1        0    0.13000000000000E+02
      g2        0    0.55000000000000E+01
      f000      0    0.00000000000000E+00
      f001      1    0.25298724728304E+01
      f100      1    0.76001446034650E+01
      ...
      end


Using the structural parameters in the ``POTEN`` section is important for PES refinements, see the corresponding section for details.
The input file example is  `h2s_step1.inp <https://raw.githubusercontent.com/Trovemaster/TROVE/develop/docs/source/input/h2s_step1.inp>`_ .




``poten_co2_ames1``
^^^^^^^^^^^^^^^^^^^


An empirical PES of CO\ :sub:`2` is from [17HuScFr]_. The input file example is  `CO2_bisect_xyz_step1.inp <https://raw.githubusercontent.com/Trovemaster/TROVE/develop/docs/source/input/CO2_bisect_xyz_step1.inp>`_  where this PES is used.



It is programmed using the ``powers`` format as follows:
::


     POTEN
     NPARAM  234
     POT_TYPE  poten_co2_ames1
     COEFF  powers  (powers or list)
     r12ref    0 0 0  1    1.16139973893
     alpha2    0 0 0  1    1.0000000000
     De1       0 0 0  1    155000.0000
     De2       0 0 0  1    40000.0000
     Ae1       0 0 0  1   40000.0000
     Ae2       0 0 0  1   20000.0000
     edamp2    0 0 0  1   -2.0000
     edamp4    0 0 0  1   -4.0000
     edamp5    0 0 0  1   -0.2500
     edamp6    0 0 0  1   -0.5000
     Emin      0 0 0  1  -0.000872131085d0
     Rmin      0 0 0  1  0.1161287540428520D+01
     rminbohr  0 0 0  1  0.2194515245360671D+01
     alpha     0 0 0  1  1.00000000000000000000
     rref      0 0 0  1  1.1613997389
     f000    0    0    0  0    0.3219238090183E+01
     f001    0    0    1  0   -0.2181831901567E+02
     f002    0    0    2  0    0.7355163655139E+02
     f003    0    0    3  0   -0.1531231748456E+03
     f004    0    0    4  0    0.2090079612238E+03
     f005    0    0    5  0   -0.1883325770080E+03
     ......
     end


The first part contains some structural parameters with 'powers' indexes filled with dummy zeros to maintain the ``powers`` format.







**Examples** and **References**




XYZ type
--------


``POTEN_XYZ_KOPUT``
^^^^^^^^^^^^^^^^^^^


The PEF is given by (see [22OwMiYu]_)

.. math::

       V =  \sum_{ijk} f_{ijk} \xi_1^{i} \xi_2^{j} \xi_3^{k},

The vibrational coordinates are

.. math::

   \begin{split}
      \xi_1 &= (r_1-r_1^{\rm eq})/r_1, \\
      \xi_2 &= (r_2-r_2^{\rm eq})/r_2, \\
     \xi_3 &= \alpha-\alpha_{\rm eq},
   \end{split}

where :math:`r_1  = r_1^{\rm eq}`  and :math:`r_2  = r_2^{\rm eq}` are the internal stretching coordinates  and  :math:`\alpha` is the interbond angle, and the equilibrium parameters are :math:`r_1^{\rm eq}`, :math:`r_2^{\rm eq}` and :math:`\alpha_{\rm eq}`. Note that the exponent :math:`k` associated with the bending coordinate :math:`\xi_3` assumes only even values because of the symmetry of the XYZ molecule.


The input file example is  `CaOH_Koput_step1.inp <https://raw.githubusercontent.com/Trovemaster/TROVE/develop/docs/source/input/CaOH_Koput_step1.inp>`_  where this PES is used.



ZXY\ :sub:`2` type
------------------


``poten_zxy2_morse_cos``
^^^^^^^^^^^^^^^^^^^^^^^^

The PEF is expressed analytically as  (see [15AlOvPo]_, [23MeOwTe]_)

.. math::

     V =  \sum_{ijklmn} a_{ijklmn} \xi_{1}^{i} \xi_{2}^{j} \xi_{3}^{k} \xi_{4}^{l} \xi_{5}^{m} \xi_6^{n},


The vibrational coordinates are

.. math::

   \begin{split}
     \xi_i &= 1 - \exp[-b_i(r_i-r_i^{\rm e})],\quad i={\rm 0, 1, 2},  \\
     \xi_4 &= \alpha_1-\alpha_{\rm e}, \\
     \xi_5 &= \alpha_2-\alpha_{\rm e}, \\
     \xi_6 &= 1+\cos\tau.
   \end{split}


where :math:`a_{ijklmn}` are the expansion parameters. Here, :math:`r_0`, :math:`r_1` and :math:`r_2` are the bond lengths, :math:`\alpha_{1}` and :math:`\alpha_{2}` are the bond angles, :math:`\tau` is the dihedral angle between two bond planes, and :math:`b_i` is a Morse-oscillator parameter.

The expansion has a symmetry adapted form, which is invariant upon interchanging the two equivalent atoms Y\ :sub:`1` and Y\ :sub:`2``, i.e. for :math:`r_1 \leftrightarrow r_2` simultaneously with  :math:`\alpha_1 \leftrightarrow \alpha_2`.

The TROVE function for this PEF is ``MLpoten_zxy2_morse_cos``, which can be found in the module ``mol_zxy2.f09``.


``poten_zxy2_morse_cos`` uses the ``powers`` form for the potential parameters. An input file example for H\ :sub:`2`\ CS is  `H2CS_zxy2_morse_cos_step1.inp <https://raw.githubusercontent.com/Trovemaster/TROVE/develop/docs/source/input/H2CS_zxy2_morse_cos_step1.inp>`_  where this PES is used.



XY\ :sub:`3` type (pyramidal)
-----------------------------


``POTEN_XY3_MORBID_10``
^^^^^^^^^^^^^^^^^^^^^^^

This PEF is for ammonia (non-rigid) molecules with an umbrella coordinate representing the 1D inversion motion and other 5 are displacements from their equilibrium values, but it can be equally used for rigid molecules like phosphine as well. As all PEFs in TROVE, ``POTEN_XY3_MORBID_10``  is a symmetry adapted to be fully symmetric for all operations (permutations) of the D\ :sub: `3h`\ (M) (and also C\ :sub:`3v`\ (M)) molecular group symmetry. It uses the following coordinates


.. math::

   \begin{split}
     \xi_i &= 1 - \exp[-a(r_i-r_i^{\rm e})],\quad i={\rm 1, 2, 3},  \\
     \xi_4 &= \frac{\sqrt{2}}{\sqrt{3}} (2\alpha_1-\alpha_2-\alpha_3), \\
     \xi_5 &= \frac{\sqrt{2}}{2} (\alpha_2-\alpha_3), \\
     \xi_6 &= \sin\rho_{\rm e}-\sin\rho,
   \end{split}

where :math:`r_i` (\ :math:`i=1,2,3,`\ ) are three bond lengths, :math:`\alpha_i` (\ :math:`i=1,2,3,`\ ) are three bond angles with :math:`\alpha_i` opposite to :math:`r_i` and :math:`\rho` is an umbrella coordinate defined as an 'average' angle between three bonds and an average symmetry axis as follows:

.. math::

   \sin\rho = \frac{2}{\sqrt{3}}\sin\left(\frac{\bar\alpha}{2}\right)

and

.. math::

     \bar\alpha = \frac{1}{3} (\alpha_1+\alpha_2+\alpha_3).


The corresponding Fortran function is ``MLpoten_xy3_morbid_10``, which can found in ``mol_xy3.f90``.  It uses the ``Coeffs LIST`` form. A TROVE input example  for NH\ :sub:`3` is  `NH3_BYTe_step1.inp <https://raw.githubusercontent.com/Trovemaster/TROVE/develop/docs/source/input/NH3_BYTe_step1.inp>`_ , see [09YuBaYa]_ where this empirical PES was used to compute the BYTe line list for ammonia.


``POTEN_XY3_MORBID_11``
^^^^^^^^^^^^^^^^^^^^^^^

This is a similar PES to ``POTEN_XY3_MORBID_10``, but with the structural parameters (\ :math:`r_{\rm e}`\ , :math:`\alpha_{\rm e}` and :math:`a` (Morse parameter)) included into the POTENTIAL block:
::

    POTEN
    NPARAM   307
    POT_TYPE  poten_xy3_morbid_11
    COEFF  list
     re           0  0.97580369000000E+00
     alphae       0  0.11195131680000E+03
     beta         0  0.21500000000000E+01
     VE           0  0.00000000000000E+00
     FA1          1  0.19897404278336E+02
     FA2          1  0.34574086467373E+06
     FA3          1 -0.45576605893993E+06
     FA4          1  0.21457903182324E+07
     .....
    end


The Fortran function is ``MLpoten_xy3_morbid_11``, which can found in ``mol_xy3.f90``.  A TROVE input example  for H\ :sub:`3`\ O\ :sup:`+` is  `1H3-16O_p__eXeL__model-TROVE.inp <https://raw.githubusercontent.com/Trovemaster/TROVE/develop/docs/source/input/1H3-16O_p__eXeL__model-TROVE.inp>`_ from [20YuTeMi]_ , where it was used to compute an ExoMol line list for this molecule.




A chain molecule of HOOH type
-----------------------------



``POTEN_H2O2_KOPUT_UNIQUE``
^^^^^^^^^^^^^^^^^^^^^^^^^^^


The HOOH PES is given by the expansion

.. math::

     V = \sum_{i_1,i_2,i_3,i_4,i_5,i_6} f_{i_1,i_2,i_3,i_4,i_5,i_6} \xi_1^{i_1} \xi_2^{i_2}  \xi_3^{i_3} \xi_4^{i_4} \xi_5^{i_5} \cos(i_6 \delta)


in terms of the following coordinates:

.. math::

   \begin{split}
     \xi_1 &= \frac{R - R_{\rm e}}{R},  \\
     \xi_2 &= \frac{r_1 - r_{\rm e}}{r_1},  \\
     \xi_3 &= \frac{r_2 - r_{\rm e}}{r_2},  \\
     \xi_4 &= \alpha_1-\alpha_{\rm e}, \\
     \xi_5 &= \alpha_2-\alpha_{\rm e}, \\
     \xi_6 &= \delta,
   \end{split}

where :math:`R, r_1, r_2`  are three bond lengths, :math:`\alpha_i` (\ :math:`i=1,2,`\ ) are two bond angles with :math:`\alpha_i` and :math:`\delta` is a dihedral angle.


This a ``powers`` type:
::

    POTEN
    NPARAM  166
    POT_TYPE  POTEN_H2O2_KOPUT_UNIQUE
    COEFF  powers  (powers or list)  (oo oh1  oh2 h1oo h2oo rho)
    f000000       0  0  0  0  0  0  1  0.14557772800000E+01
    f000000       0  0  0  0  0  0  1  0.96253006000000E+00
    f000000       0  0  0  0  0  0  1  0.10108194717000E+03
    f000000       0  0  0  0  0  0  1  0.11246000000000E+03
    f000000       0  0  0  0  0  0  1  0.39611300000000E-02
    f000001       0  0  0  0  0  1  1  0.48279121261621E-02
    f000002       0  0  0  0  0  2  1  0.31552076592194E-02
    f000003       0  0  0  0  0  3  1  0.27380895575892E-03
    f000004       0  0  0  0  0  4  1  0.53200000000000E-04



The Fortran function is ``MLpoten_h2o2_koput_unique``, which can found in ``mol_abcd.f90``.  A TROVE input example  for HOOH is  `HOOH_step1.inp <https://raw.githubusercontent.com/Trovemaster/TROVE/develop/docs/source/input/HOOH_step1.inpp>`_ from [15AlOvYu]_ , where it was used to compute an ExoMol line list for this molecule.






