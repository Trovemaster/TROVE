Symmetry
********
.. _symmetry:


TROVE uses the Molecule Symmetry Group [98BuJe]_ to classify the ro-vibrational states, motion, coordinates etc. The symmetries are defined in the ``molecule.f90`` file.
In order to specify the symmetry, the keyword ``SYMGROUP`` should be given anywhere in the input file outside any sections, but before the ``DIAGONALIZER`` section, e.g.
::

     SYMGROUP D2H(M)

Here the molecular symmetry group is D\ :sub:`2h`\ (M).

C(M)
=====

C(M) is the simplest symmetry which means no symmetry with one irreducible representation (irrep) :math:`A`.


Cs(M)
=====

``Cs(M)`` is the second simplest symmetry group with two irreducible representations :math:`A'` and :math:`A''`:
::

     SYMGROUP Cs(M)

It is usually used for non-symmetric planar molecules such as triatomics XYZ.


C\ :sub:`2v`\ (M)
=================


``C2v(M)`` is a molecular symmetry group consisting of 4 irreps: :math:`A_1`, :math:`A_2`, :math:`B_1`, :math:`B_2`. The meaning of these irreps depends on the molecule  as well as its embedding. The characters are shown in the following table

+------------+-------------+------+-------------+----------------+
|            | :math:`E`   |(12)\ | :math:`E^*` | (12)\ :sup:`*` |
+============+=============+======+=============+================+
|:math:`A_1` |      1      |  1   |       1     |       1        |
+------------+-------------+------+-------------+----------------+
|:math:`A_2` |      1      |  1   |      -1     |      -1        |
+------------+-------------+------+-------------+----------------+
|:math:`B_1` |      1      | -1   |      -1     |       1        |
+------------+-------------+------+-------------+----------------+
|:math:`B_2` |      1      | -1   |       1     |      -1        |
+------------+-------------+------+-------------+----------------+

where :math:`E`, (12), :math:`E^*` and (12)\ :sup:`*` are the 4 point group operations. Typical C\ :sub:`2v`\ (M) molecules are XY\ :sub:`2` (e.g. water, H\ :sub:`2`\ S), ZXY\ :sub:`2` (e.g. formaldehyde.)


Symmetry properties of the vibrational coordinates
--------------------------------------------------

TROVE uses the symmetry properties of the vibrational coordinates, i.e. how they transform upon applying the symmetry operations, to build symmetry adapted vibrational basis functions. The symmetrisation method is described in [17YuYaOv]_. In the following, we show how the coordinates, we show these transformation properties for the corresponding coordination frames implemented in TROVE.


``R-RHO-Z``, ``R-ALPHA-Z``
^^^^^^^^^^^^^^^^^^^^^^^^^^


+--------------+---------------+--------------+---------------+---------------+
| Coordinates  | :math:`E`     |     (12)\    |  :math:`E^*`  |(12)\ :sup:`*` |
+==============+===============+==============+===============+===============+
|:math:`r_1`   | :math:`r_1`   |  :math:`r_2` |:math:`r_1`    |  :math:`r_2`  |
+--------------+---------------+--------------+---------------+---------------+
|:math:`r_2`   | :math:`r_2`   |  :math:`r_1` |:math:`r_2`    |  :math:`r_1`  |
+--------------+---------------+--------------+---------------+---------------+
| :math:`\rho` |  :math:`\rho` | :math:`\rho` |:math:`\rho`   | :math:`\rho`  |
+--------------+---------------+--------------+---------------+---------------+
|:math:`\alpha`| :math:`\alpha`|:math:`\alpha`|:math:`\alpha` | :math:`\alpha`|
+--------------+---------------+--------------+---------------+---------------+





C3v(M)
======

``C3v(M)`` is a molecular symmetry group consisting of 3 irreps: :math:`A_1`, :math:`A_2`, :math:`E` and 6 operations:


+------------+------------+-------+----------------+
|            | :math:`E`  |(123)\ | (23)\ :sup:`*` |
|            |            +-------+----------------+
|            |            |(132)\ | (12)\ :sup:`*` |
|            |            |       +----------------+
|            |            |       | (23)\ :sup:`*` |
+============+============+=======+================+
|:math:`A_1` |      1     |   1   |       1        |
+------------+------------+-------+----------------+
|:math:`A_2` |      1     |   1   |      -1        |
+------------+------------+-------+----------------+
|:math:`E`   |      2     |  -1   |       0        |
+------------+------------+-------+----------------+

It can be used for the molecules PH\ :sub:`3` [15SoAlTe]_, SbH\ :sub:`3` [10YuCaYa]_, AsH\ :sub:`3` [19CoYuKo]_, PF\ :sub:`3` [19MaChYa]_, CH\ :sub:`3`\ Cl [18OwYaTe]_, CH\ :sub:`3`\ F, isotopologue CDH\ :sub:`3` etc.


Coordinate transformation properties
------------------------------------


``R-ALPHA``
^^^^^^^^^^^

This is a rigid frame with 6 valence coordinates :math:`r_1`, :math:`r_2`, :math:`r_3`, :math:`\alpha_1`, :math:`\alpha_2` and :math:`\alpha_3`.


+----------------+----------------+------------------+------------------+------------------+-------------------+------------------+
| Coordinates    |     :math:`E`  |    (123)\        |       (132)\     |    (23)\ :sup:`*`|    (13)\ :sup:`*` |    (12)\ :sup:`*`|
+================+================+==================+==================+==================+===================+==================+
|    :math:`r_1` |    :math:`r_1` |     :math:`r_2`  |     :math:`r_3`  |     :math:`r_1`  |      :math:`r_3`  |     :math:`r_2`  |
+----------------+----------------+------------------+------------------+------------------+-------------------+------------------+
|    :math:`r_2` |    :math:`r_2` |     :math:`r_3`  |     :math:`r_1`  |     :math:`r_3`  |      :math:`r_2`  |     :math:`r_1`  |
+----------------+----------------+------------------+------------------+------------------+-------------------+------------------+
|    :math:`r_3` |    :math:`r_3` |     :math:`r_2`  |     :math:`r_2`  |     :math:`r_2`  |      :math:`r_1`  |     :math:`r_3`  |
+----------------+----------------+------------------+------------------+------------------+-------------------+------------------+
|:math:`\alpha_1`|:math:`\alpha_1`| :math:`\alpha_2` | :math:`\alpha_3` | :math:`\alpha_1` |  :math:`\alpha_3` | :math:`\alpha_2` |
+----------------+----------------+------------------+------------------+------------------+-------------------+------------------+
|:math:`\alpha_2`|:math:`\alpha_2`| :math:`\alpha_3` | :math:`\alpha_1` | :math:`\alpha_3` |  :math:`\alpha_2` | :math:`\alpha_1` |
+----------------+----------------+------------------+------------------+------------------+-------------------+------------------+
|:math:`\alpha_3`|:math:`\alpha_3`| :math:`\alpha_2` | :math:`\alpha_2` | :math:`\alpha_2` |  :math:`\alpha_1` | :math:`\alpha_3` |
+----------------+----------------+------------------+------------------+------------------+-------------------+------------------+



D\ :sub:`3h`\ (M)
=================

``D3h(M)`` is a molecular symmetry consisting of 6 irreps: :math:`A'_1`, :math:`A'_2`, :math:`E'`, :math:`A''_1`, :math:`A''_2`, :math:`E''` and 12 operations:


+-------------+------------+-------+-------------+------------------+-----------------+----------------+
|             | :math:`E`  |(123)\ | (23)\       |:math:`E^*`       | (123)\ :sup:`*` | (23)\ :sup:`*` |
|             |            +-------+-------------+                  +-----------------+----------------+
|             |            |(132)\ | (12)\       |                  | (132)\ :sup:`*` | (12)\ :sup:`*` |
|             |            |       +-------------+                  |                 +----------------+
|             |            |       | (23)\       |                  |                 | (23)\ :sup:`*` |
+=============+============+=======+=============+==================+=================+================+
|:math:`A'_1` |      1     |   1   |       1     |      1           |   1             |       1        |
+-------------+------------+-------+-------------+------------------+-----------------+----------------+
|:math:`A'_2` |      1     |   1   |      -1     |      1           |   1             |      -1        |
+-------------+------------+-------+-------------+------------------+-----------------+----------------+
|:math:`E'`   |      2     |  -1   |       0     |      2           |  -1             |       0        |
+-------------+------------+-------+-------------+------------------+-----------------+----------------+
|:math:`A''_1`|      1     |   1   |       1     |     -1           |  -1             |      -1        |
+-------------+------------+-------+-------------+------------------+-----------------+----------------+
|:math:`A''_2`|      1     |   1   |      -1     |     -1           |  -1             |       1        |
+-------------+------------+-------+-------------+------------------+-----------------+----------------+
|:math:`E''`  |      2     |  -1   |       0     |     -2           |   1             |       0        |
+-------------+------------+-------+-------------+------------------+-----------------+----------------+


The D\ :sub:`3h`\ (M) group has been used for NH\ :sub:`3` [10CoYuTe]_, CH\ :sub:`3` [19AdJeYa]_.




``R-S-DELTA``
^^^^^^^^^^^^^

This is a non-rigid frame with 3 valence stretching coordinates :math:`r_1`, :math:`r_2`, :math:`r_3`, a symmetry adapted bending vector :math:`(S_a,S_b)` and an umbrella coordinate :math:`\delta`, where

.. math::

    \begin{split}
    S_a &= \frac{1}{\sqrt{6}} (2 \alpha_{23}-\alpha_{13}-\alpha_{12}),  \\
    S_b &= \frac{1}{\sqrt{2}} ( \alpha_{13}-\alpha_{12}).  \\
    \end{split}



The transformation properties of the stretching coordinates are given by


+----------------+----------------+------------------------+------------------------+---------------------+---------------------+---------------------+
| Coordinates    |     :math:`E`  |(123)\ , (123)\ :sup:`*`|(321)\ , (321)\ :sup:`*`|(23)\ ,(23)\ :sup:`*`|(12)\ ,(12)\ :sup:`*`|(13)\ ,(13)\ :sup:`*`|
+================+================+========================+========================+=====================+=====================+=====================+
|    :math:`r_1` |    :math:`r_1` |           :math:`r_2`  |           :math:`r_3`  |        :math:`r_1`  |        :math:`r_2`  |        :math:`r_3`  |
+----------------+----------------+------------------------+------------------------+---------------------+---------------------+---------------------+
|    :math:`r_2` |    :math:`r_2` |           :math:`r_3`  |           :math:`r_1`  |        :math:`r_3`  |        :math:`r_1`  |        :math:`r_2`  |
+----------------+----------------+------------------------+------------------------+---------------------+---------------------+---------------------+
|    :math:`r_2` |    :math:`r_3` |           :math:`r_1`  |           :math:`r_2`  |        :math:`r_2`  |        :math:`r_3`  |        :math:`r_1`  |
+----------------+----------------+------------------------+------------------------+---------------------+---------------------+---------------------+

The bending vector :math:`(S_a,S_b)` transforms as follows

.. math::

    {\bf S} = {\bf D}(G) {\bf S},

where :math:`{\bf D}(G)` are 2x2 transformation matrices given by


.. math::

    \begin{split}
    {\bf D}(E) &= \left( \begin{array}{cc}
                           1 & 0 \\
                           0 & 1 \\
                          \end{array}
                  \right) \\
    {\bf D}(123) &= \left( \begin{array}{cc}
                           -\frac{1}{2} & \frac{\sqrt{3}}{2}  \\
                           -\frac{\sqrt{3}}{2} & -\frac{1}{2}  \\
                          \end{array}
                  \right) \\
    {\bf D}(321) &= \left( \begin{array}{cc}
                           -\frac{1}{2} & -\frac{\sqrt{3}}{2}  \\
                            \frac{\sqrt{3}}{2} & -\frac{1}{2}  \\
                          \end{array}
                  \right) \\
    {\bf D}(23) &= \left( \begin{array}{cc}
                           1 &  0 \\
                           0 & -1 \\
                          \end{array}
                  \right) \\
    {\bf D}(12) &= \left( \begin{array}{cc}
                           -\frac{1}{2} &  \frac{\sqrt{3}}{2}  \\
                            \frac{\sqrt{3}}{2} &  \frac{1}{2}  \\
                          \end{array}
                  \right) \\
    {\bf D}(13) &= \left( \begin{array}{cc}
                           -\frac{1}{2} & -\frac{\sqrt{3}}{2}  \\
                           -\frac{\sqrt{3}}{2} &  \frac{1}{2}  \\
                          \end{array}
                  \right) \\
        \end{split}

The operations with inversion have the same matrices, :math:`{\bf D}(G^*) = {\bf D}(G)`.

Finally, the umbrella coordinate transform as follows

.. math::

    G \delta = \left\{ \begin{array}{cc}
                               \delta  &  G=  E ,(123) ,(321), (23), (12), (13),\\
                               -\delta &   G = E^* ,(123)^* ,(321)^*,(23)^*,(12)^*,(13)^*.\\
                             \end{array} \right.





T\ :sub:`d`\ (M)
=================

``Td(M)`` is a molecular symmetry group is used for the methane-like molecules, CH\ :sub:`4` [14YuJo]_, SiH\ :sub:`3` [17OwYuYa]_. It consists of 5 irreps and 24 symmetry operations spanning 5 classes:


+-------------+------------+-------+-------------+------------------+----------------+
|             | :math:`E`  |(123)\ | (14)\ (23)\ |(1423)\ :sup:`*`  | (23)\ :sup:`*` |
+=============+============+=======+=============+==================+================+
|Elements     |      1     |   8   |      3      |      6           |       6        |
+-------------+------------+-------+-------------+------------------+----------------+
|:math:`A_1`  |      1     |   1   |       1     |      1           |       1        |
+-------------+------------+-------+-------------+------------------+----------------+
|:math:`A_2`  |      1     |   1   |       1     |     -1           |      -1        |
+-------------+------------+-------+-------------+------------------+----------------+
|:math:`E`    |      2     |  -1   |       2     |      0           |       0        |
+-------------+------------+-------+-------------+------------------+----------------+
|:math:`F_1`  |      3     |   0   |      -1     |      1           |      -1        |
+-------------+------------+-------+-------------+------------------+----------------+
|:math:`F_2`  |      3     |   0   |      -1     |     -1           |       1        |
+-------------+------------+-------+-------------+------------------+----------------+


Artificial symmetries (AEM)
===========================

The concept of artificial molecular symmetries was introduced in [21MeYuJe]_. 

C\ :sub:`2vn`\ (AEM)
-------------------


Linear molecules usually represent a special case in rotational-vibrational calculations due to a singularity of the kinetic energy operator that arises from the rotation about the :math:`a` (the principal axis of least moment of inertia, becoming the molecular axis at the linear equilibrium geometry) being undefined. Assuming the standard ro-vibrational basis functions, in the :math:`3N-6` approach,  of the form :math:`\ket{\nu_1, \nu_2, \nu_3^{l_3}; J, k, m}`, tackling the unique difficulties of linear molecules involves constraining the vibrational and rotational functions with :math:`k=l_3`,  which are the projections, in units of :math:`\hbar`, of the corresponding angular momenta onto the molecular axis. These basis functions are assigned to irreps of the C\ :math:`_{2{\rm v}}`\ (M) molecular symmetry group. This, in turn, necessitates purpose-built codes that specifically deal with linear molecules. In the present work, we describe an alternative scheme and introduce an (artificial) group that ensures that the condition :math:`l_3 =k` is automatically applied solely through symmetry group algebra. The advantage of such an approach is that the application of symmetry group algebra in ro-vibrational calculations is ubiquitous, and so this method can be used to enable ro-vibrational calculations of linear molecules in polyatomic codes with fairly minimal modifications.

In TROVE an alternative scheme is implemented as an (artificial) group that ensures that the condition :math:`l_3 =k` is automatically applied solely through symmetry group algebra. The advantage of such an approach is that the application of symmetry group algebra in ro-vibrational calculations is ubiquitous, and so this method can be used to enable ro-vibrational calculations of linear molecules in polyatomic codes with fairly minimal modifications. To this end, we construct an artificial molecular symmetry  group C\ :sub:`2vn`\ (AEM), which consists of one-dimensional (non-degenerate) irreducible representations and use it to classify vibrational and rotational basis  functions according to :math:`l` and :math:`k`.  This extension to non-rigorous, artificial symmetry groups  is based on cyclic groups of prime-order. Opposite to the usual scenario, where the form of symmetry adapted basis sets is dictated by the symmetry group the molecule belongs to, here the symmetry  group C\ :sub:`2vn`\ (AEM) is built to satisfy properties for the convenience of the basis set construction and matrix elements calculations. We believe that the idea of  purpose-built artificial symmetry groups can be useful in other~applications.

Examples of character tables for :math:`n=4` are given in Table below 

+---------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+--------------+--------------+------------------+-------------------+
|C\ :sub:`2vn`\ (AEM) | :math:`E^0` |:math:`C_2^0` | :math:`\sigma^0` |:math:`\sigma_v^0` | :math:`E^1` |:math:`C_2^1` | :math:`\sigma^1` |:math:`\sigma_v^1` | :math:`E^2` |:math:`C_2^2` | :math:`\sigma^2` |:math:`\sigma_v^2` |  :math:`E^3` |:math:`C_2^3` | :math:`\sigma^3` |:math:`\sigma_v^3` |
+=====================+=============+==============+==================+===================+=============+==============+==================+===================+=============+==============+==================+===================+==============+==============+==================+===================+
|       :math:`A_1^0` |    1        |       1      |      1           |         1         |        1    |     1        |        1         |      1            |      1      |     1        |         1        |      1            |         1    |        1     |      1           |      1            |
+---------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+--------------+--------------+------------------+-------------------+
|       :math:`B_1^0` |    1        |      -1      |      1           |        -1         |        1    |    -1        |        1         |     -1            |      1      |    -1        |         1        |     -1            |         1    |       -1     |      1           |     -1            |
+---------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+--------------+--------------+------------------+-------------------+
|       :math:`A_2^0` |    1        |      1       |      -1          |        -1         |        1    |     1        |       -1         |     -1            |      1      |     1        |        -1        |     -1            |         1    |        1     |     -1           |     -1            |
+---------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+--------------+--------------+------------------+-------------------+
|       :math:`B_2^0` |    1        |      -1      |      -1          |         1         |        1    |    -1        |       -1         |      1            |      1      |    -1        |        -1        |      1            |         1    |       -1     |     -1           |      1            |
+---------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+--------------+--------------+------------------+-------------------+
|       :math:`A_1^1` |    1        |       1      |      1           |         1         |       -1    |    -1        |       -1         |     -1            |      1      |     1        |         1        |      1            |        -1    |       -1     |     -1           |     -1            |
+---------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+--------------+--------------+------------------+-------------------+
|       :math:`B_1^1` |    1        |      -1      |      1           |        -1         |       -1    |     1        |       -1         |      1            |      1      |    -1        |         1        |     -1            |        -1    |        1     |     -1           |      1            |
+---------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+--------------+--------------+------------------+-------------------+
|       :math:`A_2^1` |    1        |       1      |      -1          |        -1         |       -1    |    -1        |        1         |      1            |      1      |     1        |        -1        |     -1            |        -1    |       -1     |     -1           |      1            |
+---------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+--------------+--------------+------------------+-------------------+
|       :math:`B_2^1` |    1        |      -1      |      -1          |         1         |       -1    |     1        |        1         |     -1            |      1      |    -1        |        -1        |      1            |        -1    |        1     |      1           |     -1            |
+---------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+--------------+--------------+------------------+-------------------+
|       :math:`A_1^2` |    1        |       1      |      1           |         1         |        1    |     1        |        1         |      1            |      -1     |    -1        |        -1        |     -1            |      -1      |       -1     |     -1           |     -1            |
+---------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+--------------+--------------+------------------+-------------------+
|       :math:`B_1^2` |    1        |      -1      |      1           |        -1         |        1    |    -1        |        1         |     -1            |      -1     |     1        |        -1        |      1            |      -1      |        1     |     -1           |      1            |
+---------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+--------------+--------------+------------------+-------------------+
|       :math:`A_2^2` |    1        |      1       |      -1          |        -1         |        1    |     1        |       -1         |     -1            |      -1     |    -1        |         1        |      1            |      -1      |       -1     |      1           |      1            |
+---------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+--------------+--------------+------------------+-------------------+
|       :math:`B_2^2` |    1        |      -1      |      -1          |         1         |        1    |    -1        |       -1         |      1            |      -1     |    -1        |         1        |      1            |      -1      |        1     |      1           |     -1            |
+---------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+--------------+--------------+------------------+-------------------+
|       :math:`A_1^3` |    1        |       1      |      1           |         1         |       -1    |    -1        |       -1         |     -1            |      -1     |    -1        |        -1        |     -1            |       1      |        1     |      1           |      1            |
+---------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+--------------+--------------+------------------+-------------------+
|       :math:`B_1^3` |    1        |      -1      |      1           |        -1         |       -1    |     1        |       -1         |      1            |      -1     |     1        |        -1        |      1            |       1      |       -1     |      1           |     -1            |
+---------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+--------------+--------------+------------------+-------------------+
|       :math:`A_2^3` |    1        |       1      |      -1          |        -1         |       -1    |    -1        |        1         |      1            |      -1     |    -1        |         1        |      1            |       1      |        1     |     -1           |     -1            |
+---------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+--------------+--------------+------------------+-------------------+
|       :math:`B_2^3` |    1        |      -1      |      -1          |         1         |       -1    |     1        |        1         |     -1            |      -1     |     1        |         1        |     -1            |       1      |       -1     |     -1           |      1            |
+---------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+-------------+--------------+------------------+-------------------+--------------+--------------+------------------+-------------------+

How to use C\ :sub:`2vn`\ (AEM)   
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An artificial symmetry C\ :sub:`2vn`\ (AEM)    is invoked via the following card places anywhere (but before the ``diagonalizer`` section):
::

    SYMGROUP C2vn 18

Here the integer number :math:`n=18` corresponds to the maximal value of the vibrational angular momentum :math:`l_{\rm max}` and therefore to the maximal value of the rotational quantum number :math:`k` due to the constraint :math:`k=l` used for linear triatomic molecules (case 3N-6). This number muster coincide with the value of the ``k`` or ``kmax`` cards in the ``BASIS`` block (rotational basis line), e.g.
::

    BASIS
      0,'JKtau', Jrot 0, krot  18
      1,'numerov','rational', 'morse',  range 0,20, r 8, resc 3.0, points   2000, borders -0.3,0.90
      2,'numerov','rational', 'morse',  range 0,36, r 8, resc 1.5, points   3000, borders -0.4,0.90
      3,'laguerre-k','linear','linear', range 0,48, r 8, resc 1.0, points  12000, borders  0.,120.0 deg
    END



C\ :sub:`ns`\ (AEM) 
-------------------


The artificial symmetry group C\ :sub:`ns`\ (AEM) consists of one-dimensional, real irreps, which are contracted to correlate with irreps of C\ :math:`_{\infty v}`. The irreps of C\ :sub:`ns`\ (AEM) are labelled as :math:`\Gamma=A'` and :math:`A''` irreps with an extra subscript (see Table below), e.g., :math:`A'_4`.  For example, a vibrational function with :math:`l=4` and transforming as :math:`A'` in C\ :sub:`s` would be assigned the symmetry :math:`A'_4` in C\ :sub:`ns`\ (AEM). The 0-superscripted irreps are the only physical irreps, matched to :math:`A'` and :math:`A''` of \Cs\ together with the corresponding characters of each element, while all irreps of :math:`l>0` are non-physical, i.e. "artificial". The full description of this case is given in [24YuMeTe]_.

+--------------+-----------+-------------+-----------------+-------------+-----------------+--------------+-----------------+-------------+----------------+
|C\ :sub:`4s`  |C\ :sub:`s`| :math:`E^0` |:math:`\sigma^0` | :math:`E^1` |:math:`\sigma^1` |  :math:`E^2` |:math:`\sigma^2` | :math:`E^3` |:math:`\sigma^3`|
+==============+===========+=============+=================+=============+=================+==============+=================+=============+================+
|:math:`A'_0`  |:math:`A'` |   1         |   1             |       1     |        1        |       1      |        1        |          1  |   1            |
+--------------+-----------+-------------+-----------------+-------------+-----------------+--------------+-----------------+-------------+----------------+
|:math:`A''_0` |:math:`A''`|   1         |  -1             |       1     |       -1        |       1      |       -1        |          1  |  -1            |
+--------------+-----------+-------------+-----------------+-------------+-----------------+--------------+-----------------+-------------+----------------+
|:math:`A'_1`  |           |   1         |   1             |      -1     |       -1        |       1      |        1        |         -1  |  -1            |
+--------------+-----------+-------------+-----------------+-------------+-----------------+--------------+-----------------+-------------+----------------+
|:math:`A''_1` |           |   1         |  -1             |      -1     |        1        |       1      |       -1        |         -1  |   1            |
+--------------+-----------+-------------+-----------------+-------------+-----------------+--------------+-----------------+-------------+----------------+
|:math:`A'_2`  |           |   1         |   1             |       1     |        1        |      -1      |       -1        |         -1  |  -1            |
+--------------+-----------+-------------+-----------------+-------------+-----------------+--------------+-----------------+-------------+----------------+
|:math:`A''_2` |           |   1         |  -1             |       1     |       -1        |      -1      |        1        |         -1  |   1            |
+--------------+-----------+-------------+-----------------+-------------+-----------------+--------------+-----------------+-------------+----------------+
|:math:`A'_3`  |           |   1         |   1             |      -1     |       -1        |      -1      |       -1        |          1  |   1            |
+--------------+-----------+-------------+-----------------+-------------+-----------------+--------------+-----------------+-------------+----------------+
|:math:`A''_3` |           |   1         |  -1             |      -1     |        1        |      -1      |        1        |          1  |  -1            |
+--------------+-----------+-------------+-----------------+-------------+-----------------+--------------+-----------------+-------------+----------------+

The  effects of the C\ :sub:`ns`\ (AEM) group operations on the coordinates is as follows: all operations leave the vibrational coordinates invariant; the  :math:`E^a` operations (in the notation of Table above leave the rotational functions invariant while the :math:`\sigma^a` operation has the same effect as the :math:`\sigma^0` operation. 

How to use C\ :sub:`ns`\ (AEM) 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An artificial symmetry C\ :sub:`ns`\ (AEM) is invoked via the following card places anywhere (but before the ``diagonalizer`` section):
::
    
    SYMGROUP Csn 18
   
The integer number :math:`n=18` corresponds to the maximal value of the vibrational angular momentum :math:`l_{\rm max}` and therefore to the maximal value of the rotational quantum number :math:`k` due to the constraint :math:`k=l` used for linear triatomic molecules (case 3N-6). This number muster coincide with the value of the ``k`` or ``kmax`` cards in the ``BASIS`` block (rotational basis line) (see above). 





