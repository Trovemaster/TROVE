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
|:math:`\alpha1` |:math:`\alpha1` | :math:`\alpha2`  | :math:`\alpha3`  | :math:`\alpha1`  |  :math:`\alpha3`  | :math:`\alpha2`  |
+----------------+----------------+------------------+------------------+------------------+-------------------+------------------+
|:math:`\alpha2` |:math:`\alpha2` | :math:`\alpha3`  | :math:`\alpha1`  | :math:`\alpha3`  |  :math:`\alpha2`  | :math:`\alpha1`  |
+----------------+----------------+------------------+------------------+------------------+-------------------+------------------+
|:math:`\alpha3` |:math:`\alpha3` | :math:`\alpha2`  | :math:`\alpha2`  | :math:`\alpha2`  |  :math:`\alpha1`  | :math:`\alpha3`  |
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
