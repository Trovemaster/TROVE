Symmetry
********
.. _symmetry:


TROVE uses the Molecule Symmetry Group [98BuJe]_ to classify the ro-vibrational states, motion, coordinates etc. The symmetries are defined in the ``molecule.f90`` file.
In order to specify the symmetry, the keyword ``SYMGROUP`` should be given anywhere in the input file outside any sections, but before the ``DIAGONALIZER`` section, e.g.
::

     SYMGROUP D2H(M)

Here the molecular symmetry group is D\ :sub:`2h`(M). 

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


``C2v(M)`` is a molecular symmetry consisting of 4 irreps: :math:`A_1`, :math:`A_2`, :math:`B_1`, :math:`B_2`. The meaning of these irreps depends on the molecule  as well as its embedding. The characters are shown in the following table

+------------+-------------+------+-------------+----------------+
|            | :math:`E`   | (12) | :math:`E^*` | (12)\ :sup:`*` |
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

C3v(M)
======

``C3v(M)`` is a molecular symmetry consisting of 3 irreps: :math:`A_1`, :math:`A_2`, :math:`E` and 6 operations: 


============= ============ ======= ================
               :math:`E`    (123)   (23)\ :sup:`*`
 ------------              ------- ----------------
                            (132)   (12)\ :sup:`*`
-------------              ------- ----------------
                                    (23)\ :sup:`*`
============= ============ ======= ================
 :math:`A_1`        1         1           1
 :math:`A_2`        1         1          -1
 :math:`E`          2        -1           0
============= ============ ======= ================

It can be used for the molecules PH\ :sub:`3`, SbH\ :sub:`3`, SO\ :sub:`3`, CH\ :sub:`3`\ Cl, CH\ :sub:`3`\ Cl, isotopologue CDH\ :sub:`3` etc.


D\ :sub:`3h`\ (M)
=================

``D3h(M)`` is a molecular symmetry consisting of 6 irreps: :math:`A'_1`, :math:`A'_2`, :math:`E'`, :math:`A''_1`, :math:`A''_2`, :math:`E''` and 12 operations:


============= ============ ======= ================
               :math:`E`    (123)   (23)\ :sup:`*`
------------- ------------ ------- ----------------
                            (132)   (12)\ :sup:`*`
------------- ------------         ----------------
                                    (23)\ :sup:`*`
============= ============ ======= ================
 :math:`A_1`        1         1           1
 :math:`A_2`        1         1          -1
 :math:`E`          2        -1           0
============= ============ ======= ================

It can be used for the molecules PH\ :sub:`3`, SbH\ :sub:`3`, SO\ :sub:`3`, CH\ :sub:`3`Cl, CH\ :sub:`3`Cl, isotopologue CDH\ :sub:`3` etc.


============= ============ ======= ================
               :math:`E`    (123)   (23)\ :sup:`*`
                            (132)   (12)\ :sup:`*`
                                    (23)\ :sup:`*`
------------- ------------ ------- ----------------
 :math:`A_1`        1         1           1
 :math:`A_2`        1         1          -1
 :math:`E`          2        -1           0
============= ============ ======= ================

