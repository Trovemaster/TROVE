Symmetry
********
.. _symmetry:


TROVE uses the Molecule Symmetry Group [98BuJe]_ to classify the ro-vibrational states, motion, coordinates etc. The symmetries are defined in the ``molecule.f90`` file.
In order to specify the symmetry, the keyword ``SYMGROUP`` should be given anywhere in the input file outside any sections, but before the ``DIAGONALIZER`` section, e.g.
::

     SYMGROUP D2H(M)

Here the molecular symmetry group is :math:`D_{2h}`(M). 

C(M)
=====

C(M) is the simplest symmetry which means no symmetry with one irreducible representation (irrep) :math:`A`.


Cs(M)
=====

``Cs(M)`` is the second simplest symmetry group with two irreducible representations :math:`A'` and :math:`A''`:
::

     SYMGROUP Cs(M)

It is usually used for non-symmetric planar molecules such as triatomics XYZ.  


C2v(M)
=====


``C2v(M)`` is a molecular symmetry consisting of for irreps: :math:`A_1`, :math:`A_2`, :math:`B_1`, :math:`B_2`. The meaning of these irreps depends on the molecule  as well as its embedding. For an XY\ :sub:`2` molecule with the :math:`x` axis chosen as a bisector of the inter-bond angle, the characters are shown in the following table

+------------+-------------+------+-------------+----------------+
|            | :math:`E^*` | (12) | :math:`E^*` | (12)\ :sup:`*` |
+============+=============+======+=============+================+
|:math:`A_1` |      1      |  1   |       1     |       1        | 
+------------+-------------+------+-------------+----------------+
|:math:`A_2` |      1      |  1   |      -1     |      -1        | 
+------------+-------------+------+-------------+----------------+
|:math:`B_1` |      1      | -1   |      -1     |       1        | 
+------------+-------------+------+-------------+----------------+
|:math:`B_2` |      1      | -1   |       1     |      -1        | 
+------------+-------------+------+-------------+----------------+

