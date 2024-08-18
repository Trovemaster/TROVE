========
Z-matrix
========


TROVE uses the standard Z-matrix scheme to define the molecular structure and introduce the basic internal coordinates, which are Z-matrix coordinates in TROVE. For, example, the hydrogen peroxide HOOH coordinate system is defined using the following Z-matrix defined with the ``ZMAT`` card: 
::

   ZMAT
       O1  0  0  0  0  15.99491463
       O2  1  0  0  0  15.99491463
       H1  1  2  0  0   1.00782505
       H2  2  1  3  2   1.00782505
   end

.. note: Zmatrix is also used to introduce the atomic (or nuclear) masses. 

Here, the Zmatrix coordinates are as follows.


.. sidebar::

   .. figure:: img/A2B2.jpg
       :alt: A2B2

       Valence coordinates used for HOOH.



- Three valence bond length (stretching) coordinates 1,2,3:

.. math::
      
      \begin{split}
       \xi_1 &= r_{{\rm O}_1{\rm O}_2} \\
       \xi_2 &= r_{{\rm O}_1{\rm H}_1} \\
       \xi_3 &= r_{{\rm O}_2{\rm H}_2} 
     \end{split}
    

- Two valence bond-angle (bending) coordinates 4,5:

.. math::
      
      \begin{split}
       \xi_4 &= \alpha_{{\rm H}_1{\rm O}_1 {\rm O}_2} \\
       \xi_5 &= \alpha_{{\rm H}_2{\rm O}_2 {\rm O}_3}
      \end{split}
      

- One dihedral coordinate defined as a book angle between planes :math:`{\rm H}_1{\rm O}_1 {\rm O}_2` and {\rm H}_2{\rm O}_2 {\rm O}_3: 

.. math::

   \xi_6 = \delta_{{\rm H}_1{\rm O}_1 {\rm O}_2 {\rm H}_2}.
   
.. note: The order of the coordinates in TROVE is always: stretching, bending and dihedrals. 



Structure of the Zmatrix block
------------------------------

Let us consider the HOOH example above to describe the columns in the ``ZMAT`` block:


+---------+-------------+---------------+--------------+----------+------------+
|      0  |   1         |     2         |       3      |    4     |       5    |  
+---------+-------------+---------------+--------------+----------+------------+
|  atom   | connector 1 | connector  2  | connector  3 |  Type    | Mass       |
+---------+-------------+---------------+--------------+----------+------------+
|      O1 |    0        |       0       |       0      |    0     | 15.99491463|
+---------+-------------+---------------+--------------+----------+------------+
|      O2 |    1        |       0       |       0      |    0     | 15.99491463|
+---------+-------------+---------------+--------------+----------+------------+
|      H1 |    1        |       2       |       0      |    0     |  1.00782505|
+---------+-------------+---------------+--------------+----------+------------+
|      H2 |    2        |       1       |       3      |    2     |  1.00782505|
+---------+-------------+---------------+--------------+----------+------------+


Here:

- atom : A label for the name of the atom; it is not interpreted by the program and is currently only for clarity. 
- connector 1: The 1st connecting atom to form a molecular bond.
- connector 2: The 2nd connecting atom to form a bond angle.
- connector 3: The 3rd connecting atom to form a dihedral angle or other type angle, depending on the value in column Type.
- Type: A "Dihedral" angle type (see below).
- Mass: The mass of the particle; usually an atomic but sometimes a nuclear mass.




"Dihedral" types
----------------

The following "Dihedral" types are available: 

- Type 0: the "Dihedral"  angle  is defined as  a valence angle between bonds :math:`\vec{r_{01}}`  and  :math:`\vec{r_{03}}`. 
- Type 1: it is defined as  a usual dihedral angle  between two planes. 

**TBC**