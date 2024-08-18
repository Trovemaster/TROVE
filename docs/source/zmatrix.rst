========
Z-matrix
========


TROVE uses the standard Z-matrix scheme to define the molecular structure and introduce the basic internal coordinates, which are Z-matrix coordinates in TROVE. For, example, the hydrogen peroxide HOOH coordinate system is defined using the following Z-matrix:
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


   