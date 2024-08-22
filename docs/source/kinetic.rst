Kinetic energy operators
========================

There are several options for kinetic energy operators (KEO) in TROVE: 

- **Linearised KEO**: Automatically constructed KEO using the Eckart frame (rigid) or Eckart/Saywitz frame ((1D non-rigid) as a Taylor expansion in terms of linearised coordinates;
- **Curvilinear KEO**: Preprogrammed analytic KEO (exact or non-exact) using curvilinear (valence) coordinates and different frames;
- **External Numerical Taylor-type KEO**: Externally generated KEO as a Taylor-type numerical expansion in terms of any user-defined coordinates;
- **External sum-of-products KEO**: Externally generated KEO as a generic sum-of-products expansion via the TROVE KEO builder.


Linearised KEO
--------------

This is a standard KEO constructed using a general black-boxed algorithm applicable for an arbitrary polyatomic molecule. It has been introduced in the original TROVE paper [TROVE]_ and also explained in :doc:`theory`. For the following, it useful to recall that TROVE uses Z-matrix curvilinear coordinates to define the linearised coordinates for a given molecule and selected frame.

Linearised KEO is a default option in TROVE. It requires the following components to be present in the TROVE code:

- Coordinate transformation (``TRANSFORM``) between the Z-matrix curvilinear coordinates :math:`\xi_{i}^{\rm Zmat}` and the target curvilinear coordinates :math:`\xi_i` used to generate the linearised coordinates :math:`\xi_i^{\rm lin}`. For example, for the Z-matrix used for NH\ :sub:`3` (see :doc:`molecules`):
::


   ZMAT
       N   0  0  0  0  14.00307401
       H1  1  0  0  0   1.00782503223
       H2  1  2  0  0   1.00782503223
       H3  1  2  3  1   1.00782503223
   end

the Z-matrix coordinates :math:`\xi_{i}^{\rm Zmat}` are the three bond lengths :math:`r_1`, :math:`r_2` and :math:`r_3` and three bond angles :math:`\alpha_{12}`, :math:`\alpha_{13}` and :math:`\alpha_{23}`. The target :math:`\xi_{i}` curvilinear coordinates for the scheme (``TRANSFORM``) ``r-s-delta`` are

.. math::

    \begin{split}
      \xi_1 &= r_1-r_{\rm e} \\
      \xi_2 &= r_2-r_{\rm e} \\
      \xi_3 &= r_3-r_{\rm e} \\
      \xi_4 &= \frac{1}{6} (2\alpha_{23}-\alpha_{13}-\alpha_{12}) \\
      \xi_5 &= \frac{1}{2} (\alpha_{13}-\alpha_{12}) \\
      \xi_6& = \delta \\
    \end{split}


where  :math:`\delta = \frac{\pi}{2}-\rho`, where :math:`\rho` is the angle between the trisector and any of the bond vectors.  The five linearised coordinates :math:`\xi_i^{\rm lin}` for a non-rigid reference configuration around :math:`\delta` are constructed by linearising :math:`\xi_i`.


- Equilibrium or non-rigid reference structure associated with the frame (``Frame``) and/or coordinate transformation (``TRANSFORM``). For NH\ :sub:`3`, a non-rigid reference frame associated with the coordinates ``r-s-delta`` is  given by (in the centre-of-mass principle axis system):

.. math::

      \left( \begin{array}{ccc}
          0                              & 0                   & -\frac{3 m_H r_{\rm e} \sin\delta}{3 m_H+m_X} \\
         r_{\rm e} \cos\delta            & 0                   & \frac{m_X r_{\rm e} \cos\delta}{3 m_H+mX} \\
        -\frac{1}{2} r_{\rm e} \cos\delta&  \frac{\sqrt{3}}{2} r_{\rm e} \cos\delta & \frac{m_X r_{\rm e} \sin\delta}{3 m_H+mX} \\
        -\frac{1}{2} r_{\rm e} \cos\delta& -\frac{\sqrt{3}}{2} r_{\rm e} \cos\delta & \frac{m_X r_{\rm e} \sin\delta}{3 m_H+mX} \\
      \end{array}
      \right)

where :math:`\delta` is evaluated in a grid :math:`\delta_i = [-\delta_b \ldots  \delta_b]`.


.. sidebar::

    .. figure:: img/XY3_equil_delta.jpg
       :alt: XY3 equilibrium structure

       The non-rigid reference configuration for NH\ :sub:`3` in the trisector embedding and the ``r-s-delta``` coordinates/frame type.


For the XY\ :sub:`3` type (``MOLTYPE XY3``), the equilibrium structures are implemented as the subroutine ``ML_b0_XY3.f90`` as part of the module ``mol_xy3.f90``.

- The internal coordinates type needs to be set to `linear` (linearised).


Consider a KEO in the general form:

.. math::
       :label: e-H-total
           
       \begin{split}
       \hat{T}
         &= \frac{1}{2} \, \sum_{\alpha=x,y,z} \; \; \sum_{\alpha^\prime=x,y,z} \hat{J}_{\alpha}\, G_{\alpha,\alpha^\prime}(\xi)\, \hat{J}_{\alpha^\prime}   \\
         &+  \, \sum_{\alpha=x,y,z}\;\; \sum_{n=1}^{3N-6} \left[
                \hat{J}_{\alpha}\, G_{\alpha,\lambda}(\xi)\, \hat{p}_\lambda +
               \hat{p}_\lambda  \, G_{\alpha,\lambda}(\xi)\, \hat{J}_{\alpha} \right]  \\
         &+  \, \sum_{\lambda=1}^{M}\; \sum_{\lambda^\prime=1}^{M}
               \hat{p}_\lambda \, G_{\lambda,\lambda'}(\xi)\,  \hat{p}_{\lambda'} + U(\xi),
       \end{split}
       



Each term in the linearised KEO is represented by a Taylor expansion in terms of :math:`\xi_{i}^{\rm lin}` (or some 1D functions of them) around the (non-) rigid reference configuration, e.g.

.. math::
      :label: e-Taylor

      G_{n,n'} = \sum_{i,j,k,l,\ldots } G_{i,j,k,l,\ldots}^{(n,n')} \xi_1^i \xi_2^j \xi_3^k \xi_4^l \ldots

where the subscript 'lin' is omitted for compactness. This is a sum-of-products form which simplifies the integration with the basis functions.


Analytic curvilinear KEO
------------------------

A number of analytic curvilinear KEOs are available (see :doc:`frames`) and associated with the card ``kinetic_type`` as part of the ``KINETIC`` section, e.g.:
::

   KINETIC
     kinetic_type  KINETIC_XYZ_EKE_bisect
   END

Here are the implemented KEOs (quasi-linear triatomic molecules):

.. table:: t-KEOs
 
     Kinetic energy types implemented with associated basis set and frame types. 

+----------------+-------------------------------------+---------------------+--------------------------+
| ``MolType``    | ``kinetic_type``                    |  basis set          |    frames                |
+----------------+-------------------------------------+---------------------+--------------------------+
|    XY2         |   ``KINETIC_XY2_EKE_BISECT``        |``laguerre-k``       |  ``R-RHO-Z``             |
+----------------+-------------------------------------+---------------------+--------------------------+
|    XY2         |   ``KINETIC_XY2_EKE_BISECT_SINRHO`` |``sinrho-laguerre-k``|  ``R-RHO-Z``             |
+----------------+-------------------------------------+---------------------+--------------------------+
|    XYZ         |   ``KINETIC_XYZ_EKE_bisect``        |``laguerre-k``       |  ``R-RHO-Z-M2-M3-BISECT``|
+----------------+-------------------------------------+---------------------+--------------------------+
|    XYZ         |   ``KINETIC_XYZ_EKE_BOND``          |``laguerre-k``       |  ``R1-Z-R2-RHO``         |
+----------------+-------------------------------------+---------------------+--------------------------+
|    XYZ         |   ``KINETIC_XYZ_EKE_BOND_SINRHO``   |``sinrho-laguerre-k``|  ``R-RHO-Z-M2-M3-BISECT``|
+----------------+-------------------------------------+---------------------+--------------------------+
|    XYZ         |   ``KINETIC_XYZ_EKE_BOND-R2``       |``laguerre-k``       |  ``R2-Z-R1-RHO``         |
+----------------+-------------------------------------+---------------------+--------------------------+

These are included in the module ``kin_xy2.f90``.

All these analytic KEO are given in a sum-of-products of 1D terms, i.e. similar to the linearised type:

.. math::

    U = \sum_{l,m,n} u_{l,m,n}   f_{l}(r_1)f_{m}(r_2)f_{n}(\rho)

for the KEO coordinates chosen as :math:`r_1`, :math:`r_2` and :math:`\rho`. Once all the matrix elements of :math:`f_{l}(r_1)`, :math:`f_{m}(r_2)` and :math:`f_{n}(\alpha)` are computed, the same TROVE pipeline can be used for any types of  coordinates, regardless of if they are linearised or curvilinear.

In the case of analytic ``XYZ`` and ``XY2`` KEOs from :numref:`Table %s <t-KEOs>` above, the corresponding KEO terms :math:`G_{\lambda,\lambda'}` or pseudo-potential function :math:`U` depend on :math:`r_1` and :math:`r_2` as inverse expressions,  :math:`1/r_1`, :math:`1/r_1^2`, :math:`1/r_2`, :math:`1/r_2^2` and :math:`1/(r_1 r_2)`. We therefore represent them as sum-of-products expansions in terms of these 1D functions around the non-rigid configuration along the angle :math:`\rho`:

.. math::
    :label: e-G
    
    G_{\lambda,\lambda'} = \sum_{l,m} u_{l,m}\lambda,\lambda'   \frac{1}{r_1^l} \frac{1}{r_2^m} f_{n}(\rho)
    
for example:

.. math::
     
     \begin{split}
     G_{1,1}^{\rm vib} &= G_{2,2}^{\rm vib} = \frac{1}{\mu_{XY}},\\
     G_{1,2}^{\rm vib}&= G_{2,1}^{\rm vib} = \frac{\cos\alpha}{m_X},\\
     G_{1,3}^{\rm vib}&= G_{3,1}^{\rm vib} = -\frac{\sin\alpha}{r_2 m_X},\\
     G_{2,3}^{\rm vib}&= G_{3,2}^{\rm vib} = -\frac{\sin\alpha}{r_1 m_X},\\
     G_{3,3}^{\rm vib} &=  \frac{1}{\mu_{XY}} \left( \frac{1}{r_1^2} + \frac{1}{r_2^2} \right) -\frac{2 \cos\alpha}{r_1 r_2 m_X};
     \end{split}
     

These KEOs become singular at the linear geometry. This is resolved by combining them with special basis set functions that exactly cancel these singularities. Such special basis functions linked to the corresponding KEOs include: ``laguerre-k`` and ``sinrho-laguerre-k`` as listed in the table above. 

Moreover, the type of the 1D expansion (basic) terms in Eq. :eq:`e-G` should be set in the ``basis`` block for the appropriate mode. The expansion form in terms of  :math:`1/r_1` and :math:`1/r_2` are called ``rational``. Here is an example of the ``basis`` block associated with the KEO type ``KINETIC_XYZ_EKE_bisect``:
::
     
     BASIS
       0,'JKtau', Jrot 0, krot  18
       1,'numerov','rational', 'morse',  range 0,20, r 8, resc 3.0, points   2000, borders -0.3,0.90
       2,'numerov','rational', 'morse',  range 0,36, r 8, resc 1.5, points   3000, borders -0.4,0.90
       3,'laguerre-k','linear','linear', range 0,48, r 8, resc 1.0, points  12000, borders  0.,120.0 deg
     END
     
where the 3rd card  ``rational`` in the lines ``1`` and ``2``  indicate that the KEO is expanded in terms of :math:`1/r_1` and :math:`1/r_2`. Other related cards include: 
::
     
    KinOrder  2
    TRANSFORM    R-RHO-Z-M2-M3-BISECT
    MOLTYPE      XY2
    REFER-CONF   non-RIGID
    
    KINETIC
      kinetic_type  KINETIC_XYZ_EKE_bisect
    END
    
where ``KinOrder 2`` tells TROVE to expand the KEO :math:`1/r_1` and :math:`1/r_2` up to the 2nd order. 


For intensity calculations, it is also important to link these KEO to the appropriate frame. This is because the dipole moment vectors representations strongly depend on the molecular coordinate system. In the table above, the appropriate frames (``TRANSFORM``) are also listed. 

How to use analytic KEOs
************************

1. Define the KEO type using the ``kinetic`` block:
::

    KINETIC
      kinetic_type  KINETIC_XYZ_EKE_bisect
    END
    
2. Set the ``Coords`` type to ``local``, KEO expansion to 2,  choose the appropriate frame/and coordinates and molecule type, set the reference configuration to ``non-rigid``: 
::
    
    KinOrder  2    
    COORDS local (curvilinear)
    TRANSFORM   R-RHO-Z-M2-M3-BISECT (FRAME)
    MOLTYPE XY2
    REFER-CONF NON-RIGID
    
3. Choose appropriate basis set configuration, i.e. the KEO expansion type (e.g. ``rational``) and the non-rigid basis set type (e.g. ``laguerre-k``):
:: 

     BASIS
       0,'JKtau', Jrot 0, krot  18
       1,'numerov','rational', 'morse',  range 0,20, r 8, resc 3.0, points   2000, borders -0.3,0.90
       2,'numerov','rational', 'morse',  range 0,36, r 8, resc 1.5, points   3000, borders -0.4,0.90
       3,'laguerre-k','linear','linear', range 0,48, r 8, resc 1.0, points  12000, borders  0.,120.0 deg
     END



External Numerical Taylor-type KEO
----------------------------------

As was mentioned above, TROVE can work with any ro-vibrational KEOs as long as they are represented as sum-of-products of 1D functions of vibrational coordinates. It is possible to input any externally constructed (non-singular) KEO and to be used with the TROVE pipe line. One of the most robust methods is to use the TROVE checkpoints functionality.  Therefore, ``kinetic.chk`` can be used to input externally constructed KEOs in a numerical form, providing that the format (line and column orders, numbering) is preserved. 

Once produced by TROVE, KEO can be saved in an ascii file in a form of expansion coefficients, which fully define it. The structure of the kinetic.chk is explained in :doc:`checkpoints`. 
An example of externally constructed exact KEO for the H\ :sub:`2`\ CS molecule using Mathematica can be found in [23MeOwTe]_ (both as analytic formulas and a Mathematica ``.nb`` script), where it was used with TROVE to compute an ExoMol line list. It is given in a sum-of-products form:

.. math::
    :label:  e-G-H2CS

    G_{\lambda,\lambda'}(r_1,r_2,r_3,\alpha_1,\alpha_2,\tau) = \sum_{l,m} u_{l,m,n,o,p,q}\lambda,\lambda'   u_{l}(r_1) u_{m}(r_2) u_{n}(r_3) u_{o}(\alpha_1) u_{p}(\alpha_2) u_{o}(\tau).
    

.. sidebar::

    .. figure:: img/H2CS_coords.jpg
       :alt: H2CS coordinates

       The curvilinear  coordinates used for H\ :sub:`2`\ CS [23MeOwTe]_ (type ``R-THETA-TAU``).


This coordinates type (``TRANSFORM``) and associated frame is ``R-THETA-TAU`` (see :doc:`molecules`). 

The Mathemtica script ``Supp_Info_kinetic_energy_generator_h2cs_paper.nb`` (see supplementary to [23MeOwTe]_) was used to produce the ``kinetic.chk`` in a numerical format using these coordinates in a rigid configuration. This can be best explained by a example, see the TROVE input at `MOTY spectroscopic model <https://exomol.com/models/H2CS/1H2-12C-32S/MOTY/>`__.

Let us start with the ``basis`` set block: 
::
    
   BASIS
     0,'JKtau', Jrot 0
     1,'numerov','BOND-LENGTH', 'morse',  range 0,17, resc 1.0, points 1000, borders -0.35,1.00
     2,'numerov','BOND-LENGTH', 'morse',  range 0, 7, resc 2.0, points 1000, borders -0.5,0.9
     2,'numerov','BOND-LENGTH', 'morse',  range 0, 7, resc 2.0, points 1000, borders -0.5,0.9
     3,'numerov','ANGLE',       'linear', range 0,14, resc 1.0, points 1000, borders -1.2,1.2
     3,'numerov','ANGLE',       'linear', range 0,14, resc 1.0, points 1000, borders -1.2,1.2
     4,'numerov','DIHEDRAL',    'linear', range 0,14, resc 1.0, points 2000, borders -120.,  120.0 deg
   END
   
Here, the expansion types ``BOND-LENGTH``, ``ANGLE`` and ``DIHEDRAL`` are implemented in TROVE and fully define the 1D expansion terms :math:`u_{l}(\xi_l)` in Eq. :eq:`e-G-H2CS`. In TROVE, the expansion index :math:`l` is used as an analogy for the exponent in the Taylor type expansion, Eq :eq:`e-Taylor`. The correlation between these  indices and their actual forms are as follows. 



``BOND-LENGTH``
^^^^^^^^^^^^^^^

+-------+-----------------+
| Index | Term            |
+-------+-----------------+
| 0     |   1             |
+-------+-----------------+
| 1     |   :math:`1/r`   |
+-------+-----------------+
| 2     |   :math:`1/r^2` |
+-------+-----------------+


``ANGLE``
^^^^^^^^^

+-------+-------------------------------------+
| Index | Term                                |
+-------+-------------------------------------+
| 0     |   1                                 |
+-------+-------------------------------------+
| 1     |   :math:`\cos\alpha`                |
+-------+-------------------------------------+
| 2     |   :math:`1/\tan\alpha`              |
+-------+-------------------------------------+
| 3     |   :math:`1/\sin\alpha`              |
+-------+-------------------------------------+
| 4     |   :math:`\sin\alpha`                |
+-------+-------------------------------------+
| 5     |   :math:`1/\tan^2\alpha`            |
+-------+-------------------------------------+
| 6     |   :math:`1/(\sin\alpha \tan\alpha)` |
+-------+-------------------------------------+
| 5     |   :math:`1/\sin^2\alpha`            |
+-------+-------------------------------------+


``DIHEDRAL``
^^^^^^^^^^^^

+-------+------------------------------------+
| Index | Term                               |
+-------+------------------------------------+
| 0     |   1                                |
+-------+------------------------------------+
| 1     |   :math:`\cos(\tau/2)`             |
+-------+------------------------------------+
| 2     |   :math:`\sin(\tau/2)`             |
+-------+------------------------------------+
| 3     |   :math:`\cos^2(\tau/2)`           |
+-------+------------------------------------+
| 4     |   :math:`\cos(\tau/2)\sin(\tau/2)` |
+-------+------------------------------------+
| 5     |   :math:`\sin^2(\tau/2)`           |
+-------+------------------------------------+

Other relevant input cards include:
::
     
     COORDS     local
     TRANSFORM  R-THETA-TAU
     MOLTYPE    zxy2
     REFER-CONF RIGID 


The corresponding numerical KEO of H\ :sub:`2`\ CS in the form of the checkpoint file ``kinetic.chk`` can be also found at  `MOTY spectroscopic model <https://exomol.com/models/H2CS/1H2-12C-32S/MOTY/>`__.

Although efficient in interfacing TROVE with external KEO constructors, the latter scheme "External Numerical Taylor-type" has a few disadvantages:

- The 1D basic expansion types must be implemented in TROVE, which reduces its flexibility.
- When representing as sum-of-products, TROVE treats it as a Tylor-type expansion of a given order :math:`N_{\rm kin}` (defined with ``KinOrder``). The number of expansion terms grows very quickly with the varying of the expansion types. In the case of the H\ :sub:`2`\ CS , the formal expansion order of the KEO terms in Eq. :eq:`e-G-H2CS` (``KinOrder``) was 15 (i.e. :math:`n+l+m+o+p+q`). This was necessary to incorporate all the terms in the exact KEO of this molecule using the basic terms from the tables above. With this order, the number of the formal expansion terms in TROVE was 74613 (see card ``Ncoeff`` in ``kineti.chk``) fr each KEO term :math:`G_{\lambda,\lambda'}`  and TROVE needed to allocated all of them, even though the actual number of non-zero terms was much smaller, only 208. That is, it is expensive. It should be noted that once inputed into the TROVE pipeline, with the ``sparse`` representation, TROVE will only work with 208 terms, but initially it would need to assume all 74613. 

Ideally, we would like to have a similar functionality of the "External Numerical Taylor-type", but with more flexibility in constructing KEO and also in the way it is inputed.  In response to these requirements, the ``BASIC-FUNCTION`` feature has been introduced. 





Basic-Function: External sum-of-products KEO of general type
------------------------------------------------------------

This feature allows user to build their KEO from the existing "basic-functions" using the ``BASIC-FUNCTION`` builder. 

``BASIC-FUNCTION`` is a TROVE block defining the shape of the KEO expansions in Eq. :eq:`e-G-H2CS` as follows. Using the same example of H\ :sub:`2`\ CS, we now introduce the basic coordinates using the following constructor: 
::
     
     BASIC-FUNCTION
     Mode 1 2
     1 1 -1 I 1 1
     2 1 -2 I 1 1
     Mode 2 2
     1 1 -1 I 1 1
     2 1 -2 I 1 1
     Mode 3 2
     1 1 -1 I 1 1
     2 1 -2 I 1 1
     Mode 4 7
     1 1 1 Cos 1 1
     2 1 1 Cot 1 1
     3 1 1 Csc 1 1
     4 1 1 Sin 1 1
     5 1 2 Cot 1 1
     6 2 1 Cot 1 1 1 Csc 1 1
     7 1 2 Csc 1 1
     Mode 5 7
     1 1 1 Cos 1 1
     2 1 1 Cot 1 1
     3 1 1 Csc 1 1
     4 1 1 Sin 1 1
     5 1 2 Cot 1 1
     6 2 1 Cot 1 1 1 Csc 1 1
     7 1 2 Csc 1 1
     Mode 6 5
     1 1 1 Cos 0.5 1
     2 1 1 Sin 0.5 1
     3 1 2 Cos 0.5 1
     4 2 1 Cos 0.5 1 1 Sin 0.5 1
     5 1 2 Sin 0.5 1
     END
    
    
We assume that each   1D contributing term :math:`u_{l}{\xi_i}` for a given mode :math:`\xi` in the expansion of Eq. :eq:`e-G-H2CS` can be represented using the following general form:

.. math::   
            
      u_{l}(\xi_i) = f_1( a_i \xi_i^{k_1} )^{n_i} f_2( b_i \xi_i^{l_1} )^{m_i} \ldots
       

For example, for mode 1 (:math:`r_1`), two expansion types are selected: :math:`f_1{(1)}(r_1) = \frac{1}{r_1}` and :math:`f_1{(1)}(r_1) = \frac{1}{r_1^2}`, which are summarised in the following table:

+-------+-----------------------+-------------+-------+------------+-------------+   
|                               |        mode |:math:`N_{\rm func}`|             |
+-------+-----------------------+-------------+-------+------------+-------------+
|       |  ``Mode``             |           1 |       2            |             |
+-------+-----------------------+-------------+-------+------------+-------------+
|   #   |  :math:`N_{\rm contr}`| :math:`n_i` | type  |:math:`a_i` | :math:`k_i` |   
+-------+-----------------------+-------------+-------+------------+-------------+
|   1   |          1            |  -1         |  I    |  1         |     1       | 
+-------+-----------------------+-------------+-------+------------+-------------+
|   2   |          1            |  -2         |  I    |  1         |     1       |
+-------+-----------------------+-------------+-------+------------+-------------+

Here ``I`` indicates the inverse function of :math:`r_1`. 

For mode 4 (:math:`\alpha_1`) however 7 basic functions are used:

.. math:: 
    
    \begin{split}
      f^{(1)} & = \cos(\tau/2) \\
      f^{(2)} & = \sin(\tau/2) \\
      f^{(3)} & = \cos^2(\tau/2) \\
      f^{(4)} & = \cos(\tau/2)\sin(\tau/2) \\
      f^{(4)} & = \sin^2(\tau/2) 
    \end{split}



which is summarised in the following table: 

+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+
|                               |        mode |:math:`N_{\rm func}`|             |                    |                         |
+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+
|       |  ``Mode``             |           4 |     7              |             |             |      |           |             |
+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+
|   #   |  :math:`N_{\rm contr}`| :math:`n_i` | type  |:math:`a_i` | :math:`k_i` | :math:`m_i` | type |:math:`b_i`| :math:`l_i` |
+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+
|   1   |          1            |   1         |  Cos  |  1         |     1       |             |      |           |             |
+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+
|   2   |          1            |   1         |  Cot  |  1         |     1       |             |      |           |             |
+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+
|   3   |          1            |   1         |  Csc  |  1         |     1       |             |      |           |             |
+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+
|   4   |          1            |   1         |  Sin  |  1         |     1       |             |      |           |             |
+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+
|   5   |          1            |   2         |  Cot  |  1         |     1       |             |      |           |             |
+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+
|   6   |          2            |   1         |  Cot  |  1         |     1       |        1    | Csc  |    1      |     1       |
+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+
|   7   |          1            |   2         |  Csc  |  1         |     1       |             |      |           |             |
+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+


And for the last mode 6 (:math:`\tau`), the basic coordinates 

.. math::

    \begin{split}
       f^{(1)} & = \cos\alpha_1 \\
       f^{(2)} & = \cot\alpha_1 \\
       f^{(3)} & = \csc\alpha_1 \\
       f^{(4)} & = \sin\alpha_1 \\
       f^{(5)} & = \cot^2\alpha_1 \\
       f^{(6)} & = \cot\alpha_1\csc\alpha_1 \\
       f^{(7)} & = \csc^2\alpha_1
    \end{split}

which is summarised in the following table:

+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+
|                               |        mode |:math:`N_{\rm func}`|             |                    |                         |
+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+
|       |  ``Mode``             |           6 |     5              |             |             |      |           |             |
+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+
|   #   |  :math:`N_{\rm contr}`| :math:`n_i` | type  |:math:`a_i` | :math:`k_i` | :math:`m_i` | type |:math:`b_i`| :math:`l_i` |
+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+
|   1   |          1            |   1         |  Cos  |  0.5       |     1       |             |      |           |             |
+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+
|   2   |          1            |   1         |  Sin  |  0.5       |     1       |             |      |           |             |
+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+
|   3   |          1            |   1         |  Cos  |  0.5       |     1       |             |      |           |             |
+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+
|   4   |          1            |   2         |  Cos  |  0.5       |     1       |       1     |Sin   |    0.5    |     1       |
+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+
|   5   |          1            |   1         |  Sin  |  0.5       |     1       |             |      |           |             |
+-------+-----------------------+-------------+-------+------------+-------------+-------------+------+-----------+-------------+


Structure of kinetic.chk
------------------------

For the Basic-Function feature,  we use the extended format of the checkpoint file kinetic.chk, see :doc:`checkpoints. It has the following form 
::
    
       0   7   93    <--- Gvib Npoints,Norder,Ncoeff
       1     1     1   0   3.86412653967     0    0    0    0    0    0
       1     2     1   0   2.80960448197     0    0    0    1    0    0
       1     3     1   0   2.80960448197     0    0    0    0    1    0
       1     4     1   0  -2.80960448197     0    1    0    4    0    0
       1     5     1   0  -2.80960448197     0    0    1    0    4    0
       1     6     1   0               0     0    0    0    0    0    0
       2     1     1   0   2.80960448197     0    0    0    1    0    0
       2     2     1   0   36.2630831225     0    0    0    0    0    0
      .....
           

where the first three integer numbers are  :math:`N_{\rm points}`, :math:`N_{\rm order}` and :math:`N_{\rm coeff}` with :math:`N_{\rm coeff}` is the number of non-zero coefficients. The columns of the main body are as follows:
::

     ----- ----- ---- --- ---------------- ---- ---- ---- ---- ----- ----
       1     2     3   4           5         6    7    8    9    10   11
     ----- ----- ---- --- ---------------- ---- ---- ---- ---- ----- ----
       1     1     1   0   3.86412653967     0    0    0    0    0    0
       1     2     1   0   2.80960448197     0    0    0    1    0    0
       1     3     1   0   2.80960448197     0    0    0    0    1    0
       .....
       
The first 5 columns are as before:

- col 1: the value of the index :math:`\lambda` in :math:`G_{\lambda,\lambda'}^{\rm vib}`
- col 2: the value of the index :math:`\lambda'` in :math:`G_{\lambda,\lambda'}^{\rm vib}`
- col 3: expansion index :math:`n`: 

.. math::
       
       G_{\lambda,\lambda'}^{\rm vib}(\xi) = \sum_i g_{\lambda,\lambda'}[n]\, \xi_1^{i_1} \xi_1^{i_2} \ldots  \xi_{i_M}
      

- col 4: The non-rigid grid point :math:`k` (:math:`\alpha_k`), zero for the rigid and :math:`k=0\ldots N_{\rm ponits}`
- col 5: Value of the expansion parameter :math:`g_{\lambda,\lambda'}[n](\alpha_k)`.
- col 6: :math:`i_1`
- col 7: :math:`i_2`
- col 8: :math:`i_3`
- col 9: :math:`i_4`
- col 10: :math:`i_5`
- col 11: :math:`i_6` (i.e. :math:`i_M`).

Special Expansion ``Automatic`` type in the basis set block
-----------------------------------------------------------

Finally, the expansion basic functions now have to be specified in the kinetic part of the ``basis`` block as the ``automatic`` type as in the following example: 
:: 

     BASIS
       0,'JKtau', Jrot 0
       1,'numerov','AUTOMATIC', 'morse',  range 0, 8, resc 1.0, points 1000, borders -0.35,0.9
       2,'numerov','AUTOMATIC', 'morse',  range 0, 4, resc 2.0, points 1000, borders -0.5,0.9
       2,'numerov','AUTOMATIC', 'morse',  range 0, 4, resc 2.0, points 1000, borders -0.5,0.9
       3,'numerov','AUTOMATIC', 'linear', range 0, 8, resc 1.0, points 1000, borders -1.0,1.0
       3,'numerov','AUTOMATIC', 'linear', range 0, 8, resc 1.0, points 1000, borders -1.0,1.0
       4,'numerov','AUTOMATIC', 'linear', range 0, 8, resc 1.0, points 2000, borders -120.,  120.0 deg
     END
 
This will make sure that the correct expansion functions are used when computing KEO matrix elements. 

