Kinetic energy operators
========================

There are several options for kinetic energy operators (KEO) in TROVE: 

- **Linearised KEO**: Automatically constructed KEO using the Eckart frame (rigid) or Eckart/Saywitz frame ((1D non-rigid) as a Taylor expansion in terms of linearised coordinates;
- **Curvilinear KEO**: Preprogrammed analytic KEO (exact or non-exact) using curvilinear (valence) coordinates and different frames;
- **External Numerical Taylor-type KEO**: Externally generated KEO as a Taylor-type numerical expansion in terms of any user-defined coordinates;
- **External sum-of-products KEO**: Externally generated KEO as a generic sum-of-products expansion via the TROVE KEO builder.


Linearised KEO
--------------

This is a standard KEO constructed using a general black-boxed algorithm applicable for an arbitrary polyatomic molecule. It has been introduced in the original TROVE paper [TROVE]_ and also explained in Chapter `Theory <https://spectrove.readthedocs.io/en/latest/theory.html>`__. For the following, it useful to recall that TROVE uses Z-matrix curvilinear coordinates to define the linearised coordinates for a given molecule and selected frame.

Linearised KEO is a default option in TROVE. It requires the following components to be present in the TROVE code:

- Coordinate transformation (``TRANSFORM``) between the Z-matrix curvilinear coordinates :math:`\xi_{i}^{\rm Zmat}` and the target curvilinear coordinates :math:`\xi_i` used to generate the linearised coordinates :math:`\xi_i^{\rm lin}`. For example, for the Z-matrix used for NH\ :sub:`3` (see Chapter `Molecules <https://spectrove.readthedocs.io/en/latest/molecules.html>`__):
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


Each term in the linearised KEO is represented by a Taylor expansion in terms of :math:`\xi_{i}^{\rm lin}` (or some 1D functions of them) around the (non-) rigid reference configuration, e.g.
.. math::

      G_{n,n'} = \sum_{i,j,k,l,\ldots } G_{i,j,k,l,\ldots}^{(n,n')} \xi_1^i \xi_2^j \xi_3^k \xi_4^l \ldots

where the subscript 'lin' is omitted for compactness. This is a sum-of-products form which simplifies the integration with the basis functions.


Analytic curvilinear KEO
------------------------

A number of analytic curvilinear KEOs are available (see Chapter `Frames <https://spectrove.readthedocs.io/en/latest/frames.html>`__) and associated with the card ``kinetic_type`` as part of the ``KINETIC`` section, e.g.:
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

#. Define the KEO type using the ``kinetic`` block:
::

    KINETIC
      kinetic_type  KINETIC_XYZ_EKE_bisect
    END
    
#. Set the ``Coords`` type to ``local``, KEO expansion to 2,  choose the appropriate frame/and coordinates and molecule type, set the reference configuration to ``non-rigid``: 
::
    
    KinOrder  2    
    COORDS local (curvilinear)
    TRANSFORM   R-RHO-Z-M2-M3-BISECT (FRAME)
    MOLTYPE XY2
    REFER-CONF NON-RIGID
    
#. Choose appropriate basis set configuration, i.e. the KEO expansion type (e.g. ``rational``) and the non-rigid basis set type (e.g. ``laguerre-k``):
:: 

     BASIS
       0,'JKtau', Jrot 0, krot  18
       1,'numerov','rational', 'morse',  range 0,20, r 8, resc 3.0, points   2000, borders -0.3,0.90
       2,'numerov','rational', 'morse',  range 0,36, r 8, resc 1.5, points   3000, borders -0.4,0.90
       3,'laguerre-k','linear','linear', range 0,48, r 8, resc 1.0, points  12000, borders  0.,120.0 deg
     END



External Numerical Taylor-type KEO
----------------------------------

As was mentioned above, TROVE can work with any ro-vibrational KEOs as long as they are represented as sum-of-products of 1D functions of vibrational coordinates. It is possible to input any externally constructed (non-singular) KEO and to be used with the TROVE pipe line. One of the most robust methods is to use the TROVE checkpoints functionality.  


