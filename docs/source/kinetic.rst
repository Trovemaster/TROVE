Kinetic energy operators
========================

There are several options for kinetic energy operators (KEO) in TROVE: 

- **Linearised KEO**: Automatically constructed KEO using the Eckart frame (rigid) or Eckart/Saywitz frame ((1D non-rigid) as a Taylor expansion in terms of linearised coordinates;
- **Curvilinear KEO**: Preprogrammed analytic KEO (exact or non-exact) using curvilinear (valence) coordinates and different frames;
- **External Numerical Taylor-type KEO**: Externally generated KEO as a Taylor-type numerical expansion in terms of any user-defined coordinates;
- **External sum-of-products KEO**: Externally generated KEO as a generic sum-of-products expansion via the TROVE KEO builder.


Linearised KEO
--------------

This is a standard KEO constructed using a general black-boxed algorithm applicable for an arbitrary polyatomic molecule. It has been introduced in the original TROVE paper [TROVE]_ and also explained in Chapter `Theory <../theory.html>`__. For the following, it useful to recall that TROVE uses Z-matrix curvilinear coordinates to define the linearised coordinates for a given molecule and selected frame. 

Linearised KEO is a default option in TROVE. It requires the following components to be present in the TROVE code:

- Coordinate transformation (``TRANSFORM``) between the Z-matrix curvilinear coordinates :math:`\xi_{i}^{\rm Zmat}` and the target curvilinear coordinates :math:`\xi_i` used to generate the linearised coordinates :math:`\xi_i^{\rm lin}`. For example, for the Z-matrix used for NH\ :sub:`3` (see <../molecules.html>__):
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

where  :math:`\delta` is the angle between the trisector and any of the bond vectors.  The five linearised coordinates :math:`\xi_i^{\rm lin}` for a non-rigid reference configuration around :math:`\delta` are constructed by linearising :math:`\xi_i`. 


- Equilibrium or non-rigid reference structure associated with the frame (``Frame``) and/or coordinate transformation (``TRANSFORM``). For NH\ :sub:`3`, a non-rigid reference frame associated with the coordinates ``r-s-delta`` is  given by:

.. math::

      \left( \begin{array}{cc}
        0 & 0 & -\frac{3 m_H r_{\rm e} \cos\rho}{3 m_H+m_X} \\
         r_{\rm e} \sin\rho& 0 & \frac{m_X r_{\rm e} \cos\rho}{3 m_H+mX} \\
        -\frac{1}{2} r_{\rm e} \sin\rho&  \frac{\sqrt{3}}{2} r_{\rm e} \sin\rh & \frac{m_X r_{\rm e} \cos\rho}{3 m_H+mX} \\
        -\frac{1}{2} r_{\rm e} \sin\rho& -\frac{\sqrt{3}}{2} r_{\rm e} \sin\rh & \frac{m_X r_{\rm e} \cos\rho}{3 m_H+mX} \\
      \end{array} 
      \right)



