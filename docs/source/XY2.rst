Triatomics
==========

Here we introduce different ingredients available for triatomic molecules, including

- Molecular frames :math:`xyz`;
- :math:`3N-6` (:math:`3N-5`) vibrational coordinates :math:`\xi_n`;
- Kinetic energy operators (KEO);
- Potential energy functions (PEF);
- For intensity calculations, 3D dipole moment functions.

XY\ :sub:`2` type molecules
---------------------------

The main frame is the 'bisector', with the :math:`x` axis bisecting the bond angle and the :math:`z` in the plane of the molecule, but other embeddings are possible. In TROVE, the definition of the frame is combine the definition of the internal coordinates via the keyword ``transform``. In the following, these are described.


``R-RHO-Z``
^^^^^^^^^^^

- ``R-RHO-Z`` is used for (quasi-)linear molecules of the XY\ :sub:`2` type. it defined the curvilinear vibrational coordinates as the two bond angles :math:`r_1` and :math:`r_2` with  the bending mode described by the angle :math:`\rho = \pi - \alpha`, where :math:`\alpha` is the interbond angle (:math:`\rho = 0 \ldots \rho_{\rm max}`). For the rigid reference frame (``REFER-CONF RIGID``), the actual internal coordinates are the displacements of :math:`r_1`, :math:`r_2` and :math:`\rho` from the corresponding equilibrium:
.. math::

    \xi_1 = r_1 - r_{\rm e},
    \xi_2 = r_2 - r_{\rm e},
    \xi_3 = \rho,
     
where :math:`r_{\rm e}` is the equilibrium bond length. If the non-rigid reference frame is used (``REFER-CONF NON-RIGID``), the bending mode is given on an equidistant grid, typically of 1000-2000 points, while the stretching modes are the displacements from the given :math:`\rho` point along the non-rigid reference frame, the latter is usually defined as the principal axes system with the bond length fixed to the equilibrium. 
.. math::

    \xi_1 = r_1 - r_{\rm ref},
    \xi_2 = r_2 - r_{\rm ref},
    \xi_3 = \rho,

Alternatively, the reference value of the bond length :math:`r_{\rm ref}` can also vary with :math:`\rho` as e.g. in the minimum energy path (MEP) definition with :math:`r_{\rm ref}` being the optimised value at the given value of :math:`\rho` corresponding to the local energy minimum. In this case, the non-rigid frame must be defined using the ``MEP`` block (see the corresponding section). 

For the linearised coordinates type (``COORDS Linear``), the actual internal coordinates are the linearised versions of :math:`\xi_i` above. More specifically, for the ``non-rigid`` reference configuration, the bending coordinate :math`\rho` is kept curvilinear on a grid of :math`\rho_k` points as before, while the stretching coordinates are defined by linearly expanding :math`r_1` and :math`r_2` in terms of the Cartesian displacement around the corresponding reference values :math:`r_{\rm ref}`. In the ``Rigid`` case, the bending coordinate is also linearised. 

The advantage of the linearised coordinates is that the corresponding KEO can be constructed on the fly as part of the TROVE generalised procedure as a Taylor type expansion. The main disadvantage however is that the approximate linearised KEO operator is less accurate than the (exact) curvilinear EKO. Besides, the convergence of the variational solution is also poorer for the linearised case (see [15YaYu]_). 


'R-RHO-Z-ECKART'
^^^^^^^^^^^^^^^^

This ``Transform`` type is very similar to ``R-RHO-Z``, but with the molecular frame define using the Eckart conditions. 


``R-ALPHA-Z``
^^^^^^^^^^^

- ``R-ALPHA-Z`` is very similar to ``R-RHO-Z`` with the difference in the bending coordinate, which in the interbond angle :math:`\alpha` in this case. In the ```Rigid` reference configuration, it is a displacement from the equilibrium value :math:`\alpha_{\rm e}`:
.. math::

    \xi_1 = r_1 - r_{\rm ref},
    \xi_2 = r_2 - r_{\rm ref},
    \xi_3 = \alpha-\alpha_{\rm e}. 

In the ``Non-rigid`` reference configuration, :math:`\alpha` is given on a grid of points ranging from :math:`\alpha_{\rm min}` to :math:`\alpha_{\rm max}` and including the equilibrium value. In the linearised ``Rigid`` case, the bending coordinated is defined as a linear expansion of :math:`\alpha` at :math:`\alpha_{\rm eq}`  in terms of the Cartesian displacements. 





