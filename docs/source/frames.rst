Frames and vibrational coordinates 
**********************************


Here we introduce different ingredients available for triatomic molecules, including

- Molecular frames :math:`xyz`;
- :math:`3N-6` (:math:`3N-5`) vibrational coordinates :math:`\xi_n`;

For the linearised coordinates, the default frame is Eckart. The equilibrium structures, required for the definition of the linearised KEO and PEF, are chosen as the principal axis system (PAS). For the curvilinear KEOs, the frames are defined by the construction of the KEOs in their analytic representations.

Triatomics
==========

XY\ :sub:`2` type molecules
---------------------------


in the curvilinear KEO,  it is common in TROVE to use the bisector frame for the XY\ :sub:`2` molecules, with the :math:`x` axis bisecting the bond angle and the :math:`z` in the plane of the molecule, but other embeddings are possible. The PAS frame coincides with the bisector frame at the equilibrium or non-rigid reference configuration (i.e. symmetric).  In TROVE, the definition of the frame is combined with the definition of the internal coordinates via the keyword ``transform``. In the following, these are described.

There are currently at least two  exact, curvilinear KEO forms are provided for a quasi-linear XY\ :sub:`2` molecules, ``MLkinetic_xy2_bisect_EKE``, ``MLkinetic_xy2_bisect_EKE_sinrho``, see below. 


``R-RHO-Z``
^^^^^^^^^^^

- ``R-RHO-Z`` is used for (quasi-)linear molecules of the XY\ :sub:`2` type. it defined the curvilinear vibrational coordinates as the two bond angles :math:`r_1` and :math:`r_2` with  the bending mode described by the angle :math:`\rho = \pi - \alpha`, where :math:`\alpha` is the interbond angle (:math:`\rho = 0 \ldots \rho_{\rm max}`). For the rigid reference frame (``REFER-CONF RIGID``), the actual internal coordinates are the displacements of :math:`r_1`, :math:`r_2` and :math:`\rho` from the corresponding equilibrium:

.. math::

    \xi_1 = r_1 - r_{\rm e}, \\
    \xi_2 = r_2 - r_{\rm e}, \\
    \xi_3 = \rho,

where :math:`r_{\rm e}` is the equilibrium bond length. If the non-rigid reference frame is used (``REFER-CONF NON-RIGID``), the bending mode is given on an equidistant grid, typically of 1000-2000 points, while the stretching modes are the displacements from the given :math:`\rho` point along the non-rigid reference frame, the latter is usually defined as the principal axes system with the bond length fixed to the equilibrium.

.. math::

    \xi_1 = r_1 - r_{\rm ref}, \\
    \xi_2 = r_2 - r_{\rm ref}, \\
    \xi_3 = \rho,

Alternatively, the reference value of the bond length :math:`r_{\rm ref}` can also vary with :math:`\rho` as e.g. in the minimum energy path (MEP) definition with :math:`r_{\rm ref}` being the optimised value at the given value of :math:`\rho` corresponding to the local energy minimum. In this case, the non-rigid frame must be defined using the ``MEP`` block (see the corresponding section).

For the linearised coordinates type (``COORDS Linear``), the actual internal coordinates are the linearised versions of :math:`\xi_i` above. More specifically, for the ``non-rigid`` reference configuration, the bending coordinate :math`\rho` is kept curvilinear on a grid of :math`\rho_k` points as before, while the stretching coordinates are defined by linearly expanding :math`r_1` and :math`r_2` in terms of the Cartesian displacement around the corresponding reference values :math:`r_{\rm ref}`. In the ``Rigid`` case, the bending coordinate is also linearised.

The advantage of the linearised coordinates is that the corresponding KEO can be constructed on the fly as part of the TROVE generalised procedure as a Taylor type expansion. The main disadvantage however is that the approximate linearised KEO operator is less accurate than the (exact) curvilinear EKO. Besides, the convergence of the variational solution is also poorer for the linearised case (see [15YaYu]_).



``R-RHO-Z-ECKART``
^^^^^^^^^^^^^^^^^^

This ``Transform`` type is very similar to ``R-RHO-Z``, but with the molecular frame define using the Eckart conditions.


``R-ALPHA-Z``
^^^^^^^^^^^

- ``R-ALPHA-Z`` is very similar to ``R-RHO-Z`` with the difference in the bending coordinate, which in the interbond angle :math:`\alpha` in this case. In the ```Rigid` reference configuration, it is a displacement from the equilibrium value :math:`\alpha_{\rm e}`:
.. math::

    \begin{split}
    \xi_1 &= r_1 - r_{\rm e}, \\
    \xi_2 &= r_2 - r_{\rm e},\\ 
    \xi_3 &= \alpha-\alpha_{\rm e}.
    \end{split}

In the ``Non-rigid`` reference configuration, :math:`\alpha` is given on a grid of points ranging from :math:`\alpha_{\rm min}` to :math:`\alpha_{\rm max}` and including the equilibrium value. In the linearised ``Rigid`` case, the bending coordinated is defined as a linear expansion of :math:`\alpha` at :math:`\alpha_{\rm eq}`  in terms of the Cartesian displacements.


TROVE input example:
::

COORDS       local    (curvilinear coordinates)
TRANSFORM    r-rho-z  (r1, r2, rho with the x parallel to the bisector)
MOLTYPE      XY2
REFER-CONF   non-RIGID  (Reference configuration)

.. Note:: The text in brackets is used for comments.


XYZ type molecules
------------------

The main embedding here is the 'bond'-embedding, with the :math:`z` axis placed parallel to the bond Y-Z with a heavier atom Z comparing to X (second bond).
For molecules XYZ with  comparable masses X and Z (e.g. in similar isotopologues), the bisector frames and associated ``TRANSFORM`` can be used.



``R1-Z-R2-RHO``
^^^^^^^^^^^^^^^^^

This is a 'bond'-embedding with the same vibrational coordinates as in ``R-RHO-Z``.


``R1-Z-R2-ALPHA``
^^^^^^^^^^^^^^^^^

This is another 'bond'-embedding with the same vibrational coordinates as in ``R-ALPHA-Z``.


Tetratomics
===========

XY\ :sub:`3` rigid  molecules (PH\ :sub:`3` type)
-------------------------------------------------

Linearized KEOs use the Eckart frame with the PAS at the equilibrium configuration. The latter has the :math:`z` axis along the axis of symmetry :math:`C_3` with the :math:`x` axis chosen in plane containing the X-Y\ :sub:`1` bond and passing through :math:`C_3`. 


``R-ALPHA``
^^^^^^^^^^^

For the rigid XY\ :sub:`3`, like PH\ :sub:`3`, the logical coordinate choice of the valence coordinates consists of three bond lengths :math:`r_1`, :math:`r_2`, :math:`r_3`, :math:`\alpha_{23}`, :math:`\alpha_{13}` and :math:`\alpha_{12}`. For the linearised KEO, these valence are used to form the linearised coordinates in the same way as before (1st order expansion in terms of the Cartesian displacement). For the curvilinear KEO (``local``), the vibrational coordinates are then defined as displacement from the corresponding equilibrium (or non-rigid reference) values:

.. math::
    \begin{split}
    \xi_1 &= r_1 - r_{\rm e}, \\
    \xi_2 &= r_2 - r_{\rm e}, \\
    \xi_3 &= r_3 - r_{\rm e}, \\
    \xi_4 &= \alpha_{23}-\alpha_{\rm e}, \\
    \xi_5 &= \alpha_{13}-\alpha_{\rm e}, \\
    \xi_6 &= \alpha_{12}-\alpha_{\rm e}.
    \end{split}

.. sidebar::

    .. figure:: img/PH3.jpg
       :alt: PH3 equilibrium structure 
     
       PH\ :sub:`3` equilibrium structure 


This representation has been used for PH\ :sub:`3` [15SoAlTe]_, SbH\ :sub:`3` [10YuCaYa]_, AsH\ :sub:`3` [19CoYuKo]_, PF\ :sub:`3` [19MaChYa]_.


XY\ :sub:`3` non-rigid with umbrella motion (NH\ :sub:`3` type)
---------------------------------------------------------------

Consider the Ammonia molecule NH3\ :sub:`3` with a relatively small barrier to the planarity. The three bending angles are not suitable in this case  as they cannot distinguish the two opposite inversion configurations above and below the planarity. Instead, an umbrella mode has to be introduced as one of the bending modes. An example of an umbrella coordinate is an angle between the :math:`C_3` symmetry axis and the bond X-Y, see Figure. It is natural to use the non-rigid reference configuration along the umbrella, inversion motion and build the KEO as an expansion around it. For two other bending modes, in principle one can use two inter-bond angles, e.g.  :math:`\alpha_2` and :math:`\alpha_3`, two dihedral angles :math:`\phi_2` and :math:`\phi_3`. However, for symmetry reasons, TROVE employs the symmetry-adapted bending pair :math:`S_a` and :math:`S_b`, defined as follows:

.. math::

    S_a = \frac{1}{\sqrt{6}} (2 \alpha_{23}-\alpha_{13}-\alpha_{12}), \\
    S_b  = \frac{1}{\sqrt{2}} ( \alpha_{13}-\alpha_{12})
    

or 


.. math::

    S_a = \frac{1}{\sqrt{6}} (2 \phi_{23}-\phi_{13}-\phi_{12}), \\
    S_b  = \frac{1}{\sqrt{2}} ( \phi_{13}-\phi_{12})


The umbrella mode for any instantaneous configuration of the nuclei is defined in TROVE as the angle between a trisector 



.. sidebar::

   .. figure:: img/umbrella.jpg
       :alt: Umbrella motions

       NH\ :sub:`3`: umbrella modes :math:`\rho` and :math:`\delta`. 



Linearized KEOs use the Eckart frame with the PAS at the equilibrium configuration. The latter has the :math:`z` axis along the axis of symmetry :math:`C_3` with the :math:`x` axis chosen in plane containing the X-Y\ :sub:`1` bond and passing through :math:`C_3`.



``R-S-DELTA``
^^^^^^^^^^^^^

For this ``TRANSFORM`` case, the following valence-based coordinates are used: 


.. math::

    \begin{split}
    \xi_1 &= r_1 - r_{\rm e}, \\
    \xi_2 &= r_2 - r_{\rm e}, \\
    \xi_3 &= r_3 - r_{\rm e}, \\
    \xi_4 &= \frac{1}{\sqrt{6}} (2 \alpha_{23}-\alpha_{13}-\alpha_{12}),  \\
    \xi_5 &= \frac{1}{\sqrt{2}} ( \alpha_{13}-\alpha_{12}),  \\
    \xi_6 &= \delta. 
    \end{split}

The umbrella mode :math:``\delta`` is defined as an angle between the trisector and any of the bonds X-Y. The other 5 coordinates are then used to construct the corresponding linearised vibrational coordinates (see above) for the linearised (``linear``) representation. 



ZXY\ :sub:`2` (Formaldehyde type)
---------------------------------

.. sidebar::

   .. figure:: img/H2CO.jpg
       :alt: H2CO

       Valence coordinates and the bisector frame used for H\ :sub:`2`\ CO.




The common valence coordinate choice for ZXY\ :sub:`2` includes three bond lengths , two bond angles and a dihedral angle :math:`\tau`. The latter can be treated as the reference for a non-rigid reference configuration in TROVE on a grid of :math:`\tau_i` ranging from  :math`[-\tau_{0}\ldots \tau_{0}]`, while other 5 modes are treated as displacement from their equilibrium values at each grid point :math:`\tau_i`. Apart from the standard linearised KEO, a curvilinear exact KEO has been recently introduced into TROVE. This is exactly the ``R-THETA-TAU`` type, detailed below.


``R-THETA-TAU``
^^^^^^^^^^^^^^^

.. math::

    \begin{split}
    \xi_1 &= r_1 - r_{\rm e}, \\
    \xi_2 &= r_2 - r_{\rm e}, \\
    \xi_3 &= r_3 - r_{\rm e}, \\
    \xi_4 &= \theta_1,  \\
    \xi_5 &= \theta_2,  \\
    \xi_6 &= \tau.
    \end{split}





