Dipole moment functions
=======================


TROVE provides a larger number of dipole moment functions (DMFs) for different molecules already implemented. Most of these PEFs are in modules ``pot_*`` contained in  file ``pot_*.f90``.

 - ``pot_xy2.f90``
 - ``pot_xy3.f90``
 - ``pot_zxy2.f90``
 - ``pot_abcd.f90``
 - ``pot_xy4.f90``
 - ``pot_zxy3.f90``
 - ........

These are a part of the standard TROVE compilation set. Alternatively, a user-defined DMF can be included into the TROVE compilation as a generic 'user-defined' module ``pot_user``.



Dipole (External)  Block
------------------------

The DMFs are defined in the TROVE input file using the ``Dipole`` block, which is just an alias for the ``external`` input structure. A typical ``Dipole`` input is as follows:
::

      DIPOLE
      dimension 3
      NPARAM  264 0 0
      compact
      DMS_TYPE  XY3_SYMMB
      COEFF   list
      COORDS  linear linear linear linear linear linear
      Order   6  6  6
      dstep 0.005
      parameters
      charge              0.00000000
      nparamA           112.00000000
      RhoE               90.00000000
      RE14                1.01032000
      beta                1.00000000
      d0                  4.56621083
      f1a                -9.36438932
      f2a                32.96400671
      ...
      end

For an example, see :download:`14N-1H3__BYTe__TROVE__step1.inp <./input/14N-1H3__BYTe__TROVE__step1.inp>`  where this DMF is used.

- ``dimension`` (aliases ``rank``, ``dim``): This is the dimension of the external field.   An "external" is treated in TROVE as a vector of dimension :math:`D`, which in the case of dipole can be up to :math:`D=3`, but will depend on the implementation. This parameter is to help structure the input dipole parameters according to the dipole components, if necessary.
- ``NPARAM`` is used to specify the number of parameters to define the DMF and should contain :math:`D` input values.
- ``compact``: is recently implemented card which switches to the "compact" format with no "fitting-indexes" column  present.
- ``DMS_TYPE`` (``TYPE``) is the name of the DMF as implemented in ``pot_*.f90 file`` and referenced in ``mlecules.f90``.
- ``COEFF`` indicates if the DFM parameters are given as a list of parameter values (``LIST``) or  values with the corresponding expansion powers (``POWERS``), see example below.
- ``COORDS``: coordinate types used to re-expand the dipole field in terms of the internal TROVE coordinates.
- ``Order``: The corresponding expansion order.
- ``dstep``: finite difference step used in the re-expansion. The default value is 0.005 Ang.
- ``parameters``: card indicating the section with the dipole parameters entries specific for the given ``DMS_TYPE``.

The dipole moment parameters are listed after the keyword ``parameters`` and terminated with the keyword ``END``. The number of the entries should be equal exactly to the sum of  ``NPARAM`` values.


For the ``COEFF  list`` option, the meaning of the columns is as follows:


+---------+-----------------------+
| Label   | Value                 |
+---------+-----------------------+
| charge  |        0.00000000     |
+---------+-----------------------+
| nparamA |       112.00000000    |
+---------+-----------------------+
| RhoE    |     90.00000000       |
+---------+-----------------------+
| RE14    |      1.01032000       |
+---------+-----------------------+
| beta    |      1.00000000       |
+---------+-----------------------+
| d0      |      4.56621083       |
+---------+-----------------------+
| f1a     |     -9.36438932       |
+---------+-----------------------+
| f2a     |     32.96400671       |
+---------+-----------------------+




The first column with the name of the parameters, which is for clearness only. This field is only used for printing purposes and otherwise nor-referenced in the code in any way. The second column contains the actual value of the given parameter. The input is directly associated with the corresponding implementation and therefore the order is important.

An alternative, legacy, format with no ``compact`` card assumes an additional column with the so-called "fitting-indexes" indicating if the parameter was varied in the fitting to the ab initio data. Here is an example:
::

    parameters
    charge       0           0.00000000
    nparamA      1          112.00000000
    RhoE         0         90.00000000
    RE14         0          1.01032000
    beta         0          1.00000000
    d0           4          4.56621083
    f1a          8         -9.36438932
    f2a          8         32.96400671
    f3a          7        -80.82339377
    ....

where column 2 contains the "fitting-indexes". These indexes are not used by TROVE. They are kept in order to simplify the interfacing between the ab initio fitting and TROVE, but can be always omitted with the help of the card ``compact``.


Here is an example of the input format using individual expansion "powers", ``COEFF powers`` (from CO\ :sub:`2`):
::

         DIPOLE
         rank 3
         NPARAM  971 0 0
         compact
         TYPE  DIPOLE_AMES1
         COEFF   powers  (powers or list)
         COORDS  linear linear linear
         Orders   16 16 16
         threshold   1e-8
         Parameters
         re    0 0 0 0  1.15958d0
         ae    0 0 0 0  180.00
         d000    0    0    0   0      -0.4801402388843266D+00
         d001    0    0    1   0       0.1203598337496481D+00
         d002    0    0    2   0      -0.5662267278952241D-01
         d003    0    0    3   0      -0.2529381009630170D-01
         d004    0    0    4   0      -0.1271678002798687D+00
         d005    0    0    5   0       0.3033049145401118D+01
         d006    0    0    6   0      -0.1754036600894653D+02
     .....
     end

See the TROVE input :download:`CO2_bisect_xyz_step1.inp <./input/CO2_bisect_xyz_step1.inp>`.

Assuming the DMF form as an expansion

.. math::
        
   \mu_\alpha(\xi_1,\xi_2,\xi_3) =  \sum_{k,j,k} d_{i,j,k} \xi_1^i \xi_2^j  \xi_3^k,
    
    
the input card has the following format

   +---------+-----------+----------+----------+-----+-------------------------+
   | Label   | :math:`i` | :math:`j`| :math:`k`|Index| Value :math:`_{i,j,k}`  |
   +---------+-----------+----------+----------+-----+-------------------------+
   |  d000   |          0|    0     |    0     |    0|  -0.4801402388843266D+00|
   +---------+-----------+----------+----------+-----+-------------------------+
   |  d001   |          1|    0     |    0     |    0|  0.1203598337496481D+00 |
   +---------+-----------+----------+----------+-----+-------------------------+
   |  d002   |          0|    1     |    0     |    0|  -0.5662267278952241D-01|
   +---------+-----------+----------+----------+-----+-------------------------+

where

 - 'Labels' are the parameter name,  for printing purposes only;
 - :math:`i`, :math:`j`, :math:`k` are the 'powers' of an expansion term;
  - 'Index' is a switch to indicate if the corresponding parameter was fitted or can be fitted, with no impact on any evaluations of the PEF values. It is not present in the ``compact`` form.
  - 'Values' are the actual dipole parameters. For ``powers``, their order is not important.

In case the definition of DMF requires also structural parameters, such as equilibrium bond lengths :math:`r_{\rm e}`, equilibrium inter-bond angles :math:`\alpha_{\rm e}`,  in the ``COEFF  Powers`` form these parameters should be listed exactly in the order expected by the  implemented of the PEF (similar to the ``COEFF LIST`` form), but with dummy "powers" columns so that their 'values' appear in the right column, as in the example above, ``re`` and ``ae`` are two the equilibrium values and the three columns with ``0 0 0`` are given in order to parse their values using exactly column 6.



Implemented DMFs
================


XY\ :sub:`2` type
-----------------

There are several PEFs available for this molecule type.


``xy2_pq_coeff``
^^^^^^^^^^^^^^^^

.. sidebar::

    .. figure:: img/XY2_dm-pq.jpg
       :alt: PQ frame 

       The :math:`pq` bisector frame ``xy2_pq_coeff`` used for for XY\ :sub:`2`.



This is a bisector-frame DMF, given by two components, :math:`\mu^{(q)}` and :math:`\mu^{(p)}` with the :math:`q` axis being the bisector. The following expansions in terms of the coordinate displacements :math:`\Delta r_1 = r_{\rm 1} - r_{\rm e}`, :math:`\Delta r_2 = r_2 - r_{\rm e}`, and :math:`\cos\rho_{\rm e} - \cos\bar\rho`, where :math:`\bar\rho = \pi - \theta` are used, with :math:`\theta` is the bond angle, and :math:`r_1` and :math:`r_2` are the bond lengths:

.. math::
       :label: e-muQ-1
       
      \begin{split}
        \mu^{(q)} (\Delta r_1, \Delta r_2, \Delta \alpha ) &=  \sin\alpha \left[ \mu_0^{(q)}(\alpha) + \sum_{j} \mu_{j}^{(q)}(\alpha)  \Delta r_j + \sum_{j\le k} \mu_{jk}^{(q)}(\alpha)   \Delta r_j \Delta r_k \right.   \\
        &  \left . + \sum_{j\le k \le m} \mu_{jkm}^{(q)}(\alpha) \Delta r_j \Delta r_k \Delta r_m  + \sum_{j\le k \le m \le n} \mu_{jkmn}^{(q)}(\alpha)  \Delta r_j \Delta r_k \Delta r_m  \Delta r_n  + \ldots \right], \\
        \mu^{(p)} (\Delta r_1, \Delta r_2, \Delta \alpha ) &=  \mu_0^{(p)}(\alpha) + \sum_{j}^{(p)} \mu_{j}^{(p)} (\alpha) \Delta r_j   + \sum_{j\le k}  \mu_{jk}^{(p)}(\alpha) \Delta r_j \Delta r_k   \\
        &    + \sum_{j\le k \le m} \mu_{jkm}^{(p)}(\alpha)  \Delta r_j \Delta r_k \Delta r_m  + \sum_{j\le k \le m \le n} \mu_{jkmn}^{(p)}(\alpha) \Delta r_j \Delta r_k \Delta r_m  \Delta r_n  + \ldots ,
      \end{split}


where all indices :math:`j, k, m`, and :math:`n` assume the values 1 or 2,

.. math::
       :label: e-muQ-exp
       
       \begin{split}
         \mu_{jk\ldots}^{(q)}(\alpha)  =&  \sum_{i=0}^{N} Q_{ij\ldots}^{(i)} (\cos\alpha_{\rm e} - \cos\alpha )^i, \\
         \mu_{jk\ldots}^{(p)}(\alpha)  =&  \sum_{i=0}^{N} P_{ij\ldots}^{(i)} (\cos\alpha_{\rm e} - \cos\alpha )^i,
       \end{split}
        


       
and the :math:`Q_{ij\ldots}^{(i)}` and :math:`P_{ij\ldots}^{(i)}` are molecular dipole parameters. The expansion coefficients in Eqs. :eq:`e-muQ-exp` are subject to the conditions that the functions :math:`\mu^{(q)}` are unchanged under the interchange of the identical protons, whereas the function :math:`\mu^{(p)}` is antisymmetric under this operation. There are 72 and 99 paramters :math:`Q_{ij\ldots}^{(i)}` and :math:`P_{ij\ldots}^{(i)}`, respectively. An example of ``xy2_pq_coeff`` is illustrated above and can be foound in :download:`H2S_EKE_basic-functions_step1.inp <./input/H2S_EKE_basic-functions_step1.inp`.

The implementation can be found in :code:`subroutine MLdms2pqr_xy2` from the module pot_xy2.f90. The transformation between the TROVE frame and the frame of the specifc dipole of the XY\ :sub:`2` is perfomed in the :code:`subroutine MLloc2pqr_xy2`, e.g.:

.. code:: c
    
    !
    select case(trim(molec%frame))
       !
    case('R-RHO-Z','R-RHO-Z-M2-M3','R-RHO-Z-M2-M3-BISECT','BISECT-Z')
       !
       a0(2, 1) = -r(1) * cos(alpha_2)
       a0(2, 3) = -r(1) * sin(alpha_2)
       !
       a0(3, 1) = -r(2) * cos(alpha_2)
       a0(3, 3) =  r(2) * sin(alpha_2)
    case ...
    


``XY2_PQ_LINEAR``
^^^^^^^^^^^^^^^^^

This is similar to ``xy2_pq_coeff``, but with the bending expansion in Eq. :eq:` e-muQ-exp` in terms of the displacement :math:`\alpha-\alpha_{\rm e}`: 

.. math::
       :label: e-muQ-exp-2

       \begin{split}
         \mu_{jk\ldots}^{(q)}(\alpha)  =&  \sum_{i=0}^{N} Q_{ij\ldots}^{(i)} (\alpha - \alpha_{\rm e} )^i, \\
         \mu_{jk\ldots}^{(p)}(\alpha)  =&  \sum_{i=0}^{N} P_{ij\ldots}^{(i)} (\alpha - \alpha_{\rm e} )^i,
       \end{split}





``DIPOLE_AMES1``
^^^^^^^^^^^^^^^^

This DMF is of the AMES1 type represented using the point-charge molecular bond frame [14HuScLe]_ given by projections on the molecular bond vectors :math:`\vec{r}_1` and :math:`\vec{r}_2`:

.. math::

        \vec{\mu} = \mu_x \vec{i} + \mu_z \vec{k} =  \mu_1 \vec{r}_1 + \mu_2 \vec{r}_2


where :math:`\mu_x` and :math:`\mu_z` are the TROVE frame vectors and :math:`\mu_1` and :math:`\mu_2`  are the *ab initio* dipoles in the molecular bond frame (the :math:`\mu_y` component is always zero). The two point-charge dipole moment components  :math:`\mu_1` and :math:`\mu_2` are represented in terms of the vibrational coordinates as

.. math::
       :label: eq:coords_dms
       
       \begin{split}

         \zeta_1 &= r_1-r^{\rm ref}_1, \\
         \zeta_2 &= r_2-r^{\rm ref}_2, \\
         \zeta_3 &= \cos\alpha-1.
       \end{split}
       
       
with the following analytic Taylor-type expansions used  (see e.g. [14HuScLe]_):

.. math::
       
       \begin{split}
        \mu_1 &=  \sum_{ijk} F^{(1)}_{ijk} \zeta_1^{i} \zeta_2^{j} \zeta_3^{k} , \\
        \mu_2 &=  \sum_{ijk} F^{(2)}_{ijk} \zeta_1^{j} \zeta_2^{i} \zeta_3^{k} ,
       \end{split}


As an example can be found of a system where this form was used, see :download:`CO2_bisect_xyz_step1.inp <./input/CO2_bisect_xyz_step1.inp>`.


``DIPOLE_SO2_AMES1``
^^^^^^^^^^^^^^^^^^^^

This form is essentially the same as ``DIPOLE_AMES1`` but some specific characteristic used for the SO\ :sub:`2` molecule in [14HuScLe]_. 


``XY2_C3_SCHROEDER``
--------------------

.. sidebar::

    .. figure:: img/XY2_rot_from_Eckart.jpg
       :alt: Eckart frame

       Transfromation from tje Eckart to any other frames for XY\ :sub:`2` in the :math:`xz` plane.


This DMF is based on the DMF form reported by Schroeder et al. [16ScSe]_ for C\ :sub:`3`. This DMF is in the Ecakrt frmame expressed in terms of two in-plane components, :math:`\mu^{\parallel}` and :math:`\mu^{\perpl}`,  as Taylor expansions around the equilibrium geometry:

.. math::

     \begin{split}
       \mu^{\parallel} &= D^{\parallel}_{ijk} \Delta r_{1}^i \,\Delta r_{2}^j \, \Delta \alpha^k ,  \\
       \mu^{\perp} &= D^{\perp}_{ijk} \Delta r_{1}^i \,\Delta r_{2^j \, \Delta \alpha^k .
     \end{split}


Since TROVE's frame is usually different from the DMF frame (e.g. bisector) in the  ro-vibrational calculations, this dipole moments functions needs to be rotated.  This is done using the rotation angle :math:`\phi`  from an equilibrium bysector frame :math:`r_{n\alpha}^0` to the instantaneous frame  :math:`r_{n\alpha}` (:math:`n=1,2,3` and :math:`\alpha=x,y,z`) in the in the :math:`xz` plane  as given by 

.. math::
    
    \tan\phi = \frac{m_{\rm X} ( r^0_{1,z} r_{1,x}-r^0_{1,x} r_{1,z} )+m_{{\rm Y}_1} (r^0_{2,z}r_{2,x}-r^0_{2,x}r_{2,z} )+m_{{\rm Y}_2}(r^0_{3,z}r_{3,x}-r^0_{3,x}r_{3,z}) }{ m_{\rm X} (r^0_{1,z} r_{1,z}+r^0_{1,x}r_{1,x})+m_{{\rm Y}_1}(r^0_{2,z}r_{2,z}+r^0_{2,x}r_{2,x})+m_{{\rm Y}_2}(r^0_{3,z}r_{3,z}+r^0_{3,x}r_{3,x}) }




``DIPOLE_PQR_XYZ_Z-FRAME``
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. sidebar::

    .. figure:: img/XYZ-pqy_DM.jpg
       :alt: PQY frame

       The :math:`pqy` bisector frame ``DIPOLE_PQR_XYZ_Z-FRAME`` used for for XYZ.


This is frame used to represent DMF of XYZ non-symmetri molecules with the :math:`z` (:math:`p`) axis along the vecror :math:`r_1` and other two axes defined using the following conditions:

.. math::

      \begin{split}
      \vec{p} &= \frac{\vec{r}_1}{r_1} \\
      \vec{y} &= \frac{\vec{r}_1 \times \vec{r}_2 }{|\vec{r}_1 \times \vec{r}_2|} \\
      \vec{q} &= \vec{y}\times \vec{p}
      \end{split}



The corrsponding components :math:`\mu^{(q)}` and :math:`\mu^{(p)}` are expanded using the same form as in Eq. :eq:`e-muQ-1` but with no constraints on the permutations of the atoms.





``DIPOLE_PQR_XYZ_Z-FRAME_SINRHO``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The same as ``DIPOLE_PQR_XYZ_Z-FRAME`` but with :math:`\sin(\alpha_{\rm e}-\alpha)` as an expansion variable in Eq. :eq:`e-muQ-1` instead of :math:`\cos\alpha_{\rm e} - \cos\alpha`:

.. math::
       :label: e-muQ-exp-3

       \begin{split}
         \mu_{jk\ldots}^{(q)}(\alpha)  =&  \sum_{i=0}^{N} Q_{ij\ldots}^{(i)} (\sin(\alpha_{\rm e}-\alpha))^i, \\
         \mu_{jk\ldots}^{(p)}(\alpha)  =&  \sum_{i=0}^{N} P_{ij\ldots}^{(i)} (\sin(\alpha_{\rm e}-\alpha))^i,
       \end{split}

When :math:`\alpha_{\rm e} = \pi` (linear molecules),    :math:`\sin(\pi-\alpha) = \sin\rho`, which explanes the suffix ``_sinrho`` in the name of thi DMF, wich is aimed at linear molecules. 




``DIPOLE_AMES1_XYZ``
^^^^^^^^^^^^^^^^^^^^

This form is a modification of ``DIPOLE_AMES1`` for non-symmetric molecules. 

As an example can be found of a system where this form was used, see :download:`16O-12C-32S__OYT8__TROVE.model <./input/16O-12C-32S__OYT8__TROVE.model>` as well in `OYT8 spectroscopic model <https://exomol.com/models/OCS/16O-12C-32S/OYT8/>`__, where it was used to compute an ExoMol line list for OCS [24OwYuTe]_.





``XY2_SCHROEDER_XYZ_ECKART``
----------------------------

This is an XYZ version of the ``XY2_C3_SCHROEDER`` type. 


DIPOLE_XY2_LORENZO

DIPOLE_H2O_LPT2011

DIPOLE_PQR_XYZ

DIPOLE_PQR_XYZ_Z-BOND

DIPOLE_PQR_XYZ_BISECTING

DIPOLE_BISECT_S1S2T_XYZ

XY2_QMOM_SYM

XY2_ALPHA_SYM

XY2_QMOM_BISECT_FRAME

TEST_XY2_QMOM_BISECT_FRAME

XY2_SR-BISECT-NONLIN

TEST_XY2_SR-BISECT-NONLIN


XY3
---


XY3_MB

XY3_MB4

XY3_SYMMB

XY3_NSS_MB

ABCD
-----

HOOH_MB

HPPH_MB

HCCH_MB

HCCH_DMS_7D

HCCH_DMS_7D_7ORDER

HCCH_DMS_7D_7ORDER_LINEAR

HCCH_DMS_7D_LOCAL

HCCH_ALPHA_ISO_7D_LINEAR


ZXY2
----

ZXY2_SYMADAP


ZXY3
----

ZXY3_SYM

C2H4
----

DIPOLE_C2H4_4M




XY2_SR-BISECT

XY2_SS_DIPOLE_YY

XY2_G-BISECT

XY2_G-ROT-ELEC

XY2_G-COR-ELEC

XY2_G-TENS-RANK3

XY2_G-TENS-NUC

COORDINATES

