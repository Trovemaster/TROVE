Potential energy functions
**************************

TROVE provides a larger number of potential energy functions (PEFs) for different molecules already implemented. Most of these PEFs are in modules ``pot_*`` contained in  file ``pot_*.f90``.

 - ``pot_xy2.f90``
 - ``pot_xy3.f90``
 - ``pot_zxy2.f90``
 - ``pot_abcd.f90``
 - ``pot_xy4.f90``
 - ``pot_zxy3.f90``
 - ........

These are a part of the standard TROVE compilation set. Alternatively, a user-defined PEF can be included into the TROVE compilation as a generic 'user-defined' module ``pot_user``, see details below.





Potential Block
===============

The ``Potential`` (``Poten``) block used to specify a PEF, has the following generic structure (using XY\ :sub:`2` as an example)

::

      POTENTIAL
      NPARAM  99
      POT_TYPE  poten_xy2_tyuterev
      COEFF  list
      b1        0    0.80000000000000E+06
      b2        0    0.80000000000000E+05
      g1        0    0.13000000000000E+02
      g2        0    0.55000000000000E+01
      f000      0    0.00000000000000E+00
      f001      1    0.25298724728304E+01
      f100      1    0.76001446034650E+01
      ...
      end

For an example, see `h2s_step1.inp <./input/h2s_step1.inp>`_  where this PES is used.  

Here ``NPARAM`` is used to specify the number of parameters used to define the PES. ``POT_TYPE`` is the name of the potential energy surface being used which is defined in the ``pot_*.f90 file``. The keywords ``COEFF`` indicates if the potential contains a list of parameter values (``LIST``) or  values with the corresponding expansion powers (``POWERS``), e.g. (for H\ :sub:`2`\ CO):
::

     POTEN
     NPARAM  114
     POT_TYPE  poten_zxy2_mep_r_alpha_rho_powers
     COEFF  powers
     f_0_0      0 0 0 0 0 0   1  0.00000000000000E+00
     f_0_1      1 0 0 0 0 0   1  0.00000000000000E+00
     f_0_1      0 1 0 0 0 0   1  0.00000000000000E+00
     f_0_2      0 0 0 1 0 0   1  0.00000000000000E+00
     f_1_1      0 0 0 0 0 1   1  0.13239727881219E+05
     f_2_1      0 0 0 0 0 2   1  0.46279621687684E+04
     f_3_1      0 0 0 0 0 3   1  0.14394787420943E+04
     f_4_1      0 0 0 0 0 4   1  0.10067584554678E+04
     f_5_1      0 0 0 2 0 0   1  0.31402651686299E+05
     .....
     end

which is from the variational calculations of H\ :sub:`2`\ CO, see the TROVE input `1H2-12C-16O__AYTY__TROVE.inp <input/1H2-12C-16O__AYTY__TROVE.inp>`_. 


The potential parameters are listed after the keyword ``COEFF`` and terminated with the keyword ``END`` with exactly ``NPARAM``. For the ``COEFF  list`` option, the meaning of the columns is as follows:

  +---------+-----+-----------------------+
  | Label   |Index| Value                 |
  +=========+=====+=======================+
  |   b1    |   0 |  0.80000000000000E+06 |
  +---------+-----+-----------------------+
  |   b2    |   0 |  0.80000000000000E+05 |
  +---------+-----+-----------------------+
  |   g1    |   0 |  0.13000000000000E+02 |
  +---------+-----+-----------------------+
  |   g2    |   0 |  0.55000000000000E+01 |
  +---------+-----+-----------------------+
  |   f000  |   0 |  0.00000000000000E+00 |
  +---------+-----+-----------------------+
  |   f001  |   1 |  0.25298724728304E+01 |
  +---------+-----+-----------------------+
  |   f100  |   1 |  0.76001446034650E+01 |
  +---------+-----+-----------------------+


Here 'Labels' are the parameter names,  used only for printing purposes and not in any calculations. The 'Index' field can be used to as a switch to indicate if the corresponding parameter was fitted or can be fitted. Otherwise it has no impact on any evaluations of the PEF values. 'Values' are the actual potential parameters, listed in the order implemented in the corresponding PEF ``POT_TYPE``, for example

.. math:: 
   
   \begin{split}
   V(r_1,r_2,\alpha) &= f_{000} + f_{001} y_3 + f_{100} [ y_1 + y_2 ] + f_{100} [ y_1 + y_2 ] + f_{002} y_3^2 + \ldots +  \\
                     & + b_1 e^{-g_1 r_{\rm HH}} + b_2 e^{-g_2 r_{\rm HH}^2} \\
   \end{split}
   
 

where

.. math:: 

   \begin{split}
      y_1 & = 1-e^{-a (r_1 - r_{\rm e})}, \\
      y_2 & = 1-e^{-a (r_2 - r_{\rm e})}, \\
      y_3 &= \cos\alpha-\cos\alpha_{\rm e}, \\
      r_{HH}=\sqrt{r_1^2+r_2^2-2 r_1 r_2 \cos\alpha}.\\
   \end{split}
   
   
   
For the ``COEFF  Powers`` option, the meaning of the columns is as follows:

   +---------+-----+--+--+--+--+--+-----+-----------------------+
   | Label   |   n1|n2|n3|n4|n5|n6|Index| Value                 |
   +---------+-----+--+--+--+--+--+-----+-----------------------+
   |  f000000|    0| 0| 0| 0| 0| 0|    1|  0.00000000000000E+00 |
   +---------+-----+--+--+--+--+--+-----+-----------------------+
   |  f100000|    1| 0| 0| 0| 0| 0|    1|  0.00000000000000E+00 |
   +---------+-----+--+--+--+--+--+-----+-----------------------+
   |  f010000|    0| 1| 0| 0| 0| 0|    1|  0.00000000000000E+00 |
   +---------+-----+--+--+--+--+--+-----+-----------------------+
   |  f000100|    0| 0| 0| 1| 0| 0|    1|  0.00000000000000E+00 |
   +---------+-----+--+--+--+--+--+-----+-----------------------+
   |  f000001|    0| 0| 0| 0| 0| 1|    1|  0.13239727881219E+05 |
   +---------+-----+--+--+--+--+--+-----+-----------------------+
   |  f000002|    0| 0| 0| 0| 0| 2|    1|  0.46279621687684E+04 |
   +---------+-----+--+--+--+--+--+-----+-----------------------+
   |  f000003|    0| 0| 0| 0| 0| 3|    1|  0.14394787420943E+04 |
   +---------+-----+--+--+--+--+--+-----+-----------------------+

where 

 - 'Labels' are the parameter name,  for printing purposes only; 
 - 'n1', 'n2', 'n3', ... are the 'powers' of an expansion term, e.g. 
   :math:`V(r_1,r_2,r_3,r_4,r_5, r_6) = \sum_{n_1,n_2,n_3,n_4,n_5,n_1} f_{n_1,n_2,n_3,n_4,n_5,n_1} \xi_1^{n_1} \xi_2^{n_2} \xi_3^{n_3} \xi_4^{n_4} \xi_5^{n_5} \xi_6^{n_6}`
  - 'Index' is a switch to indicate if the corresponding parameter was fitted or can be fitted, with no impact on any evaluations of the PEF values. 
  - 'Values' are the actual potential parameters. Their order is not important for this implementation as long as the corresponding powers are defined. 



In case the definition of PEF requires also structural parameters, such as equilibrium bond lengths :math:`r_{\rm e}`\ , equilibrium inter-bond angles :math:`\alpha_{\rm e}`, Morse exponents :math:`a` etc., in the ``COEFF  Powers`` form these parameters should be listed exactly in the order expected by the  implemented of the PEF (similar to the ``COEFF LIST`` form), but with dummy "powers" columns so that their 'values' appear in the right column. For example: 
::

    POTEN
    NPARAM  58
    POT_TYPE  poten_C3_R_theta
    COEFF  powers  (powers or list)
    RE12          0      0      0      0     1.29397
    theta0        0      0      0      0     0.000000000000E+00
    f200          2      0      0      0        0.33240693
    f300          3      0      0      0       -0.35060064
    f400          4      0      0      0        0.22690209
    f500          5      0      0      0       -0.11822982
    .....
    

Here, ``RE12`` and ``theta0`` are two the equilibrium values and the three columns with ``0 0 0`` are given in order to parse their values using column 6. 

Implemented PEFs
================


XY\ :sub:`2` type
-----------------

There are several PEFs available for this molecule type.


``POTEN_XY2_MORBID``
^^^^^^^^^^^^^^^^^^^^

This form is given by 

.. math:: 

   \begin{split}
   V(r_1,r_2,\alpha) &= f_{000} + f_{001} y_3 + f_{002} y_3^2 +  + f_{003} y_4 + \ldots  \\
                   & + (f_{100} + f_{101} y_3 + f_{102} y_3^2 +  + f_{103} y_4 + \ldots) [y_1 + y_2] \\ 
                   & + (f_{200} + f_{201} y_3 + f_{202} y_3^2 +  + f_{203} y_4 + \ldots) [y_1^2 + y_2^2] \\
                   & + (f_{110} + f_{111} y_3 + f_{112} y_3^2 +  + f_{113} y_4 + \ldots) y_1y_2 \\
                   & + \ldots \\

   \end{split}

where

.. math::

   \begin{split}
      y_1 & = 1-e^{-a (r_1 - r_{\rm e})}, \\
      y_2 & = 1-e^{-a (r_2 - r_{\rm e})}, \\
      y_3 &= \cos\alpha-\cos\alpha_{\rm e}, \\
   \end{split}


**Examples** and **References**



