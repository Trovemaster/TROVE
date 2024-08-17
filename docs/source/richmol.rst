Interfacing with other molecular codes
**************************************

RichMol
=======

`Richmol <https://github.com/CFEL-CMI/richmol>`__ is a package aimed at providing simple and efficient platform for simulations of ro-vibrational energies, spectra, and dynamics in the presence of external (laser) and induced internal (hyperfine) electromagnetic fields for general molecule. Richmol can be interfaced with other variational codes, like TROVE and Duo.

The role of TROVE is to provide the interacting molecular properties (dipole moments, polarizabilities etc) in the so-called spectral representations, i.e.inthe representation of field-free ro-vibrational eigenfunctions. Consider, e.g., a Hamiltonian operator describing the molecule interacting with a general time-dependent electro-magnetic field :math:`\vec{F}(t)`:
.. math:: 
    
    \hat{H} = \hat{H}_0 +  (\vec{\mu} \cdot \vec{F}) + \frac{1}}{} (\vec{F} \underline(\alpha) \vec{F}) + \cdots 
    
where :math:`\vec{\mu}` is a dipole moment vector and :math:`\underline{\alpha}` is polarizability rank 2 tensor. In the spectral representation,  :math:`\hat{H}` is give by the matrix elements 
.. math:: 
    :label: e-Hpsi-matrix
    \begin{split}
    \langle \Phi_{n,m}^{J,\Gamma} |\hat{H}| \Phi_{n',m'}^{J',\Gamma'}  \rangle  & = E_{n,m}^{J,\Gamma} \delta_{n,n'}\delta_{m,m'} \delta_{J,J'}\delta_{\Gamma,\Gamma'} + \\
            & \sum_{\beta=X,Y,Z} \langle \Phi_{n,m}^{J,\Gamma} |\vec{\mu}_{\beta}| \Phi_{n',m'}^{J',\Gamma'} F_{\beta} + \\
            & \sum_{\beta,\gamma=X,Y,Z}  F_{\gamma} \langle \Phi_{n,m}^{J,\Gamma} |\vec{\alpha}_{\beta}| \Phi_{n',m'}^{J',\Gamma'} F_{\beta} + \cdots
    \end{split}
    
where :math:`\Phi_{n,m}^{J,\Gamma}` is an eigenfunction and :math:`E_{n,m}^{J,\Gamma}` is an eigenvalue of the field free Hamiltonian:
.. math:: 
      :label: e-H0psi
      \hat{H}_0 \Phi_{n,m}^{J,\Gamma} = E_{n,m}^{J,\Gamma} \Phi_{n,m}^{J,\Gamma}
       
:math:`J` is the rotational angular momentum, :math:`m` is its projection on the space-fixed axis :math:`Z`, :math:`\Gamma` is an irrep of the (field-free) molecule and :math:`n` is a general eigen-state counting number. Thus, this is  TROVE's role to solve Eq. :ref:`e-H0psi` and to compute the corresponding matrix elements of the properties in Eq. :math:`e-Hpsi-matrix`. It is then the role of RichMol to find the time dependent solution for the molecule interacting with :math:`\vec{F}(t)` and thus to solve for the molecular dynamics:

.. math:: 
    
    \hbar \frac{\partial \Psi(t)}{\partial t} = \hat{H}\Psi(t). 
    

Similarity, RichMol can be also used to solve a  general problem in a spectral representation of an unperturbed Hamiltonian. 


How to interface TROVE with RichMol
===================================

First of all, the user must provide a Fortan subroutine for the property in question  for calculating their values at any arbitrary molecular geometry :math`\vec{r}`. Consider e.g. a polarizability tensor :math:`\underline{\alpha}` as a set of surface components :math:`\alpha_{\beta,\gamma}`, where :math:`\beta,\gamma=x,y,z` are projections on the molecular fixed axes. In TROVE, the property is processed through the block ``External`` in the same way as for the dipoles or for the correction of PEF used in refinements.

Step 1
------

::

    Control
     Step 1
     external
    end


Follow the same general procedure described in `TROVE Quick start <https://spectrove.readthedocs.io/en/latest/quickstart.html>`__ and use the ``external``
(:code:`EXTERNAL`) block to define the polarizability components, e.g. 
::
      
      external
        dimension 6
        NPARAM  5 2 2 2 2 2
        DMS_TYPE ALPHA_C2H6_ZERO
        COEFF   list
        COORDS  linear
        Order   0
        compact
        parameters
        re1     1.52563613
        re2      1.09068849
        beta   111.19731929
        alphaxx0      -26.595
        alphaxx3      0.1004
        alphaxy0      0.0
        alphaxy3      0.0
        alphaxz0      0.0
        alphaxz3      0.0
        alphayy0     -26.595
        alphayy3      0.1004
        alphayz0      0.0
        alphayz3      0.0
        alphazz0      -30.336
        alphazz3      0.0824
      end

which represents a simplistic form of the polarizability tensor of C\ :sub:`2`\ H:sub:`6` using the TROVE function ``ALPHA_C2H6_ZERO``.Here, there are six independent components each of which is represented by a single  value (cards ``alpha***``) at the molecular equilibrium (cards ``re1``, ``re2`` and ``beta``).
 

Step 2 
------

Business as usual: 
::

    Control
     Step 2
     external
    end

Step 3
------

Business as usual, e.g.:
::

    Control
     Step 3
     J 0 
    end


Step 4
------

This is the main step of computing the matrix elements of :math:`\alpha_{\beta,\gamma}`, for which the `Intensity` card is used. We first define the calculation step 4 in the control block (anywhere in the input file):

::

    Control
     Step 4
     J 0,1
    end

and then define the ``intensity`` block using the RichMol-related cards ``field_me`` in conjunction with ``oper_alpha``, e.g.
::
   
   INTENSITY
     field_me
     oper alpha
     THRESH_INTES  1e-10
     THRESH_LINE   1e-10
     THRESH_COEFF  1e-20
     GNS          6.0 10.0 6.0 10.0 4.0 4.0 2.0 6.0 12.0 0 0 0 0 0 0 0 0 0
     selection (rules) 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
     J, 6,8
     freq-window  0, 10000
     energy low   -0.001, 10000.00, upper   -0.001, 10000.0
   END


The keyword ``field_me`` is to switch on the "Field's Matrix Elements". The keyword ``oper`` is to specify which type of the property to process; in this case it is the polarizability (``alpha``). 

Currently, the following properties are available in TROVE (see module :code:`extfield.f90`): 

- ``ALPHA``
- ``MU``

as well as 

- ``QUAD``
- ``SPINROT``
- ``SPINSPIN``
- ``GTENS``
- ``WIGNER``
- ``COSTHETA``
- ``J``
- ``COS2THETA``
- ``RICHMOL_LEVELS_FILE``
- ``MF_TENSOR``



