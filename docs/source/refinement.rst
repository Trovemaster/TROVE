Refinement
**********

.. _refine:


In this chapter details will be given of the refinement procedure implemented in TROVE. Refinement in this context means to adjust the *ab intio* potential energy surface by comparing computed rotational-vibrational energy levels to experimental values. Parameters of the PES are varied, energy levels re-computed and compared to experiment. This process is continued until acceptable agreement between the calculated and experimental energy levels is obtained. Usually there is relatively few experimental energies and so *ab intio* electronic energies are used to constrain the refinement to prevent over fitting. Although experimental data is usually at fairly low energies, it is often the case that correcting the lower energy  region of the PES gives more accurate values at high energies also.

Refinement with TROVE: Theory
=============================

Details of the method of refinement implemented in TROVE have been published [11YuBaTe]_ and only a brief summary  will be given here. Assuming a reasonable PES has already been obtained, a correction is added in terms of internal coordinates :math:`\xi`

.. math::
     
    \Delta V = \sum_{ijk...} \Delta f_{ijk...} \left(\xi_1^i \xi_2^j \xi_3^k ...\right)^A
     
where :math:`\left(\xi_1^i \xi_2^j \xi_3^k ... \right)^A` corresponds to totally symmetric permutation of the internal coordinates so that all symmetry properties of the molecule are properly accounted for. :math:`\Delta f_{ijk}` are the expansion coefficients which are found by refinement. The Hamiltonian is now given as

.. math::
    
    H = T + V + \Delta V = H_0 + \sum_{ijk...} \Delta f_{ijk...} \left(\xi_1^i \xi_2^j \xi_3^k \right)^A
    
where :math:`H_0` is the initial Hamiltonian with *ab initio* PES.

If the eigenvalue problem for the initial Hamiltonian has been solved,

.. math::
    
    H_0 \psi^{J,\Gamma}_{0,i} = E^{J,\Gamma}_{0,i} \psi^{J,\Gamma}_{0,i}
    
where :math:`J` is the total angular momentum quantum number and :math:`\Gamma` is a symmetry label, then matrix elements of :math:`H`, using the :math:`H_0` solutions as a basis, are

.. math::
      
      \langle  \psi^{J,\Gamma}_{0,i} | H |\psi^{J,\Gamma}_{0,i'}   \rangle = E^{J,\Gamma}_{0,i} + \sum_{ijk...} \Delta f_{ijk...} \Xi_{i,i'}^{J, \Gamma}
      
where

.. math::
      
      \Xi_{i,i'}^{J, \Gamma} = \langle  \psi^{J,\Gamma}_{0,i} | \left(\xi_1^i \xi_2^j \xi_3^k ...\right)^A | \psi^{J,\Gamma}_{0,i'} \rangle.
       

The derivatives of the energies with respect to adjustable parameters, which are required for least squares fitting, are given by the Hellman-Feynman theorem

.. math::
      
      \frac{\partial E^{J,\Gamma}_{n} }{ \partial \Delta f_{ijk...} } = \langle \psi^{J,\Gamma}_{n} \left| \frac{\partial \Delta V}{\partial \Delta f_{ijk...} }       \right |\psi^{J,\Gamma}_{n} \rangle  = \langle  \psi^{J,\Gamma}_{n} \left| \left(\xi_1^i \xi_2^j \xi_3^k ...\right)^A \right| \psi^{J,\Gamma}_{n} \rangle .
       
where :math:`E^{J,\Gamma}_{n}` and :math:`\psi^{J,\Gamma}_{n}` are eigenvalues and eigenvectors of :math:`H` respectively. In the J=0 representation :math:`\psi^{J,\Gamma}_{n}` is given by

.. math::
     
     \psi^{J,\Gamma}_{n} = \sum_i C_i^{J, \Gamma} \psi_{0,i}^{J, \gamma}
      

Since TROVE typically uses linearised coordinates :math:`\xi^{\rm lin}_i` to represent PESs internally, i.e. in most cases different from the analytic representation of the *ab initio* PES, typically given in terms of some valence coordinates :math:`r_\lambda`, the TROVE eigenfunctions :math:`\psi^{J,\Gamma}_{n}` as well as basis functions :math:`\psi^{J,\Gamma}_{0,n}` are also solved and represented in terms of :math:`\xi^{\rm lin}_i`. Therefore, in order ro evaluate the integrals :math:`\langle  \psi^{J,\Gamma}_{0,i} | \left(\xi_1^i \xi_2^j \xi_3^k ...\right)^A | \psi^{J,\Gamma}_{0,i'} \rangle`, each combination :math:`\left(\xi_1^i \xi_2^j \xi_3^k ...\right)^A` has to be Taylor re-expanded in terms of the linearised coordinates :math:`\xi^{\rm lin}_i`. 

In the fittings the derivatives of the energy levels with respect to the correction parameters are evaluated at each fitting iteration in line with the  standard least squares fitting Newton procedures. This is all implemented in TROVE.


Refinement Implementation with TROVE
====================================

Setting up Refinement
---------------------

The specific inputs and checkpoint files required to carry out refinement of a PES using TROVE are discussed in this section.

Prior to refinement, TROVE requires checkpoint files and eigenfunctions for the basis set being used (see above). If a calculation of the rotational-vibrational levels using an unrefined PES has already been carried out, then all necessary files for refinement will have been generated. Refinement can be carried out in the :math:`J=0` basis.

 As explained above, refinement in TROVE is represented as a correction :math:`\Delta V(r)` to the *ab initio* PES :math:`V(r)` an represented by refinement parameters :math:`\Delta f_{ijk...}`. In the current TROVE implementation, the refinement part :math:`\Delta V(r)` is required to have exactly the same analytic representation as :math:`V(r)`, i.e. the refined PES is represented by the expansion parameters `f'_{ijk...}` given by 
 
 .. math::  
            f'_{ijk...} =  f^{\rm ai}_{ijk...} + \Delta f_{ijk...}
            
             
 While the potential is defined in the ``POTEN`` block, the refined PES  the ``external`` block on the TROVE input file. This is the same structure as used to define the ``dipole`` moment for intensity calculations and can assume a vector structure of dimension :math:`D`, for example in the case of DMS, the dimension is 3. For the refinement, each expansion term  :math:`\left(\xi_1^i \xi_2^j \xi_3^k ...\right)^A` is treated as an independent function and thus the ``external`` field is represented as a vector of dimension :math:`N`,  where  :math:`N` is the number of expansion parameters :math:`\Delta f_{ijk...}` and each vector elements holds a combination  :math:`\left(\xi_1^i \xi_2^j \xi_3^k ...\right)^A`. 
 
Thus the structure of the ``external`` parameter section is just a repeat of the ``potential`` block.

  ..Note:: Only linear parameters like :math:`\Delta f_{ijk...}` can be fitted in TROVE. Non-linear parameters such as equilibrium positions, structural parameters currently cannot be refined in TROVE. 
 
 A typical fitting ``external`` section has the following form
::

     external
     dimension 102
     Nparam  1
     compact
     TYPE  potential
     COEFF   list  (powers or list)
     dstep   0.005
     Order   4
     COORDS  morse morse linear
     parameters
     RE13            1.5144017558        fix
     alphae          92.00507388         fix
     a               0.127050746200E+01  fix
     b1              0.500000000000E+06  fix
     b2              0.500000000000E+05  fix
     g1              0.150000000000E+02  fix
     g2              0.100000000000E+02  fix
     V0              0.0000000000000000
     F_0_0_1     0.0              fit
     F_1_0_0     0.0              fit
     F_0_0_2    -0.173956405672E+05 fit
     F_1_0_1     0.241119856834E+04 fit
     F_1_1_0     0.0
     ....
     end
     
Here 

  - ``dimension`` is the number of all parameters :math:`N` of the PEF (number of lines between the cards ``parameters`` and ``end``). 
  - ``Nparam`` tells TROVE that each component of the :math:`N` dimensional external field (:math:`\left(\xi_1^i \xi_2^j \xi_3^k ...\right)^A`) is a single parameter object. Compare this with ``DIPOLE`` which can have 3 components each of which is an :math:`N_i`-dimensional analytic expansion (:math:`i=`1..3`) represented  with :math:`N_i` parameters. 
  - ``compact`` is the format without a fitting weight in the column before the parameter values (penultimate). 
  - ``type`` in the case of the fitting must be `potential` (in the current implementation), which tells TROVE to refer to the functional type of the PEF (``POT_TYPE``), see `Potentials  <https://spectrove.readthedocs.io/en/latest/potentials.html>`__. 
  - ``coeff``  is the card specifying whether  that the parameters are given as a ``list`` with predefined order as implemented in the code or via list with ``powers``-s (exponents). Here we use the ``list`` form.
    .. Note:: Although ``Coords`` in ``external`` does not have to coincide with that in ``POTEN``, it is advised to used the same type in both fields for consistency.
  - ``dstep`` defines the derivation step size used to evaluate high order derivatives with finite differences. 
  - ``Order`` defines the re-expansion order of each :math:`\left(\xi_1^i \xi_2^j \xi_3^k ...\right)^A` term.
  - ``Coords`` defines the types of the linearised coordinates used in the re-expansion. 
  - ``parameters`` indicates the beginning of the section with parameters. 
  - ``fix`` and ``fit`` are the keywords to distinguish the parameters to fix and parameters to fit. It is important that all structural parameters are marked with the ``fix`` card. This will insure that the derivatives and expansions of :math:`\left(\xi_1^i \xi_2^j \xi_3^k ...\right)^A` are evaluated correctly. ``fit`` needs to be set only to the parameters that will need t o be varied. 
  
  .. Note:: It is important to ``fix`` all structural parameters in the ``external`` section. For the example above, the potential function it is linked to has the ``pot_type`` ``POTEN_XY2_MORSE_COS`` and is defined as follows 
:: 

    POTEN 
    NPARAM  102
    POT_TYPE  POTEN_XY2_MORSE_COS
    compact
    COEFF  list  (powers or list)
    RE13              1.5144017558
    alphae            92.00507388
    a                 0.127050746200E+01
    b1                0.500000000000E+06
    b2                0.500000000000E+05
    g1                0.150000000000E+02
    g2                0.100000000000E+02
    V0                0.000000000000E+00
    F_0_0_1           0.000000000000E+00
    F_1_0_0           0.000000000000E+00
    F_0_0_2           0.173956405672E+05
    F_1_0_1          -0.241119856834E+04
    F_1_1_0           0.223873811001E+03
    ...
    ...
    end 
    

Here, the structural parameters are :math:`r_{\rm e}`, :math:`\alpha_{\rm e}`, Morse parameter :math:`a` as well as parameters :math:`b_1`, :math:`b_2`, :math:`g_1`, :math:`g_2`. The must be fixed to their values when doing the re-expansion of the external part. 


Here is another example, where the potential function type ``poten_C3_R_theta`` was used: 
::


     POTEN
     compact
     POT_TYPE  poten_C3_R_theta
     COEFF  powers  (powers or list)
     RE12          0      0      0     1.29397
     theta0        0      0      0     0.000000000000E+00
     f000          0      0      0        0.00000000
     f100          1      0      0        0.00000000
     f200          2      0      0        0.33240693
     f300          3      0      0       -0.35060064
     f400          4      0      0        0.22690209
     f500          5      0      0       -0.11822982
     .....
     .....
     end 
      

The ``external`` field is then given by 
::
     
     external
     dimension 60
     compact
     NPARAM  1
     compact
     type potential
     order 8
     coords morse morse   linear
     COEFF  powers  (powers or list)
     parameters
     RE12          1      0      0        1.29397   fix
     theta0        0      0      0        0.0000000 fix
     f000          0      0      0        0.00000000 fit
     f100          1      0      0        0.00000000 fit
     f200          2      0      0        0.00000000 fit
     f300          3      0      0        0.00000000 fit
     f400          4      0      0        0.00000000 fit
     f500          5      0      0        0.00000000 fit
     ....
     ....
     end
     

The fitted potential parameters in the ``external`` section  can be assumed  to be zero but never actually featured until step 4, so the actual values won't matter at steps 1,2,3. 
     
Calculation steps 
-----------------

At step 1 and 2, for the ``external`` field to be processed, the ``control`` block has to include the card ``external``: 

- Step 1:
::

    Control
    Step 1
    external
    end

- Step 2:
::

    Control
    Step 2
    end

See `Quick Start  <https://spectrove.readthedocs.io/en/latest/quickstart.html>`__. 

Step 3 does not involve any operations with the external field and therefore should be processed as usual, e.g. 
::

    Control
    Step 3
    J 1
    end

After step 3, in the case of the refinements, in the control block  we skip step 4 (``intensities``) and start step 5 (``fitpot``), at which matrix elements of the expansion terms :math:`\left(\xi_1^i \xi_2^j \xi_3^k ...\right)^A` are computed on the final ro-vibrational eigenfunctions obtained at step 3 for the *ab initio* model, for all values of :math:`J` and all symmetries considered. A typical Step 5 ``control`` block has the following structure:
::

    Control
      Step 5  (FITPOT)
      external  3 60
      J  0, 1, 2, 3
      symmetries 1, 2, 3, 4
    end

Here, 
 - 3-60 in the ``external`` is the range of the expansion terms (i.e. corresponding to the expansion parameters :math:`f^{\rm ai}_{ijk...}`) to be processed. 
 - `J` card  lists all values of :math:`J` to be processed. Note that it is not a range but a list i.e. the parameters can appear in any combination or order. Alias `Jrot`. 
 - `symmetries` card (alias `gamma`) is similar to the `gamma` card used in step 3. It gives a list of symmetries to be processed, again, in any combination or order. 

  .. Note:: Alias for ``step 5`` is ``step fitpot``. 

At ``step 6`` (alias ``step refinement``), the actual fits are taking place. At this step, the control block will have similar format as for step 5:
::

    Control
      Step 5  (FITPOT)
      external  3 60
      J  0, 1, 2, 3
      symmetries 1, 2, 3, 4
    end


At this step, as an additional to the ``control`` block, the user needs to include the ``fitting`` section, which will be extensively used at step 6 ``refinement``. Let's consider the following example of the ``fitting``  block:
::
 
 FITTING
 J-LIST         0
 symmetries     1 2 3 4
 itmax   0
 fit_factor     1e-4
 THRESH_COEFF   1e-20
 lock           100
 output         f01
 TARGET_RMS     1e-18
 fit_scale      0.25
 thresh_obs-calc  10
 robust  0.0
 geometries     poten.dat
 OBS_ENERGIES  26   (J  symmetry NN )
     0    1    1      0.00000  0  0   0   0   1.00
     0    1    3   1978.1533   0  0   0   2   1.00
     0    1    4   2005.469    0  1   0   0   1.00
     0    1    5   2952.7      0  0   0   3   1.00
     0    1    6   2998.6      0  1   0   1   1.00
     0    1    7   3907.4      0  1   0   2   1.00
 .....
 end
 






The required matrix elements of the refinement parameters are computed in stages. First matrix elements of the primitive basis functions used by TROVE are calculated. This is carried out by putting
::

     extmatelem save split n n

in the TROVE input file in the checkpoint block. n is the number of a specific expansion parameter. If all parameters are required this can be ignored but this may not be possible for all molecules or if there are lots of expansion parameters. This will generate extmatelemn.chk files (where n is number of expansion parameter). These files can then be converted to the `J=0`` representation using
::

     extmatelem convert split n n.


As discussed above, the refinement procedure requires matrix elements of the :math:`H_0` Hamiltonian and so eigenfunctions for each :math:`J` of interest must  be computed. To save matrix elements of the eigenfunctions, TROVE is run for each :math:`J` with
::

    fit_poten save split

in the checkpoint block. This generates ``fitpot-J-Gamma-n.chk`` files. Since a file is
generated for each :math:`J` and symmetry :math:`\Gamma` for each expansion parameter n, many files are generated in this step.


Running Refinement
------------------

When all of the required checkpoint files described above have been generated, TROVE can be used to refine the PES. The following block should be added to the TROVE input file
::

    FITTING
    J-LIST         0 1......
    symmetries     1 2 .....n
    itmax          10
    fit_factor     100
    lock           20.0
    output         fittest
    robust         0.0
    geometries     c2h4_pes.dat
    OBS_ENERGIES   52 (J  symmetry NN Energy 12 QNs weight  )
    0 1 2  1343.31   0   0   0   0   0   0   0   1   0   0   0   0   10  e  o
    .
    .
    .

``J-LIST`` is a list of total angular momentums to be included in the refinement, all checkpoint files for :math:`J` selected must have been already computed.

``symmetries`` is a list of symmetries to be included, again all checkpoint files for each :math:`\Gamma` must have already been computed.

``itmax`` is the number of iterations of refining carried out.

``fit_factor`` is the relative weighting for the experimental data compared to *ab initio* energies. The larger this is, the more importance will be given to the experimental energies.

``output`` is a string which specifies the pre-fix for output file names.

``robust`` specifies whether Watson Robust fitting is used, for 0.0 it is not, for 2.0 it is.

``geometries`` is the name of the file which contains *ab initio* energies. This file should give geometries in the same coordinates as specified by the potential energy surface for the molecule of interest in TROVE followed by the *ab initio* energy (from Molpro for example) and a weighting.

``OBS_ENERGIES`` is the number of observed (experimental) energies used. Below this a list of energies is given in the format
::

     J \Gamma NN E_i t_1 t_2 t_3 . . .    weight e o

where :math:`J` and :math:`\Gamma` are the angular momentum and symmetry number of the energies, NN is the block number, which is the number of the energy given by TROVE. The following numbers are the TROVE assignment of the energy level, followed by a weighting.

With the fitting block added to the input, TROVE can be used to refine a PES. In the external block ``NPARAM`` should be set to the number of parameters which are to be refined. In the list of parameters, the first column of integers specifies if a parameter is to be refined. `1` will include in refinement, `0` will exclude. The next column of real numbers are the starting values of the refinement parameters and should be set to 0.0 if initial refinement.

To carry out refinement all parts of the checkpoint block should be set to `read`` or `none``. TROVE will carry
out refinement until the number of iterations specified is
reached. The first iteration is essentially a checking step and does not change the value of the parameters.



Refinement Output
-----------------

The refinement procedure produces three output files. A regular .out file with a prefix the same as the .inp file and a
.pot file and .en file with prefixes as determined by the name given in the ``output`` keyword in the Fitting block.

The main output file for refinement is straightforward. The input is repeated as with other TROVE output files and then
some information is given about the eigenfunctions which were read in, etc. After this Trove prints the iteration number
and then a list comparing the observed to calculated energies. For example
::

    ---------------------------------------------------------------------------------
    | ## |  N |  J | sym|  Obs. | Calc.| Obs.-Calc. | Weight | K     vib. quanta
    ---------------------------------------------------------------------------------
    1  2  0  Ag  1343.5400  1346.2786  -2.7386 0.51E-03 (0) ( 0 0 0 0 0 1 0 0 0 0 0 0)*
    2  3  0  Ag  1625.4000  1632.5923  -7.1923 0.26E-03 (0) ( 1 0 0 0 0 0 0 0 0 0 0 0)
    3  4  0  Ag  1662.2000  1667.4972  -5.2972 0.26E-04 (0) ( 0 0 0 0 0 1 0 1 0 0 0 0)


The first number in a row is just a label to order the output. The second is the block number which was given to a particular state in the input file in the Fitting block. For the :math:`A_g` state in the example the first energy corresponding to a
fundamental mode has a block number of 2 since 1 would correspond to the ground state with relative energy of 0. After this the angular momentum of the state, :math:`J`, is given along with the symmetry. The observed energy as given in the input file
is then given followed by the current iterations calculation of the energy using the adjusted potential parameters and the difference between them. The weighting given to the state is then given. The rotational :math:`k` quantum number and vibrational
quantum numbers are then given. If an asterisk (*) is printed at the end of the row (as in the first row of this example) it means that TROVE has assigned the state differently to how it was labelled in the input in the Fitting block.

TROVE then prints a list of corrections to the potential parameters followed by the new values for the potential parameters and the corrections rounded according to their error.

A table is then printed which gives details on the fit for this iteration.
::

    -----------------------------------------------------------------------------------
    |  Iter | Points | Params | Deviat       | ssq_ener     | ssq_pot   | Convergence |
    -----------------------------------------------------------------------------------
    |  1    | 18107  |   21   |  0.34175E-01 |  0.61230E+01 | 0.173E+03 | 0.293E+12   |
    -----------------------------------------------------------------------------------

This gives the statistics of the fit including both the experimental energies and the *ab initio* energies used to constrain the fit.

The Obs-Calc table and fit statistics is then repeated for each iteration.

The .en file gives similar information to the Obs-Calc table in the output file but gives calculated energies for all states calculated by TROVE. The .pot file is a list of the *ab initio* geometries with the observed (that is,
the energy given for that geometry in the file listed under geometries in the Fitting bock) energies. The calculated  energy is also given, which is the energy given by the potential with the corrections from refinement, along with
zero-calc and the weight for the energy.







