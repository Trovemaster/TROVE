Line lists and Spectra
**********************
.. _linelists:

This chapter will give details on how to use TROVE and associated programs to make production quality line lists and spectra. That is, large line lists involving millions or even billions of transitions between states. The programs involved are called GAIN and Exocross_. Both have been designed to interface with TROVE outputs and using them does not require much additional syntax.


Intensity calculations are a part of step 4 (see Chapter Quickstart_ for the description of TROVE steps an Chapter Theory_ for theoretical details of the intensity calculations). 


Prerequisites of intensity calculations in TROVE
================================================

We assume that all three initial steps have been accomplished. This includes: solution of the ro-vibrational problem by computing the ro-vibrational wavefunctions :math:`\Psi_{\rm rv}^{(i)}` and energy term values :math:`\tilde{E}_i`, steps 1-3. The wavefunction coefficients :math:`C_{v'' k'' }^{(i)}` are stored in the check points ``j0eigen_vectors*.chk``, while the energies are in the ascii files ``j0eigen_descr*.chk`` for a set of :math:`J` from 0 to :math:`J_{\rm max}`. We also assume that the DMS components have fully been processed in steps 1 and 2. At step 1 this incudes (i) definition of DMS in the step 1 input file, (ii) re-expansion in terms of the TROVE internal vibrational coordinates (``external.chk``) and (iii) calculation of the vibrational matrix elements of the DMS components (``extmatelem1.chk``, ``extmatelem2.chk`` and ``extmatelem3.chk``). At step, the DM matrix elements are converted to the :math:`J=0` representation (``j0_extmatelem1.chk``, ``j0_extmatelem2.chk`` and ``j0_extmatelem3.chk``). Providing that these ingredients are available, the intensities can be generated as follows. 


INTENSITY input file: step 4
============================

In order to start the intensity calculations, the following changes in the input are required. In the ``CONTROL`` section, the keyword ``step`` has to be set to  4 and the range of the values of :math:`J` needs to be defined, e.g. 

::

    Control
    Step 4
    J 0 10
    end

Besides, an ``INTENSITY`` section has to be added anywhere in the input (see Quickstart_), e.g. 
::

      INTENSITY
       absorption
       exomol
       name *filename*
       THRESH_INTES  1e-20
       THRESH_LINE   1e-20
       THRESH_COEFF  1e-18
       TEMPERATURE   300.0
       Partition     1000.0
       GNS          8.0 8.0 8.0
       selection (rules) 1 1 2  (N of irreps)
       J,  0,10
       freq-window  -0.1,   4000.0
       energy low   -0.1,  2000.0, upper   -0.1, 6000.0
     END

Here 
``absorption`` specifies that absorption intensities between states are to be calculated.


``THRESH_INTES/LINE/COEFF`` are used to control the level of print out for intensities. Very large outputs can be produced if these are set very low (as needed for 'production' quality line lists) but for quicker checks higher values should be used.


``TEMPERATURE`` is used to specify the temperature of interest. This will affect the population of states (Boltzmann population).


``Partition`` is the value of the partition function. This can be calculated from all of the ro-vibrational energy levels used. Note that at high temperatures enough energy levels must be included for accurate results. If this is not the case (for example, for a test calculation) then a literature value could be used.


``GNS`` is the spin statistical weights for each symmetry. These can be looked up for many molecules or worked out from the procedure in Bunker and Jensen,  chapter 8 [98BuJe]_. ``selection`` is used to specify which symmetries can make up the initial and final states of a transition. The product of the upper and lower eigenfunctions must contain a component of the dipole itself [98BuJe]_. Thus for the PF\ :sub:`3` example, :math:`A_1` and :math:`A_2` are grouped together while E can only go to E. Integers are used to form groups, in this case 1 1 are for :math:`A_1` and :math:`A_2` and 2 is for E.


``J,  i,j`` specifies the rotational states to be included. In the example above 0 to 10 were used. It is often better to split a calculation into 0,1-1,2-2,3, etc to fit into time allocations on computers. The vibrational states to be included can also be specified by the ``v i, lower x, y, upper x', y'``
where i is the number of a vibrational mode and x, x' and y, y' give the limits for the lower and upper states included. If this is not included then all vibrational states are considered.


``freq-window`` This specifies the frequency window (in wavenumbers) in the spectra to be used. In the example here -0.1 is used as the minimum to guarantee values from 0 are used while 4000 is the maximum considered. ``energy low`` specifies the energies of the lower and upper states to be included. In the example the highest energy lower state to include it 2000 so since the maximum frequency of light considered is 4000, the upper state needs a maximum of 6000 (energy proportional to frequency, :math:`E = h \nu`).

``exomol``is the keyword to produce a line list in the ExoMol format, which will be written into the .trans and .states files. In this case, the absorption intensities computed for the temperature and partition functions specifies are not printed out only used together with the corresponding threshold ``THRESH_INTES`` to decide to keep the line in the line list or not.


``name``  is the *filename* of the .trans and .states files.

To calculate absorption intensities the eigenfunctions and eigenvalue files of the states to be included must be included in the directory where TROVE is run. More on this will be described below.



Choosing absorption parameters
==============================

The intensity block in TROVE requires a choice for the minimum intensities to be printed out and for the range of states and frequencies to be included. The value of intensity thresholds should be set very small for production quality line lists, for example ``THRESH_INTES`` and ``THRESH_LINE`` can be set at :math:`1\times 10^{-50}` to ensure all transitions are included.

Values for ``freq-window`` and ``energy low`` and ``upper`` depend on the molecule and temperature of interest. The lower energy range required will depend on the desired temperature range. For room temperature line lists, only relatively low energy states will be significantly populated. For hot line lists, this range will be increased. The partition function for the molecule can be used to judge which states are required for coverage at a certain temperature (see below for how to calculate using ExoCross_). The frequency window (and thus upper states to include) depends on the frequency of light which of interest.
Of course, the range which is included will also be limited by practical considerations such as computational time, memory, basis set convergence, etc.

Intensities with GAIN
=====================

As discussed in Chapter Quickstart_, TROVE is capable of calculating transition intensities once the relevant eigenfunctions and dipole matrix elements have been calculated. This procedure was used in early line list papers using TROVE.

A more efficient way of calculating intensities is to make use of the GAIN program. GAIN (GPU Accelerated INtensities) is a program which was written by Ahmed Al-Rafaie (see GAIN_). It uses graphical processing units (GPUs) to calculate intensities far quicker than can be achieved using conventional TROVE.

GAIN uses the same input file as TROVE but only the ``intensity`` block is actually used to control the calculation. GAIN requires the eigenvectors, eigen description and eigen quanta files for the states of interest. It also requires the eigen descrption and eigen quanta of the :math:`J = 0` state and extfield file for the dipole matrix elements (note that currently GAIN cannot accept split dipole files, these must be stitched together).

GAIN can be found at gitGAIN_. 

Using GAIN
----------

The number of states for a polynomial molecule quickly increases with :math:`J` and energy. This leads to millions of transitions and so even with GAIN, intensity calculations scale quite drastically. There are a few ways in which calculations can be sped up however so that they can be run within wall clock limits.

The first is to increase the number of nodes used. GAIN is an mpi parallel program and can make use of multiple nodes, which themselves have multiple cores. A rule of thumb for how many cores to use is: size of eigenvectors / memory available per core.

Another speed increase is to split the intensities which are being calculated by :math:`J` and symmetry. Rotational selection rules limit transitions to :math:`J'' = J'` and :math:`J'' = J' \pm 1`. Currently GAIN does not have a rule for only computing upper or lower Q branch transitions and so these duplicates should be removed for a complete line list. Using the selection rule, intensities can be calculated by setting :math:`J` in the ``intensity`` block to 0,1 then 1,2 then 2,3, etc. Symmetry also limits transitions but these are molecule dependent. For example, for PF\ :sub:`3` transitions can only take place for :math:`A_1` :math:`\leftrightarrow` :math:`A_2` and E :math:`\leftrightarrow` E. To make use of this symmetry the nuclear statistical weights (:math:`g_{ns}`) for the symmetries which are allowed should be set to their usual values but others set to 0. For example for :math:`A_1 \leftrightarrow A_2` in PF\ :sub:`3` the :math:`g_{\rm ns}` would be set to ``8.0 8.0 0.0``. For both :math:`J` and symmetry selection rules, a separate input file and run of GAIN should be carried out for each selection rule.

GAIN produces two types of output files. The .out files begins with a repeat of the input file. Information is then given on which .chk files were opened and which GPUs are being used and their memory. Information is then given on how GAIN splits up the calculation and how many transitions are to be computed. GAIN then cycles through the energy states starting with the lowest energy and computes all transitions to higher energies. For each complete lower energy calculation the current lines per second computed (L/s) is reported along with the predicted total time required. The other output file produced is a ``__n__.out`` file. Here ``n`` is an integer starting at 0 going up to number of nodes :math:`-1`. This file(s) contain the GAIN results and lists the frequency and the Einstein A coefficient \cite{98BuJexx} for a transition. Labels are also given for which states the transition is between.

Einstein A coefficients are calculated as opposed to intensities as these are temperature (and pressure) independent. To simulate and plot a spectrum Exocross_ is used which is discussed in the next section.


Currently the format for the intensities from GAIN is not compatible with Exocross_. Programs can be used however to convert the GAIN output to the slightly more compact Exomol format. Code for doing this can be obtained from Sergey Yurchenko. In the future it may be that GAIN is modified to directly output the correct format for Exocross_.


Exocross
========

As discussed above, GAIN produces a list of temperature and pressure independent Einstein A coefficients. To simulate a spectra, these must be converted into intensities. This can be achieved using ExoCross_ [ExoCross]_, providing the data is correctly formatted. TROVE can directly produced intensities but ExoCross_ has features which make it a better choice for production quality simulations.

To run ExoCross_, two types of file are required. A TRANS file which contains information about the intensity of transitions and a STATES file which contains the energy levels of the molecule. These files should be obtained using a program to change GAIN output or from TROVE directly. This is likely to change depending on situation and will not be discussed here.

A simple but important calculation which can be performed using ExoCross_ is finding the partition function at a given temperature. This is determined from the States file only. An example input is
::

     mem 63 gb

     partfunc
      ntemps 10
      tempmax 800 (K)
     end

     NPROCS 12

     verbose 4

     States C2H4_v01.states

The keyword ``partfunc`` is used to select a partition function calculation. ``ntemps`` is the number of partitions of the temperature, ``tempmax`` which will be calculated. For this example the partition function will be calculated at 80, 160, ... and 800 K. ``verbose`` is the level of print out. ``States`` is the name of the states file.

The output for this calculation is simple. A repeat of the input is first given and then the partition function calculation for each temperature is given in columns. The running total of the partition function with :math:`J` is given in rows.



ExoCross_ can also be used to make a `stick spectrum``. This is an idealised spectrum where each absorption is only represented by a line at a given wavenumber and intensity and broadening effects (doppler, collision, etc) are ignored. An input example is
::

     mem 63.0 gb

     Temperature  296
     Range 0 9000.0

     Npoints 90001

     absorption
     stick

     mass 28
     threshold 1e-25

     pf 11000.0


     output C2H4_thr_1e-25_T296

     ncache 1000000

     NPROCS 16

     verbose 4

     States C2H4_v01.states

     Transitions
      c2h4_initial_vib_2016_intense_j0_j1__0__.out_0.-9000..trans
      c2h4_initial_vib_2016_intense_j1_j2__0__.out_0.-9000..trans
      ...
      ...
     end

``Temperature`` is the temperature of interest in Kelvin.

``Range`` specifies the wavelength range to be used,
in this case 0 to 9000 cm\ :sup:`-1`.

``Npoints`` controls the density of the grid produced. In this example there will be
10 points per cm\ :sup:`-1`.

 ``absoprtion`` specifies that a spectra is to be computed and ``stick`` indicates
that a stick spectrum is required.

``mass`` is the molecule's mass in atomic mass units.

``threshold`` is the minimum intensity of transition to be included. This is important for keeping the output file
manageable so it can be used for making plots.

``pf`` is an optional keyword which is used to give the value
of the partition function rather than calculate it from the States file (the default case). This is useful if, for example,
not all :math:`J` have been calculated but you want to check the spectrum looks reasonable.


``output`` specifies what to call the output file.

``ncache`` is how much memory will be cached on the cpu during calculations. ``nprocs`` is the number of threads
to use.

``States`` is the States file to use and ``Transitions`` is a list of Trans files to use.


ExoCross_ has other options for simulating spectra. Examples include accounting for line broadening by using Gaussian or Voigt profiles for each line. The effects of particular background gas collisions can also be taken into account. These features are fully discussed in a recent publication and manual for the ExoCross_ program and the reader is directed there for full details [ExoCross]_.





.. _ExoCross: https://github.com/Trovemaster/ExoCross


.. _Quickstart: https://spectrove.readthedocs.io/en/latest/quickstart.html

.. _Theory: https://spectrove.readthedocs.io/en/latest/theory.html#intensity-calculations-in-trove

.. _gitGAIN: https://github.com/ExoMol/GAIN-MPI



