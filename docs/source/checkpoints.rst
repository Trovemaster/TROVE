Checkpoints
***********
.. _checkpoints:

TROVE provides an extended set of the so-called 'checkpoints', with the main functionality to provide the capability of restarting the calculations at any critical point. The checkpoints are TROVE files with extension .chk that contain the necessary information for the restart. Some of the files are in the ASCII format (i.e. 'formated'), some are written in as machine readable, i.e 'unformatted', some - as direct access files. For instance, the transition between steps 1,2,3, ... is driven by the checkpoint functionality. Moreover, since step 1 contains a number of critical sub-steps, TROVE checkpoints those sub-steps as well. 

Apart from the restarting functionality, the .chk files are also used to store the eigenvectors (i.e. eigencoefficients). Although these objects are formally a the end result of the calculations, in same cases they provide support for some intermediate steps. For example, step 2 (transformation of matrix elements to the :math:`J=0` representation) used the eigenfunctions of the :math:`J=0` solution obtained at step 1 to create a basis set for step 3 (production of eigenfunctions for :math:`J>0`).

The checkpoints control section is ``CHECK_POINT``, e.g.:
::

     CHECK_POINT
      ascii
      kinetic     read
      potential   read
      external    none
      basis_set   save
      contract    save
      matelem     save    split
      extmatelem  none    split
      eigenfunc   save
     END
     

with the main control keywords ``save``, ``read`` and ``none`` and  the additional keywords ``append``, ``stitch``, ``split`` etc. The functional checkpoint keywords include 

- ``kinetic`` 
- ``potential``
- ``external``
- ``basis_set``
- ``contract``
- ``matelem``
- ``extmatelem``
- ``eigenfunc``
- ``fit_poten``

The order of these lines is unimportant. The ``CHECK_POINT`` section can appear anywhere in the input file. 



List of checkpoints
===================

- kinetic.chk (``kinetic``)  contains the expansion coefficients of the KEO as controlled by the ``Kinorder`` keyword. 

  The KEO is the first to be constructed and saved if ``kinetic`` is set to ``save`` and can be then used at any stage by switching to ``read``. One opt out from saving the KEO (or any other checkpoint) by setting it to ``none``.

- potential.chk (``potential``) is the second object created in TROVE as part of step 1. 

  It contains expansion coefficient of PEF in terms of the TROVE vibrational coordinates. It can be saved, if generated from scratched and then read after a TROVE restart.

- external.chk (``external``) is used to store expansion coefficients of any other objects apart from KEO and PEF, 

  dipole moment, polarizability, spin-rotation, a correction to the PEF used for the refinement of the PEF etc, which are called ``External`` or ``Dipole`` The same functionality as above is applied to the External field.

- hamiltoian.chk: 

  This file  contains the definition of the linearised coordinates (A and B matrices, see the TROVE paper [TROVE]_), i.e. for ``COORDS linear`` in the case of the linearised coordinate. It is saved and read as part of the ``kinetic`` switch.  hamiltoian.chk is not produced for the curvilinear coordinate type ``COORDS local``.

  kinetic.chk,  potential.chk and external.chk have the ASCII format and can be used to plot the corresponding fields or their individual components.  

- numerov_bset.chk and prim_bset.chk (``basis_set``). 

  These two unformatted files containing information on the 1D primitive basis sets and are produced after the objects described above. The same rules are applied for saving, reading or doing nothing (``none``). 


- contr_descr.chk, contr_quanta.chk, contr_vectors.chk (``contract``) are three files to store the information on the symmetry adapted basis functions for individual basis sub-sets. 

  - contr_vectors.chk is an unformatted file containing eigne-coefficients. 
  - contr_quanta.chk is a formatted (ASCII) file with the information on the basis set bookkeeping - mapping of the multimode quantum numbers into a 1D array. 
  - contr_descr.chk is a formatted file containing descriptions of the basis sets: 
  
    individual energies of the basis functions from different sub-classes together with their classifications, symmetries, TROVE quantum number, IDs, largest expansion coefficients used in their assignment as well as a placeholder for the spectroscopic quantum numbers. This file can be edited in order to include these spectroscopic quantum numbers.
  
- contr_matelem.chk (``matelem``) contains vibrational matrix of the different pars of the Hamiltonian operator (KEO and PEF)

  on the basis set produced at the ``contract`` step, i.e. contracted symmetry adapted basis functions. 

- matelem1.chk, matelem2.chk ... matelem12.chk represent a ``split`` version of contr_matelem.chk, which contain the rotational and Coriolis parts of the KEO, 
  
  while contr_matelem.chk contains the matrix elements of the pure vibrational part of the Hamiltonian operator, i.e. vibrational KEO + PEF. The ``split`` option allows one to compute all these 13 parts of ``matelems`` separately, which is especially useful for larger molecules. Historically, the ``matelem`` step the main bottleneck in the TROVE pipeline in the case of large molecules and the ``split`` feature allows one to at least partly mitigate this issue. An example of the ``split`` option is 
  
  - To process all 13 parts of the Hamiltonian: 
  :: 
   
   contract  save split 
   
  - To process individual parts of the Hamiltonian (from 0 to 6):
  ::
   
   contract  save split 0 6 
   
   here 0 stands for the pure vibrational part of the Hamiltonian operator.  
 
  - To process a single part of the Hamiltonian (from 11 to 11):
  ::
     
     contract  save split 11 11 
     
- contr_extfield.chk ``extmatelem`` contains all (vibrational) matrix elements of the external field. 

  ``extmatelem`` step is not a compulsory step in the TROVE pipeline. It is invoked when keyword it is set to ``save``. 
  
- extmatelem1.chk, extmatelem2.chk, extmatelem3.chk .... are the ``split`` analogy of contr_extfield.chk, 

  where different components are written into separate extmatelem chk-files.
  

Example of the split option include 
:
     
     extmatelem  save split 1  1
     
:

     extmatelem  save split 
     
     
   
- eigen_descr*chk, eigen_vector*chk and eigen_quant*chk (``eigenfunc``) contain the  eigencoefficients and their descriptions. 

  - eigen_descr\ :math:`J`\ _\ :math:`G`\ .chk contain the eigenvalues (energy term values in cm\ :sup:`-1`\ ):
   
   state IDs, symmetries, TROVE quantum numbers, largest coefficients as well as a placeholder for the spectroscopic quantum numbers. These files formatted (ASCII) and can be used for the analysis or postprocessing (e.g. construction of line lists). Here :math:`J` is the rotational angular momentum and :math:`G` is the symmetry (irrep), i.e. there is a description for each J/symmetry.

- eigen_quant\ :math:`J`\ .chk contain the bookkeeping information for the basis sets used: 

  the mapping between the multidimensional, multimode description of the product-form basis functions to a 1D basis set index.
  
- eigen_vectors\ :math:`J`\ _\ :math:`G`\ .chk contain the eigencoefficients written in direct unformatted form. 

  For each eigen_descr*chk there is an eigen_vector*chk file.   


