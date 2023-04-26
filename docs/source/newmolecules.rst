New Molecules
=============
.. _newmol:

In this chapter details will be given on how to set up a new molecule in TROVE. 
There are two approaches to this. If a molecule of the same symmetry has already been set up in TROVE then it is
relatively straightforward to add a similar molecule. This is one of the strengths of the TROVE approach. 
If a completely new kind of molecular structure is needed then
more work is required. The former case will be discussed first. 

If a new isotopologue is required then a combination of the two approaches will be required, for example the PES may 
be of the same form but the symmetry will change. 


Adding a Similar Symmetry Molecule
----------------------------------

Adding a molecule of similar symmetry is relatively straightforward. It is possible that no new TROVE files for defining 
the potential, etc will be required. An example of this is adding arsine (AsH:sub:`3`) when
phosphine (PH:sub:`3`) is already set up in TROVE. Both molecules have trigonal pyramidal structure and belong to the 
:math:`C_{3v}` point group. 

In the TROVE input file the molecule is defined in various places as was discussed in Chapter Quickstart_. If a 
new molecule of the same symmetry as a previously defined molecule is to be added, then the symmetry group 
defined by the keyword ``SYMGROUP`` should be kept the same.

``TRANSFORM``, which determines the subroutine to transform coordinates from Z-matrix to those for computations, 
should also be left the same. If a new type of coordinates is required then this needs to be changed either to another 
existing transform or a new one. Setting this up is described below. 

``MOLTYPE`` specifies which routines should be used to define the molecule's geometry in Cartesian coordinates,
how the rotational symmetry is implemented and other pointers to symmetry transformations. Again it should be kept the same 
as used for a previously defined molecule. 

``MOLECULE`` is an optional keyword which does not necessarily need to be changed. 

The ``ZMAT`` does need to be changed for a new molecule type but probably only the atomic masses. For example
on going from PH:sub:`3` to AsH:sub:`3`, only the central atom mass changes. Similarly the ``Equilibrium`` block needs to 
be changed with the new molecule's bond lengths and angles. The order of these will be the same as the 
previously defined molecule if the same Z-matrix is used.

The ``SPECPARAM`` block depends on the specific potential used and may need to be changed as required. 

If a new molecule is being set up with a similar structure to a previously defined molecule then it may be that the 
same form of function can be used to fit the PES. If this is the case then in the ``POTEN`` block, the 
``POT_TYPE`` should be kept the same but the list of parameters (equilibrium geometry values, Morse parameters and
expansion coefficients, etc) should be changed to the values for the new molecule. Similarly, if the previously defined
DMS is suitable tthen in the ``EXTERNAL`` block, the ``DMS_TYPE`` can be kept the same and parameters changed to those
of the new molecule. 

If new functional forms for the PES or DMS are required then these need to be defined as described below.


Setting Up a New Molecule
-------------------------

Setting up a new molecule for TROVE requires the user to modify or add new subroutines to TROVE which will be added to the
list of files to be compiled with the main program. Before the specific details of these are provided a description 
of what is required will be given.

As discussed in Chapter theory_, TROVE expands the Hamiltonian in terms of molecule-fixed internal coordinates.
To do this, transformations between various coordinates are required. The molecule's geometry (including the equilibrium
geometry) is defined in the input file via the Z-matrix. This gives a straightforward and intuitive 
way to define a molecule's geometry. Z-matrix coordinates are usually not used to expand the kinetic or potential
energy however as they do not take into account symmetry, are not of Morse form, etc. Transformation to the working coordinates
and potential is given in a mol and pot file, usually one of each is defined per molecule. 

The mol file contains one or multiple subroutines for transforming from Z-matrix coordinates to `TROVE coordinates' 
(the coordinates used to expand the kinetic energy and re-expand the potential energy) and back again. Choosing which
coordinates to use is part of the problem of setting up a new molecule. Various options are sometimes needed to assess 
the best choice. Criteria for choosing coordinates are: symmetry considerations, ease of use, numerical stability. etc.
The mol file also contains a routine for expressing the molecules equilibrium geometry in Cartesian coordinates. Another 
important part of the mol file is a routine which expresses how the TROVE coordinates transform under the symmetry which has
been chosen to work in. 

The pot file contains the potential energy surface function for the molecule of interest. This file should include the potential
function itself and how parameters are defined from the input file. TROVE will call the potential subroutine using
Z-matrix geometries and the pot file gives the transformation from these coordinates into those used for the potential.
This file typically also contains the dipole moment surface function. The dipole is set up in a similar way to the 
potential but also typically Cartesian vectors expressed in body-fixed coordinates need to be set up.

As well as the mol and .pot files, the symmetry of a new molecular structure may also need to be set up. There are already 
many molecular symmetry groups set up in TROVE but the particular one for the molecule of interest may be missing. The 
symmetries are collected together in the symmetry.f90 file. 

A final part of setting up a new molecule is describing the rotational symmetry. 

A detailed discussion of all of these files and examples are discussed below. TROVE is written in Fortran90 and in this 
chapter knowledge of that programming language shall be assumed (although the main points can be appreciated without it).

The mol File
^^^^^^^^^^^^

The mol file controls the coordinates used for a molecule in TROVE. Many such files have been written and when implementing
a new mol file it is recommended to follow a similar style to these. Specific details such as which routines need to be 
set to ``Public`` and which modules are to be used can be obtained from previous mol files. The PF:sub:`3` molecule will be used
 to illustrate  each section of the mol file.

The details of the TROVE coordinates and transformations to Z-matrix coordinates are given in the 
subroutine ``ML_coordinate_transform_XXX((src,ndst,direct) result (dst)`` (where XXX is some molecule or molecule type). 
The arguments to this subroutine are as follows.
``src`` (`source') is a vector of coordinates. These can either be Z-matrix or TROVE coordinates depending on 
which part of the main program is calling the subroutine. ``dst`` (`destination') are again either the Z-matrix or 
TROVE coordinates but will be the result of a transform from one to the other (this is why ``dst`` is
labelled ``result``). ``direct`` is a logical (true of false) variable to used to choose if the subroutine is 
going from Z-matrix to TROVE coordinates or vice-versa.

As is standard in Fortran, any variables required in the subroutine are then declared. A typical mol file has 
multiple IF statements to choose which transform and coordinates to use. As mentioned, many may be set up as 
some will work better for specific applications/molecules. 

The PF:sub:`3` molecule is of the generic type XY:sub:`3` and the mol file used is  ``mol_xy3.f90``. The first part of the
transform subroutine is 
::
     
     if (verbose>=5) write(out,"('ML_coordinate_transform_XY3/start')")
     !
     if (direct) then
     !
     dsrc(:) = src(:) - molec%local_eq(:)
     !
     else
     !
     dsrc(:) = src(:)
     !
     endif
     !
     nsrc = size(src)
     
This will print out the message if the ``verbose`` value is $>5$. Next the value of ``direct`` is checked. If true
then the molecule's equilibrium parameters (defined in a global vector from the input file) are subtracted from the
``src``. This is for Z-matrix to TROVE. Otherwise, the ``src`` vector is transferred to ``dsrc``. 

After this initial step many different choices of coordinates and transforms are defined. From Chapter Quickstart_ the PF:sub:`3` example was defined using
::
     
     dstep            0.01
     COORDS           linear
     TRANSFORM        r-alpha
     MOLTYPE          XY3
     MOLECULE         PF3
     REFER-CONF       RIGID
     
The ``MOLTYPE`` keyword selected the  ``mol_xy3.f90`` file. The specific coordinate transform to use is given by the 
``TRANSFORM`` keyword and is ``r-alpha``. This corresponds to one of the options in the mol file. The option is 
selected as
::
     
     case('R-ALPHA')
     !
     if (size(src)/=6) then
       write(out,"('MLcoordinate_transform_func: r-alpha  works only with 6 coords')")
       stop 'MLcoordinate_transform_func: r-alpha  works only with 6 coords'
       endif
       !
     if (direct) then
       !
       dst(1:3) = dsrc(1:3)
       dst(6) = dsrc(4)
       dst(5) = dsrc(5)
       dst(4) = dsrc(6)
       !
     else ! not direct
       !
       dst(1:3) = dsrc(1:3)+molec%local_eq(1:3)
       dst(6) = dsrc(4)+molec%local_eq(4)
       dst(5) = dsrc(5)+molec%local_eq(5)
       dst(4) = dsrc(6)+molec%local_eq(6)
       !
     endif
     
``case`` chooses the transform to use. There is then a check of how many coordinates are used. This routine only
works with 6 (other choices make use of extra redundant coordinates). ``direct`` is then used to check to which 
coordinates are being transformed. For Z-matrix to TROVE, the coordinates are taken directly from ``dsrc`` 
(as the equilibrium coordinates
were already subtracted at the start of the routine). If TROVE to Z-matrix, equilibrium coordinates are added to the TROVE
coordinates to get back to the Z-matrix values.

This is a very simple transformation but illustrates the idea. Other molecules have more complicated coordinates which
usually requires the application of more geometry transforms/trigonometry etc and symmetrised coordinates may be introduced.


The symmetry properties of the TROVE coordinates used is defined in the subroutine 
``ML_symmetry_transformation_XXX(ioper,nmodes,src,dst)``. The subroutine is used to define how the coordinates of the 
molecule permute into each other with a given symmetry operation.  The arguments to this subroutine are: ``ioper`` which is an integer do choose a symmetry operation, ``nmodes`` which is the number of vibrational modes and ``src``
and ``dst`` which are the coordinates before and after the symmetry operation. 

The symmetry group and coordinates used are chosen using ``case`` statements similar to the transform subroutine. These
are defined in the input file. For each symmetry operation the ``dst`` coordinates should be defined in terms of the 
initial ``src`` coordinates. This may involve introducing normalisation constants or other variables as needed. 

For PF:sub:`3` the symmetry transforms are defined in  ``ML_symmetry_transformation_XY3(ioper,nmodes,src,dst)``. The subroutine
starts by performing checks on the number of modes. The symmetry group is then chosen as
::
     
     select case(trim(molec%symmetry))
     case default
        write (out,"('ML_symmetry_transformation_XY3: symmetry ',a,' unknown')") 
        trim(molec%symmetry)
       stop 'ML_symmetry_transformation_XY3 - bad symm. type'
     case('C3V','C3V(M)')
     
where both ``C3V`` and ``C3V(M)`` can be used in the input file. As there are many TROVE coordinates defined for  XY:sub:`3` molecules, further ``case`` selections are required (if for a given molecule only one type of TROVE coordinates
has been set up then no further selects are necessary). For the ``r-alpha`` example the symmetry is defined by
::
     
     select case(trim(molec%coords_transform))
     !
     !
     case('R-ALPHA')
     !
     select case(ioper)
     !
     case (1) ! identity
     !
     dst = src
     !
     case (3) ! (132)
     !
     !dst(1) = src(2)
     !dst(2) = src(3)
     !dst(3) = src(1)
     !dst(4) = src(5)
     !dst(5) = src(6)
     !dst(6) = src(4)
     ...
     ...
     
Once the ``R-ALPHA`` coordinates are chosen, further ``case`` selects each symmetry operation. For the identity, $E$ 
operation, no change is required and so ``dst`` = ``src``. Here, case 3 corresponds to the operation (132) and the
bond lengths and angles are changed accordingly. The 4 other operations for this group have similar transforms. 


The centre of mass of the molecule in Cartesian coordinates is defined in the subroutine 
`` ML_b0_XXX(Npoints,Natoms,b0,rho_i,rho_ref,rho_borders)``. ``Natoms`` is the number of atoms and
``b0`` is a matrix containing the Cartesian coordinates of the atoms at the molecule's equilibrium geometry. The 
other subroutine arguments are optional and are for defining multiple geometries. This is needed if HBJ theory
is being used for a large amplitude coordinate. 

For PF:sub:`3` the subroutine is ``ML_b0_XY3``. This routine starts by performing checks to see if the number of 
atoms, equilibrium coordinates and atomic masses are consistent for an XY:sub:`3` molecule. Coordinates are then defined from
the input file equilibrium block as
::
     
     re14 = molec%req(1)
     alpha = molec%alphaeq(1)
     rho = pi-asin(2.0_ark/sqrt(3.0_ark)*sin(alpha/2.0_ark))
     
Using these coordinates the ``b0`` matrix is filled in with the Cartesian coordinates of the atoms
::
     
     cosr = cos(rho)
     sinr = sin(rho)
     !
     b0(2,1,0) = re14*sinr
     b0(2,2,0) = 0
     b0(2,3,0) = mX*re14*cosr/(Mtotal+mX)
     b0(3,1,0) = -re14*sinr/2.0_ark
     b0(3,2,0) = sqrt(3.0_ark)*re14*sinr/2.0_ark
     b0(3,3,0) = mX*re14*cosr/(Mtotal+mX)
     b0(4,1,0) = -re14*sinr/2.0_ark
     b0(4,2,0) = -sqrt(3.0_ark)*re14*sinr/2.0_ark
     b0(4,3,0) = mX*re14*cosr/(Mtotal+mX)
     b0(1,1,0) = 0
     b0(1,2,0) = 0
     b0(1,3,0) = -Mtotal*re14*cosr/(Mtotal+mX)
     

In this case ``b0`` has been defined explicitly with respect to the centre of mass of the molecule. If this is 
not the case then the centre of mass can be found using a subroutine. This step is part of the XY:sub:`3` subroutine as
::
     
     if (any(molec%AtomMasses(2:4)/=mH1)) then
     !
     do n = 1,3
     CM_shift = sum(b0(:,n,0)*molec%AtomMasses(:))/sum(molec%AtomMasses(:))
     b0(:,n,0) = b0(:,n,0) - CM_shift
     enddo
     


If the molecule contains a non-rigid degree of freedom (for example, the umbrella motion in NH:sub:`3`) then HBJ theory is used
as discussed in Chapter theory_. In this case TROVE expands the Hamiltonian on a grid of geometries along the 
non-rigid degree of freedom. The other arguments to the subroutine then come into play. ``Npoints`` is the number of 
points the non-rigid degree of freedom is split into, chosen in the ``BASIS`` block of the input file. ``rho_i`` 
is the value of the non-rigid coordinate for that ``npoint``. ``rho_ref`` and ``rho_borders`` are the reference
geometry (usually at equilibrium) and the ends of the grid along the non-rigid coordinate.

The array which contains the Cartesian coordinates, ``b0`` is of size ``(Natoms,3,Npoints)``. For rigid molecules, 
``Npoints`` = 0 and only the equilibrium geometry is necessary. For non-rigid, the coordinates of each atom are required
at each point along the non-rigid coordinate. A loop over  ``Npoints`` is required and the way the other rigid 
coordinates change at each ``rho_i`` is given. The mol file for NH:sub:`3` or H$_2$O$_2$ shows examples of this. 
Ideally the rigid coordinates should be set to change along the least energy path. Quantum chemistry programs such as 
MOLPRO can be used to find this where a geometry optimisation is carried out at each step. Alternatively
it can be done `by hand' from the PES.



A final part of the mol file which needs to be set up is the ``ML_rotsymmetry_XXX`` subroutine which defines
the rotational symmetry.

The pot File
^^^^^^^^^^^^

The pot file is used to define potential energy surfaces in TROVE. Although TROVE re-expands the PES in whichever 
coordinates have been chosen in the mol file (see Chapter theory_, the program needs the potential energy function as part 
of this process. As with the mol file the pot file can make use of parameters defined in the input file.

A typical pot file contains multiple PES functions which return the energy for a given geometry. For a given molecule
class many functions may be implemented to test different PESs or compare against functions given 
in the literature. The choice of PES is defined in the input file.

Each PES function is initiated by 
::
     
     function MLpoten_xxx(ncoords,natoms,local,xyz,force) result(f). 
     
The function arguments are as follows. ``ncoords`` and ``natoms`` are the number of vibrational coordinates and atoms respectively.
``local`` is the molecule's coordinates given in Z-matrix form as defined in the input file. ``xyz`` is a matrix
of atomic positions in Cartesian coordinates. ``force`` is a list of parameters for the PES defined in the input. The
energy at a given coordinate is the output (result) of the function, ``f``.  

For the PF:sub:`3` molecule the pot file is ``pot_xy3.f90``. This file contains multiple PES and DMS functions. From the PF:sub:`3`
example the PES is chosen in the input file as `` MLpoten_xy3_morbid_10``. This function starts by defining equilibrium
parameters from the input file and coordinates from ``local``. The specific choice for the ``r-alpha`` coordinate
transform is not given by a ``case`` (unlike others in the function) but instead by the specifics of the coordinates
::
     

     elseif (size(local)==6.and.molec%Ndihedrals==0) then
     !
     alpha3 = local(4)
     alpha2 = local(5)
     alpha1 = local(6)
     !
     tau = sqrt(1.0_ark-cos(alpha1)**2-cos(alpha2)**2-cos(alpha3)**2 &
                        +2.0_ark*cos(alpha1)*cos(alpha2)*cos(alpha3) )
     
as there is no dihedral angles for the ``r-alpha`` choice. After this the coordinates are transformed into those of the
PES used and a separate function for the PES called. Up to this point the function has been to transform to these coordinates
from whichever Z-matrix coordinates were specified.
::
     
     y1=1.0_ark-exp(-aa1*(r14-re14))
     y2=1.0_ark-exp(-aa1*(r24-re14))
     y3=1.0_ark-exp(-aa1*(r34-re14))
     !
     y4=(2.0_ark*alpha1-alpha2-alpha3)/sqrt(6.0_ark)
     y5=(alpha2-alpha3)/sqrt(2.0_ark)
     !
     f = poten_xy3_morbid_10(y1,y2,y3,y4,y5,coro,force)
     

The function ``poten_xy3_morbid_10`` itself is the PES function and uses the coordinates ``y1-y5`` along 
with the parameters in ``force``. The function is rather large and can be viewed in the pot file. 
The function is a sum of symmetrised combinations of the coordinates raised to powers
and multiplied by the relevant expansion parameters. These expansion are usually not all programmed by hand but 
obtained from symbolic mathematical software such as Mathematica or Python.

Rather than explicitly give all the symmetrised expansion coordinates in a PES routine, another approach is to 
do the symmetry `on the fly'. This means to apply the symmetry operations to coordinates by making use of the 
symmetry operation matrices for the group. This method is used in TROVE for the C:sub:`2`H:sub:`4` molecule. In the pot file this
is specified as
::
     
     f = 0
     !
     do ioper = 1,12
     !
       call ML_symmetry_transformation_XY3_II(ioper,xi,chi(:,ioper),18)
     !
     enddo
     !
     do i = 6, molec%parmax
       ipower(1:18) = molec%pot_ind(1:18,i)
       term = 0
         do ioper = 1,12
           term = term + product(chi(1:18,ioper)**ipower(1:18))
         end do
       term = term/12.0_ark
       f = f + term*force(i)
     end do
     
This starts by calling a symmetry transform subroutine (similar to that in the mol file discussed above) for each 
symmetry operation (12 in this case). All permutations are stored in the ``chi`` matrix.  The parameters 
of the potential are then looped over. The power to which each coordinate is raise is extracted from the 
list given in the input file (recall that parameters can be given as a simple list or including the powers, see Chapter Quickstart_. The symmetries 
are then looped over and each permutation raised by that power. The division by 12 is then applied to match how the 
PES was fit. Finally the relevant parameter multiplies the geometry term and then another loop over then next parameter
is started.

This approach guarantees that the symmetry of the molecule is taken into account. For example, if a C-H bond length was varied
then all other permutations are taken into account so that all C-H stretches are equivalent. 


The best way of setting up the pot file is molecule dependent. Many options are possible, as long as the energy is returned 
for a certain geometry. Many pot files have already been set up in TROVE, some with multiple choices per molecule type. These
can be referred to for more details of the procedure or used as a starting point for new potentials.








