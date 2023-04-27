Molecules
#########
.. _molecules:

In this chapter a list is given of all of the molecules which have been studied using TROVE. This serves two purposes: as a summary of the development of TROVE over the years and as a reference for adding further molecules. As discussed in Chapter newmol_, setting up a new molecule in TROVE is fairly straightforward if a molecule of the same symmetry and structure has already been implemented.

The list of molecules are given in roughly chronological order with the relevant references. Details are given of the symmetry, basis sets used, coordinates used, the analytical forms of the PES and DMS and the frequency ranges of linelists if they were calculated.



Hydrogen sulfide, H\ :sub:`2` S
----------------------------

Symmetry: :math:`C_{2v}`

Coordinates: Linearized coordinates. :math:`\xi_1 = r_1^l - r_e`, :math:`\xi_2 = r_2^l - r_e` and :math:`\xi = \alpha^l - \alpha_e` .

Coordinate to expand kinetic energy: :math:`g_n = \xi_n (n=1,2,3)`.

Coordinates to expand Potential energy: :math:`f_n = 1 - \exp(-a(r_1^l - r_e))` :math:`(n = 1, 2)`, :math:`f_3 = \cos(\alpha^l) - \cos(\alpha_e)`

Primitive basis set: Numerov generated for all coordinates

Kinetic energy expansion order: 8

Potential expansion order: 10

Polyad scheme: :math:`P = 2(v_1 + v_2) + v_3 \leq 20`

Potential energy function: Morbid expansion of published potential.

Dipole moment surface expansion: N/A

Results: :math:`J = 0` vibrational energy levels up to 8000 cm\ :sup:`-1`.

.. Note:: Usually TROVE would not be used for a triatomic molecule as programs exist which have the 'exact' triatomic Hamiltonian (DVR3D for example [DVR3D]_). This molecule was used as a test case to demonstrate TROVE's capabilities in the original TROVE paper. The calculation was also repeated using HBJ theory, taking the bending coordinate to be the large amplitude coordinate and using normal coordinates.

Reference: [TROVE]_


Methyl cation, :math:`{\rm CH}_3^+`
-----------------------------------



Symmetry: :math:`D_{3h}`

Coordinates: Linearized coordinates. :math:`\xi_k = r_k^l - r_e` :math:`k = 1,2,3` for vibrational coordinates,
two symmetrized bending variables :math:`\xi_4 = S_{4a} = \frac{1}{\sqrt{6}} (2 \alpha_1^l - \alpha_2^l - \alpha_3^l)` and
:math:`\xi_5 = S_{4b} = \frac{1}{\sqrt{2}}(\alpha_2^l - \alpha_3^l)` and an out of plane vibration coordinate
:math:`\xi_6 = \rho = \frac{\mathbf{r_1} \cdot (\mathbf{r_2} \times \mathbf{r_3} )} {r_{1}r_{2}r_{3}} `. See paper for details.


Coordinate to expand kinetic energy: :math:`g_n = \xi_n (n=1-6)`

Coordinates to expand Potential energy: :math:`f_n = 1 - \exp(-a(\xi_n))` :math:`(n = 1, 3)` for stretching coordinates, :math:`f_4 = S_{4a}`, :math:`f_5 = S_{4b}` for two bending modes and :math:`f_6 = \rho`.

Primitive basis set: Numerov generated for all coordinates.

Kinetic energy expansion order: 6

Potential expansion order: 6

Polyad scheme: :math:`P = 1.5(v_1 + v_2 + v_3) + v_4 + v_5 + v_6 \leq 18`

Potential energy function: Published potential.

Dipole moment surface expansion: N/A

Results: :math:`J = 0` vibrational energy levels up to 6000 cm\ :sup:`-1`.

.. Note:: This was also used as a test example in the original TROVE paper. The coordinate scheme employed is similar to that for Ammonia (see below).

Reference: [TROVE]_



Carbon monoxide, CO
-------------------



Symmetry: :math:`C_{\infty V}`

Coordinates: r, bond coordinate between C and O.


Coordinate to expand kinetic energy: :math:`g_n = r`

Coordinates to expand Potential energy: Analytical and Morse (See paper).

Primitive basis set: Numerov generated for all coordinates.

Kinetic energy expansion order:

Potential expansion order:

Polyad scheme:

Potential energy function: Published empirical PEC. (REF)

Dipole moment surface expansion: N/A

Results: :math:`J = 0` vibrational energy levels up to 43000 cm\ :sup:`-1` (corresponding to :math:`v = 22`).

.. Note:: This was also used as a test example in the original TROVE paper. For diatomic molecules specialist programs
are of course recommended such as Duo [Duo]_. CO, like H\ :sub:`2`, CO is included only as a test case.

Reference: [TROVE]_


Ammonia, NH\ :sub:`3`
-------------------

Symmetry: :math:`C_{3v}`

Coordinates: Similar to those for :math:`{\rm CH}_3^+` but for sixth coordinate, :math:`xi_6 = \sin \rho_e - \sin \rho` where
:math:`\sin \rho = \frac{2}{\sqrt{3}} \sin\left[ (\alpha_1 + \alpha_2 + \alpha_3)/6) \right]`.

Coordinate to expand kinetic energy: :math:`g_n = \xi_n (n=1-6)`

Coordinates to expand Potential energy: Morse for stretching coordinates, angles themselves for bends.

Primitive basis set: Numerov generated for all coordinates.

Kinetic energy expansion order: 6

Potential expansion order: 8

Polyad scheme: For BYTe line list it is :math:`P = 2(v_1 + v_2 + v_3) + v_4 + v_5 + \frac{v_6}{2} \leq 14`

Potential energy function: Refinement of published potential [09YuBaYa]_.

Dipole moment surface expansion: For BYTe line list, an *ab initio* DMS was computed at the CCSD(T)/aug-cc-pVQZ level of
 theory [09YuBaYa]_.

Results:  Hot line list called BYTe. BYTe is applicable for temperatures up to 1500 K. It Comprises of 1138 323 351 transitions in the frequency range from 0 to 12 000 wavenumbers, constructed from 1373 897 energy levels below 18 000 wavenumbers having J values :math:`\le` 36.

.. Note:: Apart from BYTe, ammonia was used to develop TROVE itself, specifically for the J=0 contraction and refinement methods. The BYTe line list remains important for astronomical applications but will also soon be joined by an even more accurate line list from the work of Coles \textit{et al.} [10CoYuTe]_.

Reference:  [09YuBaYa]_, [11YuBaTe]_, [10CoYuTe]_.


Methane, CH\ :sub:`4`
--------------------

Symmetry: :math:`{T}_d`

Coordinates: Linearised coordinates. :math:`\xi_i = (r_i - r_e) \exp(-\beta(r_i - r_e)^2)` :math:`i = 1,4` for stretching coordinates.
 :math:`\xi_5 = \frac{1}{12}(2\alpha_{12} - \alpha_{13} - \alpha_{14} - \alpha_{23} - \alpha_{24} + 2\alpha_{34}`),
 :math:`\xi_6 = \frac{1}{2}(\alpha_{13} - \alpha_{14} - \alpha_{24} + \alpha_{24})`, :math:`\xi_7 = \frac{1}{\sqrt{2}}(\alpha_{24}
 - \alpha_{23})`, :math:`\xi_8 = \frac{1}{\sqrt{2}}(\alpha_{23} - \alpha_{14})` and :math:`\xi_9 = \frac{1}{\sqrt{2}}(\alpha_{34}
 - \alpha_{12})`. Where :math:`\alpha_{ij}` is the interbond angles. Also complimented by redundancy conditions (see paper).

Coordinate to expand kinetic energy: :math:`g_n = \xi_n (n=1-9)`, linearised coordinates.

Coordinates to expand Potential energy: :math:`f_n = 1 - \exp(-a(\xi_i^l))` :math:`(i = 1, 4)` for stretching coordinates and :math:`f_n = \xi_i^l` :math:`(i = 5, 9)` for bending coordinates.

Primitive basis set: Numerov generated for stretching coordinates, harmonic oscillator basis for bends.

Kinetic energy expansion order: 6

Potential expansion order: 8

Polyad scheme: :math:`P = 2(v_1 + v_2 + v_3 + v_4) + v_5 + v_6 + v_7 + v_8 + v_9 \leq 20` with caveats, see paper.

Potential energy function:  CCSD(T)-F12c/aug-cc-pVQZ-F12 + DK relativistic corrections *ab initio* data fit using polynomial of symmetrised coordinates given above. Refined using experimental :math:`J = 0, 4` data from HITRAN 2008 database.

Dipole moment surface expansion: CCSD(T)-F12c/aug-cc-pVTZ-F12 *ab initio* points fit using polynomial of symmetrised coordinates which is then expressed in symmetrised molecular bond (SMB) representation.

Results: 10to10 linelist complete for up to 1500 K. All states up to 18000 cm\ :sup:`-1` included, up to `J = 39`.

.. Note:: This describes the 10to10 calculation which was based on a previous calculation for lower frequencies. The high symmetry of methane meant special symmetry considerations are required. Details of this are given in the papers.

Reference: [13YuTeBa]_, [14YuJo]_.


Sulfur trioxide, SO\ :sub:`3`
----------------------------

Symmetry: :math:`D_{3h}`

Coordinates: As for ammonia.

Coordinate to expand kinetic energy: As for ammonia.

Coordinates to expand Potential energy: As for ammonia.

Primitive basis set: As for ammonia.

Kinetic energy expansion order: 6

Potential expansion order: 8

Polyad scheme: :math:`P = 2(n_1 + n_2 + n_3) + n_4 + n_5 + \frac{n_6}{2} \leq 18 `

Potential energy function: CCSD(T)-F12b/aug-cc-pVTZ-F12 + scalar relativistic corrections and DBOCs *ab initio* energies fitted to polynomial expansion of symmetrised coordinates. Refined using :math:`J \leq 5` experimental energies.

Dipole moment surface expansion: *ab initio* calculations at the same levels as for PES. Fitted using SMB
representation.

Results: Linelist complete up to 5000 cm\ :sup:`-1` for temperatures up to 800 K.

.. Note:: As SO\ :sub:`3` has a large moment of inertia, many :math:`J`s need to be included. Up to :math:`J = 130` was included for
a complete linelist at 800 K. For calculating :math:`J` this large, special procedures were used as discussed in the paper.

Reference: [16UnTeYu]_.


Hydrogen peroxide, H\ :sub:`2` O\ :sub:`2`
----------------------------------------

Symmetry: :math:`D_{2h}`. This is not the same as the point group of the molecule which is C\ :sub:`2`.

Coordinates: :math:`\xi_i = (x_i^l - x_i^e)` where :math:`i = 1, 6` are :math:`R`, :math:`r_1`, :math:`r_2`, :math:`\theta_1`, :math:`\theta
_2` and :math:`\tau`.

Coordinate to expand kinetic energy: :math:`g_n = \xi_n (n=1-6)`, linearised coordinates

Coordinates to expand Potential energy: :math:`f_n = 1 - \exp(-a_i(\xi_i^l))` :math:`(i = 1, 3)` for stretches and
:math:`f_n = \xi_i^l` :math:`(i = 4, 6)` for bending coordinates.

Primitive basis set: Numerov generated for all coordinates.

Kinetic energy expansion order: 6

Potential expansion order: 8

Polyad scheme: :math:`P = 4n_1 + 8(n_2 + n_3 + n_4 + n_5) +n_6 \leq 42`

Potential energy function: *ab initio* energies using CCSD(T)-F12b/aug-cc-pVNZ for N up to 7
for different parts of surface including DBO, relativistic, core-valence corrections fit to polynomial function
of coordinates. Refined to experimental energies for :math:`J \leq 4`.

Dipole moment surface expansion:  CCSD(T)-F12b/aug-cc-pV(T+d)Z fittied to polynomial function.


Results:  Linelist complete up to 6000 cm\ :sup:`-1`. Extended linelist up to 8000 cm\ :sup:`-1` with reduced completeness
at high temperatures.

.. Note:: The :math:`\tau` coordinate for this molecule adds complications to expansion of dipole, etc. See papers for details.

Reference: [15AlOvYu]_, [16AlPoOv]_.



Phosphine, PH\ :sub:`3`
----------------------

Symmetry: :math:`C_{3v}`

Coordinates: As for ammonia

Coordinate to expand kinetic energy: As for ammonia

Coordinates to expand Potential energy: As for ammonia

Primitive basis set: Numerov generated for all coordinates.

Kinetic energy expansion order: 6

Potential expansion order: 8

Polyad scheme: :math:`P = 2(s_1 + s_2 + s_3) + b_1 + b_2 + b_3 \leq 16` plus some additions, see paper.

Potential energy function:  CCSD(T)/aug-cc-pV(Q+d)Z) *ab initio* energies fitted to polynomial expansion.
Refined using HITRAN data up to :math:`J = 4`.

Dipole moment surface expansion: CCSD(T)/aug-cc-pVTZ *ab initio* dipole data fitted to polynomial expansion.


Results: SAlTY linelist, complete for up to 1500 K. All states up to 18000 cm\ :sup:`-1` included, up to :math:`J = 46`

.. Note:: For PH\ :sub:`3`, tunneling splitting via the umbrella motion may exist (as for NH\ :sub:`3`) may exist  but has yet to be detected [16SoYuTe]_.


References: [13SoYuTe]_, [15SoAlTe]_.



Formaldehyde, H\ :sub:`2`CO
--------------------------

Symmetry: :math:`C_{2v}`

Coordinates: :math:`\xi_i = (x_i^l - x_i^e)` where :math:`i = 1, 6` are :math:`r_1^l`, :math:`r_2^l`, :math:`r_3^l`, :math:`\theta_1^l`, :math:`\theta_2^l` and :math:`\tau`.

Coordinate to expand kinetic energy: :math:`g_n = \xi_n`, linearised.

Coordinates to expand Potential energy: :math:`f_n = 1 - \exp(-a_i(\xi_i^l))` :math:`(i = 1, 3)` for stretches, :math:`f_n = xi_i` :math:`(i = 4, 6)`
for bends.

Primitive basis set: Numerov generated for all coordinates.

Kinetic energy expansion order: 6

Potential expansion order: 8

Polyad scheme: :math:`P = 2(n_2 + n_3) + n_1 + n_4 + n_5 \leq 16` plus some additions, see paper.

Potential energy function:  CCSD(T)/aug-cc-pVQZ) *ab initio* energies fitted to polynomial expansion.
Refined using HITRAN data up to :math:`J = 5`.

Dipole moment surface expansion: CCSD(T)/aug-cc-pVQZ *ab initio* dipole data fitted to polynomial expansion.

Results: Linelist for temperatures up to 1500 K for transitions up to 10,000 cm\ :sup:`-1` and :math:`J = 70`.


Reference: [15AlOvPo]_.


Silane, SiH\ :sub:`4`
--------------------

Symmetry: :math:`T_d`

Coordinates: Linearised coordinates. As for methane.

Coordinate to expand kinetic energy: As for methane but with curvilinear coordinates.

Coordinates to expand Potential energy: As for methane.

Primitive basis set: Numerov generated for all coordinates.

Kinetic energy expansion order: 6

Potential expansion order: 8

Polyad scheme: :math:`P = 2(n_1 + n_2 n_3 + n_4) + n_5 + n_6 + n_7 + n_8 + n_9 \leq 12` plus some additions, see paper.

Potential energy function: CBS-F12 PES including extensive corrections, see paper. Fitted to polynomial expansion.
Refined using 1452 experimental energies up to :math:`J = 6`.

Dipole moment surface expansion: CCSD(T)/aug-cc-pVT(+d for Si)Z *ab initio* dipole data fitted to polynomial expansion.

Results: Linelist for temperatures up to 1200 K for transitions up to 5000 cm\ :sup:`-1` and :math:`J = 42`.


Reference: [17OwYuYa]_.



Methyl chloride, CH\ :sub:`3`Cl
-----------------------------

Symmetry: :math:`C_{3v}`

Coordinates:  :math:`\xi_k = r_k^l - r_e` :math:`k = 1,2,3,4` for vibrational coordinates,
:math:`\xi_i = \beta_i - \beta_e` :math:`i = 5,6,7` for bending coordinates, :math:`\xi_8 = \frac{1}{\sqrt{6}} (2 \tau_{23} -
 \tau_{13} - \tau_{12})` and :math:`\xi_9 = \frac{1}{2}(\tau_{13} - \tau_{12})`.

Coordinate to expand kinetic energy: :math:`g_n = \xi_n`, curvilinear coordinates used.

Coordinates to expand Potential energy: :math:`f_n = 1 - \exp(-a_i(\xi_i^l))` :math:`(i = 1, 4)` for stretches and
:math:`f_n = \xi_i^l` :math:`(i = 4, 9)` for bending coordinates.

Primitive basis set: Numerov generated for all coordinates.

Kinetic energy expansion order: 6

Potential expansion order: 8

Polyad scheme: :math:`P = n_1 + 2(n_2 + n_3 + n_4) + n_5 + n_6 + n_7 + n_8 + n_9 \leq 14` plus some additions, see paper.

Potential energy function: CBS-F12 PES including extensive corrections, see paper. Fitted to polynomial form.

Dipole moment surface expansion: CCSD(T)/aug-cc-pVQZ(+d for Cl) level of theory. Fitted to polynomial form.

Results: Line list applicable up to 1200 K.

.. Note:: Data for :math:`^{35}{\rm Cl}` and :math:`^{37}{\rm Cl} isotopologues.

Reference: [15OwYuTa]_,[18OwYaTe]_.


Ethylene, C\ :sub:`2`\ H\ :sub:`4`
----------------------------------

Symmetry: :math:`D_2h`

Coordinates: :math:`\xi_n = r_i-r_e` :math:`i=1,5` for stretches, :math:`\xi_n = \theta_i - \theta_e`  :math:`i = 1, 4` for bends,
:math:`\xi_10 = \pi - \beta_1`, :math:`\xi_11 = \beta_2 - \pi` for two :math:`\beta` H-C-H 'book type' angles and
:math:`\xi_12 = 2 \tau - \beta_1 + \beta_2` where :math:`\tau` is H-C-C-H dihedral angle.

Coordinate to expand kinetic energy: :math:`g_n = \xi_n`. Curvilinear coordinates.

Coordinates to expand Potential energy: Morse coordinates for stretches, other coordinates expanded as :math:`\xi` themselves.

Primitive basis set: Numerov generated for all coordinates.

Kinetic energy expansion order: 6

Potential expansion order: 8

Polyad scheme: :math:`P = n_1 + 2(n_2 + n_3 + n_4 + n_5) + n_6 + n_7 + n_8 + n_9 + n_{10} + n_{11} + n_{12} \leq 10` plus additions,
see paper.

Potential energy function: *ab initio* PES calculated at  CCSD(T)-F12b/cc-pVTZ-F12 level of theory. Fit to polynomial
form. Refined PES using HITRAN data for :math:`J=1-4` and other sources for vibrational band centres.

Dipole moment surface expansion: DMS calculated at CCSD(T)-F12b/aug-cc-pVTZ level of theory and fit to polynomial form with
appropriate axis system.

Results: Line list for 0-7000 cm\ :sup:`-1` up to :math:`J=78`. Applicable up to 700 K.

.. Note:: Largest molecule in TROVE so far. Special techniques developed to cope with such a large molecule.

Reference: [18MaYaTe]_.


Phosphorus trifluoride, PF\ :sub:`3`
----------------------------------

Symmetry: :math:`C_{3v}`

Coordinates: :math:`\xi_n = r_i - r_e` :math:`i=1,3` for stretching coordinates and :math:`\xi_n = \alpha_i - \alpha_e` :math:`i=1,3` for bends.

Coordinate to expand kinetic energy: :math:`g_n = \xi_n`. Linearised expansion.

Coordinates to expand Potential energy: Morse coordinates for stretches, bends expanded as :math:`\xi` themselves.

Primitive basis set: Numerov generated for all coordinates.

Kinetic energy expansion order: 6

Potential expansion order: 8

Polyad scheme: :math:`P = 2(n_1 + n_2 + n_3) + n_4 + n_5 + n_6 \leq 14`.

Potential energy function:  *ab initio* PES calculated at CCSD(T)-F12b/cc-pVTZ-f12 level of theory fitted using
polynomial expansion of symmetrized coordinates.

Dipole moment surface expansion: CCSD(T)/aug-cc-pVTZ *ab initio* dipole data fitted to polynomial expansion.

Results: Room temperature line list for up to :math:`J = 60`.

.. Note:: The room temperature line list for this molecule is not complete but could be easily extended using the methods applied
to SO\ :sub:`3` and :math:`{\rm C}_2{\rm H}_6`.

Reference: 19MaChYa]_.


References
----------



.. [DVR3D] J. Tennyson, M. A. Kostin, P. Barletta, G. J. Harris, O. L. Polyansky, J. Ramanlal, N. F. Zobov, Comput. Phys. Commun. 163, 85 (2004), DVR3D: a program suite for the calculation of rotation-vibration spectra of triatomic molecules.


.. [TROVE] S. N. Yurchenko, W. Thiel, P. Jensen, J. Mol. Spectrosc. 245, 126 (2007), Theoretical ROVibrational Energies (TROVE): A robust numerical approach to the calculation of rovibrational energies for polyatomic molecules.


.. [Duo] S. N. Yurchenko, L. Lodi, J. Tennyson, A. V. Stolyarov, Comput. Phys. Commun. 202, 262 (2016), Duo: A general program for calculating spectra of diatomic molecules.

.. [09YuBaYa] S. N. Yurchenko, R. J. Barber, A. Yachmenev, W. Thiel, P. Jensen, J. Tennyson, J. Phys. Chem. A 113, 11845 (2009), A variationally computed T=300 K line list for NH\ :sub:`3`.

.. [10CoYuTe] Phillip A Coles, Sergei N Yurchenko, and Jonathan Tennyson. Mon. Not. R. Astron. Soc., 490, 4638 - 4647 (2019), ExoMol molecular line lists XXXV: a rotation-vibration line list for hot ammonia.

.. [11YuBaTe]  S. N. Yurchenko, R. J. Barber, J. Tennyson, Mon. Not. R. Astron. Soc. 413, 1828 (2011), A variationally computed hot line list for NH\ :sub:`3`.

.. [13YuTeBa] S. N. Yurchenko, J. Tennyson, R. J. Barber, W. Thiel, J. Mol. Spectrosc. 291, 69 (2013), Vibrational transition moments of CH\ :sub:`4` from the first principles.

.. [14YuJo] Sergei N. Yurchenko and Jonathan Tennyson.  Mon. Not. R. Astron. Soc., 440, 1649-1661 (2014), ExoMol line lists - IV. The rotation-vibration spectrum of methane up to 1500 K.

.. [16UnTeYu] D. S. Underwood, J. Tennyson, S. N. Yurchenko, S. Clausen, A. Fateev, Mon. Not. R. Astron. Soc. 462, 4300 (2016), ExoMol line lists XVII: A line list for hot SO\ :sub:`3`.

.. [15AlOvYu] A. F. Al-Refaie, R. I. Ovsyannikov, O. L. Polyansky, S. N. Yurchenko, J. Tennyson, J. Mol. Spectrosc. 318, 84 (2015), A variationally calculated room temperature line-list for H\ :sub:`2`O\ :sub:`2`.

.. [16AlPoOv] A. F. Al-Refaie, O. L. Polyansky, R. I. Ovsyannikov, J. Tennyson, S. N. Yurchenko, Mon. Not. R. Astron. Soc. 461, 1012 (2016), ExoMol line lists XV: A hot line-list for hydrogen peroxide.

.. [13SoYuTe]  C. Sousa-Silva, S. N. Yurchenko, J. Tennyson, J. Mol. Spectrosc. 288, 28 (2013), A computed room temperature line list for phosphine.

.. [15SoAlTe] C. Sousa-Silva, A. F. Al-Refaie, J. Tennyson, S. N. Yurchenko, Mon. Not. R. Astron. Soc. 446, 2337 (2015), ExoMol line lists - VII. the rotation-vibration spectrum of phosphine up to 1500 K.

.. [16SoYuTe] C. Sousa-Silva, S. N. Yurchenko, J. Tennyson, J. Chem. Phys. 145, 091102 (2016), Tunelling splitting in the phosphine molecule.

.. [15AlOvPo] A. F. Al-Refaie, R. I. Ovsyannikov, O. L. Polyansky, S. N. Yurchenko, J. Tennyson, J. Mol. Spectrosc. 318, 84 (2015), A variationally calculated room temperature line-list for H\ :sub:`2`CO.

.. [17OwYuYa] A. Owens, S. N. Yurchenko, A. Yachmenev, W. Thiel, J. Tennyson, Mon. Not. R. Astron. Soc. 471, 5025 (2017), ExoMol molecular line lists XXII. The rotation-vibration spectrum of silane up to 1200K.

.. [18OwYaTe] A. Owens, A. Yachmenev, J. Tennyson, W. Thiel, S. N. Yurchenko, Mon. Not. R. Astron. Soc. 479, 3002 (2018), ExoMol Molecular line lists XXIX: The rotation-vibration spectrum of methyl chloride up to 1200 K.

.. [15OwYuTa] A. Owens, S. N. Yurchenko, A. Yachmenev, J. Tennyson, W. Thiel, J. Chem. Phys. 142, 244306 (2015), Accurate ab initio vibrational energies of methyl chloride.

.. [18MaYaTe] B. P. Mant, A. Yachmenev, J. Tennyson, S. N. Yurchenko, Mon. Not. R. Astron. Soc. 478, 3220 (2018), ExoMol molecular line lists - XXVII: spectra of C\ :sub:`2`H\ :sub:`4`.

.. [19MaChYa] Barry P. Mant, Katy L. Chubb, Andrey Yachmenev, Jonathan Tennyson, and Sergei N. Yurchenko, Mol. Phys., 118, e1581951 (2019), The infrared spectrum of PF\ :sub:`3` and analysis of rotational energy clustering effect.
