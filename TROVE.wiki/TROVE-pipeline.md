
 ## The typical TROVE intensity project consists of the following steps:

1. Expansion of the Hamiltonian operator (generating kinetic and potential energy expansion coefficients numerically on-the-fly) as well as of any “external” function (e.g., dipole moment, polarizability, spin–spin coupling, or any other property; PES correction used in the refinement process);
2. Numerov–Cooley solution of the 1D Schrodinger equations;
3. Eigen-solutions of the reduced Hamiltonian problems;
4. Symmetrization of the contracted eigenfunctions from Step 3 and construction of the symmetry-adapted vibrational basis set;
5. Calculation of the vibrational matrix of the Hamiltonian operator as well as external functions (e.g., dipole) when required;
6. Diagonalizaitons of the J=0 Hamiltonian matrices for each irreducible representation in question;
7. Conversion of the primitive basis set representation (vibrational matrix elements from Step 5) to the J=0 representation;
8. Construction of the symmetry-adapted ro-vibrational basis set as a direct product of the J=0 eigenfunctions and rigid rotor wavefunctions;
9. Construction of the ro-vibrational Hamiltonian matrices for each J>0 and irreducible representation Gamma;
10. Diagonalization of the Hamiltonian matrices and storing eigenvectors for the postprocessing (e.g., intensity calculations) if necessary;
11. For the intensity calculations (line list production), all pairs of the ro-vibrational eigenvectors (bra and ket) from Step 10 (subject to the selection rules as well as to the energy, frequency and J thresholds) are cross-correlated with the dipole moment XYZ components in the laboratory-fixed frame via a vector-matrix-vector product, where the body-fixed xyz components of the dipole moment from Step 5 are transformed to the XYZ-frame using the Wigner-matrices.
