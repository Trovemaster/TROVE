Kinetic energy operators
========================

There are several options for kinetic energy operators (KEO) in TROVE: 

- **Linearised KEO**: Automatically constructed KEO using the Eckart frame (rigid) or Eckart/Saywitz frame ((1D non-rigid) as a Taylor expansion in terms of linearised coordinates;
- **Curvilinear KEO**: Preprogrammed analytic KEO (exact or non-exact) using curvilinear (valence) coordinates and different frames;
- **External Numerical Taylor-type KEO**: Externally generated KEO as a Taylor-type numerical expansion in terms of any user-defined coordinates;
- **External sum-of-products KEO**: Externally generated KEO as a generic sum-of-products expansion via the TROVE KEO builder.


Linearised KEO
--------------

This is a standard KEO constructed using a general black-boxed algorithm applicable for an arbitrary polyatomic molecule. It has been introduced in the original TROVE paper [TROVE]_ and also explained in Chapter `Theory <../theory.html>`_. 




