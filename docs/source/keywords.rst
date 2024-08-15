TROVE Keywords
==============


POTENTIAL_SIMPLE

``IRON-OUT``: the card to switch on an automatic smoothing of all expansion terms of the PEF, DMF, KEO and external field when expanded around a non-rigid reference configuration. TROVE does not use this feature by default. It can however recommend to use it in the case of to large errors in the derivatives of these fields. The card needs to be place anywhere in the main body of the step 1 input. 

``verbose``: Verbosity level


``sparse``:


``exp_coeff_thresh``
^^^^^^^^^^^^^^^^^^^^

Expansion coefficients of  the field checkpoints ``kinetic.chk``, ``potential.chk`` and ``external.chk`` that are smaller by magnitude than this threshold are not included in the corresponding checkpoint.

