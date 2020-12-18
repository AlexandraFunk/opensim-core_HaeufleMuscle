NOTES to be read first:

This folder shall implement the Rockenfeller activation dynamics into the 
Opensim Standard Millard 2012 Equilibrium muscle. Therefore the calcMuscleDynamicsInfo
function from the millard implementation is overridden. 

To achieve this private functions from Millard 2012 needs to be copied to the new class, they are:
- calcFiberForce
- calcFiberStiffness
- calc_DFiberForceAT_DFiberLength
- calc_DFiberForceAT_DFiberLengthAT

The results of this muscle need to be rechecked at least against a Matlab implementation,
a second implementation into demoa would be even better!

Mike
